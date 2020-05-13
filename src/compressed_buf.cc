//-----------------------------------------------------------------------------
/*!
 * \file   compressed_buf.cc
 * \author Wayne Joubert
 * \date   Tue May  5 20:11:12 EDT 2020
 * \brief  Mirrored buffer allowing for compression.
 * \note   Copyright (C) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#if defined COMET_USE_CUDA
#  include "cub.cuh"
#elif defined COMET_USE_HIP
// TODO
#endif

#include "env.hh"
#include "mirrored_buf.hh"
#include "compressed_buf.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Constructor for CompressedBuf class.

CompressedBuf::CompressedBuf(MirroredBuf& buf, CEnv& env)
  : env_(env)
  , buf_(&buf)
  , num_nonzeros_buf_(env)
  , num_runs_buf_(env)
  , keys_buf_(env)
  , lengths_buf_(env)
  , keys_alias_buf_(env)
  , lengths_alias_buf_(env)
  , reduce_workspace_buf_(env)
  , rle_workspace_buf_(env)
    // NOTE: assume this is the largest size to ever be used.
  , length_max_(buf_->dim0 * buf_->dim1 * NUM_VALUES_PER_METRIC)
  , num_nonzeros_approx_(0)
  , do_compress_(false)
  , num_runs_(0)
  , state_(State::IDLE) {
  //, is_open_(false)
  //, read_ptr_(0) {

  if (!can_compress_())
    return;

  COMET_INSIST(env_.metric_format() == METRIC_FORMAT);

  num_nonzeros_buf_.allocate(1, 1, sizeof(MFTTypeIn));

  keys_buf_.allocate(length_max_, 1, sizeof(MFTTypeIn));

  lengths_buf_.allocate(length_max_, 1, sizeof(size_t));

  num_runs_buf_.allocate(1, 1, sizeof(size_t));

  // Anticipate size--usually seems to be = 1.
  reduce_workspace_buf_.allocate(1, 1, sizeof(char));
}

//-----------------------------------------------------------------------------
/// \brief Replace the base underlying MirroredBuf object.

void CompressedBuf::attach(MirroredBuf& buf) {
  COMET_INSIST(State::IDLE == state_);

  buf_ = &buf;
  COMET_INSIST(length_() <= length_max_);
}

//-----------------------------------------------------------------------------
/// \brief Compute number of nonzero MFTTypeIn values in buffer device memory.

void CompressedBuf::compute_num_nonzeros_() {
  COMET_INSIST(can_compress_());

# if defined COMET_USE_CUDA

    // see https://nvlabs.github.io/cub/structcub_1_1_device_reduce.html

    ReductionOp reduction_op;
    size_t temp_storage_bytes = 0;

    // Get amount of temp storage needed.

    const MFTTypeIn initial_value = ReductionOp::INITIAL_VALUE;

    cub::DeviceReduce::Reduce(NULL,
      temp_storage_bytes, (MFTTypeIn*)buf_->d, (MFTTypeIn*)num_nonzeros_buf_.d,
      length_(), reduction_op, initial_value, env_.stream_compute());
    num_nonzeros_buf_.from_accel(env_.stream_compute());

    // Re/allocate temp storage.

    if (temp_storage_bytes > reduce_workspace_buf_.size()) {
      reduce_workspace_buf_.deallocate();
      // Allocate slightly extra to allow some future growth.
      const double fuzz = 0;
      const size_t alloc_size = temp_storage_bytes * (1+fuzz);
      reduce_workspace_buf_.allocate(alloc_size, 1, sizeof(char));
    }

    // Perform reduction.

    cub::DeviceReduce::Reduce((Workspace_t*)reduce_workspace_buf_.d,
      temp_storage_bytes, (MFTTypeIn*)buf_->d, (MFTTypeIn*)num_nonzeros_buf_.d,
      length_(), reduction_op, initial_value, env_.stream_compute());

    // Retrieve count.
    // NOTE: this value may be approximate, since reduction is over floats.

    num_nonzeros_buf_.from_accel(env_.stream_compute());

    num_nonzeros_approx_ = (size_t)
      -(num_nonzeros_buf_.elt_const<MFTTypeIn>(0, 0) - initial_value);

# elif defined COMET_USE_HIP

  // TODO

# endif // COMET_USE_CUDA || COMET_USE_HIP
}

//-----------------------------------------------------------------------------
/// \brief Compress the underlying buffer if it meets conditions.

void CompressedBuf::compress() {

  if (!can_compress_())
    return;

  compute_num_nonzeros_();

  do_compress_ = num_nonzeros_approx_ <= compress_threshold_() * length_();

  if (do_compress_) {

    COMET_INSIST(State::IDLE == state_);

#   if defined COMET_USE_CUDA
      //https://nvlabs.github.io/cub/structcub_1_1_device_run_length_encode.html
      size_t temp_storage_bytes = 0;

      // Get amount of temp storage needed.

      cub::DeviceRunLengthEncode::Encode(NULL,
        temp_storage_bytes, (MFTTypeIn*)buf_->d, (MFTTypeIn*)keys_buf_.d,
        (size_t*)lengths_buf_.d, (size_t*)num_runs_buf_.d, length_(),
        env_.stream_compute());
      num_runs_buf_.from_accel(env_.stream_compute());

      // Re/allocate temp storage.

      if (temp_storage_bytes > rle_workspace_buf_.size()) {
        rle_workspace_buf_.deallocate();
        // Allocate slightly extra to allow some future growth.
        const double fuzz = 0;
        const size_t alloc_size = temp_storage_bytes * (1+fuzz);
        rle_workspace_buf_.allocate(alloc_size, 1, sizeof(char));
      }

      // Perform RLE.

      cub::DeviceRunLengthEncode::Encode((Workspace_t*)rle_workspace_buf_.d,
        temp_storage_bytes, (MFTTypeIn*)buf_->d, (MFTTypeIn*)keys_buf_.d,
        (size_t*)lengths_buf_.d, (size_t*)num_runs_buf_.d, length_(),
        env_.stream_compute());

      // Retrieve size.

      num_runs_buf_.from_accel(env_.stream_compute());

      num_runs_ = num_runs_buf_.elt_const<size_t>(0, 0);

      // Create right-sized buf aliases to later capture run lengths/keys.

      keys_alias_buf_.allocate(keys_buf_, num_runs_);
      lengths_alias_buf_.allocate(lengths_buf_, num_runs_);

#   elif defined COMET_USE_HIP

  // TODO

#   endif // COMET_USE_CUDA || COMET_USE_HIP

    state_ = State::COMPRESSED;

  } // if (do_compress_)
}

//-----------------------------------------------------------------------------

void CompressedBuf::from_accel_start() {

  if (do_compress_) {
    COMET_INSIST(State::COMPRESSED == state_);
    keys_alias_buf_.from_accel_start();
    lengths_alias_buf_.from_accel_start();
    state_ = State::TRANSFER_STARTED;
    //buf_->from_accel_start();

  } else {
    buf_->from_accel_start();
  }
}

//-----------------------------------------------------------------------------

void CompressedBuf::from_accel_wait() {

  if (do_compress_) {
    COMET_INSIST(State::TRANSFER_STARTED == state_);
    keys_alias_buf_.from_accel_wait();
    lengths_alias_buf_.from_accel_wait();

#if 0
    for (size_t i=0; i<num_runs_; ++i) {
      printf("%i %i %f\n", (int)i,
        (int)lengths_alias_buf_.elt_const<size_t>(i, 0),
        (double)keys_alias_buf_.elt_const<MFTTypeIn>(i, 0) ); 
    }
#endif

    reader_.init();

    state_ = State::IDLE;
    //buf_->from_accel_wait();

#if 0
    for (size_t i=0; i<64; ++i) {
      if ((double)(((MFTTypeIn*)(buf_->h))[i]))
        printf("%i %f\n", (int)i,
          (double)(((MFTTypeIn*)(buf_->h))[i]));
    }
#endif

  } else {

    buf_->from_accel_wait();

  }
}

//-----------------------------------------------------------------------------

void CompressedBuf::lock_h() {

  if (do_compress_) {
    keys_alias_buf_.lock_h();
    lengths_alias_buf_.lock_h();
  } else {
    buf_->lock_h();
  }
}

//-----------------------------------------------------------------------------

void CompressedBuf::unlock_h() {

  if (do_compress_) {
    keys_alias_buf_.unlock_h();
    lengths_alias_buf_.unlock_h();
  } else {
    buf_->unlock_h();
  }
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
