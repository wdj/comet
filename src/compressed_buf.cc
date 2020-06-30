//-----------------------------------------------------------------------------
/*!
 * \file   compressed_buf.cc
 * \author Wayne Joubert
 * \date   Tue May  5 20:11:12 EDT 2020
 * \brief  Mirrored buffer allowing for compression.
 */
//-----------------------------------------------------------------------------
/*-----------------------------------------------------------------------------

Copyright 2020, UT-Battelle, LLC

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-----------------------------------------------------------------------------*/

#if defined COMET_USE_CUDA
#  include "cub.cuh"
#elif defined COMET_USE_HIP
#  include "rocprim/rocprim.hpp"
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
  , buf_length_max_(buf_->dim0 * buf_->dim1 * NUM_VALUES_PER_METRIC)
  , num_nonzeros_approx_(0)
  , do_compress_(false)
  , num_runs_(0)
  , state_(State::IDLE)
  , num_entries_(buf_length_()) {
  //, is_open_(false)
  //, read_ptr_(0) {

  if (!can_compress_())
    return;

  COMET_INSIST(env_.metric_format() == METRIC_FORMAT);

  num_nonzeros_buf_.allocate(1, 1, sizeof(MFTTypeIn));

  keys_buf_.allocate(buf_length_max_, 1, sizeof(MFTTypeIn));

  lengths_buf_.allocate(buf_length_max_, 1, sizeof(Lengths_t));

  num_runs_buf_.allocate(1, 1, sizeof(size_t));

  // Anticipate size--usually seems to be = 1.
  reduce_workspace_buf_.allocate(1, 1, sizeof(char));
}

//-----------------------------------------------------------------------------
/// \brief Replace the base underlying MirroredBuf object.

void CompressedBuf::attach(MirroredBuf& buf) {
  COMET_INSIST(State::IDLE == state_);

  buf_ = &buf;
  COMET_INSIST(buf_length_() <= buf_length_max_);

  num_entries_ = buf_length_();
}

//-----------------------------------------------------------------------------
/// \brief Compute number of nonzero MFTTypeIn values in buffer device memory.

void CompressedBuf::compute_num_nonzeros_() {
  COMET_INSIST(can_compress_());

# if defined COMET_USE_ACCEL

    ReductionOp reduction_op;
    size_t temp_storage_bytes = 0;

    // Get amount of temp storage needed.

    const MFTTypeIn initial_value = ReductionOp::INITIAL_VALUE;

# if defined COMET_USE_CUDA

    // see https://nvlabs.github.io/cub/structcub_1_1_device_reduce.html

    cub::DeviceReduce::Reduce(NULL,
      temp_storage_bytes, (MFTTypeIn*)buf_->d, (MFTTypeIn*)num_nonzeros_buf_.d,
      buf_length_(), reduction_op, initial_value, env_.stream_compute());

# elif defined COMET_USE_HIP

    // see https://github.com/ROCmSoftwarePlatform/rocPRIM/blob/develop/rocprim/include/rocprim/device/device_reduce.hpp

    rocprim::reduce(NULL,
      temp_storage_bytes, (MFTTypeIn*)buf_->d, (MFTTypeIn*)num_nonzeros_buf_.d,
      initial_value, buf_length_(), reduction_op, env_.stream_compute());

# endif // COMET_USE_CUDA || COMET_USE_HIP

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

# if defined COMET_USE_CUDA

    cub::DeviceReduce::Reduce((Workspace_t*)reduce_workspace_buf_.d,
      temp_storage_bytes, (MFTTypeIn*)buf_->d, (MFTTypeIn*)num_nonzeros_buf_.d,
      buf_length_(), reduction_op, initial_value, env_.stream_compute());

# elif defined COMET_USE_HIP

    rocprim::reduce((Workspace_t*)reduce_workspace_buf_.d,
      temp_storage_bytes, (MFTTypeIn*)buf_->d, (MFTTypeIn*)num_nonzeros_buf_.d,
      initial_value, buf_length_(), reduction_op, env_.stream_compute());

# endif // COMET_USE_CUDA || COMET_USE_HIP

    // Retrieve count.
    // NOTE: this value may be approximate, since reduction is over floats.

    num_nonzeros_buf_.from_accel(env_.stream_compute());

    num_nonzeros_approx_ = (size_t)
      -(num_nonzeros_buf_.elt_const<MFTTypeIn>(0, 0) - initial_value);

#endif // COMET_USE_ACCEL
}

//-----------------------------------------------------------------------------
/// \brief Compress the underlying buffer if it meets conditions.

void CompressedBuf::compress() {

  if (!can_compress_())
    return;

  compute_num_nonzeros_();

  const size_t num_runs_max = 2 * num_nonzeros_approx_;
  const size_t estimated_storage_compressed = num_runs_max *
    (sizeof(MFTTypeIn) + sizeof(Lengths_t));
  const size_t storage_uncompressed = buf_length_() * sizeof(MFTTypeIn);

  do_compress_ = estimated_storage_compressed <
    compression_factor_required_() * storage_uncompressed &&
    env_.try_compress();

  if (do_compress_) {

    COMET_INSIST(State::IDLE == state_);

#   if defined COMET_USE_ACCEL


      size_t temp_storage_bytes = 0;

      // Get amount of temp storage needed.

#   if defined COMET_USE_CUDA

      // see https://nvlabs.github.io/cub/structcub_1_1_device_run_length_encode.html

      cub::DeviceRunLengthEncode::Encode(NULL,
        temp_storage_bytes, (MFTTypeIn*)buf_->d, (MFTTypeIn*)keys_buf_.d,
        (Lengths_t*)lengths_buf_.d, (size_t*)num_runs_buf_.d, buf_length_(),
        env_.stream_compute());

# elif defined COMET_USE_HIP

      // see https://github.com/ROCmSoftwarePlatform/rocPRIM/blob/develop/rocprim/include/rocprim/device/device_run_length_encode.hpp

      rocprim::run_length_encode(NULL,
        temp_storage_bytes, (MFTTypeIn*)buf_->d, buf_length_(),
        (MFTTypeIn*)keys_buf_.d, (Lengths_t*)lengths_buf_.d,
        (size_t*)num_runs_buf_.d,
        env_.stream_compute());

# endif // COMET_USE_CUDA || COMET_USE_HIP

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

#   if defined COMET_USE_CUDA

      cub::DeviceRunLengthEncode::Encode((Workspace_t*)rle_workspace_buf_.d,
        temp_storage_bytes, (MFTTypeIn*)buf_->d, (MFTTypeIn*)keys_buf_.d,
        (Lengths_t*)lengths_buf_.d, (size_t*)num_runs_buf_.d, buf_length_(),
        env_.stream_compute());

# elif defined COMET_USE_HIP

      rocprim::run_length_encode((Workspace_t*)rle_workspace_buf_.d,
        temp_storage_bytes, (MFTTypeIn*)buf_->d, buf_length_(),
        (MFTTypeIn*)keys_buf_.d, (Lengths_t*)lengths_buf_.d,
        (size_t*)num_runs_buf_.d,
        env_.stream_compute());

# endif // COMET_USE_CUDA || COMET_USE_HIP

      // Retrieve size.

      num_runs_buf_.from_accel(env_.stream_compute());

      num_runs_ = num_runs_buf_.elt_const<size_t>(0, 0);

      // Create right-sized buf aliases to later capture run lengths/keys.

      keys_alias_buf_.allocate(keys_buf_, num_runs_);
      lengths_alias_buf_.allocate(lengths_buf_, num_runs_);

#endif // COMET_USE_ACCEL

    state_ = State::COMPRESSED;

  } // if (do_compress_)
}

//-----------------------------------------------------------------------------
/// \brief CompressedBuf copy back from accel: start.

void CompressedBuf::from_accel_start() {

  if (do_compress_) {
    COMET_INSIST(State::COMPRESSED == state_);
    keys_alias_buf_.from_accel_start();
    lengths_alias_buf_.from_accel_start();
    state_ = State::TRANSFER_STARTED;
#ifdef COMET_ASSERTIONS_ON
    //buf_->from_accel_start();
#endif
  } else {
    buf_->from_accel_start();
  }
}

//-----------------------------------------------------------------------------
/// \brief CompressedBuf copy back from accel: wait.

void CompressedBuf::from_accel_wait() {

  if (do_compress_) {
    COMET_INSIST(State::TRANSFER_STARTED == state_);
    keys_alias_buf_.from_accel_wait();
    lengths_alias_buf_.from_accel_wait();

#if 0
    for (size_t i=0; i<num_runs_; ++i) {
      printf("cbuf %i %i %f\n", (int)i,
        (int)lengths_alias_buf_.elt_const<size_t>(i, 0),
        (double)keys_alias_buf_.elt_const<MFTTypeIn>(i, 0) ); 
    }
#endif

    num_entries_ = 0;
    for (size_t i = 0; i < num_runs_; ++i) {
      if (keys_alias_buf_.elt_const<MFTTypeIn>(i, 0))
        num_entries_ += lengths_alias_buf_.elt_const<Lengths_t>(i, 0);
    } // i
//printf("%i\n", (int)num_entries_);

#if 0
double mysum = 0;
    for (size_t i = 0; i < num_runs_; ++i) {
mysum += 
keys_alias_buf_.elt_const<MFTTypeIn>(i, 0)
* lengths_alias_buf_.elt_const<Lengths_t>(i, 0);
    } // i
printf("%.20e\n", mysum);
#endif

    elt_read_start();

    state_ = State::IDLE;
#ifdef COMET_ASSERTIONS_ON
    //buf_->from_accel_wait();
#endif

#if 0
    for (size_t i=0; i<64; ++i) {
      if ((double)(((MFTTypeIn*)(buf_->h))[i]))
        printf("buf %i %f\n", (int)i,
          (double)(((MFTTypeIn*)(buf_->h))[i]));
    }
#endif

  } else {

    buf_->from_accel_wait();

  }
}

//-----------------------------------------------------------------------------
/// \brief CompressedBuf lock buffer for exclusive use.

void CompressedBuf::lock_h() {

  if (do_compress_) {
    keys_alias_buf_.lock_h();
    lengths_alias_buf_.lock_h();
  } else {
    buf_->lock_h();
  }
}

//-----------------------------------------------------------------------------
/// \brief CompressedBuf unlock buffer.

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
