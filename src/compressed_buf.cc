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

CompressedBuf::CompressedBuf(MirroredBuf& buf, CEnv& env)
  : env_(env)
  , buf_(buf)
  , num_nonzeros_buf_(env)
  , num_runs_buf_(env)
  , keys_buf_(env)
  , lengths_buf_(env)
  , keys_alias_buf_(env)
  , lengths_alias_buf_(env)
  , reduce_workspace_buf_(env)
  , rle_workspace_buf_(env)
  , length_max_(buf_.dim0 * buf_.dim1 * NUM_VALUES_PER_METRIC)
  , num_nonzeros_approx_(0)
  , do_compress_(false)
  , num_runs_(0)
  , is_open_(false)
  , read_ptr_(0) {

  if (!try_compress())
    return;

  COMET_INSIST(env_.metric_format() == METRIC_FORMAT);

  num_nonzeros_buf_.allocate(1, 1, sizeof(MFTTypeIn));

  keys_buf_.allocate(length_max_, 1, sizeof(MFTTypeIn));

  lengths_buf_.allocate(length_max_, 1, sizeof(size_t));

  num_runs_buf_.allocate(1, 1, sizeof(size_t));

}

//-----------------------------------------------------------------------------

void CompressedBuf::compute_num_nonzeros_() {

# if defined COMET_USE_CUDA

    // see https://nvlabs.github.io/cub/structcub_1_1_device_reduce.html

    ReductionOp reduction_op;
    size_t temp_storage_bytes = 0;

    cub::DeviceReduce::Reduce((Workspace_t*)reduce_workspace_buf_.d,
      temp_storage_bytes, (MFTTypeIn*)buf_.d, (MFTTypeIn*)num_nonzeros_buf_.d,
      length_(), reduction_op, (MFTTypeIn)0, env_.stream_compute());

    if (temp_storage_bytes > reduce_workspace_buf_.size()) {
      reduce_workspace_buf_.deallocate();
      reduce_workspace_buf_.allocate(temp_storage_bytes, 1, sizeof(char));
    }

    cub::DeviceReduce::Reduce((Workspace_t*)reduce_workspace_buf_.d,
      temp_storage_bytes, (MFTTypeIn*)buf_.d, (MFTTypeIn*)num_nonzeros_buf_.d,
      length_(), reduction_op, (MFTTypeIn)0, env_.stream_compute());

    num_nonzeros_buf_.from_accel(env_.stream_compute());

    // NOTE: this may be approximate, since reduction is over floats.
    num_nonzeros_approx_ = (size_t)
      -num_nonzeros_buf_.elt_const<MFTTypeIn>(0, 0);

# elif defined COMET_USE_HIP

  // TODO

# endif // COMET_USE_CUDA || COMET_USE_HIP
}

//-----------------------------------------------------------------------------

void CompressedBuf::compress() {

  if (!try_compress())
    return;

  compute_num_nonzeros_();

  do_compress_ = num_nonzeros_approx_ <= compress_multiplier_() * length_();

  if (do_compress_) {

#   if defined COMET_USE_CUDA
      // see https://nvlabs.github.io/cub/structcub_1_1_device_run_length_encode.html
      size_t temp_storage_bytes = 0;

      cub::DeviceRunLengthEncode::Encode((Workspace_t*)rle_workspace_buf_.d,
        temp_storage_bytes, (MFTTypeIn*)buf_.d, (MFTTypeIn*)keys_buf_.d,
        (size_t*)lengths_buf_.d, (size_t*)num_runs_buf_.d, length_(),
        env_.stream_compute());

      if (temp_storage_bytes > rle_workspace_buf_.size()) {
        rle_workspace_buf_.deallocate();
        rle_workspace_buf_.allocate(temp_storage_bytes, 1, sizeof(char));
      }

      cub::DeviceRunLengthEncode::Encode((Workspace_t*)rle_workspace_buf_.d,
        temp_storage_bytes, (MFTTypeIn*)buf_.d, (MFTTypeIn*)keys_buf_.d,
        (size_t*)lengths_buf_.d, (size_t*)num_runs_buf_.d, length_(),
        env_.stream_compute());

      num_runs_buf_.from_accel(env_.stream_compute());

      num_runs_ = num_runs_buf_.elt_const<size_t>(0, 0);

      keys_alias_buf_.allocate(keys_buf_, num_runs_);
      lengths_alias_buf_.allocate(lengths_buf_, num_runs_);

      keys_alias_buf_.from_accel(env_.stream_compute());
      lengths_alias_buf_.from_accel(env_.stream_compute());

    } // if (do_compress)

#   elif defined COMET_USE_HIP

  // TODO

#   endif // COMET_USE_CUDA || COMET_USE_HIP
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
