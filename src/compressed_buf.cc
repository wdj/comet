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
  , is_open_(false)
  , read_ptr_(0) {

  if (!try_compress())
    return;

  COMET_INSIST(env_.metric_format() == METRIC_FORMAT);

  num_nonzeros_buf_.allocate(1, 1, sizeof(MFTTypeIn));



  keys_buf_.allocate(length_max_, 1, sizeof(MFTTypeIn));

  lengths_buf_.allocate(length_max_, 1, sizeof(MFTTypeIn));



}

//-----------------------------------------------------------------------------

void CompressedBuf::compute_nonzeros_() {

  if (!try_compress())
    return;

# if defined COMET_USE_CUDA

    ReductionOp reduction_op;
    const size_t length = buf_.dim0 * buf_.dim1 * NUM_VALUES_PER_METRIC;
    COMET_INSIST(length <= length_max_);
    size_t temp_storage_bytes = 0;

    cub::DeviceReduce::Reduce((MFTTypeIn*)reduce_workspace_buf_.d,
      temp_storage_bytes, (MFTTypeIn*)buf_.d, (MFTTypeIn*)num_nonzeros_buf_.d,
      length, reduction_op, (MFTTypeIn)0);





# elif defined COMET_USE_HIP

  // TODO

# endif // COMET_USE_CUDA || COMET_USE_HIP





#if 0
see https://nvlabs.github.io/cub/structcub_1_1_device_reduce.html#aa4adabeb841b852a7a5ecf4f99a2daeb

- if defined COMET_USE_CUDA ...
- call cub to get temp array size required
- if larger than workspace then
  - deallocate workspace if allocated
  - allocate to new size - perhaps mke larger than requirement by small fraction





#endif

}



//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
