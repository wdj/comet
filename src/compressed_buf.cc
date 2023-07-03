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
#  define CUB_NS_QUALIFIER cub
#  include "cub.cuh"
#elif defined COMET_USE_HIP && defined COMET_USE_ROCPRIM
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

  if (!is_compress_enabled_())
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

#if defined COMET_USE_HIP && ! defined COMET_USE_ROCPRIM

//-----------------------------------------------------------------------------

template<class MFTTypeIn>
__global__ static void reduce_kernel_(MFTTypeIn* buf_d, size_t buf_length,
  MFTTypeIn* num_nonzeros_buf_d) {

  const int threadindex = threadIdx_x_();
  const int blockindex = blockIdx_x_();
  const int blocksize = blockDim_x_();
  const int gridsize = gridDim_x_();

  const size_t thread = threadindex + blocksize * static_cast<size_t>(
                        blockindex);

  typedef int Dummy_t;
  extern __shared__ Dummy_t sdata_[];
  auto sdata = reinterpret_cast<MFTTypeIn*>(&(sdata_[0]));
  sdata[threadindex] = 0;

  size_t index = thread;

 // Serially reduce across (threadblock*grid)-sized chunks.

  while (index < buf_length) {
    sdata[threadindex] += static_cast<MFTTypeIn>(buf_d[index] != static_cast<MFTTypeIn>(0));
    index += blocksize * static_cast<size_t>(gridsize);
  } // while

  // Reduce within threadblock.

  for (unsigned int s = blocksize/2; s > 0; s >>= 1) {
    __syncthreads();
    if (threadindex < s) { // parallel sweep reduction
      sdata[threadindex] += sdata[threadindex + s];
    }
  }

  // First thread of each threadblock adds in its contribution.

  if (0 == threadindex)
    atomicAdd(&(num_nonzeros_buf_d[0]), sdata[0]);
}

//-----------------------------------------------------------------------------

template<class MFTTypeIn>
static void reduce_(MFTTypeIn* buf_d, size_t buf_length,
  MFTTypeIn* num_nonzeros_buf_d, AccelStream_t accel_stream) {

  const auto num_thread = static_cast<NML_t>(buf_length);

  hipMemsetAsync(num_nonzeros_buf_d, 0, 1*sizeof(MFTTypeIn), accel_stream);

  // ISSUE: this may need to be a power of 2.
  const int threadblocksize = 256;
  COMET_INSIST((threadblocksize <= 256 || ! BuildHas::HIP) &&
               "Current HIP limitation.");
  const int num_threadblocks_0 = utils::ceil(num_thread,
                                   static_cast<NML_t>(threadblocksize));

  COMET_LAUNCH_KERNEL((reduce_kernel_<MFTTypeIn>),
    dim3(num_threadblocks_0, 1, 1),
    dim3(threadblocksize, 1, 1),
    threadblocksize * sizeof(MFTTypeIn),
    accel_stream,
    buf_d,
    buf_length,
    num_nonzeros_buf_d);

  System::accel_last_call_succeeded();
} // reduce_

//-----------------------------------------------------------------------------

#endif

//-----------------------------------------------------------------------------
/// \brief Compute number of nonzero MFTTypeIn values in buffer device memory.

void CompressedBuf::compute_num_nonzeros_() {
  COMET_INSIST(is_compress_enabled_());

# if defined COMET_USE_ACCEL

    // Initializations.

    ReductionOp reduction_op;
    size_t temp_storage_bytes = 0;

    // First get amount of temp storage needed for reduction.

    const MFTTypeIn initial_value = ReductionOp::INITIAL_VALUE;

#   if defined COMET_USE_CUDA

      // see https://nvlabs.github.io/cub/structcub_1_1_device_reduce.html

      cub::DeviceReduce::Reduce(NULL,
        temp_storage_bytes, (MFTTypeIn*)buf_->d, (MFTTypeIn*)num_nonzeros_buf_.d,
        buf_length_(), reduction_op, initial_value, env_.stream_fromgpu());

#   elif defined COMET_USE_HIP && defined COMET_USE_ROCPRIM

      // see https://github.com/ROCmSoftwarePlatform/rocPRIM/blob/develop/rocprim/include/rocprim/device/device_reduce.hpp

      rocprim::reduce(NULL,
        temp_storage_bytes, (MFTTypeIn*)buf_->d, (MFTTypeIn*)num_nonzeros_buf_.d,
        initial_value, buf_length_(), reduction_op, env_.stream_fromgpu());

#   elif defined COMET_USE_HIP // && ! defined COMET_USE_ROCPRIM

      // no temp storage needed //

#   endif // COMET_USE_CUDA || COMET_USE_HIP

    num_nonzeros_buf_.from_accel(env_.stream_fromgpu());

    // Re/allocate temp storage as needed.

    if (temp_storage_bytes > reduce_workspace_buf_.size()) {
      reduce_workspace_buf_.deallocate();
      // Possibly allocate slightly extra to allow some future growth
      // without need to resize.
      const double fuzz = 0;
      const size_t alloc_size = temp_storage_bytes * (1+fuzz);
      COMET_INSIST(alloc_size >= temp_storage_bytes);
      reduce_workspace_buf_.allocate(alloc_size, 1, sizeof(char));
    }

    // Now perform reduction.

#   if defined COMET_USE_CUDA

      cub::DeviceReduce::Reduce((Workspace_t*)reduce_workspace_buf_.d,
        temp_storage_bytes, (MFTTypeIn*)buf_->d, (MFTTypeIn*)num_nonzeros_buf_.d,
        buf_length_(), reduction_op, initial_value, env_.stream_fromgpu());

#   elif defined COMET_USE_HIP && defined COMET_USE_ROCPRIM

      rocprim::reduce((Workspace_t*)reduce_workspace_buf_.d,
        temp_storage_bytes, (MFTTypeIn*)buf_->d, (MFTTypeIn*)num_nonzeros_buf_.d,
        initial_value, buf_length_(), reduction_op, env_.stream_fromgpu());

#   elif defined COMET_USE_HIP // && ! defined COMET_USE_ROCPRIM

      reduce_(reinterpret_cast<MFTTypeIn*>(buf_->d), buf_length_(),
              reinterpret_cast<MFTTypeIn*>(num_nonzeros_buf_.d),
              env_.stream_fromgpu());

#   endif // COMET_USE_CUDA || COMET_USE_HIP

    // Retrieve nonzeros count.
    // NOTE: this value may be approximate, since reduction is over floats.

    num_nonzeros_buf_.from_accel(env_.stream_fromgpu());

    num_nonzeros_approx_ = (size_t)
      -(num_nonzeros_buf_.elt_const<MFTTypeIn>(0, 0) - initial_value);

    const size_t roundoff_limit = 4e6;
    COMET_INSIST((size_t)(float)roundoff_limit == roundoff_limit);
    COMET_INSIST(num_nonzeros_approx_ <= buf_length_() ||
                 buf_length_() > roundoff_limit);

# endif // COMET_USE_ACCEL
}

//-----------------------------------------------------------------------------

#if defined COMET_USE_HIP && ! defined COMET_USE_ROCPRIM

//-----------------------------------------------------------------------------

  // TODO - need prefix scan code

//-----------------------------------------------------------------------------

#endif

//-----------------------------------------------------------------------------
/// \brief Compress the underlying buffer if needed.

void CompressedBuf::compress() {

  if (!is_compress_enabled_())
    return;

 // Determine whether it is worth the effort to compress.

  compute_num_nonzeros_();

  const size_t num_runs_max = 2 * num_nonzeros_approx_;
  const size_t estimated_storage_compressed = num_runs_max *
    (sizeof(MFTTypeIn) + sizeof(Lengths_t));
  const size_t storage_uncompressed = buf_length_() * sizeof(MFTTypeIn);

  do_compress_ = estimated_storage_compressed <
    compression_factor_required_() * storage_uncompressed;

  if (!do_compress_)
    return;

  COMET_INSIST(State::IDLE == state_);

# if defined COMET_USE_ACCEL

    size_t temp_storage_bytes = 0;

    // First get amount of temp storage needed for rle.

#   if defined COMET_USE_CUDA

      // see https://nvlabs.github.io/cub/structcub_1_1_device_run_length_encode.html

      cub::DeviceRunLengthEncode::Encode(NULL,
        temp_storage_bytes, (MFTTypeIn*)buf_->d, (MFTTypeIn*)keys_buf_.d,
        (Lengths_t*)lengths_buf_.d, (size_t*)num_runs_buf_.d, buf_length_(),
        env_.stream_fromgpu());

#   elif defined COMET_USE_HIP && defined COMET_USE_ROCPRIM

      // see https://github.com/ROCmSoftwarePlatform/rocPRIM/blob/develop/rocprim/include/rocprim/device/device_run_length_encode.hpp

      rocprim::run_length_encode(NULL,
        temp_storage_bytes, (MFTTypeIn*)buf_->d, buf_length_(),
        (MFTTypeIn*)keys_buf_.d, (Lengths_t*)lengths_buf_.d,
        (size_t*)num_runs_buf_.d,
        env_.stream_fromgpu());

#   elif defined COMET_USE_HIP // && ! defined COMET_USE_ROCPRIM

  // TODO

  // https://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_website/samples.html
  // https://github.com/ROCmSoftwarePlatform/rocPRIM/blob/develop/rocprim/include/rocprim/device/device_run_length_encode.hpp
  // https://erkaman.github.io/posts/cuda_rle.html
  // https://www.eecs.umich.edu/courses/eecs570/hw/parprefix.pdf
  // https://github.com/NVIDIA/cuda-samples/tree/master/Samples/2_Concepts_and_Techniques/scan

#   endif // COMET_USE_CUDA || COMET_USE_HIP

    num_runs_buf_.from_accel(env_.stream_fromgpu());

    // Re/allocate temp storage as needed.

    if (temp_storage_bytes > rle_workspace_buf_.size()) {
      rle_workspace_buf_.deallocate();
      // Possibly allocate slightly extra to allow some future growth
      // without need to resize.
      const double fuzz = 0;
      const size_t alloc_size = temp_storage_bytes * (1+fuzz);
      COMET_INSIST(alloc_size >= temp_storage_bytes);
      rle_workspace_buf_.allocate(alloc_size, 1, sizeof(char));
    }

    // Now perform RLE.

#   if defined COMET_USE_CUDA

      cub::DeviceRunLengthEncode::Encode((Workspace_t*)rle_workspace_buf_.d,
        temp_storage_bytes, (MFTTypeIn*)buf_->d, (MFTTypeIn*)keys_buf_.d,
        (Lengths_t*)lengths_buf_.d, (size_t*)num_runs_buf_.d, buf_length_(),
        env_.stream_fromgpu());

#   elif defined COMET_USE_HIP && defined COMET_USE_ROCPRIM

      rocprim::run_length_encode((Workspace_t*)rle_workspace_buf_.d,
        temp_storage_bytes, (MFTTypeIn*)buf_->d, buf_length_(),
        (MFTTypeIn*)keys_buf_.d, (Lengths_t*)lengths_buf_.d,
        (size_t*)num_runs_buf_.d,
        env_.stream_fromgpu());

#   elif defined COMET_USE_HIP // && ! defined COMET_USE_ROCPRIM

  // TODO

#   endif // COMET_USE_CUDA || COMET_USE_HIP

    // Retrieve number of runs.

    num_runs_buf_.from_accel(env_.stream_fromgpu());

    num_runs_ = num_runs_buf_.elt_const<size_t>(0, 0);

    // Create smaller-sized buf aliases, for later retrieval
    // of run lengths and keys.

    keys_alias_buf_.allocate(keys_buf_, num_runs_);
    lengths_alias_buf_.allocate(lengths_buf_, num_runs_);

# endif // COMET_USE_ACCEL

  state_ = State::COMPRESSED;
}

//-----------------------------------------------------------------------------
/// \brief CompressedBuf copy back from accel: start.

void CompressedBuf::from_accel_start() {

  if (do_compress_) {

    COMET_INSIST(State::COMPRESSED == state_);

    // Begin retrieval of rle data.

    keys_alias_buf_.from_accel_start();
    lengths_alias_buf_.from_accel_start();

    state_ = State::TRANSFER_STARTED;

#   ifdef _COMET_COMPRESSED_BUF_CHECK_RESULT
        buf_->from_accel_start();
#   endif

  } else { // ! do_compress_

    buf_->from_accel_start();

  } // if (do_compress_)
}

//-----------------------------------------------------------------------------
/// \brief CompressedBuf copy back from accel: wait.

void CompressedBuf::from_accel_wait() {

  if (do_compress_) {

    COMET_INSIST(State::TRANSFER_STARTED == state_);

    // Finish retrieval of rle data.

    keys_alias_buf_.from_accel_wait();
    lengths_alias_buf_.from_accel_wait();

    // Count total number of nonzero values present after uncompression.

    num_entries_ = 0;
    for (size_t i = 0; i < num_runs_; ++i) {
      if (keys_alias_buf_.elt_const<MFTTypeIn>(i, 0))
        num_entries_ += lengths_alias_buf_.elt_const<Lengths_t>(i, 0);
    } // i

    // Prepare for calls to read rle data.

    elt_read_start();

    state_ = State::IDLE;

#   ifdef _COMET_COMPRESSED_BUF_CHECK_RESULT
      buf_->from_accel_wait();
#   endif

  } else { // ! do_compress_

    buf_->from_accel_wait();


  } // if (do_compress_)
}

//-----------------------------------------------------------------------------
/// \brief Lock CompressedBuf buffer for exclusive use.

void CompressedBuf::lock_h() {

  if (do_compress_) {
    keys_alias_buf_.lock_h();
    lengths_alias_buf_.lock_h();
  } else { // ! do_compress_
    buf_->lock_h();
  } // if (do_compress_)
}

//-----------------------------------------------------------------------------
/// \brief Unlock CompressedBuf buffer.

void CompressedBuf::unlock_h() {

  if (do_compress_) {
    keys_alias_buf_.unlock_h();
    lengths_alias_buf_.unlock_h();
  } else { // ! do_compress_
    buf_->unlock_h();
  } // if (do_compress_)
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
