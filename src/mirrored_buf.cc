//-----------------------------------------------------------------------------
/*!
 * \file   mirrored_buf.cc
 * \author Wayne Joubert
 * \date   Thu Aug  3 15:04:05 EDT 2017
 * \brief  Struct/code to manage dual CPU/GPU reflected arrays.
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

#include "magma_wrapper.hh"
#include "mirrored_buf.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief MirroredBuf contructor: create empty with no allocations.

MirroredBuf::MirroredBuf(CEnv& env)
  : env_(env)
  , h(NULL)
  , d(NULL)
  , active(NULL)
  , dim0(0)
  , dim1(0)
  , num_elts_(0)
  , elt_size_(0)
  , size_allocated_(0)
  , is_alias(false)
  , is_allocated(false)
  , is_locked_h_(false)
  , is_locked_d_(false)
  , use_linalg_(false) {
  }

//-----------------------------------------------------------------------------
/// \brief MirroredBuf: constructor: specify elt size.

MirroredBuf::MirroredBuf(size_t dim0_, size_t dim1_, int elt_size, CEnv& env)
  : env_(env)
  , h(NULL)
  , d(NULL)
  , active(NULL)
  , dim0(dim0_)
  , dim1(dim1_)
  , num_elts_(dim0_ * dim1_)
  , elt_size_(elt_size)
  , size_allocated_(num_elts_ * elt_size_)
  , is_alias(false)
  , is_allocated(false)
  , is_locked_h_(false)
  , is_locked_d_(false)
  , use_linalg_(false) {
  allocate(dim0_, dim1_, elt_size);
}

//-----------------------------------------------------------------------------
/// \brief MirroredBuf: constructor: infer elt size from metric_type.

MirroredBuf::MirroredBuf(size_t dim0_, size_t dim1_, CEnv& env)
  : env_(env)
  , h(NULL)
  , d(NULL)
  , active(NULL)
  , dim0(dim0_)
  , dim1(dim1_)
  , num_elts_(dim0_ * dim1_)
  , elt_size_(env_.matrix_buf_elt_size())
  , size_allocated_(num_elts_ * elt_size_)
  , is_alias(false)
  , is_allocated(false)
  , is_locked_h_(false)
  , is_locked_d_(false)
  , use_linalg_(false) {
  allocate(dim0_, dim1_);
}

//-----------------------------------------------------------------------------
/// \brief MirroredBuf: constructor: alias, with possible redimension.

MirroredBuf::MirroredBuf(MirroredBuf& buf, size_t dim0_, CEnv& env)
  : env_(env)
  , h(buf.h)
  , d(buf.d)
  , active(buf.active)
  , dim0(dim0_)
  , dim1(buf.dim1)
  , num_elts_(dim0 * dim1)
  , elt_size_(buf.elt_size_)
  , size_allocated_(buf.size_allocated_)
  , is_alias(true)
  , is_allocated(true)
  , is_locked_h_(false)
  , is_locked_d_(false)
  //, use_linalg_(BuildHas::MAGMA) {
  , use_linalg_(false) {
  COMET_INSIST(dim0_ <= buf.dim0);
  COMET_INSIST(buf.is_allocated);
}

//-----------------------------------------------------------------------------

MirroredBuf::~MirroredBuf() {
  deallocate();
}

//-----------------------------------------------------------------------------
/// \brief MirroredBuf: allocate memory: specify elt size.

void MirroredBuf::allocate(size_t dim0_, size_t dim1_, int elt_size) {
  COMET_INSIST(is_alias || !is_allocated);

  use_linalg_ = false;

  dim0 = dim0_;
  dim1 = dim1_;
  num_elts_ = dim0 * dim1;
  elt_size_ = elt_size;
  size_allocated_ = num_elts_ * elt_size_;
  size_allocated_ = size_allocated_ ? size_allocated_ : 1;
# if defined COMET_USE_CUDA
//#if defined(COMET_PLATFORM_CORI_GPU) || defined(COMET_PLATFORM_JUWELS_BOOSTER)
    if (System::compute_capability() >= 800 ||
        !env_.is_compute_method_gpu()) {
      // WORKAROUND
      h = malloc(size_allocated_);
    } else {
      cudaMallocHost((void**)&h, size_allocated_);
    }
    COMET_INSIST(System::accel_last_call_succeeded());
    if (env_.is_compute_method_gpu()) {
//printf("1 %zu %i\n", size_allocated_, System::proc_num());
      cudaMalloc((void**)&d, size_allocated_);
      COMET_INSIST(System::accel_last_call_succeeded());
    }
# elif defined COMET_USE_HIP
    hipHostMalloc((void**)&h, size_allocated_);
    COMET_INSIST(System::accel_last_call_succeeded());
    if (env_.is_compute_method_gpu()) {
      hipMalloc((void**)&d, size_allocated_);
      COMET_INSIST(System::accel_last_call_succeeded());
    }
# else
    h = malloc(size_allocated_);
    COMET_INSIST(!env_.is_compute_method_gpu() &&
      "GPU not supported for this build.");
# endif

  env_.cpu_mem_local_inc(size_allocated_);
  if (env_.is_compute_method_gpu())
    env_.gpu_mem_local_inc(size_allocated_);

  COMET_INSIST(h &&
    "Invalid host pointer created, possibly due to insufficient memory.");
  COMET_INSIST((d || !env_.is_compute_method_gpu()) &&
    "Invalid device pointer created, possibly due to insufficient memory.");

  active = env_.is_compute_method_gpu() ? d : h;
  is_alias = false;
  is_allocated = true;
}

//-----------------------------------------------------------------------------
/// \brief MirroredBuf: allocate memory: infer elt size from metric_type.

void MirroredBuf::allocate(size_t dim0_, size_t dim1_) {
  COMET_INSIST(is_alias || !is_allocated);

  dim0 = dim0_;
  dim1 = dim1_;
  num_elts_ = dim0 * dim1;
  elt_size_ = env_.matrix_buf_elt_size();
  size_allocated_ = num_elts_ * elt_size_;
  size_allocated_ = size_allocated_ ? size_allocated_ : 1;

  use_linalg_ = BuildHas::MAGMA;

  if (use_linalg_) {

    MagmaWrapper::malloc(this, dim0_, dim1_, env_);

    env_.cpu_mem_local_inc(size_allocated_);
    if (env_.is_compute_method_gpu())
      env_.gpu_mem_local_inc(size_allocated_);

    active = env_.is_compute_method_gpu() ? d : h;
    is_alias = false;
    is_allocated = true;

  } else {

    MirroredBuf::allocate(dim0_, dim1_, elt_size_);

  } // if (use_linalg_)
}

//-----------------------------------------------------------------------------
/// \brief MirroredBuf: allocate memory: alias, with possible redimension.

void MirroredBuf::allocate(MirroredBuf& buf, size_t dim0_) {
  COMET_INSIST(is_alias || !is_allocated);
  COMET_INSIST(dim0_ <= buf.dim0);
  COMET_INSIST(!buf.is_alias);
  COMET_INSIST(buf.is_allocated);

  use_linalg_ = buf.use_linalg_;
  h = buf.h;
  d = buf.d;
  active = buf.active;
  dim0 = dim0_;
  dim1 = buf.dim1;
  num_elts_ = dim0 * dim1;
  elt_size_ = buf.elt_size_;
  size_allocated_ = buf.size_allocated_;
  is_alias = true;
  is_allocated = true;
}

//-----------------------------------------------------------------------------

void MirroredBuf::set_zero_h() {
  COMET_INSIST(is_allocated);

//# pragma omp parallel for schedule(dynamic,1000)
  for (size_t i=0; i<size_allocated_; ++i)
    ((char*)h)[i] = 0;
}

//-----------------------------------------------------------------------------

void MirroredBuf::deallocate() {
  COMET_INSIST(!is_locked_h_ && !is_locked_d_);

  if (is_allocated && !is_alias) {

    if (use_linalg_) {

      MagmaWrapper::free(this, env_);

      env_.cpu_mem_local_dec(size_allocated_);
      if (env_.is_compute_method_gpu())
        env_.gpu_mem_local_dec(size_allocated_);

    } else {

#     if defined COMET_USE_CUDA
//#if defined(COMET_PLATFORM_CORI_GPU) || defined(COMET_PLATFORM_JUWELS_BOOSTER)
        if (System::compute_capability() >= 800 ||
            !env_.is_compute_method_gpu()) {
          // WORKAROUND
          free(h);
        } else {
          cudaFreeHost(h);
        }
        COMET_INSIST(System::accel_last_call_succeeded());
        if (env_.is_compute_method_gpu()) {
          cudaFree(d);
          COMET_INSIST(System::accel_last_call_succeeded());
        }
#     elif defined COMET_USE_HIP
        hipHostFree(h);
        COMET_INSIST(System::accel_last_call_succeeded());
        if (env_.is_compute_method_gpu()) {
          hipFree(d);
          COMET_INSIST(System::accel_last_call_succeeded());
        }
#     else
        free(h);
        COMET_INSIST(!env_.is_compute_method_gpu() &&
          "GPU not supported for this build.");
#     endif

      env_.cpu_mem_local_dec(size_allocated_);
      if (env_.is_compute_method_gpu())
        env_.gpu_mem_local_dec(size_allocated_);

    } // if (use_linalg_)

    h = NULL;
    d = NULL;

    active = NULL;
    is_allocated = false;

  } // if

  use_linalg_ = false;
}

//-----------------------------------------------------------------------------

void MirroredBuf::to_accel_start() {

  if (!env_.is_compute_method_gpu())
    return;

  lock();

  if (use_linalg_) {
    MagmaWrapper::set_matrix_start(this, env_);
  } else {
#   if defined COMET_USE_CUDA
      cudaMemcpyAsync(d, h, size(), cudaMemcpyHostToDevice, env_.stream_togpu());
#   elif defined COMET_USE_HIP
      hipMemcpyAsync(d, h, size(), hipMemcpyHostToDevice, env_.stream_togpu());
#   endif
    COMET_INSIST(System::accel_last_call_succeeded());
  } // if (use_linalg_)
}

//-----------------------------------------------------------------------------

void MirroredBuf::to_accel_wait() {

  if (!env_.is_compute_method_gpu())
    return;

  if (use_linalg_)
    MagmaWrapper::set_matrix_wait(env_);
  else
    env_.stream_synchronize(env_.stream_togpu());

  unlock();
}

//-----------------------------------------------------------------------------

void MirroredBuf::to_accel() {
  to_accel_start();
  to_accel_wait();
}

//-----------------------------------------------------------------------------

void MirroredBuf::from_accel_start(AccelStream_t stream) {

  if (!env_.is_compute_method_gpu())
    return;

  lock();

  if (use_linalg_) {
    COMET_INSIST(env_.stream_fromgpu() == stream);
    MagmaWrapper::get_matrix_start(this, env_);
  } else {
#   if defined COMET_USE_CUDA
      cudaMemcpyAsync(h, d, size(), cudaMemcpyDeviceToHost, stream);
#   elif defined COMET_USE_HIP
      hipMemcpyAsync(h, d, size(), hipMemcpyDeviceToHost, stream);
#   endif
    COMET_INSIST(System::accel_last_call_succeeded());
  } // if (use_linalg_)
}

//-----------------------------------------------------------------------------

void MirroredBuf::from_accel_wait(AccelStream_t stream) {

  if (!env_.is_compute_method_gpu())
    return;

  if (use_linalg_) {
    COMET_INSIST(env_.stream_fromgpu() == stream);
    MagmaWrapper::get_matrix_wait(env_);
  } else
    env_.stream_synchronize(stream);

  unlock();
}

//-----------------------------------------------------------------------------

void MirroredBuf::from_accel(AccelStream_t stream) {
  from_accel_start(stream);
  from_accel_wait(stream);
}

//-----------------------------------------------------------------------------

void MirroredBuf::from_accel_start() {

  from_accel_start(env_.stream_fromgpu());
}

//-----------------------------------------------------------------------------

void MirroredBuf::from_accel_wait() {

  from_accel_wait(env_.stream_fromgpu());
}

//-----------------------------------------------------------------------------

void MirroredBuf::from_accel() {
  from_accel_start();
  from_accel_wait();
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
