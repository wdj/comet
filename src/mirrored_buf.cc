//-----------------------------------------------------------------------------
/*!
 * \file   mirrored_buf.cc
 * \author Wayne Joubert
 * \date   Thu Aug  3 15:04:05 EDT 2017
 * \brief  Struct/code to manage dual CPU/GPU reflected arrays.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "linalg.hh"
#include "mirrored_buf.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

MirroredBuf::MirroredBuf(CEnv& env)
  : h(NULL)
  , d(NULL)
  , active(NULL)
  , dim0(0)
  , dim1(0)
  , size(0)
  , is_alias(false)
  , is_allocated(false)
  , env_(env)
  , is_locked_h_(false)
  , is_locked_d_(false)
  , use_linalg_(BuildHas::MAGMA) {
  }

//-----------------------------------------------------------------------------

MirroredBuf::MirroredBuf(size_t dim0_, size_t dim1_, int elt_size, CEnv& env)
  : h(NULL)
  , d(NULL)
  , active(NULL)
  , dim0(dim0_)
  , dim1(dim1_)
  , size(dim0_ * dim1_)
  , is_alias(false)
  , is_allocated(false)
  , env_(env)
  , is_locked_h_(false)
  , is_locked_d_(false)
  , use_linalg_(false) {
  allocate(dim0_, dim1_, elt_size);
}

//-----------------------------------------------------------------------------

MirroredBuf::MirroredBuf(size_t dim0_, size_t dim1_, CEnv& env)
  : h(NULL)
  , d(NULL)
  , active(NULL)
  , dim0(dim0_)
  , dim1(dim1_)
  , size(dim0_ * dim1_)
  , is_alias(false)
  , is_allocated(false)
  , env_(env)
  , is_locked_h_(false)
  , is_locked_d_(false)
  , use_linalg_(BuildHas::MAGMA) {
  allocate(dim0_, dim1_);
}

//-----------------------------------------------------------------------------

MirroredBuf::MirroredBuf(MirroredBuf& buf, size_t dim0_, CEnv& env)
  : h(buf.h)
  , d(buf.d)
  , active(buf.active)
  , dim0(dim0_)
  , dim1(buf.dim1)
  , size(buf.size)
  , is_alias(true)
  , is_allocated(true)
  , env_(env)
  , is_locked_h_(false)
  , is_locked_d_(false)
  , use_linalg_(BuildHas::MAGMA) {
  COMET_INSIST(dim0_ <= buf.dim0);
  COMET_INSIST(buf.is_allocated);
}

//-----------------------------------------------------------------------------

MirroredBuf::~MirroredBuf() {
  deallocate();
}

//-----------------------------------------------------------------------------

void MirroredBuf::allocate(size_t dim0_, size_t dim1_, int elt_size) {
  COMET_INSIST(is_alias || !is_allocated);

  dim0 = dim0_;
  dim1 = dim1_;
  size = dim0 * dim1 * elt_size;
# if defined COMET_USE_CUDA
    cudaMallocHost((void**)&h, size);
    if (env_.is_compute_method_gpu())
      cudaMalloc((void**)&d, size);
# elif defined COMET_USE_HIP
    hipHostMalloc((void**)&h, size);
    if (env_.is_compute_method_gpu())
      hipMalloc((void**)&d, size);
# else
    h = malloc(size);
    COMET_INSIST(!env_.is_compute_method_gpu() &&
      "GPU not supported for this build.");
# endif

  env_.cpu_mem_local_inc(size);
  if (env_.is_compute_method_gpu())
    env_.gpu_mem_local_inc(size);

  COMET_INSIST(h &&
    "Invalid host pointer created, possibly due to insufficient memory.");
  COMET_INSIST((d || !env_.is_compute_method_gpu()) &&
    "Invalid device pointer created, possibly due to insufficient memory.");

  active = env_.is_compute_method_gpu() ? d : h;
  is_alias = false;
  is_allocated = true;
}

//-----------------------------------------------------------------------------

void MirroredBuf::allocate(size_t dim0_, size_t dim1_) {
  COMET_INSIST(is_alias || !is_allocated);

  if (use_linalg_) {

    gm_linalg_malloc(this, dim0_, dim1_, &env_);

    active = env_.is_compute_method_gpu() ? d : h;
    is_alias = false;
    is_allocated = true;

  } else {

    MirroredBuf::allocate(dim0_, dim1_, env_.matrix_buf_elt_size());

  } // if (use_linalg_)
}

//-----------------------------------------------------------------------------

void MirroredBuf::allocate(MirroredBuf& buf, size_t dim0_) {
  COMET_INSIST(is_alias || !is_allocated);
  COMET_INSIST(dim0_ <= buf.dim0);
  COMET_INSIST(!buf.is_alias);
  COMET_INSIST(buf.is_allocated);

  h = buf.h;
  d = buf.d;
  active = buf.active;
  dim0 = dim0_;
  dim1 = buf.dim1;
  size = buf.size;
  is_alias = true;
  is_allocated = true;
}

//-----------------------------------------------------------------------------

void MirroredBuf::set_zero_h() {
  COMET_INSIST(is_allocated);

  for (size_t i=0; i<size; ++i) {
    ((char*)h)[i] = 0;
  }
}

//-----------------------------------------------------------------------------

void MirroredBuf::deallocate() {
  COMET_INSIST(!is_locked_h_ && !is_locked_d_);

  if (is_allocated && !is_alias) {

    if (use_linalg_) {

      gm_linalg_free(this, &env_);

    } else {

#     if defined COMET_USE_CUDA
        cudaFreeHost(h);
        if (env_.is_compute_method_gpu())
          cudaFree(d);
#     elif defined COMET_USE_HIP
        hipHostFree(h);
        if (env_.is_compute_method_gpu())
          hipFree(d);
#     else
        free(h);
        COMET_INSIST(!env_.is_compute_method_gpu() &&
          "GPU not supported for this build.");
#     endif

      env_.cpu_mem_local_dec(size);
      if (env_.is_compute_method_gpu())
        env_.gpu_mem_local_dec(size);

      h = NULL;
      d = NULL;

    } // if (use_linalg_)

    active = NULL;
    is_allocated = false;

  } // if
}

//-----------------------------------------------------------------------------

void MirroredBuf::to_accel_start() {

//FIX ...
  if (env_.is_compute_method_gpu())
    lock();
  else
    return;

  if (use_linalg_) {

    gm_linalg_set_matrix_start(this, &env_);

  } else {

#   if defined COMET_USE_CUDA
      cudaMemcpyAsync(d, h, size, cudaMemcpyHostToDevice, env_.stream_togpu());
#   elif defined COMET_USE_HIP
      hipMemcpyAsync(d, h, size, hipMemcpyHostToDevice, env_.stream_togpu());
#   endif

  } // if (use_linalg_)
}

//-----------------------------------------------------------------------------

void MirroredBuf::to_accel_wait() {

  if (use_linalg_)
    gm_linalg_set_matrix_wait(&env_);
  else
    env_.stream_synchronize(env_.stream_togpu());

  if (env_.is_compute_method_gpu())
    unlock();
  else
    return;
}

//-----------------------------------------------------------------------------

void MirroredBuf::to_accel() {
  to_accel_start();
  to_accel_wait();
}

//-----------------------------------------------------------------------------

void MirroredBuf::from_accel_start() {

  if (env_.is_compute_method_gpu())
    lock();
  else
    return;

  if (use_linalg_) {

    gm_linalg_get_matrix_start(this, &env_);

  } else {

#   if defined COMET_USE_CUDA
      cudaMemcpyAsync(h, d, size, cudaMemcpyDeviceToHost,
        env_.stream_fromgpu());
#   elif defined COMET_USE_HIP
      hipMemcpyAsync(h, d, size, hipMemcpyDeviceToHost, env_.stream_fromgpu());
#   endif

  } // if (use_linalg_)
}

//-----------------------------------------------------------------------------

void MirroredBuf::from_accel_wait() {

  if (use_linalg_)
    gm_linalg_get_matrix_wait(&env_);
  else
    env_.stream_synchronize(env_.stream_fromgpu());

  if (env_.is_compute_method_gpu())
    unlock();
  else
    return;
}

//-----------------------------------------------------------------------------

void MirroredBuf::from_accel() {
  from_accel_start();
  from_accel_wait();
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
