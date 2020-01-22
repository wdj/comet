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

GMMirroredBuf::GMMirroredBuf(Env& env)
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

GMMirroredBuf::GMMirroredBuf(size_t dim0_, size_t dim1_, Env& env)
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

GMMirroredBuf::GMMirroredBuf(GMMirroredBuf& buf, size_t dim0_, Env& env)
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

GMMirroredBuf::~GMMirroredBuf() {
  deallocate();
}

//-----------------------------------------------------------------------------

void GMMirroredBuf::allocate(size_t dim0, size_t dim1) {
  COMET_INSIST(is_alias || !is_allocated);

  if (use_linalg_) {

    gm_linalg_malloc(this, dim0, dim1, &env_);

  } else {

    this->dim0 = dim0;
    this->dim1 = dim1;
    size = dim0 * dim1 * env_.matrix_buf_elt_size();
#   if defined COMET_USE_CUDA
      cudaMallocHost((void**)&h, size);
      if (env_.is_compute_method_gpu())
        cudaMalloc((void**)&d, size);
#   elif defined COMET_USE_HIP
      hipHostMalloc((void**)&h, size);
      if (env_.is_compute_method_gpu())
        hipMalloc((void**)&d, size);
#   else
      h = malloc(size);
      COMET_INSIST(!env_.is_compute_method_gpu() &&
        "GPU not supported for this build.");
#   endif

    env_.cpu_mem_local_inc(size);
    if (env_.is_compute_method_gpu())
      env_.gpu_mem_local_inc(size);

    COMET_INSIST(h &&
      "Invalid host pointer created, possibly due to insufficient memory.");
    COMET_INSIST((d || !env_.is_compute_method_gpu()) &&
      "Invalid device pointer created, possibly due to insufficient memory.");

  } // if (use_linalg_)

  active = env_.is_compute_method_gpu() ? d : h;
  is_alias = false;
  is_allocated = true;
}

//-----------------------------------------------------------------------------

void GMMirroredBuf::allocate(GMMirroredBuf& buf, size_t dim0_) {
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

void GMMirroredBuf::deallocate() {

  COMET_INSIST(!is_locked_h_ && !is_locked_d_);

  if (is_allocated && !is_alias)
    gm_linalg_free(this, &env_);

  is_allocated = false;
}

//-----------------------------------------------------------------------------

void GMMirroredBuf::to_accel_start() {
  if (env_.is_compute_method_gpu())
    lock();
  gm_linalg_set_matrix_start(this, &env_);
}

//-----------------------------------------------------------------------------

void GMMirroredBuf::to_accel_wait() {
  gm_linalg_set_matrix_wait(&env_);
  if (env_.is_compute_method_gpu())
    unlock();
}

//-----------------------------------------------------------------------------

void GMMirroredBuf::to_accel() {
  to_accel_start();
  to_accel_wait();
}

//-----------------------------------------------------------------------------

void GMMirroredBuf::from_accel_start() {
  if (env_.is_compute_method_gpu())
    lock();
  gm_linalg_get_matrix_start(this, &env_);
}

//-----------------------------------------------------------------------------

void GMMirroredBuf::from_accel_wait() {
  gm_linalg_get_matrix_wait(&env_);
  if (env_.is_compute_method_gpu())
    unlock();
}

//-----------------------------------------------------------------------------

void GMMirroredBuf::from_accel() {
  from_accel_start();
  from_accel_wait();
}

//-----------------------------------------------------------------------------








#if 0
//-----------------------------------------------------------------------------
/// \brief Mirrored buf pseudo-constructor: allocate CPU and GPU arrays.

void GMMirroredBuf_create(GMMirroredBuf* p, size_t dim0, size_t dim1, 
                          GMEnv* env) {
  COMET_INSIST(p && env);
  COMET_INSIST(dim0 + 1 >= 1 && dim1 + 1 >= 1);

  gm_linalg_malloc(p, dim0, dim1, env);
  p->active = env->compute_method() == ComputeMethod::GPU ? p->d : p->h;
  p->is_allocated = true;
}

//-----------------------------------------------------------------------------
/// \brief Mirrored buf pseudo-constructor: create alias to existing mirror.

void GMMirroredBuf_create(GMMirroredBuf* p, GMMirroredBuf* p_old, size_t dim0,
                          GMEnv* env) {
  COMET_INSIST(p && p_old && env);
  COMET_INSIST(dim0 <= p_old->dim0);
  COMET_INSIST(p_old->is_allocated);

  p->h = p_old->h;
  p->d = p_old->d;
  p->size = p_old->size;
  p->dim0 = dim0;
  p->dim1 = p_old->dim1;
  p->is_alias = true;
  p->active = p_old->active;
  p->is_allocated = true;
}

//-----------------------------------------------------------------------------
/// \brief Mirrored buf pseudo-destructor

void GMMirroredBuf_destroy(GMMirroredBuf* p, GMEnv* env) {
  COMET_INSIST(p && env);

  if (! p->is_alias) {
    gm_linalg_free(p, env);
  }

  p->is_allocated = false;
}
#endif

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
