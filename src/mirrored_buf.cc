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
  , env_(env) {
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
  , env_(env) {

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
  , env_(env) {
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

  gm_linalg_malloc(this, dim0, dim1, &env_);
  active = env_.compute_method() == ComputeMethod::GPU ? d : h;
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

  if (is_allocated && !is_alias)
    gm_linalg_free(this, &env_);

  is_allocated = false;
}

//-----------------------------------------------------------------------------






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

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
