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

MirroredBuf::MirroredBuf(size_t dim0_, size_t dim1_, Env& env)
  : h(NULL)
  , d(NULL)
  , active(NULL)
  , dim0(dim0_)
  , dim1(dim1_)
  , size(dim0_ * dim1_)
  , is_alias(false)
  , env_(env) {

#if 0
  gm_linalg_malloc(this, dim0, dim1, &env_);
#endif
  active = env_.compute_method() == ComputeMethod::GPU ? d : h;
}

//-----------------------------------------------------------------------------

MirroredBuf::MirroredBuf(MirroredBuf& b_old, size_t dim0_, Env& env)
  : h(b_old.h)
  , d(b_old.d)
  , active(b_old.active)
  , dim0(dim0_)
  , dim1(b_old.dim1)
  , size(b_old.size) //FIX
  , is_alias(true)
  , env_(env) {
}

//-----------------------------------------------------------------------------

MirroredBuf::~MirroredBuf() {

#if 0
  if (!is_alias)
    gm_linalg_free(this, &env_);
#endif
}

//-----------------------------------------------------------------------------








//-----------------------------------------------------------------------------

GMMirroredBuf GMMirroredBuf_null(void) {
  GMMirroredBuf p;
  p.h = NULL;
  p.d = NULL;
  p.active = NULL;
  p.size = 0;
  p.dim0 = 0;
  p.dim1 = 0;
  p.is_alias = false;
  return p;
}

//-----------------------------------------------------------------------------
/// \brief Mirrored buf pseudo-constructor: allocate CPU and GPU arrays.

void GMMirroredBuf_create(GMMirroredBuf* p, size_t dim0, size_t dim1, 
                          GMEnv* env) {
  COMET_INSIST(p && env);
  COMET_INSIST(dim0 + 1 >= 1 && dim1 + 1 >= 1);

  gm_linalg_malloc(p, dim0, dim1, env);
  p->active = env->compute_method() == ComputeMethod::GPU ? p->d : p->h;
}

//-----------------------------------------------------------------------------
/// \brief Mirrored buf pseudo-constructor: create alias to existing mirror.

void GMMirroredBuf_create(GMMirroredBuf* p, GMMirroredBuf* p_old, size_t dim0,
                          GMEnv* env) {

  COMET_INSIST(p && p_old && env);
  COMET_INSIST(dim0 <= p_old->dim0);

  p->h = p_old->h;
  p->d = p_old->d;
  p->size = p_old->size;
  p->dim0 = dim0;
  p->dim1 = p_old->dim1;
  p->is_alias = true;
  p->active = p_old->active;
}

//-----------------------------------------------------------------------------
/// \brief Mirrored buf pseudo-destructor

void GMMirroredBuf_destroy(GMMirroredBuf* p, GMEnv* env) {
  COMET_INSIST(p && env);

  if (! p->is_alias) {
    gm_linalg_free(p, env);
  }

  *p = GMMirroredBuf_null();
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
