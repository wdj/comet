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

GMMirroredBuf GMMirroredBuf_null(void) {
  GMMirroredBuf p;
  p.h = NULL;
  p.d = NULL;
  p.size = 0;
  p.dim0 = 0;
  p.dim1 = 0;
  p.is_alias = false;
  return p;
}

//-----------------------------------------------------------------------------

void GMMirroredBuf_create(GMMirroredBuf* p, size_t dim0, size_t dim1, 
                          GMEnv* env) {
  GMInsist(p && env);
  GMInsist(dim0 + 1 >= 1 && dim1 + 1 >= 1);

  gm_linalg_malloc(p, dim0, dim1, env);
}

//-----------------------------------------------------------------------------

void GMMirroredBuf_create(GMMirroredBuf* p, GMMirroredBuf* p_old, size_t dim0,
                          GMEnv* env) {

  GMInsist(p && p_old && env);
  GMInsist(dim0 <= p_old->dim0);

  p->h = p_old->h;
  p->d = p_old->d;
  p->size = p_old->size;
  p->dim0 = dim0;
  p->dim1 = p_old->dim1;
  p->is_alias = true;
}

//-----------------------------------------------------------------------------

void GMMirroredBuf_destroy(GMMirroredBuf* p, GMEnv* env) {
  GMInsist(p && env);

  if (! p->is_alias) {
    gm_linalg_free(p, env);
  }

  *p = GMMirroredBuf_null();
}

//=============================================================================

//-----------------------------------------------------------------------------
