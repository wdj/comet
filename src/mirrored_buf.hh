//-----------------------------------------------------------------------------
/*!
 * \file   mirrored_buf.hh
 * \author Wayne Joubert
 * \date   Thu Aug  3 15:04:05 EDT 2017
 * \brief  CPU/GPU mirrored buffer, header.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_mirrored_buf_hh_
#define _gm_mirrored_buf_hh_

#include "env.hh"

//=============================================================================

typedef struct {
  void* __restrict__ h;
  void* __restrict__ d;
  size_t size;
  size_t dim0;
  size_t dim1;
  bool is_alias;
} GMMirroredBuf;

//-----------------------------------------------------------------------------

GMMirroredBuf GMMirroredBuf_null(void);

void GMMirroredBuf_create(GMMirroredBuf* p, size_t dim0, size_t dim1,
                          GMEnv* env);

void GMMirroredBuf_create(GMMirroredBuf* p, GMMirroredBuf* p_old, size_t dim0,
                          GMEnv* env);

void GMMirroredBuf_destroy(GMMirroredBuf* p, GMEnv* env);

//-----------------------------------------------------------------------------

template<typename T>
static T& GMMirroredBuf_elt(GMMirroredBuf* p, int i0, int i1) {
  GMAssert(p);
  GMAssert(i0 >= 0 && (size_t)i0 < p->dim0);
  GMAssert(i1 >= 0 && (size_t)i1 < p->dim1);

  return ((T*)(p->h))[i0 + p->dim0 * i1];
}

//-----------------------------------------------------------------------------

template<typename T>
static T GMMirroredBuf_elt_const(const GMMirroredBuf* p, int i0, int i1) {
  GMAssert(p);
  GMAssert(i0 >= 0 && (size_t)i0 < p->dim0);
  GMAssert(i1 >= 0 && (size_t)i1 < p->dim1);

  return ((T*)(p->h))[i0 + p->dim0 * i1];
}

//=============================================================================

#endif // _gm_mirrored_buf_hh_

//-----------------------------------------------------------------------------
