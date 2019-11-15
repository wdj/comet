//-----------------------------------------------------------------------------
/*!
 * \file   mirrored_buf.hh
 * \author Wayne Joubert
 * \date   Thu Aug  3 15:04:05 EDT 2017
 * \brief  Struct/code to manage dual CPU/GPU reflected arrays.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_mirrored_buf_hh_
#define _gm_mirrored_buf_hh_

#include "env.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

typedef struct {
  void* __restrict__ h;
  void* __restrict__ d;
  void* __restrict__ active;
  size_t size;
  size_t dim0;
  size_t dim1;
  bool is_alias;
} GMMirroredBuf;

// TODO: is it appropriate to put copy to host / copy to device fns here

//-----------------------------------------------------------------------------

GMMirroredBuf GMMirroredBuf_null(void);

void GMMirroredBuf_create(GMMirroredBuf* p, size_t dim0, size_t dim1,
                          GMEnv* env);

void GMMirroredBuf_create(GMMirroredBuf* p, GMMirroredBuf* p_old, size_t dim0,
                          GMEnv* env);

void GMMirroredBuf_destroy(GMMirroredBuf* p, GMEnv* env);

//-----------------------------------------------------------------------------
/// \brief Mirrored buf element accessor.

template<typename T>
static T& GMMirroredBuf_elt(GMMirroredBuf* p, int i0, int i1) {
  COMET_ASSERT(p);
  COMET_ASSERT(i0 >= 0 && (size_t)i0 < p->dim0);
  COMET_ASSERT(i1 >= 0 && (size_t)i1 < p->dim1);

  return ((T*)(p->h))[i0 + p->dim0 * i1];
}

//-----------------------------------------------------------------------------
/// \brief Mirrored buf const element accessor.

template<typename T>
static T GMMirroredBuf_elt_const(const GMMirroredBuf* p, int i0, int i1) {
  COMET_ASSERT(p);
  COMET_ASSERT(i0 >= 0 && (size_t)i0 < p->dim0);
  COMET_ASSERT(i1 >= 0 && (size_t)i1 < p->dim1);

  return ((T*)(p->h))[i0 + p->dim0 * i1];
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _gm_mirrored_buf_hh_

//-----------------------------------------------------------------------------
