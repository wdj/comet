//-----------------------------------------------------------------------------
/*!
 * \file   mirrored_buf.hh
 * \author Wayne Joubert
 * \date   Thu Aug  3 15:04:05 EDT 2017
 * \brief  Struct/code to manage dual CPU/GPU reflected arrays.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_mirrored_buf_hh_
#define _comet_mirrored_buf_hh_

#include "env.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

class GMMirroredBuf {
public:

  GMMirroredBuf(Env& env);
  GMMirroredBuf(size_t dim0, size_t dim1, Env& env);
  GMMirroredBuf(GMMirroredBuf& buf, size_t dim0, Env& env);

  ~GMMirroredBuf();

  void allocate(size_t dim0, size_t dim1);
  void allocate(GMMirroredBuf& buf, size_t dim0);
  void deallocate();

  template<typename T>
  T& elt(int i0, int i1) {
    COMET_ASSERT(is_allocated);
    COMET_ASSERT(i0 >= 0 && (size_t)i0 < dim0);
    COMET_ASSERT(i1 >= 0 && (size_t)i1 < dim1);

    return ((T*)(h))[i0 + dim0 * i1];
  }

  template<typename T>
  T elt_const(int i0, int i1) const {
    COMET_ASSERT(is_allocated);
    COMET_ASSERT(i0 >= 0 && (size_t)i0 < dim0);
    COMET_ASSERT(i1 >= 0 && (size_t)i1 < dim1);

    return ((T*)(h))[i0 + dim0 * i1];
  }

  void to_accel_start();
  void to_accel_wait();
  void to_accel();

  void from_accel_start();
  void from_accel_wait();
  void from_accel();

  void* __restrict__ h;
  void* __restrict__ d;
  void* __restrict__ active;
  size_t dim0;
  size_t dim1;
  size_t size;
  bool is_alias;
  bool is_allocated;

private:

  Env& env_;
};

// TODO: put copy to host / copy to device fns here

//-----------------------------------------------------------------------------

//GMMirroredBuf GMMirroredBuf_null(void);

void GMMirroredBuf_create(GMMirroredBuf* p, size_t dim0, size_t dim1,
                          GMEnv* env);

void GMMirroredBuf_create(GMMirroredBuf* p, GMMirroredBuf* p_old, size_t dim0,
                          GMEnv* env);

void GMMirroredBuf_destroy(GMMirroredBuf* p, GMEnv* env);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_mirrored_buf_hh_

//-----------------------------------------------------------------------------
