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

class MirroredBuf {
private:

  CEnv& env_;

public:

  MirroredBuf(CEnv& env);
  MirroredBuf(size_t dim0, size_t dim1, int elt_size, CEnv& env);
  MirroredBuf(size_t dim0, size_t dim1, CEnv& env);
  MirroredBuf(MirroredBuf& buf, size_t dim0, CEnv& env);

  ~MirroredBuf();

  void allocate(size_t dim0, size_t dim1, int elt_size);
  void allocate(size_t dim0, size_t dim1);
  void allocate(MirroredBuf& buf, size_t dim0);
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

  void set_zero_h();

  void to_accel_start();
  void to_accel_wait();
  void to_accel();

  void from_accel_start();
  void from_accel_wait();
  void from_accel();

  // TODO: make private; ? change active to a
  void* __restrict__ h;
  void* __restrict__ d;
  void* __restrict__ active;
  size_t dim0;
  size_t dim1;
  size_t num_elts_;
  int elt_size_;
  size_t size_;
  bool is_alias;
  bool is_allocated;

  size_t num_elts() const {return num_elts_;}

  void lock_h() const {
    COMET_INSIST(!is_locked_h_);
    is_locked_h_ = true;
  }

  void lock_d() const {
    COMET_INSIST(!is_locked_d_);
    is_locked_d_ = true;
  }

  void unlock_h() const {
    COMET_INSIST(is_locked_h_);
    is_locked_h_ = false;
  }

  void unlock_d() const {
    COMET_INSIST(is_locked_d_);
    is_locked_d_ = false;
  }

  void lock() const {
    lock_h();
    lock_d();
  }

  void unlock() const {
    unlock_h();
    unlock_d();
  }

private:

  mutable bool is_locked_h_;
  mutable bool is_locked_d_;
  bool use_linalg_;

  // Disallowed methods.
  MirroredBuf(const MirroredBuf&);
  void operator=(const MirroredBuf&);

  friend void gm_linalg_malloc(MirroredBuf* p, size_t dim0, size_t dim1, CEnv* env);
  friend void gm_linalg_free(MirroredBuf* p, CEnv* env);

};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_mirrored_buf_hh_

//-----------------------------------------------------------------------------
