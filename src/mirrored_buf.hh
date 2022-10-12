//-----------------------------------------------------------------------------
/*!
 * \file   mirrored_buf.hh
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

#ifndef _COMET_MIRRORED_BUF_HH_
#define _COMET_MIRRORED_BUF_HH_

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
  void terminate() {
    deallocate();
  }

  template<typename T>
  T& elt(size_t ind0, size_t ind1) {
    COMET_ASSERT(is_allocated);
    COMET_ASSERT(ind0+1 >= 0+1 && ind0 < dim0);
    COMET_ASSERT(ind1+1 >= 0+1 && ind1 < dim1);

    return ((T*)(h))[ind0 + dim0 * ind1];
  }

  template<typename T>
  T elt_const(size_t ind0, size_t ind1) const {
    COMET_ASSERT(is_allocated);
    COMET_ASSERT(ind0+1 >= 0+1 && ind0 < dim0);
    COMET_ASSERT(ind1+1 >= 0+1 && ind1 < dim1);

    return ((T*)(h))[ind0 + dim0 * ind1];
  }

  void set_zero_h();

  void to_accel_start();
  void to_accel_wait();
  void to_accel();

  void from_accel_start(AccelStream_t stream);
  void from_accel_wait(AccelStream_t stream);
  void from_accel(AccelStream_t stream);

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
  size_t size_allocated_;
  bool is_alias;
  bool is_allocated;

  size_t num_elts() const {return num_elts_;}
  size_t size() const {return num_elts_ * elt_size_;}

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

  bool is_compute_method_gpu() const {
    return env_.is_compute_method_gpu();
  }

private:

  mutable bool is_locked_h_;
  mutable bool is_locked_d_;
  bool use_linalg_;

  friend class MagmaWrapper;

  // Disallowed methods.
  MirroredBuf(const MirroredBuf&);
  void operator=(const MirroredBuf&);
};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_MIRRORED_BUF_HH_

//-----------------------------------------------------------------------------
