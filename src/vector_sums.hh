//-----------------------------------------------------------------------------
/*!
 * \file   vector_sums.hh
 * \author Wayne Joubert
 * \date   Sat Oct 24 10:48:32 EDT 2015
 * \brief  Compute the denominators needed by the methods.
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

#ifndef _COMET_VECTOR_SUMS_HH_
#define _COMET_VECTOR_SUMS_HH_

#include "env.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

class VectorSums {

public:

  typedef GMFloat Float_t;

  VectorSums(size_t num_vector_local, CEnv& env_);
  VectorSums(const Vectors& vectors, CEnv& env);
  ~VectorSums() {}

  void compute(const Vectors& vectors);

  void compute_accel(const Vectors& vectors, AccelStream_t accel_stream);

  void to_accel() {
    sums_.to_accel();
    if (need_counts_())
      counts_.to_accel();
  }

  void from_accel() {
    sums_.from_accel();
    if (need_counts_())
      counts_.from_accel();
  }

  Float_t sum(size_t i) const {
    COMET_ASSERT(i < sums_.dim0); // && i >= 0
    return ((Float_t*)sums_.h)[i];
  }

  Float_t count(size_t i) const {
    COMET_ASSERT(need_counts_());
    COMET_ASSERT(i < counts_.dim0); // && i >= 0
    return ((Float_t*)counts_.h)[i];
  }

  MirroredBuf* sums() {return &sums_;}
  MirroredBuf* counts() {return &counts_;}

private:

  CEnv& env_;

  const int num_vector_local_;

  MirroredBuf sums_;
  MirroredBuf sums_tmp_;
  MirroredBuf counts_;
  MirroredBuf counts_tmp_;

  void compute_float_(const Vectors& vectors);
  void compute_bits2_(const Vectors& vectors);
  void compute_float_accel_(const Vectors& vectors,
    AccelStream_t accel_stream);
  void compute_bits2_accel_(const Vectors& vectors,
    AccelStream_t accel_stream);

  Float_t& elt_ref_(MirroredBuf& b, size_t i) {
    COMET_ASSERT(i < b.dim0); // && i >= 0
    return (((Float_t*)b.h)[i]);
  }

  void allocate_();

  bool need_counts_() const {
    return env_.sparse() && env_.is_metric_type_bitwise();
  }

  // Disallowed methods.

  VectorSums(  const VectorSums&);
  void operator=(const VectorSums&);

}; // VectorSums

//-----------------------------------------------------------------------------
/// \brief Utility class for aggregating vector-related objects.

// TODO: put in separate file.

struct VData {
  Vectors* vectors;
  MirroredBuf* buf;
  VectorSums* sums;

  VData(Vectors* vectors_in, MirroredBuf* buf_in,
    VectorSums* sums_in = NULL)
    : vectors(vectors_in), buf(buf_in), sums(sums_in) {}
};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_VECTOR_SUMS_HH_

//-----------------------------------------------------------------------------
