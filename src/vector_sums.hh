//-----------------------------------------------------------------------------
/*!
 * \file   vector_sums.hh
 * \author Wayne Joubert
 * \date   Sat Oct 24 10:48:32 EDT 2015
 * \brief  Compute the denominators needed by the methods.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_vector_sums_hh_
#define _comet_vector_sums_hh_

#include "env.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

class VectorSums {

  typedef GMFloat Float_t;

public:

  VectorSums(size_t num_vector_local, Env& env_);
  VectorSums(const GMVectors& vectors, Env& env);
  ~VectorSums() {}

  void compute(const GMVectors& vectors);

  Float_t sum_const(size_t i) const {
    COMET_ASSERT(i < sums_.dim0); // && i >= 0
    return ((Float_t*)sums_.h)[i];
  }

  Float_t count_const(size_t i) const {
    COMET_ASSERT(i < counts_.dim0); // && i >= 0
    return ((Float_t*)counts_.h)[i];
  }

private:

  Env& env_;

  const int num_vector_local_;

  GMMirroredBuf sums_;
  GMMirroredBuf sums_tmp_;
  GMMirroredBuf counts_;
  GMMirroredBuf counts_tmp_;

  void compute_float(const GMVectors& vectors);
  void compute_bits2(const GMVectors& vectors);

  Float_t& sum(size_t i) {
    COMET_ASSERT(i < sums_.dim0); // && i >= 0
    return (((Float_t*)sums_.h)[i]);
  }

  Float_t& count(size_t i) {
    COMET_ASSERT(i < counts_.dim0); // && i >= 0
    return (((Float_t*)counts_.h)[i]);
  }

  void allocate_();

  // Disallowed methods.

  VectorSums(  const VectorSums&);
  void operator=(const VectorSums&);

}; // VectorSums

//-----------------------------------------------------------------------------







//-----------------------------------------------------------------------------
// Struct declaration

typedef struct {
  GMFloat* __restrict__ sums;
  GMFloat* __restrict__ counts;
  GMFloat* __restrict__ sums_tmp_;
  GMFloat* __restrict__ counts_tmp_;
  size_t size_;
  size_t num_field_;
} GMVectorSums;

//=============================================================================
// Null object

GMVectorSums GMVectorSums_null(void);

//=============================================================================
// Pseudo-constructor

void GMVectorSums_create(GMVectorSums* this_,
                         int num_vector_local,
                         GMEnv* env);

//=============================================================================
// Pseudo-destructor

void GMVectorSums_destroy(GMVectorSums* this_, GMEnv* env);

//=============================================================================
// Compute

void GMVectorSums_compute(GMVectorSums* this_, GMVectors* vectors, GMEnv* env);

//=============================================================================
// Accessors

GMFloat GMVectorSums_sum(const GMVectorSums* this_, int i,  GMEnv* env);

GMFloat GMVectorSums_count(const GMVectorSums* this_, int i,  GMEnv* env);





//-----------------------------------------------------------------------------
/// \brief Utility class for aggregating vector-related objects.

// TODO: put in separate file.

struct VData {
  GMVectors* vectors;
  GMMirroredBuf* buf;
  GMVectorSums* sums;

  VData(GMVectors* vectors_in, GMMirroredBuf* buf_in,
    GMVectorSums* sums_in = NULL)
    : vectors(vectors_in), buf(buf_in), sums(sums_in) {}
};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_vector_sums_hh_

//-----------------------------------------------------------------------------
