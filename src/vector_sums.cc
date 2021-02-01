//-----------------------------------------------------------------------------
/*!
 * \file   vector_sums.cc
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

#include "cstdint"
#include "cstdlib"

#include "env.hh"
#include "vectors.hh"
#include "vector_sums.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

VectorSums::VectorSums(size_t num_vector_local, CEnv& env)
  : env_(env)
  , num_vector_local_(num_vector_local)
  , sums_(env)
  , sums_tmp_(env)
  , counts_(env)
  , counts_tmp_(env) {

  allocate_();
}

//-----------------------------------------------------------------------------

VectorSums::VectorSums(const GMVectors& vectors, CEnv& env)
  : env_(env)
  , num_vector_local_(vectors.num_vector_local)
  , sums_(env)
  , sums_tmp_(env)
  , counts_(env)
  , counts_tmp_(env) {

  allocate_();
  compute(vectors);
}

//-----------------------------------------------------------------------------

void VectorSums::allocate_() {

  if (!num_vector_local_)
    return;

  sums_.allocate(num_vector_local_, 1, sizeof(Float_t));
  if (env_.do_reduce())
    sums_tmp_.allocate(num_vector_local_, 1, sizeof(Float_t));

  if (env_.sparse() && env_.is_metric_type_bitwise()) {
    counts_.allocate(num_vector_local_, 1, sizeof(Float_t));
    if (env_.do_reduce())
      counts_tmp_.allocate(num_vector_local_, 1, sizeof(Float_t));
  }
}

//-----------------------------------------------------------------------------

void VectorSums::compute(const GMVectors& vectors) {
  COMET_INSIST(num_vector_local_ == vectors.num_vector_local);

  // Compute the vector sums.

  if (env_.is_metric_type_bitwise()) {
    compute_bits2_(vectors);
  } else {
    compute_float_(vectors);
  }

  // Put on accelerator in case needed later.

  sums_.to_accel();

  if (env_.sparse() && env_.is_metric_type_bitwise())
    counts_.to_accel();
}

//-----------------------------------------------------------------------------

void VectorSums::compute_float_(const GMVectors& vectors) {
  COMET_INSIST(num_vector_local_ == vectors.num_vector_local);

  MirroredBuf& sums_local = env_.do_reduce() ? sums_tmp_ : sums_;

  // Sum up all values in each vector.

# pragma omp parallel for schedule(dynamic,1000)
  for (int i = 0; i < num_vector_local_; ++i) {
    GMFloat sum = 0;
    //#pragma omp parallel for reduction(+:sum)
    for (int f = 0; f < vectors.num_field_local; ++f) {
      const Float_t value = GMVectors_float_get(&vectors, f, i, &env_);
      sum += value;
    }
    elt_ref_(sums_local, i) = sum;
  }

  // Do reduction across field procs if needed.

  if (env_.do_reduce())
    COMET_MPI_SAFE_CALL(MPI_Allreduce(sums_local.h, sums_.h,
      num_vector_local_, COMET_MPI_FLOAT, MPI_SUM, env_.comm_field()));

  env_.ops_local_inc(2 * num_vector_local_ * (double)vectors.num_field_local);
}

//-----------------------------------------------------------------------------

void VectorSums::compute_bits2_(const GMVectors& vectors) {
  COMET_INSIST(num_vector_local_ == vectors.num_vector_local);

  MirroredBuf& sums_local = env_.do_reduce() ? sums_tmp_ : sums_;
  MirroredBuf& counts_local = env_.do_reduce() ? counts_tmp_ : counts_;

  // Count number of 1-bits in each vector

  const int cbpe = env_.counted_bits_per_elt();
  const bool is_cbpe_2 = cbpe == 2;

  //----------
  if (env_.compute_method() == ComputeMethod::REF) {
  //----------

#   pragma omp parallel for schedule(dynamic,1000)
    for (int i = 0; i < num_vector_local_; ++i) {
      Float_t sum = 0;
      if (env_.sparse()) {
        Float_t count = 0;
        for (int f = 0; f < (int)vectors.dm->num_field_active_local; ++f) {
          // Slow way: sum each seminibble individually
          const GMBits2 v = GMVectors_bits2_get(&vectors, f, i, &env_);
          if (GM_2BIT_UNKNOWN != v){
            sum += is_cbpe_2 ? ((v & 1) != 0) + ((v & 2) != 0)
                             : ((v & 1) != 0);
            count++;
          }
        }
        COMET_ASSERT(count >= 0 && count <= vectors.dm->num_field_active_local);
        COMET_ASSERT(sum >= 0 && sum <= cbpe * count);
        elt_ref_(sums_local, i) = sum;
        elt_ref_(counts_local, i) = count;
      } else { // ! sparse
        //#pragma omp parallel for reduction(+:sum)
        for (int f = 0; f < vectors.num_field_local; ++f) {
          // Slow way: sum each seminibble individually
          const GMBits2 v = GMVectors_bits2_get(&vectors, f, i, &env_);
          sum += is_cbpe_2 ? ((v & 1) != 0) + ((v & 2) != 0)
                           : ((v & 1) != 0);
        }
        COMET_ASSERT(sum >= 0 &&
                     sum <= cbpe * vectors.dm->num_field_active_local);
        elt_ref_(sums_local, i) = sum;
      } // if sparse
    } // for i

    //----------
  } else { // env_.compute_method() != ComputeMethod::REF
    //----------

    const int num_pad_field_local = vectors.dm->num_pad_field_local;

#   pragma omp parallel for schedule(dynamic,1000)
    for (int i = 0; i < num_vector_local_; ++i) {
      Float_t sum = 0;
      if (env_.sparse()) {
        const uint64_t oddbits = 0x5555555555555555;
        Float_t count = 0;
        for (int f = 0; f < vectors.num_packedfield_local; ++f) {
          // Fast way: sum all 64 bits of each word immediately
          const GMBits2x64 v = GMVectors_bits2x64_get(&vectors, f, i, &env_);
          const uint64_t data0 = v.data[0];
          const uint64_t data1 = v.data[1];
          const uint64_t v10_oddmask0 = (data0 | ~(data0 >> 1)) & oddbits;
          const uint64_t v10_oddmask1 = (data1 | ~(data1 >> 1)) & oddbits;
          const uint64_t v10_mask0 = v10_oddmask0 | (v10_oddmask0 << 1);
          const uint64_t v10_mask1 = v10_oddmask1 | (v10_oddmask1 << 1);
          sum += is_cbpe_2 ? utils::popc64(data0 & v10_mask0)
                           : utils::popc64(data0 & oddbits & v10_mask0);
          sum += is_cbpe_2 ? utils::popc64(data1 & v10_mask1)
                           : utils::popc64(data1 & oddbits & v10_mask1);
          // NOTE: the code below interlaces half the bits of each of the two
          // 64-bit words being processed here.
          // In fact, "count" counts the VECTOR ELEMENTS that are defined, not
          // the number of BITS for all the defined elements.
          count += utils::popc64(v10_oddmask0 | (v10_oddmask1 << 1));
        }
        // Adjust for end pad
        count -= num_pad_field_local;
        // Finish
        COMET_ASSERT(count >= 0 && count <= vectors.dm->num_field_active_local);
        COMET_ASSERT(sum >= 0 && sum <= cbpe * count);
        elt_ref_(sums_local, i) = sum;
        elt_ref_(counts_local, i) = count;
      } else { // ! sparse
        const uint64_t oddbits = 0x5555555555555555;
        for (int f = 0; f < vectors.num_packedfield_local; ++f) {
          // Fast way: sum all 64 bits of each word immediately
          const GMBits2x64 v = GMVectors_bits2x64_get(&vectors, f, i, &env_);
          sum += is_cbpe_2 ? utils::popc64(v.data[0])
                            : utils::popc64(v.data[0] & oddbits);
          sum += is_cbpe_2 ? utils::popc64(v.data[1])
                            : utils::popc64(v.data[1] & oddbits);
          // NOTE: for this case pad entries of vec all zero so no effect on sum
        }
        COMET_ASSERT(sum >= 0 &&
                     sum <= cbpe * vectors.dm->num_field_active_local);
        elt_ref_(sums_local, i) = sum;
      } // if sparse
    } // for i

    //----------
  } // if
  //----------

  // Do reduction across field procs if needed

  if (env_.do_reduce()) {
    COMET_MPI_SAFE_CALL(MPI_Allreduce(sums_local.h, sums_.h,
      num_vector_local_, COMET_MPI_FLOAT, MPI_SUM, env_.comm_field()));
    if (env_.sparse())
      COMET_MPI_SAFE_CALL(MPI_Allreduce(counts_local.h, counts_.h,
        num_vector_local_, COMET_MPI_FLOAT, MPI_SUM, env_.comm_field()));
  }
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
