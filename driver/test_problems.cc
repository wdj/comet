//-----------------------------------------------------------------------------
/*!
 * \file   test_problems.cc
 * \author Wayne Joubert
 * \date   Mon Aug  7 17:02:51 EDT 2017
 * \brief  Generator for synthetic test problems.
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

#include "cstdio"
#include "tuple"

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"

#include "driver.hh"
#include "vectors_io.hh"
#include "test_problems.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Set the entries of the vectors.

void set_vectors_random_(GMVectors* vectors, int verbosity, CEnv* env) {
  COMET_INSIST(vectors && env);

  if (! env->is_proc_active()) {
    return;
  }

  const size_t nva = vectors->dm->num_vector_active;
  const size_t nfa = vectors->dm->num_field_active;

  switch (env->data_type_vectors()) {
    //--------------------
    case GM_DATA_TYPE_FLOAT: {
    //--------------------
#pragma omp parallel for
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        size_t vector = vl +
            vectors->num_vector_local * (size_t)env->proc_num_vector();
        // Fill pad vectors with copies of the last vector.
        // By construction, active vectors are packed for lower procs.
        const size_t vector_capped = utils::min(vector, nva);
        int fl = 0;
        for (fl = 0; fl < vectors->num_field_local; ++fl) {
          size_t field = fl +
              vectors->num_field_local * (size_t)env->proc_num_field();
          if (field >= vectors->num_field_active) {
            continue; // These vector entries will be padded to zero elsewhere.
          }
          // Compute element unique id.
          const size_t uid = field + nfa * vector_capped;
          // Generate large random number.
          size_t rand1 = uid;
          rand1 = utils::randomize(rand1);
          rand1 = utils::randomize(rand1);
          size_t rand2 = uid;
          rand2 = utils::randomize(rand2);
          rand2 = utils::randomize(rand2);
          rand2 = utils::randomize(rand2);
          const size_t rand_max = utils::randomize_max();
          size_t rand_value = rand1 + rand_max * rand2;
          /*---Reduce so that after summing num_field times the integer
               still exactly representable by floating point type---*/
          const size_t rand_max2 = rand_max * rand_max;
          const int log2_num_summands_3way_numer = 2;
          const int shift_amount1 = utils::max(0,
             utils::log2(log2_num_summands_3way_numer * rand_max2 * nfa)
             - mantissa_digits<GMFloat>() + 1);
          // Account for cast to float in magma Volta version.
          const int shift_amount2 = utils::max(0,
                             utils::log2(rand_max2) - mantissa_digits<float>() + 1);
          const int shift_amount = utils::max(shift_amount1, shift_amount2);
          rand_value >>= shift_amount > 0 ? shift_amount : 0;
          // Store.
          GMFloat float_value = (GMFloat)rand_value;
          COMET_INSIST((size_t)float_value == rand_value);
          COMET_INSIST(float_value * vectors->num_field_active <
                         ((size_t)1)<<mantissa_digits<GMFloat>());
          GMVectors_float_set(vectors, fl, vl, float_value, env);
        } // field_local
      }   // vector_local
      // Print.
//TODO: move this
      //if (verbosity > 2) {
      //  VectorsIO::print(*vectors, *env);
      //}
    } break;
    //--------------------
    case GM_DATA_TYPE_BITS2: {
    //--------------------
#pragma omp parallel for
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {

        size_t vector = vl +
            vectors->num_vector_local * (size_t)env->proc_num_vector();
        // Fill pad vectors with copies of the last vector.
        const size_t vector_capped = utils::min(vector, nva);

        for (int fl = 0; fl < vectors->num_field_local; ++fl) {
          size_t field = fl +
              vectors->num_field_local * (size_t)env->proc_num_field();
          if (field >= vectors->num_field_active) {
            continue; // These vector entries will be padded to zero elsewhere.
          }
          // Compute element unique id.
          const size_t uid = field + vectors->num_field_active * vector_capped;
          size_t index = uid;
          // Randomize.
          index = utils::randomize(index);
          index = utils::randomize(index);
          // Calculate random number between 0 and 3.
          const float float_rand_value = index / (float)utils::randomize_max();
          // Create 2-bit value - make extra sure less than 4.
          GMBits2 value = (int)((4. - 1e-5) * float_rand_value);
          // Store.
          GMVectors_bits2_set(vectors, fl, vl, value, env);
        } // fl
      }   // vl
      // Print.
//TODO: move this
      //if (verbosity > 2) {
      //  VectorsIO::print(*vectors, *env);
      //}
    } break;
    //--------------------
    default:
    //--------------------
      COMET_INSIST(false && "Invalid data type.");
  } // switch
}

//-----------------------------------------------------------------------------

static size_t perm_shuffle(size_t key, size_t i, size_t n) {
  COMET_ASSERT((key & (~(size_t)1)) == 0);
  COMET_ASSERT(i>=0 && i<n);
  COMET_ASSERT(n>=0);

  // For an integer between 0 and n-1, permute it to another such integer.
  // The permutation choice is specified by 1 bit.
  // For an ascending sequence of integers, first output the even values,
  // then the odd values, for key=0. If key=1, same with even/odd reversed.

  const size_t nhalf = (n+1-key)/2;
  const size_t result = i < nhalf ? 2*i + key : 2*(i-nhalf) + 1 - key;
  COMET_ASSERT(result>=0 && result<n);
  return result;
}

//-----------------------------------------------------------------------------

// This is not a totally (pseudo)random permutation.  However it does
// have the advantage that it can be computed formulaically and quickly.

#ifndef __clang__
#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
#endif

static size_t perm(size_t key, size_t i, size_t n) {
  COMET_ASSERT((key & (~(size_t)((1<<TestProblemInfo::NUM_SHUFFLE)-1))) == 0);
  COMET_ASSERT(i>=0 && i<n);
  COMET_ASSERT(n>=0);

  // For an integer between 0 and n-1, permute it to another such integer.
  // The permutation choice is specified by NUM_SHUFFLE bits.

  size_t result = i;
  size_t key_resid = key;
  for (int shuffle_num = 0; shuffle_num < TestProblemInfo::NUM_SHUFFLE; ++shuffle_num) {
    result = perm_shuffle(key_resid&1, result, n);
    key_resid >>= 1;
  }
  COMET_ASSERT(result>=0 && result<n);
  return result;
}

#ifndef __clang__
#pragma GCC pop_options
#endif

//-----------------------------------------------------------------------------

void set_vectors_analytic_(GMVectors* vectors, int verbosity, CEnv* env) {
  COMET_INSIST(vectors && env);

  if (! env->is_proc_active())
    return;

  const size_t nfa = vectors->num_field_active;
  const size_t nva = vectors->dm->num_vector_active;

  const auto tpi = TestProblemInfo(nva, nfa, *env);

  // The elements of a single permuted vector are partitioned into
  // "groups", with all elements in a group contiguous and having
  // the same value.
  // By keeping the number of groups (here = 8) much smaller than
  // the vector length, the calculation of the exact comparisons
  // is much cheaper -- the comparison of 2 or 3 vectors by element
  // is the same across all elements of the group.

  switch (env->data_type_vectors()) {
    //--------------------
    case GM_DATA_TYPE_FLOAT: {
    //--------------------
#pragma omp parallel for
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        size_t vector = vl +
            vectors->num_vector_local * (size_t)env->proc_num_vector();
        // Fill pad vectors with copies of the last vector.
        const size_t vector_capped = utils::min(vector, nva-1);
        for (int fl = 0; fl < vectors->num_field_local; ++fl) {
          size_t field = fl +
              vectors->num_field_local * (size_t)env->proc_num_field();
          if (field >= nfa) {
            continue; // These vector entries will be padded to zero elsewhere.
          }
          const size_t f = field; // field number
          const size_t v = vector_capped; // vector number

          const size_t pf = perm(0, f, nfa); // permuted field number
          const size_t g = pf / tpi.group_size_max_; // group number
          COMET_ASSERT(g>=0 && g<tpi.num_group_);

          const size_t pv = perm(g, v, nva); // permuted vector number

          // Linearly map pv to small interval.
          const size_t value = tpi.value_min_ + (pv * tpi.value_max_) / (tpi.value_min_+nva);

          const GMFloat float_value = value;

          // Store.
          COMET_INSIST(float_value * nfa >= 1);
          COMET_INSIST(float_value * nfa < tpi.max_float_);
          GMVectors_float_set(vectors, fl, vl, float_value, env);

        } // field_local
      }   // vector_local
    } break;
    //--------------------
    case GM_DATA_TYPE_BITS2: {
    //--------------------

#pragma omp parallel for
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        size_t vector = vl +
            vectors->num_vector_local * (size_t)env->proc_num_vector();
        // Fill pad vectors with copies of the last vector.
        const size_t vector_capped = utils::min(vector, nva-1);
        for (int fl = 0; fl < vectors->num_field_local; ++fl) {
          size_t field = fl +
              vectors->num_field_local * (size_t)env->proc_num_field();
          if (field >= nfa) {
            continue; // These vector entries will be padded to zero elsewhere.
          }
          // Create 2-bit value - make extra sure less than 4.

          const size_t f = field;
          const size_t v = vector_capped;

          const size_t pf = perm(0, f, nfa);
          const size_t g = pf / tpi.group_size_max_;
          COMET_ASSERT(g>=0 && g<tpi.num_group_);

          const size_t pv = perm(g, v, nva);

          const size_t value = tpi.value_min_ + ( pv * tpi.value_max_ ) / (nva+tpi.value_min_);

          const GMBits2 bval = ((size_t)3) & (value - tpi.value_min_);

          // Store.
          GMVectors_bits2_set(vectors, fl, vl, bval, env);

        } // field_local
      }   // vector_local
    } break;
    //--------------------
    default:
    //--------------------
      COMET_INSIST(false && "Invalid data type.");
  } // switch
}

//=============================================================================

void TestProblem::set_vectors_synthetic(GMVectors* vectors, int problem_type,
                                        int verbosity, CEnv* env) {
  COMET_INSIST(vectors && env);

  if (problem_type == ProblemType::RANDOM) {
    set_vectors_random_(vectors, verbosity, env);
  } else if (problem_type == ProblemType::ANALYTIC) {
    set_vectors_analytic_(vectors, verbosity, env);
  } else {
    COMET_INSIST(false && "Invalid problem_type");
  }
}

//=============================================================================
// Compute metric value, analytic case, czek 2-way.

static GMFloat metric_value_analytic_(size_t vi,
  size_t vj, const TestProblemInfo& tpi) {

  GMFloat float_n = 0;
  GMFloat float_d = 0;

  size_t n = 0;
  size_t d = 0;

  // For each comparison of vectors, the compared/summed
  // elements are treated as num_group groups.  All element
  // comparisons in the group have the same value, so we just
  // compute once and multiply that by the group size.

  for (size_t g=0; g<tpi.num_group_; ++g) {

    const size_t pf_min = g * tpi.group_size_max_;
    const size_t pf_max = utils::min((g+1) * tpi.group_size_max_, tpi.nfa_);
    const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

    const size_t pvi = perm(g, vi, tpi.nva_);
    const size_t pvj = perm(g, vj, tpi.nva_);

    const size_t value_i = tpi.value_min_ + ( pvi * tpi.value_max_ ) /
                                       (tpi.value_min_+tpi.nva_);
    const size_t value_j = tpi.value_min_ + ( pvj * tpi.value_max_ ) /
                                       (tpi.value_min_+tpi.nva_);
    float_n += utils::min(value_i, value_j) * gs_this;
    float_d += (value_i + value_j) * gs_this;
    n += utils::min(value_i, value_j) * gs_this;
    d += (value_i + value_j) * gs_this;

  } //---g

  COMET_INSIST(n == (size_t)float_n);
  COMET_INSIST(d == (size_t)float_d);

  const GMFloat multiplier = (GMFloat)2;

  const GMFloat value = (multiplier * float_n) / float_d;

  return value;
}

//=============================================================================
// Compute metric value, analytic case, czek 3-way.

static GMFloat metric_value_analytic_(size_t vi,
  size_t vj, size_t vk, const TestProblemInfo& tpi) {

  GMFloat float_n = 0;
  GMFloat float_d = 0;

  size_t n = 0;
  size_t d = 0;

  // For each comparison of vectors, the compared/summed
  // elements are treated as num_group groups.  All element
  // comparisons in the group have the same value, so we just
  // compute once and multiply that by the group size.

  for (size_t g=0; g<tpi.num_group_; ++g) {

    const size_t pf_min = g * tpi.group_size_max_;
    const size_t pf_max = utils::min((g+1) * tpi.group_size_max_, tpi.nfa_);
    const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

    const size_t pvi = perm(g, vi, tpi.nva_);
    const size_t pvj = perm(g, vj, tpi.nva_);
    const size_t pvk = perm(g, vk, tpi.nva_);

    const size_t value_i = tpi.value_min_ + ( pvi * tpi.value_max_ ) /
                                            (tpi.nva_+tpi.value_min_);
    const size_t value_j = tpi.value_min_ + ( pvj * tpi.value_max_ ) /
                                            (tpi.nva_+tpi.value_min_);
    const size_t value_k = tpi.value_min_ + ( pvk * tpi.value_max_ ) /
                                            (tpi.nva_+tpi.value_min_);

    float_n += utils::min(value_i, value_j) * gs_this;
    float_n += utils::min(value_i, value_k) * gs_this;
    float_n += utils::min(value_j, value_k) * gs_this;
    float_n -= utils::min(value_i, utils::min(value_j, value_k)) * gs_this;
    float_d += (value_i + value_j + value_k) * gs_this;

    n += utils::min(value_i, value_j) * gs_this;
    n += utils::min(value_i, value_k) * gs_this;
    n += utils::min(value_j, value_k) * gs_this;
    n -= utils::min(value_i, utils::min(value_j, value_k)) * gs_this;
    d += (value_i + value_j + value_k) * gs_this;

  } //---g

  COMET_INSIST(n == (size_t)float_n);
  COMET_INSIST(d == (size_t)float_d);

  const GMFloat multiplier = (GMFloat)1.5;

  const GMFloat value = (multiplier * float_n) / float_d;

  return value;
}

//=============================================================================
// Compute metric value, analytic case, ccc/duo 2-way.

static GMFloat metric_value_analytic_(size_t vi,
  size_t vj, int iE, int jE, const TestProblemInfo& tpi, CEnv& env) {

  const int cbpe = env.counted_bits_per_elt();
  const double recip_m = 1. / tpi.nfa_;

  const size_t iG = vi;
  const size_t jG = vj;

  GMTally1 rij = 0;
  GMTally1 si = 0;
  GMTally1 sj = 0;
  GMTally1 ci = 0;
  GMTally1 cj = 0;
  GMTally1 cij = 0;

  for (size_t g=0; g<tpi.num_group_; ++g) {

    const size_t pf_min = g * tpi.group_size_max_;
    const size_t pf_max = utils::min((g+1) * tpi.group_size_max_, tpi.nfa_);
    const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

    const size_t piG = perm(g, iG, tpi.nva_);
    const size_t pjG = perm(g, jG, tpi.nva_);

    const size_t value_i = tpi.value_min_ + (piG * tpi.value_max_) /
                                            (tpi.nva_+tpi.value_min_);
    const size_t value_j = tpi.value_min_ + (pjG * tpi.value_max_ ) /
                                            (tpi.nva_+tpi.value_min_);

    const GMBits2 bval_i = ((size_t)3) & (value_i - tpi.value_min_);
    const GMBits2 bval_j = ((size_t)3) & (value_j - tpi.value_min_);

    const int bval_i_0 = !!(bval_i&1);
    const int bval_i_1 = !!(bval_i&2);
    const int bval_j_0 = !!(bval_j&1);
    const int bval_j_1 = !!(bval_j&2);

    const bool unk_i = env.sparse() && bval_i == GM_2BIT_UNKNOWN;
    const bool unk_j = env.sparse() && bval_j == GM_2BIT_UNKNOWN;
    const bool unk_ij = unk_i || unk_j;

    if (! unk_i) {
      ci += gs_this;
      si += cbpe == 2 ?
        ((bval_i_0 == iE) + (bval_i_1 == iE)) * gs_this :
        (bval_i_0 == iE) * gs_this;
    }

    if (! unk_j) {
      cj += gs_this;
      sj += cbpe == 2 ?
        ((bval_j_0 == jE) + (bval_j_1 == jE)) * gs_this :
        (bval_j_0 == jE) * gs_this;
    }

    if (! unk_ij) {
      cij += cbpe * cbpe * gs_this;
      rij += cbpe == 2 ?
             (((bval_i_0 == iE) && (bval_j_0 == jE)) +
              ((bval_i_0 == iE) && (bval_j_1 == jE)) +
              ((bval_i_1 == iE) && (bval_j_0 == jE)) +
              ((bval_i_1 == iE) && (bval_j_1 == jE))) *
             gs_this :
             ((bval_i_0 == iE) && (bval_j_0 == jE)) *
             gs_this;
    }
  } //---g

  double value = 0;
  const bool is_zero_denom = ci == 0 || cj == 0 || cij == 0;
  if (!is_zero_denom) {
    // CHECK typing here
    const double f_one = 1;

    const double f_ci = (double) ci;
    const double f_cj = (double) cj;

    const double f_cicj_min = f_ci < f_cj ? f_ci : f_cj;
    const double f_cicj_max = f_ci > f_cj ? f_ci : f_cj;

    const double f_cij = (double) cij;
    const double recip_cicjcij = f_one /
                                  (f_cicj_min * f_cicj_max * f_cij);

    const double recip_ci = env.sparse() ?
      f_cj * f_cij * recip_cicjcij : recip_m;
    const double recip_cj = env.sparse() ?
      f_ci * f_cij * recip_cicjcij : recip_m;

    const double recip_sumcij = env.sparse() ?
      f_cicj_min * f_cicj_max * recip_cicjcij :
      (f_one / (cbpe * cbpe)) * recip_m;

    value = cbpe == CBPE::CCC ?
      ccc_duo_value<CBPE::CCC>(rij, si, sj,
          recip_ci, recip_cj, recip_sumcij,
          env_ccc_duo_multiplier<CBPE::CCC>(env), env.ccc_param()) :
      ccc_duo_value<CBPE::DUO>(rij, si, sj,
          recip_ci, recip_cj, recip_sumcij,
          env_ccc_duo_multiplier<CBPE::DUO>(env), env.ccc_param());
  } // is_zero_denom

  return value;
}

//=============================================================================
// Compute metric value, analytic case, ccc/duo 3-way.

static GMFloat metric_value_analytic_(size_t vi,
  size_t vj, size_t vk, int iE, int jE, int kE, const TestProblemInfo& tpi,
  CEnv& env) {

  const int cbpe = env.counted_bits_per_elt();
  const double recip_m = 1. / tpi.nfa_;

  const size_t iG = vi;
  const size_t jG = vj;
  const size_t kG = vk;

    GMTally1 rijk = 0;
    GMTally1 si = 0;
    GMTally1 sj = 0;
    GMTally1 sk = 0;
    GMTally1 ci = 0;
    GMTally1 cj = 0;
    GMTally1 ck = 0;
    GMTally1 cijk = 0;

    for (size_t g=0; g<tpi.num_group_; ++g) {

      const size_t pf_min = g * tpi.group_size_max_;
      const size_t pf_max = utils::min((g+1) * tpi.group_size_max_, tpi.nfa_);
      const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

      const size_t piG = perm(g, iG, tpi.nva_);
      const size_t pjG = perm(g, jG, tpi.nva_);
      const size_t pkG = perm(g, kG, tpi.nva_);

      const size_t value_i = tpi.value_min_ + (piG * tpi.value_max_ ) /
                                         (tpi.nva_+tpi.value_min_);
      const size_t value_j = tpi.value_min_ + (pjG * tpi.value_max_ ) /
                                         (tpi.nva_+tpi.value_min_);
      const size_t value_k = tpi.value_min_ + (pkG * tpi.value_max_ ) /
                                         (tpi.nva_+tpi.value_min_);

      const GMBits2 bval_i = ((size_t)3) & (value_i - tpi.value_min_);
      const GMBits2 bval_j = ((size_t)3) & (value_j - tpi.value_min_);
      const GMBits2 bval_k = ((size_t)3) & (value_k - tpi.value_min_);

      const int bval_i_0 = !!(bval_i&1);
      const int bval_i_1 = !!(bval_i&2);
      const int bval_j_0 = !!(bval_j&1);
      const int bval_j_1 = !!(bval_j&2);
      const int bval_k_0 = !!(bval_k&1);
      const int bval_k_1 = !!(bval_k&2);


      const bool unk_i = env.sparse() && bval_i == GM_2BIT_UNKNOWN;
      const bool unk_j = env.sparse() && bval_j == GM_2BIT_UNKNOWN;
      const bool unk_k = env.sparse() && bval_k == GM_2BIT_UNKNOWN;
      const bool unk_ijk = unk_i || unk_j || unk_k;

      if (! unk_i) {
        ci += gs_this;
        si += cbpe == 2 ?
          ((bval_i_0 == iE) + (bval_i_1 == iE)) * gs_this :
          (bval_i_0 == iE) * gs_this;
      }

      if (! unk_j) {
        cj += gs_this;
        sj += cbpe == 2 ?
          ((bval_j_0 == jE) + (bval_j_1 == jE)) * gs_this :
          (bval_j_0 == jE) * gs_this;
      }

      if (! unk_k) {
        ck += gs_this;
        sk += cbpe == 2 ?
          ((bval_k_0 == kE) + (bval_k_1 == kE)) * gs_this :
          (bval_k_0 == kE) * gs_this;
      }

      if (! unk_ijk) {
        cijk += cbpe * cbpe * cbpe * gs_this;
        rijk += cbpe == 2 ?
                (((bval_i_0==iE) && (bval_j_0==jE) && (bval_k_0==kE))+
                 ((bval_i_1==iE) && (bval_j_0==jE) && (bval_k_0==kE))+
                 ((bval_i_0==iE) && (bval_j_1==jE) && (bval_k_0==kE))+
                 ((bval_i_1==iE) && (bval_j_1==jE) && (bval_k_0==kE))+
                 ((bval_i_0==iE) && (bval_j_0==jE) && (bval_k_1==kE))+
                 ((bval_i_1==iE) && (bval_j_0==jE) && (bval_k_1==kE))+
                 ((bval_i_0==iE) && (bval_j_1==jE) && (bval_k_1==kE))+
                 ((bval_i_1==iE) && (bval_j_1==jE) && (bval_k_1==kE)))
               * gs_this :
               ((bval_i_0 == iE) && (bval_j_0 == jE) &&
                (bval_k_0 == kE)) * gs_this;
      }
    } //---g

    double value = 0;
    const bool is_zero_denom = ci == 0 || cj == 0 || ck == 0 || cijk == 0;

    if (!is_zero_denom) {
      // CHECK typing here
      const double f_one = 1;

      const double recip_ci = env.sparse() ? f_one/ci : recip_m;
      const double recip_cj = env.sparse() ? f_one/cj : recip_m;
      const double recip_ck = env.sparse() ? f_one/ck : recip_m;

      const double recip_sumcijk = env.sparse() ? f_one/cijk :
                                     (f_one / 8) * recip_m;

      value = cbpe == CBPE::CCC ?
        ccc_duo_value<CBPE::CCC>(rijk, si, sj, sk,
                 recip_ci, recip_cj, recip_ck, recip_sumcijk, env) :
        ccc_duo_value<CBPE::DUO>(rijk, si, sj, sk,
                 recip_ci, recip_cj, recip_ck, recip_sumcijk, env);

    } // is_zero_denom

  return value;
}

//=============================================================================
// Check correctness of metrics, if possible.

void TestProblem::check_metrics_analytic_(GMMetrics* metrics, Driver& driver,
                                          CEnv* env) {
  COMET_INSIST(metrics && env);
  COMET_INSIST(ProblemType::ANALYTIC == driver.options_.problem_type);
  COMET_INSIST(NULL == driver.options_.input_file);

  if (! env->is_proc_active())
    return;

  const size_t nfa = metrics->num_field_active;
  const size_t nva = metrics->num_vector_active;

  const auto tpi = TestProblemInfo(nva, nfa, *env);

  size_t num_incorrect = 0;
  const size_t max_to_print = 10;
  double max_incorrect_diff = 0.;

  const size_t hnlen = 256;
  char hn[hnlen];
  gethostname(hn, hnlen);
  int rank = 0;
  COMET_MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

  switch (env->data_type_metrics()) {
    //--------------------
    case GM_DATA_TYPE_FLOAT: {
    //--------------------
      if (env->num_way() == NumWay::_2) {

#       pragma omp parallel for reduction(+:num_incorrect) reduction(max:max_incorrect_diff)
        for (size_t index = 0; index < metrics->num_metrics_local; ++index) {
          const size_t vi = Metrics_coords_getG(*metrics, index, 0, *env);
          const size_t vj = Metrics_coords_getG(*metrics, index, 1, *env);
          if (vi >= nva || vj >= nva)
            continue;

          const auto value = Metrics_elt_const<GMFloat>(*metrics, index, *env);

          GMFloat value_expected = metric_value_analytic_(vi, vj, tpi);

          const bool is_incorrect = value_expected != value;
          if (is_incorrect) {
            const double diff = value - value_expected;
            max_incorrect_diff = utils::max(fabs(diff), max_incorrect_diff);
            if (num_incorrect < max_to_print)
              fprintf(stderr, "Error: incorrect result detected.  "
                     "coords %zu %zu  "
                     "expected %.20e  actual %.20e  diff %.20e\n", vi, vj,
                     (double)value_expected, (double)value, diff);
          }

          num_incorrect += is_incorrect;
        } //---for index

      } //---if
      if (env->num_way() == NumWay::_3) {

#       pragma omp parallel for reduction(+:num_incorrect) reduction(max:max_incorrect_diff)
        for (size_t index = 0; index < metrics->num_metrics_local; ++index) {
          const size_t vi = Metrics_coords_getG(*metrics, index, 0, *env);
          const size_t vj = Metrics_coords_getG(*metrics, index, 1, *env);
          const size_t vk = Metrics_coords_getG(*metrics, index, 2, *env);
          if (vi >= nva || vj >= nva || vk >= nva)
            continue;

          const auto value = Metrics_elt_const<GMFloat>(*metrics, index, *env);

          GMFloat value_expected = metric_value_analytic_(vi, vj, vk, tpi);

          const bool is_incorrect = value_expected != value;
          if (is_incorrect) {
            const double diff = value - value_expected;
            max_incorrect_diff = utils::max(fabs(diff), max_incorrect_diff);
            if (num_incorrect < max_to_print)
              fprintf(stderr, "Error: incorrect result detected.  "
                      "coords %zu %zu %zu  "
                      "expected %.20e  actual %.20e  diff %.20e\n", vi, vj, vk,
                      (double)value_expected, (double)value, diff);
          }

          num_incorrect += is_incorrect;
        } //---for index
      } //---if

    } break;
    //--------------------
    case GM_DATA_TYPE_TALLY2X2: {
    //--------------------

#     pragma omp parallel for reduction(+:num_incorrect) reduction(max:max_incorrect_diff)
      for (size_t index = 0; index < metrics->num_metric_items_local_computed;
           ++index) {
        const MetricItemCoords_t coords = metrics->coords_value(index);
        const size_t iG = CoordsInfo::getiG(coords, *metrics, *env);
        const size_t jG = CoordsInfo::getjG(coords, *metrics, *env);
        if (iG >= nva || jG >= nva)
          continue;
        for (int entry_num = 0; entry_num < env->num_entries_per_metric_item();
             ++entry_num) {
          const int iE = CoordsInfo::getiE(coords, entry_num, *metrics, *env);
          const int jE = CoordsInfo::getjE(coords, entry_num, *metrics, *env);

          // Get actual.
          const GMFloat value = Metrics_ccc_duo_get_2(*metrics, index,
            entry_num, *env);
          //const bool is_pass_threshold
          //  = Metrics_is_pass_threshold(*metrics, index, iE, jE, *env);

          // Get expected.
          const GMFloat value_expected_nothreshold
            = metric_value_analytic_(iG, jG, iE, jE, tpi, *env);

          // if metric may be thresholded to zero before now, then check.
          // NOTE: will never be here if is_shrink and entry failed threshold.
          // NOTE: if not is_thresold_tc, disregard thresholding issue
          // since actual is not thresholded yet.

          bool do_set_zero = false;
          if (env->is_threshold_tc()) {

            // NOTE: using MF::SINGLE if and only if is_threshold_tc()
            typedef Tally2x2<MetricFormat::SINGLE> TTable_t;
            TTable_t ttable = TTable_t::null();

            for (int iE_ = 0; iE_ < 2; ++iE_) {
              for (int jE_ = 0; jE_ < 2; ++jE_) {
                const GMFloat mva
                  = metric_value_analytic_(iG, jG, iE_, jE_, tpi, *env);
                // NOTE: GMFloat == float if we are here.
                TTable_t::set(ttable, iE_, jE_, mva);
              }
            }
            do_set_zero = ! env->thresholds().is_pass(ttable, iE, jE);

          } // if (env->is_threshold_tc())

          const GMFloat value_expected = do_set_zero ?
            (GMFloat)0 : value_expected_nothreshold;

#if 0
//#ifdef COMET_USE_INT128
          if (env->are_ccc_params_default()) {
          if (!(0 == ci || 0 == cj || 0 == cij)) {
            value_expected = GMMetrics_ccc_value_nofp_2(metrics,
              rij, si, sj, ci, cj, cij, env); 
          }
          }
#endif

          const bool is_incorrect = value_expected != value;
          if (is_incorrect) {
            const double diff = value - value_expected;
            max_incorrect_diff = utils::max(fabs(diff), max_incorrect_diff);
            if (num_incorrect < max_to_print)
              fprintf(stderr, "Error: incorrect result detected.  "
                     "coords %zu %zu  %i %i  "
                     "expected %.20e  actual %.20e  diff %.20e\n",
                     iG, jG, iE, jE,
                     (double)value_expected, (double)value, diff);
          } // is_correct

          num_incorrect += is_incorrect;
        } // for entry_num
      } // for index

    } break;
    //--------------------
    case GM_DATA_TYPE_TALLY4X2: {
    //--------------------

#     pragma omp parallel for reduction(+:num_incorrect) reduction(max:max_incorrect_diff)
      for (size_t index = 0; index < metrics->num_metric_items_local_computed;
           ++index) {
        const MetricItemCoords_t coords = metrics->coords_value(index);
        const size_t iG = CoordsInfo::getiG(coords, *metrics, *env);
        const size_t jG = CoordsInfo::getjG(coords, *metrics, *env);
        const size_t kG = CoordsInfo::getkG(coords, *metrics, *env);
        if (iG >= nva || jG >= nva || kG >= nva)
          continue;
        for (int entry_num = 0; entry_num < env->num_entries_per_metric_item();
             ++entry_num) {
          const int iE = CoordsInfo::getiE(coords, entry_num, *metrics, *env);
          const int jE = CoordsInfo::getjE(coords, entry_num, *metrics, *env);
          const int kE = CoordsInfo::getkE(coords, entry_num, *metrics, *env);

          // Get actual.
          const GMFloat value = Metrics_ccc_duo_get_3(*metrics, index,
            entry_num, *env);

          // Get expected.
          GMFloat value_expected_nothreshold
            = metric_value_analytic_(iG, jG, kG, iE, jE, kE, tpi, *env);

          // if metric may be thresholded to zero before now, then check.
          // NOTE: will never be here if is_shrink and entry failed threshold.
          // NOTE: if not is_thresold_tc, disregard thresholding issue
          // since actual is not thresholded yet.

          bool do_set_zero = false;
          if (env->is_threshold_tc()) {

            // NOTE: using MF::SINGLE if and only if is_threshold_tc()
            typedef Tally4x2<MetricFormat::SINGLE> TTable_t;
            TTable_t ttable = TTable_t::null();

            for (int iE_ = 0; iE_ < 2; ++iE_) {
              for (int jE_ = 0; jE_ < 2; ++jE_) {
                for (int kE_ = 0; kE_ < 2; ++kE_) {
                  const GMFloat mva
                    = metric_value_analytic_(iG, jG, kG, iE_, jE_, kE_, tpi, *env);
                  // NOTE: GMFloat == float if we are here.
                  TTable_t::set(ttable, iE_, jE_, kE_, mva);
                }
              }
            }
            do_set_zero = ! env->thresholds().is_pass(ttable, iE, jE, kE);

          } // if (env->is_threshold_tc())

          const GMFloat value_expected = do_set_zero ?
            (GMFloat)0 : value_expected_nothreshold;

//          // If threshold_tc, threshold this to match computed result.
//          const bool do_set_zero = env->is_threshold_tc() &&
//            //!env->pass_threshold((double)(float)value_expected_nothreshold);
//            !env->thresholds().is_pass((double)(float)value_expected_nothreshold);
//
//          GMFloat value_expected = do_set_zero ? 0. : value_expected_nothreshold;

#if 0
//#ifdef COMET_USE_INT128
         if (env->are_ccc_params_default()) {
          if (!(0 == ci || 0 == cj || 0 == ck || 0 == cijk)) {
            value_expected = GMMetrics_ccc_value_nofp_3(metrics,
              rijk, si, sj, sk, ci, cj, ck, cijk, env); 
          }
          }
#endif

          const bool is_incorrect = value_expected != value;
          if (is_incorrect) {
            const double diff = value - value_expected;
            max_incorrect_diff = utils::max(fabs(diff), max_incorrect_diff);
            if (num_incorrect < max_to_print) {
              fprintf(stderr, "Error: incorrect result detected.  "
                     "coords %zu %zu %zu  %i %i %i  "
                     "expected %.20e  actual %.20e  diff %.20e  exp2 %.20e  hostname  %s  rank  %i\n",
                     iG, jG, kG, iE, jE, kE,
                     (double)value_expected, (double)value, diff,
                     value_expected_nothreshold, hn, rank);

            }

          } // is_correct

          num_incorrect += is_incorrect;
        } // for entry_num
      } // for index

    } break;
    //--------------------
    default:
      COMET_INSIST(false && "Invalid data type.");
  } // switch
  driver.counters_.num_incorrect += num_incorrect;
  driver.counters_.max_incorrect_diff = utils::max(max_incorrect_diff,
                                       driver.counters_.max_incorrect_diff);
}

//=============================================================================

void TestProblem::check_metrics(GMMetrics* metrics, Driver& driver,
                                CEnv* env) {
  COMET_INSIST(metrics && env);

  if (driver.options_.input_file)
    return;

  if (ProblemType::ANALYTIC == driver.options_.problem_type)
    check_metrics_analytic_(metrics, driver, env);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
