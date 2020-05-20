//-----------------------------------------------------------------------------
/*!
 * \file   test_problems.cc
 * \author Wayne Joubert
 * \date   Mon Aug  7 17:02:51 EDT 2017
 * \brief  Generator for synthetic test problems.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "cstdio"
//#include "stdlib.h"
//#include "stddef.h"
//#include "string.h"
//#include "float.h"
//#include "errno.h"

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
            continue; // These entries will be padded to zero elsewhere.
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
      if (verbosity > 2) {
        VectorsIO::print(*vectors, *env);
      }
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
            continue; // These entries will be padded to zero elsewhere.
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
      if (verbosity > 2) {
        VectorsIO::print(*vectors, *env);
      }
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

enum {NUM_SHUFFLE = 3};

// This is not a totally (pseudo)random permutation.  However it does
// have the advantage that it can be computed formulaically and quickly.

#ifndef __clang__
#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
#endif

static size_t perm(size_t key, size_t i, size_t n) {
  COMET_ASSERT((key & (~(size_t)((1<<NUM_SHUFFLE)-1))) == 0);
  COMET_ASSERT(i>=0 && i<n);
  COMET_ASSERT(n>=0);

  // For an integer between 0 and n-1, permute it to another such integer.
  // The permutation choice is specified by NUM_SHUFFLE bits.

  size_t result = i;
  size_t key_resid = key;
  for (int shuffle_num = 0; shuffle_num < NUM_SHUFFLE; ++shuffle_num) {
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

  if (! env->is_proc_active()) {
    return;
  }

  const size_t nfa = vectors->num_field_active;
  const size_t nva = vectors->dm->num_vector_active;

  // Upper bound on integer representable exactly by floating point type.
  // Account for cast to float in magma Volta version.
  const size_t max_float = ((size_t)1) <<
    (env->data_type_vectors() == GM_DATA_TYPE_FLOAT ?
     mantissa_digits<float>() : mantissa_digits<GMFloat>());
  // Czek account for number of terms summed in denom or num
  const size_t overflow_limit =
    env->data_type_vectors() != GM_DATA_TYPE_FLOAT ? 1 :
    env->num_way() == NUM_WAY::_2 ? 2 : 4;
  // Sum nfa times down the vector, is it still exact.
  const size_t value_limit = (max_float - 1) / (overflow_limit * nfa);

  const size_t value_min = 1;
  const size_t value_max = utils::min(value_min+nva, value_limit);

  // The elements of a single permuted vector are partitioned into
  // "groups", with all elements in a group contiguous and having
  // the same value.
  // By keeping the number of groups (here = 8) much smaller than
  // the vector length, the calculation of the exact comparisons
  // is much cheaper -- the comparison of 2 or 3 vectors by element
  // is the same across all elements of the group.

  const size_t num_group = 1 << NUM_SHUFFLE;
  const size_t group_size_max = utils::ceil(nfa, (size_t)num_group);

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
            continue; // These entries will be padded to zero elsewhere.
          }
          const size_t f = field; // field number
          const size_t v = vector_capped; // vector number

          const size_t pf = perm(0, f, nfa); // permuted field number
          const size_t g = pf / group_size_max; // group number
          COMET_ASSERT(g>=0 && g<num_group);

          const size_t pv = perm(g, v, nva); // permuted vector number

          // Linearly map pv to small interval.
          const size_t value = value_min + (pv * value_max) / (value_min+nva);

          const GMFloat float_value = value;

          // Store.
          COMET_INSIST(float_value * nfa >= 1);
          COMET_INSIST(float_value * nfa < max_float);
          GMVectors_float_set(vectors, fl, vl, float_value, env);

        } // field_local
      }   // vector_local
      // Print.
//TODO: move this
      if (verbosity > 2) {
        VectorsIO::print(*vectors, *env);
      }
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
            continue; // These entries will be padded to zero elsewhere.
          }
          // Create 2-bit value - make extra sure less than 4.

          const size_t f = field;
          const size_t v = vector_capped;

          const size_t pf = perm(0, f, nfa);
          const size_t g = pf / group_size_max;
          COMET_ASSERT(g>=0 && g<num_group);

          const size_t pv = perm(g, v, nva);

          const size_t value = value_min + ( pv * value_max ) / (nva+value_min);

          const GMBits2 bval = ((size_t)3) & (value - value_min);

          // Store.
          GMVectors_bits2_set(vectors, fl, vl, bval, env);

        } // field_local
      }   // vector_local
//TODO: move this
      if (verbosity > 2) {
        VectorsIO::print(*vectors, *env);
      }
    } break;
    //--------------------
    default:
    //--------------------
      COMET_INSIST(false && "Invalid data type.");
  } // switch
}

//=============================================================================

void set_vectors_synthetic(GMVectors* vectors, int problem_type, int verbosity,
                           CEnv* env) {
  COMET_INSIST(vectors && env);

  if (problem_type == GM_PROBLEM_TYPE_RANDOM) {
    set_vectors_random_(vectors, verbosity, env);
  } else if (problem_type == GM_PROBLEM_TYPE_ANALYTIC) {
    set_vectors_analytic_(vectors, verbosity, env);
  } else {
    COMET_INSIST(false && "Invalid problem_type");
  }
}

//=============================================================================
// Check correctness of metrics, if possible.

void check_metrics_analytic_(GMMetrics* metrics, DriverOptions* do_,
                             CEnv* env) {
  COMET_INSIST(metrics && do_ && env);
  COMET_INSIST(GM_PROBLEM_TYPE_ANALYTIC == do_->problem_type);
  COMET_INSIST(NULL == do_->input_file_path);

  if (! env->is_proc_active()) {
    return;
  }

  const size_t nfa = metrics->num_field_active;
  const size_t nva = metrics->num_vector_active;

  // Upper bound on integer representable exactly by floating point type.
  // Account for cast to float in magma Volta version.
  const size_t max_float = ((size_t)1) <<
    (env->data_type_vectors() == GM_DATA_TYPE_FLOAT ?
     mantissa_digits<float>() : mantissa_digits<GMFloat>());
  // Czek account for number of terms summed in denom or num
  const size_t overflow_limit =
    env->data_type_vectors() != GM_DATA_TYPE_FLOAT ? 1 :
    env->num_way() == NUM_WAY::_2 ? 2 : 4;
  // Sum nfa times down the vector, is it still exact.
  const size_t value_limit = (max_float - 1) / (overflow_limit * nfa);

  const size_t value_min = 1;
  const size_t value_max = utils::min(value_min+nva, value_limit);

  const size_t num_group = 1 << NUM_SHUFFLE;
  const size_t group_size_max = utils::ceil(nfa, (size_t)num_group);

  size_t num_incorrect = 0;
  const size_t max_to_print = 10;
  double max_incorrect_diff = 0.;

  switch (env->data_type_metrics()) {
    //--------------------
    case GM_DATA_TYPE_FLOAT: {
    //--------------------
      if (env->num_way() == NUM_WAY::_2) {
#pragma omp parallel for reduction(+:num_incorrect) reduction(max:max_incorrect_diff)
        for (size_t index = 0; index < metrics->num_elts_local; ++index) {
          const size_t vi =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          const size_t vj =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          if (vi >= nva || vj >= nva) {
            continue;
          }
          const auto value = Metrics_elt_const<GMFloat>(*metrics, index, *env);

          GMFloat float_n = 0;
          GMFloat float_d = 0;

          size_t n = 0;
          size_t d = 0;

          // For each comparison of vectors, the compared/summed
          // elements are treated as num_group groups.  All element
          // comparisons in the group have the same value, so we just
          // compute once and multiply that by the group size.

          for (size_t g=0; g<num_group; ++g) {

            const size_t pf_min = g * group_size_max;
            const size_t pf_max = utils::min((g+1) * group_size_max, nfa);
            const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

            const size_t pvi = perm(g, vi, nva);
            const size_t pvj = perm(g, vj, nva);

            const size_t value_i = value_min + ( pvi * value_max ) /
                                               (value_min+nva);
            const size_t value_j = value_min + ( pvj * value_max ) /
                                               (value_min+nva);
            float_n += utils::min(value_i, value_j) * gs_this;
            float_d += (value_i + value_j) * gs_this;
            n += utils::min(value_i, value_j) * gs_this;
            d += (value_i + value_j) * gs_this;

          } //---g

          COMET_INSIST(n == (size_t)float_n);
          COMET_INSIST(d == (size_t)float_d);

          const GMFloat multiplier = (GMFloat)2;

          const GMFloat value_expected = (multiplier * float_n) / float_d;

          const bool is_incorrect = value_expected != value;
          if (is_incorrect) {
            const double diff = fabs(value - value_expected);
            max_incorrect_diff = diff > max_incorrect_diff ? diff : max_incorrect_diff;
            if (num_incorrect < max_to_print) {
              fprintf(stderr, "Error: incorrect result detected.  coords %zu %zu  "
                     "expected %.20e  actual %.20e  diff %.20e\n", vi, vj,
                     (double)value_expected, (double)value,
                     (double)value-(double)value_expected);
            }
          }

          num_incorrect += is_incorrect;
        } //---for index
      } //---if
      if (env->num_way() == NUM_WAY::_3) {
#pragma omp parallel for reduction(+:num_incorrect) reduction(max:max_incorrect_diff)
        for (size_t index = 0; index < metrics->num_elts_local; ++index) {
          const size_t vi =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          const size_t vj =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          const size_t vk =
            GMMetrics_coord_global_from_index(metrics, index, 2, env);
          if (vi >= nva || vj >= nva || vk >= nva) {
            continue;
          }
          const auto value = Metrics_elt_const<GMFloat>(*metrics, index, *env);

          GMFloat float_n = 0;
          GMFloat float_d = 0;

          for (size_t g=0; g<num_group; ++g) {

            const size_t pf_min = g * group_size_max;
            const size_t pf_max = utils::min((g+1) * group_size_max, nfa);
            const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

            const size_t pvi = perm(g, vi, nva);
            const size_t pvj = perm(g, vj, nva);
            const size_t pvk = perm(g, vk, nva);

            const size_t value_i = value_min + ( pvi * value_max ) /
                                               (nva+value_min);
            const size_t value_j = value_min + ( pvj * value_max ) /
                                               (nva+value_min);
            const size_t value_k = value_min + ( pvk * value_max ) /
                                               (nva+value_min);

            float_n += utils::min(value_i, value_j) * gs_this;
            float_n += utils::min(value_i, value_k) * gs_this;
            float_n += utils::min(value_j, value_k) * gs_this;

            float_n -= utils::min(value_i, utils::min(value_j, value_k)) * gs_this;

            float_d += (value_i + value_j + value_k) * gs_this;

          } //---g

          const GMFloat multiplier = (GMFloat)1.5;

          const GMFloat value_expected = (multiplier * float_n) / float_d;

          const bool is_incorrect = value_expected != value;
          if (is_incorrect) {
            const double diff = fabs(value - value_expected);
            max_incorrect_diff = diff > max_incorrect_diff ? diff : max_incorrect_diff;
            if (num_incorrect < max_to_print) {
              fprintf(stderr, "Error: incorrect result detected.  coords %zu %zu %zu  "
                     "expected %.20e  actual %.20e  diff %.20e\n", vi, vj, vk,
                     (double)value_expected, (double)value,
                     (double)value-(double)value_expected);
            }
          }

          num_incorrect += is_incorrect;
        } //---for index
      } //---if
    } break;
    //--------------------
    case GM_DATA_TYPE_TALLY2X2: {
    //--------------------

    const int cbpe = env->counted_bits_per_elt();

#pragma omp parallel for reduction(+:num_incorrect) reduction(max:max_incorrect_diff)
      for (size_t index = 0; index < metrics->num_elts_local; ++index) {
        const size_t vi =
          GMMetrics_coord_global_from_index(metrics, index, 0, env);
        const size_t vj =
          GMMetrics_coord_global_from_index(metrics, index, 1, env);
        if (vi >= nva || vj >= nva) {
          continue;
        }
        for (int iE = 0; iE < 2; ++iE) {
          for (int jE = 0; jE < 2; ++jE) {
            const GMFloat value = Metrics_ccc_duo_get_2(*metrics,
              index, iE, jE, *env);

            GMTally1 rij = 0;
            GMTally1 si = 0;
            GMTally1 sj = 0;
            GMTally1 ci = 0;
            GMTally1 cj = 0;
            GMTally1 cij = 0;

            for (size_t g=0; g<num_group; ++g) {

              const size_t pf_min = g * group_size_max;
              const size_t pf_max = utils::min((g+1) * group_size_max, nfa);
              const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

              const size_t pvi = perm(g, vi, nva);
              const size_t pvj = perm(g, vj, nva);

              const size_t value_i = value_min + ( pvi * value_max ) /
                                                 (nva+value_min);
              const size_t value_j = value_min + ( pvj * value_max ) /
                                                 (nva+value_min);

              const GMBits2 bval_i = ((size_t)3) & (value_i - value_min);
              const GMBits2 bval_j = ((size_t)3) & (value_j - value_min);

              const int bval_i_0 = !!(bval_i&1);
              const int bval_i_1 = !!(bval_i&2);
              const int bval_j_0 = !!(bval_j&1);
              const int bval_j_1 = !!(bval_j&2);

              const bool unknown_i = env->sparse() && bval_i == GM_2BIT_UNKNOWN;
              const bool unknown_j = env->sparse() && bval_j == GM_2BIT_UNKNOWN;
              const bool unknown_ij = unknown_i || unknown_j;

              if (! unknown_i) {
                ci += gs_this;
                si += cbpe == 2 ?
                  ((bval_i_0 == iE) + (bval_i_1 == iE)) * gs_this :
                  (bval_i_0 == iE) * gs_this;
              }

              if (! unknown_j) {
                cj += gs_this;
                sj += cbpe == 2 ?
                  ((bval_j_0 == jE) + (bval_j_1 == jE)) * gs_this :
                  (bval_j_0 == jE) * gs_this;
              }

              if (! unknown_ij) {
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

            GMFloat value_expected_floatcalc = 0;
            if (!(ci == 0 || cj == 0 || cij == 0)) {
              // FIX typing here
              const double f_one = 1;

              const double f_ci = (double) ci;
              const double f_cj = (double) cj;

              const double f_cicj_min = f_ci < f_cj ? f_ci : f_cj;
              const double f_cicj_max = f_ci > f_cj ? f_ci : f_cj;

              const double f_cij = (double) cij;
              const double recip_cicjcij = f_one /
                                            (f_cicj_min * f_cicj_max * f_cij);

              const double recip_ci = env->sparse() ?
                f_cj * f_cij * recip_cicjcij : metrics->recip_m;
              const double recip_cj = env->sparse() ?
                f_ci * f_cij * recip_cicjcij : metrics->recip_m;

              const double recip_sumcij = env->sparse() ?
                f_cicj_min * f_cicj_max * recip_cicjcij :
                (f_one / (cbpe * cbpe)) * metrics->recip_m;

              value_expected_floatcalc = cbpe == 2 ?
                ccc_duo_value<CBPE::CCC>(rij, si, sj,
                    recip_ci, recip_cj, recip_sumcij,
                    env_ccc_duo_multiplier<CBPE::CCC>(*env), env->ccc_param()) :
                ccc_duo_value<CBPE::DUO>(rij, si, sj,
                    recip_ci, recip_cj, recip_sumcij,
                    env_ccc_duo_multiplier<CBPE::DUO>(*env), env->ccc_param());
            }

            GMFloat value_expected = value_expected_floatcalc;

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
              const double diff = fabs(value - value_expected);
              max_incorrect_diff = diff > max_incorrect_diff ?
                                   diff : max_incorrect_diff;
              if (num_incorrect < max_to_print) {
                fprintf(stderr, "Error: incorrect result detected.  coords %zu %zu  "
                       "expected %.20e  actual %.20e  diff %.20e\n", vi, vj,
                       (double)value_expected, (double)value,
                       (double)value-(double)value_expected);
              }
            }

            num_incorrect += is_incorrect;
          } //---j
        } //---i
      } //---for index
    } break;
    //--------------------
    case GM_DATA_TYPE_TALLY4X2: {
    //--------------------

      const int cbpe = env->counted_bits_per_elt();

#pragma omp parallel for reduction(+:num_incorrect) reduction(max:max_incorrect_diff)
      for (size_t index = 0; index < metrics->num_elts_local; ++index) {
        const size_t vi =
          GMMetrics_coord_global_from_index(metrics, index, 0, env);
        const size_t vj =
          GMMetrics_coord_global_from_index(metrics, index, 1, env);
        const size_t vk =
          GMMetrics_coord_global_from_index(metrics, index, 2, env);
        if (vi >= nva || vj >= nva || vk >= nva) {
          continue;
        }
        for (int iE = 0; iE < 2; ++iE) {
          for (int jE = 0; jE < 2; ++jE) {
            for (int kE = 0; kE < 2; ++kE) {
              const GMFloat value = Metrics_ccc_duo_get_3(*metrics,
                index, iE, jE, kE, *env);

              GMTally1 rijk = 0;
              GMTally1 si = 0;
              GMTally1 sj = 0;
              GMTally1 sk = 0;
              GMTally1 ci = 0;
              GMTally1 cj = 0;
              GMTally1 ck = 0;
              GMTally1 cijk = 0;

              for (size_t g=0; g<num_group; ++g) {

                const size_t pf_min = g * group_size_max;
                const size_t pf_max = utils::min((g+1) * group_size_max, nfa);
                const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

                const size_t pvi = perm(g, vi, nva);
                const size_t pvj = perm(g, vj, nva);
                const size_t pvk = perm(g, vk, nva);

                const size_t value_i = value_min + ( pvi * value_max ) /
                                                   (nva+value_min);
                const size_t value_j = value_min + ( pvj * value_max ) /
                                                   (nva+value_min);
                const size_t value_k = value_min + ( pvk * value_max ) /
                                                   (nva+value_min);

                const GMBits2 bval_i = ((size_t)3) & (value_i - value_min);
                const GMBits2 bval_j = ((size_t)3) & (value_j - value_min);
                const GMBits2 bval_k = ((size_t)3) & (value_k - value_min);

                const int bval_i_0 = !!(bval_i&1);
                const int bval_i_1 = !!(bval_i&2);
                const int bval_j_0 = !!(bval_j&1);
                const int bval_j_1 = !!(bval_j&2);
                const int bval_k_0 = !!(bval_k&1);
                const int bval_k_1 = !!(bval_k&2);


                const bool unknown_i = env->sparse() && bval_i == GM_2BIT_UNKNOWN;
                const bool unknown_j = env->sparse() && bval_j == GM_2BIT_UNKNOWN;
                const bool unknown_k = env->sparse() && bval_k == GM_2BIT_UNKNOWN;
                const bool unknown_ijk = unknown_i || unknown_j || unknown_k;

                if (! unknown_i) {
                  ci += gs_this;
                  si += cbpe == 2 ?
                    ((bval_i_0 == iE) + (bval_i_1 == iE)) * gs_this :
                    (bval_i_0 == iE) * gs_this;
                }

                if (! unknown_j) {
                  cj += gs_this;
                  sj += cbpe == 2 ?
                    ((bval_j_0 == jE) + (bval_j_1 == jE)) * gs_this :
                    (bval_j_0 == jE) * gs_this;
                }

                if (! unknown_k) {
                  ck += gs_this;
                  sk += cbpe == 2 ?
                    ((bval_k_0 == kE) + (bval_k_1 == kE)) * gs_this :
                    (bval_k_0 == kE) * gs_this;
                }

                if (! unknown_ijk) {
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

              GMFloat value_expected_floatcalc = 0;
              if (!(ci == 0 || cj == 0 || ck == 0 || cijk == 0)) {
                // FIX typing here
                const double f_one = 1;
  
                const double recip_ci = env->sparse() ? f_one/ci
                                                      : metrics->recip_m;
                const double recip_cj = env->sparse() ? f_one/cj
                                                      : metrics->recip_m;
                const double recip_ck = env->sparse() ? f_one/ck
                                                      : metrics->recip_m;
  
                const double recip_sumcijk = env->sparse() ? f_one/cijk :
                                               (f_one / 8) * metrics->recip_m;
  
                value_expected_floatcalc = cbpe == CBPE::CCC ?
                  Metrics_ccc_duo_value<CBPE::CCC>(*metrics, rijk, si, sj, sk,
                           recip_ci, recip_cj, recip_ck, recip_sumcijk, *env) :
                  Metrics_ccc_duo_value<CBPE::DUO>(*metrics, rijk, si, sj, sk,
                           recip_ci, recip_cj, recip_ck, recip_sumcijk, *env);
              }

              GMFloat value_expected = value_expected_floatcalc;

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
                const double diff = fabs(value - value_expected);
                max_incorrect_diff = diff > max_incorrect_diff ? diff : max_incorrect_diff;
                if (num_incorrect < max_to_print) {
                  fprintf(stderr, "Error: incorrect result detected.  coords %zu %zu %zu  "
                         "expected %.20e  actual %.20e  diff %.20e\n", vi, vj, vk,
                         (double)value_expected, (double)value,
                         (double)value-(double)value_expected);
                }
              }

              num_incorrect += is_incorrect;
            } //---k
          } //---j
        } //---i
      } //---for index
    } break;
    //--------------------
    default:
      COMET_INSIST(false && "Invalid data type.");
  } // switch
  do_->num_incorrect += num_incorrect;
  do_->max_incorrect_diff = max_incorrect_diff > do_->max_incorrect_diff ?
                            max_incorrect_diff : do_->max_incorrect_diff;
}

//=============================================================================

void check_metrics(GMMetrics* metrics, DriverOptions* do_, CEnv* env) {
  COMET_INSIST(metrics && do_ && env);

  if (NULL != do_->input_file_path) {
    return;
  }

  if (GM_PROBLEM_TYPE_ANALYTIC == do_->problem_type) {
    check_metrics_analytic_(metrics, do_, env);
  }
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
