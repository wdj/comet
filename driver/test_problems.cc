//-----------------------------------------------------------------------------
/*!
 * \file   test_problems.cc
 * \author Wayne Joubert
 * \date   Mon Aug  7 17:02:51 EDT 2017
 * \brief  Generator for synthetic test problems.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "stdio.h"
//#include "stdlib.h"
//#include "stddef.h"
//#include "string.h"
//#include "float.h"
//#include "errno.h"

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"

#include "driver.hh"
#include "test_problems.hh"

//=============================================================================
/*---Set the entries of the vectors---*/

void set_vectors_random(GMVectors* vectors, DriverOptions* do_, GMEnv* env) {
  GMInsist(vectors && do_ && env);
  GMInsist(do_->problem_type == GM_PROBLEM_TYPE_RANDOM);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  switch (GMEnv_data_type_vectors(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/
#pragma omp parallel for
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        size_t vector = vl +
            vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
        /*---Fill pad vectors with copies of the last vector---*/
        const size_t vector_capped = vector <= do_->num_vector_active-1 ?
                                     vector : do_->num_vector_active-1;
        int fl = 0;
        for (fl = 0; fl < vectors->num_field_local; ++fl) {
          size_t field = fl +
              vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
          if (field >= vectors->num_field_active) {
            continue;
          }
          /*---Compute element unique id---*/
          const size_t uid = field + vectors->num_field_active*vector_capped;
          /*---Generate large random number---*/
          size_t rand1 = uid;
          rand1 = gm_randomize(rand1);
          rand1 = gm_randomize(rand1);
          size_t rand2 = uid;
          rand2 = gm_randomize(rand2);
          rand2 = gm_randomize(rand2);
          rand2 = gm_randomize(rand2);
          const size_t randomize_max = gm_randomize_max();
          size_t rand_value = rand1 + randomize_max * rand2;
          /*---Reduce so that after summing num_field times the integer
               still exactly representable by floating point type---*/
          const size_t rand_max = randomize_max * randomize_max;
          const int log2_num_summands_3way_numerator = 2;
          const int shift_amount = gm_log2(log2_num_summands_3way_numerator*
                                           rand_max*vectors->num_field_active)
                                   - GMFloat_mant_dig();
          rand_value >>= shift_amount > 0 ? shift_amount : 0;
          /*---Store---*/
          GMFloat float_value = (GMFloat)rand_value;
          GMInsist((size_t)float_value == rand_value);
          GMInsist(float_value * vectors->num_field_active <
                         ((size_t)1)<<GMFloat_mant_dig());
          GMVectors_float_set(vectors, fl, vl, float_value, env);
          /*---Print---*/
          if (do_->verbosity > 2) {
            printf("vec_proc %i vec %i field_proc %i field %i value %e\n",
                   GMEnv_proc_num_vector_i(env), vl,
                   GMEnv_proc_num_field(env), fl, float_value);
          }
        } /*---field_local---*/
      }   /*---vector_local---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
    /*--------------------*/
#pragma omp parallel for
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        size_t vector = vl +
            vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
        /*---Fill pad vectors with copies of the last vector---*/
        const size_t vector_capped = vector <= do_->num_vector_active-1 ?
                                     vector : do_->num_vector_active-1;
        int fl;
        for (fl = 0; fl < vectors->num_field_local; ++fl) {
          size_t field = fl +
              vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
          if (field >= vectors->num_field_active) {
            continue;
          }
          /*---Compute element unique id---*/
          const size_t uid = field + vectors->num_field_active*vector_capped;
          size_t index = uid;
          /*---Randomize---*/
          index = gm_randomize(index);
          index = gm_randomize(index);
          /*---Calculate random number between 0 and 3---*/
          const float float_rand_value = index / (float)gm_randomize_max();
          /*---Create 2-bit value - make extra sure less than 4---*/
          GMBits2 value = (int)((4. - 1e-5) * float_rand_value);
          /*---Store---*/
          GMVectors_bits2_set(vectors, fl, vl, value, env);
          /*---Print---*/
          if (do_->verbosity > 2) {
            printf("vec_proc %i vec %i "
                   "field_proc %i field %i value %.1i%.1i\n",
                   GMEnv_proc_num_vector_i(env), vl,
                   GMEnv_proc_num_field(env), fl, value / 2, value % 2);
          }
        } /*---fl---*/
      }   /*---vl---*/
    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      GMInsist(false && "Invalid data type.");
  } /*---switch---*/
}

//-----------------------------------------------------------------------------

static size_t perm_shuffle(size_t key, size_t i, size_t n) {
  GMAssert((key & (~(size_t)1)) == 0);
  GMAssert(i>=0 && i<n);
  GMAssert(n>=0);

  // For an integer between 0 and n-1, permute it to another such integer.
  // The permutation choice is specified by 1 bit.
  // For an ascending sequence of integers, first output the even values,
  // then the odd values, for key=0. If key=1, same with even/odd reversed.

  const size_t nhalf = (n+1-key)/2;
  const size_t result = i < nhalf ? 2*i + key : 2*(i-nhalf) + 1 - key;
  GMAssert(result>=0 && result<n);
  return result;
}

//-----------------------------------------------------------------------------

enum {NUM_SHUFFLE = 3};

// This is not a totally (pseudo)random permutation.  However it does
// have the advantage that it can be computed formulaically and quickly.

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")

static size_t perm(size_t key, size_t i, size_t n) {
  GMAssert((key & (~(size_t)((1<<NUM_SHUFFLE)-1))) == 0);
  GMAssert(i>=0 && i<n);
  GMAssert(n>=0);

  // For an integer between 0 and n-1, permute it to another such integer.
  // The permutation choice is specified by NUM_SHUFFLE bits.

  size_t result = i;
  size_t key_resid = key;
  for (int shuffle_num = 0; shuffle_num < NUM_SHUFFLE; ++shuffle_num) {
    result = perm_shuffle(key_resid&1, result, n);
    key_resid >>= 1;
  }
  GMAssert(result>=0 && result<n);
  return result;
}

#pragma GCC pop_options

//-----------------------------------------------------------------------------

void set_vectors_analytic(GMVectors* vectors, DriverOptions* do_, GMEnv* env) {
  GMInsist(vectors && do_ && env);
  GMInsist(do_->problem_type == GM_PROBLEM_TYPE_ANALYTIC);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  const size_t nf = vectors->num_field_active;
  const size_t nv = do_->num_vector_active;

  // Upper bound on integer representable exactly by GMFloat.
  const size_t max_float = ((size_t)1) << GMFloat_mant_dig();
  // Sum nf times down the vector, is it still exact.
  const size_t value_limit = (max_float - 1) / nf;

  const size_t value_min = 1;
  const size_t value_max = (nv+value_min) < value_limit ?
                           (nv+value_min) : value_limit;

  // The elements of a single permuted vector are partitioned into
  // "groups", with all elements in a group contiguous and having
  // the same value.
  // By keeping the number of groups (here = 8) much smaller than
  // the vector length, the calculation of the exact comparisons
  // is much cheaper -- the comparison of 2 or 3 vectors by element
  // is the same across all elements of the group.

  const size_t num_group = 1 << NUM_SHUFFLE;
  const size_t group_size_max = (nf+num_group-1) / num_group;

  switch (GMEnv_data_type_vectors(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/
#pragma omp parallel for
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        size_t vector = vl +
            vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
        /*---Fill pad vectors with copies of the last vector---*/
        const size_t vector_capped = vector <= do_->num_vector_active-1 ?
                                     vector : do_->num_vector_active-1;
        int fl = 0;
        for (fl = 0; fl < vectors->num_field_local; ++fl) {
          size_t field = fl +
              vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
          if (field >= vectors->num_field_active) {
            continue;
          }
          const size_t f = field; // field number
          const size_t v = vector_capped; // vector number

          const size_t pf = perm(0, f, nf); // permuted field number
          const size_t g = pf / group_size_max; // group number
          GMAssert(g>=0 && g<num_group);

          const size_t pv = perm(g, v, nv); // permuted vector number

          // Linearly map pv to small interval.
          const size_t value = value_min + ( pv * value_max ) / (nv+value_min);

          const GMFloat float_value = value;

          /*---Store---*/
          GMInsist(float_value * vectors->num_field_active < max_float);
          GMVectors_float_set(vectors, fl, vl, float_value, env);

          /*---Print---*/
          if (do_->verbosity > 2) {
            printf("vec_proc %i vec %i field_proc %i field %i value %e\n",
                   GMEnv_proc_num_vector_i(env), vl,
                   GMEnv_proc_num_field(env), fl, float_value);
          }
        } /*---field_local---*/
      }   /*---vector_local---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
    /*--------------------*/
#pragma omp parallel for
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        size_t vector = vl +
            vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
        /*---Fill pad vectors with copies of the last vector---*/
        const size_t vector_capped = vector <= do_->num_vector_active-1 ?
                                     vector : do_->num_vector_active-1;
        int fl = 0;
        for (fl = 0; fl < vectors->num_field_local; ++fl) {
          size_t field = fl +
              vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
          if (field >= vectors->num_field_active) {
            continue;
          }
          /*---Create 2-bit value - make extra sure less than 4---*/

          const size_t f = field;
          const size_t v = vector_capped;

          const size_t pf = perm(0, f, nf);
          const size_t g = pf / group_size_max;
          GMAssert(g>=0 && g<num_group);

          const size_t pv = perm(g, v, nv);

          const size_t value = value_min + ( pv * value_max ) / (nv+value_min);

          const GMBits2 bval = ((size_t)3) & (value - value_min);

          /*---Store---*/
          GMVectors_bits2_set(vectors, fl, vl, bval, env);

          /*---Print---*/
          if (do_->verbosity > 2) {
            printf("vec_proc %i vec %i "
                   "field_proc %i field %i value %.1i%.1i\n",
                   GMEnv_proc_num_vector_i(env), vl,
                   GMEnv_proc_num_field(env), fl,
                   (int)bval / 2, (int)bval % 2);
          }
        } /*---field_local---*/
      }   /*---vector_local---*/
    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      GMInsist(false && "Invalid data type.");
  } /*---switch---*/
}

//=============================================================================
/*---Check correctness of metrics, if possible---*/

void check_metrics(GMMetrics* metrics, DriverOptions* do_, GMEnv* env) {
  GMInsist(metrics && do_ && env);

  if (GM_PROBLEM_TYPE_ANALYTIC != do_->problem_type ||
      NULL != do_->input_file_path) {
    return;
  }

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  const size_t nf = metrics->num_field_active;
  const size_t nv = metrics->num_vector_active;

  // Upper bound on integer representable exactly by GMFloat.
  const size_t max_float = ((size_t)1) << GMFloat_mant_dig();
  // Sum nf times down the vector, is it still exact.
  const size_t value_limit = (max_float - 1) / nf;

  const size_t value_min = 1;
  const size_t value_max = (nv+value_min) < value_limit ?
                           (nv+value_min) : value_limit;

  const size_t num_group = 1 << NUM_SHUFFLE;
  const size_t group_size_max = (nf+num_group-1) / num_group;

  size_t num_incorrect = 0;

  switch (GMEnv_data_type_metrics(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/
      if (GMEnv_num_way(env) == GM_NUM_WAY_2) {
#pragma omp parallel for reduction(+:num_incorrect)
        for (size_t index = 0; index < metrics->num_elts_local; ++index) {
          const size_t vi =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          const size_t vj =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          if (vi >= metrics->num_vector_active ||
              vj >= metrics->num_vector_active) {
            continue;
          }
          const GMFloat value
            = GMMetrics_czek_get_from_index(metrics, index, env);

          GMFloat n = 0;
          GMFloat d = 0;

          // For each comparison of vectors, the compared/summed
          // elements are treated as num_group groups.  All element
          // comparisons in the group have the same value, so we just
          // compute once and multiply that by the group size.

          for (size_t g=0; g<num_group; ++g) {

            const size_t pf_min = g * group_size_max;
            const size_t pf_max = gm_min_i8((g+1) * group_size_max, nf);
            const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

            const size_t pvi = perm(g, vi, nv);
            const size_t pvj = perm(g, vj, nv);

            const size_t value_i = value_min + ( pvi * value_max ) /
                                               (nv+value_min);
            const size_t value_j = value_min + ( pvj * value_max ) /
                                               (nv+value_min);
            n += gm_min_i8(value_i, value_j) * gs_this;
            d += (value_i + value_j) * gs_this;

          } //---g

          const GMFloat value_expected = (((GMFloat)2) * n) / d;

          num_incorrect += value_expected != value;
        } //---for index
      } //---if
      if (GMEnv_num_way(env) == GM_NUM_WAY_3) {
#pragma omp parallel for reduction(+:num_incorrect)
        for (size_t index = 0; index < metrics->num_elts_local; ++index) {
          const size_t vi =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          const size_t vj =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          const size_t vk =
            GMMetrics_coord_global_from_index(metrics, index, 2, env);
          if (vi >= metrics->num_vector_active ||
              vj >= metrics->num_vector_active ||
              vk >= metrics->num_vector_active) {
            continue;
          }
          const GMFloat value
            = GMMetrics_czek_get_from_index(metrics, index, env);

          GMFloat n = 0;
          GMFloat d = 0;

          for (size_t g=0; g<num_group; ++g) {

            const size_t pf_min = g * group_size_max;
            const size_t pf_max = gm_min_i8((g+1) * group_size_max, nf);
            const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

            const size_t pvi = perm(g, vi, nv);
            const size_t pvj = perm(g, vj, nv);
            const size_t pvk = perm(g, vk, nv);

            const size_t value_i = value_min + ( pvi * value_max ) /
                                               (nv+value_min);
            const size_t value_j = value_min + ( pvj * value_max ) /
                                               (nv+value_min);
            const size_t value_k = value_min + ( pvk * value_max ) /
                                               (nv+value_min);

            n += gm_min_i8(value_i, value_j) * gs_this;
            n += gm_min_i8(value_i, value_k) * gs_this;
            n += gm_min_i8(value_j, value_k) * gs_this;

            n -= gm_min_i8(value_i, gm_min_i8(value_j, value_k)) * gs_this;

            d += (value_i + value_j + value_k) * gs_this;

          } //---g

          const GMFloat value_expected = (((GMFloat)1.5) * n) / d;

          num_incorrect += value_expected != value;
        } //---for index
      } //---if
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY2X2: {
    /*--------------------*/
#pragma omp parallel for reduction(+:num_incorrect)
      for (size_t index = 0; index < metrics->num_elts_local; ++index) {
        const size_t vi =
          GMMetrics_coord_global_from_index(metrics, index, 0, env);
        const size_t vj =
          GMMetrics_coord_global_from_index(metrics, index, 1, env);
        if (vi >= metrics->num_vector_active ||
            vj >= metrics->num_vector_active) {
          continue;
        }
        for (int i0 = 0; i0 < 2; ++i0) {
          for (int i1 = 0; i1 < 2; ++i1) {
            const GMFloat value =
                GMMetrics_ccc_get_from_index_2(metrics, index, i0, i1, env);

            GMTally1 rij = 0;
            GMTally1 si = 0;
            GMTally1 sj = 0;
            GMTally1 ci = 0;
            GMTally1 cj = 0;
            GMTally1 cij = 0;

            for (size_t g=0; g<num_group; ++g) {

              const size_t pf_min = g * group_size_max;
              const size_t pf_max = gm_min_i8((g+1) * group_size_max, nf);
              const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

              const size_t pvi = perm(g, vi, nv);
              const size_t pvj = perm(g, vj, nv);

              const size_t value_i = value_min + ( pvi * value_max ) /
                                                 (nv+value_min);
              const size_t value_j = value_min + ( pvj * value_max ) /
                                                 (nv+value_min);

              const GMBits2 bval_i = ((size_t)3) & (value_i - value_min);
              const GMBits2 bval_j = ((size_t)3) & (value_j - value_min);

              const int bval_i_0 = !!(bval_i&1);
              const int bval_i_1 = !!(bval_i&2);
              const int bval_j_0 = !!(bval_j&1);
              const int bval_j_1 = !!(bval_j&2);

              const bool unknown_i = env->sparse && bval_i == GM_2BIT_UNKNOWN;
              const bool unknown_j = env->sparse && bval_j == GM_2BIT_UNKNOWN;
              const bool unknown_ij = unknown_i || unknown_j;

              if (! unknown_i) {
                ci += gs_this;
                si += ((bval_i_0 == i0) + (bval_i_1 == i0)) * gs_this;
              }

              if (! unknown_j) {
                cj += gs_this;
                sj += ((bval_j_0 == i1) + (bval_j_1 == i1)) * gs_this;
              }

              if (! unknown_ij) {
                cij += 4 * gs_this;
                rij += (((bval_i_0 == i0) && (bval_j_0 == i1)) +
                        ((bval_i_0 == i0) && (bval_j_1 == i1)) +
                        ((bval_i_1 == i0) && (bval_j_0 == i1)) +
                        ((bval_i_1 == i0) && (bval_j_1 == i1))) *
                       gs_this;
              }
            } //---g

            GMFloat value_expected = 0;
            if (!(ci == 0 || cj == 0 || cij == 0)) {
              const GMFloat one = 1;

              const GMFloat recip_ci = env->sparse ? one/ci : metrics->recip_m;
              const GMFloat recip_cj = env->sparse ? one/cj : metrics->recip_m;

              const GMFloat recip_sumcij = env->sparse ? one/cij :
                                             (one / 4) * metrics->recip_m;

              value_expected =
                GMMetrics_ccc_value_2(metrics, rij, si, sj,
                                    recip_ci, recip_cj, recip_sumcij, env);
            }
            num_incorrect += value_expected != value;
          } //---j
        } //---i
      } //---for index
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY4X2: {
    /*--------------------*/
#pragma omp parallel for reduction(+:num_incorrect)
      for (size_t index = 0; index < metrics->num_elts_local; ++index) {
        const size_t vi =
          GMMetrics_coord_global_from_index(metrics, index, 0, env);
        const size_t vj =
          GMMetrics_coord_global_from_index(metrics, index, 1, env);
        const size_t vk =
          GMMetrics_coord_global_from_index(metrics, index, 2, env);
        if (vi >= metrics->num_vector_active ||
            vj >= metrics->num_vector_active ||
            vk >= metrics->num_vector_active) {
          continue;
        }
        for (int i0 = 0; i0 < 2; ++i0) {
          for (int i1 = 0; i1 < 2; ++i1) {
            for (int i2 = 0; i2 < 2; ++i2) {
              const GMFloat value =
               GMMetrics_ccc_get_from_index_3( metrics, index, i0, i1, i2, env);

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
                const size_t pf_max = gm_min_i8((g+1) * group_size_max, nf);
                const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

                const size_t pvi = perm(g, vi, nv);
                const size_t pvj = perm(g, vj, nv);
                const size_t pvk = perm(g, vk, nv);

                const size_t value_i = value_min + ( pvi * value_max ) /
                                                   (nv+value_min);
                const size_t value_j = value_min + ( pvj * value_max ) /
                                                   (nv+value_min);
                const size_t value_k = value_min + ( pvk * value_max ) /
                                                   (nv+value_min);

                const GMBits2 bval_i = ((size_t)3) & (value_i - value_min);
                const GMBits2 bval_j = ((size_t)3) & (value_j - value_min);
                const GMBits2 bval_k = ((size_t)3) & (value_k - value_min);

                const int bval_i_0 = !!(bval_i&1);
                const int bval_i_1 = !!(bval_i&2);
                const int bval_j_0 = !!(bval_j&1);
                const int bval_j_1 = !!(bval_j&2);
                const int bval_k_0 = !!(bval_k&1);
                const int bval_k_1 = !!(bval_k&2);


                const bool unknown_i = env->sparse && bval_i == GM_2BIT_UNKNOWN;
                const bool unknown_j = env->sparse && bval_j == GM_2BIT_UNKNOWN;
                const bool unknown_k = env->sparse && bval_k == GM_2BIT_UNKNOWN;
                const bool unknown_ijk = unknown_i || unknown_j || unknown_k;

                if (! unknown_i) {
                  ci += gs_this;
                  si += ((bval_i_0 == i0) + (bval_i_1 == i0)) * gs_this;
                }

                if (! unknown_j) {
                  cj += gs_this;
                  sj += ((bval_j_0 == i1) + (bval_j_1 == i1)) * gs_this;
                }

                if (! unknown_k) {
                  ck += gs_this;
                  sk += ((bval_k_0 == i2) + (bval_k_1 == i2)) * gs_this;
                }

                if (! unknown_ijk) {
                  cijk += 8 * gs_this;
                  rijk += (((bval_i_0==i0) && (bval_j_0==i1) && (bval_k_0==i2))+
                           ((bval_i_1==i0) && (bval_j_0==i1) && (bval_k_0==i2))+
                           ((bval_i_0==i0) && (bval_j_1==i1) && (bval_k_0==i2))+
                           ((bval_i_1==i0) && (bval_j_1==i1) && (bval_k_0==i2))+
                           ((bval_i_0==i0) && (bval_j_0==i1) && (bval_k_1==i2))+
                           ((bval_i_1==i0) && (bval_j_0==i1) && (bval_k_1==i2))+
                           ((bval_i_0==i0) && (bval_j_1==i1) && (bval_k_1==i2))+
                           ((bval_i_1==i0) && (bval_j_1==i1) && (bval_k_1==i2)))
                          * gs_this;
                }
              } //---g

              GMFloat value_expected = 0;
              if (!(ci == 0 || cj == 0 || ck == 0 || cijk == 0)) {
                const GMFloat one = 1;
  
                const GMFloat recip_ci = env->sparse ? one/ci
                                                     : metrics->recip_m;
                const GMFloat recip_cj = env->sparse ? one/cj
                                                     : metrics->recip_m;
                const GMFloat recip_ck = env->sparse ? one/ck
                                                     : metrics->recip_m;
  
                const GMFloat recip_sumcijk = env->sparse ? one/cijk :
                                               (one / 8) * metrics->recip_m;
  
                value_expected =
                  GMMetrics_ccc_value_3(metrics, rijk, si, sj, sk, recip_ci,
                                        recip_cj, recip_ck, recip_sumcijk, env);
              }
              num_incorrect += value_expected != value;
            } //---k
          } //---j
        } //---i
      } //---for index
    } break;
    /*--------------------*/
    default:
      GMInsist(false && "Invalid data type.");
  } /*---switch---*/
  do_->num_incorrect += num_incorrect;
}

//=============================================================================
