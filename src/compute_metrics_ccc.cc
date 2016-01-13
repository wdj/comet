/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_ccc.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Functions for computing CCC metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>

#include "env.h"
#include "vectors.h"
#include "vector_sums.h"
#include "metrics.h"
#include "compute_metrics_2way.h"
#include "compute_metrics_ccc.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

#if 0
void gm_compute_metrics_ccc_2way_cpu(GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, Env_num_proc_field(env) == 1
    ? "num_proc_field>1 for CPU case not supported" : 0);

  if (Env_all2all(env)) {
    gm_compute_metrics_2way_all2all(metrics, vectors, env);
    return;
  }

  /*---Compute sums---*/

  GMFloat* vector_sums = GMFloat_malloc(metrics->num_vector_local);
  GMFloat* vector_sums_tmp = GMFloat_malloc(metrics->num_vector_local);

  gm_compute_bits2_vector_sums(vectors, vector_sums, vector_sums_tmp, env);

  /*---Compute R (up to scaling)---*/

  int i = 0;
  int j = 0;
  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      GMTally2x2 sum = GMTally2x2_null();
      int f = 0;
      for (f = 0; f < vectors->num_field_local; ++f) {
        const GMBits2 value_i = GMVectors_bits2_get(vectors, f, i, env);
        const GMBits2 value_j = GMVectors_bits2_get(vectors, f, j, env);

        /* clang-format off */
        const int r00 = ( ( !(value_i & 1) ) && ( !(value_j & 1) ) ) +
                        ( ( !(value_i & 1) ) && ( !(value_j & 2) ) ) +
                        ( ( !(value_i & 2) ) && ( !(value_j & 1) ) ) +
                        ( ( !(value_i & 2) ) && ( !(value_j & 2) ) );
        const int r01 = ( ( !(value_i & 1) ) && (  (value_j & 1) ) ) +
                        ( ( !(value_i & 1) ) && (  (value_j & 2) ) ) +
                        ( ( !(value_i & 2) ) && (  (value_j & 1) ) ) +
                        ( ( !(value_i & 2) ) && (  (value_j & 2) ) );
        const int r10 = ( (  (value_i & 1) ) && ( !(value_j & 1) ) ) +
                        ( (  (value_i & 1) ) && ( !(value_j & 2) ) ) +
                        ( (  (value_i & 2) ) && ( !(value_j & 1) ) ) +
                        ( (  (value_i & 2) ) && ( !(value_j & 2) ) );
        const int r11 = ( (  (value_i & 1) ) && (  (value_j & 1) ) ) +
                        ( (  (value_i & 1) ) && (  (value_j & 2) ) ) +
                        ( (  (value_i & 2) ) && (  (value_j & 1) ) ) +
                        ( (  (value_i & 2) ) && (  (value_j & 2) ) );
        /* clang-format on */

// NOTE: "since the sum of all 4 of these relative co-occurences is 1 we really only need to compute 3 of them. Then the last one is just 1 minus the rest."

        sum.data[0] += GMTally1_encode(r00, r01);
        sum.data[1] += GMTally1_encode(r10, r11);
      } /*---for field---*/
      GMMetrics_tally2x2_set_2(metrics, i, j, sum, env);
    } /*---for j---*/
  }   /*---for i---*/

  /*---Compute multipliers---*/
 
  //const GMFloat one = 1;
  //const GMFloat m = vectors->num_field;
  //const GMFloat recip_m = one / m;

  for (i = 0; i < metrics->num_vector_local; ++i) {
    const GMTally1 si_1 = (GMTally1)(vector_sums[i]);
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      const GMTally1 sj_1 = (GMTally1)(vector_sums[j]);
      const GMFloat si1_sj1 = GMTally1_encode(si_1, sj_1);
      GMMetrics_float_M_set_2(metrics, i, j, si1_sj1, env);
    } /*---for j---*/
  }   /*---for i---*/

  /*---Deallocations---*/

  free(vector_sums);
  free(vector_sums_tmp);
}
#endif

/*===========================================================================*/

void gm_compute_metrics_ccc_3way_cpu(GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, Env_num_proc_field(env) == 1
    ? "num_proc_field>1 for CPU case not supported" : 0);

  if (Env_all2all(env)) {
//    gm_compute_metrics_czekanowski_2way_all2all(metrics, vectors, env);
    GMInsist(env, (!Env_all2all(env)) ? "Unimplemented." : 0);
    return;
  }

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*---Compute sums---*/

  GMFloat* vector_sums = GMFloat_malloc(metrics->num_vector_local);
  GMFloat* vector_sums_tmp = GMFloat_malloc(metrics->num_vector_local);

  gm_compute_bits2_vector_sums(vectors, vector_sums, vector_sums_tmp, env);

  /*---Compute R (up to scaling)---*/

//---NOTE: The following is presently a GUESS as to what the algorithm is.

  int i = 0;
  int j = 0;
  int k = 0;
  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      for (k = j + 1; k < metrics->num_vector_local; ++k) {
        GMTally4x2 sum = GMTally4x2_null();
        int field_local = 0;
        for (field_local = 0; field_local < vectors->num_field_local;
             ++field_local) {
          const GMBits2 value_i = GMVectors_bits2_get(vectors,
                                                      field_local, i, env);
          const GMBits2 value_j = GMVectors_bits2_get(vectors,
                                                      field_local, j, env);
          const GMBits2 value_k = GMVectors_bits2_get(vectors,
                                                      field_local, k, env);
          /* clang-format off */
          const int r000 =
            ( ( !(value_i & 1) ) && ( !(value_j & 1) ) && ( !(value_k & 1) ) ) +
            ( ( !(value_i & 1) ) && ( !(value_j & 1) ) && ( !(value_k & 2) ) ) +
            ( ( !(value_i & 1) ) && ( !(value_j & 2) ) && ( !(value_k & 1) ) ) +
            ( ( !(value_i & 1) ) && ( !(value_j & 2) ) && ( !(value_k & 2) ) ) +
            ( ( !(value_i & 2) ) && ( !(value_j & 1) ) && ( !(value_k & 1) ) ) +
            ( ( !(value_i & 2) ) && ( !(value_j & 1) ) && ( !(value_k & 2) ) ) +
            ( ( !(value_i & 2) ) && ( !(value_j & 2) ) && ( !(value_k & 1) ) ) +
            ( ( !(value_i & 2) ) && ( !(value_j & 2) ) && ( !(value_k & 2) ) );
          const int r001 =
            ( ( !(value_i & 1) ) && ( !(value_j & 1) ) && (  (value_k & 1) ) ) +
            ( ( !(value_i & 1) ) && ( !(value_j & 1) ) && (  (value_k & 2) ) ) +
            ( ( !(value_i & 1) ) && ( !(value_j & 2) ) && (  (value_k & 1) ) ) +
            ( ( !(value_i & 1) ) && ( !(value_j & 2) ) && (  (value_k & 2) ) ) +
            ( ( !(value_i & 2) ) && ( !(value_j & 1) ) && (  (value_k & 1) ) ) +
            ( ( !(value_i & 2) ) && ( !(value_j & 1) ) && (  (value_k & 2) ) ) +
            ( ( !(value_i & 2) ) && ( !(value_j & 2) ) && (  (value_k & 1) ) ) +
            ( ( !(value_i & 2) ) && ( !(value_j & 2) ) && (  (value_k & 2) ) );
          const int r010 =
            ( ( !(value_i & 1) ) && (  (value_j & 1) ) && ( !(value_k & 1) ) ) +
            ( ( !(value_i & 1) ) && (  (value_j & 1) ) && ( !(value_k & 2) ) ) +
            ( ( !(value_i & 1) ) && (  (value_j & 2) ) && ( !(value_k & 1) ) ) +
            ( ( !(value_i & 1) ) && (  (value_j & 2) ) && ( !(value_k & 2) ) ) +
            ( ( !(value_i & 2) ) && (  (value_j & 1) ) && ( !(value_k & 1) ) ) +
            ( ( !(value_i & 2) ) && (  (value_j & 1) ) && ( !(value_k & 2) ) ) +
            ( ( !(value_i & 2) ) && (  (value_j & 2) ) && ( !(value_k & 1) ) ) +
            ( ( !(value_i & 2) ) && (  (value_j & 2) ) && ( !(value_k & 2) ) );
          const int r011 =
            ( ( !(value_i & 1) ) && (  (value_j & 1) ) && (  (value_k & 1) ) ) +
            ( ( !(value_i & 1) ) && (  (value_j & 1) ) && (  (value_k & 2) ) ) +
            ( ( !(value_i & 1) ) && (  (value_j & 2) ) && (  (value_k & 1) ) ) +
            ( ( !(value_i & 1) ) && (  (value_j & 2) ) && (  (value_k & 2) ) ) +
            ( ( !(value_i & 2) ) && (  (value_j & 1) ) && (  (value_k & 1) ) ) +
            ( ( !(value_i & 2) ) && (  (value_j & 1) ) && (  (value_k & 2) ) ) +
            ( ( !(value_i & 2) ) && (  (value_j & 2) ) && (  (value_k & 1) ) ) +
            ( ( !(value_i & 2) ) && (  (value_j & 2) ) && (  (value_k & 2) ) );
          const int r100 =
            ( (  (value_i & 1) ) && ( !(value_j & 1) ) && ( !(value_k & 1) ) ) +
            ( (  (value_i & 1) ) && ( !(value_j & 1) ) && ( !(value_k & 2) ) ) +
            ( (  (value_i & 1) ) && ( !(value_j & 2) ) && ( !(value_k & 1) ) ) +
            ( (  (value_i & 1) ) && ( !(value_j & 2) ) && ( !(value_k & 2) ) ) +
            ( (  (value_i & 2) ) && ( !(value_j & 1) ) && ( !(value_k & 1) ) ) +
            ( (  (value_i & 2) ) && ( !(value_j & 1) ) && ( !(value_k & 2) ) ) +
            ( (  (value_i & 2) ) && ( !(value_j & 2) ) && ( !(value_k & 1) ) ) +
            ( (  (value_i & 2) ) && ( !(value_j & 2) ) && ( !(value_k & 2) ) );
          const int r101 =
            ( (  (value_i & 1) ) && ( !(value_j & 1) ) && (  (value_k & 1) ) ) +
            ( (  (value_i & 1) ) && ( !(value_j & 1) ) && (  (value_k & 2) ) ) +
            ( (  (value_i & 1) ) && ( !(value_j & 2) ) && (  (value_k & 1) ) ) +
            ( (  (value_i & 1) ) && ( !(value_j & 2) ) && (  (value_k & 2) ) ) +
            ( (  (value_i & 2) ) && ( !(value_j & 1) ) && (  (value_k & 1) ) ) +
            ( (  (value_i & 2) ) && ( !(value_j & 1) ) && (  (value_k & 2) ) ) +
            ( (  (value_i & 2) ) && ( !(value_j & 2) ) && (  (value_k & 1) ) ) +
            ( (  (value_i & 2) ) && ( !(value_j & 2) ) && (  (value_k & 2) ) );
          const int r110 =
            ( (  (value_i & 1) ) && (  (value_j & 1) ) && ( !(value_k & 1) ) ) +
            ( (  (value_i & 1) ) && (  (value_j & 1) ) && ( !(value_k & 2) ) ) +
            ( (  (value_i & 1) ) && (  (value_j & 2) ) && ( !(value_k & 1) ) ) +
            ( (  (value_i & 1) ) && (  (value_j & 2) ) && ( !(value_k & 2) ) ) +
            ( (  (value_i & 2) ) && (  (value_j & 1) ) && ( !(value_k & 1) ) ) +
            ( (  (value_i & 2) ) && (  (value_j & 1) ) && ( !(value_k & 2) ) ) +
            ( (  (value_i & 2) ) && (  (value_j & 2) ) && ( !(value_k & 1) ) ) +
            ( (  (value_i & 2) ) && (  (value_j & 2) ) && ( !(value_k & 2) ) );
          const int r111 =
            ( (  (value_i & 1) ) && (  (value_j & 1) ) && (  (value_k & 1) ) ) +
            ( (  (value_i & 1) ) && (  (value_j & 1) ) && (  (value_k & 2) ) ) +
            ( (  (value_i & 1) ) && (  (value_j & 2) ) && (  (value_k & 1) ) ) +
            ( (  (value_i & 1) ) && (  (value_j & 2) ) && (  (value_k & 2) ) ) +
            ( (  (value_i & 2) ) && (  (value_j & 1) ) && (  (value_k & 1) ) ) +
            ( (  (value_i & 2) ) && (  (value_j & 1) ) && (  (value_k & 2) ) ) +
            ( (  (value_i & 2) ) && (  (value_j & 2) ) && (  (value_k & 1) ) ) +
            ( (  (value_i & 2) ) && (  (value_j & 2) ) && (  (value_k & 2) ) );
          /* clang-format on */

//---NOTE: must check that order is correct.

          sum.data[0] += GMTally1_encode(r000, r001);
          sum.data[1] += GMTally1_encode(r010, r011);
          sum.data[2] += GMTally1_encode(r100, r101);
          sum.data[3] += GMTally1_encode(r110, r111);
        } /*---for field---*/
        GMMetrics_tally4x2_set_3(metrics, i, j, k, sum, env);
      } /*---for k---*/
    } /*---for j---*/
  }   /*---for i---*/

  /*---Compute multipliers---*/
 
  const GMFloat one = 1;

  const GMFloat m = vectors->num_field;
  const GMFloat recip_m = one / m;

//---NOTE: fix this
  for (i = 0; i < metrics->num_vector_local; ++i) {
    //const GMFloat f_i = ((one / 2 ) * recip_m) * vector_sums[i];
    const GMFloat s_i = vector_sums[i];
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      //const GMFloat f_j = ((one / 2 ) * recip_m) * vector_sums[j];
      const GMFloat s_j = vector_sums[j];
       for (k = j + 1; j < metrics->num_vector_local; ++j) {
         //const GMFloat f_k = ((one / 2 ) * recip_m) * vector_sums[k];
         const GMFloat s_k = vector_sums[k];
        //const GMFloat result = ((one / 2 ) * recip_m) *
        //                       (one - (2*one/3) * f_i) *
        //                       (one - (2*one/3) * f_j) *
        //                       (one - (2*one/3) * f_k);
        /*---Arrange so as to guarantee each factor nonnegative---*/
        const GMFloat result = ((one / 2 ) * recip_m) *
                               (3 * m - s_i) * (one/3) * recip_m *
                               (3 * m - s_j) * (one/3) * recip_m *
                               (3 * m - s_k) * (one/3) * recip_m;
        GMMetrics_float_M_set_3(metrics, i, j, k, result, env);
      } /*---for k---*/
    } /*---for j---*/
  }   /*---for i---*/

  /*---Deallocations---*/

  free(vector_sums);
  free(vector_sums_tmp);
}

/*===========================================================================*/

#if 0
void gm_compute_metrics_ccc_2way_gpu(GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, (!Env_all2all(env)) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}
#endif

/*===========================================================================*/

void gm_compute_metrics_ccc_3way_gpu(GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, (!Env_all2all(env)) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
