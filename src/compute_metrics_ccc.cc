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
#include "metrics.h"
#include "compute_metrics_ccc.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

void gm_compute_metrics_ccc_2way_cpu(GMMetrics* metrics,
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



  /*---Compute sums for denominators---*/

//  GMFloat* vector_sums = GMFloat_malloc(metrics->num_vector_local);
//  GMFloat* vector_sums_tmp = GMFloat_malloc(metrics->num_vector_local);

//  gm_compute_float_vector_sums(vectors, vector_sums, vector_sums_tmp, env);

  /*---Compute numerators---*/

  int i = 0;
  int j = 0;
  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      GMTally2x2 sum = GMTally2x2_null();
      int field_local = 0;
      for (field_local = 0; field_local < vectors->num_field_local;
           ++field_local) {
        const GMBits2 value_i = GMVectors_bits2_get(vectors,
                                                    field_local, i, env);
        const GMBits2 value_j = GMVectors_bits2_get(vectors,
                                                    field_local, j, env);

        const int r00 = ( ( ! (value_i & 1) ) && ( ! (value_j & 1) ) ) +
                        ( ( ! (value_i & 1) ) && ( ! (value_j & 2) ) ) +
                        ( ( ! (value_i & 2) ) && ( ! (value_j & 1) ) ) +
                        ( ( ! (value_i & 2) ) && ( ! (value_j & 2) ) );
        const int r01 = ( ( ! (value_i & 1) ) && (   (value_j & 1) ) ) +
                        ( ( ! (value_i & 1) ) && (   (value_j & 2) ) ) +
                        ( ( ! (value_i & 2) ) && (   (value_j & 1) ) ) +
                        ( ( ! (value_i & 2) ) && (   (value_j & 2) ) );
        const int r10 = ( (   (value_i & 1) ) && ( ! (value_j & 1) ) ) +
                        ( (   (value_i & 1) ) && ( ! (value_j & 2) ) ) +
                        ( (   (value_i & 2) ) && ( ! (value_j & 1) ) ) +
                        ( (   (value_i & 2) ) && ( ! (value_j & 2) ) );
        const int r11 = ( (   (value_i & 1) ) && (   (value_j & 1) ) ) +
                        ( (   (value_i & 1) ) && (   (value_j & 2) ) ) +
                        ( (   (value_i & 2) ) && (   (value_j & 1) ) ) +
                        ( (   (value_i & 2) ) && (   (value_j & 2) ) );

// NOTE: "since the sum of all 4 of these relative co-occurences is 1 we really only need to compute 3 of them. Then the last one is just 1 minus the rest.

        sum.data[0] += r00 + (1<<GM_TALLY1_MAX_VALUE_BITS) * r01;
        sum.data[1] += r10 + (1<<GM_TALLY1_MAX_VALUE_BITS) * r11;
      } /*---for k---*/
      GMMetrics_tally2x2_set_2(metrics, i, j, sum, env);
    } /*---for j---*/
  }   /*---for i---*/

  /*---Combine---*/

//  for (i = 0; i < metrics->num_vector_local; ++i) {
//    for (j = i + 1; j < metrics->num_vector_local; ++j) {
//      const GMFloat numerator = GMMetrics_float_get_2(metrics, i, j, env);
//      const GMFloat denominator = vector_sums[i] + vector_sums[j];
//      GMMetrics_float_set_2(metrics, i, j, 2 * numerator / denominator, env);
//    } /*---for j---*/
//  }   /*---for i---*/

  /*---Deallocations---*/

//  free(vector_sums);
//  free(vector_sums_tmp);
}

/*===========================================================================*/

void gm_compute_metrics_ccc_2way_gpu(GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, (!Env_all2all(env)) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

void gm_compute_metrics_ccc_3way_cpu(GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, (!Env_all2all(env)) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

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
