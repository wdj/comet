/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_sorenson.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Functions for computing Sorenson metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>

#include "env.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_metrics_sorenson.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

void gm_compute_metrics_sorenson_2way_cpu(GMMetrics* metrics,
                                          GMVectors* vectors,
                                          GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);

  GMInsist(env, (!Env_all2all(env)) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

void gm_compute_metrics_sorenson_2way_gpu(GMMetrics* metrics,
                                          GMVectors* vectors,
                                          GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);

  GMInsist(env, (!Env_all2all(env)) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

void gm_compute_metrics_sorenson_2way_ref(GMMetrics* metrics,
                                          GMVectors* vectors,
                                          GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);

  GMInsist(env, (!Env_all2all(env)) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

void gm_compute_metrics_sorenson_3way_cpu(GMMetrics* metrics,
                                          GMVectors* vectors,
                                          GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);

  GMInsist(env, (!Env_all2all(env)) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

void gm_compute_metrics_sorenson_3way_gpu(GMMetrics* metrics,
                                          GMVectors* vectors,
                                          GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);

  GMInsist(env, (!Env_all2all(env)) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

void gm_compute_metrics_sorenson_3way_ref(GMMetrics* metrics,
                                          GMVectors* vectors,
                                          GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);

  GMInsist(env, (!Env_all2all(env)) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
