/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Functions for computing metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h>

#include <stdlib.h>

#include "env.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_metrics.h"
#include "compute_metrics_sorenson.h"
#include "compute_metrics_czekanowski_2way.h"
#include "compute_metrics_czekanowski_3way.h"
#include "compute_metrics_ccc.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

void gm_compute_metrics(GMMetrics* metrics, GMVectors* vectors, GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  switch (env->metric_type +
          GM_NUM_METRIC_TYPE *
              (env->compute_method + GM_NUM_COMPUTE_METHOD * (env->num_way))) {
    /*--------------------*/

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                            GM_NUM_COMPUTE_METHOD * (2)):
      gm_compute_metrics_sorenson_2way_cpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                            GM_NUM_COMPUTE_METHOD * (2)):
      gm_compute_metrics_sorenson_2way_gpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                            GM_NUM_COMPUTE_METHOD * (2)):
      gm_compute_metrics_sorenson_2way_ref(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                            GM_NUM_COMPUTE_METHOD * (3)):
      gm_compute_metrics_sorenson_3way_cpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                            GM_NUM_COMPUTE_METHOD * (3)):
      gm_compute_metrics_sorenson_3way_gpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                            GM_NUM_COMPUTE_METHOD * (3)):
      gm_compute_metrics_sorenson_3way_ref(metrics, vectors, env);
      break;

    /*--------------------*/

    case GM_METRIC_TYPE_CZEKANOWSKI +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                            GM_NUM_COMPUTE_METHOD * (2)):
      gm_compute_metrics_czekanowski_2way_cpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_CZEKANOWSKI +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                            GM_NUM_COMPUTE_METHOD * (2)):
      gm_compute_metrics_czekanowski_2way_gpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_CZEKANOWSKI +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                            GM_NUM_COMPUTE_METHOD * (2)):
      gm_compute_metrics_czekanowski_2way_cpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_CZEKANOWSKI +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                            GM_NUM_COMPUTE_METHOD * (3)):
      gm_compute_metrics_czekanowski_3way_cpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_CZEKANOWSKI +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                            GM_NUM_COMPUTE_METHOD * (3)):
      gm_compute_metrics_czekanowski_3way_gpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_CZEKANOWSKI +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                            GM_NUM_COMPUTE_METHOD * (3)):
      gm_compute_metrics_czekanowski_3way_cpu(metrics, vectors, env);
      break;

    /*--------------------*/

    case GM_METRIC_TYPE_CCC + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                                                  GM_NUM_COMPUTE_METHOD * (2)):
      gm_compute_metrics_ccc_2way_cpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_CCC + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                                                  GM_NUM_COMPUTE_METHOD * (2)):
      gm_compute_metrics_ccc_2way_gpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_CCC + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                                                  GM_NUM_COMPUTE_METHOD * (2)):
      gm_compute_metrics_ccc_2way_cpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_CCC + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                                                  GM_NUM_COMPUTE_METHOD * (3)):
      gm_compute_metrics_ccc_3way_cpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_CCC + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                                                  GM_NUM_COMPUTE_METHOD * (3)):
      gm_compute_metrics_ccc_3way_gpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_CCC + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                                                  GM_NUM_COMPUTE_METHOD * (3)):
      gm_compute_metrics_ccc_3way_cpu(metrics, vectors, env);
      break;

    default:
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---switch---*/
}

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
