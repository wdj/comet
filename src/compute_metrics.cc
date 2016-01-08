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
#include "compute_metrics_2way.cc"
#include "compute_metrics_sorenson.h"
#include "compute_metrics_czekanowski_2way.h"
#include "compute_metrics_czekanowski_3way.h"
#include "compute_metrics_ccc.h"
#include "compute_metrics.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

void gm_compute_metrics(GMMetrics* metrics, GMVectors* vectors, GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  if (!Env_is_proc_active(env)) {
    return;
  }

  double time_begin = GMEnv_get_synced_time(env);

  switch (Env_metric_type(env) + GM_NUM_METRIC_TYPE * (
          Env_compute_method(env) + GM_NUM_COMPUTE_METHOD * (
          Env_num_way(env)))) {
    /*====================*/

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                            GM_NUM_COMPUTE_METHOD * (2)):
      gm_compute_metrics_sorenson_2way_ref(metrics, vectors, env);
      break;

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

    /*--------------------*/

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                            GM_NUM_COMPUTE_METHOD * (3)):
      gm_compute_metrics_sorenson_3way_ref(metrics, vectors, env);
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

    /*====================*/

    case GM_METRIC_TYPE_CZEKANOWSKI + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF+
                            GM_NUM_COMPUTE_METHOD * (2)): {
        if (Env_all2all(env)) {
          gm_compute_metrics_2way_all2all(metrics, vectors, env);
        } else {
          gm_compute_metrics_czekanowski_2way_cpu(metrics, vectors, env);
        }
      } break;

    case GM_METRIC_TYPE_CZEKANOWSKI + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU+
                            GM_NUM_COMPUTE_METHOD * (2)): {
        if (Env_all2all(env)) {
          gm_compute_metrics_2way_all2all(metrics, vectors, env);
        } else {
          gm_compute_metrics_czekanowski_2way_cpu(metrics, vectors, env);
        }
      } break;

    case GM_METRIC_TYPE_CZEKANOWSKI + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU+
                            GM_NUM_COMPUTE_METHOD * (2)): {
        if (Env_all2all(env)) {
          gm_compute_metrics_2way_all2all(metrics, vectors, env);
        } else {
          gm_compute_metrics_czekanowski_2way_gpu(metrics, vectors, env);
        }
      } break;

    /*--------------------*/

    case GM_METRIC_TYPE_CZEKANOWSKI + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF+
                            GM_NUM_COMPUTE_METHOD * (3)):
      gm_compute_metrics_czekanowski_3way_cpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_CZEKANOWSKI + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU+
                            GM_NUM_COMPUTE_METHOD * (3)):
      gm_compute_metrics_czekanowski_3way_cpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_CZEKANOWSKI + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU+
                            GM_NUM_COMPUTE_METHOD * (3)):
      gm_compute_metrics_czekanowski_3way_gpu(metrics, vectors, env);
      break;

    /*====================*/

    case GM_METRIC_TYPE_CCC + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                              GM_NUM_COMPUTE_METHOD * (2)): {
        if (Env_all2all(env)) {
          gm_compute_metrics_2way_all2all(metrics, vectors, env);
        } else {
          gm_compute_metrics_ccc_2way_cpu(metrics, vectors, env);
        }
      } break;

    case GM_METRIC_TYPE_CCC + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                              GM_NUM_COMPUTE_METHOD * (2)): {
        if (Env_all2all(env)) {
          gm_compute_metrics_2way_all2all(metrics, vectors, env);
        } else {
          gm_compute_metrics_ccc_2way_cpu(metrics, vectors, env);
        }
      } break;

    case GM_METRIC_TYPE_CCC + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                              GM_NUM_COMPUTE_METHOD * (2)): {
        if (Env_all2all(env)) {
          gm_compute_metrics_2way_all2all(metrics, vectors, env);
        } else {
          gm_compute_metrics_ccc_2way_gpu(metrics, vectors, env);
        }
      } break;

    /*--------------------*/

    case GM_METRIC_TYPE_CCC + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                                                  GM_NUM_COMPUTE_METHOD * (3)):
      gm_compute_metrics_ccc_3way_cpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_CCC + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                                                  GM_NUM_COMPUTE_METHOD * (3)):
      gm_compute_metrics_ccc_3way_cpu(metrics, vectors, env);
      break;

    case GM_METRIC_TYPE_CCC + GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                                                  GM_NUM_COMPUTE_METHOD * (3)):
      gm_compute_metrics_ccc_3way_gpu(metrics, vectors, env);
      break;

    /*====================*/

    default:
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---switch---*/

  double time_end = GMEnv_get_synced_time(env);

  env->time += time_end - time_begin;

  env->ops += (Env_num_way(env) == 2 && ! Env_all2all(env))
            ?    Env_num_proc_vector(env) * 1. *
                 vectors->num_vector_local * 1. *
                 (vectors->num_vector_local - 1) * (1./2.) *
                 vectors->num_field
            : (Env_num_way(env) == 2 && Env_all2all(env)) 
            ?    vectors->num_vector * 1. *
                 (vectors->num_vector - 1) * (1./2.) *
                 vectors->num_field
            : (Env_num_way(env) == 3 && ! Env_all2all(env)) 
            ?    Env_num_proc_vector(env) * 1. *
                 vectors->num_vector_local * 1. *
                 (vectors->num_vector_local - 1) * 1. *
                 (vectors->num_vector_local - 2) * (1./6.) *
                 vectors->num_field
            :    vectors->num_vector * 1. *
                 (vectors->num_vector - 1) * 1. *
                 (vectors->num_vector - 2) * (1./6.) *
                 vectors->num_field;
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
