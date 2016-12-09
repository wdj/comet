/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Functions for computing metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "compute_metrics_2way.cc"
#include "compute_metrics_3way.cc"
#include "compute_metrics.hh"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

void gm_compute_metrics(GMMetrics* metrics, GMVectors* vectors, GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);

  if (!GM_FP_PRECISION_DOUBLE) {
    GMInsist(env, GMEnv_metric_type(env) != GM_METRIC_TYPE_CCC ? 
    "CCC metric currently not functional under single precision build." : 0);
  }

  if (!GMEnv_is_proc_active(env)) {
    return;
  }

  double time_begin = GMEnv_get_synced_time(env);

  switch (GMEnv_metric_type(env) +
          GM_NUM_METRIC_TYPE * (GMEnv_compute_method(env) +
                                GM_NUM_COMPUTE_METHOD * (GMEnv_num_way(env)))) {
    /*====================*/
    /*---Sorenson---*/
    /*====================*/

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                            GM_NUM_COMPUTE_METHOD * (2)):
      //gm_compute_metrics_sorenson_2way_ref(metrics, vectors, env);
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
      break;

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                            GM_NUM_COMPUTE_METHOD * (2)):
      //gm_compute_metrics_sorenson_2way_cpu(metrics, vectors, env);
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
      break;

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                            GM_NUM_COMPUTE_METHOD * (2)):
      //gm_compute_metrics_sorenson_2way_gpu(metrics, vectors, env);
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
      break;

    /*--------------------*/

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                            GM_NUM_COMPUTE_METHOD * (3)):
      //gm_compute_metrics_sorenson_3way_ref(metrics, vectors, env);
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
      break;

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                            GM_NUM_COMPUTE_METHOD * (3)):
      //gm_compute_metrics_sorenson_3way_cpu(metrics, vectors, env);
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
      break;

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                            GM_NUM_COMPUTE_METHOD * (3)):
      //gm_compute_metrics_sorenson_3way_gpu(metrics, vectors, env);
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
      break;

    /*====================*/
    /*---Czekanowski---*/
    /*====================*/

    case GM_METRIC_TYPE_CZEKANOWSKI +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                            GM_NUM_COMPUTE_METHOD * (2)): {
      if (GMEnv_all2all(env)) {
        gm_compute_metrics_2way_all2all(metrics, vectors, env);
      } else {
        gm_compute_metrics_2way_notall2all(metrics, vectors, env);
      }
    } break;

    case GM_METRIC_TYPE_CZEKANOWSKI +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                            GM_NUM_COMPUTE_METHOD * (2)): {
      if (GMEnv_all2all(env)) {
        gm_compute_metrics_2way_all2all(metrics, vectors, env);
      } else {
        gm_compute_metrics_2way_notall2all(metrics, vectors, env);
      }
    } break;

    case GM_METRIC_TYPE_CZEKANOWSKI +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                            GM_NUM_COMPUTE_METHOD * (2)): {
      if (GMEnv_all2all(env)) {
        gm_compute_metrics_2way_all2all(metrics, vectors, env);
      } else {
        gm_compute_metrics_2way_notall2all(metrics, vectors, env);
      }
    } break;

    /*--------------------*/

    case GM_METRIC_TYPE_CZEKANOWSKI +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                            GM_NUM_COMPUTE_METHOD * (3)): {
      if (GMEnv_all2all(env)) {
        gm_compute_metrics_3way_all2all(metrics, vectors, env);
      } else {
        gm_compute_metrics_3way_notall2all(metrics, vectors, env);
      }
    } break;

    case GM_METRIC_TYPE_CZEKANOWSKI +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                            GM_NUM_COMPUTE_METHOD * (3)): {
      if (GMEnv_all2all(env)) {
        gm_compute_metrics_3way_all2all(metrics, vectors, env);
      } else {
        gm_compute_metrics_3way_notall2all(metrics, vectors, env);
      }
    } break;

    case GM_METRIC_TYPE_CZEKANOWSKI +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                            GM_NUM_COMPUTE_METHOD * (3)): {
      if (GMEnv_all2all(env)) {
        gm_compute_metrics_3way_all2all(metrics, vectors, env);
      } else {
        gm_compute_metrics_3way_notall2all(metrics, vectors, env);
      }
    } break;

    /*====================*/
    /*---CCC---*/
    /*====================*/

    case GM_METRIC_TYPE_CCC +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                            GM_NUM_COMPUTE_METHOD * (2)): {
      if (GMEnv_all2all(env)) {
        gm_compute_metrics_2way_all2all(metrics, vectors, env);
      } else {
        gm_compute_metrics_2way_notall2all(metrics, vectors, env);
        // gm_compute_metrics_ccc_2way_cpu(metrics, vectors, env);
      }
    } break;

    case GM_METRIC_TYPE_CCC +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                            GM_NUM_COMPUTE_METHOD * (2)): {
      if (GMEnv_all2all(env)) {
        gm_compute_metrics_2way_all2all(metrics, vectors, env);
      } else {
        gm_compute_metrics_2way_notall2all(metrics, vectors, env);
        // gm_compute_metrics_ccc_2way_cpu(metrics, vectors, env);
      }
    } break;

    case GM_METRIC_TYPE_CCC +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                            GM_NUM_COMPUTE_METHOD * (2)): {
      if (GMEnv_all2all(env)) {
        gm_compute_metrics_2way_all2all(metrics, vectors, env);
      } else {
        gm_compute_metrics_2way_notall2all(metrics, vectors, env);
        // gm_compute_metrics_ccc_2way_gpu(metrics, vectors, env);
      }
    } break;

    /*--------------------*/

    case GM_METRIC_TYPE_CCC +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                            GM_NUM_COMPUTE_METHOD * (3)): {
      if (GMEnv_all2all(env)) {
        gm_compute_metrics_3way_all2all(metrics, vectors, env);
      } else {
        gm_compute_metrics_3way_notall2all(metrics, vectors, env);
      }
    } break;

    case GM_METRIC_TYPE_CCC +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                            GM_NUM_COMPUTE_METHOD * (3)): {
      if (GMEnv_all2all(env)) {
        gm_compute_metrics_3way_all2all(metrics, vectors, env);
      } else {
        gm_compute_metrics_3way_notall2all(metrics, vectors, env);
      }
    } break;

    case GM_METRIC_TYPE_CCC +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                            GM_NUM_COMPUTE_METHOD * (3)): {
      if (GMEnv_all2all(env)) {
        gm_compute_metrics_3way_all2all(metrics, vectors, env);
      } else {
        gm_compute_metrics_3way_notall2all(metrics, vectors, env);
        // gm_compute_metrics_ccc_3way_gpu(metrics, vectors, env);
      }
    } break;

    /*====================*/

    default:
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---switch---*/

  double time_end = GMEnv_get_synced_time(env);

  env->time += time_end - time_begin;

  /*---Compute a (global) operation count - total work done on all procs---*/

  /* clang-format off */
  env->ops += (GMEnv_num_way(env) == GM_NUM_WAY_2 && ! GMEnv_all2all(env))
            ?    GMEnv_num_proc_vector_i(env) * 1. *
                 vectors->num_vector_local * 1. *
                 (vectors->num_vector_local - 1) * (1./2.) *
                 vectors->num_field
            : (GMEnv_num_way(env) == GM_NUM_WAY_3 && ! GMEnv_all2all(env))
            ?    GMEnv_num_proc_vector_i(env) * 1. *
                 vectors->num_vector_local * 1. *
                 (vectors->num_vector_local - 1) * 1. *
                 (vectors->num_vector_local - 2) * (1./6.) *
                 vectors->num_field
            : (GMEnv_num_way(env) == GM_NUM_WAY_2 && GMEnv_all2all(env))
            ?    vectors->num_vector * 1. *
                 (vectors->num_vector - 1) * (1./2.) *
                 vectors->num_field
            :    vectors->num_vector * 1. *
                 (vectors->num_vector - 1) * 1. *
                 (vectors->num_vector - 2) * (1./6.) *
                 vectors->num_field;
  /* clang-format on */
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
