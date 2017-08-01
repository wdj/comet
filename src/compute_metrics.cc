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

//TODOTODO: helper function for case statement
// collapse some cases

/*===========================================================================*/

void gm_compute_metrics(GMMetrics* metrics, GMVectors* vectors, GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);

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
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
      break;

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                            GM_NUM_COMPUTE_METHOD * (2)):
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
      break;

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                            GM_NUM_COMPUTE_METHOD * (2)):
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
      break;

    /*--------------------*/

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_REF +
                            GM_NUM_COMPUTE_METHOD * (3)):
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
      break;

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                            GM_NUM_COMPUTE_METHOD * (3)):
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
      break;

    case GM_METRIC_TYPE_SORENSON +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                            GM_NUM_COMPUTE_METHOD * (3)):
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
      }
    } break;

    case GM_METRIC_TYPE_CCC +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_CPU +
                            GM_NUM_COMPUTE_METHOD * (2)): {
      if (GMEnv_all2all(env)) {
        gm_compute_metrics_2way_all2all(metrics, vectors, env);
      } else {
        gm_compute_metrics_2way_notall2all(metrics, vectors, env);
      }
    } break;

    case GM_METRIC_TYPE_CCC +
        GM_NUM_METRIC_TYPE*(GM_COMPUTE_METHOD_GPU +
                            GM_NUM_COMPUTE_METHOD * (2)): {
      if (GMEnv_all2all(env)) {
        gm_compute_metrics_2way_all2all(metrics, vectors, env);
      } else {
        gm_compute_metrics_2way_notall2all(metrics, vectors, env);
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
      }
    } break;

    /*====================*/

    default:
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } // switch

  // Stop timer.

  double time_end = GMEnv_get_synced_time(env);
  env->time += time_end - time_begin;

  // Check computed element count.

  GMAssertAlways(metrics->num_elts_local == metrics->num_elts_local_computed);

  // Compute global counts of compares and operations.

  double num_elts_local = metrics->num_elts_local;
  double num_elts = 0;

  int mpi_code = 0;
  mpi_code *= 1; // Avoid unused variable warning.

  mpi_code = MPI_Allreduce(&num_elts_local, &num_elts, 1,
                           MPI_DOUBLE, MPI_SUM, GMEnv_mpi_comm_vector(env));
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  env->compares += metrics->num_field*num_elts*metrics->data_type_num_values;

  mpi_code = MPI_Allreduce(&env->ops_local, &env->ops, 1, MPI_DOUBLE, MPI_SUM,
                           GMEnv_mpi_comm_vector(env));
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  // Compute global CPU, GPU memory high water marks.

  const size_t cpu_mem_max_local = env->cpu_mem_max;
  mpi_code = MPI_Allreduce(&cpu_mem_max_local, &env->cpu_mem_max, 1,
                           MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                           GMEnv_mpi_comm_vector(env));
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  const size_t gpu_mem_max_local = env->gpu_mem_max;
  mpi_code = MPI_Allreduce(&gpu_mem_max_local, &env->gpu_mem_max, 1,
                           MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                           GMEnv_mpi_comm_vector(env));
  GMAssertAlways(mpi_code == MPI_SUCCESS);
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
