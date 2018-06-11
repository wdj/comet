//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Functions for computing metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "stdio.h"

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "compute_metrics.hh"

//=============================================================================

void GMComputeMetrics_create(
    GMComputeMetrics* this_,
    GMDecompMgr* dm,
    GMEnv* env) {
  GMInsist(this_ && dm && env);

  if (!(GMEnv_num_way(env) == 2 && GMEnv_all2all(env))) {
    return;
  }

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  double time_begin = GMEnv_get_synced_time(env);

  *this_ = {0};

  GMVectorSums_create(&this_->vector_sums_onproc, dm->num_vector_local, env);
  GMVectorSums_create(&this_->vector_sums_offproc, dm->num_vector_local, env);

  for (int i = 0; i < 2; ++i) {
    GMVectors_create_with_buf(&this_->vectors_01[i],
                              GMEnv_data_type_vectors(env), dm, env);
  }

  for (int i = 0; i < 2; ++i) {
    GMMirroredBuf_create(&this_->metrics_buf_01[i],
                         dm->num_vector_local, dm->num_vector_local, env);
  }

  GMMirroredBuf_create(&this_->vectors_buf, dm->num_packedfield_local,
                       dm->num_vector_local, env);

  if (env->do_reduce) {
    GMMirroredBuf_create(&this_->metrics_tmp_buf,
                         dm->num_vector_local, dm->num_vector_local, env);
  }

  double time_end = GMEnv_get_synced_time(env);
  env->time += time_end - time_begin;
}

//-----------------------------------------------------------------------------

void GMComputeMetrics_destroy(
    GMComputeMetrics* this_,
    GMEnv* env) {
  GMInsist(this_ && env);

  if (!(GMEnv_num_way(env) == 2 && GMEnv_all2all(env))) {
    return;
  }

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  double time_begin = GMEnv_get_synced_time(env);

  GMVectorSums_destroy(&this_->vector_sums_onproc, env);
  GMVectorSums_destroy(&this_->vector_sums_offproc, env);

  for (int i = 0; i < 2; ++i) {
    GMVectors_destroy(&this_->vectors_01[i], env);
  }

  for (int i = 0; i < 2; ++i) {
    GMMirroredBuf_destroy(&this_->metrics_buf_01[i], env);
  }

  GMMirroredBuf_destroy(&this_->vectors_buf, env);

  if (env->do_reduce) {
    GMMirroredBuf_destroy(&this_->metrics_tmp_buf, env);
  }

  double time_end = GMEnv_get_synced_time(env);
  env->time += time_end - time_begin;
}

//=============================================================================

void gm_compute_metrics(GMMetrics* metrics, GMVectors* vectors, GMEnv* env) {
  GMInsist(metrics && vectors && env);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  GMComputeMetrics compute_metrics_value = {0},
                  *compute_metrics = &compute_metrics_value;

  GMComputeMetrics_create(compute_metrics, vectors->dm, env);

  gm_compute_metrics(compute_metrics, metrics, vectors, env);

  GMComputeMetrics_destroy(compute_metrics, env);
}

//=============================================================================

void gm_compute_metrics(GMComputeMetrics* compute_metrics, GMMetrics* metrics,
                        GMVectors* vectors, GMEnv* env) {
  GMInsist(metrics && vectors && env);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  // Start timer.

  double time_begin = GMEnv_get_synced_time(env);

  // Perform metrics computation.

  if (GMEnv_num_way(env) == 2 && ! GMEnv_all2all(env)) {

    gm_compute_metrics_2way_notall2all(compute_metrics, metrics, vectors, env);

  } else if (GMEnv_num_way(env) == 2 && GMEnv_all2all(env)) {

    gm_compute_metrics_2way_all2all(compute_metrics, metrics, vectors, env);

  } else if (GMEnv_num_way(env) == 3 && ! GMEnv_all2all(env)) {

    gm_compute_metrics_3way_notall2all(compute_metrics, metrics, vectors, env);

  } else if (GMEnv_num_way(env) == 3 && GMEnv_all2all(env)) {

    gm_compute_metrics_3way_all2all(compute_metrics, metrics, vectors, env);

  } else {

    GMInsistInterface(env, false && "Unimplemented.");

  }

  // Stop timer.

  double time_end = GMEnv_get_synced_time(env);
  env->time += time_end - time_begin;

  // Check computed element count.

  GMInsist(metrics->num_elts_local == metrics->num_elts_local_computed);

  // Compute global counts of compares and operations.

  double num_elts_local = metrics->num_elts_local;
  double num_elts = 0;

  int mpi_code = MPI_Allreduce(&num_elts_local, &num_elts, 1,
                           MPI_DOUBLE, MPI_SUM, GMEnv_mpi_comm_repl_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS);

  env->compares += metrics->num_field * num_elts *
                   metrics->data_type_num_values;

  env->eltcompares += metrics->num_field * num_elts;

  mpi_code = MPI_Allreduce(&env->ops_local, &env->ops, 1, MPI_DOUBLE, MPI_SUM,
                           GMEnv_mpi_comm_repl_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS);

  // Compute global CPU, GPU memory high water marks.

  const size_t cpu_mem_max_local = env->cpu_mem_max;
  mpi_code = MPI_Allreduce(&cpu_mem_max_local, &env->cpu_mem_max, 1,
                           MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                           GMEnv_mpi_comm_repl_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS);

  const size_t gpu_mem_max_local = env->gpu_mem_max;
  mpi_code = MPI_Allreduce(&gpu_mem_max_local, &env->gpu_mem_max, 1,
                           MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                           GMEnv_mpi_comm_repl_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS);
}

//-----------------------------------------------------------------------------
