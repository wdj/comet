//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Top-level function to calculate metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "stdio.h"

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "compute_metrics_2way.hh"
#include "compute_metrics_3way.hh"
#include "compute_metrics.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

/// Constructor for ComputeMetrics class.

ComputeMetrics::ComputeMetrics(GMDecompMgr& dm, GMEnv& env)
  : env_(&env) {

  if (! GMEnv_is_proc_active(env_))
    return;

  double time_begin = GMEnv_get_synced_time(env_);

  if (GMEnv_num_way(env_) == GM_NUM_WAY_2) {
    GMComputeMetrics2Way_create(&compute_metrics_2way_, &dm, env_);
  } else {
    GMComputeMetrics3Way_create(&compute_metrics_3way_, &dm, env_);
  }

  env_->time += GMEnv_get_synced_time(env_) - time_begin;
}

//-----------------------------------------------------------------------------

/// Destructor for ComputeMetrics class.

ComputeMetrics::~ComputeMetrics() {

  if (! GMEnv_is_proc_active(env_))
    return;

  double time_begin = GMEnv_get_synced_time(env_);

  if (GMEnv_num_way(env_) == GM_NUM_WAY_2) {
    GMComputeMetrics2Way_destroy(&compute_metrics_2way_, env_);
  } else {
    GMComputeMetrics3Way_destroy(&compute_metrics_3way_, env_);
  }

  env_->time += GMEnv_get_synced_time(env_) - time_begin;
}

//-----------------------------------------------------------------------------

/// \brief Function to compute 2-way or 3-way metrics.

void ComputeMetrics::compute_metrics(GMMetrics& metrics, GMVectors& vectors) {

  if (! GMEnv_is_proc_active(env_))
    return;

  double time_begin = GMEnv_get_synced_time(env_);

  //--------------------

  if (GMEnv_num_way(env_) == 2 && ! GMEnv_all2all(env_)) {
  
    gm_compute_metrics_2way_notall2all(&compute_metrics_2way_,
      &metrics, &vectors, env_);

  } else if (GMEnv_num_way(env_) == 2 && GMEnv_all2all(env_)) {

    gm_compute_metrics_2way_all2all(&compute_metrics_2way_,
      &metrics, &vectors, env_);

  } else if (GMEnv_num_way(env_) == 3 && ! GMEnv_all2all(env_)) {

    gm_compute_metrics_3way_notall2all(&compute_metrics_3way_,
      &metrics, &vectors, env_);

  } else { // (GMEnv_num_way(env_) == 3 && GMEnv_all2all(env_))

    gm_compute_metrics_3way_all2all(&compute_metrics_3way_,
      &metrics, &vectors, env_);

  }

  //--------------------

  env_->time += GMEnv_get_synced_time(env_) - time_begin;

  compute_stats_(metrics);
}

//-----------------------------------------------------------------------------

/// \brief Internal function to calculate misc. statistics for the run.

void ComputeMetrics::compute_stats_(GMMetrics& metrics) {

  // Check computed element count.

  GMInsist(metrics.num_elts_local == metrics.num_elts_local_computed &&
           "Failure to compute all requested metrics.");

  // Compute global counts of compares and operations.

  size_t metrics_num_elts = 0;

  // NOTE: metrics elts have no field axis so just sum across repl/vector procs.

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&metrics.num_elts_local, &metrics_num_elts,
    1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, GMEnv_mpi_comm_repl_vector(env_)));

  env_->eltcompares += metrics.num_field_active * metrics_num_elts;
  env_->compares += env_->eltcompares * metrics.data_type_num_values;
  env_->veccompares += num_elts;

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&env_->ops_local, &env_->ops, 1,
    MPI_DOUBLE, MPI_SUM, GMEnv_mpi_comm_repl_vector(env_)));

  // Compute global CPU, GPU memory high water marks.

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&env_->cpu_mem_max_local,
    &env_->cpu_mem_max, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX,
    GMEnv_mpi_comm_repl_vector(env_)));

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&env_->gpu_mem_max_local,
    &env_->gpu_mem_max, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX,
    GMEnv_mpi_comm_repl_vector(env_)));
}

//-----------------------------------------------------------------------------

/// \brief Convenience function for single call to compute the metrics.

void ComputeMetrics::compute_metrics(GMMetrics& metrics, GMVectors& vectors,
  GMEnv& env) {

  ComputeMetrics compute_metrics(*vectors.dm, env);
  compute_metrics.compute_metrics(metrics, vectors);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
