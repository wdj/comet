//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Top-level function to calculate metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "compute_metrics_2way.hh"
#include "compute_metrics_3way.hh"
#include "compute_metrics.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Constructor for ComputeMetrics class.

ComputeMetrics::ComputeMetrics(GMDecompMgr& dm, CEnv& env)
  : env_(env)
  , compute_metrics_2way_(NULL)
  , compute_metrics_3way_(NULL) {

  if (! can_run_())
    return;

  CodeBlockTimer timer(env_);

  if (env_.num_way() == NUM_WAY::_2) {
    compute_metrics_2way_ = new ComputeMetrics2Way(dm, env_);
  } else {
    compute_metrics_3way_ = new ComputeMetrics3Way(dm, env);
  }
}

//-----------------------------------------------------------------------------
/// \brief Destructor for ComputeMetrics class.

ComputeMetrics::~ComputeMetrics() {

  if (! can_run_())
    return;

  CodeBlockTimer timer(env_);

  delete compute_metrics_2way_;
  delete compute_metrics_3way_;
}

//-----------------------------------------------------------------------------
/// \brief Function to compute 2-way or 3-way metrics.

void ComputeMetrics::compute(GMMetrics& metrics, GMVectors& vectors) {

  if (! can_run_())
    return;

  {
    CodeBlockTimer timer(env_);

    if (env_.num_way() == 2) {
    compute_metrics_2way_->compute(metrics, vectors);
    } else { // (env_.num_way() == NUM_WAY::_3)
      compute_metrics_3way_->compute(metrics, vectors);
    }

  }

  compute_stats_(metrics);
}

//-----------------------------------------------------------------------------
/// \brief Internal function to calculate misc. statistics for the run.

void ComputeMetrics::compute_stats_(GMMetrics& metrics) {

  // Check computed element count.

  COMET_INSIST(metrics.num_elts_local == metrics.num_elts_local_computed &&
           "Failure to compute all requested metrics.");

  // Compute global counts of compares.

  size_t metrics_num_elts = 0;

  // NOTE: metrics elts have no field axis so just sum across repl/vector procs.

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&metrics.num_elts_local, &metrics_num_elts,
    1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, env_.comm_repl_vector()));

  env_.compares_inc(metrics.num_field_active * metrics_num_elts *
                    metrics.data_type_num_values);
  env_.eltcompares_inc(metrics.num_field_active * metrics_num_elts);
  env_.veccompares_inc(metrics_num_elts);
}

//-----------------------------------------------------------------------------
/// \brief Convenience function for constructor, compute, destructor.

void ComputeMetrics::compute(GMMetrics& metrics, GMVectors& vectors,
  CEnv& env) {

  ComputeMetrics compute_metrics(*vectors.dm, env);
  compute_metrics.compute(metrics, vectors);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
