//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Top-level function to calculate metrics.
 */
//-----------------------------------------------------------------------------
/*-----------------------------------------------------------------------------

Copyright 2020, UT-Battelle, LLC

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-----------------------------------------------------------------------------*/

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

  if (env_.num_way() == NumWay::_2) {
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
    } else { // (env_.num_way() == NumWay::_3)
      compute_metrics_3way_->compute(metrics, vectors);
    }

  }

  compute_stats_(metrics);
}

//-----------------------------------------------------------------------------
/// \brief Internal function to calculate misc. statistics for the run.

void ComputeMetrics::compute_stats_(GMMetrics& metrics) {

  //--------------------
  // Check computed element count.
  //--------------------

  COMET_INSIST((metrics.num_metric_items_local ==
                metrics.num_metric_items_local_computed ||
                env_.is_shrink()) &&
           "Failure to compute all requested metrics.");

  // NOTE: below: metrics elts have no field axis so just sum across
  // repl/vector procs.

  //--------------------
  // Compute counter values: compares.
  //--------------------

  const double num_metrics_local = metrics.num_metrics_local;
  double num_metrics = 0;

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&num_metrics_local,
    &num_metrics, 1, MPI_DOUBLE, MPI_SUM, env_.comm_repl_vector()));

  const double num_metric_compares = num_metrics * metrics.num_field_active;
  const double num_entry_compares = num_metric_compares *
                                    env_.num_entries_per_metric();

  env_.vec_compares_inc(num_metrics);
  env_.metric_compares_inc(num_metric_compares);
  env_.entry_compares_inc(num_entry_compares);

  //--------------------
  // Compute counter values: active compares.
  //--------------------

  const double num_metrics_active_local = metrics.num_metrics_active_local;
  double num_metrics_active = 0;

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&num_metrics_active_local,
    &num_metrics_active, 1, MPI_DOUBLE, MPI_SUM, env_.comm_repl_vector()));

  const double num_metric_active_compares = num_metrics_active *
                                            metrics.num_field_active;
  const double num_entry_active_compares = num_metric_active_compares *
                                    env_.num_entries_per_metric();

  env_.vec_active_compares_inc(num_metrics_active);
  env_.metric_active_compares_inc(num_metric_active_compares);
  env_.entry_active_compares_inc(num_entry_active_compares);

  //--------------------
  // Compute counter values: computed metric entries.
  //--------------------

  const size_t metric_entries_local = metrics.num_metric_items_local *
    env_.num_entries_per_metric_item();
  const size_t metric_entries_local_computed =
    metrics.num_metric_items_local_computed *
    env_.num_entries_per_metric_item();
  size_t metric_entries = 0;
  size_t metric_entries_computed = 0;

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&metric_entries_local, &metric_entries,
    1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, env_.comm_repl_vector()));
  COMET_MPI_SAFE_CALL(MPI_Allreduce(&metric_entries_local_computed,
    &metric_entries_computed,
    1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, env_.comm_repl_vector()));

  env_.metric_entries_inc(metric_entries);
  env_.metric_entries_computed_inc(metric_entries_computed);

  //--------------------
  // Compute counter values: shrink_achieved.
  //--------------------

  const double shrink_achieved_tmp = utils::min(env_.shrink_achieved(),
    metrics.shrink_achieved_local());

  double shrink_achieved = 0;

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&shrink_achieved_tmp, &shrink_achieved,
    1, MPI_DOUBLE, MPI_MIN, env_.comm_repl_vector()));

  env_.shrink_achieved_set(shrink_achieved);
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
