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

  // Check computed element count.

  COMET_INSIST((metrics.num_metric_items_local ==
                metrics.num_metric_items_local_computed ||
                env_.is_shrink()) &&
           "Failure to compute all requested metrics.");

  // Compute counter values: compares.

  size_t num_metrics = 0;

  // NOTE: metrics elts have no field axis so just sum across repl/vector procs.

  COMET_MPI_SAFE_CALL(MPI_Allreduce(&metrics.num_metrics_local, &num_metrics,
    1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, env_.comm_repl_vector()));

  env_.metriccompares_inc(metrics.num_field_active * num_metrics *
                    env_.num_entries_per_metric());
  env_.entrycompares_inc(metrics.num_field_active * num_metrics);
  env_.veccompares_inc(num_metrics);

  // Compute counter values: shrink_achieved.

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
