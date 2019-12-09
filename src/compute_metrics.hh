//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Top-level function to calculate metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_compute_metrics_hh_
#define _comet_compute_metrics_hh_

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "compute_metrics_2way.hh"
#include "compute_metrics_3way.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

class ComputeMetrics {

public:

  ComputeMetrics(GMDecompMgr& dm, Env& env);
  ~ComputeMetrics();

  void compute(GMMetrics& metrics, GMVectors& vectors);

  static void compute(GMMetrics& metrics, GMVectors& vectors,
    Env& env);

private:

  Env& env_;

  ComputeMetrics2Way* compute_metrics_2way_;
  ComputeMetrics3Way* compute_metrics_3way_;

  void compute_stats_(GMMetrics& metrics);

  bool can_run_() const {return env_.is_proc_active();}

  // Convenience class for timing code blocks.

  class CodeBlockTimer {
  public:
    CodeBlockTimer(Env& env) : env_(env), time_begin_(0) {
      time_begin_ = env_.synced_time();
    }
    ~CodeBlockTimer() {
      env_.ctime_inc(env_.synced_time() - time_begin_);
    }
  private:
    Env& env_;
    double time_begin_;
  };

  // Disallowed methods.

  ComputeMetrics(  const ComputeMetrics&);
  void operator=(const ComputeMetrics&);

}; // ComputeMetrics

//=============================================================================

} // namespace comet

#endif // _comet_compute_metrics_hh_

//-----------------------------------------------------------------------------
