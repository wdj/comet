//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics.hh
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

#ifndef _COMET_COMPUTE_METRICS_HH_
#define _COMET_COMPUTE_METRICS_HH_

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

  ComputeMetrics(GMDecompMgr& dm, CEnv& env);
  ~ComputeMetrics();
  void terminate();

  void compute(GMMetrics& metrics, GMVectors& vectors);

  static void compute(GMMetrics& metrics, GMVectors& vectors,
    CEnv& env);

private:

  CEnv& env_;
  bool is_active_;

  ComputeMetrics2Way* compute_metrics_2way_;
  ComputeMetrics3Way* compute_metrics_3way_;

  void compute_stats_(GMMetrics& metrics);

  bool can_run_() const {return env_.is_proc_active();}

  // Convenience class for timing code blocks.

  class CodeBlockTimer {
  public:
    CodeBlockTimer(CEnv& env) : env_(env), time_begin_(0) {
      time_begin_ = env_.synced_time();
    }
    ~CodeBlockTimer() {
      env_.ctime_inc(env_.synced_time() - time_begin_);
    }
  private:
    CEnv& env_;
    double time_begin_;
  };

  // Disallowed methods.

  ComputeMetrics(  const ComputeMetrics&);
  void operator=(const ComputeMetrics&);

}; // ComputeMetrics

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_COMPUTE_METRICS_HH_

//-----------------------------------------------------------------------------
