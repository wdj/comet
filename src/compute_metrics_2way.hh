//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_2way.hh
 * \author Wayne Joubert
 * \date   Thu Jan  7 10:21:09 EST 2016
 * \brief  Calculate metrics, 2-way.
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

#ifndef _COMET_COMPUTE_METRICS_2WAY_HH_
#define _COMET_COMPUTE_METRICS_2WAY_HH_

#include "env.hh"
#include "decomp_mgr.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

class ComputeMetrics2Way {

  enum {NUM_BUF = 2};

public:

  ComputeMetrics2Way(GMDecompMgr& dm, CEnv& env);
  ~ComputeMetrics2Way();

  void compute(GMMetrics& metrics, Vectors& vectors);

private:

  CEnv& env_;

  Vectors vectors_0_;
  Vectors vectors_1_;
  Vectors* vectors_01_[NUM_BUF];
  Vectors vectors_left_alt_;
  MirroredBuf metrics_buf_0_;
  MirroredBuf metrics_buf_1_;
  MirroredBuf* metrics_buf_01_[NUM_BUF];
  MirroredBuf metrics_tmp_buf_;
  VectorSums vector_sums_onproc_;
  VectorSums vector_sums_offproc_0_;
  VectorSums vector_sums_offproc_1_;
  VectorSums* vector_sums_offproc_01_[NUM_BUF];

  void compute_notall2all_(GMMetrics& metrics, Vectors& vectors);
  void compute_all2all_(GMMetrics& metrics, Vectors& vectors);

  void lock(bool& lock_val) {
    COMET_INSIST(! lock_val);
    lock_val = true;
  };

  void unlock(bool& lock_val) {
    COMET_INSIST(lock_val);
    lock_val = false;
  };

  // Disallowed methods.

  ComputeMetrics2Way(const ComputeMetrics2Way&);
  void operator=(const ComputeMetrics2Way&);

};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_COMPUTE_METRICS_2WAY_HH_

//-----------------------------------------------------------------------------
