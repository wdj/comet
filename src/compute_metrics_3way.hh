//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way.hh
 * \author Wayne Joubert
 * \date   Thu Jan 21 19:07:47 EST 2016
 * \brief  Calculate metrics, 3-way.
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

#ifndef _COMET_COMPUTE_METRICS_3WAY_HH_
#define _COMET_COMPUTE_METRICS_3WAY_HH_

#include "env.hh"
#include "decomp_mgr.hh"
#include "vectors.hh"
#include "metrics.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

class ComputeMetrics3Way {

public:

  ComputeMetrics3Way(GMDecompMgr& dm, CEnv& env);
  ~ComputeMetrics3Way();

  void compute(GMMetrics& metrics, GMVectors& vectors);

private:

  CEnv& env_;

  void compute_notall2all_(GMMetrics& metrics, GMVectors& vectors);
  void compute_all2all_(GMMetrics& metrics, GMVectors& vectors);

  // Disallowed methods.

  ComputeMetrics3Way(const ComputeMetrics3Way&);
  void operator=(const ComputeMetrics3Way&);

}; // ComputeMetrics3Way

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_COMPUTE_METRICS_3WAY_HH_

//-----------------------------------------------------------------------------
