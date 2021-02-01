//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way_block.hh
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Calculate numerators, 3-way, for a single block.
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

#ifndef _COMET_COMPUTE_METRICS_3WAY_BLOCK_HH_
#define _COMET_COMPUTE_METRICS_3WAY_BLOCK_HH_

#include "env.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"

//=============================================================================

namespace comet {

// Forward declaration.
class MagmaWrapper;

//-----------------------------------------------------------------------------
/// \brief Class for computing numerators for 3-way methods.

class ComputeMetrics3WayBlock {

  enum {NUM_BUF = 2};

public:

  ComputeMetrics3WayBlock(int nvl, int npfl, CEnv& env);
  ~ComputeMetrics3WayBlock();

  void compute(VData vdata_i, VData vdata_j, VData vdata_k,
    GMMetrics& numerators, int j_block, int k_block, int section_step,
    MagmaWrapper& magma_wrapper);

private:

  CEnv& env_;

  MirroredBuf matM_ij_buf_;
  MirroredBuf matM_jk_buf_;
  MirroredBuf matM_kik_buf_;

  MirroredBuf tmp_buf0_;
  MirroredBuf tmp_buf1_;
  MirroredBuf matXitem_buf0_;
  MirroredBuf matXitem_buf1_;
  MirroredBuf matB_buf0_;
  MirroredBuf matB_buf1_;
  MirroredBuf* tmp_buf_[NUM_BUF];
  MirroredBuf* matXitem_buf_[NUM_BUF];
  MirroredBuf* matB_buf_[NUM_BUF];

  void compute_linalg_(VData vdata_i, VData vdata_j, VData vdata_k,
    GMMetrics& metrics, int j_block, int k_block, int section_step,
    MagmaWrapper& magma_wrapper);

  void compute_nonlinalg_(VData vdata_i, VData vdata_j, VData vdata_k,
    GMMetrics& metrics, int j_block, int k_block, int section_step);

  // Disallowed methods.
  ComputeMetrics3WayBlock(const ComputeMetrics3WayBlock&);
  void operator=(const ComputeMetrics3WayBlock&);
};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_COMPUTE_METRICS_3WAY_BLOCK_HH_

//-----------------------------------------------------------------------------
