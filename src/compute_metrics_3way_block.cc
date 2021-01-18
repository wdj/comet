//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way_block.cc
 * \author Wayne Joubert
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

#include "string.h"

#include "env.hh"
#include "linalg.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "comm_xfer_utils.hh"
#include "compute_metrics_3way_block.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

ComputeMetrics3WayBlock::ComputeMetrics3WayBlock(int nvl, int npfl, CEnv& env)
  : env_(env)
  , matM_ij_buf_(env)
  , matM_jk_buf_(env)
  , matM_kik_buf_(env)
  , tmp_buf0_(env)
  , tmp_buf1_(env)
  , matXitem_buf0_(env)
  , matXitem_buf1_(env)
  , matB_buf0_(env)
  , matB_buf1_(env)
  , tmp_buf_{&tmp_buf0_, &tmp_buf1_}
  , matXitem_buf_{&matXitem_buf0_, &matXitem_buf1_}
  , matB_buf_{&matB_buf0_, &matB_buf1_}

 {
  COMET_INSIST(nvl >= 0 && npfl >= 0);

  if (!env_.is_using_linalg())
    return;

  for (int i=0; i<NUM_BUF; ++i) {
    if (env_.do_reduce())
      tmp_buf_[i]->allocate(nvl, nvl);
    const int matXitem_buf_num_cols = env_.form_matX_tc() ? 1 : nvl;
    matXitem_buf_[i]->allocate(npfl, matXitem_buf_num_cols);
    matB_buf_[i]->allocate(nvl, nvl);
    // ensure determinism
    matXitem_buf_[i]->set_zero_h();
    matXitem_buf_[i]->to_accel();
    matB_buf_[i]->set_zero_h();
    matB_buf_[i]->to_accel();
  }
  if (env_.does_3way_need_2way()) {
    matM_ij_buf_.allocate(nvl, nvl);
    matM_jk_buf_.allocate(nvl, nvl);
    matM_kik_buf_.allocate(nvl, nvl);
  }
}

//-----------------------------------------------------------------------------

ComputeMetrics3WayBlock::~ComputeMetrics3WayBlock() {

  if (!env_.is_using_linalg())
    return;

}

//-----------------------------------------------------------------------------

void ComputeMetrics3WayBlock::compute(
  VData vdata_i, VData vdata_j, VData vdata_k, 
  GMMetrics& numerators, int j_block, int k_block, int section_step,
  MagmaWrapper& magma_wrapper) {
  COMET_INSIST(j_block >= 0 && j_block < env_.num_block_vector());
  COMET_INSIST(k_block >= 0 && k_block < env_.num_block_vector());
  COMET_INSIST(! (env_.proc_num_vector() == j_block &&
                  env_.proc_num_vector() != k_block));
  COMET_INSIST(! (env_.proc_num_vector() == k_block &&
                  env_.proc_num_vector() != j_block));

  if (env_.is_using_linalg()) {

    compute_linalg_(vdata_i, vdata_j, vdata_k, numerators, j_block, k_block,
      section_step, magma_wrapper);

  } else {

    compute_nonlinalg_(vdata_i, vdata_j, vdata_k, numerators, j_block, k_block,
      section_step);
  }
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
