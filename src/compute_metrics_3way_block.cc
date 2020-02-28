//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way_block.cc
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Calculate numerators, 3-way, for a single block.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

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

ComputeMetrics3WayBlock::ComputeMetrics3WayBlock(int nvl, int npvfl, CEnv& env)
  : env_(env)
  , tmp_buf_{MirroredBuf(env), MirroredBuf(env)}
  , matM_ij_buf_(env)//{MirroredBuf(env)}
  , matM_jk_buf_(env)//{MirroredBuf(env)}
  , matM_kik_buf_(env)//{MirroredBuf(env)}
  , matXitem_buf_{MirroredBuf(env), MirroredBuf(env)}
  , matB_buf_{MirroredBuf(env), MirroredBuf(env)} {
  COMET_INSIST(nvl >= 0 && npvfl >= 0);

  if (!env_.is_using_linalg())
    return;

  for (int i=0; i<NUM_BUF; ++i) {
    if (env_.do_reduce()) {
      tmp_buf_[i].allocate(nvl, nvl);
    }
    const int matXitem_buf_num_cols = env_.form_matX_on_accel() ? 1 : nvl;
    matXitem_buf_[i].allocate(npvfl, matXitem_buf_num_cols);
    matB_buf_[i].allocate(nvl, nvl);
    matB_buf_[i].set_zero_h(); // ensure determinism
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
  GMMetrics& numerators, int j_block, int k_block, int section_step) {
  COMET_INSIST(j_block >= 0 && j_block < env_.num_block_vector());
  COMET_INSIST(k_block >= 0 && k_block < env_.num_block_vector());
  COMET_INSIST(! (env_.proc_num_vector() == j_block &&
                  env_.proc_num_vector() != k_block));
  COMET_INSIST(! (env_.proc_num_vector() == k_block &&
                  env_.proc_num_vector() != j_block));

  if (env_.is_using_linalg()) {

    compute_linalg_(vdata_i, vdata_j, vdata_k, numerators, j_block, k_block,
      section_step);

  } else if (env_.metric_type() == MetricType::CZEK)  {

    compute_czek_(vdata_i, vdata_j, vdata_k, numerators, j_block, k_block,
      section_step);

  } else if (env_.is_metric_type_bitwise())  {

    compute_ccc_duo_(vdata_i, vdata_j, vdata_k, numerators, j_block, k_block,
      section_step);

  } else {

    COMET_INSIST_INTERFACE(&env_, false &&
      "Selected metric_type unimplemented.");

  }
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
