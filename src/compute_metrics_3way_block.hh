//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way_block.hh
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Calculate numerators, 3-way, for a single block.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_compute_metrics_3way_block_hh_
#define _comet_compute_metrics_3way_block_hh_

#include "env.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Class for computing numerators for 3-way methods.

class ComputeMetrics3WayBlock {

  enum {NUM_BUF = 2};

public:

  ComputeMetrics3WayBlock(int nvl, int npvfl, CEnv& env);
  ~ComputeMetrics3WayBlock();

  void compute(VData vdata_i, VData vdata_j, VData vdata_k,
    GMMetrics& numerators, int j_block, int k_block, int section_step);

private:

  CEnv& env_;

  MirroredBuf tmp_buf_[NUM_BUF];
  MirroredBuf matM_ij_buf_;
  MirroredBuf matM_jk_buf_;
  MirroredBuf matM_kik_buf_;
  MirroredBuf matXitem_buf_[NUM_BUF];
  MirroredBuf matB_buf_[NUM_BUF];

  void compute_linalg_(VData vdata_i, VData vdata_j, VData vdata_k,
    GMMetrics& metrics, int j_block, int k_block, int section_step);

  void compute_czek_(VData vdata_i, VData vdata_j, VData vdata_k,
    GMMetrics& metrics, int j_block, int k_block, int section_step);

  void compute_ccc_duo_(VData vdata_i, VData vdata_j, VData vdata_k,
    GMMetrics& metrics, int j_block, int k_block, int section_step);

  // Disallowed methods.
  ComputeMetrics3WayBlock(const ComputeMetrics3WayBlock&);
  void operator=(const ComputeMetrics3WayBlock&);
};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_compute_metrics_3way_block_hh_

//-----------------------------------------------------------------------------
