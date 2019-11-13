//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_2way_block_nums.hh
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Calculate metrics numerators, 2-way, for a single block.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_compute_metrics_2way_block_nums_hh_
#define _gm_compute_metrics_2way_block_nums_hh_

#include "env.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

void gm_compute_2way_proc_nums_start(
  GMVectors* vectors_left,
  GMVectors* vectors_right,
  GMMetrics* metrics,
  GMMirroredBuf* vectors_left_buf,
  GMMirroredBuf* vectors_right_buf,
  GMMirroredBuf* metrics_buf,
  int j_proc,
  bool compute_triang_only,
  GMEnv* env);

void gm_compute_2way_proc_nums_wait(
  GMVectors* vectors_left,
  GMVectors* vectors_right,
  GMMetrics* metrics,
  GMMirroredBuf* vectors_left_buf,
  GMMirroredBuf* vectors_right_buf,
  GMMirroredBuf* metrics_buf,
  int j_proc,
  bool compute_triang_only,
  GMEnv* env);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _gm_compute_metrics_2way_block_nums_hh_

//-----------------------------------------------------------------------------
