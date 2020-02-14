//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_2way_block_combine.hh
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Combine numerators and denominators, 2-way, for a single block.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_compute_metrics_2way_block_combine_hh_
#define _comet_compute_metrics_2way_block_combine_hh_

#include "env.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

void gm_compute_2way_proc_combine(
  GMMetrics* metrics,
  MirroredBuf* metrics_buf,
  const VectorSums* const vector_sums_left,
  const VectorSums* const vector_sums_right,
  int j_proc,
  bool do_compute_triang_only,
  GMEnv* env);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_compute_metrics_2way_block_combine_hh_

//-----------------------------------------------------------------------------
