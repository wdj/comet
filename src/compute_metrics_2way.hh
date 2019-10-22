//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_2way.hh
 * \author Wayne Joubert
 * \date   Thu Jan  7 10:21:09 EST 2016
 * \brief  Calculate metrics, 2-way.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_compute_metrics_2way_hh_
#define _gm_compute_metrics_2way_hh_

#include "env.hh"
#include "decomp_mgr.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

typedef struct {
  GMVectorSums vector_sums_onproc;
  GMVectorSums vector_sums_offproc;
  GMVectors vectors_01[2];
  GMMirroredBuf metrics_buf_01[2];
  GMMirroredBuf vectors_buf;
  GMMirroredBuf metrics_tmp_buf;
} GMComputeMetrics2Way;
    
//=============================================================================

void GMComputeMetrics2Way_create(
  GMComputeMetrics2Way* this_,
  GMDecompMgr* dm,
  GMEnv* env);
  
void GMComputeMetrics2Way_destroy(
  GMComputeMetrics2Way* this_,
  GMEnv* env);

//-----------------------------------------------------------------------------

void gm_compute_metrics_2way_notall2all(
  GMComputeMetrics2Way* this_,
  GMMetrics* metrics,
  GMVectors* vectors,
  GMEnv* env);

void gm_compute_metrics_2way_all2all(
  GMComputeMetrics2Way* this_,
  GMMetrics* metrics,
  GMVectors* vectors,
  GMEnv* env);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _gm_compute_metrics_2way_hh_

//-----------------------------------------------------------------------------
