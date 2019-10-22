//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way.hh
 * \author Wayne Joubert
 * \date   Thu Jan 21 19:07:47 EST 2016
 * \brief  Calculate metrics, 3-way.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_compute_metrics_3way_hh_
#define _gm_compute_metrics_3way_hh_

#include "env.hh"
#include "decomp_mgr.hh"
#include "vectors.hh"
#include "metrics.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

typedef struct {
} GMComputeMetrics3Way;
    
//=============================================================================

void GMComputeMetrics3Way_create(
  GMComputeMetrics3Way* this_,
  GMDecompMgr* dm,
  GMEnv* env);

void GMComputeMetrics3Way_destroy(
  GMComputeMetrics3Way* this_,
  GMEnv* env);                
  
//-----------------------------------------------------------------------------

void gm_compute_metrics_3way_notall2all(GMComputeMetrics3Way* this_,
                                        GMMetrics* metrics,
                                        GMVectors* vectors,
                                        GMEnv* env);

void gm_compute_metrics_3way_all2all(GMComputeMetrics3Way* this_,
                                     GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _gm_compute_metrics_3way_hh_

//-----------------------------------------------------------------------------
