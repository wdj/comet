//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_2way.hh
 * \author Wayne Joubert
 * \date   Thu Jan  7 10:21:09 EST 2016
 * \brief  Compute metrics, 2-way, headers.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_compute_metrics_2way_hh_
#define _gm_compute_metrics_2way_hh_

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "compute_metrics.hh"

//=============================================================================

void gm_compute_metrics_2way_notall2all(GMComputeMetrics* compute_metrics,
                                        GMMetrics* metrics,
                                        GMVectors* vectors,
                                        GMEnv* env);

void gm_compute_metrics_2way_all2all(GMComputeMetrics* compute_metrics,
                                     GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env);

//=============================================================================

#endif // _gm_compute_metrics_2way_hh_

//-----------------------------------------------------------------------------
