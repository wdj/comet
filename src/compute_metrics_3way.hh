//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way.hh
 * \author Wayne Joubert
 * \date   Thu Jan 21 19:07:47 EST 2016
 * \brief  Functions for computing 3-way metrics, header.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_compute_metrics_3way_hh_
#define _gm_compute_metrics_3way_hh_

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"

#ifdef __cplusplus
extern "C" {
#endif

//=============================================================================

void gm_compute_metrics_3way_notall2all(GMMetrics* metrics,
                                        GMVectors* vectors,
                                        GMEnv* env);

void gm_compute_metrics_3way_all2all(GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env);

//=============================================================================

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

//=============================================================================

#endif /*---_gm_compute_metrics_3way_hh_---*/

//-----------------------------------------------------------------------------
