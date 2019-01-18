//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Top-level function to calculate metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_compute_metrics_hh_
#define _gm_compute_metrics_hh_

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "compute_metrics_2way.hh"
#include "compute_metrics_3way.hh"

//=============================================================================

typedef struct {
  GMComputeMetrics2Way compute_metrics_2way;
  GMComputeMetrics3Way compute_metrics_3way;
} GMComputeMetrics;

//=============================================================================

void GMComputeMetrics_create(
    GMComputeMetrics* this_,
    GMDecompMgr* dm,
    GMEnv* env);

void GMComputeMetrics_destroy(
    GMComputeMetrics* this_,
    GMEnv* env);

//-----------------------------------------------------------------------------

void gm_compute_metrics(GMComputeMetrics* compute_metrics, GMMetrics* metrics,
                        GMVectors* vectors, GMEnv* env);

void gm_compute_metrics(GMMetrics* metrics, GMVectors* vectors, GMEnv* env);

//=============================================================================

#endif // _gm_compute_metrics_hh_

//-----------------------------------------------------------------------------
