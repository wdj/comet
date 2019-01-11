//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Functions for computing metrics, headers.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_compute_metrics_hh_
#define _gm_compute_metrics_hh_

#include "env.hh"
#include "vectors.hh"
#include "vector_sums.hh"
#include "metrics.hh"

//=============================================================================

typedef struct {
  GMVectorSums vector_sums_onproc;
  GMVectorSums vector_sums_offproc;
  GMVectors vectors_01[2];
  GMMirroredBuf metrics_buf_01[2];
  GMMirroredBuf vectors_buf;
  GMMirroredBuf metrics_tmp_buf;
} GMComputeMetrics;

//=============================================================================

// TODO: fix this by C++ forward declarqtion of GMComputeMetrics class in these include files.

#include "compute_metrics_2way.hh"
#include "compute_metrics_3way.hh"

//-----------------------------------------------------------------------------

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

//-----------------------------------------------------------------------------

#if 0
void gm_compute_metrics_2way_notall2all(GMComputeMetrics* compute_metrics,
                                        GMMetrics* metrics,
                                        GMVectors* vectors,
                                        GMEnv* env);

void gm_compute_metrics_2way_all2all(GMComputeMetrics* compute_metrics,
                                     GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env);

//-----------------------------------------------------------------------------

void gm_compute_metrics_3way_notall2all(GMComputeMetrics* compute_metrics,
                                        GMMetrics* metrics,
                                        GMVectors* vectors,
                                        GMEnv* env);

void gm_compute_metrics_3way_all2all(GMComputeMetrics* compute_metrics,
                                     GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env);
#endif

//=============================================================================

#endif // _gm_compute_metrics_hh_

//-----------------------------------------------------------------------------
