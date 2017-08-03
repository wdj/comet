/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_utils_2way.hh
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing metrics, utilities, for 2-way case, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _gm_compute_metrics_utils_2way_hh_
#define _gm_compute_metrics_utils_2way_hh_

#include "env.hh"
#include "mirrored_buf.hh"
#include "vector_sums.hh"
#include "vectors.hh"
#include "metrics.hh"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

void gm_compute_numerators_2way_start(GMVectors* vectors_left,
                                      GMVectors* vectors_right,
                                      GMMetrics* metrics,
                                      GMMirroredBuf* vectors_left_buf,
                                      GMMirroredBuf* vectors_right_buf,
                                      GMMirroredBuf* metrics_buf,
                                      int j_proc,
                                      bool compute_triang_only,
                                      GMEnv* env);

void gm_compute_2way_combine(GMMetrics* metrics,
                             GMMirroredBuf* metrics_buf,
                             const GMVectorSums* vector_sums_left,
                             const GMVectorSums* vector_sums_right,
                             int j_proc,
                             bool do_compute_triang_only,
                             GMEnv* env);

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_gm_compute_metrics_utils_2way_hh_---*/

/*---------------------------------------------------------------------------*/
