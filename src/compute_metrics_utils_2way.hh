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
#include "vector_sums.hh"
#include "vectors.hh"
#include "metrics.hh"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

#if 0
void gm_compute_czekanowski_numerators_2way_start(
    GMVectors* vectors_left,
    GMVectors* vectors_right,
    GMMetrics* numerators,
    GMMirroredPointer* vectors_left_buf,
    GMMirroredPointer* vectors_right_buf,
    GMMirroredPointer* numerators_buf,
    int j_proc,
    _Bool do_compute_triang_only,
    GMEnv* env);

void gm_compute_ccc_numerators_2way_start(GMVectors* vectors_left,
                                          GMVectors* vectors_right,
                                          GMMetrics* numerators,
                                          GMMirroredPointer* vectors_left_buf,
                                          GMMirroredPointer* vectors_right_buf,
                                          GMMirroredPointer* numerators_buf,
                                          int j_proc,
                                          _Bool do_compute_triang_only,
                                          GMEnv* env);
#endif

void gm_compute_numerators_2way_start(GMVectors* vectors_left,
                                      GMVectors* vectors_right,
                                      GMMetrics* numerators,
                                      GMMirroredPointer* vectors_left_buf,
                                      GMMirroredPointer* vectors_right_buf,
                                      GMMirroredPointer* numerators_buf,
                                      int j_proc,
                                      _Bool compute_triang_only,
                                      GMEnv* env);

/*----------*/

#if 0
void gm_compute_czekanowski_2way_combine(
    GMMetrics* metrics,
    GMMirroredPointer* metrics_buf,
    GMFloat* __restrict__ vector_sums_left,
    GMFloat* __restrict__ vector_sums_right,
    int j_proc,
    _Bool compute_triang_only,
    GMEnv* env);

void gm_compute_ccc_2way_combine(GMMetrics* metrics,
                                 GMMirroredPointer* metrics_buf,
                                 GMFloat* vector_sums_left,
                                 GMFloat* vector_sums_right,
                                 int j_proc,
                                 _Bool compute_triang_only,
                                 GMEnv* env);
#endif

void gm_compute_2way_combine(GMMetrics* metrics,
                             GMMirroredPointer* metrics_buf,
                             GMVectorSums* vector_sums_left,
                             GMVectorSums* vector_sums_right,
                             int j_proc,
                             _Bool do_compute_triang_only,
                             GMEnv* env);

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_gm_compute_metrics_utils_2way_hh_---*/

/*---------------------------------------------------------------------------*/
