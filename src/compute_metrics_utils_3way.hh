/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_utils_3way.hh
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing metrics, utilities, 3-way case, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _gm_compute_metrics_utils_3way_hh_
#define _gm_compute_metrics_utils_3way_hh_

#include "env.hh"
#include "vector_sums.hh"
#include "vectors.hh"
#include "metrics.hh"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

typedef struct {
  GMMirroredPointer mat_buf_tmp[2];
  GMMirroredPointer matM_ij_buf;
  GMMirroredPointer matM_jk_buf;
  GMMirroredPointer matM_kik_buf;
  GMMirroredPointer matV_buf[2];
  GMMirroredPointer matB_buf[2];
} GMComputeNumerators3Way;

void GMComputeNumerators3Way_create(
    GMComputeNumerators3Way* this_,
    int nvl,
    int npvfl,
    GMEnv* env);

void GMComputeNumerators3Way_destroy(
    GMComputeNumerators3Way* this_,
    GMEnv* env);

void GMComputeNumerators3Way_start(
    GMComputeNumerators3Way* this_,
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* numerators,
    GMMirroredPointer* vectors_i_buf,
    GMMirroredPointer* vectors_j_buf,
    GMMirroredPointer* vectors_k_buf,
    int j_proc,
    int k_proc,
    const GMVectorSums* vector_sums_i,
    const GMVectorSums* vector_sums_j,
    const GMVectorSums* vector_sums_k,
    int section_step,
    GMEnv* env);

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_gm_compute_metrics_utils_3way_hh_---*/

/*---------------------------------------------------------------------------*/
