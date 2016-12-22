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

#if 0
void gm_compute_czekanowski_numerators_3way_nongpu_start(
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* numerators,
    GMMirroredPointer* vectors_i_buf,
    GMMirroredPointer* vectors_j_buf,
    GMMirroredPointer* vectors_k_buf,
    int j_proc,
    int k_proc,
    GMVectorSums* vector_sums_i,
    GMVectorSums* vector_sums_j,
    GMVectorSums* vector_sums_k,
    int section_step,
    GMEnv* env);

void gm_compute_ccc_numerators_3way_nongpu_start(
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* numerators,
    GMMirroredPointer* vectors_i_buf,
    GMMirroredPointer* vectors_j_buf,
    GMMirroredPointer* vectors_k_buf,
    int j_proc,
    int k_proc,
    GMVectorSums* vector_sums_i,
    GMVectorSums* vector_sums_j,
    GMVectorSums* vector_sums_k,
    int section_step,
    GMEnv* env);

void gm_compute_numerators_3way_gpu_start(
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* numerators,
    GMMirroredPointer* vectors_i_buf,
    GMMirroredPointer* vectors_j_buf,
    GMMirroredPointer* vectors_k_buf,
    int j_proc,
    int k_proc,
    GMVectorSums* vector_sums_i,
    GMVectorSums* vector_sums_j,
    GMVectorSums* vector_sums_k,
    int section_step,
    GMEnv* env);
#endif

void gm_compute_numerators_3way_start(
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* numerators,
    GMMirroredPointer* vectors_i_buf,
    GMMirroredPointer* vectors_j_buf,
    GMMirroredPointer* vectors_k_buf,
    int j_proc,
    int k_proc,
    GMVectorSums* vector_sums_i,
    GMVectorSums* vector_sums_j,
    GMVectorSums* vector_sums_k,
    int section_step,
    GMEnv* env);

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_gm_compute_metrics_utils_3way_hh_---*/

/*---------------------------------------------------------------------------*/
