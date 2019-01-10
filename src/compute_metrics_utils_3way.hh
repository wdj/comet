//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_utils_3way.hh
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing metrics, utilities, 3-way case, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_compute_metrics_utils_3way_hh_
#define _gm_compute_metrics_utils_3way_hh_

#include "env.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"

//=============================================================================

typedef struct {
  GMMirroredBuf tmp_buf[2];
  GMMirroredBuf matM_ij_buf;
  GMMirroredBuf matM_jk_buf;
  GMMirroredBuf matM_kik_buf;
  GMMirroredBuf matX_buf[2];
  GMMirroredBuf matB_buf[2];
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
    GMMirroredBuf* vectors_i_buf,
    GMMirroredBuf* vectors_j_buf,
    GMMirroredBuf* vectors_k_buf,
    int j_proc,
    int k_proc,
    const GMVectorSums* vector_sums_i,
    const GMVectorSums* vector_sums_j,
    const GMVectorSums* vector_sums_k,
    int section_step,
    GMEnv* env);

#if 0
void gm_compute_numerators_3way_gpu_form_matX_(
  const GMVectors* vectors_i,
  const GMMirroredBuf* vectors_I_buf,
  const GMMirroredBuf* vectors_J_buf,
  GMMirroredBuf* const matX_buf,
  const int J,
  const int step_2way,
  const int I_min,
  const int I_max,
  GMEnv* const env);

void gm_compute_numerators_3way_gpu_form_metrics_(
  GMMirroredBuf* const matM_IJ_buf,
  GMMirroredBuf* const matM_JK_buf,
  GMMirroredBuf* const matM_KIK_buf,
  GMMirroredBuf* const matB_buf,
  GMMetrics* metrics,
  const int nvl,
  const int J,
  const int step_2way,
  const int I_min,
  const int I_max,
  const int K_min,
  const int K_max,
  const int j_block,
  const int k_block,
  const GMSectionInfo* const si,
  const GMVectorSums* vector_sums_i,
  const GMVectorSums* vector_sums_j,
  const GMVectorSums* vector_sums_k,
  GMEnv* const env);
#endif

void gm_compute_numerators_3way_gpu_start_(
    GMComputeNumerators3Way* this_,
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* metrics,
    GMMirroredBuf* vectors_i_buf,
    GMMirroredBuf* vectors_j_buf,
    GMMirroredBuf* vectors_k_buf,
    int j_block,
    int k_block,
    const GMVectorSums* vector_sums_i,
    const GMVectorSums* vector_sums_j,
    const GMVectorSums* vector_sums_k,
    int section_step,
    GMEnv* env);

//=============================================================================

#endif /*---_gm_compute_metrics_utils_3way_hh_---*/

//-----------------------------------------------------------------------------
