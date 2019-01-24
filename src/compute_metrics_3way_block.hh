//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way_block.hh
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Calculate numerators, 3-way, for a single block.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_compute_metrics_3way_block_hh_
#define _gm_compute_metrics_3way_block_hh_

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

//=============================================================================

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

//=============================================================================

#endif // _gm_compute_metrics_3way_block_hh_

//-----------------------------------------------------------------------------
