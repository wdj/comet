//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way_block_nongpu.hh
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Calculate numerators, 3-way, for a single block, non-GPU case.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_compute_metrics_3way_block_nongpu_hh_
#define _gm_compute_metrics_3way_block_nongpu_hh_

#include "env.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "compute_metrics_3way_block.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

void gm_compute_3way_nums_nongpu_czek_start_(
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

//-----------------------------------------------------------------------------

void gm_compute_3way_nums_nongpu_ccc_start_(
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

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _gm_compute_metrics_3way_block_nongpu_hh_

//-----------------------------------------------------------------------------
