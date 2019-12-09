//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way_block.hh
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Calculate numerators, 3-way, for a single block.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_compute_metrics_3way_block_hh_
#define _comet_compute_metrics_3way_block_hh_

#include "env.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Class for computing numerators for 3-way methods.

class ComputeNumerators3Way {

  enum {NUM_BUF = 2};

public:

  ComputeNumerators3Way(int nvl, int npvfl, Env& env);
  ~ComputeNumerators3Way();

  void compute(VData vdata_i, VData vdata_j, VData vdata_k,
    GMMetrics& numerators, int j_block, int k_block, int section_step);

private:

  Env& env_;

  GMMirroredBuf tmp_buf_[NUM_BUF];
  GMMirroredBuf matM_ij_buf_;
  GMMirroredBuf matM_jk_buf_;
  GMMirroredBuf matM_kik_buf_;
  GMMirroredBuf matX_buf_[NUM_BUF];
  GMMirroredBuf matB_buf_[NUM_BUF];

  void compute_linalg_(VData vdata_i, VData vdata_j, VData vdata_k,
    GMMetrics& numerators, int j_block, int k_block, int section_step);


  void compute_linalg_matX_(VData vdata_i, VData vdata_j,
    const GMMirroredBuf& matX_buf,
    const int J, const int step_2way, const int I_min, const int I_max);

  void compute_linalg_metrics_(
    const GMMirroredBuf& matM_IJ_buf,
    const GMMirroredBuf& matM_JK_buf,
    const GMMirroredBuf& matM_KIK_buf,
    const GMMirroredBuf& matB_buf,
    GMMetrics& metrics,
    const int nvl, const int J, const int step_2way,
    const int I_min, const int I_max, const int K_min, const int K_max,
    const int j_block, const int k_block,
    const GMSectionInfo& si,
    const GMVectorSums& vector_sums_i,
    const GMVectorSums& vector_sums_j,
    const GMVectorSums& vector_sums_k);


  void compute_czek_(VData vdata_i, VData vdata_j, VData vdata_k,
    GMMetrics& numerators, int j_block, int k_block, int section_step);

  void compute_ccc_(VData vdata_i, VData vdata_j, VData vdata_k,
    GMMetrics& numerators, int j_block, int k_block, int section_step);

  // Disallowed methods.
  ComputeNumerators3Way(const ComputeNumerators3Way&);
  void operator=(const ComputeNumerators3Way&);
};

//-----------------------------------------------------------------------------










//-----------------------------------------------------------------------------

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

#endif // _comet_compute_metrics_3way_block_hh_

//-----------------------------------------------------------------------------
