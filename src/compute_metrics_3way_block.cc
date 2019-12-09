//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way_block.cc
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Calculate numerators, 3-way, for a single block.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "string.h"

#include "env.hh"
#include "linalg.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "comm_xfer_utils.hh"
#include "compute_metrics_3way_block_gpu.hh"
#include "compute_metrics_3way_block_nongpu.hh"
#include "compute_metrics_3way_block.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

ComputeNumerators3Way::ComputeNumerators3Way(int nvl, int npvfl, Env& env)
  : env_(env)
  , tmp_buf_{}
  , matM_ij_buf_{}
  , matM_jk_buf_{}
  , matM_kik_buf_{}
  , matX_buf_{}
  , matB_buf_{} {
  COMET_INSIST(nvl >= 0 && npvfl >= 0);

  if (!env_.is_using_linalg())
    return;

  for (int i=0; i<NUM_BUF; ++i) {
    if (env_.do_reduce()) {
      GMMirroredBuf_create(&tmp_buf_[i], nvl, nvl, &env_);
    }
    GMMirroredBuf_create(&matX_buf_[i], npvfl, nvl, &env_);
    GMMirroredBuf_create(&matB_buf_[i], nvl, nvl, &env_);
  }
  if (env_.does_3way_need_2way()) {
    GMMirroredBuf_create(&matM_ij_buf_, nvl, nvl, &env_);
    GMMirroredBuf_create(&matM_jk_buf_, nvl, nvl, &env_);
    GMMirroredBuf_create(&matM_kik_buf_, nvl, nvl, &env_);
  }
}

//-----------------------------------------------------------------------------

ComputeNumerators3Way::~ComputeNumerators3Way() {

  if (!env_.is_using_linalg())
    return;

  for (int i=0; i<NUM_BUF; ++i) {
    if (env_.do_reduce()) {
      GMMirroredBuf_destroy(&tmp_buf_[i], &env_);
    }
    GMMirroredBuf_destroy(&matX_buf_[i], &env_);
    GMMirroredBuf_destroy(&matB_buf_[i], &env_);
  }
  if (env_.does_3way_need_2way()) {
    GMMirroredBuf_destroy(&matM_ij_buf_, &env_);
    GMMirroredBuf_destroy(&matM_jk_buf_, &env_);
    GMMirroredBuf_destroy(&matM_kik_buf_, &env_);
  }
}

//-----------------------------------------------------------------------------

void ComputeNumerators3Way::compute(
  VData vdata_i, VData vdata_j, VData vdata_k, 
  GMMetrics& numerators, int j_block, int k_block, int section_step) {
  COMET_INSIST(j_block >= 0 && j_block < env_.num_block_vector());
  COMET_INSIST(k_block >= 0 && k_block < env_.num_block_vector());
  COMET_INSIST(! (env_.proc_num_vector() == j_block &&
                  env_.proc_num_vector() != k_block));
  COMET_INSIST(! (env_.proc_num_vector() == k_block &&
                  env_.proc_num_vector() != j_block));

  if (env_.is_using_linalg()) {

    compute_linalg_(vdata_i, vdata_j, vdata_k, numerators, j_block, k_block,
      section_step);

  } else if (env_.metric_type() == MetricType::CZEK)  {

    compute_czek_(vdata_i, vdata_j, vdata_k, numerators, j_block, k_block,
      section_step);

  } else if (env_.metric_type() == MetricType::CCC)  {

    compute_ccc_(vdata_i, vdata_j, vdata_k, numerators, j_block, k_block,
      section_step);

  } else {

    COMET_INSIST_INTERFACE(&env_, false &&
      "Selected metric_type unimplemented.");

  }
}

//-----------------------------------------------------------------------------










//-----------------------------------------------------------------------------

void GMComputeNumerators3Way_create(GMComputeNumerators3Way* this_,
                                    int nvl, int npvfl, GMEnv* env) {
  COMET_INSIST(this_ && env);
  COMET_INSIST(nvl >= 0 && npvfl >= 0);

  this_->matM_ij_buf = GMMirroredBuf_null();
  this_->matM_jk_buf = GMMirroredBuf_null();
  this_->matM_kik_buf = GMMirroredBuf_null();

  for (int i=0; i<2; ++i) {
    this_->tmp_buf[i] = GMMirroredBuf_null();
    this_->matX_buf[i] = GMMirroredBuf_null();
    this_->matB_buf[i] = GMMirroredBuf_null();
  }

  if (env->is_using_linalg()) {
    for (int i=0; i<2; ++i) {
      if (env->do_reduce()) {
        GMMirroredBuf_create(&(this_->tmp_buf[i]), nvl, nvl, env);
      }
      GMMirroredBuf_create(&(this_->matX_buf[i]), npvfl, nvl, env);
      GMMirroredBuf_create(&(this_->matB_buf[i]), nvl, nvl, env);
    }
    if (env->does_3way_need_2way()) {
      GMMirroredBuf_create(&(this_->matM_ij_buf), nvl, nvl, env);
      GMMirroredBuf_create(&(this_->matM_jk_buf), nvl, nvl, env);
      GMMirroredBuf_create(&(this_->matM_kik_buf), nvl, nvl, env);
    }
  }
}

//=============================================================================

void GMComputeNumerators3Way_destroy(GMComputeNumerators3Way* this_,
                                     GMEnv* env) {
  COMET_INSIST(this_ && env);

  if (env->is_using_linalg()) {
    for (int i=0; i<2; ++i) {
      if (env->do_reduce()) {
        GMMirroredBuf_destroy(&this_->tmp_buf[i], env);
      }
      GMMirroredBuf_destroy(&this_->matX_buf[i], env);
      GMMirroredBuf_destroy(&this_->matB_buf[i], env);
    }
    if (env->does_3way_need_2way()) {
      GMMirroredBuf_destroy(&this_->matM_ij_buf, env);
      GMMirroredBuf_destroy(&this_->matM_jk_buf, env);
      GMMirroredBuf_destroy(&this_->matM_kik_buf, env);
    }
  }
}

//=============================================================================
/*---Start calculation of numerators, 3-way generic---*/

void GMComputeNumerators3Way_start(
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
    GMEnv* env) {
  COMET_INSIST(this_ && metrics && env);
  COMET_INSIST(vectors_i && vectors_j && vectors_k);
  COMET_INSIST(vectors_i_buf && vectors_j_buf && vectors_k_buf);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(k_block >= 0 && k_block < env->num_block_vector());
  COMET_INSIST(! (env->proc_num_vector() == j_block &&
              env->proc_num_vector() != k_block));
  COMET_INSIST(! (env->proc_num_vector() == k_block &&
              env->proc_num_vector() != j_block));
  COMET_INSIST(env->num_way() == NUM_WAY::_3);
  COMET_INSIST(vector_sums_i && vector_sums_j && vector_sums_k);

  /*----------------------------------------*/
  if (env->is_using_linalg()) {
    /*----------------------------------------*/
    gm_compute_3way_nums_gpu_start_(this_,
                                    vectors_i, vectors_j, vectors_k,
                                    metrics, vectors_i_buf, vectors_j_buf,
                                    vectors_k_buf, j_block, k_block,
                                    vector_sums_i, vector_sums_j,
                                    vector_sums_k,
                                    section_step, env);
    /*----------------------------------------*/
  } else /*---(!env->is_using_linalg())---*/ {
    /*----------------------------------------*/
    switch (env->metric_type()) {
      case MetricType::CZEK: {
        gm_compute_3way_nums_nongpu_czek_start_(this_,
            vectors_i, vectors_j, vectors_k, metrics, vectors_i_buf,
            vectors_j_buf, vectors_k_buf, j_block, k_block,
            vector_sums_i, vector_sums_j, vector_sums_k,
            section_step, env);
      } break;
      case MetricType::CCC: {
        gm_compute_3way_nums_nongpu_ccc_start_(this_,
            vectors_i, vectors_j, vectors_k, metrics, vectors_i_buf,
            vectors_j_buf, vectors_k_buf, j_block, k_block,
            vector_sums_i, vector_sums_j, vector_sums_k,
            section_step, env);
      } break;
      default:
        COMET_INSIST_INTERFACE(env, false && "Selected metric_type unimplemented.");
    } /*---case---*/
    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
