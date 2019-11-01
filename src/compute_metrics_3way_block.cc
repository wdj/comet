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

void GMComputeNumerators3Way_create(GMComputeNumerators3Way* this_,
                                    int nvl, int npvfl, GMEnv* env) {
  GMInsist(this_ && env);
  GMInsist(nvl >= 0 && npvfl >= 0);

  this_->matM_ij_buf = GMMirroredBuf_null();
  this_->matM_jk_buf = GMMirroredBuf_null();
  this_->matM_kik_buf = GMMirroredBuf_null();

  for (int i=0; i<2; ++i) {
    this_->tmp_buf[i] = GMMirroredBuf_null();
    this_->matX_buf[i] = GMMirroredBuf_null();
    this_->matB_buf[i] = GMMirroredBuf_null();
  }

  if (env->compute_method() == ComputeMethod::GPU) {
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
  GMInsist(this_ && env);

  if (env->compute_method() == ComputeMethod::GPU) {
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
  GMInsist(this_ && metrics && env);
  GMInsist(vectors_i && vectors_j && vectors_k);
  GMInsist(vectors_i_buf && vectors_j_buf && vectors_k_buf);
  GMInsist(j_block >= 0 && j_block < env->num_block_vector());
  GMInsist(k_block >= 0 && k_block < env->num_block_vector());
  GMInsist(! (env->proc_num_vector() == j_block &&
              env->proc_num_vector() != k_block));
  GMInsist(! (env->proc_num_vector() == k_block &&
              env->proc_num_vector() != j_block));
  GMInsist(env->num_way() == NUM_WAY::_3);
  GMInsist(vector_sums_i && vector_sums_j && vector_sums_k);

  /*----------------------------------------*/
  if (env->compute_method() == ComputeMethod::GPU) {
    /*----------------------------------------*/
    gm_compute_3way_nums_gpu_start_(this_,
                                    vectors_i, vectors_j, vectors_k,
                                    metrics, vectors_i_buf, vectors_j_buf,
                                    vectors_k_buf, j_block, k_block,
                                    vector_sums_i, vector_sums_j,
                                    vector_sums_k,
                                    section_step, env);
    /*----------------------------------------*/
  } else /*---(env->compute_method() != ComputeMethod::GPU)---*/ {
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
        GMInsistInterface(env, false && "Selected metric_type unimplemented.");
    } /*---case---*/
    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
