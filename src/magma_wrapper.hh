//-----------------------------------------------------------------------------
/*!
 * \file   magma_wrapper.hh
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Interface to generalized linear algebra functions, e.g. MAGMA.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_magma_wrapper_hh_
#define _comet_magma_wrapper_hh_

#include "env.hh"
#include "mirrored_buf.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Helpers

static bool use_minproduct(CEnv* env) {
  return env->metric_type() == MetricType::CZEK;
}

static bool use_mgemm2(CEnv* env) {
  return env->metric_type() == MetricType::CCC &&
         env->num_way() == NumWay::_2 && ! env->sparse();
}

static bool use_mgemm3(CEnv* env) {
  return env->metric_type() == MetricType::CCC &&
         env->num_way() == NumWay::_3 && ! env->sparse();
}

static bool use_mgemm4(CEnv* env) {
  return env->metric_type() == MetricType::CCC && env->sparse();
}

static bool use_mgemm5(CEnv* env) {
  return env->metric_type() == MetricType::DUO;
}

//-----------------------------------------------------------------------------

void gm_linalg_initialize(CEnv* env);

void gm_linalg_finalize(CEnv* env);

//----------

void gm_linalg_malloc(MirroredBuf* p, size_t dim0, size_t dim1, CEnv* env);

void gm_linalg_free(MirroredBuf* p, CEnv* env);

//----------

void gm_linalg_set_matrix_zero_start_(MirroredBuf* matrix_buf, CEnv* env);

void gm_linalg_gemm_magma_start(size_t m, size_t n, size_t k,
  const void* matA, size_t ldda, const void* matB, size_t lddb,
  void* matC, size_t lddc, CEnv* env);

//----------

void gm_linalg_set_matrix_start(MirroredBuf* matrix_buf, CEnv* env);

void gm_linalg_set_matrix_wait(CEnv* env);

void gm_linalg_get_matrix_start(MirroredBuf* matrix_buf, CEnv* env);

void gm_linalg_get_matrix_wait(CEnv* env);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_magma_wrapper_hh_

//-----------------------------------------------------------------------------
