//-----------------------------------------------------------------------------
/*!
 * \file   linalg.cc
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Perform (actual or modified) GEMM operations.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifdef COMET_USE_MAGMA
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "magma_mgemm2.h"
#include "magma_mgemm2_lapack.h"
#include "magma_mgemm3.h"
#include "magma_mgemm3_lapack.h"
#include "magma_mgemm4.h"
#include "magma_mgemm4_lapack.h"
#include "magma_mgemm5.h"
#include "magma_mgemm5_lapack.h"
//#elif defined COMET_USE_CUDA
//  #include "cublas_v2.h"
#elif defined COMET_USE_HIP
  #include "hip/hip_runtime_api.h"
//  #include "hip/hip_runtime.h"
#endif

#include "env.hh"
#include "assertions.hh"
#include "tc.hh"
#include "decomp_mgr.hh"
#include "magma_wrapper.hh"
#include "linalg.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

void gm_linalg_gemm_start(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA1, const MirroredBuf* matA2,
  const MirroredBuf* matB, MirroredBuf* matC,
  MirroredBuf* sums_I, MirroredBuf* sums_J, MirroredBuf* sums_K,
  MirroredBuf* counts_I, MirroredBuf* counts_J, MirroredBuf* counts_K, int J,
  int step_2way, GMDecompMgr* dm, CEnv* env) {
  COMET_INSIST(matA1 && matA2 && matB && matC && env);
  COMET_ASSERT(sums_I && sums_J && sums_K && counts_I && counts_J && counts_K);

  if (m==0 || n==0 || k==0)
    return;

  if (env->is_compute_method_gpu()) {
    matA1->lock_d();
    if (matB != matA1) {
      matB->lock_d();
    }
    matC->lock_d();
  }

  if (env->is_using_tc()) {
    if (env->is_compute_method_gpu()) {
      tc_gemm_start(m, n, k,
        matA1->active, matA1->dim0, matA2->active, matA2->dim0,
        matB->active, matB->dim0, matC->active, matC->dim0,
//        0,0,0,0,0,0,0,
        (GMFloat*)sums_I->active, (GMFloat*)sums_J->active, (GMFloat*)sums_K->active,
        (GMFloat*)counts_I->active, (GMFloat*)counts_J->active, (GMFloat*)counts_K->active, J,
        dm->num_field_active_local, step_2way, dm->tc_bufs, *env);
    }
  } else {
    gm_linalg_set_matrix_zero_start_(matC, env); // apparently needed by magma.
    gm_linalg_gemm_magma_start(m, n, k, matA1->active, matA1->dim0,
      matB->active, matB->dim0, matC->active, matC->dim0, env);
  }
}

//-----------------------------------------------------------------------------

void gm_linalg_gemm_wait(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA1, const MirroredBuf* matA2,
  const MirroredBuf* matB, MirroredBuf* matC,
  MirroredBuf* sums_I, MirroredBuf* sums_J, MirroredBuf* sums_K,
  MirroredBuf* counts_I, MirroredBuf* counts_J, MirroredBuf* counts_K, int J,
  int step_2way, GMDecompMgr* dm, CEnv* env) {
  COMET_INSIST(matA1 && matA2 && matB && matC && env);
  COMET_ASSERT(sums_I && sums_J && sums_K && counts_I && counts_J && counts_K);

  if (m==0 || n==0 || k==0)
    return;

  if (env->is_using_tc()) {
    if (!env->is_compute_method_gpu()) {
      matA1->lock_h();
      if (matA2 != matA1 && matA2 != matB) {
        matA2->lock_h();
      }
      if (matB != matA1) {
        matB->lock_h();
      }
      matC->lock_h();
      tc_gemm_start(m, n, k,
        matA1->active, matA1->dim0, matA2->active, matA2->dim0,
        matB->active, matB->dim0, matC->active, matC->dim0,
//        0,0,0,0,0,0,0,
        (GMFloat*)sums_I->active, (GMFloat*)sums_J->active, (GMFloat*)sums_K->active,
        (GMFloat*)counts_I->active, (GMFloat*)counts_J->active, (GMFloat*)counts_K->active, J,
        dm->num_field_active_local, step_2way, dm->tc_bufs, *env);
      matA1->unlock_h();
      if (matA2 != matA1 && matA2 != matB) {
        matA2->unlock_h();
      }
      if (matB != matA1) {
        matB->unlock_h();
      }
      matC->unlock_h();
    }
  }

  env->stream_synchronize(env->stream_compute());

  if (env->is_compute_method_gpu()) {
    matA1->unlock_d();
    if (matB != matA1) {
      matB->unlock_d();
    }
    matC->unlock_d();
  }
}

//-----------------------------------------------------------------------------

void gm_linalg_gemm_start(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA, const MirroredBuf* matB, MirroredBuf* matC,
  MirroredBuf* sums_I, MirroredBuf* sums_J,
  MirroredBuf* counts_I, MirroredBuf* counts_J,
  GMDecompMgr* dm, CEnv* env) {

  gm_linalg_gemm_start(m, n, k, matA, matA, matB, matC,
    sums_I, sums_J, sums_J, counts_I, counts_J, counts_J, 0, 0, dm, env);
}

//-----------------------------------------------------------------------------

void gm_linalg_gemm_wait(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA, const MirroredBuf* matB, MirroredBuf* matC,
  MirroredBuf* sums_I, MirroredBuf* sums_J,
  MirroredBuf* counts_I, MirroredBuf* counts_J,
  GMDecompMgr* dm, CEnv* env) {

  gm_linalg_gemm_wait(m, n, k, matA, matA, matB, matC,
    sums_I, sums_J, sums_J, counts_I, counts_J, counts_J, 0, 0, dm, env);
}

//-----------------------------------------------------------------------------

void gm_linalg_gemm(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA, const MirroredBuf* matB, MirroredBuf* matC,
  MirroredBuf* sums_I, MirroredBuf* sums_J,
  MirroredBuf* counts_I, MirroredBuf* counts_J,
  GMDecompMgr* dm, CEnv* env) {
  COMET_INSIST(matA && matB && matC && env);

  gm_linalg_gemm_start(m, n, k, matA, matB, matC,
    sums_I, sums_J, counts_I, counts_J, dm, env);
  gm_linalg_gemm_wait(m, n, k, matA, matB, matC,
    sums_I, sums_J, counts_I, counts_J, dm, env);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
