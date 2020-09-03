//-----------------------------------------------------------------------------
/*!
 * \file   linalg.cc
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Perform (actual or modified) GEMM operations.
 */
//-----------------------------------------------------------------------------
/*-----------------------------------------------------------------------------

Copyright 2020, UT-Battelle, LLC

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-----------------------------------------------------------------------------*/

#ifdef COMET_USE_MAGMA
//#include "magma_minproduct.h"
//#include "magma_minproduct_lapack.h"
//#include "magma_mgemm2.h"
//#include "magma_mgemm2_lapack.h"
//#include "magma_mgemm3.h"
//#include "magma_mgemm3_lapack.h"
//#include "magma_mgemm4.h"
//#include "magma_mgemm4_lapack.h"
//#include "magma_mgemm5.h"
//#include "magma_mgemm5_lapack.h"
//#elif defined COMET_USE_CUDA
//  #include "cublas_v2.h"
#elif defined COMET_USE_HIP
//  #include "hip/hip_runtime_api.h"
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
/// \brief GEMM start, A matrix+column case.

void LinAlg::gemm_start(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA1, const MirroredBuf* matA2,
  const MirroredBuf* matB, MirroredBuf* matC,
  MirroredBuf* sums_I, MirroredBuf* sums_J, MirroredBuf* sums_K,
  MirroredBuf* counts_I, MirroredBuf* counts_J, MirroredBuf* counts_K, int J,
  int step_2way, GMDecompMgr& dm, MagmaWrapper& magma_wrapper, CEnv& env) {
  COMET_INSIST(matA1 && matA2 && matB && matC);
  COMET_ASSERT(sums_I && sums_J && sums_K && counts_I && counts_J && counts_K);

  if (!m || !n || !k)
    return;

  // Lock.

  if (env.is_compute_method_gpu()) {
    matA1->lock_d();
    if (matB != matA1)
      matB->lock_d();
    matC->lock_d();
  }

  // GPU case.

  if (env.is_using_tc()) {
    if (env.is_compute_method_gpu())
      // GEMM call, tc case.
      tc_gemm_start(m, n, k,
        matA1->active, matA1->dim0, matA2->active, matA2->dim0,
        matB->active, matB->dim0, matC->active, matC->dim0,
        (GMFloat*)sums_I->active, (GMFloat*)sums_J->active, (GMFloat*)sums_K->active,
        (GMFloat*)counts_I->active, (GMFloat*)counts_J->active, (GMFloat*)counts_K->active, J,
        dm.num_field_active_local, step_2way, dm.tc_bufs, env);
  } else {
    // apparently needed by magma.
    MagmaWrapper::set_matrix_zero_start(matC, env);
    // GEMM call, non-tc case.
    MagmaWrapper::gemm_start(m, n, k, matA1->active, matA1->dim0,
      matB->active, matB->dim0, matC->active, matC->dim0, env);
  }
}

//-----------------------------------------------------------------------------
/// \brief GEMM wait, A matrix+column case.

void LinAlg::gemm_wait(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA1, const MirroredBuf* matA2,
  const MirroredBuf* matB, MirroredBuf* matC,
  MirroredBuf* sums_I, MirroredBuf* sums_J, MirroredBuf* sums_K,
  MirroredBuf* counts_I, MirroredBuf* counts_J, MirroredBuf* counts_K, int J,
  int step_2way, GMDecompMgr& dm, CEnv& env) {
  COMET_INSIST(matA1 && matA2 && matB && matC);
  COMET_ASSERT(sums_I && sums_J && sums_K && counts_I && counts_J && counts_K);

  if (!m || !n || !k)
    return;

  // CPU case.

  if (env.is_using_tc()) {
    if (!env.is_compute_method_gpu()) {
      // Lock
      matA1->lock_h();
      if (matA2 != matA1 && matA2 != matB)
        matA2->lock_h();
      if (matB != matA1)
        matB->lock_h();
      matC->lock_h();
      // GEMM call, tc case.
      tc_gemm_start(m, n, k,
        matA1->active, matA1->dim0, matA2->active, matA2->dim0,
        matB->active, matB->dim0, matC->active, matC->dim0,
        (GMFloat*)sums_I->active, (GMFloat*)sums_J->active, (GMFloat*)sums_K->active,
        (GMFloat*)counts_I->active, (GMFloat*)counts_J->active, (GMFloat*)counts_K->active, J,
        dm.num_field_active_local, step_2way, dm.tc_bufs, env);
      // Unlock
      matA1->unlock_h();
      if (matA2 != matA1 && matA2 != matB)
        matA2->unlock_h();
      if (matB != matA1)
        matB->unlock_h();
      matC->unlock_h();
    }
  }

  env.stream_synchronize(env.stream_compute());

  // Unlock.

  if (env.is_compute_method_gpu()) {
    matA1->unlock_d();
    if (matB != matA1)
      matB->unlock_d();
    matC->unlock_d();
  }
}

//-----------------------------------------------------------------------------
/// \brief GEMM start, standard A case.

void LinAlg::gemm_start(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA, const MirroredBuf* matB, MirroredBuf* matC,
  MirroredBuf* sums_I, MirroredBuf* sums_J,
  MirroredBuf* counts_I, MirroredBuf* counts_J,
  GMDecompMgr& dm, MagmaWrapper& magma_wrapper, CEnv& env) {

  gemm_start(m, n, k, matA, matA, matB, matC,
    sums_I, sums_J, sums_J, counts_I, counts_J, counts_J, 0, 0, dm,
    magma_wrapper, env);
}

//-----------------------------------------------------------------------------
/// \brief GEMM wait, standard A case.

void LinAlg::gemm_wait(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA, const MirroredBuf* matB, MirroredBuf* matC,
  MirroredBuf* sums_I, MirroredBuf* sums_J,
  MirroredBuf* counts_I, MirroredBuf* counts_J,
  GMDecompMgr& dm, CEnv& env) {

  gemm_wait(m, n, k, matA, matA, matB, matC,
    sums_I, sums_J, sums_J, counts_I, counts_J, counts_J, 0, 0, dm, env);
}

//-----------------------------------------------------------------------------
/// \brief GEMM start and wait, standard A case.

void LinAlg::gemm(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA, const MirroredBuf* matB, MirroredBuf* matC,
  MirroredBuf* sums_I, MirroredBuf* sums_J,
  MirroredBuf* counts_I, MirroredBuf* counts_J,
  GMDecompMgr& dm, MagmaWrapper& magma_wrapper, CEnv& env) {
  COMET_INSIST(matA && matB && matC);

  gemm_start(m, n, k, matA, matB, matC,
    sums_I, sums_J, counts_I, counts_J, dm, magma_wrapper, env);
  gemm_wait(m, n, k, matA, matB, matC,
    sums_I, sums_J, counts_I, counts_J, dm, env);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
