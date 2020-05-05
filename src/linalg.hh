//-----------------------------------------------------------------------------
/*!
 * \file   linalg.hh
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Perform (actual or modified) GEMM operations.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_linalg_hh_
#define _comet_linalg_hh_

#include "env.hh"
#include "decomp_mgr.hh"
#include "mirrored_buf.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

#if 0
void gm_linalg_initialize(CEnv* env);

void gm_linalg_finalize(CEnv* env);

/*----------*/

void gm_linalg_malloc(MirroredBuf* p, size_t dim0, size_t dim1, CEnv* env);

void gm_linalg_free(MirroredBuf* p, CEnv* env);

/*----------*/
#endif

void gm_linalg_gemm_start(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA1, const MirroredBuf* matA2,
  const MirroredBuf* matB, MirroredBuf* matC,
  MirroredBuf* sums_I, MirroredBuf* sums_J, MirroredBuf* sums_K,
  MirroredBuf* counts_I, MirroredBuf* counts_J, MirroredBuf* counts_K, int J,
  int step_2way, GMDecompMgr* dm, CEnv* env);

void gm_linalg_gemm_wait(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA1, const MirroredBuf* matA2,
  const MirroredBuf* matB, MirroredBuf* matC,
  MirroredBuf* sums_I, MirroredBuf* sums_J, MirroredBuf* sums_K,
  MirroredBuf* counts_I, MirroredBuf* counts_J, MirroredBuf* counts_K, int J,
  int step_2way, GMDecompMgr* dm, CEnv* env);

void gm_linalg_gemm_start(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA, const MirroredBuf* matB, MirroredBuf* matC,
  MirroredBuf* sums_I, MirroredBuf* sums_J,
  MirroredBuf* counts_I, MirroredBuf* counts_J,
  GMDecompMgr* dm, CEnv* env);

void gm_linalg_gemm_wait(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA, const MirroredBuf* matB, MirroredBuf* matC,
  MirroredBuf* sums_I, MirroredBuf* sums_J,
  MirroredBuf* counts_I, MirroredBuf* counts_J,
  GMDecompMgr* dm, CEnv* env);

void gm_linalg_gemm(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA, const MirroredBuf* matB, MirroredBuf* matC,
  MirroredBuf* sums_I, MirroredBuf* sums_J,
  MirroredBuf* counts_I, MirroredBuf* counts_J,
  GMDecompMgr* dm, CEnv* env);

#if 0
/*----------*/

void gm_linalg_set_matrix_start(MirroredBuf* matrix_buf, CEnv* env);

void gm_linalg_set_matrix_wait(CEnv* env);

void gm_linalg_get_matrix_start(MirroredBuf* matrix_buf, CEnv* env);

void gm_linalg_get_matrix_wait(CEnv* env);
#endif

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_linalg_hh_

//-----------------------------------------------------------------------------
