//-----------------------------------------------------------------------------
/*!
 * \file   linalg.hh
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Interface to generalized linear algebra functions, e.g. MAGMA.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_linalg_hh_
#define _gm_linalg_hh_

#include "env.hh"
#include "decomp_mgr.hh"
#include "mirrored_buf.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

void gm_linalg_initialize(GMEnv* env);

void gm_linalg_finalize(GMEnv* env);

/*----------*/

void gm_linalg_malloc(GMMirroredBuf* p, size_t dim0, size_t dim1, GMEnv* env);

void gm_linalg_free(GMMirroredBuf* p, GMEnv* env);

void gm_linalg_set_matrix_zero_start(GMMirroredBuf* matrix_buf,
                                     GMEnv* env);

/*----------*/

void gm_linalg_gemm_start(size_t m,
                          size_t n,
                          size_t k,
                          void* dA,
                          size_t ldda,
                          void* dB,
                          size_t lddb,
                          void* dC,
                          size_t lddc,
                          GMDecompMgr* dm,
                          GMEnv* env);

void gm_compute_wait(GMEnv* env);

/*----------*/

void gm_linalg_set_matrix_start(GMMirroredBuf* matrix_buf, GMEnv* env);

void gm_linalg_set_matrix_wait(GMEnv* env);

void gm_linalg_get_matrix_start(GMMirroredBuf* matrix_buf, GMEnv* env);

void gm_linalg_get_matrix_wait(GMEnv* env);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _gm_linalg_hh_

//-----------------------------------------------------------------------------
