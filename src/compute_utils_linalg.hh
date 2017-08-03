/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_utils_linalg.hh
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Magma utilities, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _gm_compute_utils_magma_hh_
#define _gm_compute_utils_magma_hh_

#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"

#include "env.hh"
#include "mirrored_buf.hh"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

void gm_linalg_initialize(GMEnv* env);

void gm_linalg_finalize(GMEnv* env);

/*----------*/

void gm_linalg_malloc(GMMirroredBuf* p, size_t dim0, size_t dim1, GMEnv* env);

void gm_linalg_free(GMMirroredBuf* p, GMEnv* env);

void gm_linalg_set_matrix_zero_start(GMMirroredBuf* matrix_buf,
                                     GMEnv* env);

/*----------*/

void gm_linalg_gemm_start(magma_minproduct_int_t m,
                          magma_minproduct_int_t n,
                          magma_minproduct_int_t k,
                          void* dA,
                          magma_minproduct_int_t ldda,
                          void* dB,
                          magma_minproduct_int_t lddb,
                          void* dC,
                          magma_minproduct_int_t lddc,
                          GMEnv* env);

void gm_compute_wait(GMEnv* env);

/*----------*/

void gm_linalg_set_matrix_start(GMMirroredBuf* matrix_buf, GMEnv* env);

void gm_linalg_set_matrix_wait(GMEnv* env);

void gm_linalg_get_matrix_start(GMMirroredBuf* matrix_buf, GMEnv* env);

void gm_linalg_get_matrix_wait(GMEnv* env);

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_gm_compute_utils_magma_hh_---*/

/*---------------------------------------------------------------------------*/
