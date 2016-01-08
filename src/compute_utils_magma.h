/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_utils_magma.h
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Magma utilities, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _compute_utils_magma_h_
#define _compute_utils_magma_h_

#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"

#include "env.h"
#include "vector_sums.h"
#include "vectors.h"
#include "metrics.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

void gm_magma_initialize(GMEnv* env);

void gm_magma_finalize(GMEnv* env);

/*----------*/

GMMirroredPointer gm_malloc_magma(size_t n, GMEnv* env);

void gm_free_magma(GMMirroredPointer* p, GMEnv* env);

void gm_magma_set_matrix_zero_start(GMMirroredPointer* matrix_buf,
                                    int mat_dim1,
                                    int mat_dim2,
                                    GMEnv* env);

void gm_magma_gemm_start(magma_minproduct_int_t m,
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

void gm_set_matrix_start(GMMirroredPointer* matrix_buf,
                         int mat_dim1,
                         int mat_dim2,
                         GMEnv* env);

void gm_set_matrix_wait(GMEnv* env);

void gm_get_matrix_start(GMMirroredPointer* matrix_buf,
                         int mat_dim1,
                         int mat_dim2,
                         GMEnv* env);

void gm_get_matrix_wait(GMEnv* env);

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_compute_utils_magma_h---*/

/*---------------------------------------------------------------------------*/
