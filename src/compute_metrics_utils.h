/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_utils.h
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing metrics, utilities, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _compute_metrics_utils_h_
#define _compute_metrics_utils_h_

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

void magma_gemm_start(magma_minproduct_int_t m,
                      magma_minproduct_int_t n,
                      magma_minproduct_int_t k,
                      void* dA,
                      magma_minproduct_int_t ldda,
                      void* dB,
                      magma_minproduct_int_t lddb,
                      void* dC,
                      magma_minproduct_int_t lddc);

void gm_compute_wait(GMEnv* env);

/*----------*/

MPI_Request gm_send_vectors_start(GMVectors* vectors, int proc_num, GMEnv* env);

MPI_Request gm_recv_vectors_start(GMVectors* vectors, int proc_num, GMEnv* env);

void gm_send_vectors_wait(MPI_Request* mpi_request, GMEnv* env);

void gm_recv_vectors_wait(MPI_Request* mpi_request, GMEnv* env);

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

/*----------*/

void gm_set_vectors_start(GMVectors* vectors,
                          GMMirroredPointer* vectors_buf,
                          GMEnv* env);

void gm_set_vectors_wait(GMEnv* env);

void gm_get_metrics_start(GMMetrics* metrics,
                          GMMirroredPointer* metrics_buf,
                          GMEnv* env);

void gm_get_metrics_wait(GMEnv* env);

/*----------*/

void gm_vectors_to_buf(GMVectors* vectors,
                       GMMirroredPointer* vectors_buf,
                       GMEnv* env);

/*----------*/

void gm_compute_numerators_2way_start(GMVectors* vectors_left,
                                      GMVectors* vectors_right,
                                      GMMetrics* numerators,
                                      GMMirroredPointer* vectors_left_buf,
                                      GMMirroredPointer* vectors_right_buf,
                                      GMMirroredPointer* numerators_buf,
                                      int j_proc,
                                      _Bool compute_triang_only,
                                      GMEnv* env);

void gm_compute_2way_combine(GMMetrics* metrics,
                             GMMirroredPointer* metrics_buf,
                             GMVectorSums* vector_sums_left,
                             GMVectorSums* vector_sums_right,
                             int j_proc,
                             _Bool do_compute_triang_only,
                             GMEnv* env);

/*----------*/

void gm_compute_czekanowski_numerators_2way_start(
    GMVectors* vectors_left,
    GMVectors* vectors_right,
    GMMetrics* numerators,
    GMMirroredPointer* vectors_left_buf,
    GMMirroredPointer* vectors_right_buf,
    GMMirroredPointer* numerators_buf,
    int j_proc,
    _Bool do_compute_triang_only,
    GMEnv* env);

void gm_compute_czekanowski_2way_combine(
    GMMetrics* metrics,
    GMMirroredPointer* metrics_buf,
    GMFloat* __restrict__ vector_sums_left,
    GMFloat* __restrict__ vector_sums_right,
    int j_proc,
    _Bool compute_triang_only,
    GMEnv* env);

void gm_compute_czekanowski_numerators_3way_start(
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* numerators,
    GMMirroredPointer* vectors_i_buf,
    GMMirroredPointer* vectors_j_buf,
    GMMirroredPointer* vectors_k_buf,
    int j_proc,
    int k_proc,
    GMEnv* env);

void gm_compute_czekanowski_3way_combine(GMMetrics* metrics,
                                         GMFloat* __restrict__ vector_sums_i,
                                         GMFloat* __restrict__ vector_sums_j,
                                         GMFloat* __restrict__ vector_sums_k,
                                         int j_proc,
                                         int k_proc,
                                         GMEnv* env);

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_compute_metrics_utils_h---*/

/*---------------------------------------------------------------------------*/
