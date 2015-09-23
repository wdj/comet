/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproductblas_z.h normal z -> c, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_minproductBLAS_C_H
#define MAGMA_minproductBLAS_C_H

#include "magma_minproduct_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_minproductblas_ctranspose_inplace(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda );

void
magma_minproductblas_ctranspose(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr       dAT, magma_minproduct_int_t lddat );

void
magma_minproductblas_cgetmatrix_transpose(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dAT,   magma_minproduct_int_t ldda,
    magma_minproductFloatComplex          *hA,    magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr       dwork, magma_minproduct_int_t lddwork, magma_minproduct_int_t nb );

void
magma_minproductblas_csetmatrix_transpose(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const magma_minproductFloatComplex *hA,    magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr    dAT,   magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr    dwork, magma_minproduct_int_t lddwork, magma_minproduct_int_t nb );

  /*
   * RBT-related functions
   */
void
magma_minproductblas_cprbt(
    magma_minproduct_int_t n, 
    magma_minproductFloatComplex *dA, magma_minproduct_int_t ldda, 
    magma_minproductFloatComplex *du, magma_minproductFloatComplex *dv);

void
magma_minproductblas_cprbt_mv(
    magma_minproduct_int_t n, 
    magma_minproductFloatComplex *dv, magma_minproductFloatComplex *db);

void
magma_minproductblas_cprbt_mtv(
    magma_minproduct_int_t n, 
    magma_minproductFloatComplex *du, magma_minproductFloatComplex *db);

void
magma_minproductblas_caxpycp2(
    magma_minproduct_int_t m, magma_minproductFloatComplex *r, magma_minproductFloatComplex *x,
    const magma_minproductFloatComplex *b);

  /*
   * Multi-GPU copy functions
   */
void
magma_minproductblas_cgetmatrix_transpose_mgpu(
    magma_minproduct_int_t ngpu, magma_minproduct_queue_t queues[][2],
    magma_minproductFloatComplex_const_ptr const dAT[],    magma_minproduct_int_t ldda,
    magma_minproductFloatComplex                *hA,       magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr             dwork[],  magma_minproduct_int_t lddw,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb );

void
magma_minproductblas_csetmatrix_transpose_mgpu(
    magma_minproduct_int_t ngpu, magma_minproduct_queue_t queues[][2],
    const magma_minproductFloatComplex *hA,      magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr    dAT[],   magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr    dwork[], magma_minproduct_int_t lddw,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb );

void
magma_minproduct_cgetmatrix_1D_col_bcyclic(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr const dA[], magma_minproduct_int_t ldda,
    magma_minproductFloatComplex                *hA,   magma_minproduct_int_t lda,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb );

void
magma_minproduct_csetmatrix_1D_col_bcyclic(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const magma_minproductFloatComplex *hA,   magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr    dA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb );

void
magma_minproduct_cgetmatrix_1D_row_bcyclic(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr const dA[], magma_minproduct_int_t ldda,
    magma_minproductFloatComplex                *hA,   magma_minproduct_int_t lda,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb );

void
magma_minproduct_csetmatrix_1D_row_bcyclic(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const magma_minproductFloatComplex *hA,   magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr    dA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb );

// in src/chetrd_mgpu.cpp
// TODO rename csetmatrix_sy or similar
magma_minproduct_int_t
magma_minproduct_chtodhe(
    magma_minproduct_int_t ngpu, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductFloatComplex     *A,   magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr dA[], magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queues[][10],
    magma_minproduct_int_t *info );

// in src/cpotrf3_mgpu.cpp
// TODO same as magma_minproduct_chtodhe?
magma_minproduct_int_t
magma_minproduct_chtodpo(
    magma_minproduct_int_t ngpu, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t off_i, magma_minproduct_int_t off_j, magma_minproduct_int_t nb,
    magma_minproductFloatComplex     *A,   magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr dA[], magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queues[][3],
    magma_minproduct_int_t *info );

// in src/cpotrf3_mgpu.cpp
// TODO rename cgetmatrix_sy or similar
magma_minproduct_int_t
magma_minproduct_cdtohpo(
    magma_minproduct_int_t ngpu, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t off_i, magma_minproduct_int_t off_j, magma_minproduct_int_t nb, magma_minproduct_int_t NB,
    magma_minproductFloatComplex     *A,   magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr dA[], magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queues[][3],
    magma_minproduct_int_t *info );


  /*
   * Multi-GPU BLAS functions (alphabetical order)
   */
void
magma_minproductblas_chemm_mgpu(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dA[],    magma_minproduct_int_t ldda,  magma_minproduct_int_t offset,
    magma_minproductFloatComplex_ptr dB[],    magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dC[],    magma_minproduct_int_t lddc,
    magma_minproductFloatComplex_ptr dwork[], magma_minproduct_int_t dworksiz,
    magma_minproductFloatComplex    *C,       magma_minproduct_int_t ldc,
    magma_minproductFloatComplex    *work[],  magma_minproduct_int_t worksiz,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][20], magma_minproduct_int_t nbevents );

void
magma_minproductblas_chemm_mgpu_com(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dA[],    magma_minproduct_int_t ldda,  magma_minproduct_int_t offset,
    magma_minproductFloatComplex_ptr dB[],    magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dC[],    magma_minproduct_int_t lddc,
    magma_minproductFloatComplex_ptr dwork[], magma_minproduct_int_t dworksiz,
    magma_minproductFloatComplex    *C,       magma_minproduct_int_t ldc,
    magma_minproductFloatComplex    *work[],  magma_minproduct_int_t worksiz,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents,
    magma_minproduct_int_t gnode[Magma_minproductMaxGPUs][Magma_minproductMaxGPUs+2], magma_minproduct_int_t nbcmplx );

void
magma_minproductblas_chemm_mgpu_spec(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dA[],    magma_minproduct_int_t ldda,  magma_minproduct_int_t offset,
    magma_minproductFloatComplex_ptr dB[],    magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dC[],    magma_minproduct_int_t lddc,
    magma_minproductFloatComplex_ptr dwork[], magma_minproduct_int_t dworksiz,
    magma_minproductFloatComplex    *C,       magma_minproduct_int_t ldc,
    magma_minproductFloatComplex    *work[],  magma_minproduct_int_t worksiz,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents,
    magma_minproduct_int_t gnode[Magma_minproductMaxGPUs][Magma_minproductMaxGPUs+2], magma_minproduct_int_t nbcmplx );

void
magma_minproductblas_chemm_mgpu_spec33(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dA[],    magma_minproduct_int_t ldda,  magma_minproduct_int_t offset,
    magma_minproductFloatComplex_ptr dB[],    magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dC[],    magma_minproduct_int_t lddc,
    magma_minproductFloatComplex_ptr dVIN[],  magma_minproduct_int_t lddv, magma_minproduct_int_t v_offset,
    magma_minproductFloatComplex_ptr dwork[], magma_minproduct_int_t dworksiz,
    magma_minproductFloatComplex *C,       magma_minproduct_int_t ldc,
    magma_minproductFloatComplex *work[],  magma_minproduct_int_t worksiz,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents,
    magma_minproduct_int_t gnode[Magma_minproductMaxGPUs][Magma_minproductMaxGPUs+2], magma_minproduct_int_t nbcmplx );

magma_minproduct_int_t
magma_minproductblas_chemv_mgpu(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr const d_lA[], magma_minproduct_int_t ldda, magma_minproduct_int_t offset,
    magma_minproductFloatComplex const *x,         magma_minproduct_int_t incx,
    magma_minproductFloatComplex beta,             
    magma_minproductFloatComplex       *y,         magma_minproduct_int_t incy,
    magma_minproductFloatComplex       *hwork,     magma_minproduct_int_t lhwork,
    magma_minproductFloatComplex_ptr    dwork[],   magma_minproduct_int_t ldwork,
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[] );

magma_minproduct_int_t
magma_minproductblas_chemv_mgpu_sync(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr const d_lA[], magma_minproduct_int_t ldda, magma_minproduct_int_t offset,
    magma_minproductFloatComplex const *x,         magma_minproduct_int_t incx,
    magma_minproductFloatComplex beta,             
    magma_minproductFloatComplex       *y,         magma_minproduct_int_t incy,
    magma_minproductFloatComplex       *hwork,     magma_minproduct_int_t lhwork,
    magma_minproductFloatComplex_ptr    dwork[],   magma_minproduct_int_t ldwork,
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[] );

// Ichi's version, in src/chetrd_mgpu.cpp
void
magma_minproduct_cher2k_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t nb, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    float beta,
    magma_minproductFloatComplex_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t nqueue, magma_minproduct_queue_t queues[][10] );

void
magma_minproductblas_cher2k_mgpu_spec(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dA[], magma_minproduct_int_t ldda, magma_minproduct_int_t a_offset,
    magma_minproductFloatComplex_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    float beta,
    magma_minproductFloatComplex_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb, magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue );

void
magma_minproductblas_cher2k_mgpu_spec324(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dVIN[], magma_minproduct_int_t lddv, magma_minproduct_int_t v_offset,
    magma_minproductFloatComplex_ptr dWIN[], magma_minproduct_int_t lddw, magma_minproduct_int_t w_offset,
    float beta,
    magma_minproductFloatComplex_ptr dC[],   magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproductFloatComplex_ptr dwork[], magma_minproduct_int_t lndwork,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents );

void
magma_minproductblas_cher2k_mgpu_spec325(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dVIN[], magma_minproduct_int_t lddv, magma_minproduct_int_t v_offset,
    magma_minproductFloatComplex_ptr dWIN[], magma_minproduct_int_t lddw, magma_minproduct_int_t w_offset,
    float beta,
    magma_minproductFloatComplex_ptr dC[],   magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproductFloatComplex_ptr dwork[], magma_minproduct_int_t lndwork,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproductFloatComplex *harray[],
    magma_minproductFloatComplex_ptr *darray[],
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents );

void
magma_minproductblas_cher2k_mgpu2(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dA[], magma_minproduct_int_t ldda, magma_minproduct_int_t a_offset,
    magma_minproductFloatComplex_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    float beta,
    magma_minproductFloatComplex_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue );

// in src/cpotrf_mgpu_right.cpp
void
magma_minproduct_cherk_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t nb, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloatComplex_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    float beta,
    magma_minproductFloatComplex_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t nqueue, magma_minproduct_queue_t queues[][10]);

// in src/cpotrf_mgpu_right.cpp
void
magma_minproduct_cherk_mgpu2(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t nb, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloatComplex_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    float beta,
    magma_minproductFloatComplex_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t nqueue, magma_minproduct_queue_t queues[][10]);


  /*
   * LAPACK auxiliary functions (alphabetical order)
   */
void
magma_minproductblas_cgeadd(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr       dB, magma_minproduct_int_t lddb );

void
magma_minproductblas_clacpy(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr       dB, magma_minproduct_int_t lddb );

void
magma_minproductblas_clacpy_cnjg(
    magma_minproduct_int_t n, magma_minproductFloatComplex *dA1, magma_minproduct_int_t lda1,
    magma_minproductFloatComplex *dA2, magma_minproduct_int_t lda2);

float
magma_minproductblas_clange(
    magma_minproduct_norm_t norm,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dwork );

float
magma_minproductblas_clanhe(
    magma_minproduct_norm_t norm, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dwork );

float
magma_minproductblas_clansy(
    magma_minproduct_norm_t norm, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dwork );

void
magma_minproductblas_clarfg(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dalpha, magma_minproductFloatComplex_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr dtau );

void
magma_minproductblas_clarfg_work(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dalpha, magma_minproductFloatComplex_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr dtau, magma_minproductFloatComplex_ptr dwork );

void
magma_minproductblas_clascl(
    magma_minproduct_type_t type, magma_minproduct_int_t kl, magma_minproduct_int_t ku,
    float cfrom, float cto,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_clascl_2x2(
    magma_minproduct_type_t type, magma_minproduct_int_t m,
    magma_minproductFloatComplex *dW, magma_minproduct_int_t lddw,
    magma_minproductFloatComplex *dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_clascl2(
    magma_minproduct_type_t type,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dD,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_clascl_diag(
    magma_minproduct_type_t type, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dD, magma_minproduct_int_t lddd,
          magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_claset(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex offdiag, magma_minproductFloatComplex diag,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda );

void
magma_minproductblas_claset_band(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex offdiag, magma_minproductFloatComplex diag,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda);

void
magma_minproductblas_claswp(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dAT, magma_minproduct_int_t ldda,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci );

void
magma_minproductblas_claswp2(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dAT, magma_minproduct_int_t ldda,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    magma_minproductInt_const_ptr d_ipiv, magma_minproduct_int_t inci );

void
magma_minproductblas_claswpx(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldx, magma_minproduct_int_t ldy,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci );

void
magma_minproductblas_csymmetrize(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda );

void
magma_minproductblas_csymmetrize_tiles(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t ntile, magma_minproduct_int_t mstride, magma_minproduct_int_t nstride );

void
magma_minproductblas_ctrtri_diag(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr d_dinvA);

  /*
   * to cleanup (alphabetical order)
   */
void
magma_minproductblas_scnrm2_adjust(
    magma_minproduct_int_t k,
    magma_minproductFloat_ptr dxnorm,
    magma_minproductFloatComplex_ptr dc);

void
magma_minproductblas_scnrm2_check(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dxnorm,
    magma_minproductFloat_ptr dlsticc);

void
magma_minproductblas_scnrm2_cols(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dxnorm);

void
magma_minproductblas_scnrm2_row_check_adjust(
    magma_minproduct_int_t k, float tol,
    magma_minproductFloat_ptr dxnorm,
    magma_minproductFloat_ptr dxnorm2,
    magma_minproductFloatComplex_ptr dC, magma_minproduct_int_t lddc,
    magma_minproductFloat_ptr dlsticc);

void
magma_minproduct_clarfbx_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t k,
    magma_minproductFloatComplex_ptr V,  magma_minproduct_int_t ldv,
    magma_minproductFloatComplex_ptr dT, magma_minproduct_int_t ldt,
    magma_minproductFloatComplex_ptr c,
    magma_minproductFloatComplex_ptr dwork);

void
magma_minproduct_clarfg_gpu(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dx0,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex_ptr dtau,
    magma_minproductFloat_ptr        dxnorm,
    magma_minproductFloatComplex_ptr dAkk );

void
magma_minproduct_clarfgtx_gpu(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dx0,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex_ptr dtau,
    magma_minproductFloat_ptr        dxnorm,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t iter,
    magma_minproductFloatComplex_ptr V,  magma_minproduct_int_t ldv,
    magma_minproductFloatComplex_ptr T,  magma_minproduct_int_t ldt,
    magma_minproductFloatComplex_ptr dwork);

void
magma_minproduct_clarfgx_gpu(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dx0,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex_ptr dtau,
    magma_minproductFloat_ptr        dxnorm,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t iter);

void
magma_minproduct_clarfx_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr v,
    magma_minproductFloatComplex_ptr tau,
    magma_minproductFloatComplex_ptr C,  magma_minproduct_int_t ldc,
    magma_minproductFloat_ptr        xnorm,
    magma_minproductFloatComplex_ptr dT, magma_minproduct_int_t iter,
    magma_minproductFloatComplex_ptr work);


  /*
   * Level 1 BLAS (alphabetical order)
   */
void
magma_minproductblas_cswap(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr dy, magma_minproduct_int_t incy );

void
magma_minproductblas_cswapblk(
    magma_minproduct_order_t order,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t i1, magma_minproduct_int_t i2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci,
    magma_minproduct_int_t offset );

void
magma_minproductblas_cswapdblk(
    magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t inca,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb, magma_minproduct_int_t incb );

  /*
   * Level 2 BLAS (alphabetical order)
   */
void
magma_minproductblas_cgemv(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dy, magma_minproduct_int_t incy );


void
magma_minproductblas_cgemv_conjv(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dy, magma_minproduct_int_t incy);

magma_minproduct_int_t
magma_minproductblas_chemv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dy, magma_minproduct_int_t incy );

magma_minproduct_int_t
magma_minproductblas_csymv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dy, magma_minproduct_int_t incy );

magma_minproduct_int_t
magma_minproductblas_chemv_work(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dy, magma_minproduct_int_t incy,
    magma_minproductFloatComplex_ptr       dwork, magma_minproduct_int_t lwork,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproductblas_csymv_work(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dy, magma_minproduct_int_t incy,
    magma_minproductFloatComplex_ptr       dwork, magma_minproduct_int_t lwork,
    magma_minproduct_queue_t queue );

  /*
   * Level 3 BLAS (alphabetical order)
   */
void
magma_minproductblas_cgemm(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_cgemm_reduce(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_chemm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_csymm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_csyr2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_cher2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dB, magma_minproduct_int_t lddb,
    float  beta,
    magma_minproductFloatComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_csyrk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_cherk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    float  alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    float  beta,
    magma_minproductFloatComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_ctrsm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr       dB, magma_minproduct_int_t lddb );

void
magma_minproductblas_ctrsm_outofplace(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr       dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t flag, magma_minproductFloatComplex_ptr d_dinvA, magma_minproductFloatComplex_ptr dX );

void
magma_minproductblas_ctrsm_work(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr       dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t flag, magma_minproductFloatComplex_ptr d_dinvA, magma_minproductFloatComplex_ptr dX );


  /*
   * Wrappers for platform independence.
   * These wrap CUBLAS or AMD OpenCL BLAS functions.
   */

// ========================================
// copying vectors
// set copies host to device
// get copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

#define magma_minproduct_csetvector(           n, hx_src, incx, dy_dst, incy ) \
        magma_minproduct_csetvector_internal(  n, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_minproduct_cgetvector(           n, dx_src, incx, hy_dst, incy ) \
        magma_minproduct_cgetvector_internal(  n, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_minproduct_ccopyvector(          n, dx_src, incx, dy_dst, incy ) \
        magma_minproduct_ccopyvector_internal( n, dx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_minproduct_csetvector_async(           n, hx_src, incx, dy_dst, incy, queue ) \
        magma_minproduct_csetvector_async_internal(  n, hx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_cgetvector_async(           n, dx_src, incx, hy_dst, incy, queue ) \
        magma_minproduct_cgetvector_async_internal(  n, dx_src, incx, hy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_ccopyvector_async(          n, dx_src, incx, dy_dst, incy, queue ) \
        magma_minproduct_ccopyvector_async_internal( n, dx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

void
magma_minproduct_csetvector_internal(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex const    *hx_src, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr       dy_dst, magma_minproduct_int_t incy,
    const char* func, const char* file, int line );

void
magma_minproduct_cgetvector_internal(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dx_src, magma_minproduct_int_t incx,
    magma_minproductFloatComplex          *hy_dst, magma_minproduct_int_t incy,
    const char* func, const char* file, int line );

void
magma_minproduct_ccopyvector_internal(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dx_src, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr       dy_dst, magma_minproduct_int_t incy,
    const char* func, const char* file, int line );

void
magma_minproduct_csetvector_async_internal(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex const    *hx_src, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr       dy_dst, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void
magma_minproduct_cgetvector_async_internal(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dx_src, magma_minproduct_int_t incx,
    magma_minproductFloatComplex          *hy_dst, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void
magma_minproduct_ccopyvector_async_internal(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dx_src, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr       dy_dst, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// copying sub-matrices (contiguous columns)
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

#define magma_minproduct_csetmatrix(           m, n, hA_src, lda, dB_dst, lddb ) \
        magma_minproduct_csetmatrix_internal(  m, n, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_minproduct_cgetmatrix(           m, n, dA_src, ldda, hB_dst, ldb ) \
        magma_minproduct_cgetmatrix_internal(  m, n, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )

#define magma_minproduct_ccopymatrix(          m, n, dA_src, ldda, dB_dst, lddb ) \
        magma_minproduct_ccopymatrix_internal( m, n, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_minproduct_csetmatrix_async(           m, n, hA_src, lda, dB_dst, lddb, queue ) \
        magma_minproduct_csetmatrix_async_internal(  m, n, hA_src, lda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_cgetmatrix_async(           m, n, dA_src, ldda, hB_dst, ldb, queue ) \
        magma_minproduct_cgetmatrix_async_internal(  m, n, dA_src, ldda, hB_dst, ldb, queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_ccopymatrix_async(          m, n, dA_src, ldda, dB_dst, lddb, queue ) \
        magma_minproduct_ccopymatrix_async_internal( m, n, dA_src, ldda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

void
magma_minproduct_csetmatrix_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex const    *hA_src, magma_minproduct_int_t ldha,
    magma_minproductFloatComplex_ptr       dB_dst, magma_minproduct_int_t lddb,
    const char* func, const char* file, int line );

void
magma_minproduct_cgetmatrix_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA_src, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex          *hB_dst, magma_minproduct_int_t ldhb,
    const char* func, const char* file, int line );

void
magma_minproduct_ccopymatrix_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA_src, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr       dB_dst, magma_minproduct_int_t lddb,
    const char* func, const char* file, int line );

void
magma_minproduct_csetmatrix_async_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex const    *hA_src, magma_minproduct_int_t ldha,
    magma_minproductFloatComplex_ptr       dB_dst, magma_minproduct_int_t lddb,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void
magma_minproduct_cgetmatrix_async_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA_src, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex          *hB_dst, magma_minproduct_int_t ldhb,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void
magma_minproduct_ccopymatrix_async_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA_src, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr       dB_dst, magma_minproduct_int_t lddb,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// Level 1 BLAS (alphabetical order)

// in cublas_v2, result returned through output argument
magma_minproduct_int_t
magma_minproduct_icamax(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx );

// in cublas_v2, result returned through output argument
magma_minproduct_int_t
magma_minproduct_icamin(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx );

// in cublas_v2, result returned through output argument
float
magma_minproduct_scasum(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx );

void
magma_minproduct_caxpy(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr       dy, magma_minproduct_int_t incy );

void
magma_minproduct_ccopy(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr       dy, magma_minproduct_int_t incy );

// in cublas_v2, result returned through output argument
magma_minproductFloatComplex
magma_minproduct_cdotc(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_const_ptr dy, magma_minproduct_int_t incy );

// in cublas_v2, result returned through output argument
magma_minproductFloatComplex
magma_minproduct_cdotu(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_const_ptr dy, magma_minproduct_int_t incy );

// in cublas_v2, result returned through output argument
float
magma_minproduct_scnrm2(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx );

void
magma_minproduct_crot(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr dy, magma_minproduct_int_t incy,
    float dc, magma_minproductFloatComplex ds );

void
magma_minproduct_csrot(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr dy, magma_minproduct_int_t incy,
    float dc, float ds );

#ifdef REAL
void
magma_minproduct_crotm(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr dy, magma_minproduct_int_t incy,
    magma_minproductFloat_const_ptr param );

void
magma_minproduct_crotmg(
    magma_minproductFloat_ptr d1, magma_minproductFloat_ptr       d2,
    magma_minproductFloat_ptr x1, magma_minproductFloat_const_ptr y1,
    magma_minproductFloat_ptr param );
#endif

void
magma_minproduct_cscal(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dx, magma_minproduct_int_t incx );

void
magma_minproduct_csscal(
    magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloatComplex_ptr dx, magma_minproduct_int_t incx );

void
magma_minproduct_cswap(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr dy, magma_minproduct_int_t incy );

// ========================================
// Level 2 BLAS (alphabetical order)

void
magma_minproduct_cgemv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dy, magma_minproduct_int_t incy );

void
magma_minproduct_cgerc(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_const_ptr dy, magma_minproduct_int_t incy,
    magma_minproductFloatComplex_ptr       dA, magma_minproduct_int_t ldda );

void
magma_minproduct_cgeru(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_const_ptr dy, magma_minproduct_int_t incy,
    magma_minproductFloatComplex_ptr       dA, magma_minproduct_int_t ldda );

void
magma_minproduct_chemv(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dy, magma_minproduct_int_t incy );

void
magma_minproduct_cher(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr       dA, magma_minproduct_int_t ldda );

void
magma_minproduct_cher2(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_const_ptr dy, magma_minproduct_int_t incy,
    magma_minproductFloatComplex_ptr       dA, magma_minproduct_int_t ldda );

void
magma_minproduct_ctrmv(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr       dx, magma_minproduct_int_t incx );

void
magma_minproduct_ctrsv(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr       dx, magma_minproduct_int_t incx );

// ========================================
// Level 3 BLAS (alphabetical order)

void
magma_minproduct_cgemm(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_csymm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_chemm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_csyr2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_cher2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_const_ptr dB, magma_minproduct_int_t lddb,
    float beta,
    magma_minproductFloatComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_csyrk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_cherk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    float beta,
    magma_minproductFloatComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_ctrmm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr       dB, magma_minproduct_int_t lddb );

void
magma_minproduct_ctrsm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr       dB, magma_minproduct_int_t lddb );


#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif  /* MAGMA_minproductBLAS_C_H */
