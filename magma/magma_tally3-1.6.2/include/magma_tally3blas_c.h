/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally3blas_z.h normal z -> c, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally3BLAS_C_H
#define MAGMA_tally3BLAS_C_H

#include "magma_tally3_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_tally3blas_ctranspose_inplace(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda );

void
magma_tally3blas_ctranspose(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr       dAT, magma_tally3_int_t lddat );

void
magma_tally3blas_cgetmatrix_transpose(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dAT,   magma_tally3_int_t ldda,
    magma_tally3FloatComplex          *hA,    magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr       dwork, magma_tally3_int_t lddwork, magma_tally3_int_t nb );

void
magma_tally3blas_csetmatrix_transpose(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3FloatComplex *hA,    magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr    dAT,   magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr    dwork, magma_tally3_int_t lddwork, magma_tally3_int_t nb );

  /*
   * RBT-related functions
   */
void
magma_tally3blas_cprbt(
    magma_tally3_int_t n, 
    magma_tally3FloatComplex *dA, magma_tally3_int_t ldda, 
    magma_tally3FloatComplex *du, magma_tally3FloatComplex *dv);

void
magma_tally3blas_cprbt_mv(
    magma_tally3_int_t n, 
    magma_tally3FloatComplex *dv, magma_tally3FloatComplex *db);

void
magma_tally3blas_cprbt_mtv(
    magma_tally3_int_t n, 
    magma_tally3FloatComplex *du, magma_tally3FloatComplex *db);

void
magma_tally3blas_caxpycp2(
    magma_tally3_int_t m, magma_tally3FloatComplex *r, magma_tally3FloatComplex *x,
    const magma_tally3FloatComplex *b);

  /*
   * Multi-GPU copy functions
   */
void
magma_tally3blas_cgetmatrix_transpose_mgpu(
    magma_tally3_int_t ngpu, magma_tally3_queue_t queues[][2],
    magma_tally3FloatComplex_const_ptr const dAT[],    magma_tally3_int_t ldda,
    magma_tally3FloatComplex                *hA,       magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr             dwork[],  magma_tally3_int_t lddw,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb );

void
magma_tally3blas_csetmatrix_transpose_mgpu(
    magma_tally3_int_t ngpu, magma_tally3_queue_t queues[][2],
    const magma_tally3FloatComplex *hA,      magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr    dAT[],   magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr    dwork[], magma_tally3_int_t lddw,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb );

void
magma_tally3_cgetmatrix_1D_col_bcyclic(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr const dA[], magma_tally3_int_t ldda,
    magma_tally3FloatComplex                *hA,   magma_tally3_int_t lda,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb );

void
magma_tally3_csetmatrix_1D_col_bcyclic(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3FloatComplex *hA,   magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr    dA[], magma_tally3_int_t ldda,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb );

void
magma_tally3_cgetmatrix_1D_row_bcyclic(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr const dA[], magma_tally3_int_t ldda,
    magma_tally3FloatComplex                *hA,   magma_tally3_int_t lda,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb );

void
magma_tally3_csetmatrix_1D_row_bcyclic(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3FloatComplex *hA,   magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr    dA[], magma_tally3_int_t ldda,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb );

// in src/chetrd_mgpu.cpp
// TODO rename csetmatrix_sy or similar
magma_tally3_int_t
magma_tally3_chtodhe(
    magma_tally3_int_t ngpu, magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3FloatComplex     *A,   magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr dA[], magma_tally3_int_t ldda,
    magma_tally3_queue_t queues[][10],
    magma_tally3_int_t *info );

// in src/cpotrf3_mgpu.cpp
// TODO same as magma_tally3_chtodhe?
magma_tally3_int_t
magma_tally3_chtodpo(
    magma_tally3_int_t ngpu, magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t off_i, magma_tally3_int_t off_j, magma_tally3_int_t nb,
    magma_tally3FloatComplex     *A,   magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr dA[], magma_tally3_int_t ldda,
    magma_tally3_queue_t queues[][3],
    magma_tally3_int_t *info );

// in src/cpotrf3_mgpu.cpp
// TODO rename cgetmatrix_sy or similar
magma_tally3_int_t
magma_tally3_cdtohpo(
    magma_tally3_int_t ngpu, magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t off_i, magma_tally3_int_t off_j, magma_tally3_int_t nb, magma_tally3_int_t NB,
    magma_tally3FloatComplex     *A,   magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr dA[], magma_tally3_int_t ldda,
    magma_tally3_queue_t queues[][3],
    magma_tally3_int_t *info );


  /*
   * Multi-GPU BLAS functions (alphabetical order)
   */
void
magma_tally3blas_chemm_mgpu(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dA[],    magma_tally3_int_t ldda,  magma_tally3_int_t offset,
    magma_tally3FloatComplex_ptr dB[],    magma_tally3_int_t lddb,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dC[],    magma_tally3_int_t lddc,
    magma_tally3FloatComplex_ptr dwork[], magma_tally3_int_t dworksiz,
    magma_tally3FloatComplex    *C,       magma_tally3_int_t ldc,
    magma_tally3FloatComplex    *work[],  magma_tally3_int_t worksiz,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_event_t redevents[][20], magma_tally3_int_t nbevents );

void
magma_tally3blas_chemm_mgpu_com(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dA[],    magma_tally3_int_t ldda,  magma_tally3_int_t offset,
    magma_tally3FloatComplex_ptr dB[],    magma_tally3_int_t lddb,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dC[],    magma_tally3_int_t lddc,
    magma_tally3FloatComplex_ptr dwork[], magma_tally3_int_t dworksiz,
    magma_tally3FloatComplex    *C,       magma_tally3_int_t ldc,
    magma_tally3FloatComplex    *work[],  magma_tally3_int_t worksiz,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_event_t redevents[][Magma_tally3MaxGPUs*Magma_tally3MaxGPUs+10], magma_tally3_int_t nbevents,
    magma_tally3_int_t gnode[Magma_tally3MaxGPUs][Magma_tally3MaxGPUs+2], magma_tally3_int_t nbcmplx );

void
magma_tally3blas_chemm_mgpu_spec(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dA[],    magma_tally3_int_t ldda,  magma_tally3_int_t offset,
    magma_tally3FloatComplex_ptr dB[],    magma_tally3_int_t lddb,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dC[],    magma_tally3_int_t lddc,
    magma_tally3FloatComplex_ptr dwork[], magma_tally3_int_t dworksiz,
    magma_tally3FloatComplex    *C,       magma_tally3_int_t ldc,
    magma_tally3FloatComplex    *work[],  magma_tally3_int_t worksiz,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_event_t redevents[][Magma_tally3MaxGPUs*Magma_tally3MaxGPUs+10], magma_tally3_int_t nbevents,
    magma_tally3_int_t gnode[Magma_tally3MaxGPUs][Magma_tally3MaxGPUs+2], magma_tally3_int_t nbcmplx );

void
magma_tally3blas_chemm_mgpu_spec33(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dA[],    magma_tally3_int_t ldda,  magma_tally3_int_t offset,
    magma_tally3FloatComplex_ptr dB[],    magma_tally3_int_t lddb,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dC[],    magma_tally3_int_t lddc,
    magma_tally3FloatComplex_ptr dVIN[],  magma_tally3_int_t lddv, magma_tally3_int_t v_offset,
    magma_tally3FloatComplex_ptr dwork[], magma_tally3_int_t dworksiz,
    magma_tally3FloatComplex *C,       magma_tally3_int_t ldc,
    magma_tally3FloatComplex *work[],  magma_tally3_int_t worksiz,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_event_t redevents[][Magma_tally3MaxGPUs*Magma_tally3MaxGPUs+10], magma_tally3_int_t nbevents,
    magma_tally3_int_t gnode[Magma_tally3MaxGPUs][Magma_tally3MaxGPUs+2], magma_tally3_int_t nbcmplx );

magma_tally3_int_t
magma_tally3blas_chemv_mgpu(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr const d_lA[], magma_tally3_int_t ldda, magma_tally3_int_t offset,
    magma_tally3FloatComplex const *x,         magma_tally3_int_t incx,
    magma_tally3FloatComplex beta,             
    magma_tally3FloatComplex       *y,         magma_tally3_int_t incy,
    magma_tally3FloatComplex       *hwork,     magma_tally3_int_t lhwork,
    magma_tally3FloatComplex_ptr    dwork[],   magma_tally3_int_t ldwork,
    magma_tally3_int_t ngpu,
    magma_tally3_int_t nb,
    magma_tally3_queue_t queues[] );

magma_tally3_int_t
magma_tally3blas_chemv_mgpu_sync(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr const d_lA[], magma_tally3_int_t ldda, magma_tally3_int_t offset,
    magma_tally3FloatComplex const *x,         magma_tally3_int_t incx,
    magma_tally3FloatComplex beta,             
    magma_tally3FloatComplex       *y,         magma_tally3_int_t incy,
    magma_tally3FloatComplex       *hwork,     magma_tally3_int_t lhwork,
    magma_tally3FloatComplex_ptr    dwork[],   magma_tally3_int_t ldwork,
    magma_tally3_int_t ngpu,
    magma_tally3_int_t nb,
    magma_tally3_queue_t queues[] );

// Ichi's version, in src/chetrd_mgpu.cpp
void
magma_tally3_cher2k_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t nb, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dB[], magma_tally3_int_t lddb, magma_tally3_int_t b_offset,
    float beta,
    magma_tally3FloatComplex_ptr dC[], magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3_int_t nqueue, magma_tally3_queue_t queues[][10] );

void
magma_tally3blas_cher2k_mgpu_spec(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dA[], magma_tally3_int_t ldda, magma_tally3_int_t a_offset,
    magma_tally3FloatComplex_ptr dB[], magma_tally3_int_t lddb, magma_tally3_int_t b_offset,
    float beta,
    magma_tally3FloatComplex_ptr dC[], magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb, magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue );

void
magma_tally3blas_cher2k_mgpu_spec324(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dVIN[], magma_tally3_int_t lddv, magma_tally3_int_t v_offset,
    magma_tally3FloatComplex_ptr dWIN[], magma_tally3_int_t lddw, magma_tally3_int_t w_offset,
    float beta,
    magma_tally3FloatComplex_ptr dC[],   magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3FloatComplex_ptr dwork[], magma_tally3_int_t lndwork,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_event_t redevents[][Magma_tally3MaxGPUs*Magma_tally3MaxGPUs+10], magma_tally3_int_t nbevents );

void
magma_tally3blas_cher2k_mgpu_spec325(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dVIN[], magma_tally3_int_t lddv, magma_tally3_int_t v_offset,
    magma_tally3FloatComplex_ptr dWIN[], magma_tally3_int_t lddw, magma_tally3_int_t w_offset,
    float beta,
    magma_tally3FloatComplex_ptr dC[],   magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3FloatComplex_ptr dwork[], magma_tally3_int_t lndwork,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb,
    magma_tally3FloatComplex *harray[],
    magma_tally3FloatComplex_ptr *darray[],
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_event_t redevents[][Magma_tally3MaxGPUs*Magma_tally3MaxGPUs+10], magma_tally3_int_t nbevents );

void
magma_tally3blas_cher2k_mgpu2(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dA[], magma_tally3_int_t ldda, magma_tally3_int_t a_offset,
    magma_tally3FloatComplex_ptr dB[], magma_tally3_int_t lddb, magma_tally3_int_t b_offset,
    float beta,
    magma_tally3FloatComplex_ptr dC[], magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue );

// in src/cpotrf_mgpu_right.cpp
void
magma_tally3_cherk_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t nb, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3FloatComplex_ptr dB[], magma_tally3_int_t lddb, magma_tally3_int_t b_offset,
    float beta,
    magma_tally3FloatComplex_ptr dC[], magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3_int_t nqueue, magma_tally3_queue_t queues[][10]);

// in src/cpotrf_mgpu_right.cpp
void
magma_tally3_cherk_mgpu2(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t nb, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3FloatComplex_ptr dB[], magma_tally3_int_t lddb, magma_tally3_int_t b_offset,
    float beta,
    magma_tally3FloatComplex_ptr dC[], magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3_int_t nqueue, magma_tally3_queue_t queues[][10]);


  /*
   * LAPACK auxiliary functions (alphabetical order)
   */
void
magma_tally3blas_cgeadd(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr       dB, magma_tally3_int_t lddb );

void
magma_tally3blas_clacpy(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr       dB, magma_tally3_int_t lddb );

void
magma_tally3blas_clacpy_cnjg(
    magma_tally3_int_t n, magma_tally3FloatComplex *dA1, magma_tally3_int_t lda1,
    magma_tally3FloatComplex *dA2, magma_tally3_int_t lda2);

float
magma_tally3blas_clange(
    magma_tally3_norm_t norm,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dwork );

float
magma_tally3blas_clanhe(
    magma_tally3_norm_t norm, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dwork );

float
magma_tally3blas_clansy(
    magma_tally3_norm_t norm, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dwork );

void
magma_tally3blas_clarfg(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dalpha, magma_tally3FloatComplex_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr dtau );

void
magma_tally3blas_clarfg_work(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dalpha, magma_tally3FloatComplex_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr dtau, magma_tally3FloatComplex_ptr dwork );

void
magma_tally3blas_clascl(
    magma_tally3_type_t type, magma_tally3_int_t kl, magma_tally3_int_t ku,
    float cfrom, float cto,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info );

void
magma_tally3blas_clascl_2x2(
    magma_tally3_type_t type, magma_tally3_int_t m,
    magma_tally3FloatComplex *dW, magma_tally3_int_t lddw,
    magma_tally3FloatComplex *dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info );

void
magma_tally3blas_clascl2(
    magma_tally3_type_t type,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dD,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info );

void
magma_tally3blas_clascl_diag(
    magma_tally3_type_t type, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dD, magma_tally3_int_t lddd,
          magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info );

void
magma_tally3blas_claset(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex offdiag, magma_tally3FloatComplex diag,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda );

void
magma_tally3blas_claset_band(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex offdiag, magma_tally3FloatComplex diag,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda);

void
magma_tally3blas_claswp(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dAT, magma_tally3_int_t ldda,
    magma_tally3_int_t k1, magma_tally3_int_t k2,
    const magma_tally3_int_t *ipiv, magma_tally3_int_t inci );

void
magma_tally3blas_claswp2(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dAT, magma_tally3_int_t ldda,
    magma_tally3_int_t k1, magma_tally3_int_t k2,
    magma_tally3Int_const_ptr d_ipiv, magma_tally3_int_t inci );

void
magma_tally3blas_claswpx(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldx, magma_tally3_int_t ldy,
    magma_tally3_int_t k1, magma_tally3_int_t k2,
    const magma_tally3_int_t *ipiv, magma_tally3_int_t inci );

void
magma_tally3blas_csymmetrize(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda );

void
magma_tally3blas_csymmetrize_tiles(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t ntile, magma_tally3_int_t mstride, magma_tally3_int_t nstride );

void
magma_tally3blas_ctrtri_diag(
    magma_tally3_uplo_t uplo, magma_tally3_diag_t diag, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr d_dinvA);

  /*
   * to cleanup (alphabetical order)
   */
void
magma_tally3blas_scnrm2_adjust(
    magma_tally3_int_t k,
    magma_tally3Float_ptr dxnorm,
    magma_tally3FloatComplex_ptr dc);

void
magma_tally3blas_scnrm2_check(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dxnorm,
    magma_tally3Float_ptr dlsticc);

void
magma_tally3blas_scnrm2_cols(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dxnorm);

void
magma_tally3blas_scnrm2_row_check_adjust(
    magma_tally3_int_t k, float tol,
    magma_tally3Float_ptr dxnorm,
    magma_tally3Float_ptr dxnorm2,
    magma_tally3FloatComplex_ptr dC, magma_tally3_int_t lddc,
    magma_tally3Float_ptr dlsticc);

void
magma_tally3_clarfbx_gpu(
    magma_tally3_int_t m, magma_tally3_int_t k,
    magma_tally3FloatComplex_ptr V,  magma_tally3_int_t ldv,
    magma_tally3FloatComplex_ptr dT, magma_tally3_int_t ldt,
    magma_tally3FloatComplex_ptr c,
    magma_tally3FloatComplex_ptr dwork);

void
magma_tally3_clarfg_gpu(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dx0,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex_ptr dtau,
    magma_tally3Float_ptr        dxnorm,
    magma_tally3FloatComplex_ptr dAkk );

void
magma_tally3_clarfgtx_gpu(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dx0,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex_ptr dtau,
    magma_tally3Float_ptr        dxnorm,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t iter,
    magma_tally3FloatComplex_ptr V,  magma_tally3_int_t ldv,
    magma_tally3FloatComplex_ptr T,  magma_tally3_int_t ldt,
    magma_tally3FloatComplex_ptr dwork);

void
magma_tally3_clarfgx_gpu(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dx0,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex_ptr dtau,
    magma_tally3Float_ptr        dxnorm,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t iter);

void
magma_tally3_clarfx_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr v,
    magma_tally3FloatComplex_ptr tau,
    magma_tally3FloatComplex_ptr C,  magma_tally3_int_t ldc,
    magma_tally3Float_ptr        xnorm,
    magma_tally3FloatComplex_ptr dT, magma_tally3_int_t iter,
    magma_tally3FloatComplex_ptr work);


  /*
   * Level 1 BLAS (alphabetical order)
   */
void
magma_tally3blas_cswap(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr dy, magma_tally3_int_t incy );

void
magma_tally3blas_cswapblk(
    magma_tally3_order_t order,
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t i1, magma_tally3_int_t i2,
    const magma_tally3_int_t *ipiv, magma_tally3_int_t inci,
    magma_tally3_int_t offset );

void
magma_tally3blas_cswapdblk(
    magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t inca,
    magma_tally3FloatComplex_ptr dB, magma_tally3_int_t lddb, magma_tally3_int_t incb );

  /*
   * Level 2 BLAS (alphabetical order)
   */
void
magma_tally3blas_cgemv(
    magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dy, magma_tally3_int_t incy );


void
magma_tally3blas_cgemv_conjv(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dy, magma_tally3_int_t incy);

magma_tally3_int_t
magma_tally3blas_chemv(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dy, magma_tally3_int_t incy );

magma_tally3_int_t
magma_tally3blas_csymv(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dy, magma_tally3_int_t incy );

magma_tally3_int_t
magma_tally3blas_chemv_work(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dy, magma_tally3_int_t incy,
    magma_tally3FloatComplex_ptr       dwork, magma_tally3_int_t lwork,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3blas_csymv_work(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dy, magma_tally3_int_t incy,
    magma_tally3FloatComplex_ptr       dwork, magma_tally3_int_t lwork,
    magma_tally3_queue_t queue );

  /*
   * Level 3 BLAS (alphabetical order)
   */
void
magma_tally3blas_cgemm(
    magma_tally3_trans_t transA, magma_tally3_trans_t transB,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dB, magma_tally3_int_t lddb,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3blas_cgemm_reduce(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dB, magma_tally3_int_t lddb,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3blas_chemm(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dB, magma_tally3_int_t lddb,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3blas_csymm(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dB, magma_tally3_int_t lddb,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3blas_csyr2k(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dB, magma_tally3_int_t lddb,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3blas_cher2k(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dB, magma_tally3_int_t lddb,
    float  beta,
    magma_tally3FloatComplex_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3blas_csyrk(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3blas_cherk(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t n, magma_tally3_int_t k,
    float  alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    float  beta,
    magma_tally3FloatComplex_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3blas_ctrsm(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transA, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr       dB, magma_tally3_int_t lddb );

void
magma_tally3blas_ctrsm_outofplace(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transA, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr       dB, magma_tally3_int_t lddb,
    magma_tally3_int_t flag, magma_tally3FloatComplex_ptr d_dinvA, magma_tally3FloatComplex_ptr dX );

void
magma_tally3blas_ctrsm_work(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transA, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr       dB, magma_tally3_int_t lddb,
    magma_tally3_int_t flag, magma_tally3FloatComplex_ptr d_dinvA, magma_tally3FloatComplex_ptr dX );


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

#define magma_tally3_csetvector(           n, hx_src, incx, dy_dst, incy ) \
        magma_tally3_csetvector_internal(  n, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally3_cgetvector(           n, dx_src, incx, hy_dst, incy ) \
        magma_tally3_cgetvector_internal(  n, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally3_ccopyvector(          n, dx_src, incx, dy_dst, incy ) \
        magma_tally3_ccopyvector_internal( n, dx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally3_csetvector_async(           n, hx_src, incx, dy_dst, incy, queue ) \
        magma_tally3_csetvector_async_internal(  n, hx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_tally3_cgetvector_async(           n, dx_src, incx, hy_dst, incy, queue ) \
        magma_tally3_cgetvector_async_internal(  n, dx_src, incx, hy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_tally3_ccopyvector_async(          n, dx_src, incx, dy_dst, incy, queue ) \
        magma_tally3_ccopyvector_async_internal( n, dx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

void
magma_tally3_csetvector_internal(
    magma_tally3_int_t n,
    magma_tally3FloatComplex const    *hx_src, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr       dy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally3_cgetvector_internal(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dx_src, magma_tally3_int_t incx,
    magma_tally3FloatComplex          *hy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally3_ccopyvector_internal(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dx_src, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr       dy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally3_csetvector_async_internal(
    magma_tally3_int_t n,
    magma_tally3FloatComplex const    *hx_src, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr       dy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally3_cgetvector_async_internal(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dx_src, magma_tally3_int_t incx,
    magma_tally3FloatComplex          *hy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally3_ccopyvector_async_internal(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dx_src, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr       dy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// copying sub-matrices (contiguous columns)
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

#define magma_tally3_csetmatrix(           m, n, hA_src, lda, dB_dst, lddb ) \
        magma_tally3_csetmatrix_internal(  m, n, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_tally3_cgetmatrix(           m, n, dA_src, ldda, hB_dst, ldb ) \
        magma_tally3_cgetmatrix_internal(  m, n, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )

#define magma_tally3_ccopymatrix(          m, n, dA_src, ldda, dB_dst, lddb ) \
        magma_tally3_ccopymatrix_internal( m, n, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_tally3_csetmatrix_async(           m, n, hA_src, lda, dB_dst, lddb, queue ) \
        magma_tally3_csetmatrix_async_internal(  m, n, hA_src, lda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

#define magma_tally3_cgetmatrix_async(           m, n, dA_src, ldda, hB_dst, ldb, queue ) \
        magma_tally3_cgetmatrix_async_internal(  m, n, dA_src, ldda, hB_dst, ldb, queue, __func__, __FILE__, __LINE__ )

#define magma_tally3_ccopymatrix_async(          m, n, dA_src, ldda, dB_dst, lddb, queue ) \
        magma_tally3_ccopymatrix_async_internal( m, n, dA_src, ldda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

void
magma_tally3_csetmatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex const    *hA_src, magma_tally3_int_t ldha,
    magma_tally3FloatComplex_ptr       dB_dst, magma_tally3_int_t lddb,
    const char* func, const char* file, int line );

void
magma_tally3_cgetmatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA_src, magma_tally3_int_t ldda,
    magma_tally3FloatComplex          *hB_dst, magma_tally3_int_t ldhb,
    const char* func, const char* file, int line );

void
magma_tally3_ccopymatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA_src, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr       dB_dst, magma_tally3_int_t lddb,
    const char* func, const char* file, int line );

void
magma_tally3_csetmatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex const    *hA_src, magma_tally3_int_t ldha,
    magma_tally3FloatComplex_ptr       dB_dst, magma_tally3_int_t lddb,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally3_cgetmatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA_src, magma_tally3_int_t ldda,
    magma_tally3FloatComplex          *hB_dst, magma_tally3_int_t ldhb,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally3_ccopymatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA_src, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr       dB_dst, magma_tally3_int_t lddb,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// Level 1 BLAS (alphabetical order)

// in cublas_v2, result returned through output argument
magma_tally3_int_t
magma_tally3_icamax(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx );

// in cublas_v2, result returned through output argument
magma_tally3_int_t
magma_tally3_icamin(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx );

// in cublas_v2, result returned through output argument
float
magma_tally3_scasum(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx );

void
magma_tally3_caxpy(
    magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr       dy, magma_tally3_int_t incy );

void
magma_tally3_ccopy(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr       dy, magma_tally3_int_t incy );

// in cublas_v2, result returned through output argument
magma_tally3FloatComplex
magma_tally3_cdotc(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex_const_ptr dy, magma_tally3_int_t incy );

// in cublas_v2, result returned through output argument
magma_tally3FloatComplex
magma_tally3_cdotu(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex_const_ptr dy, magma_tally3_int_t incy );

// in cublas_v2, result returned through output argument
float
magma_tally3_scnrm2(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx );

void
magma_tally3_crot(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr dy, magma_tally3_int_t incy,
    float dc, magma_tally3FloatComplex ds );

void
magma_tally3_csrot(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr dy, magma_tally3_int_t incy,
    float dc, float ds );

#ifdef REAL
void
magma_tally3_crotm(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_ptr dy, magma_tally3_int_t incy,
    magma_tally3Float_const_ptr param );

void
magma_tally3_crotmg(
    magma_tally3Float_ptr d1, magma_tally3Float_ptr       d2,
    magma_tally3Float_ptr x1, magma_tally3Float_const_ptr y1,
    magma_tally3Float_ptr param );
#endif

void
magma_tally3_cscal(
    magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dx, magma_tally3_int_t incx );

void
magma_tally3_csscal(
    magma_tally3_int_t n,
    float alpha,
    magma_tally3FloatComplex_ptr dx, magma_tally3_int_t incx );

void
magma_tally3_cswap(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr dy, magma_tally3_int_t incy );

// ========================================
// Level 2 BLAS (alphabetical order)

void
magma_tally3_cgemv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dy, magma_tally3_int_t incy );

void
magma_tally3_cgerc(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex_const_ptr dy, magma_tally3_int_t incy,
    magma_tally3FloatComplex_ptr       dA, magma_tally3_int_t ldda );

void
magma_tally3_cgeru(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex_const_ptr dy, magma_tally3_int_t incy,
    magma_tally3FloatComplex_ptr       dA, magma_tally3_int_t ldda );

void
magma_tally3_chemv(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dy, magma_tally3_int_t incy );

void
magma_tally3_cher(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    float alpha,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex_ptr       dA, magma_tally3_int_t ldda );

void
magma_tally3_cher2(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3FloatComplex_const_ptr dy, magma_tally3_int_t incy,
    magma_tally3FloatComplex_ptr       dA, magma_tally3_int_t ldda );

void
magma_tally3_ctrmv(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_diag_t diag,
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr       dx, magma_tally3_int_t incx );

void
magma_tally3_ctrsv(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_diag_t diag,
    magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr       dx, magma_tally3_int_t incx );

// ========================================
// Level 3 BLAS (alphabetical order)

void
magma_tally3_cgemm(
    magma_tally3_trans_t transA, magma_tally3_trans_t transB,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dB, magma_tally3_int_t lddb,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3_csymm(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dB, magma_tally3_int_t lddb,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3_chemm(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dB, magma_tally3_int_t lddb,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3_csyr2k(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dB, magma_tally3_int_t lddb,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3_cher2k(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_const_ptr dB, magma_tally3_int_t lddb,
    float beta,
    magma_tally3FloatComplex_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3_csyrk(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3_cherk(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    float beta,
    magma_tally3FloatComplex_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3_ctrmm(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr       dB, magma_tally3_int_t lddb );

void
magma_tally3_ctrsm(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr       dB, magma_tally3_int_t lddb );


#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif  /* MAGMA_tally3BLAS_C_H */
