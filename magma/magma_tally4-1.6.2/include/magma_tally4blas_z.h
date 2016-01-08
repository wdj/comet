/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
*/

#ifndef MAGMA_tally4BLAS_Z_H
#define MAGMA_tally4BLAS_Z_H

#include "magma_tally4_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_tally4blas_ztranspose_inplace(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda );

void
magma_tally4blas_ztranspose(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr       dAT, magma_tally4_int_t lddat );

void
magma_tally4blas_zgetmatrix_transpose(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dAT,   magma_tally4_int_t ldda,
    magma_tally4DoubleComplex          *hA,    magma_tally4_int_t lda,
    magma_tally4DoubleComplex_ptr       dwork, magma_tally4_int_t lddwork, magma_tally4_int_t nb );

void
magma_tally4blas_zsetmatrix_transpose(
    magma_tally4_int_t m, magma_tally4_int_t n,
    const magma_tally4DoubleComplex *hA,    magma_tally4_int_t lda,
    magma_tally4DoubleComplex_ptr    dAT,   magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr    dwork, magma_tally4_int_t lddwork, magma_tally4_int_t nb );

  /*
   * RBT-related functions
   */
void
magma_tally4blas_zprbt(
    magma_tally4_int_t n, 
    magma_tally4DoubleComplex *dA, magma_tally4_int_t ldda, 
    magma_tally4DoubleComplex *du, magma_tally4DoubleComplex *dv);

void
magma_tally4blas_zprbt_mv(
    magma_tally4_int_t n, 
    magma_tally4DoubleComplex *dv, magma_tally4DoubleComplex *db);

void
magma_tally4blas_zprbt_mtv(
    magma_tally4_int_t n, 
    magma_tally4DoubleComplex *du, magma_tally4DoubleComplex *db);

void
magma_tally4blas_zaxpycp2(
    magma_tally4_int_t m, magma_tally4DoubleComplex *r, magma_tally4DoubleComplex *x,
    const magma_tally4DoubleComplex *b);

  /*
   * Multi-GPU copy functions
   */
void
magma_tally4blas_zgetmatrix_transpose_mgpu(
    magma_tally4_int_t ngpu, magma_tally4_queue_t queues[][2],
    magma_tally4DoubleComplex_const_ptr const dAT[],    magma_tally4_int_t ldda,
    magma_tally4DoubleComplex                *hA,       magma_tally4_int_t lda,
    magma_tally4DoubleComplex_ptr             dwork[],  magma_tally4_int_t lddw,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nb );

void
magma_tally4blas_zsetmatrix_transpose_mgpu(
    magma_tally4_int_t ngpu, magma_tally4_queue_t queues[][2],
    const magma_tally4DoubleComplex *hA,      magma_tally4_int_t lda,
    magma_tally4DoubleComplex_ptr    dAT[],   magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr    dwork[], magma_tally4_int_t lddw,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nb );

void
magma_tally4_zgetmatrix_1D_col_bcyclic(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr const dA[], magma_tally4_int_t ldda,
    magma_tally4DoubleComplex                *hA,   magma_tally4_int_t lda,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb );

void
magma_tally4_zsetmatrix_1D_col_bcyclic(
    magma_tally4_int_t m, magma_tally4_int_t n,
    const magma_tally4DoubleComplex *hA,   magma_tally4_int_t lda,
    magma_tally4DoubleComplex_ptr    dA[], magma_tally4_int_t ldda,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb );

void
magma_tally4_zgetmatrix_1D_row_bcyclic(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr const dA[], magma_tally4_int_t ldda,
    magma_tally4DoubleComplex                *hA,   magma_tally4_int_t lda,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb );

void
magma_tally4_zsetmatrix_1D_row_bcyclic(
    magma_tally4_int_t m, magma_tally4_int_t n,
    const magma_tally4DoubleComplex *hA,   magma_tally4_int_t lda,
    magma_tally4DoubleComplex_ptr    dA[], magma_tally4_int_t ldda,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb );

// in src/zhetrd_mgpu.cpp
// TODO rename zsetmatrix_sy or similar
magma_tally4_int_t
magma_tally4_zhtodhe(
    magma_tally4_int_t ngpu, magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4DoubleComplex     *A,   magma_tally4_int_t lda,
    magma_tally4DoubleComplex_ptr dA[], magma_tally4_int_t ldda,
    magma_tally4_queue_t queues[][10],
    magma_tally4_int_t *info );

// in src/zpotrf3_mgpu.cpp
// TODO same as magma_tally4_zhtodhe?
magma_tally4_int_t
magma_tally4_zhtodpo(
    magma_tally4_int_t ngpu, magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t off_i, magma_tally4_int_t off_j, magma_tally4_int_t nb,
    magma_tally4DoubleComplex     *A,   magma_tally4_int_t lda,
    magma_tally4DoubleComplex_ptr dA[], magma_tally4_int_t ldda,
    magma_tally4_queue_t queues[][3],
    magma_tally4_int_t *info );

// in src/zpotrf3_mgpu.cpp
// TODO rename zgetmatrix_sy or similar
magma_tally4_int_t
magma_tally4_zdtohpo(
    magma_tally4_int_t ngpu, magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t off_i, magma_tally4_int_t off_j, magma_tally4_int_t nb, magma_tally4_int_t NB,
    magma_tally4DoubleComplex     *A,   magma_tally4_int_t lda,
    magma_tally4DoubleComplex_ptr dA[], magma_tally4_int_t ldda,
    magma_tally4_queue_t queues[][3],
    magma_tally4_int_t *info );


  /*
   * Multi-GPU BLAS functions (alphabetical order)
   */
void
magma_tally4blas_zhemm_mgpu(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dA[],    magma_tally4_int_t ldda,  magma_tally4_int_t offset,
    magma_tally4DoubleComplex_ptr dB[],    magma_tally4_int_t lddb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dC[],    magma_tally4_int_t lddc,
    magma_tally4DoubleComplex_ptr dwork[], magma_tally4_int_t dworksiz,
    magma_tally4DoubleComplex    *C,       magma_tally4_int_t ldc,
    magma_tally4DoubleComplex    *work[],  magma_tally4_int_t worksiz,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_event_t redevents[][20], magma_tally4_int_t nbevents );

void
magma_tally4blas_zhemm_mgpu_com(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dA[],    magma_tally4_int_t ldda,  magma_tally4_int_t offset,
    magma_tally4DoubleComplex_ptr dB[],    magma_tally4_int_t lddb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dC[],    magma_tally4_int_t lddc,
    magma_tally4DoubleComplex_ptr dwork[], magma_tally4_int_t dworksiz,
    magma_tally4DoubleComplex    *C,       magma_tally4_int_t ldc,
    magma_tally4DoubleComplex    *work[],  magma_tally4_int_t worksiz,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_event_t redevents[][Magma_tally4MaxGPUs*Magma_tally4MaxGPUs+10], magma_tally4_int_t nbevents,
    magma_tally4_int_t gnode[Magma_tally4MaxGPUs][Magma_tally4MaxGPUs+2], magma_tally4_int_t nbcmplx );

void
magma_tally4blas_zhemm_mgpu_spec(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dA[],    magma_tally4_int_t ldda,  magma_tally4_int_t offset,
    magma_tally4DoubleComplex_ptr dB[],    magma_tally4_int_t lddb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dC[],    magma_tally4_int_t lddc,
    magma_tally4DoubleComplex_ptr dwork[], magma_tally4_int_t dworksiz,
    magma_tally4DoubleComplex    *C,       magma_tally4_int_t ldc,
    magma_tally4DoubleComplex    *work[],  magma_tally4_int_t worksiz,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_event_t redevents[][Magma_tally4MaxGPUs*Magma_tally4MaxGPUs+10], magma_tally4_int_t nbevents,
    magma_tally4_int_t gnode[Magma_tally4MaxGPUs][Magma_tally4MaxGPUs+2], magma_tally4_int_t nbcmplx );

void
magma_tally4blas_zhemm_mgpu_spec33(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dA[],    magma_tally4_int_t ldda,  magma_tally4_int_t offset,
    magma_tally4DoubleComplex_ptr dB[],    magma_tally4_int_t lddb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dC[],    magma_tally4_int_t lddc,
    magma_tally4DoubleComplex_ptr dVIN[],  magma_tally4_int_t lddv, magma_tally4_int_t v_offset,
    magma_tally4DoubleComplex_ptr dwork[], magma_tally4_int_t dworksiz,
    magma_tally4DoubleComplex *C,       magma_tally4_int_t ldc,
    magma_tally4DoubleComplex *work[],  magma_tally4_int_t worksiz,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_event_t redevents[][Magma_tally4MaxGPUs*Magma_tally4MaxGPUs+10], magma_tally4_int_t nbevents,
    magma_tally4_int_t gnode[Magma_tally4MaxGPUs][Magma_tally4MaxGPUs+2], magma_tally4_int_t nbcmplx );

magma_tally4_int_t
magma_tally4blas_zhemv_mgpu(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr const d_lA[], magma_tally4_int_t ldda, magma_tally4_int_t offset,
    magma_tally4DoubleComplex const *x,         magma_tally4_int_t incx,
    magma_tally4DoubleComplex beta,             
    magma_tally4DoubleComplex       *y,         magma_tally4_int_t incy,
    magma_tally4DoubleComplex       *hwork,     magma_tally4_int_t lhwork,
    magma_tally4DoubleComplex_ptr    dwork[],   magma_tally4_int_t ldwork,
    magma_tally4_int_t ngpu,
    magma_tally4_int_t nb,
    magma_tally4_queue_t queues[] );

magma_tally4_int_t
magma_tally4blas_zhemv_mgpu_sync(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr const d_lA[], magma_tally4_int_t ldda, magma_tally4_int_t offset,
    magma_tally4DoubleComplex const *x,         magma_tally4_int_t incx,
    magma_tally4DoubleComplex beta,             
    magma_tally4DoubleComplex       *y,         magma_tally4_int_t incy,
    magma_tally4DoubleComplex       *hwork,     magma_tally4_int_t lhwork,
    magma_tally4DoubleComplex_ptr    dwork[],   magma_tally4_int_t ldwork,
    magma_tally4_int_t ngpu,
    magma_tally4_int_t nb,
    magma_tally4_queue_t queues[] );

// Ichi's version, in src/zhetrd_mgpu.cpp
void
magma_tally4_zher2k_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t nb, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dB[], magma_tally4_int_t lddb, magma_tally4_int_t b_offset,
    double beta,
    magma_tally4DoubleComplex_ptr dC[], magma_tally4_int_t lddc, magma_tally4_int_t c_offset,
    magma_tally4_int_t nqueue, magma_tally4_queue_t queues[][10] );

void
magma_tally4blas_zher2k_mgpu_spec(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dA[], magma_tally4_int_t ldda, magma_tally4_int_t a_offset,
    magma_tally4DoubleComplex_ptr dB[], magma_tally4_int_t lddb, magma_tally4_int_t b_offset,
    double beta,
    magma_tally4DoubleComplex_ptr dC[], magma_tally4_int_t lddc, magma_tally4_int_t c_offset,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb, magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue );

void
magma_tally4blas_zher2k_mgpu_spec324(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dVIN[], magma_tally4_int_t lddv, magma_tally4_int_t v_offset,
    magma_tally4DoubleComplex_ptr dWIN[], magma_tally4_int_t lddw, magma_tally4_int_t w_offset,
    double beta,
    magma_tally4DoubleComplex_ptr dC[],   magma_tally4_int_t lddc, magma_tally4_int_t c_offset,
    magma_tally4DoubleComplex_ptr dwork[], magma_tally4_int_t lndwork,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_event_t redevents[][Magma_tally4MaxGPUs*Magma_tally4MaxGPUs+10], magma_tally4_int_t nbevents );

void
magma_tally4blas_zher2k_mgpu_spec325(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dVIN[], magma_tally4_int_t lddv, magma_tally4_int_t v_offset,
    magma_tally4DoubleComplex_ptr dWIN[], magma_tally4_int_t lddw, magma_tally4_int_t w_offset,
    double beta,
    magma_tally4DoubleComplex_ptr dC[],   magma_tally4_int_t lddc, magma_tally4_int_t c_offset,
    magma_tally4DoubleComplex_ptr dwork[], magma_tally4_int_t lndwork,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb,
    magma_tally4DoubleComplex *harray[],
    magma_tally4DoubleComplex_ptr *darray[],
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_event_t redevents[][Magma_tally4MaxGPUs*Magma_tally4MaxGPUs+10], magma_tally4_int_t nbevents );

void
magma_tally4blas_zher2k_mgpu2(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dA[], magma_tally4_int_t ldda, magma_tally4_int_t a_offset,
    magma_tally4DoubleComplex_ptr dB[], magma_tally4_int_t lddb, magma_tally4_int_t b_offset,
    double beta,
    magma_tally4DoubleComplex_ptr dC[], magma_tally4_int_t lddc, magma_tally4_int_t c_offset,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue );

// in src/zpotrf_mgpu_right.cpp
void
magma_tally4_zherk_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t nb, magma_tally4_int_t n, magma_tally4_int_t k,
    double alpha,
    magma_tally4DoubleComplex_ptr dB[], magma_tally4_int_t lddb, magma_tally4_int_t b_offset,
    double beta,
    magma_tally4DoubleComplex_ptr dC[], magma_tally4_int_t lddc, magma_tally4_int_t c_offset,
    magma_tally4_int_t nqueue, magma_tally4_queue_t queues[][10]);

// in src/zpotrf_mgpu_right.cpp
void
magma_tally4_zherk_mgpu2(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t nb, magma_tally4_int_t n, magma_tally4_int_t k,
    double alpha,
    magma_tally4DoubleComplex_ptr dB[], magma_tally4_int_t lddb, magma_tally4_int_t b_offset,
    double beta,
    magma_tally4DoubleComplex_ptr dC[], magma_tally4_int_t lddc, magma_tally4_int_t c_offset,
    magma_tally4_int_t nqueue, magma_tally4_queue_t queues[][10]);


  /*
   * LAPACK auxiliary functions (alphabetical order)
   */
void
magma_tally4blas_zgeadd(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr       dB, magma_tally4_int_t lddb );

void
magma_tally4blas_zlacpy(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr       dB, magma_tally4_int_t lddb );

void
magma_tally4blas_zlacpy_cnjg(
    magma_tally4_int_t n, magma_tally4DoubleComplex *dA1, magma_tally4_int_t lda1,
    magma_tally4DoubleComplex *dA2, magma_tally4_int_t lda2);

double
magma_tally4blas_zlange(
    magma_tally4_norm_t norm,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dwork );

double
magma_tally4blas_zlanhe(
    magma_tally4_norm_t norm, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dwork );

double
magma_tally4blas_zlansy(
    magma_tally4_norm_t norm, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dwork );

void
magma_tally4blas_zlarfg(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dalpha, magma_tally4DoubleComplex_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_ptr dtau );

void
magma_tally4blas_zlarfg_work(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dalpha, magma_tally4DoubleComplex_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_ptr dtau, magma_tally4DoubleComplex_ptr dwork );

void
magma_tally4blas_zlascl(
    magma_tally4_type_t type, magma_tally4_int_t kl, magma_tally4_int_t ku,
    double cfrom, double cto,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info );

void
magma_tally4blas_zlascl_2x2(
    magma_tally4_type_t type, magma_tally4_int_t m,
    magma_tally4DoubleComplex *dW, magma_tally4_int_t lddw,
    magma_tally4DoubleComplex *dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info );

void
magma_tally4blas_zlascl2(
    magma_tally4_type_t type,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_const_ptr dD,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info );

void
magma_tally4blas_zlascl_diag(
    magma_tally4_type_t type, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dD, magma_tally4_int_t lddd,
          magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info );

void
magma_tally4blas_zlaset(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex offdiag, magma_tally4DoubleComplex diag,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda );

void
magma_tally4blas_zlaset_band(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex offdiag, magma_tally4DoubleComplex diag,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda);

void
magma_tally4blas_zlaswp(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dAT, magma_tally4_int_t ldda,
    magma_tally4_int_t k1, magma_tally4_int_t k2,
    const magma_tally4_int_t *ipiv, magma_tally4_int_t inci );

void
magma_tally4blas_zlaswp2(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dAT, magma_tally4_int_t ldda,
    magma_tally4_int_t k1, magma_tally4_int_t k2,
    magma_tally4Int_const_ptr d_ipiv, magma_tally4_int_t inci );

void
magma_tally4blas_zlaswpx(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldx, magma_tally4_int_t ldy,
    magma_tally4_int_t k1, magma_tally4_int_t k2,
    const magma_tally4_int_t *ipiv, magma_tally4_int_t inci );

void
magma_tally4blas_zsymmetrize(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda );

void
magma_tally4blas_zsymmetrize_tiles(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t ntile, magma_tally4_int_t mstride, magma_tally4_int_t nstride );

void
magma_tally4blas_ztrtri_diag(
    magma_tally4_uplo_t uplo, magma_tally4_diag_t diag, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr d_dinvA);

  /*
   * to cleanup (alphabetical order)
   */
void
magma_tally4blas_dznrm2_adjust(
    magma_tally4_int_t k,
    magma_tally4Double_ptr dxnorm,
    magma_tally4DoubleComplex_ptr dc);

void
magma_tally4blas_dznrm2_check(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dxnorm,
    magma_tally4Double_ptr dlsticc);

void
magma_tally4blas_dznrm2_cols(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dxnorm);

void
magma_tally4blas_dznrm2_row_check_adjust(
    magma_tally4_int_t k, double tol,
    magma_tally4Double_ptr dxnorm,
    magma_tally4Double_ptr dxnorm2,
    magma_tally4DoubleComplex_ptr dC, magma_tally4_int_t lddc,
    magma_tally4Double_ptr dlsticc);

void
magma_tally4_zlarfbx_gpu(
    magma_tally4_int_t m, magma_tally4_int_t k,
    magma_tally4DoubleComplex_ptr V,  magma_tally4_int_t ldv,
    magma_tally4DoubleComplex_ptr dT, magma_tally4_int_t ldt,
    magma_tally4DoubleComplex_ptr c,
    magma_tally4DoubleComplex_ptr dwork);

void
magma_tally4_zlarfg_gpu(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dx0,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex_ptr dtau,
    magma_tally4Double_ptr        dxnorm,
    magma_tally4DoubleComplex_ptr dAkk );

void
magma_tally4_zlarfgtx_gpu(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dx0,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex_ptr dtau,
    magma_tally4Double_ptr        dxnorm,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t iter,
    magma_tally4DoubleComplex_ptr V,  magma_tally4_int_t ldv,
    magma_tally4DoubleComplex_ptr T,  magma_tally4_int_t ldt,
    magma_tally4DoubleComplex_ptr dwork);

void
magma_tally4_zlarfgx_gpu(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dx0,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex_ptr dtau,
    magma_tally4Double_ptr        dxnorm,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t iter);

void
magma_tally4_zlarfx_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr v,
    magma_tally4DoubleComplex_ptr tau,
    magma_tally4DoubleComplex_ptr C,  magma_tally4_int_t ldc,
    magma_tally4Double_ptr        xnorm,
    magma_tally4DoubleComplex_ptr dT, magma_tally4_int_t iter,
    magma_tally4DoubleComplex_ptr work);


  /*
   * Level 1 BLAS (alphabetical order)
   */
void
magma_tally4blas_zswap(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_ptr dy, magma_tally4_int_t incy );

void
magma_tally4blas_zswapblk(
    magma_tally4_order_t order,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t i1, magma_tally4_int_t i2,
    const magma_tally4_int_t *ipiv, magma_tally4_int_t inci,
    magma_tally4_int_t offset );

void
magma_tally4blas_zswapdblk(
    magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t inca,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb, magma_tally4_int_t incb );

  /*
   * Level 2 BLAS (alphabetical order)
   */
void
magma_tally4blas_zgemv(
    magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dy, magma_tally4_int_t incy );


void
magma_tally4blas_zgemv_conjv(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dy, magma_tally4_int_t incy);

magma_tally4_int_t
magma_tally4blas_zhemv(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dy, magma_tally4_int_t incy );

magma_tally4_int_t
magma_tally4blas_zsymv(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dy, magma_tally4_int_t incy );

magma_tally4_int_t
magma_tally4blas_zhemv_work(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dy, magma_tally4_int_t incy,
    magma_tally4DoubleComplex_ptr       dwork, magma_tally4_int_t lwork,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4blas_zsymv_work(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dy, magma_tally4_int_t incy,
    magma_tally4DoubleComplex_ptr       dwork, magma_tally4_int_t lwork,
    magma_tally4_queue_t queue );

  /*
   * Level 3 BLAS (alphabetical order)
   */
void
magma_tally4blas_zgemm(
    magma_tally4_trans_t transA, magma_tally4_trans_t transB,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4blas_zgemm_reduce(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4blas_zhemm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4blas_zsymm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4blas_zsyr2k(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4blas_zher2k(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dB, magma_tally4_int_t lddb,
    double  beta,
    magma_tally4DoubleComplex_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4blas_zsyrk(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4blas_zherk(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    double  alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    double  beta,
    magma_tally4DoubleComplex_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4blas_ztrsm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t transA, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr       dB, magma_tally4_int_t lddb );

void
magma_tally4blas_ztrsm_outofplace(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t transA, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr       dB, magma_tally4_int_t lddb,
    magma_tally4_int_t flag, magma_tally4DoubleComplex_ptr d_dinvA, magma_tally4DoubleComplex_ptr dX );

void
magma_tally4blas_ztrsm_work(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t transA, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr       dB, magma_tally4_int_t lddb,
    magma_tally4_int_t flag, magma_tally4DoubleComplex_ptr d_dinvA, magma_tally4DoubleComplex_ptr dX );


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

#define magma_tally4_zsetvector(           n, hx_src, incx, dy_dst, incy ) \
        magma_tally4_zsetvector_internal(  n, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally4_zgetvector(           n, dx_src, incx, hy_dst, incy ) \
        magma_tally4_zgetvector_internal(  n, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally4_zcopyvector(          n, dx_src, incx, dy_dst, incy ) \
        magma_tally4_zcopyvector_internal( n, dx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally4_zsetvector_async(           n, hx_src, incx, dy_dst, incy, queue ) \
        magma_tally4_zsetvector_async_internal(  n, hx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_tally4_zgetvector_async(           n, dx_src, incx, hy_dst, incy, queue ) \
        magma_tally4_zgetvector_async_internal(  n, dx_src, incx, hy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_tally4_zcopyvector_async(          n, dx_src, incx, dy_dst, incy, queue ) \
        magma_tally4_zcopyvector_async_internal( n, dx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

void
magma_tally4_zsetvector_internal(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex const    *hx_src, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_ptr       dy_dst, magma_tally4_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally4_zgetvector_internal(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dx_src, magma_tally4_int_t incx,
    magma_tally4DoubleComplex          *hy_dst, magma_tally4_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally4_zcopyvector_internal(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dx_src, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_ptr       dy_dst, magma_tally4_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally4_zsetvector_async_internal(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex const    *hx_src, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_ptr       dy_dst, magma_tally4_int_t incy,
    magma_tally4_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally4_zgetvector_async_internal(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dx_src, magma_tally4_int_t incx,
    magma_tally4DoubleComplex          *hy_dst, magma_tally4_int_t incy,
    magma_tally4_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally4_zcopyvector_async_internal(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dx_src, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_ptr       dy_dst, magma_tally4_int_t incy,
    magma_tally4_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// copying sub-matrices (contiguous columns)
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

#define magma_tally4_zsetmatrix(           m, n, hA_src, lda, dB_dst, lddb ) \
        magma_tally4_zsetmatrix_internal(  m, n, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_tally4_zgetmatrix(           m, n, dA_src, ldda, hB_dst, ldb ) \
        magma_tally4_zgetmatrix_internal(  m, n, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )

#define magma_tally4_zcopymatrix(          m, n, dA_src, ldda, dB_dst, lddb ) \
        magma_tally4_zcopymatrix_internal( m, n, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_tally4_zsetmatrix_async(           m, n, hA_src, lda, dB_dst, lddb, queue ) \
        magma_tally4_zsetmatrix_async_internal(  m, n, hA_src, lda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

#define magma_tally4_zgetmatrix_async(           m, n, dA_src, ldda, hB_dst, ldb, queue ) \
        magma_tally4_zgetmatrix_async_internal(  m, n, dA_src, ldda, hB_dst, ldb, queue, __func__, __FILE__, __LINE__ )

#define magma_tally4_zcopymatrix_async(          m, n, dA_src, ldda, dB_dst, lddb, queue ) \
        magma_tally4_zcopymatrix_async_internal( m, n, dA_src, ldda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

void
magma_tally4_zsetmatrix_internal(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex const    *hA_src, magma_tally4_int_t ldha,
    magma_tally4DoubleComplex_ptr       dB_dst, magma_tally4_int_t lddb,
    const char* func, const char* file, int line );

void
magma_tally4_zgetmatrix_internal(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dA_src, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex          *hB_dst, magma_tally4_int_t ldhb,
    const char* func, const char* file, int line );

void
magma_tally4_zcopymatrix_internal(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dA_src, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr       dB_dst, magma_tally4_int_t lddb,
    const char* func, const char* file, int line );

void
magma_tally4_zsetmatrix_async_internal(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex const    *hA_src, magma_tally4_int_t ldha,
    magma_tally4DoubleComplex_ptr       dB_dst, magma_tally4_int_t lddb,
    magma_tally4_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally4_zgetmatrix_async_internal(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dA_src, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex          *hB_dst, magma_tally4_int_t ldhb,
    magma_tally4_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally4_zcopymatrix_async_internal(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dA_src, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr       dB_dst, magma_tally4_int_t lddb,
    magma_tally4_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// Level 1 BLAS (alphabetical order)

// in cublas_v2, result returned through output argument
magma_tally4_int_t
magma_tally4_izamax(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx );

// in cublas_v2, result returned through output argument
magma_tally4_int_t
magma_tally4_izamin(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx );

// in cublas_v2, result returned through output argument
double
magma_tally4_dzasum(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx );

void
magma_tally4_zaxpy(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_ptr       dy, magma_tally4_int_t incy );

void
magma_tally4_zcopy(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_ptr       dy, magma_tally4_int_t incy );

// in cublas_v2, result returned through output argument
magma_tally4DoubleComplex
magma_tally4_zdotc(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_const_ptr dy, magma_tally4_int_t incy );

// in cublas_v2, result returned through output argument
magma_tally4DoubleComplex
magma_tally4_zdotu(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_const_ptr dy, magma_tally4_int_t incy );

// in cublas_v2, result returned through output argument
double
magma_tally4_dznrm2(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx );

void
magma_tally4_zrot(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_ptr dy, magma_tally4_int_t incy,
    double dc, magma_tally4DoubleComplex ds );

void
magma_tally4_zdrot(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_ptr dy, magma_tally4_int_t incy,
    double dc, double ds );

#ifdef REAL
void
magma_tally4_zrotm(
    magma_tally4_int_t n,
    magma_tally4Double_ptr dx, magma_tally4_int_t incx,
    magma_tally4Double_ptr dy, magma_tally4_int_t incy,
    magma_tally4Double_const_ptr param );

void
magma_tally4_zrotmg(
    magma_tally4Double_ptr d1, magma_tally4Double_ptr       d2,
    magma_tally4Double_ptr x1, magma_tally4Double_const_ptr y1,
    magma_tally4Double_ptr param );
#endif

void
magma_tally4_zscal(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dx, magma_tally4_int_t incx );

void
magma_tally4_zdscal(
    magma_tally4_int_t n,
    double alpha,
    magma_tally4DoubleComplex_ptr dx, magma_tally4_int_t incx );

void
magma_tally4_zswap(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_ptr dy, magma_tally4_int_t incy );

// ========================================
// Level 2 BLAS (alphabetical order)

void
magma_tally4_zgemv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dy, magma_tally4_int_t incy );

void
magma_tally4_zgerc(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_const_ptr dy, magma_tally4_int_t incy,
    magma_tally4DoubleComplex_ptr       dA, magma_tally4_int_t ldda );

void
magma_tally4_zgeru(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_const_ptr dy, magma_tally4_int_t incy,
    magma_tally4DoubleComplex_ptr       dA, magma_tally4_int_t ldda );

void
magma_tally4_zhemv(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dy, magma_tally4_int_t incy );

void
magma_tally4_zher(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    double alpha,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_ptr       dA, magma_tally4_int_t ldda );

void
magma_tally4_zher2(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4DoubleComplex_const_ptr dy, magma_tally4_int_t incy,
    magma_tally4DoubleComplex_ptr       dA, magma_tally4_int_t ldda );

void
magma_tally4_ztrmv(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_diag_t diag,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr       dx, magma_tally4_int_t incx );

void
magma_tally4_ztrsv(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_diag_t diag,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr       dx, magma_tally4_int_t incx );

// ========================================
// Level 3 BLAS (alphabetical order)

void
magma_tally4_zgemm(
    magma_tally4_trans_t transA, magma_tally4_trans_t transB,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4_zsymm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4_zhemm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4_zsyr2k(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4_zher2k(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_const_ptr dB, magma_tally4_int_t lddb,
    double beta,
    magma_tally4DoubleComplex_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4_zsyrk(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4_zherk(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    double alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    double beta,
    magma_tally4DoubleComplex_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4_ztrmm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr       dB, magma_tally4_int_t lddb );

void
magma_tally4_ztrsm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr       dB, magma_tally4_int_t lddb );


#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif  /* MAGMA_tally4BLAS_Z_H */
