/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
*/

#ifndef MAGMA_tally2BLAS_Z_H
#define MAGMA_tally2BLAS_Z_H

#include "magma_tally2_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_tally2blas_ztranspose_inplace(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda );

void
magma_tally2blas_ztranspose(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr       dAT, magma_tally2_int_t lddat );

void
magma_tally2blas_zgetmatrix_transpose(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dAT,   magma_tally2_int_t ldda,
    magma_tally2DoubleComplex          *hA,    magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr       dwork, magma_tally2_int_t lddwork, magma_tally2_int_t nb );

void
magma_tally2blas_zsetmatrix_transpose(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const magma_tally2DoubleComplex *hA,    magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr    dAT,   magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr    dwork, magma_tally2_int_t lddwork, magma_tally2_int_t nb );

  /*
   * RBT-related functions
   */
void
magma_tally2blas_zprbt(
    magma_tally2_int_t n, 
    magma_tally2DoubleComplex *dA, magma_tally2_int_t ldda, 
    magma_tally2DoubleComplex *du, magma_tally2DoubleComplex *dv);

void
magma_tally2blas_zprbt_mv(
    magma_tally2_int_t n, 
    magma_tally2DoubleComplex *dv, magma_tally2DoubleComplex *db);

void
magma_tally2blas_zprbt_mtv(
    magma_tally2_int_t n, 
    magma_tally2DoubleComplex *du, magma_tally2DoubleComplex *db);

void
magma_tally2blas_zaxpycp2(
    magma_tally2_int_t m, magma_tally2DoubleComplex *r, magma_tally2DoubleComplex *x,
    const magma_tally2DoubleComplex *b);

  /*
   * Multi-GPU copy functions
   */
void
magma_tally2blas_zgetmatrix_transpose_mgpu(
    magma_tally2_int_t ngpu, magma_tally2_queue_t queues[][2],
    magma_tally2DoubleComplex_const_ptr const dAT[],    magma_tally2_int_t ldda,
    magma_tally2DoubleComplex                *hA,       magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr             dwork[],  magma_tally2_int_t lddw,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb );

void
magma_tally2blas_zsetmatrix_transpose_mgpu(
    magma_tally2_int_t ngpu, magma_tally2_queue_t queues[][2],
    const magma_tally2DoubleComplex *hA,      magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr    dAT[],   magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr    dwork[], magma_tally2_int_t lddw,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb );

void
magma_tally2_zgetmatrix_1D_col_bcyclic(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr const dA[], magma_tally2_int_t ldda,
    magma_tally2DoubleComplex                *hA,   magma_tally2_int_t lda,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb );

void
magma_tally2_zsetmatrix_1D_col_bcyclic(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const magma_tally2DoubleComplex *hA,   magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr    dA[], magma_tally2_int_t ldda,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb );

void
magma_tally2_zgetmatrix_1D_row_bcyclic(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr const dA[], magma_tally2_int_t ldda,
    magma_tally2DoubleComplex                *hA,   magma_tally2_int_t lda,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb );

void
magma_tally2_zsetmatrix_1D_row_bcyclic(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const magma_tally2DoubleComplex *hA,   magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr    dA[], magma_tally2_int_t ldda,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb );

// in src/zhetrd_mgpu.cpp
// TODO rename zsetmatrix_sy or similar
magma_tally2_int_t
magma_tally2_zhtodhe(
    magma_tally2_int_t ngpu, magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2DoubleComplex     *A,   magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr dA[], magma_tally2_int_t ldda,
    magma_tally2_queue_t queues[][10],
    magma_tally2_int_t *info );

// in src/zpotrf3_mgpu.cpp
// TODO same as magma_tally2_zhtodhe?
magma_tally2_int_t
magma_tally2_zhtodpo(
    magma_tally2_int_t ngpu, magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t off_i, magma_tally2_int_t off_j, magma_tally2_int_t nb,
    magma_tally2DoubleComplex     *A,   magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr dA[], magma_tally2_int_t ldda,
    magma_tally2_queue_t queues[][3],
    magma_tally2_int_t *info );

// in src/zpotrf3_mgpu.cpp
// TODO rename zgetmatrix_sy or similar
magma_tally2_int_t
magma_tally2_zdtohpo(
    magma_tally2_int_t ngpu, magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t off_i, magma_tally2_int_t off_j, magma_tally2_int_t nb, magma_tally2_int_t NB,
    magma_tally2DoubleComplex     *A,   magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr dA[], magma_tally2_int_t ldda,
    magma_tally2_queue_t queues[][3],
    magma_tally2_int_t *info );


  /*
   * Multi-GPU BLAS functions (alphabetical order)
   */
void
magma_tally2blas_zhemm_mgpu(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dA[],    magma_tally2_int_t ldda,  magma_tally2_int_t offset,
    magma_tally2DoubleComplex_ptr dB[],    magma_tally2_int_t lddb,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dC[],    magma_tally2_int_t lddc,
    magma_tally2DoubleComplex_ptr dwork[], magma_tally2_int_t dworksiz,
    magma_tally2DoubleComplex    *C,       magma_tally2_int_t ldc,
    magma_tally2DoubleComplex    *work[],  magma_tally2_int_t worksiz,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_event_t redevents[][20], magma_tally2_int_t nbevents );

void
magma_tally2blas_zhemm_mgpu_com(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dA[],    magma_tally2_int_t ldda,  magma_tally2_int_t offset,
    magma_tally2DoubleComplex_ptr dB[],    magma_tally2_int_t lddb,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dC[],    magma_tally2_int_t lddc,
    magma_tally2DoubleComplex_ptr dwork[], magma_tally2_int_t dworksiz,
    magma_tally2DoubleComplex    *C,       magma_tally2_int_t ldc,
    magma_tally2DoubleComplex    *work[],  magma_tally2_int_t worksiz,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_event_t redevents[][Magma_tally2MaxGPUs*Magma_tally2MaxGPUs+10], magma_tally2_int_t nbevents,
    magma_tally2_int_t gnode[Magma_tally2MaxGPUs][Magma_tally2MaxGPUs+2], magma_tally2_int_t nbcmplx );

void
magma_tally2blas_zhemm_mgpu_spec(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dA[],    magma_tally2_int_t ldda,  magma_tally2_int_t offset,
    magma_tally2DoubleComplex_ptr dB[],    magma_tally2_int_t lddb,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dC[],    magma_tally2_int_t lddc,
    magma_tally2DoubleComplex_ptr dwork[], magma_tally2_int_t dworksiz,
    magma_tally2DoubleComplex    *C,       magma_tally2_int_t ldc,
    magma_tally2DoubleComplex    *work[],  magma_tally2_int_t worksiz,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_event_t redevents[][Magma_tally2MaxGPUs*Magma_tally2MaxGPUs+10], magma_tally2_int_t nbevents,
    magma_tally2_int_t gnode[Magma_tally2MaxGPUs][Magma_tally2MaxGPUs+2], magma_tally2_int_t nbcmplx );

void
magma_tally2blas_zhemm_mgpu_spec33(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dA[],    magma_tally2_int_t ldda,  magma_tally2_int_t offset,
    magma_tally2DoubleComplex_ptr dB[],    magma_tally2_int_t lddb,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dC[],    magma_tally2_int_t lddc,
    magma_tally2DoubleComplex_ptr dVIN[],  magma_tally2_int_t lddv, magma_tally2_int_t v_offset,
    magma_tally2DoubleComplex_ptr dwork[], magma_tally2_int_t dworksiz,
    magma_tally2DoubleComplex *C,       magma_tally2_int_t ldc,
    magma_tally2DoubleComplex *work[],  magma_tally2_int_t worksiz,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_event_t redevents[][Magma_tally2MaxGPUs*Magma_tally2MaxGPUs+10], magma_tally2_int_t nbevents,
    magma_tally2_int_t gnode[Magma_tally2MaxGPUs][Magma_tally2MaxGPUs+2], magma_tally2_int_t nbcmplx );

magma_tally2_int_t
magma_tally2blas_zhemv_mgpu(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr const d_lA[], magma_tally2_int_t ldda, magma_tally2_int_t offset,
    magma_tally2DoubleComplex const *x,         magma_tally2_int_t incx,
    magma_tally2DoubleComplex beta,             
    magma_tally2DoubleComplex       *y,         magma_tally2_int_t incy,
    magma_tally2DoubleComplex       *hwork,     magma_tally2_int_t lhwork,
    magma_tally2DoubleComplex_ptr    dwork[],   magma_tally2_int_t ldwork,
    magma_tally2_int_t ngpu,
    magma_tally2_int_t nb,
    magma_tally2_queue_t queues[] );

magma_tally2_int_t
magma_tally2blas_zhemv_mgpu_sync(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr const d_lA[], magma_tally2_int_t ldda, magma_tally2_int_t offset,
    magma_tally2DoubleComplex const *x,         magma_tally2_int_t incx,
    magma_tally2DoubleComplex beta,             
    magma_tally2DoubleComplex       *y,         magma_tally2_int_t incy,
    magma_tally2DoubleComplex       *hwork,     magma_tally2_int_t lhwork,
    magma_tally2DoubleComplex_ptr    dwork[],   magma_tally2_int_t ldwork,
    magma_tally2_int_t ngpu,
    magma_tally2_int_t nb,
    magma_tally2_queue_t queues[] );

// Ichi's version, in src/zhetrd_mgpu.cpp
void
magma_tally2_zher2k_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t nb, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dB[], magma_tally2_int_t lddb, magma_tally2_int_t b_offset,
    double beta,
    magma_tally2DoubleComplex_ptr dC[], magma_tally2_int_t lddc, magma_tally2_int_t c_offset,
    magma_tally2_int_t nqueue, magma_tally2_queue_t queues[][10] );

void
magma_tally2blas_zher2k_mgpu_spec(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dA[], magma_tally2_int_t ldda, magma_tally2_int_t a_offset,
    magma_tally2DoubleComplex_ptr dB[], magma_tally2_int_t lddb, magma_tally2_int_t b_offset,
    double beta,
    magma_tally2DoubleComplex_ptr dC[], magma_tally2_int_t lddc, magma_tally2_int_t c_offset,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb, magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue );

void
magma_tally2blas_zher2k_mgpu_spec324(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dVIN[], magma_tally2_int_t lddv, magma_tally2_int_t v_offset,
    magma_tally2DoubleComplex_ptr dWIN[], magma_tally2_int_t lddw, magma_tally2_int_t w_offset,
    double beta,
    magma_tally2DoubleComplex_ptr dC[],   magma_tally2_int_t lddc, magma_tally2_int_t c_offset,
    magma_tally2DoubleComplex_ptr dwork[], magma_tally2_int_t lndwork,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_event_t redevents[][Magma_tally2MaxGPUs*Magma_tally2MaxGPUs+10], magma_tally2_int_t nbevents );

void
magma_tally2blas_zher2k_mgpu_spec325(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dVIN[], magma_tally2_int_t lddv, magma_tally2_int_t v_offset,
    magma_tally2DoubleComplex_ptr dWIN[], magma_tally2_int_t lddw, magma_tally2_int_t w_offset,
    double beta,
    magma_tally2DoubleComplex_ptr dC[],   magma_tally2_int_t lddc, magma_tally2_int_t c_offset,
    magma_tally2DoubleComplex_ptr dwork[], magma_tally2_int_t lndwork,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb,
    magma_tally2DoubleComplex *harray[],
    magma_tally2DoubleComplex_ptr *darray[],
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_event_t redevents[][Magma_tally2MaxGPUs*Magma_tally2MaxGPUs+10], magma_tally2_int_t nbevents );

void
magma_tally2blas_zher2k_mgpu2(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dA[], magma_tally2_int_t ldda, magma_tally2_int_t a_offset,
    magma_tally2DoubleComplex_ptr dB[], magma_tally2_int_t lddb, magma_tally2_int_t b_offset,
    double beta,
    magma_tally2DoubleComplex_ptr dC[], magma_tally2_int_t lddc, magma_tally2_int_t c_offset,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue );

// in src/zpotrf_mgpu_right.cpp
void
magma_tally2_zherk_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t nb, magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2DoubleComplex_ptr dB[], magma_tally2_int_t lddb, magma_tally2_int_t b_offset,
    double beta,
    magma_tally2DoubleComplex_ptr dC[], magma_tally2_int_t lddc, magma_tally2_int_t c_offset,
    magma_tally2_int_t nqueue, magma_tally2_queue_t queues[][10]);

// in src/zpotrf_mgpu_right.cpp
void
magma_tally2_zherk_mgpu2(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t nb, magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2DoubleComplex_ptr dB[], magma_tally2_int_t lddb, magma_tally2_int_t b_offset,
    double beta,
    magma_tally2DoubleComplex_ptr dC[], magma_tally2_int_t lddc, magma_tally2_int_t c_offset,
    magma_tally2_int_t nqueue, magma_tally2_queue_t queues[][10]);


  /*
   * LAPACK auxiliary functions (alphabetical order)
   */
void
magma_tally2blas_zgeadd(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr       dB, magma_tally2_int_t lddb );

void
magma_tally2blas_zlacpy(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr       dB, magma_tally2_int_t lddb );

void
magma_tally2blas_zlacpy_cnjg(
    magma_tally2_int_t n, magma_tally2DoubleComplex *dA1, magma_tally2_int_t lda1,
    magma_tally2DoubleComplex *dA2, magma_tally2_int_t lda2);

double
magma_tally2blas_zlange(
    magma_tally2_norm_t norm,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dwork );

double
magma_tally2blas_zlanhe(
    magma_tally2_norm_t norm, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dwork );

double
magma_tally2blas_zlansy(
    magma_tally2_norm_t norm, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dwork );

void
magma_tally2blas_zlarfg(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dalpha, magma_tally2DoubleComplex_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr dtau );

void
magma_tally2blas_zlarfg_work(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dalpha, magma_tally2DoubleComplex_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr dtau, magma_tally2DoubleComplex_ptr dwork );

void
magma_tally2blas_zlascl(
    magma_tally2_type_t type, magma_tally2_int_t kl, magma_tally2_int_t ku,
    double cfrom, double cto,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info );

void
magma_tally2blas_zlascl_2x2(
    magma_tally2_type_t type, magma_tally2_int_t m,
    magma_tally2DoubleComplex *dW, magma_tally2_int_t lddw,
    magma_tally2DoubleComplex *dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info );

void
magma_tally2blas_zlascl2(
    magma_tally2_type_t type,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dD,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info );

void
magma_tally2blas_zlascl_diag(
    magma_tally2_type_t type, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dD, magma_tally2_int_t lddd,
          magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info );

void
magma_tally2blas_zlaset(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex offdiag, magma_tally2DoubleComplex diag,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda );

void
magma_tally2blas_zlaset_band(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex offdiag, magma_tally2DoubleComplex diag,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda);

void
magma_tally2blas_zlaswp(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dAT, magma_tally2_int_t ldda,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci );

void
magma_tally2blas_zlaswp2(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dAT, magma_tally2_int_t ldda,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    magma_tally2Int_const_ptr d_ipiv, magma_tally2_int_t inci );

void
magma_tally2blas_zlaswpx(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldx, magma_tally2_int_t ldy,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci );

void
magma_tally2blas_zsymmetrize(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda );

void
magma_tally2blas_zsymmetrize_tiles(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t ntile, magma_tally2_int_t mstride, magma_tally2_int_t nstride );

void
magma_tally2blas_ztrtri_diag(
    magma_tally2_uplo_t uplo, magma_tally2_diag_t diag, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr d_dinvA);

  /*
   * to cleanup (alphabetical order)
   */
void
magma_tally2blas_dznrm2_adjust(
    magma_tally2_int_t k,
    magma_tally2Double_ptr dxnorm,
    magma_tally2DoubleComplex_ptr dc);

void
magma_tally2blas_dznrm2_check(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dxnorm,
    magma_tally2Double_ptr dlsticc);

void
magma_tally2blas_dznrm2_cols(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dxnorm);

void
magma_tally2blas_dznrm2_row_check_adjust(
    magma_tally2_int_t k, double tol,
    magma_tally2Double_ptr dxnorm,
    magma_tally2Double_ptr dxnorm2,
    magma_tally2DoubleComplex_ptr dC, magma_tally2_int_t lddc,
    magma_tally2Double_ptr dlsticc);

void
magma_tally2_zlarfbx_gpu(
    magma_tally2_int_t m, magma_tally2_int_t k,
    magma_tally2DoubleComplex_ptr V,  magma_tally2_int_t ldv,
    magma_tally2DoubleComplex_ptr dT, magma_tally2_int_t ldt,
    magma_tally2DoubleComplex_ptr c,
    magma_tally2DoubleComplex_ptr dwork);

void
magma_tally2_zlarfg_gpu(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dx0,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex_ptr dtau,
    magma_tally2Double_ptr        dxnorm,
    magma_tally2DoubleComplex_ptr dAkk );

void
magma_tally2_zlarfgtx_gpu(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dx0,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex_ptr dtau,
    magma_tally2Double_ptr        dxnorm,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t iter,
    magma_tally2DoubleComplex_ptr V,  magma_tally2_int_t ldv,
    magma_tally2DoubleComplex_ptr T,  magma_tally2_int_t ldt,
    magma_tally2DoubleComplex_ptr dwork);

void
magma_tally2_zlarfgx_gpu(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dx0,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex_ptr dtau,
    magma_tally2Double_ptr        dxnorm,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t iter);

void
magma_tally2_zlarfx_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr v,
    magma_tally2DoubleComplex_ptr tau,
    magma_tally2DoubleComplex_ptr C,  magma_tally2_int_t ldc,
    magma_tally2Double_ptr        xnorm,
    magma_tally2DoubleComplex_ptr dT, magma_tally2_int_t iter,
    magma_tally2DoubleComplex_ptr work);


  /*
   * Level 1 BLAS (alphabetical order)
   */
void
magma_tally2blas_zswap(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr dy, magma_tally2_int_t incy );

void
magma_tally2blas_zswapblk(
    magma_tally2_order_t order,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t i1, magma_tally2_int_t i2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci,
    magma_tally2_int_t offset );

void
magma_tally2blas_zswapdblk(
    magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t inca,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb, magma_tally2_int_t incb );

  /*
   * Level 2 BLAS (alphabetical order)
   */
void
magma_tally2blas_zgemv(
    magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dy, magma_tally2_int_t incy );


void
magma_tally2blas_zgemv_conjv(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dy, magma_tally2_int_t incy);

magma_tally2_int_t
magma_tally2blas_zhemv(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dy, magma_tally2_int_t incy );

magma_tally2_int_t
magma_tally2blas_zsymv(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dy, magma_tally2_int_t incy );

magma_tally2_int_t
magma_tally2blas_zhemv_work(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dy, magma_tally2_int_t incy,
    magma_tally2DoubleComplex_ptr       dwork, magma_tally2_int_t lwork,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2blas_zsymv_work(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dy, magma_tally2_int_t incy,
    magma_tally2DoubleComplex_ptr       dwork, magma_tally2_int_t lwork,
    magma_tally2_queue_t queue );

  /*
   * Level 3 BLAS (alphabetical order)
   */
void
magma_tally2blas_zgemm(
    magma_tally2_trans_t transA, magma_tally2_trans_t transB,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2blas_zgemm_reduce(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2blas_zhemm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2blas_zsymm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2blas_zsyr2k(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2blas_zher2k(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dB, magma_tally2_int_t lddb,
    double  beta,
    magma_tally2DoubleComplex_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2blas_zsyrk(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2blas_zherk(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    double  alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    double  beta,
    magma_tally2DoubleComplex_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2blas_ztrsm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t transA, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr       dB, magma_tally2_int_t lddb );

void
magma_tally2blas_ztrsm_outofplace(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t transA, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr       dB, magma_tally2_int_t lddb,
    magma_tally2_int_t flag, magma_tally2DoubleComplex_ptr d_dinvA, magma_tally2DoubleComplex_ptr dX );

void
magma_tally2blas_ztrsm_work(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t transA, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr       dB, magma_tally2_int_t lddb,
    magma_tally2_int_t flag, magma_tally2DoubleComplex_ptr d_dinvA, magma_tally2DoubleComplex_ptr dX );


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

#define magma_tally2_zsetvector(           n, hx_src, incx, dy_dst, incy ) \
        magma_tally2_zsetvector_internal(  n, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally2_zgetvector(           n, dx_src, incx, hy_dst, incy ) \
        magma_tally2_zgetvector_internal(  n, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally2_zcopyvector(          n, dx_src, incx, dy_dst, incy ) \
        magma_tally2_zcopyvector_internal( n, dx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally2_zsetvector_async(           n, hx_src, incx, dy_dst, incy, queue ) \
        magma_tally2_zsetvector_async_internal(  n, hx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_tally2_zgetvector_async(           n, dx_src, incx, hy_dst, incy, queue ) \
        magma_tally2_zgetvector_async_internal(  n, dx_src, incx, hy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_tally2_zcopyvector_async(          n, dx_src, incx, dy_dst, incy, queue ) \
        magma_tally2_zcopyvector_async_internal( n, dx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

void
magma_tally2_zsetvector_internal(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex const    *hx_src, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr       dy_dst, magma_tally2_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally2_zgetvector_internal(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dx_src, magma_tally2_int_t incx,
    magma_tally2DoubleComplex          *hy_dst, magma_tally2_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally2_zcopyvector_internal(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dx_src, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr       dy_dst, magma_tally2_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally2_zsetvector_async_internal(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex const    *hx_src, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr       dy_dst, magma_tally2_int_t incy,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally2_zgetvector_async_internal(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dx_src, magma_tally2_int_t incx,
    magma_tally2DoubleComplex          *hy_dst, magma_tally2_int_t incy,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally2_zcopyvector_async_internal(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dx_src, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr       dy_dst, magma_tally2_int_t incy,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// copying sub-matrices (contiguous columns)
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

#define magma_tally2_zsetmatrix(           m, n, hA_src, lda, dB_dst, lddb ) \
        magma_tally2_zsetmatrix_internal(  m, n, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_tally2_zgetmatrix(           m, n, dA_src, ldda, hB_dst, ldb ) \
        magma_tally2_zgetmatrix_internal(  m, n, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )

#define magma_tally2_zcopymatrix(          m, n, dA_src, ldda, dB_dst, lddb ) \
        magma_tally2_zcopymatrix_internal( m, n, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_tally2_zsetmatrix_async(           m, n, hA_src, lda, dB_dst, lddb, queue ) \
        magma_tally2_zsetmatrix_async_internal(  m, n, hA_src, lda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

#define magma_tally2_zgetmatrix_async(           m, n, dA_src, ldda, hB_dst, ldb, queue ) \
        magma_tally2_zgetmatrix_async_internal(  m, n, dA_src, ldda, hB_dst, ldb, queue, __func__, __FILE__, __LINE__ )

#define magma_tally2_zcopymatrix_async(          m, n, dA_src, ldda, dB_dst, lddb, queue ) \
        magma_tally2_zcopymatrix_async_internal( m, n, dA_src, ldda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

void
magma_tally2_zsetmatrix_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex const    *hA_src, magma_tally2_int_t ldha,
    magma_tally2DoubleComplex_ptr       dB_dst, magma_tally2_int_t lddb,
    const char* func, const char* file, int line );

void
magma_tally2_zgetmatrix_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA_src, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex          *hB_dst, magma_tally2_int_t ldhb,
    const char* func, const char* file, int line );

void
magma_tally2_zcopymatrix_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA_src, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr       dB_dst, magma_tally2_int_t lddb,
    const char* func, const char* file, int line );

void
magma_tally2_zsetmatrix_async_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex const    *hA_src, magma_tally2_int_t ldha,
    magma_tally2DoubleComplex_ptr       dB_dst, magma_tally2_int_t lddb,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally2_zgetmatrix_async_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA_src, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex          *hB_dst, magma_tally2_int_t ldhb,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally2_zcopymatrix_async_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA_src, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr       dB_dst, magma_tally2_int_t lddb,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// Level 1 BLAS (alphabetical order)

// in cublas_v2, result returned through output argument
magma_tally2_int_t
magma_tally2_izamax(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx );

// in cublas_v2, result returned through output argument
magma_tally2_int_t
magma_tally2_izamin(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx );

// in cublas_v2, result returned through output argument
double
magma_tally2_dzasum(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx );

void
magma_tally2_zaxpy(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr       dy, magma_tally2_int_t incy );

void
magma_tally2_zcopy(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr       dy, magma_tally2_int_t incy );

// in cublas_v2, result returned through output argument
magma_tally2DoubleComplex
magma_tally2_zdotc(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_const_ptr dy, magma_tally2_int_t incy );

// in cublas_v2, result returned through output argument
magma_tally2DoubleComplex
magma_tally2_zdotu(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_const_ptr dy, magma_tally2_int_t incy );

// in cublas_v2, result returned through output argument
double
magma_tally2_dznrm2(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx );

void
magma_tally2_zrot(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr dy, magma_tally2_int_t incy,
    double dc, magma_tally2DoubleComplex ds );

void
magma_tally2_zdrot(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr dy, magma_tally2_int_t incy,
    double dc, double ds );

#ifdef REAL
void
magma_tally2_zrotm(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_ptr dy, magma_tally2_int_t incy,
    magma_tally2Double_const_ptr param );

void
magma_tally2_zrotmg(
    magma_tally2Double_ptr d1, magma_tally2Double_ptr       d2,
    magma_tally2Double_ptr x1, magma_tally2Double_const_ptr y1,
    magma_tally2Double_ptr param );
#endif

void
magma_tally2_zscal(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dx, magma_tally2_int_t incx );

void
magma_tally2_zdscal(
    magma_tally2_int_t n,
    double alpha,
    magma_tally2DoubleComplex_ptr dx, magma_tally2_int_t incx );

void
magma_tally2_zswap(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr dy, magma_tally2_int_t incy );

// ========================================
// Level 2 BLAS (alphabetical order)

void
magma_tally2_zgemv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dy, magma_tally2_int_t incy );

void
magma_tally2_zgerc(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_const_ptr dy, magma_tally2_int_t incy,
    magma_tally2DoubleComplex_ptr       dA, magma_tally2_int_t ldda );

void
magma_tally2_zgeru(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_const_ptr dy, magma_tally2_int_t incy,
    magma_tally2DoubleComplex_ptr       dA, magma_tally2_int_t ldda );

void
magma_tally2_zhemv(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dy, magma_tally2_int_t incy );

void
magma_tally2_zher(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    double alpha,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr       dA, magma_tally2_int_t ldda );

void
magma_tally2_zher2(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_const_ptr dy, magma_tally2_int_t incy,
    magma_tally2DoubleComplex_ptr       dA, magma_tally2_int_t ldda );

void
magma_tally2_ztrmv(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_diag_t diag,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr       dx, magma_tally2_int_t incx );

void
magma_tally2_ztrsv(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_diag_t diag,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr       dx, magma_tally2_int_t incx );

// ========================================
// Level 3 BLAS (alphabetical order)

void
magma_tally2_zgemm(
    magma_tally2_trans_t transA, magma_tally2_trans_t transB,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2_zsymm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2_zhemm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2_zsyr2k(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2_zher2k(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_const_ptr dB, magma_tally2_int_t lddb,
    double beta,
    magma_tally2DoubleComplex_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2_zsyrk(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2_zherk(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    double beta,
    magma_tally2DoubleComplex_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2_ztrmm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr       dB, magma_tally2_int_t lddb );

void
magma_tally2_ztrsm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr       dB, magma_tally2_int_t lddb );


#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif  /* MAGMA_tally2BLAS_Z_H */
