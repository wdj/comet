/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally2blas_z.h normal z -> d, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally2BLAS_D_H
#define MAGMA_tally2BLAS_D_H

#include "magma_tally2_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_tally2blas_dtranspose_inplace(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda );

void
magma_tally2blas_dtranspose(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2Double_ptr       dAT, magma_tally2_int_t lddat );

void
magma_tally2blas_dgetmatrix_transpose(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dAT,   magma_tally2_int_t ldda,
    double          *hA,    magma_tally2_int_t lda,
    magma_tally2Double_ptr       dwork, magma_tally2_int_t lddwork, magma_tally2_int_t nb );

void
magma_tally2blas_dsetmatrix_transpose(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const double *hA,    magma_tally2_int_t lda,
    magma_tally2Double_ptr    dAT,   magma_tally2_int_t ldda,
    magma_tally2Double_ptr    dwork, magma_tally2_int_t lddwork, magma_tally2_int_t nb );

  /*
   * RBT-related functions
   */
void
magma_tally2blas_dprbt(
    magma_tally2_int_t n, 
    double *dA, magma_tally2_int_t ldda, 
    double *du, double *dv);

void
magma_tally2blas_dprbt_mv(
    magma_tally2_int_t n, 
    double *dv, double *db);

void
magma_tally2blas_dprbt_mtv(
    magma_tally2_int_t n, 
    double *du, double *db);

void
magma_tally2blas_daxpycp2(
    magma_tally2_int_t m, double *r, double *x,
    const double *b);

  /*
   * Multi-GPU copy functions
   */
void
magma_tally2blas_dgetmatrix_transpose_mgpu(
    magma_tally2_int_t ngpu, magma_tally2_queue_t queues[][2],
    magma_tally2Double_const_ptr const dAT[],    magma_tally2_int_t ldda,
    double                *hA,       magma_tally2_int_t lda,
    magma_tally2Double_ptr             dwork[],  magma_tally2_int_t lddw,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb );

void
magma_tally2blas_dsetmatrix_transpose_mgpu(
    magma_tally2_int_t ngpu, magma_tally2_queue_t queues[][2],
    const double *hA,      magma_tally2_int_t lda,
    magma_tally2Double_ptr    dAT[],   magma_tally2_int_t ldda,
    magma_tally2Double_ptr    dwork[], magma_tally2_int_t lddw,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb );

void
magma_tally2_dgetmatrix_1D_col_bcyclic(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr const dA[], magma_tally2_int_t ldda,
    double                *hA,   magma_tally2_int_t lda,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb );

void
magma_tally2_dsetmatrix_1D_col_bcyclic(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const double *hA,   magma_tally2_int_t lda,
    magma_tally2Double_ptr    dA[], magma_tally2_int_t ldda,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb );

void
magma_tally2_dgetmatrix_1D_row_bcyclic(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr const dA[], magma_tally2_int_t ldda,
    double                *hA,   magma_tally2_int_t lda,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb );

void
magma_tally2_dsetmatrix_1D_row_bcyclic(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const double *hA,   magma_tally2_int_t lda,
    magma_tally2Double_ptr    dA[], magma_tally2_int_t ldda,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb );

// in src/dsytrd_mgpu.cpp
// TODO rename dsetmatrix_sy or similar
magma_tally2_int_t
magma_tally2_dhtodhe(
    magma_tally2_int_t ngpu, magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    double     *A,   magma_tally2_int_t lda,
    magma_tally2Double_ptr dA[], magma_tally2_int_t ldda,
    magma_tally2_queue_t queues[][10],
    magma_tally2_int_t *info );

// in src/dpotrf3_mgpu.cpp
// TODO same as magma_tally2_dhtodhe?
magma_tally2_int_t
magma_tally2_dhtodpo(
    magma_tally2_int_t ngpu, magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t off_i, magma_tally2_int_t off_j, magma_tally2_int_t nb,
    double     *A,   magma_tally2_int_t lda,
    magma_tally2Double_ptr dA[], magma_tally2_int_t ldda,
    magma_tally2_queue_t queues[][3],
    magma_tally2_int_t *info );

// in src/dpotrf3_mgpu.cpp
// TODO rename dgetmatrix_sy or similar
magma_tally2_int_t
magma_tally2_ddtohpo(
    magma_tally2_int_t ngpu, magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t off_i, magma_tally2_int_t off_j, magma_tally2_int_t nb, magma_tally2_int_t NB,
    double     *A,   magma_tally2_int_t lda,
    magma_tally2Double_ptr dA[], magma_tally2_int_t ldda,
    magma_tally2_queue_t queues[][3],
    magma_tally2_int_t *info );


  /*
   * Multi-GPU BLAS functions (alphabetical order)
   */
void
magma_tally2blas_dsymm_mgpu(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_ptr dA[],    magma_tally2_int_t ldda,  magma_tally2_int_t offset,
    magma_tally2Double_ptr dB[],    magma_tally2_int_t lddb,
    double beta,
    magma_tally2Double_ptr dC[],    magma_tally2_int_t lddc,
    magma_tally2Double_ptr dwork[], magma_tally2_int_t dworksiz,
    double    *C,       magma_tally2_int_t ldc,
    double    *work[],  magma_tally2_int_t worksiz,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_event_t redevents[][20], magma_tally2_int_t nbevents );

void
magma_tally2blas_dsymm_mgpu_com(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_ptr dA[],    magma_tally2_int_t ldda,  magma_tally2_int_t offset,
    magma_tally2Double_ptr dB[],    magma_tally2_int_t lddb,
    double beta,
    magma_tally2Double_ptr dC[],    magma_tally2_int_t lddc,
    magma_tally2Double_ptr dwork[], magma_tally2_int_t dworksiz,
    double    *C,       magma_tally2_int_t ldc,
    double    *work[],  magma_tally2_int_t worksiz,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_event_t redevents[][Magma_tally2MaxGPUs*Magma_tally2MaxGPUs+10], magma_tally2_int_t nbevents,
    magma_tally2_int_t gnode[Magma_tally2MaxGPUs][Magma_tally2MaxGPUs+2], magma_tally2_int_t nbcmplx );

void
magma_tally2blas_dsymm_mgpu_spec(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_ptr dA[],    magma_tally2_int_t ldda,  magma_tally2_int_t offset,
    magma_tally2Double_ptr dB[],    magma_tally2_int_t lddb,
    double beta,
    magma_tally2Double_ptr dC[],    magma_tally2_int_t lddc,
    magma_tally2Double_ptr dwork[], magma_tally2_int_t dworksiz,
    double    *C,       magma_tally2_int_t ldc,
    double    *work[],  magma_tally2_int_t worksiz,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_event_t redevents[][Magma_tally2MaxGPUs*Magma_tally2MaxGPUs+10], magma_tally2_int_t nbevents,
    magma_tally2_int_t gnode[Magma_tally2MaxGPUs][Magma_tally2MaxGPUs+2], magma_tally2_int_t nbcmplx );

void
magma_tally2blas_dsymm_mgpu_spec33(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_ptr dA[],    magma_tally2_int_t ldda,  magma_tally2_int_t offset,
    magma_tally2Double_ptr dB[],    magma_tally2_int_t lddb,
    double beta,
    magma_tally2Double_ptr dC[],    magma_tally2_int_t lddc,
    magma_tally2Double_ptr dVIN[],  magma_tally2_int_t lddv, magma_tally2_int_t v_offset,
    magma_tally2Double_ptr dwork[], magma_tally2_int_t dworksiz,
    double *C,       magma_tally2_int_t ldc,
    double *work[],  magma_tally2_int_t worksiz,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_event_t redevents[][Magma_tally2MaxGPUs*Magma_tally2MaxGPUs+10], magma_tally2_int_t nbevents,
    magma_tally2_int_t gnode[Magma_tally2MaxGPUs][Magma_tally2MaxGPUs+2], magma_tally2_int_t nbcmplx );

magma_tally2_int_t
magma_tally2blas_dsymv_mgpu(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr const d_lA[], magma_tally2_int_t ldda, magma_tally2_int_t offset,
    double const *x,         magma_tally2_int_t incx,
    double beta,             
    double       *y,         magma_tally2_int_t incy,
    double       *hwork,     magma_tally2_int_t lhwork,
    magma_tally2Double_ptr    dwork[],   magma_tally2_int_t ldwork,
    magma_tally2_int_t ngpu,
    magma_tally2_int_t nb,
    magma_tally2_queue_t queues[] );

magma_tally2_int_t
magma_tally2blas_dsymv_mgpu_sync(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr const d_lA[], magma_tally2_int_t ldda, magma_tally2_int_t offset,
    double const *x,         magma_tally2_int_t incx,
    double beta,             
    double       *y,         magma_tally2_int_t incy,
    double       *hwork,     magma_tally2_int_t lhwork,
    magma_tally2Double_ptr    dwork[],   magma_tally2_int_t ldwork,
    magma_tally2_int_t ngpu,
    magma_tally2_int_t nb,
    magma_tally2_queue_t queues[] );

// Ichi's version, in src/dsytrd_mgpu.cpp
void
magma_tally2_dsyr2k_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t nb, magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2Double_ptr dB[], magma_tally2_int_t lddb, magma_tally2_int_t b_offset,
    double beta,
    magma_tally2Double_ptr dC[], magma_tally2_int_t lddc, magma_tally2_int_t c_offset,
    magma_tally2_int_t nqueue, magma_tally2_queue_t queues[][10] );

void
magma_tally2blas_dsyr2k_mgpu_spec(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2Double_ptr dA[], magma_tally2_int_t ldda, magma_tally2_int_t a_offset,
    magma_tally2Double_ptr dB[], magma_tally2_int_t lddb, magma_tally2_int_t b_offset,
    double beta,
    magma_tally2Double_ptr dC[], magma_tally2_int_t lddc, magma_tally2_int_t c_offset,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb, magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue );

void
magma_tally2blas_dsyr2k_mgpu_spec324(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_ptr dVIN[], magma_tally2_int_t lddv, magma_tally2_int_t v_offset,
    magma_tally2Double_ptr dWIN[], magma_tally2_int_t lddw, magma_tally2_int_t w_offset,
    double beta,
    magma_tally2Double_ptr dC[],   magma_tally2_int_t lddc, magma_tally2_int_t c_offset,
    magma_tally2Double_ptr dwork[], magma_tally2_int_t lndwork,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_event_t redevents[][Magma_tally2MaxGPUs*Magma_tally2MaxGPUs+10], magma_tally2_int_t nbevents );

void
magma_tally2blas_dsyr2k_mgpu_spec325(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_ptr dVIN[], magma_tally2_int_t lddv, magma_tally2_int_t v_offset,
    magma_tally2Double_ptr dWIN[], magma_tally2_int_t lddw, magma_tally2_int_t w_offset,
    double beta,
    magma_tally2Double_ptr dC[],   magma_tally2_int_t lddc, magma_tally2_int_t c_offset,
    magma_tally2Double_ptr dwork[], magma_tally2_int_t lndwork,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb,
    double *harray[],
    magma_tally2Double_ptr *darray[],
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_event_t redevents[][Magma_tally2MaxGPUs*Magma_tally2MaxGPUs+10], magma_tally2_int_t nbevents );

void
magma_tally2blas_dsyr2k_mgpu2(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2Double_ptr dA[], magma_tally2_int_t ldda, magma_tally2_int_t a_offset,
    magma_tally2Double_ptr dB[], magma_tally2_int_t lddb, magma_tally2_int_t b_offset,
    double beta,
    magma_tally2Double_ptr dC[], magma_tally2_int_t lddc, magma_tally2_int_t c_offset,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue );

// in src/dpotrf_mgpu_right.cpp
void
magma_tally2_dsyrk_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t nb, magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2Double_ptr dB[], magma_tally2_int_t lddb, magma_tally2_int_t b_offset,
    double beta,
    magma_tally2Double_ptr dC[], magma_tally2_int_t lddc, magma_tally2_int_t c_offset,
    magma_tally2_int_t nqueue, magma_tally2_queue_t queues[][10]);

// in src/dpotrf_mgpu_right.cpp
void
magma_tally2_dsyrk_mgpu2(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t nb, magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2Double_ptr dB[], magma_tally2_int_t lddb, magma_tally2_int_t b_offset,
    double beta,
    magma_tally2Double_ptr dC[], magma_tally2_int_t lddc, magma_tally2_int_t c_offset,
    magma_tally2_int_t nqueue, magma_tally2_queue_t queues[][10]);


  /*
   * LAPACK auxiliary functions (alphabetical order)
   */
void
magma_tally2blas_dgeadd(
    magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr       dB, magma_tally2_int_t lddb );

void
magma_tally2blas_dlacpy(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr       dB, magma_tally2_int_t lddb );

void
magma_tally2blas_dlacpy_cnjg(
    magma_tally2_int_t n, double *dA1, magma_tally2_int_t lda1,
    double *dA2, magma_tally2_int_t lda2);

double
magma_tally2blas_dlange(
    magma_tally2_norm_t norm,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dwork );

double
magma_tally2blas_dlansy(
    magma_tally2_norm_t norm, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dwork );

double
magma_tally2blas_dlansy(
    magma_tally2_norm_t norm, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dwork );

void
magma_tally2blas_dlarfg(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dalpha, magma_tally2Double_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_ptr dtau );

void
magma_tally2blas_dlarfg_work(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dalpha, magma_tally2Double_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_ptr dtau, magma_tally2Double_ptr dwork );

void
magma_tally2blas_dlascl(
    magma_tally2_type_t type, magma_tally2_int_t kl, magma_tally2_int_t ku,
    double cfrom, double cto,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info );

void
magma_tally2blas_dlascl_2x2(
    magma_tally2_type_t type, magma_tally2_int_t m,
    double *dW, magma_tally2_int_t lddw,
    double *dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info );

void
magma_tally2blas_dlascl2(
    magma_tally2_type_t type,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dD,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info );

void
magma_tally2blas_dlascl_diag(
    magma_tally2_type_t type, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dD, magma_tally2_int_t lddd,
          magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info );

void
magma_tally2blas_dlaset(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    double offdiag, double diag,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda );

void
magma_tally2blas_dlaset_band(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    double offdiag, double diag,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda);

void
magma_tally2blas_dlaswp(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dAT, magma_tally2_int_t ldda,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci );

void
magma_tally2blas_dlaswp2(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dAT, magma_tally2_int_t ldda,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    magma_tally2Int_const_ptr d_ipiv, magma_tally2_int_t inci );

void
magma_tally2blas_dlaswpx(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldx, magma_tally2_int_t ldy,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci );

void
magma_tally2blas_dsymmetrize(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda );

void
magma_tally2blas_dsymmetrize_tiles(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t ntile, magma_tally2_int_t mstride, magma_tally2_int_t nstride );

void
magma_tally2blas_dtrtri_diag(
    magma_tally2_uplo_t uplo, magma_tally2_diag_t diag, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr d_dinvA);

  /*
   * to cleanup (alphabetical order)
   */
void
magma_tally2blas_dnrm2_adjust(
    magma_tally2_int_t k,
    magma_tally2Double_ptr dxnorm,
    magma_tally2Double_ptr dc);

void
magma_tally2blas_dnrm2_check(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dxnorm,
    magma_tally2Double_ptr dlsticc);

void
magma_tally2blas_dnrm2_cols(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dxnorm);

void
magma_tally2blas_dnrm2_row_check_adjust(
    magma_tally2_int_t k, double tol,
    magma_tally2Double_ptr dxnorm,
    magma_tally2Double_ptr dxnorm2,
    magma_tally2Double_ptr dC, magma_tally2_int_t lddc,
    magma_tally2Double_ptr dlsticc);

void
magma_tally2_dlarfbx_gpu(
    magma_tally2_int_t m, magma_tally2_int_t k,
    magma_tally2Double_ptr V,  magma_tally2_int_t ldv,
    magma_tally2Double_ptr dT, magma_tally2_int_t ldt,
    magma_tally2Double_ptr c,
    magma_tally2Double_ptr dwork);

void
magma_tally2_dlarfg_gpu(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dx0,
    magma_tally2Double_ptr dx,
    magma_tally2Double_ptr dtau,
    magma_tally2Double_ptr        dxnorm,
    magma_tally2Double_ptr dAkk );

void
magma_tally2_dlarfgtx_gpu(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dx0,
    magma_tally2Double_ptr dx,
    magma_tally2Double_ptr dtau,
    magma_tally2Double_ptr        dxnorm,
    magma_tally2Double_ptr dA, magma_tally2_int_t iter,
    magma_tally2Double_ptr V,  magma_tally2_int_t ldv,
    magma_tally2Double_ptr T,  magma_tally2_int_t ldt,
    magma_tally2Double_ptr dwork);

void
magma_tally2_dlarfgx_gpu(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dx0,
    magma_tally2Double_ptr dx,
    magma_tally2Double_ptr dtau,
    magma_tally2Double_ptr        dxnorm,
    magma_tally2Double_ptr dA, magma_tally2_int_t iter);

void
magma_tally2_dlarfx_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr v,
    magma_tally2Double_ptr tau,
    magma_tally2Double_ptr C,  magma_tally2_int_t ldc,
    magma_tally2Double_ptr        xnorm,
    magma_tally2Double_ptr dT, magma_tally2_int_t iter,
    magma_tally2Double_ptr work);


  /*
   * Level 1 BLAS (alphabetical order)
   */
void
magma_tally2blas_dswap(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_ptr dy, magma_tally2_int_t incy );

void
magma_tally2blas_dswapblk(
    magma_tally2_order_t order,
    magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t i1, magma_tally2_int_t i2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci,
    magma_tally2_int_t offset );

void
magma_tally2blas_dswapdblk(
    magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t inca,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb, magma_tally2_int_t incb );

  /*
   * Level 2 BLAS (alphabetical order)
   */
void
magma_tally2blas_dgemv(
    magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx,
    double beta,
    magma_tally2Double_ptr       dy, magma_tally2_int_t incy );


void
magma_tally2blas_dgemv_conjv(
    magma_tally2_int_t m, magma_tally2_int_t n, double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx,
    double beta,
    magma_tally2Double_ptr dy, magma_tally2_int_t incy);

magma_tally2_int_t
magma_tally2blas_dsymv(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx,
    double beta,
    magma_tally2Double_ptr       dy, magma_tally2_int_t incy );

magma_tally2_int_t
magma_tally2blas_dsymv(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx,
    double beta,
    magma_tally2Double_ptr       dy, magma_tally2_int_t incy );

magma_tally2_int_t
magma_tally2blas_dsymv_work(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx,
    double beta,
    magma_tally2Double_ptr       dy, magma_tally2_int_t incy,
    magma_tally2Double_ptr       dwork, magma_tally2_int_t lwork,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2blas_dsymv_work(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx,
    double beta,
    magma_tally2Double_ptr       dy, magma_tally2_int_t incy,
    magma_tally2Double_ptr       dwork, magma_tally2_int_t lwork,
    magma_tally2_queue_t queue );

  /*
   * Level 3 BLAS (alphabetical order)
   */
void
magma_tally2blas_dgemm(
    magma_tally2_trans_t transA, magma_tally2_trans_t transB,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dB, magma_tally2_int_t lddb,
    double beta,
    magma_tally2Double_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2blas_dgemm_reduce(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dB, magma_tally2_int_t lddb,
    double beta,
    magma_tally2Double_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2blas_dsymm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dB, magma_tally2_int_t lddb,
    double beta,
    magma_tally2Double_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2blas_dsymm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dB, magma_tally2_int_t lddb,
    double beta,
    magma_tally2Double_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2blas_dsyr2k(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dB, magma_tally2_int_t lddb,
    double beta,
    magma_tally2Double_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2blas_dsyr2k(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dB, magma_tally2_int_t lddb,
    double  beta,
    magma_tally2Double_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2blas_dsyrk(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    double beta,
    magma_tally2Double_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2blas_dsyrk(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    double  alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    double  beta,
    magma_tally2Double_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2blas_dtrsm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t transA, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr       dB, magma_tally2_int_t lddb );

void
magma_tally2blas_dtrsm_outofplace(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t transA, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr       dB, magma_tally2_int_t lddb,
    magma_tally2_int_t flag, magma_tally2Double_ptr d_dinvA, magma_tally2Double_ptr dX );

void
magma_tally2blas_dtrsm_work(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t transA, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr       dB, magma_tally2_int_t lddb,
    magma_tally2_int_t flag, magma_tally2Double_ptr d_dinvA, magma_tally2Double_ptr dX );


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

#define magma_tally2_dsetvector(           n, hx_src, incx, dy_dst, incy ) \
        magma_tally2_dsetvector_internal(  n, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally2_dgetvector(           n, dx_src, incx, hy_dst, incy ) \
        magma_tally2_dgetvector_internal(  n, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally2_dcopyvector(          n, dx_src, incx, dy_dst, incy ) \
        magma_tally2_dcopyvector_internal( n, dx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally2_dsetvector_async(           n, hx_src, incx, dy_dst, incy, queue ) \
        magma_tally2_dsetvector_async_internal(  n, hx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_tally2_dgetvector_async(           n, dx_src, incx, hy_dst, incy, queue ) \
        magma_tally2_dgetvector_async_internal(  n, dx_src, incx, hy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_tally2_dcopyvector_async(          n, dx_src, incx, dy_dst, incy, queue ) \
        magma_tally2_dcopyvector_async_internal( n, dx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

void
magma_tally2_dsetvector_internal(
    magma_tally2_int_t n,
    double const    *hx_src, magma_tally2_int_t incx,
    magma_tally2Double_ptr       dy_dst, magma_tally2_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally2_dgetvector_internal(
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dx_src, magma_tally2_int_t incx,
    double          *hy_dst, magma_tally2_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally2_dcopyvector_internal(
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dx_src, magma_tally2_int_t incx,
    magma_tally2Double_ptr       dy_dst, magma_tally2_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally2_dsetvector_async_internal(
    magma_tally2_int_t n,
    double const    *hx_src, magma_tally2_int_t incx,
    magma_tally2Double_ptr       dy_dst, magma_tally2_int_t incy,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally2_dgetvector_async_internal(
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dx_src, magma_tally2_int_t incx,
    double          *hy_dst, magma_tally2_int_t incy,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally2_dcopyvector_async_internal(
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dx_src, magma_tally2_int_t incx,
    magma_tally2Double_ptr       dy_dst, magma_tally2_int_t incy,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// copying sub-matrices (contiguous columns)
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

#define magma_tally2_dsetmatrix(           m, n, hA_src, lda, dB_dst, lddb ) \
        magma_tally2_dsetmatrix_internal(  m, n, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_tally2_dgetmatrix(           m, n, dA_src, ldda, hB_dst, ldb ) \
        magma_tally2_dgetmatrix_internal(  m, n, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )

#define magma_tally2_dcopymatrix(          m, n, dA_src, ldda, dB_dst, lddb ) \
        magma_tally2_dcopymatrix_internal( m, n, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_tally2_dsetmatrix_async(           m, n, hA_src, lda, dB_dst, lddb, queue ) \
        magma_tally2_dsetmatrix_async_internal(  m, n, hA_src, lda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

#define magma_tally2_dgetmatrix_async(           m, n, dA_src, ldda, hB_dst, ldb, queue ) \
        magma_tally2_dgetmatrix_async_internal(  m, n, dA_src, ldda, hB_dst, ldb, queue, __func__, __FILE__, __LINE__ )

#define magma_tally2_dcopymatrix_async(          m, n, dA_src, ldda, dB_dst, lddb, queue ) \
        magma_tally2_dcopymatrix_async_internal( m, n, dA_src, ldda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

void
magma_tally2_dsetmatrix_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    double const    *hA_src, magma_tally2_int_t ldha,
    magma_tally2Double_ptr       dB_dst, magma_tally2_int_t lddb,
    const char* func, const char* file, int line );

void
magma_tally2_dgetmatrix_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA_src, magma_tally2_int_t ldda,
    double          *hB_dst, magma_tally2_int_t ldhb,
    const char* func, const char* file, int line );

void
magma_tally2_dcopymatrix_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA_src, magma_tally2_int_t ldda,
    magma_tally2Double_ptr       dB_dst, magma_tally2_int_t lddb,
    const char* func, const char* file, int line );

void
magma_tally2_dsetmatrix_async_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    double const    *hA_src, magma_tally2_int_t ldha,
    magma_tally2Double_ptr       dB_dst, magma_tally2_int_t lddb,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally2_dgetmatrix_async_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA_src, magma_tally2_int_t ldda,
    double          *hB_dst, magma_tally2_int_t ldhb,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally2_dcopymatrix_async_internal(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA_src, magma_tally2_int_t ldda,
    magma_tally2Double_ptr       dB_dst, magma_tally2_int_t lddb,
    magma_tally2_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// Level 1 BLAS (alphabetical order)

// in cublas_v2, result returned through output argument
magma_tally2_int_t
magma_tally2_idamax(
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx );

// in cublas_v2, result returned through output argument
magma_tally2_int_t
magma_tally2_idamin(
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx );

// in cublas_v2, result returned through output argument
double
magma_tally2_dasum(
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx );

void
magma_tally2_daxpy(
    magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_ptr       dy, magma_tally2_int_t incy );

void
magma_tally2_dcopy(
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_ptr       dy, magma_tally2_int_t incy );

// in cublas_v2, result returned through output argument
double
magma_tally2_ddot(
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_const_ptr dy, magma_tally2_int_t incy );

// in cublas_v2, result returned through output argument
double
magma_tally2_ddot(
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_const_ptr dy, magma_tally2_int_t incy );

// in cublas_v2, result returned through output argument
double
magma_tally2_dnrm2(
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx );

void
magma_tally2_drot(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_ptr dy, magma_tally2_int_t incy,
    double dc, double ds );

void
magma_tally2_drot(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_ptr dy, magma_tally2_int_t incy,
    double dc, double ds );

#ifdef REAL
void
magma_tally2_drotm(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_ptr dy, magma_tally2_int_t incy,
    magma_tally2Double_const_ptr param );

void
magma_tally2_drotmg(
    magma_tally2Double_ptr d1, magma_tally2Double_ptr       d2,
    magma_tally2Double_ptr x1, magma_tally2Double_const_ptr y1,
    magma_tally2Double_ptr param );
#endif

void
magma_tally2_dscal(
    magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_ptr dx, magma_tally2_int_t incx );

void
magma_tally2_dscal(
    magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_ptr dx, magma_tally2_int_t incx );

void
magma_tally2_dswap(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_ptr dy, magma_tally2_int_t incy );

// ========================================
// Level 2 BLAS (alphabetical order)

void
magma_tally2_dgemv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx,
    double beta,
    magma_tally2Double_ptr       dy, magma_tally2_int_t incy );

void
magma_tally2_dger(
    magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_const_ptr dy, magma_tally2_int_t incy,
    magma_tally2Double_ptr       dA, magma_tally2_int_t ldda );

void
magma_tally2_dger(
    magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_const_ptr dy, magma_tally2_int_t incy,
    magma_tally2Double_ptr       dA, magma_tally2_int_t ldda );

void
magma_tally2_dsymv(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx,
    double beta,
    magma_tally2Double_ptr       dy, magma_tally2_int_t incy );

void
magma_tally2_dsyr(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_ptr       dA, magma_tally2_int_t ldda );

void
magma_tally2_dsyr2(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_const_ptr dy, magma_tally2_int_t incy,
    magma_tally2Double_ptr       dA, magma_tally2_int_t ldda );

void
magma_tally2_dtrmv(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_diag_t diag,
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr       dx, magma_tally2_int_t incx );

void
magma_tally2_dtrsv(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_diag_t diag,
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr       dx, magma_tally2_int_t incx );

// ========================================
// Level 3 BLAS (alphabetical order)

void
magma_tally2_dgemm(
    magma_tally2_trans_t transA, magma_tally2_trans_t transB,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dB, magma_tally2_int_t lddb,
    double beta,
    magma_tally2Double_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2_dsymm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dB, magma_tally2_int_t lddb,
    double beta,
    magma_tally2Double_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2_dsymm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dB, magma_tally2_int_t lddb,
    double beta,
    magma_tally2Double_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2_dsyr2k(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dB, magma_tally2_int_t lddb,
    double beta,
    magma_tally2Double_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2_dsyr2k(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_const_ptr dB, magma_tally2_int_t lddb,
    double beta,
    magma_tally2Double_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2_dsyrk(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    double beta,
    magma_tally2Double_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2_dsyrk(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    double beta,
    magma_tally2Double_ptr       dC, magma_tally2_int_t lddc );

void
magma_tally2_dtrmm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr       dB, magma_tally2_int_t lddb );

void
magma_tally2_dtrsm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr       dB, magma_tally2_int_t lddb );


#ifdef __cplusplus
}
#endif

#undef REAL

#endif  /* MAGMA_tally2BLAS_D_H */
