/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproductblas_z.h normal z -> d, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_minproductBLAS_D_H
#define MAGMA_minproductBLAS_D_H

#include "magma_minproduct_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_minproductblas_dtranspose_inplace(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda );

void
magma_minproductblas_dtranspose(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dAT, magma_minproduct_int_t lddat );

void
magma_minproductblas_dgetmatrix_transpose(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dAT,   magma_minproduct_int_t ldda,
    double          *hA,    magma_minproduct_int_t lda,
    magma_minproductDouble_ptr       dwork, magma_minproduct_int_t lddwork, magma_minproduct_int_t nb );

void
magma_minproductblas_dsetmatrix_transpose(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const double *hA,    magma_minproduct_int_t lda,
    magma_minproductDouble_ptr    dAT,   magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr    dwork, magma_minproduct_int_t lddwork, magma_minproduct_int_t nb );

  /*
   * RBT-related functions
   */
void
magma_minproductblas_dprbt(
    magma_minproduct_int_t n, 
    double *dA, magma_minproduct_int_t ldda, 
    double *du, double *dv);

void
magma_minproductblas_dprbt_mv(
    magma_minproduct_int_t n, 
    double *dv, double *db);

void
magma_minproductblas_dprbt_mtv(
    magma_minproduct_int_t n, 
    double *du, double *db);

void
magma_minproductblas_daxpycp2(
    magma_minproduct_int_t m, double *r, double *x,
    const double *b);

  /*
   * Multi-GPU copy functions
   */
void
magma_minproductblas_dgetmatrix_transpose_mgpu(
    magma_minproduct_int_t ngpu, magma_minproduct_queue_t queues[][2],
    magma_minproductDouble_const_ptr const dAT[],    magma_minproduct_int_t ldda,
    double                *hA,       magma_minproduct_int_t lda,
    magma_minproductDouble_ptr             dwork[],  magma_minproduct_int_t lddw,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb );

void
magma_minproductblas_dsetmatrix_transpose_mgpu(
    magma_minproduct_int_t ngpu, magma_minproduct_queue_t queues[][2],
    const double *hA,      magma_minproduct_int_t lda,
    magma_minproductDouble_ptr    dAT[],   magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr    dwork[], magma_minproduct_int_t lddw,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb );

void
magma_minproduct_dgetmatrix_1D_col_bcyclic(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr const dA[], magma_minproduct_int_t ldda,
    double                *hA,   magma_minproduct_int_t lda,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb );

void
magma_minproduct_dsetmatrix_1D_col_bcyclic(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const double *hA,   magma_minproduct_int_t lda,
    magma_minproductDouble_ptr    dA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb );

void
magma_minproduct_dgetmatrix_1D_row_bcyclic(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr const dA[], magma_minproduct_int_t ldda,
    double                *hA,   magma_minproduct_int_t lda,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb );

void
magma_minproduct_dsetmatrix_1D_row_bcyclic(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const double *hA,   magma_minproduct_int_t lda,
    magma_minproductDouble_ptr    dA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb );

// in src/dsytrd_mgpu.cpp
// TODO rename dsetmatrix_sy or similar
magma_minproduct_int_t
magma_minproduct_dhtodhe(
    magma_minproduct_int_t ngpu, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    double     *A,   magma_minproduct_int_t lda,
    magma_minproductDouble_ptr dA[], magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queues[][10],
    magma_minproduct_int_t *info );

// in src/dpotrf3_mgpu.cpp
// TODO same as magma_minproduct_dhtodhe?
magma_minproduct_int_t
magma_minproduct_dhtodpo(
    magma_minproduct_int_t ngpu, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t off_i, magma_minproduct_int_t off_j, magma_minproduct_int_t nb,
    double     *A,   magma_minproduct_int_t lda,
    magma_minproductDouble_ptr dA[], magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queues[][3],
    magma_minproduct_int_t *info );

// in src/dpotrf3_mgpu.cpp
// TODO rename dgetmatrix_sy or similar
magma_minproduct_int_t
magma_minproduct_ddtohpo(
    magma_minproduct_int_t ngpu, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t off_i, magma_minproduct_int_t off_j, magma_minproduct_int_t nb, magma_minproduct_int_t NB,
    double     *A,   magma_minproduct_int_t lda,
    magma_minproductDouble_ptr dA[], magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queues[][3],
    magma_minproduct_int_t *info );


  /*
   * Multi-GPU BLAS functions (alphabetical order)
   */
void
magma_minproductblas_dsymm_mgpu(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_ptr dA[],    magma_minproduct_int_t ldda,  magma_minproduct_int_t offset,
    magma_minproductDouble_ptr dB[],    magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr dC[],    magma_minproduct_int_t lddc,
    magma_minproductDouble_ptr dwork[], magma_minproduct_int_t dworksiz,
    double    *C,       magma_minproduct_int_t ldc,
    double    *work[],  magma_minproduct_int_t worksiz,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][20], magma_minproduct_int_t nbevents );

void
magma_minproductblas_dsymm_mgpu_com(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_ptr dA[],    magma_minproduct_int_t ldda,  magma_minproduct_int_t offset,
    magma_minproductDouble_ptr dB[],    magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr dC[],    magma_minproduct_int_t lddc,
    magma_minproductDouble_ptr dwork[], magma_minproduct_int_t dworksiz,
    double    *C,       magma_minproduct_int_t ldc,
    double    *work[],  magma_minproduct_int_t worksiz,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents,
    magma_minproduct_int_t gnode[Magma_minproductMaxGPUs][Magma_minproductMaxGPUs+2], magma_minproduct_int_t nbcmplx );

void
magma_minproductblas_dsymm_mgpu_spec(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_ptr dA[],    magma_minproduct_int_t ldda,  magma_minproduct_int_t offset,
    magma_minproductDouble_ptr dB[],    magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr dC[],    magma_minproduct_int_t lddc,
    magma_minproductDouble_ptr dwork[], magma_minproduct_int_t dworksiz,
    double    *C,       magma_minproduct_int_t ldc,
    double    *work[],  magma_minproduct_int_t worksiz,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents,
    magma_minproduct_int_t gnode[Magma_minproductMaxGPUs][Magma_minproductMaxGPUs+2], magma_minproduct_int_t nbcmplx );

void
magma_minproductblas_dsymm_mgpu_spec33(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_ptr dA[],    magma_minproduct_int_t ldda,  magma_minproduct_int_t offset,
    magma_minproductDouble_ptr dB[],    magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr dC[],    magma_minproduct_int_t lddc,
    magma_minproductDouble_ptr dVIN[],  magma_minproduct_int_t lddv, magma_minproduct_int_t v_offset,
    magma_minproductDouble_ptr dwork[], magma_minproduct_int_t dworksiz,
    double *C,       magma_minproduct_int_t ldc,
    double *work[],  magma_minproduct_int_t worksiz,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents,
    magma_minproduct_int_t gnode[Magma_minproductMaxGPUs][Magma_minproductMaxGPUs+2], magma_minproduct_int_t nbcmplx );

magma_minproduct_int_t
magma_minproductblas_dsymv_mgpu(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr const d_lA[], magma_minproduct_int_t ldda, magma_minproduct_int_t offset,
    double const *x,         magma_minproduct_int_t incx,
    double beta,             
    double       *y,         magma_minproduct_int_t incy,
    double       *hwork,     magma_minproduct_int_t lhwork,
    magma_minproductDouble_ptr    dwork[],   magma_minproduct_int_t ldwork,
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[] );

magma_minproduct_int_t
magma_minproductblas_dsymv_mgpu_sync(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr const d_lA[], magma_minproduct_int_t ldda, magma_minproduct_int_t offset,
    double const *x,         magma_minproduct_int_t incx,
    double beta,             
    double       *y,         magma_minproduct_int_t incy,
    double       *hwork,     magma_minproduct_int_t lhwork,
    magma_minproductDouble_ptr    dwork[],   magma_minproduct_int_t ldwork,
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[] );

// Ichi's version, in src/dsytrd_mgpu.cpp
void
magma_minproduct_dsyr2k_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t nb, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    double beta,
    magma_minproductDouble_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t nqueue, magma_minproduct_queue_t queues[][10] );

void
magma_minproductblas_dsyr2k_mgpu_spec(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_ptr dA[], magma_minproduct_int_t ldda, magma_minproduct_int_t a_offset,
    magma_minproductDouble_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    double beta,
    magma_minproductDouble_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb, magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue );

void
magma_minproductblas_dsyr2k_mgpu_spec324(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_ptr dVIN[], magma_minproduct_int_t lddv, magma_minproduct_int_t v_offset,
    magma_minproductDouble_ptr dWIN[], magma_minproduct_int_t lddw, magma_minproduct_int_t w_offset,
    double beta,
    magma_minproductDouble_ptr dC[],   magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproductDouble_ptr dwork[], magma_minproduct_int_t lndwork,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents );

void
magma_minproductblas_dsyr2k_mgpu_spec325(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_ptr dVIN[], magma_minproduct_int_t lddv, magma_minproduct_int_t v_offset,
    magma_minproductDouble_ptr dWIN[], magma_minproduct_int_t lddw, magma_minproduct_int_t w_offset,
    double beta,
    magma_minproductDouble_ptr dC[],   magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproductDouble_ptr dwork[], magma_minproduct_int_t lndwork,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    double *harray[],
    magma_minproductDouble_ptr *darray[],
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents );

void
magma_minproductblas_dsyr2k_mgpu2(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_ptr dA[], magma_minproduct_int_t ldda, magma_minproduct_int_t a_offset,
    magma_minproductDouble_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    double beta,
    magma_minproductDouble_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue );

// in src/dpotrf_mgpu_right.cpp
void
magma_minproduct_dsyrk_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t nb, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    double beta,
    magma_minproductDouble_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t nqueue, magma_minproduct_queue_t queues[][10]);

// in src/dpotrf_mgpu_right.cpp
void
magma_minproduct_dsyrk_mgpu2(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t nb, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    double beta,
    magma_minproductDouble_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t nqueue, magma_minproduct_queue_t queues[][10]);


  /*
   * LAPACK auxiliary functions (alphabetical order)
   */
void
magma_minproductblas_dgeadd(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dB, magma_minproduct_int_t lddb );

void
magma_minproductblas_dlacpy(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dB, magma_minproduct_int_t lddb );

void
magma_minproductblas_dlacpy_cnjg(
    magma_minproduct_int_t n, double *dA1, magma_minproduct_int_t lda1,
    double *dA2, magma_minproduct_int_t lda2);

double
magma_minproductblas_dlange(
    magma_minproduct_norm_t norm,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dwork );

double
magma_minproductblas_dlansy(
    magma_minproduct_norm_t norm, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dwork );

double
magma_minproductblas_dlansy(
    magma_minproduct_norm_t norm, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dwork );

void
magma_minproductblas_dlarfg(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dalpha, magma_minproductDouble_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr dtau );

void
magma_minproductblas_dlarfg_work(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dalpha, magma_minproductDouble_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr dtau, magma_minproductDouble_ptr dwork );

void
magma_minproductblas_dlascl(
    magma_minproduct_type_t type, magma_minproduct_int_t kl, magma_minproduct_int_t ku,
    double cfrom, double cto,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_dlascl_2x2(
    magma_minproduct_type_t type, magma_minproduct_int_t m,
    double *dW, magma_minproduct_int_t lddw,
    double *dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_dlascl2(
    magma_minproduct_type_t type,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dD,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_dlascl_diag(
    magma_minproduct_type_t type, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dD, magma_minproduct_int_t lddd,
          magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_dlaset(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double offdiag, double diag,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda );

void
magma_minproductblas_dlaset_band(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double offdiag, double diag,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda);

void
magma_minproductblas_dlaswp(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dAT, magma_minproduct_int_t ldda,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci );

void
magma_minproductblas_dlaswp2(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dAT, magma_minproduct_int_t ldda,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    magma_minproductInt_const_ptr d_ipiv, magma_minproduct_int_t inci );

void
magma_minproductblas_dlaswpx(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldx, magma_minproduct_int_t ldy,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci );

void
magma_minproductblas_dsymmetrize(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda );

void
magma_minproductblas_dsymmetrize_tiles(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t ntile, magma_minproduct_int_t mstride, magma_minproduct_int_t nstride );

void
magma_minproductblas_dtrtri_diag(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr d_dinvA);

  /*
   * to cleanup (alphabetical order)
   */
void
magma_minproductblas_dnrm2_adjust(
    magma_minproduct_int_t k,
    magma_minproductDouble_ptr dxnorm,
    magma_minproductDouble_ptr dc);

void
magma_minproductblas_dnrm2_check(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dxnorm,
    magma_minproductDouble_ptr dlsticc);

void
magma_minproductblas_dnrm2_cols(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dxnorm);

void
magma_minproductblas_dnrm2_row_check_adjust(
    magma_minproduct_int_t k, double tol,
    magma_minproductDouble_ptr dxnorm,
    magma_minproductDouble_ptr dxnorm2,
    magma_minproductDouble_ptr dC, magma_minproduct_int_t lddc,
    magma_minproductDouble_ptr dlsticc);

void
magma_minproduct_dlarfbx_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t k,
    magma_minproductDouble_ptr V,  magma_minproduct_int_t ldv,
    magma_minproductDouble_ptr dT, magma_minproduct_int_t ldt,
    magma_minproductDouble_ptr c,
    magma_minproductDouble_ptr dwork);

void
magma_minproduct_dlarfg_gpu(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dx0,
    magma_minproductDouble_ptr dx,
    magma_minproductDouble_ptr dtau,
    magma_minproductDouble_ptr        dxnorm,
    magma_minproductDouble_ptr dAkk );

void
magma_minproduct_dlarfgtx_gpu(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dx0,
    magma_minproductDouble_ptr dx,
    magma_minproductDouble_ptr dtau,
    magma_minproductDouble_ptr        dxnorm,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t iter,
    magma_minproductDouble_ptr V,  magma_minproduct_int_t ldv,
    magma_minproductDouble_ptr T,  magma_minproduct_int_t ldt,
    magma_minproductDouble_ptr dwork);

void
magma_minproduct_dlarfgx_gpu(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dx0,
    magma_minproductDouble_ptr dx,
    magma_minproductDouble_ptr dtau,
    magma_minproductDouble_ptr        dxnorm,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t iter);

void
magma_minproduct_dlarfx_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr v,
    magma_minproductDouble_ptr tau,
    magma_minproductDouble_ptr C,  magma_minproduct_int_t ldc,
    magma_minproductDouble_ptr        xnorm,
    magma_minproductDouble_ptr dT, magma_minproduct_int_t iter,
    magma_minproductDouble_ptr work);


  /*
   * Level 1 BLAS (alphabetical order)
   */
void
magma_minproductblas_dswap(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr dy, magma_minproduct_int_t incy );

void
magma_minproductblas_dswapblk(
    magma_minproduct_order_t order,
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t i1, magma_minproduct_int_t i2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci,
    magma_minproduct_int_t offset );

void
magma_minproductblas_dswapdblk(
    magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t inca,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb, magma_minproduct_int_t incb );

  /*
   * Level 2 BLAS (alphabetical order)
   */
void
magma_minproductblas_dgemv(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    double beta,
    magma_minproductDouble_ptr       dy, magma_minproduct_int_t incy );


void
magma_minproductblas_dgemv_conjv(
    magma_minproduct_int_t m, magma_minproduct_int_t n, double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    double beta,
    magma_minproductDouble_ptr dy, magma_minproduct_int_t incy);

magma_minproduct_int_t
magma_minproductblas_dsymv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    double beta,
    magma_minproductDouble_ptr       dy, magma_minproduct_int_t incy );

magma_minproduct_int_t
magma_minproductblas_dsymv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    double beta,
    magma_minproductDouble_ptr       dy, magma_minproduct_int_t incy );

magma_minproduct_int_t
magma_minproductblas_dsymv_work(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    double beta,
    magma_minproductDouble_ptr       dy, magma_minproduct_int_t incy,
    magma_minproductDouble_ptr       dwork, magma_minproduct_int_t lwork,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproductblas_dsymv_work(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    double beta,
    magma_minproductDouble_ptr       dy, magma_minproduct_int_t incy,
    magma_minproductDouble_ptr       dwork, magma_minproduct_int_t lwork,
    magma_minproduct_queue_t queue );

  /*
   * Level 3 BLAS (alphabetical order)
   */
void
magma_minproductblas_dgemm(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dB, magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_dgemm_reduce(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dB, magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_dsymm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dB, magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_dsymm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dB, magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_dsyr2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dB, magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_dsyr2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dB, magma_minproduct_int_t lddb,
    double  beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_dsyrk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_dsyrk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    double  alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    double  beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_dtrsm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dB, magma_minproduct_int_t lddb );

void
magma_minproductblas_dtrsm_outofplace(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t flag, magma_minproductDouble_ptr d_dinvA, magma_minproductDouble_ptr dX );

void
magma_minproductblas_dtrsm_work(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t flag, magma_minproductDouble_ptr d_dinvA, magma_minproductDouble_ptr dX );


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

#define magma_minproduct_dsetvector(           n, hx_src, incx, dy_dst, incy ) \
        magma_minproduct_dsetvector_internal(  n, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_minproduct_dgetvector(           n, dx_src, incx, hy_dst, incy ) \
        magma_minproduct_dgetvector_internal(  n, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_minproduct_dcopyvector(          n, dx_src, incx, dy_dst, incy ) \
        magma_minproduct_dcopyvector_internal( n, dx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_minproduct_dsetvector_async(           n, hx_src, incx, dy_dst, incy, queue ) \
        magma_minproduct_dsetvector_async_internal(  n, hx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_dgetvector_async(           n, dx_src, incx, hy_dst, incy, queue ) \
        magma_minproduct_dgetvector_async_internal(  n, dx_src, incx, hy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_dcopyvector_async(          n, dx_src, incx, dy_dst, incy, queue ) \
        magma_minproduct_dcopyvector_async_internal( n, dx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

void
magma_minproduct_dsetvector_internal(
    magma_minproduct_int_t n,
    double const    *hx_src, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr       dy_dst, magma_minproduct_int_t incy,
    const char* func, const char* file, int line );

void
magma_minproduct_dgetvector_internal(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx_src, magma_minproduct_int_t incx,
    double          *hy_dst, magma_minproduct_int_t incy,
    const char* func, const char* file, int line );

void
magma_minproduct_dcopyvector_internal(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx_src, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr       dy_dst, magma_minproduct_int_t incy,
    const char* func, const char* file, int line );

void
magma_minproduct_dsetvector_async_internal(
    magma_minproduct_int_t n,
    double const    *hx_src, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr       dy_dst, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void
magma_minproduct_dgetvector_async_internal(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx_src, magma_minproduct_int_t incx,
    double          *hy_dst, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void
magma_minproduct_dcopyvector_async_internal(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx_src, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr       dy_dst, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// copying sub-matrices (contiguous columns)
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

#define magma_minproduct_dsetmatrix(           m, n, hA_src, lda, dB_dst, lddb ) \
        magma_minproduct_dsetmatrix_internal(  m, n, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_minproduct_dgetmatrix(           m, n, dA_src, ldda, hB_dst, ldb ) \
        magma_minproduct_dgetmatrix_internal(  m, n, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )

#define magma_minproduct_dcopymatrix(          m, n, dA_src, ldda, dB_dst, lddb ) \
        magma_minproduct_dcopymatrix_internal( m, n, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_minproduct_dsetmatrix_async(           m, n, hA_src, lda, dB_dst, lddb, queue ) \
        magma_minproduct_dsetmatrix_async_internal(  m, n, hA_src, lda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_dgetmatrix_async(           m, n, dA_src, ldda, hB_dst, ldb, queue ) \
        magma_minproduct_dgetmatrix_async_internal(  m, n, dA_src, ldda, hB_dst, ldb, queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_dcopymatrix_async(          m, n, dA_src, ldda, dB_dst, lddb, queue ) \
        magma_minproduct_dcopymatrix_async_internal( m, n, dA_src, ldda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

void
magma_minproduct_dsetmatrix_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double const    *hA_src, magma_minproduct_int_t ldha,
    magma_minproductDouble_ptr       dB_dst, magma_minproduct_int_t lddb,
    const char* func, const char* file, int line );

void
magma_minproduct_dgetmatrix_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA_src, magma_minproduct_int_t ldda,
    double          *hB_dst, magma_minproduct_int_t ldhb,
    const char* func, const char* file, int line );

void
magma_minproduct_dcopymatrix_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA_src, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dB_dst, magma_minproduct_int_t lddb,
    const char* func, const char* file, int line );

void
magma_minproduct_dsetmatrix_async_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double const    *hA_src, magma_minproduct_int_t ldha,
    magma_minproductDouble_ptr       dB_dst, magma_minproduct_int_t lddb,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void
magma_minproduct_dgetmatrix_async_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA_src, magma_minproduct_int_t ldda,
    double          *hB_dst, magma_minproduct_int_t ldhb,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void
magma_minproduct_dcopymatrix_async_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA_src, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dB_dst, magma_minproduct_int_t lddb,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// Level 1 BLAS (alphabetical order)

// in cublas_v2, result returned through output argument
magma_minproduct_int_t
magma_minproduct_idamax(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx );

// in cublas_v2, result returned through output argument
magma_minproduct_int_t
magma_minproduct_idamin(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx );

// in cublas_v2, result returned through output argument
double
magma_minproduct_dasum(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx );

void
magma_minproduct_daxpy(
    magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr       dy, magma_minproduct_int_t incy );

void
magma_minproduct_dcopy(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr       dy, magma_minproduct_int_t incy );

// in cublas_v2, result returned through output argument
double
magma_minproduct_ddot(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_const_ptr dy, magma_minproduct_int_t incy );

// in cublas_v2, result returned through output argument
double
magma_minproduct_ddot(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_const_ptr dy, magma_minproduct_int_t incy );

// in cublas_v2, result returned through output argument
double
magma_minproduct_dnrm2(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx );

void
magma_minproduct_drot(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr dy, magma_minproduct_int_t incy,
    double dc, double ds );

void
magma_minproduct_drot(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr dy, magma_minproduct_int_t incy,
    double dc, double ds );

#ifdef REAL
void
magma_minproduct_drotm(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr dy, magma_minproduct_int_t incy,
    magma_minproductDouble_const_ptr param );

void
magma_minproduct_drotmg(
    magma_minproductDouble_ptr d1, magma_minproductDouble_ptr       d2,
    magma_minproductDouble_ptr x1, magma_minproductDouble_const_ptr y1,
    magma_minproductDouble_ptr param );
#endif

void
magma_minproduct_dscal(
    magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_ptr dx, magma_minproduct_int_t incx );

void
magma_minproduct_dscal(
    magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_ptr dx, magma_minproduct_int_t incx );

void
magma_minproduct_dswap(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr dy, magma_minproduct_int_t incy );

// ========================================
// Level 2 BLAS (alphabetical order)

void
magma_minproduct_dgemv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    double beta,
    magma_minproductDouble_ptr       dy, magma_minproduct_int_t incy );

void
magma_minproduct_dger(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_const_ptr dy, magma_minproduct_int_t incy,
    magma_minproductDouble_ptr       dA, magma_minproduct_int_t ldda );

void
magma_minproduct_dger(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_const_ptr dy, magma_minproduct_int_t incy,
    magma_minproductDouble_ptr       dA, magma_minproduct_int_t ldda );

void
magma_minproduct_dsymv(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    double beta,
    magma_minproductDouble_ptr       dy, magma_minproduct_int_t incy );

void
magma_minproduct_dsyr(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr       dA, magma_minproduct_int_t ldda );

void
magma_minproduct_dsyr2(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_const_ptr dy, magma_minproduct_int_t incy,
    magma_minproductDouble_ptr       dA, magma_minproduct_int_t ldda );

void
magma_minproduct_dtrmv(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dx, magma_minproduct_int_t incx );

void
magma_minproduct_dtrsv(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dx, magma_minproduct_int_t incx );

// ========================================
// Level 3 BLAS (alphabetical order)

void
magma_minproduct_dgemm(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dB, magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_dsymm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dB, magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_dsymm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dB, magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_dsyr2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dB, magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_dsyr2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dB, magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_dsyrk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_dsyrk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_dtrmm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dB, magma_minproduct_int_t lddb );

void
magma_minproduct_dtrsm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dB, magma_minproduct_int_t lddb );


#ifdef __cplusplus
}
#endif

#undef REAL

#endif  /* MAGMA_minproductBLAS_D_H */
