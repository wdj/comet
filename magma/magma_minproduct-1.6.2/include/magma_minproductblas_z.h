/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
*/

#ifndef MAGMA_minproductBLAS_Z_H
#define MAGMA_minproductBLAS_Z_H

#include "magma_minproduct_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_minproductblas_ztranspose_inplace(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda );

void
magma_minproductblas_ztranspose(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr       dAT, magma_minproduct_int_t lddat );

void
magma_minproductblas_zgetmatrix_transpose(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dAT,   magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex          *hA,    magma_minproduct_int_t lda,
    magma_minproductDoubleComplex_ptr       dwork, magma_minproduct_int_t lddwork, magma_minproduct_int_t nb );

void
magma_minproductblas_zsetmatrix_transpose(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const magma_minproductDoubleComplex *hA,    magma_minproduct_int_t lda,
    magma_minproductDoubleComplex_ptr    dAT,   magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr    dwork, magma_minproduct_int_t lddwork, magma_minproduct_int_t nb );

  /*
   * RBT-related functions
   */
void
magma_minproductblas_zprbt(
    magma_minproduct_int_t n, 
    magma_minproductDoubleComplex *dA, magma_minproduct_int_t ldda, 
    magma_minproductDoubleComplex *du, magma_minproductDoubleComplex *dv);

void
magma_minproductblas_zprbt_mv(
    magma_minproduct_int_t n, 
    magma_minproductDoubleComplex *dv, magma_minproductDoubleComplex *db);

void
magma_minproductblas_zprbt_mtv(
    magma_minproduct_int_t n, 
    magma_minproductDoubleComplex *du, magma_minproductDoubleComplex *db);

void
magma_minproductblas_zaxpycp2(
    magma_minproduct_int_t m, magma_minproductDoubleComplex *r, magma_minproductDoubleComplex *x,
    const magma_minproductDoubleComplex *b);

  /*
   * Multi-GPU copy functions
   */
void
magma_minproductblas_zgetmatrix_transpose_mgpu(
    magma_minproduct_int_t ngpu, magma_minproduct_queue_t queues[][2],
    magma_minproductDoubleComplex_const_ptr const dAT[],    magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex                *hA,       magma_minproduct_int_t lda,
    magma_minproductDoubleComplex_ptr             dwork[],  magma_minproduct_int_t lddw,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb );

void
magma_minproductblas_zsetmatrix_transpose_mgpu(
    magma_minproduct_int_t ngpu, magma_minproduct_queue_t queues[][2],
    const magma_minproductDoubleComplex *hA,      magma_minproduct_int_t lda,
    magma_minproductDoubleComplex_ptr    dAT[],   magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr    dwork[], magma_minproduct_int_t lddw,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb );

void
magma_minproduct_zgetmatrix_1D_col_bcyclic(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr const dA[], magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex                *hA,   magma_minproduct_int_t lda,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb );

void
magma_minproduct_zsetmatrix_1D_col_bcyclic(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const magma_minproductDoubleComplex *hA,   magma_minproduct_int_t lda,
    magma_minproductDoubleComplex_ptr    dA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb );

void
magma_minproduct_zgetmatrix_1D_row_bcyclic(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr const dA[], magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex                *hA,   magma_minproduct_int_t lda,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb );

void
magma_minproduct_zsetmatrix_1D_row_bcyclic(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const magma_minproductDoubleComplex *hA,   magma_minproduct_int_t lda,
    magma_minproductDoubleComplex_ptr    dA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb );

// in src/zhetrd_mgpu.cpp
// TODO rename zsetmatrix_sy or similar
magma_minproduct_int_t
magma_minproduct_zhtodhe(
    magma_minproduct_int_t ngpu, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex     *A,   magma_minproduct_int_t lda,
    magma_minproductDoubleComplex_ptr dA[], magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queues[][10],
    magma_minproduct_int_t *info );

// in src/zpotrf3_mgpu.cpp
// TODO same as magma_minproduct_zhtodhe?
magma_minproduct_int_t
magma_minproduct_zhtodpo(
    magma_minproduct_int_t ngpu, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t off_i, magma_minproduct_int_t off_j, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex     *A,   magma_minproduct_int_t lda,
    magma_minproductDoubleComplex_ptr dA[], magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queues[][3],
    magma_minproduct_int_t *info );

// in src/zpotrf3_mgpu.cpp
// TODO rename zgetmatrix_sy or similar
magma_minproduct_int_t
magma_minproduct_zdtohpo(
    magma_minproduct_int_t ngpu, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t off_i, magma_minproduct_int_t off_j, magma_minproduct_int_t nb, magma_minproduct_int_t NB,
    magma_minproductDoubleComplex     *A,   magma_minproduct_int_t lda,
    magma_minproductDoubleComplex_ptr dA[], magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queues[][3],
    magma_minproduct_int_t *info );


  /*
   * Multi-GPU BLAS functions (alphabetical order)
   */
void
magma_minproductblas_zhemm_mgpu(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dA[],    magma_minproduct_int_t ldda,  magma_minproduct_int_t offset,
    magma_minproductDoubleComplex_ptr dB[],    magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dC[],    magma_minproduct_int_t lddc,
    magma_minproductDoubleComplex_ptr dwork[], magma_minproduct_int_t dworksiz,
    magma_minproductDoubleComplex    *C,       magma_minproduct_int_t ldc,
    magma_minproductDoubleComplex    *work[],  magma_minproduct_int_t worksiz,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][20], magma_minproduct_int_t nbevents );

void
magma_minproductblas_zhemm_mgpu_com(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dA[],    magma_minproduct_int_t ldda,  magma_minproduct_int_t offset,
    magma_minproductDoubleComplex_ptr dB[],    magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dC[],    magma_minproduct_int_t lddc,
    magma_minproductDoubleComplex_ptr dwork[], magma_minproduct_int_t dworksiz,
    magma_minproductDoubleComplex    *C,       magma_minproduct_int_t ldc,
    magma_minproductDoubleComplex    *work[],  magma_minproduct_int_t worksiz,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents,
    magma_minproduct_int_t gnode[Magma_minproductMaxGPUs][Magma_minproductMaxGPUs+2], magma_minproduct_int_t nbcmplx );

void
magma_minproductblas_zhemm_mgpu_spec(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dA[],    magma_minproduct_int_t ldda,  magma_minproduct_int_t offset,
    magma_minproductDoubleComplex_ptr dB[],    magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dC[],    magma_minproduct_int_t lddc,
    magma_minproductDoubleComplex_ptr dwork[], magma_minproduct_int_t dworksiz,
    magma_minproductDoubleComplex    *C,       magma_minproduct_int_t ldc,
    magma_minproductDoubleComplex    *work[],  magma_minproduct_int_t worksiz,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents,
    magma_minproduct_int_t gnode[Magma_minproductMaxGPUs][Magma_minproductMaxGPUs+2], magma_minproduct_int_t nbcmplx );

void
magma_minproductblas_zhemm_mgpu_spec33(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dA[],    magma_minproduct_int_t ldda,  magma_minproduct_int_t offset,
    magma_minproductDoubleComplex_ptr dB[],    magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dC[],    magma_minproduct_int_t lddc,
    magma_minproductDoubleComplex_ptr dVIN[],  magma_minproduct_int_t lddv, magma_minproduct_int_t v_offset,
    magma_minproductDoubleComplex_ptr dwork[], magma_minproduct_int_t dworksiz,
    magma_minproductDoubleComplex *C,       magma_minproduct_int_t ldc,
    magma_minproductDoubleComplex *work[],  magma_minproduct_int_t worksiz,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents,
    magma_minproduct_int_t gnode[Magma_minproductMaxGPUs][Magma_minproductMaxGPUs+2], magma_minproduct_int_t nbcmplx );

magma_minproduct_int_t
magma_minproductblas_zhemv_mgpu(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr const d_lA[], magma_minproduct_int_t ldda, magma_minproduct_int_t offset,
    magma_minproductDoubleComplex const *x,         magma_minproduct_int_t incx,
    magma_minproductDoubleComplex beta,             
    magma_minproductDoubleComplex       *y,         magma_minproduct_int_t incy,
    magma_minproductDoubleComplex       *hwork,     magma_minproduct_int_t lhwork,
    magma_minproductDoubleComplex_ptr    dwork[],   magma_minproduct_int_t ldwork,
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[] );

magma_minproduct_int_t
magma_minproductblas_zhemv_mgpu_sync(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr const d_lA[], magma_minproduct_int_t ldda, magma_minproduct_int_t offset,
    magma_minproductDoubleComplex const *x,         magma_minproduct_int_t incx,
    magma_minproductDoubleComplex beta,             
    magma_minproductDoubleComplex       *y,         magma_minproduct_int_t incy,
    magma_minproductDoubleComplex       *hwork,     magma_minproduct_int_t lhwork,
    magma_minproductDoubleComplex_ptr    dwork[],   magma_minproduct_int_t ldwork,
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[] );

// Ichi's version, in src/zhetrd_mgpu.cpp
void
magma_minproduct_zher2k_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t nb, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    double beta,
    magma_minproductDoubleComplex_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t nqueue, magma_minproduct_queue_t queues[][10] );

void
magma_minproductblas_zher2k_mgpu_spec(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dA[], magma_minproduct_int_t ldda, magma_minproduct_int_t a_offset,
    magma_minproductDoubleComplex_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    double beta,
    magma_minproductDoubleComplex_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb, magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue );

void
magma_minproductblas_zher2k_mgpu_spec324(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dVIN[], magma_minproduct_int_t lddv, magma_minproduct_int_t v_offset,
    magma_minproductDoubleComplex_ptr dWIN[], magma_minproduct_int_t lddw, magma_minproduct_int_t w_offset,
    double beta,
    magma_minproductDoubleComplex_ptr dC[],   magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproductDoubleComplex_ptr dwork[], magma_minproduct_int_t lndwork,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents );

void
magma_minproductblas_zher2k_mgpu_spec325(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dVIN[], magma_minproduct_int_t lddv, magma_minproduct_int_t v_offset,
    magma_minproductDoubleComplex_ptr dWIN[], magma_minproduct_int_t lddw, magma_minproduct_int_t w_offset,
    double beta,
    magma_minproductDoubleComplex_ptr dC[],   magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproductDoubleComplex_ptr dwork[], magma_minproduct_int_t lndwork,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex *harray[],
    magma_minproductDoubleComplex_ptr *darray[],
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents );

void
magma_minproductblas_zher2k_mgpu2(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dA[], magma_minproduct_int_t ldda, magma_minproduct_int_t a_offset,
    magma_minproductDoubleComplex_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    double beta,
    magma_minproductDoubleComplex_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue );

// in src/zpotrf_mgpu_right.cpp
void
magma_minproduct_zherk_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t nb, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDoubleComplex_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    double beta,
    magma_minproductDoubleComplex_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t nqueue, magma_minproduct_queue_t queues[][10]);

// in src/zpotrf_mgpu_right.cpp
void
magma_minproduct_zherk_mgpu2(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t nb, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDoubleComplex_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    double beta,
    magma_minproductDoubleComplex_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t nqueue, magma_minproduct_queue_t queues[][10]);


  /*
   * LAPACK auxiliary functions (alphabetical order)
   */
void
magma_minproductblas_zgeadd(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr       dB, magma_minproduct_int_t lddb );

void
magma_minproductblas_zlacpy(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr       dB, magma_minproduct_int_t lddb );

void
magma_minproductblas_zlacpy_cnjg(
    magma_minproduct_int_t n, magma_minproductDoubleComplex *dA1, magma_minproduct_int_t lda1,
    magma_minproductDoubleComplex *dA2, magma_minproduct_int_t lda2);

double
magma_minproductblas_zlange(
    magma_minproduct_norm_t norm,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dwork );

double
magma_minproductblas_zlanhe(
    magma_minproduct_norm_t norm, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dwork );

double
magma_minproductblas_zlansy(
    magma_minproduct_norm_t norm, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dwork );

void
magma_minproductblas_zlarfg(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dalpha, magma_minproductDoubleComplex_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_ptr dtau );

void
magma_minproductblas_zlarfg_work(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dalpha, magma_minproductDoubleComplex_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_ptr dtau, magma_minproductDoubleComplex_ptr dwork );

void
magma_minproductblas_zlascl(
    magma_minproduct_type_t type, magma_minproduct_int_t kl, magma_minproduct_int_t ku,
    double cfrom, double cto,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_zlascl_2x2(
    magma_minproduct_type_t type, magma_minproduct_int_t m,
    magma_minproductDoubleComplex *dW, magma_minproduct_int_t lddw,
    magma_minproductDoubleComplex *dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_zlascl2(
    magma_minproduct_type_t type,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dD,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_zlascl_diag(
    magma_minproduct_type_t type, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dD, magma_minproduct_int_t lddd,
          magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_zlaset(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex offdiag, magma_minproductDoubleComplex diag,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda );

void
magma_minproductblas_zlaset_band(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex offdiag, magma_minproductDoubleComplex diag,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda);

void
magma_minproductblas_zlaswp(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dAT, magma_minproduct_int_t ldda,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci );

void
magma_minproductblas_zlaswp2(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dAT, magma_minproduct_int_t ldda,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    magma_minproductInt_const_ptr d_ipiv, magma_minproduct_int_t inci );

void
magma_minproductblas_zlaswpx(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldx, magma_minproduct_int_t ldy,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci );

void
magma_minproductblas_zsymmetrize(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda );

void
magma_minproductblas_zsymmetrize_tiles(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t ntile, magma_minproduct_int_t mstride, magma_minproduct_int_t nstride );

void
magma_minproductblas_ztrtri_diag(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr d_dinvA);

  /*
   * to cleanup (alphabetical order)
   */
void
magma_minproductblas_dznrm2_adjust(
    magma_minproduct_int_t k,
    magma_minproductDouble_ptr dxnorm,
    magma_minproductDoubleComplex_ptr dc);

void
magma_minproductblas_dznrm2_check(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dxnorm,
    magma_minproductDouble_ptr dlsticc);

void
magma_minproductblas_dznrm2_cols(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dxnorm);

void
magma_minproductblas_dznrm2_row_check_adjust(
    magma_minproduct_int_t k, double tol,
    magma_minproductDouble_ptr dxnorm,
    magma_minproductDouble_ptr dxnorm2,
    magma_minproductDoubleComplex_ptr dC, magma_minproduct_int_t lddc,
    magma_minproductDouble_ptr dlsticc);

void
magma_minproduct_zlarfbx_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t k,
    magma_minproductDoubleComplex_ptr V,  magma_minproduct_int_t ldv,
    magma_minproductDoubleComplex_ptr dT, magma_minproduct_int_t ldt,
    magma_minproductDoubleComplex_ptr c,
    magma_minproductDoubleComplex_ptr dwork);

void
magma_minproduct_zlarfg_gpu(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dx0,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex_ptr dtau,
    magma_minproductDouble_ptr        dxnorm,
    magma_minproductDoubleComplex_ptr dAkk );

void
magma_minproduct_zlarfgtx_gpu(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dx0,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex_ptr dtau,
    magma_minproductDouble_ptr        dxnorm,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t iter,
    magma_minproductDoubleComplex_ptr V,  magma_minproduct_int_t ldv,
    magma_minproductDoubleComplex_ptr T,  magma_minproduct_int_t ldt,
    magma_minproductDoubleComplex_ptr dwork);

void
magma_minproduct_zlarfgx_gpu(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dx0,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex_ptr dtau,
    magma_minproductDouble_ptr        dxnorm,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t iter);

void
magma_minproduct_zlarfx_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr v,
    magma_minproductDoubleComplex_ptr tau,
    magma_minproductDoubleComplex_ptr C,  magma_minproduct_int_t ldc,
    magma_minproductDouble_ptr        xnorm,
    magma_minproductDoubleComplex_ptr dT, magma_minproduct_int_t iter,
    magma_minproductDoubleComplex_ptr work);


  /*
   * Level 1 BLAS (alphabetical order)
   */
void
magma_minproductblas_zswap(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_ptr dy, magma_minproduct_int_t incy );

void
magma_minproductblas_zswapblk(
    magma_minproduct_order_t order,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t i1, magma_minproduct_int_t i2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci,
    magma_minproduct_int_t offset );

void
magma_minproductblas_zswapdblk(
    magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t inca,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb, magma_minproduct_int_t incb );

  /*
   * Level 2 BLAS (alphabetical order)
   */
void
magma_minproductblas_zgemv(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dy, magma_minproduct_int_t incy );


void
magma_minproductblas_zgemv_conjv(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dy, magma_minproduct_int_t incy);

magma_minproduct_int_t
magma_minproductblas_zhemv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dy, magma_minproduct_int_t incy );

magma_minproduct_int_t
magma_minproductblas_zsymv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dy, magma_minproduct_int_t incy );

magma_minproduct_int_t
magma_minproductblas_zhemv_work(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dy, magma_minproduct_int_t incy,
    magma_minproductDoubleComplex_ptr       dwork, magma_minproduct_int_t lwork,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproductblas_zsymv_work(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dy, magma_minproduct_int_t incy,
    magma_minproductDoubleComplex_ptr       dwork, magma_minproduct_int_t lwork,
    magma_minproduct_queue_t queue );

  /*
   * Level 3 BLAS (alphabetical order)
   */
void
magma_minproductblas_zgemm(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_zgemm_reduce(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_zhemm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_zsymm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_zsyr2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_zher2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dB, magma_minproduct_int_t lddb,
    double  beta,
    magma_minproductDoubleComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_zsyrk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_zherk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    double  alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    double  beta,
    magma_minproductDoubleComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_ztrsm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr       dB, magma_minproduct_int_t lddb );

void
magma_minproductblas_ztrsm_outofplace(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr       dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t flag, magma_minproductDoubleComplex_ptr d_dinvA, magma_minproductDoubleComplex_ptr dX );

void
magma_minproductblas_ztrsm_work(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr       dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t flag, magma_minproductDoubleComplex_ptr d_dinvA, magma_minproductDoubleComplex_ptr dX );


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

#define magma_minproduct_zsetvector(           n, hx_src, incx, dy_dst, incy ) \
        magma_minproduct_zsetvector_internal(  n, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_minproduct_zgetvector(           n, dx_src, incx, hy_dst, incy ) \
        magma_minproduct_zgetvector_internal(  n, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_minproduct_zcopyvector(          n, dx_src, incx, dy_dst, incy ) \
        magma_minproduct_zcopyvector_internal( n, dx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_minproduct_zsetvector_async(           n, hx_src, incx, dy_dst, incy, queue ) \
        magma_minproduct_zsetvector_async_internal(  n, hx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_zgetvector_async(           n, dx_src, incx, hy_dst, incy, queue ) \
        magma_minproduct_zgetvector_async_internal(  n, dx_src, incx, hy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_zcopyvector_async(          n, dx_src, incx, dy_dst, incy, queue ) \
        magma_minproduct_zcopyvector_async_internal( n, dx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

void
magma_minproduct_zsetvector_internal(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex const    *hx_src, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_ptr       dy_dst, magma_minproduct_int_t incy,
    const char* func, const char* file, int line );

void
magma_minproduct_zgetvector_internal(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dx_src, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex          *hy_dst, magma_minproduct_int_t incy,
    const char* func, const char* file, int line );

void
magma_minproduct_zcopyvector_internal(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dx_src, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_ptr       dy_dst, magma_minproduct_int_t incy,
    const char* func, const char* file, int line );

void
magma_minproduct_zsetvector_async_internal(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex const    *hx_src, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_ptr       dy_dst, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void
magma_minproduct_zgetvector_async_internal(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dx_src, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex          *hy_dst, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void
magma_minproduct_zcopyvector_async_internal(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dx_src, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_ptr       dy_dst, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// copying sub-matrices (contiguous columns)
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

#define magma_minproduct_zsetmatrix(           m, n, hA_src, lda, dB_dst, lddb ) \
        magma_minproduct_zsetmatrix_internal(  m, n, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_minproduct_zgetmatrix(           m, n, dA_src, ldda, hB_dst, ldb ) \
        magma_minproduct_zgetmatrix_internal(  m, n, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )

#define magma_minproduct_zcopymatrix(          m, n, dA_src, ldda, dB_dst, lddb ) \
        magma_minproduct_zcopymatrix_internal( m, n, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_minproduct_zsetmatrix_async(           m, n, hA_src, lda, dB_dst, lddb, queue ) \
        magma_minproduct_zsetmatrix_async_internal(  m, n, hA_src, lda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_zgetmatrix_async(           m, n, dA_src, ldda, hB_dst, ldb, queue ) \
        magma_minproduct_zgetmatrix_async_internal(  m, n, dA_src, ldda, hB_dst, ldb, queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_zcopymatrix_async(          m, n, dA_src, ldda, dB_dst, lddb, queue ) \
        magma_minproduct_zcopymatrix_async_internal( m, n, dA_src, ldda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

void
magma_minproduct_zsetmatrix_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex const    *hA_src, magma_minproduct_int_t ldha,
    magma_minproductDoubleComplex_ptr       dB_dst, magma_minproduct_int_t lddb,
    const char* func, const char* file, int line );

void
magma_minproduct_zgetmatrix_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA_src, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex          *hB_dst, magma_minproduct_int_t ldhb,
    const char* func, const char* file, int line );

void
magma_minproduct_zcopymatrix_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA_src, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr       dB_dst, magma_minproduct_int_t lddb,
    const char* func, const char* file, int line );

void
magma_minproduct_zsetmatrix_async_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex const    *hA_src, magma_minproduct_int_t ldha,
    magma_minproductDoubleComplex_ptr       dB_dst, magma_minproduct_int_t lddb,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void
magma_minproduct_zgetmatrix_async_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA_src, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex          *hB_dst, magma_minproduct_int_t ldhb,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void
magma_minproduct_zcopymatrix_async_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA_src, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr       dB_dst, magma_minproduct_int_t lddb,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// Level 1 BLAS (alphabetical order)

// in cublas_v2, result returned through output argument
magma_minproduct_int_t
magma_minproduct_izamax(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx );

// in cublas_v2, result returned through output argument
magma_minproduct_int_t
magma_minproduct_izamin(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx );

// in cublas_v2, result returned through output argument
double
magma_minproduct_dzasum(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx );

void
magma_minproduct_zaxpy(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_ptr       dy, magma_minproduct_int_t incy );

void
magma_minproduct_zcopy(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_ptr       dy, magma_minproduct_int_t incy );

// in cublas_v2, result returned through output argument
magma_minproductDoubleComplex
magma_minproduct_zdotc(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_const_ptr dy, magma_minproduct_int_t incy );

// in cublas_v2, result returned through output argument
magma_minproductDoubleComplex
magma_minproduct_zdotu(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_const_ptr dy, magma_minproduct_int_t incy );

// in cublas_v2, result returned through output argument
double
magma_minproduct_dznrm2(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx );

void
magma_minproduct_zrot(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_ptr dy, magma_minproduct_int_t incy,
    double dc, magma_minproductDoubleComplex ds );

void
magma_minproduct_zdrot(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_ptr dy, magma_minproduct_int_t incy,
    double dc, double ds );

#ifdef REAL
void
magma_minproduct_zrotm(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr dy, magma_minproduct_int_t incy,
    magma_minproductDouble_const_ptr param );

void
magma_minproduct_zrotmg(
    magma_minproductDouble_ptr d1, magma_minproductDouble_ptr       d2,
    magma_minproductDouble_ptr x1, magma_minproductDouble_const_ptr y1,
    magma_minproductDouble_ptr param );
#endif

void
magma_minproduct_zscal(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dx, magma_minproduct_int_t incx );

void
magma_minproduct_zdscal(
    magma_minproduct_int_t n,
    double alpha,
    magma_minproductDoubleComplex_ptr dx, magma_minproduct_int_t incx );

void
magma_minproduct_zswap(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_ptr dy, magma_minproduct_int_t incy );

// ========================================
// Level 2 BLAS (alphabetical order)

void
magma_minproduct_zgemv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dy, magma_minproduct_int_t incy );

void
magma_minproduct_zgerc(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_const_ptr dy, magma_minproduct_int_t incy,
    magma_minproductDoubleComplex_ptr       dA, magma_minproduct_int_t ldda );

void
magma_minproduct_zgeru(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_const_ptr dy, magma_minproduct_int_t incy,
    magma_minproductDoubleComplex_ptr       dA, magma_minproduct_int_t ldda );

void
magma_minproduct_zhemv(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dy, magma_minproduct_int_t incy );

void
magma_minproduct_zher(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    double alpha,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_ptr       dA, magma_minproduct_int_t ldda );

void
magma_minproduct_zher2(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_const_ptr dy, magma_minproduct_int_t incy,
    magma_minproductDoubleComplex_ptr       dA, magma_minproduct_int_t ldda );

void
magma_minproduct_ztrmv(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr       dx, magma_minproduct_int_t incx );

void
magma_minproduct_ztrsv(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr       dx, magma_minproduct_int_t incx );

// ========================================
// Level 3 BLAS (alphabetical order)

void
magma_minproduct_zgemm(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_zsymm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_zhemm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_zsyr2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_zher2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_const_ptr dB, magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDoubleComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_zsyrk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_zherk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    double beta,
    magma_minproductDoubleComplex_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_ztrmm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr       dB, magma_minproduct_int_t lddb );

void
magma_minproduct_ztrsm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr       dB, magma_minproduct_int_t lddb );


#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif  /* MAGMA_minproductBLAS_Z_H */
