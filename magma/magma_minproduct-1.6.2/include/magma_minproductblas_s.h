/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproductblas_z.h normal z -> s, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_minproductBLAS_S_H
#define MAGMA_minproductBLAS_S_H

#include "magma_minproduct_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_minproductblas_stranspose_inplace(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda );

void
magma_minproductblas_stranspose(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr       dAT, magma_minproduct_int_t lddat );

void
magma_minproductblas_sgetmatrix_transpose(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dAT,   magma_minproduct_int_t ldda,
    float          *hA,    magma_minproduct_int_t lda,
    magma_minproductFloat_ptr       dwork, magma_minproduct_int_t lddwork, magma_minproduct_int_t nb );

void
magma_minproductblas_ssetmatrix_transpose(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const float *hA,    magma_minproduct_int_t lda,
    magma_minproductFloat_ptr    dAT,   magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr    dwork, magma_minproduct_int_t lddwork, magma_minproduct_int_t nb );

  /*
   * RBT-related functions
   */
void
magma_minproductblas_sprbt(
    magma_minproduct_int_t n, 
    float *dA, magma_minproduct_int_t ldda, 
    float *du, float *dv);

void
magma_minproductblas_sprbt_mv(
    magma_minproduct_int_t n, 
    float *dv, float *db);

void
magma_minproductblas_sprbt_mtv(
    magma_minproduct_int_t n, 
    float *du, float *db);

void
magma_minproductblas_saxpycp2(
    magma_minproduct_int_t m, float *r, float *x,
    const float *b);

  /*
   * Multi-GPU copy functions
   */
void
magma_minproductblas_sgetmatrix_transpose_mgpu(
    magma_minproduct_int_t ngpu, magma_minproduct_queue_t queues[][2],
    magma_minproductFloat_const_ptr const dAT[],    magma_minproduct_int_t ldda,
    float                *hA,       magma_minproduct_int_t lda,
    magma_minproductFloat_ptr             dwork[],  magma_minproduct_int_t lddw,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb );

void
magma_minproductblas_ssetmatrix_transpose_mgpu(
    magma_minproduct_int_t ngpu, magma_minproduct_queue_t queues[][2],
    const float *hA,      magma_minproduct_int_t lda,
    magma_minproductFloat_ptr    dAT[],   magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr    dwork[], magma_minproduct_int_t lddw,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb );

void
magma_minproduct_sgetmatrix_1D_col_bcyclic(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr const dA[], magma_minproduct_int_t ldda,
    float                *hA,   magma_minproduct_int_t lda,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb );

void
magma_minproduct_ssetmatrix_1D_col_bcyclic(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const float *hA,   magma_minproduct_int_t lda,
    magma_minproductFloat_ptr    dA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb );

void
magma_minproduct_sgetmatrix_1D_row_bcyclic(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr const dA[], magma_minproduct_int_t ldda,
    float                *hA,   magma_minproduct_int_t lda,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb );

void
magma_minproduct_ssetmatrix_1D_row_bcyclic(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const float *hA,   magma_minproduct_int_t lda,
    magma_minproductFloat_ptr    dA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb );

// in src/ssytrd_mgpu.cpp
// TODO rename ssetmatrix_sy or similar
magma_minproduct_int_t
magma_minproduct_shtodhe(
    magma_minproduct_int_t ngpu, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    float     *A,   magma_minproduct_int_t lda,
    magma_minproductFloat_ptr dA[], magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queues[][10],
    magma_minproduct_int_t *info );

// in src/spotrf3_mgpu.cpp
// TODO same as magma_minproduct_shtodhe?
magma_minproduct_int_t
magma_minproduct_shtodpo(
    magma_minproduct_int_t ngpu, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t off_i, magma_minproduct_int_t off_j, magma_minproduct_int_t nb,
    float     *A,   magma_minproduct_int_t lda,
    magma_minproductFloat_ptr dA[], magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queues[][3],
    magma_minproduct_int_t *info );

// in src/spotrf3_mgpu.cpp
// TODO rename sgetmatrix_sy or similar
magma_minproduct_int_t
magma_minproduct_sdtohpo(
    magma_minproduct_int_t ngpu, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t off_i, magma_minproduct_int_t off_j, magma_minproduct_int_t nb, magma_minproduct_int_t NB,
    float     *A,   magma_minproduct_int_t lda,
    magma_minproductFloat_ptr dA[], magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queues[][3],
    magma_minproduct_int_t *info );


  /*
   * Multi-GPU BLAS functions (alphabetical order)
   */
void
magma_minproductblas_ssymm_mgpu(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_ptr dA[],    magma_minproduct_int_t ldda,  magma_minproduct_int_t offset,
    magma_minproductFloat_ptr dB[],    magma_minproduct_int_t lddb,
    float beta,
    magma_minproductFloat_ptr dC[],    magma_minproduct_int_t lddc,
    magma_minproductFloat_ptr dwork[], magma_minproduct_int_t dworksiz,
    float    *C,       magma_minproduct_int_t ldc,
    float    *work[],  magma_minproduct_int_t worksiz,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][20], magma_minproduct_int_t nbevents );

void
magma_minproductblas_ssymm_mgpu_com(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_ptr dA[],    magma_minproduct_int_t ldda,  magma_minproduct_int_t offset,
    magma_minproductFloat_ptr dB[],    magma_minproduct_int_t lddb,
    float beta,
    magma_minproductFloat_ptr dC[],    magma_minproduct_int_t lddc,
    magma_minproductFloat_ptr dwork[], magma_minproduct_int_t dworksiz,
    float    *C,       magma_minproduct_int_t ldc,
    float    *work[],  magma_minproduct_int_t worksiz,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents,
    magma_minproduct_int_t gnode[Magma_minproductMaxGPUs][Magma_minproductMaxGPUs+2], magma_minproduct_int_t nbcmplx );

void
magma_minproductblas_ssymm_mgpu_spec(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_ptr dA[],    magma_minproduct_int_t ldda,  magma_minproduct_int_t offset,
    magma_minproductFloat_ptr dB[],    magma_minproduct_int_t lddb,
    float beta,
    magma_minproductFloat_ptr dC[],    magma_minproduct_int_t lddc,
    magma_minproductFloat_ptr dwork[], magma_minproduct_int_t dworksiz,
    float    *C,       magma_minproduct_int_t ldc,
    float    *work[],  magma_minproduct_int_t worksiz,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents,
    magma_minproduct_int_t gnode[Magma_minproductMaxGPUs][Magma_minproductMaxGPUs+2], magma_minproduct_int_t nbcmplx );

void
magma_minproductblas_ssymm_mgpu_spec33(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_ptr dA[],    magma_minproduct_int_t ldda,  magma_minproduct_int_t offset,
    magma_minproductFloat_ptr dB[],    magma_minproduct_int_t lddb,
    float beta,
    magma_minproductFloat_ptr dC[],    magma_minproduct_int_t lddc,
    magma_minproductFloat_ptr dVIN[],  magma_minproduct_int_t lddv, magma_minproduct_int_t v_offset,
    magma_minproductFloat_ptr dwork[], magma_minproduct_int_t dworksiz,
    float *C,       magma_minproduct_int_t ldc,
    float *work[],  magma_minproduct_int_t worksiz,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents,
    magma_minproduct_int_t gnode[Magma_minproductMaxGPUs][Magma_minproductMaxGPUs+2], magma_minproduct_int_t nbcmplx );

magma_minproduct_int_t
magma_minproductblas_ssymv_mgpu(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr const d_lA[], magma_minproduct_int_t ldda, magma_minproduct_int_t offset,
    float const *x,         magma_minproduct_int_t incx,
    float beta,             
    float       *y,         magma_minproduct_int_t incy,
    float       *hwork,     magma_minproduct_int_t lhwork,
    magma_minproductFloat_ptr    dwork[],   magma_minproduct_int_t ldwork,
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[] );

magma_minproduct_int_t
magma_minproductblas_ssymv_mgpu_sync(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr const d_lA[], magma_minproduct_int_t ldda, magma_minproduct_int_t offset,
    float const *x,         magma_minproduct_int_t incx,
    float beta,             
    float       *y,         magma_minproduct_int_t incy,
    float       *hwork,     magma_minproduct_int_t lhwork,
    magma_minproductFloat_ptr    dwork[],   magma_minproduct_int_t ldwork,
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[] );

// Ichi's version, in src/ssytrd_mgpu.cpp
void
magma_minproduct_ssyr2k_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t nb, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloat_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    float beta,
    magma_minproductFloat_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t nqueue, magma_minproduct_queue_t queues[][10] );

void
magma_minproductblas_ssyr2k_mgpu_spec(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloat_ptr dA[], magma_minproduct_int_t ldda, magma_minproduct_int_t a_offset,
    magma_minproductFloat_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    float beta,
    magma_minproductFloat_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb, magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue );

void
magma_minproductblas_ssyr2k_mgpu_spec324(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_ptr dVIN[], magma_minproduct_int_t lddv, magma_minproduct_int_t v_offset,
    magma_minproductFloat_ptr dWIN[], magma_minproduct_int_t lddw, magma_minproduct_int_t w_offset,
    float beta,
    magma_minproductFloat_ptr dC[],   magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproductFloat_ptr dwork[], magma_minproduct_int_t lndwork,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents );

void
magma_minproductblas_ssyr2k_mgpu_spec325(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_ptr dVIN[], magma_minproduct_int_t lddv, magma_minproduct_int_t v_offset,
    magma_minproductFloat_ptr dWIN[], magma_minproduct_int_t lddw, magma_minproduct_int_t w_offset,
    float beta,
    magma_minproductFloat_ptr dC[],   magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproductFloat_ptr dwork[], magma_minproduct_int_t lndwork,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    float *harray[],
    magma_minproductFloat_ptr *darray[],
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_event_t redevents[][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10], magma_minproduct_int_t nbevents );

void
magma_minproductblas_ssyr2k_mgpu2(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloat_ptr dA[], magma_minproduct_int_t ldda, magma_minproduct_int_t a_offset,
    magma_minproductFloat_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    float beta,
    magma_minproductFloat_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue );

// in src/spotrf_mgpu_right.cpp
void
magma_minproduct_ssyrk_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t nb, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloat_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    float beta,
    magma_minproductFloat_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t nqueue, magma_minproduct_queue_t queues[][10]);

// in src/spotrf_mgpu_right.cpp
void
magma_minproduct_ssyrk_mgpu2(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t nb, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloat_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    float beta,
    magma_minproductFloat_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t nqueue, magma_minproduct_queue_t queues[][10]);


  /*
   * LAPACK auxiliary functions (alphabetical order)
   */
void
magma_minproductblas_sgeadd(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr       dB, magma_minproduct_int_t lddb );

void
magma_minproductblas_slacpy(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr       dB, magma_minproduct_int_t lddb );

void
magma_minproductblas_slacpy_cnjg(
    magma_minproduct_int_t n, float *dA1, magma_minproduct_int_t lda1,
    float *dA2, magma_minproduct_int_t lda2);

float
magma_minproductblas_slange(
    magma_minproduct_norm_t norm,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dwork );

float
magma_minproductblas_slansy(
    magma_minproduct_norm_t norm, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dwork );

float
magma_minproductblas_slansy(
    magma_minproduct_norm_t norm, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dwork );

void
magma_minproductblas_slarfg(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dalpha, magma_minproductFloat_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr dtau );

void
magma_minproductblas_slarfg_work(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dalpha, magma_minproductFloat_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr dtau, magma_minproductFloat_ptr dwork );

void
magma_minproductblas_slascl(
    magma_minproduct_type_t type, magma_minproduct_int_t kl, magma_minproduct_int_t ku,
    float cfrom, float cto,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_slascl_2x2(
    magma_minproduct_type_t type, magma_minproduct_int_t m,
    float *dW, magma_minproduct_int_t lddw,
    float *dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_slascl2(
    magma_minproduct_type_t type,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dD,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_slascl_diag(
    magma_minproduct_type_t type, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dD, magma_minproduct_int_t lddd,
          magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info );

void
magma_minproductblas_slaset(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    float offdiag, float diag,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda );

void
magma_minproductblas_slaset_band(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float offdiag, float diag,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda);

void
magma_minproductblas_slaswp(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dAT, magma_minproduct_int_t ldda,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci );

void
magma_minproductblas_slaswp2(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dAT, magma_minproduct_int_t ldda,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    magma_minproductInt_const_ptr d_ipiv, magma_minproduct_int_t inci );

void
magma_minproductblas_slaswpx(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldx, magma_minproduct_int_t ldy,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci );

void
magma_minproductblas_ssymmetrize(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda );

void
magma_minproductblas_ssymmetrize_tiles(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t ntile, magma_minproduct_int_t mstride, magma_minproduct_int_t nstride );

void
magma_minproductblas_strtri_diag(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr d_dinvA);

  /*
   * to cleanup (alphabetical order)
   */
void
magma_minproductblas_snrm2_adjust(
    magma_minproduct_int_t k,
    magma_minproductFloat_ptr dxnorm,
    magma_minproductFloat_ptr dc);

void
magma_minproductblas_snrm2_check(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dxnorm,
    magma_minproductFloat_ptr dlsticc);

void
magma_minproductblas_snrm2_cols(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dxnorm);

void
magma_minproductblas_snrm2_row_check_adjust(
    magma_minproduct_int_t k, float tol,
    magma_minproductFloat_ptr dxnorm,
    magma_minproductFloat_ptr dxnorm2,
    magma_minproductFloat_ptr dC, magma_minproduct_int_t lddc,
    magma_minproductFloat_ptr dlsticc);

void
magma_minproduct_slarfbx_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t k,
    magma_minproductFloat_ptr V,  magma_minproduct_int_t ldv,
    magma_minproductFloat_ptr dT, magma_minproduct_int_t ldt,
    magma_minproductFloat_ptr c,
    magma_minproductFloat_ptr dwork);

void
magma_minproduct_slarfg_gpu(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dx0,
    magma_minproductFloat_ptr dx,
    magma_minproductFloat_ptr dtau,
    magma_minproductFloat_ptr        dxnorm,
    magma_minproductFloat_ptr dAkk );

void
magma_minproduct_slarfgtx_gpu(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dx0,
    magma_minproductFloat_ptr dx,
    magma_minproductFloat_ptr dtau,
    magma_minproductFloat_ptr        dxnorm,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t iter,
    magma_minproductFloat_ptr V,  magma_minproduct_int_t ldv,
    magma_minproductFloat_ptr T,  magma_minproduct_int_t ldt,
    magma_minproductFloat_ptr dwork);

void
magma_minproduct_slarfgx_gpu(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dx0,
    magma_minproductFloat_ptr dx,
    magma_minproductFloat_ptr dtau,
    magma_minproductFloat_ptr        dxnorm,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t iter);

void
magma_minproduct_slarfx_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr v,
    magma_minproductFloat_ptr tau,
    magma_minproductFloat_ptr C,  magma_minproduct_int_t ldc,
    magma_minproductFloat_ptr        xnorm,
    magma_minproductFloat_ptr dT, magma_minproduct_int_t iter,
    magma_minproductFloat_ptr work);


  /*
   * Level 1 BLAS (alphabetical order)
   */
void
magma_minproductblas_sswap(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr dy, magma_minproduct_int_t incy );

void
magma_minproductblas_sswapblk(
    magma_minproduct_order_t order,
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t i1, magma_minproduct_int_t i2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci,
    magma_minproduct_int_t offset );

void
magma_minproductblas_sswapdblk(
    magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t inca,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb, magma_minproduct_int_t incb );

  /*
   * Level 2 BLAS (alphabetical order)
   */
void
magma_minproductblas_sgemv(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx,
    float beta,
    magma_minproductFloat_ptr       dy, magma_minproduct_int_t incy );


void
magma_minproductblas_sgemv_conjv(
    magma_minproduct_int_t m, magma_minproduct_int_t n, float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx,
    float beta,
    magma_minproductFloat_ptr dy, magma_minproduct_int_t incy);

magma_minproduct_int_t
magma_minproductblas_ssymv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx,
    float beta,
    magma_minproductFloat_ptr       dy, magma_minproduct_int_t incy );

magma_minproduct_int_t
magma_minproductblas_ssymv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx,
    float beta,
    magma_minproductFloat_ptr       dy, magma_minproduct_int_t incy );

magma_minproduct_int_t
magma_minproductblas_ssymv_work(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx,
    float beta,
    magma_minproductFloat_ptr       dy, magma_minproduct_int_t incy,
    magma_minproductFloat_ptr       dwork, magma_minproduct_int_t lwork,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproductblas_ssymv_work(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx,
    float beta,
    magma_minproductFloat_ptr       dy, magma_minproduct_int_t incy,
    magma_minproductFloat_ptr       dwork, magma_minproduct_int_t lwork,
    magma_minproduct_queue_t queue );

  /*
   * Level 3 BLAS (alphabetical order)
   */
void
magma_minproductblas_sgemm(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dB, magma_minproduct_int_t lddb,
    float beta,
    magma_minproductFloat_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_sgemm_reduce(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dB, magma_minproduct_int_t lddb,
    float beta,
    magma_minproductFloat_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_ssymm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dB, magma_minproduct_int_t lddb,
    float beta,
    magma_minproductFloat_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_ssymm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dB, magma_minproduct_int_t lddb,
    float beta,
    magma_minproductFloat_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_ssyr2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dB, magma_minproduct_int_t lddb,
    float beta,
    magma_minproductFloat_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_ssyr2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dB, magma_minproduct_int_t lddb,
    float  beta,
    magma_minproductFloat_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_ssyrk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    float beta,
    magma_minproductFloat_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_ssyrk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    float  alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    float  beta,
    magma_minproductFloat_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproductblas_strsm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr       dB, magma_minproduct_int_t lddb );

void
magma_minproductblas_strsm_outofplace(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr       dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t flag, magma_minproductFloat_ptr d_dinvA, magma_minproductFloat_ptr dX );

void
magma_minproductblas_strsm_work(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr       dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t flag, magma_minproductFloat_ptr d_dinvA, magma_minproductFloat_ptr dX );


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

#define magma_minproduct_ssetvector(           n, hx_src, incx, dy_dst, incy ) \
        magma_minproduct_ssetvector_internal(  n, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_minproduct_sgetvector(           n, dx_src, incx, hy_dst, incy ) \
        magma_minproduct_sgetvector_internal(  n, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_minproduct_scopyvector(          n, dx_src, incx, dy_dst, incy ) \
        magma_minproduct_scopyvector_internal( n, dx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_minproduct_ssetvector_async(           n, hx_src, incx, dy_dst, incy, queue ) \
        magma_minproduct_ssetvector_async_internal(  n, hx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_sgetvector_async(           n, dx_src, incx, hy_dst, incy, queue ) \
        magma_minproduct_sgetvector_async_internal(  n, dx_src, incx, hy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_scopyvector_async(          n, dx_src, incx, dy_dst, incy, queue ) \
        magma_minproduct_scopyvector_async_internal( n, dx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

void
magma_minproduct_ssetvector_internal(
    magma_minproduct_int_t n,
    float const    *hx_src, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr       dy_dst, magma_minproduct_int_t incy,
    const char* func, const char* file, int line );

void
magma_minproduct_sgetvector_internal(
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dx_src, magma_minproduct_int_t incx,
    float          *hy_dst, magma_minproduct_int_t incy,
    const char* func, const char* file, int line );

void
magma_minproduct_scopyvector_internal(
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dx_src, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr       dy_dst, magma_minproduct_int_t incy,
    const char* func, const char* file, int line );

void
magma_minproduct_ssetvector_async_internal(
    magma_minproduct_int_t n,
    float const    *hx_src, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr       dy_dst, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void
magma_minproduct_sgetvector_async_internal(
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dx_src, magma_minproduct_int_t incx,
    float          *hy_dst, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void
magma_minproduct_scopyvector_async_internal(
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dx_src, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr       dy_dst, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// copying sub-matrices (contiguous columns)
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

#define magma_minproduct_ssetmatrix(           m, n, hA_src, lda, dB_dst, lddb ) \
        magma_minproduct_ssetmatrix_internal(  m, n, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_minproduct_sgetmatrix(           m, n, dA_src, ldda, hB_dst, ldb ) \
        magma_minproduct_sgetmatrix_internal(  m, n, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )

#define magma_minproduct_scopymatrix(          m, n, dA_src, ldda, dB_dst, lddb ) \
        magma_minproduct_scopymatrix_internal( m, n, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_minproduct_ssetmatrix_async(           m, n, hA_src, lda, dB_dst, lddb, queue ) \
        magma_minproduct_ssetmatrix_async_internal(  m, n, hA_src, lda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_sgetmatrix_async(           m, n, dA_src, ldda, hB_dst, ldb, queue ) \
        magma_minproduct_sgetmatrix_async_internal(  m, n, dA_src, ldda, hB_dst, ldb, queue, __func__, __FILE__, __LINE__ )

#define magma_minproduct_scopymatrix_async(          m, n, dA_src, ldda, dB_dst, lddb, queue ) \
        magma_minproduct_scopymatrix_async_internal( m, n, dA_src, ldda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

void
magma_minproduct_ssetmatrix_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float const    *hA_src, magma_minproduct_int_t ldha,
    magma_minproductFloat_ptr       dB_dst, magma_minproduct_int_t lddb,
    const char* func, const char* file, int line );

void
magma_minproduct_sgetmatrix_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA_src, magma_minproduct_int_t ldda,
    float          *hB_dst, magma_minproduct_int_t ldhb,
    const char* func, const char* file, int line );

void
magma_minproduct_scopymatrix_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA_src, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr       dB_dst, magma_minproduct_int_t lddb,
    const char* func, const char* file, int line );

void
magma_minproduct_ssetmatrix_async_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float const    *hA_src, magma_minproduct_int_t ldha,
    magma_minproductFloat_ptr       dB_dst, magma_minproduct_int_t lddb,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void
magma_minproduct_sgetmatrix_async_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA_src, magma_minproduct_int_t ldda,
    float          *hB_dst, magma_minproduct_int_t ldhb,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );

void
magma_minproduct_scopymatrix_async_internal(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA_src, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr       dB_dst, magma_minproduct_int_t lddb,
    magma_minproduct_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// Level 1 BLAS (alphabetical order)

// in cublas_v2, result returned through output argument
magma_minproduct_int_t
magma_minproduct_isamax(
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx );

// in cublas_v2, result returned through output argument
magma_minproduct_int_t
magma_minproduct_isamin(
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx );

// in cublas_v2, result returned through output argument
float
magma_minproduct_sasum(
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx );

void
magma_minproduct_saxpy(
    magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr       dy, magma_minproduct_int_t incy );

void
magma_minproduct_scopy(
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr       dy, magma_minproduct_int_t incy );

// in cublas_v2, result returned through output argument
float
magma_minproduct_sdot(
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_const_ptr dy, magma_minproduct_int_t incy );

// in cublas_v2, result returned through output argument
float
magma_minproduct_sdot(
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_const_ptr dy, magma_minproduct_int_t incy );

// in cublas_v2, result returned through output argument
float
magma_minproduct_snrm2(
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx );

void
magma_minproduct_srot(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr dy, magma_minproduct_int_t incy,
    float dc, float ds );

void
magma_minproduct_srot(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr dy, magma_minproduct_int_t incy,
    float dc, float ds );

#ifdef REAL
void
magma_minproduct_srotm(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr dy, magma_minproduct_int_t incy,
    magma_minproductFloat_const_ptr param );

void
magma_minproduct_srotmg(
    magma_minproductFloat_ptr d1, magma_minproductFloat_ptr       d2,
    magma_minproductFloat_ptr x1, magma_minproductFloat_const_ptr y1,
    magma_minproductFloat_ptr param );
#endif

void
magma_minproduct_sscal(
    magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_ptr dx, magma_minproduct_int_t incx );

void
magma_minproduct_sscal(
    magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_ptr dx, magma_minproduct_int_t incx );

void
magma_minproduct_sswap(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr dy, magma_minproduct_int_t incy );

// ========================================
// Level 2 BLAS (alphabetical order)

void
magma_minproduct_sgemv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx,
    float beta,
    magma_minproductFloat_ptr       dy, magma_minproduct_int_t incy );

void
magma_minproduct_sger(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_const_ptr dy, magma_minproduct_int_t incy,
    magma_minproductFloat_ptr       dA, magma_minproduct_int_t ldda );

void
magma_minproduct_sger(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_const_ptr dy, magma_minproduct_int_t incy,
    magma_minproductFloat_ptr       dA, magma_minproduct_int_t ldda );

void
magma_minproduct_ssymv(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx,
    float beta,
    magma_minproductFloat_ptr       dy, magma_minproduct_int_t incy );

void
magma_minproduct_ssyr(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr       dA, magma_minproduct_int_t ldda );

void
magma_minproduct_ssyr2(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_const_ptr dy, magma_minproduct_int_t incy,
    magma_minproductFloat_ptr       dA, magma_minproduct_int_t ldda );

void
magma_minproduct_strmv(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr       dx, magma_minproduct_int_t incx );

void
magma_minproduct_strsv(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr       dx, magma_minproduct_int_t incx );

// ========================================
// Level 3 BLAS (alphabetical order)

void
magma_minproduct_sgemm(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dB, magma_minproduct_int_t lddb,
    float beta,
    magma_minproductFloat_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_ssymm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dB, magma_minproduct_int_t lddb,
    float beta,
    magma_minproductFloat_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_ssymm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dB, magma_minproduct_int_t lddb,
    float beta,
    magma_minproductFloat_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_ssyr2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dB, magma_minproduct_int_t lddb,
    float beta,
    magma_minproductFloat_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_ssyr2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_const_ptr dB, magma_minproduct_int_t lddb,
    float beta,
    magma_minproductFloat_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_ssyrk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    float beta,
    magma_minproductFloat_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_ssyrk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    float beta,
    magma_minproductFloat_ptr       dC, magma_minproduct_int_t lddc );

void
magma_minproduct_strmm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr       dB, magma_minproduct_int_t lddb );

void
magma_minproduct_strsm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr       dB, magma_minproduct_int_t lddb );


#ifdef __cplusplus
}
#endif

#undef REAL

#endif  /* MAGMA_minproductBLAS_S_H */
