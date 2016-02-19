/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally3blas_z.h normal z -> s, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally3BLAS_S_H
#define MAGMA_tally3BLAS_S_H

#include "magma_tally3_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_tally3blas_stranspose_inplace(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda );

void
magma_tally3blas_stranspose(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3Float_ptr       dAT, magma_tally3_int_t lddat );

void
magma_tally3blas_sgetmatrix_transpose(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dAT,   magma_tally3_int_t ldda,
    float          *hA,    magma_tally3_int_t lda,
    magma_tally3Float_ptr       dwork, magma_tally3_int_t lddwork, magma_tally3_int_t nb );

void
magma_tally3blas_ssetmatrix_transpose(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const float *hA,    magma_tally3_int_t lda,
    magma_tally3Float_ptr    dAT,   magma_tally3_int_t ldda,
    magma_tally3Float_ptr    dwork, magma_tally3_int_t lddwork, magma_tally3_int_t nb );

  /*
   * RBT-related functions
   */
void
magma_tally3blas_sprbt(
    magma_tally3_int_t n, 
    float *dA, magma_tally3_int_t ldda, 
    float *du, float *dv);

void
magma_tally3blas_sprbt_mv(
    magma_tally3_int_t n, 
    float *dv, float *db);

void
magma_tally3blas_sprbt_mtv(
    magma_tally3_int_t n, 
    float *du, float *db);

void
magma_tally3blas_saxpycp2(
    magma_tally3_int_t m, float *r, float *x,
    const float *b);

  /*
   * Multi-GPU copy functions
   */
void
magma_tally3blas_sgetmatrix_transpose_mgpu(
    magma_tally3_int_t ngpu, magma_tally3_queue_t queues[][2],
    magma_tally3Float_const_ptr const dAT[],    magma_tally3_int_t ldda,
    float                *hA,       magma_tally3_int_t lda,
    magma_tally3Float_ptr             dwork[],  magma_tally3_int_t lddw,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb );

void
magma_tally3blas_ssetmatrix_transpose_mgpu(
    magma_tally3_int_t ngpu, magma_tally3_queue_t queues[][2],
    const float *hA,      magma_tally3_int_t lda,
    magma_tally3Float_ptr    dAT[],   magma_tally3_int_t ldda,
    magma_tally3Float_ptr    dwork[], magma_tally3_int_t lddw,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb );

void
magma_tally3_sgetmatrix_1D_col_bcyclic(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr const dA[], magma_tally3_int_t ldda,
    float                *hA,   magma_tally3_int_t lda,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb );

void
magma_tally3_ssetmatrix_1D_col_bcyclic(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const float *hA,   magma_tally3_int_t lda,
    magma_tally3Float_ptr    dA[], magma_tally3_int_t ldda,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb );

void
magma_tally3_sgetmatrix_1D_row_bcyclic(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr const dA[], magma_tally3_int_t ldda,
    float                *hA,   magma_tally3_int_t lda,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb );

void
magma_tally3_ssetmatrix_1D_row_bcyclic(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const float *hA,   magma_tally3_int_t lda,
    magma_tally3Float_ptr    dA[], magma_tally3_int_t ldda,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb );

// in src/ssytrd_mgpu.cpp
// TODO rename ssetmatrix_sy or similar
magma_tally3_int_t
magma_tally3_shtodhe(
    magma_tally3_int_t ngpu, magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    float     *A,   magma_tally3_int_t lda,
    magma_tally3Float_ptr dA[], magma_tally3_int_t ldda,
    magma_tally3_queue_t queues[][10],
    magma_tally3_int_t *info );

// in src/spotrf3_mgpu.cpp
// TODO same as magma_tally3_shtodhe?
magma_tally3_int_t
magma_tally3_shtodpo(
    magma_tally3_int_t ngpu, magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t off_i, magma_tally3_int_t off_j, magma_tally3_int_t nb,
    float     *A,   magma_tally3_int_t lda,
    magma_tally3Float_ptr dA[], magma_tally3_int_t ldda,
    magma_tally3_queue_t queues[][3],
    magma_tally3_int_t *info );

// in src/spotrf3_mgpu.cpp
// TODO rename sgetmatrix_sy or similar
magma_tally3_int_t
magma_tally3_sdtohpo(
    magma_tally3_int_t ngpu, magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t off_i, magma_tally3_int_t off_j, magma_tally3_int_t nb, magma_tally3_int_t NB,
    float     *A,   magma_tally3_int_t lda,
    magma_tally3Float_ptr dA[], magma_tally3_int_t ldda,
    magma_tally3_queue_t queues[][3],
    magma_tally3_int_t *info );


  /*
   * Multi-GPU BLAS functions (alphabetical order)
   */
void
magma_tally3blas_ssymm_mgpu(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_ptr dA[],    magma_tally3_int_t ldda,  magma_tally3_int_t offset,
    magma_tally3Float_ptr dB[],    magma_tally3_int_t lddb,
    float beta,
    magma_tally3Float_ptr dC[],    magma_tally3_int_t lddc,
    magma_tally3Float_ptr dwork[], magma_tally3_int_t dworksiz,
    float    *C,       magma_tally3_int_t ldc,
    float    *work[],  magma_tally3_int_t worksiz,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_event_t redevents[][20], magma_tally3_int_t nbevents );

void
magma_tally3blas_ssymm_mgpu_com(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_ptr dA[],    magma_tally3_int_t ldda,  magma_tally3_int_t offset,
    magma_tally3Float_ptr dB[],    magma_tally3_int_t lddb,
    float beta,
    magma_tally3Float_ptr dC[],    magma_tally3_int_t lddc,
    magma_tally3Float_ptr dwork[], magma_tally3_int_t dworksiz,
    float    *C,       magma_tally3_int_t ldc,
    float    *work[],  magma_tally3_int_t worksiz,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_event_t redevents[][Magma_tally3MaxGPUs*Magma_tally3MaxGPUs+10], magma_tally3_int_t nbevents,
    magma_tally3_int_t gnode[Magma_tally3MaxGPUs][Magma_tally3MaxGPUs+2], magma_tally3_int_t nbcmplx );

void
magma_tally3blas_ssymm_mgpu_spec(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_ptr dA[],    magma_tally3_int_t ldda,  magma_tally3_int_t offset,
    magma_tally3Float_ptr dB[],    magma_tally3_int_t lddb,
    float beta,
    magma_tally3Float_ptr dC[],    magma_tally3_int_t lddc,
    magma_tally3Float_ptr dwork[], magma_tally3_int_t dworksiz,
    float    *C,       magma_tally3_int_t ldc,
    float    *work[],  magma_tally3_int_t worksiz,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_event_t redevents[][Magma_tally3MaxGPUs*Magma_tally3MaxGPUs+10], magma_tally3_int_t nbevents,
    magma_tally3_int_t gnode[Magma_tally3MaxGPUs][Magma_tally3MaxGPUs+2], magma_tally3_int_t nbcmplx );

void
magma_tally3blas_ssymm_mgpu_spec33(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_ptr dA[],    magma_tally3_int_t ldda,  magma_tally3_int_t offset,
    magma_tally3Float_ptr dB[],    magma_tally3_int_t lddb,
    float beta,
    magma_tally3Float_ptr dC[],    magma_tally3_int_t lddc,
    magma_tally3Float_ptr dVIN[],  magma_tally3_int_t lddv, magma_tally3_int_t v_offset,
    magma_tally3Float_ptr dwork[], magma_tally3_int_t dworksiz,
    float *C,       magma_tally3_int_t ldc,
    float *work[],  magma_tally3_int_t worksiz,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_event_t redevents[][Magma_tally3MaxGPUs*Magma_tally3MaxGPUs+10], magma_tally3_int_t nbevents,
    magma_tally3_int_t gnode[Magma_tally3MaxGPUs][Magma_tally3MaxGPUs+2], magma_tally3_int_t nbcmplx );

magma_tally3_int_t
magma_tally3blas_ssymv_mgpu(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr const d_lA[], magma_tally3_int_t ldda, magma_tally3_int_t offset,
    float const *x,         magma_tally3_int_t incx,
    float beta,             
    float       *y,         magma_tally3_int_t incy,
    float       *hwork,     magma_tally3_int_t lhwork,
    magma_tally3Float_ptr    dwork[],   magma_tally3_int_t ldwork,
    magma_tally3_int_t ngpu,
    magma_tally3_int_t nb,
    magma_tally3_queue_t queues[] );

magma_tally3_int_t
magma_tally3blas_ssymv_mgpu_sync(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr const d_lA[], magma_tally3_int_t ldda, magma_tally3_int_t offset,
    float const *x,         magma_tally3_int_t incx,
    float beta,             
    float       *y,         magma_tally3_int_t incy,
    float       *hwork,     magma_tally3_int_t lhwork,
    magma_tally3Float_ptr    dwork[],   magma_tally3_int_t ldwork,
    magma_tally3_int_t ngpu,
    magma_tally3_int_t nb,
    magma_tally3_queue_t queues[] );

// Ichi's version, in src/ssytrd_mgpu.cpp
void
magma_tally3_ssyr2k_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t nb, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_ptr dB[], magma_tally3_int_t lddb, magma_tally3_int_t b_offset,
    float beta,
    magma_tally3Float_ptr dC[], magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3_int_t nqueue, magma_tally3_queue_t queues[][10] );

void
magma_tally3blas_ssyr2k_mgpu_spec(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_ptr dA[], magma_tally3_int_t ldda, magma_tally3_int_t a_offset,
    magma_tally3Float_ptr dB[], magma_tally3_int_t lddb, magma_tally3_int_t b_offset,
    float beta,
    magma_tally3Float_ptr dC[], magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb, magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue );

void
magma_tally3blas_ssyr2k_mgpu_spec324(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_ptr dVIN[], magma_tally3_int_t lddv, magma_tally3_int_t v_offset,
    magma_tally3Float_ptr dWIN[], magma_tally3_int_t lddw, magma_tally3_int_t w_offset,
    float beta,
    magma_tally3Float_ptr dC[],   magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3Float_ptr dwork[], magma_tally3_int_t lndwork,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_event_t redevents[][Magma_tally3MaxGPUs*Magma_tally3MaxGPUs+10], magma_tally3_int_t nbevents );

void
magma_tally3blas_ssyr2k_mgpu_spec325(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_ptr dVIN[], magma_tally3_int_t lddv, magma_tally3_int_t v_offset,
    magma_tally3Float_ptr dWIN[], magma_tally3_int_t lddw, magma_tally3_int_t w_offset,
    float beta,
    magma_tally3Float_ptr dC[],   magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3Float_ptr dwork[], magma_tally3_int_t lndwork,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb,
    float *harray[],
    magma_tally3Float_ptr *darray[],
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_event_t redevents[][Magma_tally3MaxGPUs*Magma_tally3MaxGPUs+10], magma_tally3_int_t nbevents );

void
magma_tally3blas_ssyr2k_mgpu2(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_ptr dA[], magma_tally3_int_t ldda, magma_tally3_int_t a_offset,
    magma_tally3Float_ptr dB[], magma_tally3_int_t lddb, magma_tally3_int_t b_offset,
    float beta,
    magma_tally3Float_ptr dC[], magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue );

// in src/spotrf_mgpu_right.cpp
void
magma_tally3_ssyrk_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t nb, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_ptr dB[], magma_tally3_int_t lddb, magma_tally3_int_t b_offset,
    float beta,
    magma_tally3Float_ptr dC[], magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3_int_t nqueue, magma_tally3_queue_t queues[][10]);

// in src/spotrf_mgpu_right.cpp
void
magma_tally3_ssyrk_mgpu2(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t nb, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_ptr dB[], magma_tally3_int_t lddb, magma_tally3_int_t b_offset,
    float beta,
    magma_tally3Float_ptr dC[], magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3_int_t nqueue, magma_tally3_queue_t queues[][10]);


  /*
   * LAPACK auxiliary functions (alphabetical order)
   */
void
magma_tally3blas_sgeadd(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr       dB, magma_tally3_int_t lddb );

void
magma_tally3blas_slacpy(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr       dB, magma_tally3_int_t lddb );

void
magma_tally3blas_slacpy_cnjg(
    magma_tally3_int_t n, float *dA1, magma_tally3_int_t lda1,
    float *dA2, magma_tally3_int_t lda2);

float
magma_tally3blas_slange(
    magma_tally3_norm_t norm,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dwork );

float
magma_tally3blas_slansy(
    magma_tally3_norm_t norm, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dwork );

float
magma_tally3blas_slansy(
    magma_tally3_norm_t norm, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dwork );

void
magma_tally3blas_slarfg(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dalpha, magma_tally3Float_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_ptr dtau );

void
magma_tally3blas_slarfg_work(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dalpha, magma_tally3Float_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_ptr dtau, magma_tally3Float_ptr dwork );

void
magma_tally3blas_slascl(
    magma_tally3_type_t type, magma_tally3_int_t kl, magma_tally3_int_t ku,
    float cfrom, float cto,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info );

void
magma_tally3blas_slascl_2x2(
    magma_tally3_type_t type, magma_tally3_int_t m,
    float *dW, magma_tally3_int_t lddw,
    float *dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info );

void
magma_tally3blas_slascl2(
    magma_tally3_type_t type,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dD,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info );

void
magma_tally3blas_slascl_diag(
    magma_tally3_type_t type, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dD, magma_tally3_int_t lddd,
          magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info );

void
magma_tally3blas_slaset(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    float offdiag, float diag,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda );

void
magma_tally3blas_slaset_band(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    float offdiag, float diag,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda);

void
magma_tally3blas_slaswp(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dAT, magma_tally3_int_t ldda,
    magma_tally3_int_t k1, magma_tally3_int_t k2,
    const magma_tally3_int_t *ipiv, magma_tally3_int_t inci );

void
magma_tally3blas_slaswp2(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dAT, magma_tally3_int_t ldda,
    magma_tally3_int_t k1, magma_tally3_int_t k2,
    magma_tally3Int_const_ptr d_ipiv, magma_tally3_int_t inci );

void
magma_tally3blas_slaswpx(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldx, magma_tally3_int_t ldy,
    magma_tally3_int_t k1, magma_tally3_int_t k2,
    const magma_tally3_int_t *ipiv, magma_tally3_int_t inci );

void
magma_tally3blas_ssymmetrize(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda );

void
magma_tally3blas_ssymmetrize_tiles(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t ntile, magma_tally3_int_t mstride, magma_tally3_int_t nstride );

void
magma_tally3blas_strtri_diag(
    magma_tally3_uplo_t uplo, magma_tally3_diag_t diag, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr d_dinvA);

  /*
   * to cleanup (alphabetical order)
   */
void
magma_tally3blas_snrm2_adjust(
    magma_tally3_int_t k,
    magma_tally3Float_ptr dxnorm,
    magma_tally3Float_ptr dc);

void
magma_tally3blas_snrm2_check(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dxnorm,
    magma_tally3Float_ptr dlsticc);

void
magma_tally3blas_snrm2_cols(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dxnorm);

void
magma_tally3blas_snrm2_row_check_adjust(
    magma_tally3_int_t k, float tol,
    magma_tally3Float_ptr dxnorm,
    magma_tally3Float_ptr dxnorm2,
    magma_tally3Float_ptr dC, magma_tally3_int_t lddc,
    magma_tally3Float_ptr dlsticc);

void
magma_tally3_slarfbx_gpu(
    magma_tally3_int_t m, magma_tally3_int_t k,
    magma_tally3Float_ptr V,  magma_tally3_int_t ldv,
    magma_tally3Float_ptr dT, magma_tally3_int_t ldt,
    magma_tally3Float_ptr c,
    magma_tally3Float_ptr dwork);

void
magma_tally3_slarfg_gpu(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dx0,
    magma_tally3Float_ptr dx,
    magma_tally3Float_ptr dtau,
    magma_tally3Float_ptr        dxnorm,
    magma_tally3Float_ptr dAkk );

void
magma_tally3_slarfgtx_gpu(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dx0,
    magma_tally3Float_ptr dx,
    magma_tally3Float_ptr dtau,
    magma_tally3Float_ptr        dxnorm,
    magma_tally3Float_ptr dA, magma_tally3_int_t iter,
    magma_tally3Float_ptr V,  magma_tally3_int_t ldv,
    magma_tally3Float_ptr T,  magma_tally3_int_t ldt,
    magma_tally3Float_ptr dwork);

void
magma_tally3_slarfgx_gpu(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dx0,
    magma_tally3Float_ptr dx,
    magma_tally3Float_ptr dtau,
    magma_tally3Float_ptr        dxnorm,
    magma_tally3Float_ptr dA, magma_tally3_int_t iter);

void
magma_tally3_slarfx_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr v,
    magma_tally3Float_ptr tau,
    magma_tally3Float_ptr C,  magma_tally3_int_t ldc,
    magma_tally3Float_ptr        xnorm,
    magma_tally3Float_ptr dT, magma_tally3_int_t iter,
    magma_tally3Float_ptr work);


  /*
   * Level 1 BLAS (alphabetical order)
   */
void
magma_tally3blas_sswap(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_ptr dy, magma_tally3_int_t incy );

void
magma_tally3blas_sswapblk(
    magma_tally3_order_t order,
    magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t i1, magma_tally3_int_t i2,
    const magma_tally3_int_t *ipiv, magma_tally3_int_t inci,
    magma_tally3_int_t offset );

void
magma_tally3blas_sswapdblk(
    magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t inca,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb, magma_tally3_int_t incb );

  /*
   * Level 2 BLAS (alphabetical order)
   */
void
magma_tally3blas_sgemv(
    magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx,
    float beta,
    magma_tally3Float_ptr       dy, magma_tally3_int_t incy );


void
magma_tally3blas_sgemv_conjv(
    magma_tally3_int_t m, magma_tally3_int_t n, float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx,
    float beta,
    magma_tally3Float_ptr dy, magma_tally3_int_t incy);

magma_tally3_int_t
magma_tally3blas_ssymv(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx,
    float beta,
    magma_tally3Float_ptr       dy, magma_tally3_int_t incy );

magma_tally3_int_t
magma_tally3blas_ssymv(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx,
    float beta,
    magma_tally3Float_ptr       dy, magma_tally3_int_t incy );

magma_tally3_int_t
magma_tally3blas_ssymv_work(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx,
    float beta,
    magma_tally3Float_ptr       dy, magma_tally3_int_t incy,
    magma_tally3Float_ptr       dwork, magma_tally3_int_t lwork,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3blas_ssymv_work(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx,
    float beta,
    magma_tally3Float_ptr       dy, magma_tally3_int_t incy,
    magma_tally3Float_ptr       dwork, magma_tally3_int_t lwork,
    magma_tally3_queue_t queue );

  /*
   * Level 3 BLAS (alphabetical order)
   */
void
magma_tally3blas_sgemm(
    magma_tally3_trans_t transA, magma_tally3_trans_t transB,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dB, magma_tally3_int_t lddb,
    float beta,
    magma_tally3Float_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3blas_sgemm_reduce(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dB, magma_tally3_int_t lddb,
    float beta,
    magma_tally3Float_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3blas_ssymm(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dB, magma_tally3_int_t lddb,
    float beta,
    magma_tally3Float_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3blas_ssymm(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dB, magma_tally3_int_t lddb,
    float beta,
    magma_tally3Float_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3blas_ssyr2k(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dB, magma_tally3_int_t lddb,
    float beta,
    magma_tally3Float_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3blas_ssyr2k(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dB, magma_tally3_int_t lddb,
    float  beta,
    magma_tally3Float_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3blas_ssyrk(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    float beta,
    magma_tally3Float_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3blas_ssyrk(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t n, magma_tally3_int_t k,
    float  alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    float  beta,
    magma_tally3Float_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3blas_strsm(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transA, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr       dB, magma_tally3_int_t lddb );

void
magma_tally3blas_strsm_outofplace(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transA, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr       dB, magma_tally3_int_t lddb,
    magma_tally3_int_t flag, magma_tally3Float_ptr d_dinvA, magma_tally3Float_ptr dX );

void
magma_tally3blas_strsm_work(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transA, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr       dB, magma_tally3_int_t lddb,
    magma_tally3_int_t flag, magma_tally3Float_ptr d_dinvA, magma_tally3Float_ptr dX );


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

#define magma_tally3_ssetvector(           n, hx_src, incx, dy_dst, incy ) \
        magma_tally3_ssetvector_internal(  n, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally3_sgetvector(           n, dx_src, incx, hy_dst, incy ) \
        magma_tally3_sgetvector_internal(  n, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally3_scopyvector(          n, dx_src, incx, dy_dst, incy ) \
        magma_tally3_scopyvector_internal( n, dx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally3_ssetvector_async(           n, hx_src, incx, dy_dst, incy, queue ) \
        magma_tally3_ssetvector_async_internal(  n, hx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_tally3_sgetvector_async(           n, dx_src, incx, hy_dst, incy, queue ) \
        magma_tally3_sgetvector_async_internal(  n, dx_src, incx, hy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_tally3_scopyvector_async(          n, dx_src, incx, dy_dst, incy, queue ) \
        magma_tally3_scopyvector_async_internal( n, dx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

void
magma_tally3_ssetvector_internal(
    magma_tally3_int_t n,
    float const    *hx_src, magma_tally3_int_t incx,
    magma_tally3Float_ptr       dy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally3_sgetvector_internal(
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dx_src, magma_tally3_int_t incx,
    float          *hy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally3_scopyvector_internal(
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dx_src, magma_tally3_int_t incx,
    magma_tally3Float_ptr       dy_dst, magma_tally3_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally3_ssetvector_async_internal(
    magma_tally3_int_t n,
    float const    *hx_src, magma_tally3_int_t incx,
    magma_tally3Float_ptr       dy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally3_sgetvector_async_internal(
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dx_src, magma_tally3_int_t incx,
    float          *hy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally3_scopyvector_async_internal(
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dx_src, magma_tally3_int_t incx,
    magma_tally3Float_ptr       dy_dst, magma_tally3_int_t incy,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// copying sub-matrices (contiguous columns)
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

#define magma_tally3_ssetmatrix(           m, n, hA_src, lda, dB_dst, lddb ) \
        magma_tally3_ssetmatrix_internal(  m, n, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_tally3_sgetmatrix(           m, n, dA_src, ldda, hB_dst, ldb ) \
        magma_tally3_sgetmatrix_internal(  m, n, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )

#define magma_tally3_scopymatrix(          m, n, dA_src, ldda, dB_dst, lddb ) \
        magma_tally3_scopymatrix_internal( m, n, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_tally3_ssetmatrix_async(           m, n, hA_src, lda, dB_dst, lddb, queue ) \
        magma_tally3_ssetmatrix_async_internal(  m, n, hA_src, lda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

#define magma_tally3_sgetmatrix_async(           m, n, dA_src, ldda, hB_dst, ldb, queue ) \
        magma_tally3_sgetmatrix_async_internal(  m, n, dA_src, ldda, hB_dst, ldb, queue, __func__, __FILE__, __LINE__ )

#define magma_tally3_scopymatrix_async(          m, n, dA_src, ldda, dB_dst, lddb, queue ) \
        magma_tally3_scopymatrix_async_internal( m, n, dA_src, ldda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

void
magma_tally3_ssetmatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float const    *hA_src, magma_tally3_int_t ldha,
    magma_tally3Float_ptr       dB_dst, magma_tally3_int_t lddb,
    const char* func, const char* file, int line );

void
magma_tally3_sgetmatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA_src, magma_tally3_int_t ldda,
    float          *hB_dst, magma_tally3_int_t ldhb,
    const char* func, const char* file, int line );

void
magma_tally3_scopymatrix_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA_src, magma_tally3_int_t ldda,
    magma_tally3Float_ptr       dB_dst, magma_tally3_int_t lddb,
    const char* func, const char* file, int line );

void
magma_tally3_ssetmatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float const    *hA_src, magma_tally3_int_t ldha,
    magma_tally3Float_ptr       dB_dst, magma_tally3_int_t lddb,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally3_sgetmatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA_src, magma_tally3_int_t ldda,
    float          *hB_dst, magma_tally3_int_t ldhb,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally3_scopymatrix_async_internal(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA_src, magma_tally3_int_t ldda,
    magma_tally3Float_ptr       dB_dst, magma_tally3_int_t lddb,
    magma_tally3_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// Level 1 BLAS (alphabetical order)

// in cublas_v2, result returned through output argument
magma_tally3_int_t
magma_tally3_isamax(
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx );

// in cublas_v2, result returned through output argument
magma_tally3_int_t
magma_tally3_isamin(
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx );

// in cublas_v2, result returned through output argument
float
magma_tally3_sasum(
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx );

void
magma_tally3_saxpy(
    magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_ptr       dy, magma_tally3_int_t incy );

void
magma_tally3_scopy(
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_ptr       dy, magma_tally3_int_t incy );

// in cublas_v2, result returned through output argument
float
magma_tally3_sdot(
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_const_ptr dy, magma_tally3_int_t incy );

// in cublas_v2, result returned through output argument
float
magma_tally3_sdot(
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_const_ptr dy, magma_tally3_int_t incy );

// in cublas_v2, result returned through output argument
float
magma_tally3_snrm2(
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx );

void
magma_tally3_srot(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_ptr dy, magma_tally3_int_t incy,
    float dc, float ds );

void
magma_tally3_srot(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_ptr dy, magma_tally3_int_t incy,
    float dc, float ds );

#ifdef REAL
void
magma_tally3_srotm(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_ptr dy, magma_tally3_int_t incy,
    magma_tally3Float_const_ptr param );

void
magma_tally3_srotmg(
    magma_tally3Float_ptr d1, magma_tally3Float_ptr       d2,
    magma_tally3Float_ptr x1, magma_tally3Float_const_ptr y1,
    magma_tally3Float_ptr param );
#endif

void
magma_tally3_sscal(
    magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_ptr dx, magma_tally3_int_t incx );

void
magma_tally3_sscal(
    magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_ptr dx, magma_tally3_int_t incx );

void
magma_tally3_sswap(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_ptr dy, magma_tally3_int_t incy );

// ========================================
// Level 2 BLAS (alphabetical order)

void
magma_tally3_sgemv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx,
    float beta,
    magma_tally3Float_ptr       dy, magma_tally3_int_t incy );

void
magma_tally3_sger(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_const_ptr dy, magma_tally3_int_t incy,
    magma_tally3Float_ptr       dA, magma_tally3_int_t ldda );

void
magma_tally3_sger(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_const_ptr dy, magma_tally3_int_t incy,
    magma_tally3Float_ptr       dA, magma_tally3_int_t ldda );

void
magma_tally3_ssymv(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx,
    float beta,
    magma_tally3Float_ptr       dy, magma_tally3_int_t incy );

void
magma_tally3_ssyr(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_ptr       dA, magma_tally3_int_t ldda );

void
magma_tally3_ssyr2(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_const_ptr dy, magma_tally3_int_t incy,
    magma_tally3Float_ptr       dA, magma_tally3_int_t ldda );

void
magma_tally3_strmv(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_diag_t diag,
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr       dx, magma_tally3_int_t incx );

void
magma_tally3_strsv(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_diag_t diag,
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr       dx, magma_tally3_int_t incx );

// ========================================
// Level 3 BLAS (alphabetical order)

void
magma_tally3_sgemm(
    magma_tally3_trans_t transA, magma_tally3_trans_t transB,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dB, magma_tally3_int_t lddb,
    float beta,
    magma_tally3Float_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3_ssymm(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dB, magma_tally3_int_t lddb,
    float beta,
    magma_tally3Float_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3_ssymm(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dB, magma_tally3_int_t lddb,
    float beta,
    magma_tally3Float_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3_ssyr2k(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dB, magma_tally3_int_t lddb,
    float beta,
    magma_tally3Float_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3_ssyr2k(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_const_ptr dB, magma_tally3_int_t lddb,
    float beta,
    magma_tally3Float_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3_ssyrk(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    float beta,
    magma_tally3Float_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3_ssyrk(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    float beta,
    magma_tally3Float_ptr       dC, magma_tally3_int_t lddc );

void
magma_tally3_strmm(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr       dB, magma_tally3_int_t lddb );

void
magma_tally3_strsm(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr       dB, magma_tally3_int_t lddb );


#ifdef __cplusplus
}
#endif

#undef REAL

#endif  /* MAGMA_tally3BLAS_S_H */
