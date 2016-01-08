/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally4blas_z.h normal z -> s, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally4BLAS_S_H
#define MAGMA_tally4BLAS_S_H

#include "magma_tally4_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_tally4blas_stranspose_inplace(
    magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda );

void
magma_tally4blas_stranspose(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4Float_ptr       dAT, magma_tally4_int_t lddat );

void
magma_tally4blas_sgetmatrix_transpose(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr dAT,   magma_tally4_int_t ldda,
    float          *hA,    magma_tally4_int_t lda,
    magma_tally4Float_ptr       dwork, magma_tally4_int_t lddwork, magma_tally4_int_t nb );

void
magma_tally4blas_ssetmatrix_transpose(
    magma_tally4_int_t m, magma_tally4_int_t n,
    const float *hA,    magma_tally4_int_t lda,
    magma_tally4Float_ptr    dAT,   magma_tally4_int_t ldda,
    magma_tally4Float_ptr    dwork, magma_tally4_int_t lddwork, magma_tally4_int_t nb );

  /*
   * RBT-related functions
   */
void
magma_tally4blas_sprbt(
    magma_tally4_int_t n, 
    float *dA, magma_tally4_int_t ldda, 
    float *du, float *dv);

void
magma_tally4blas_sprbt_mv(
    magma_tally4_int_t n, 
    float *dv, float *db);

void
magma_tally4blas_sprbt_mtv(
    magma_tally4_int_t n, 
    float *du, float *db);

void
magma_tally4blas_saxpycp2(
    magma_tally4_int_t m, float *r, float *x,
    const float *b);

  /*
   * Multi-GPU copy functions
   */
void
magma_tally4blas_sgetmatrix_transpose_mgpu(
    magma_tally4_int_t ngpu, magma_tally4_queue_t queues[][2],
    magma_tally4Float_const_ptr const dAT[],    magma_tally4_int_t ldda,
    float                *hA,       magma_tally4_int_t lda,
    magma_tally4Float_ptr             dwork[],  magma_tally4_int_t lddw,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nb );

void
magma_tally4blas_ssetmatrix_transpose_mgpu(
    magma_tally4_int_t ngpu, magma_tally4_queue_t queues[][2],
    const float *hA,      magma_tally4_int_t lda,
    magma_tally4Float_ptr    dAT[],   magma_tally4_int_t ldda,
    magma_tally4Float_ptr    dwork[], magma_tally4_int_t lddw,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nb );

void
magma_tally4_sgetmatrix_1D_col_bcyclic(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr const dA[], magma_tally4_int_t ldda,
    float                *hA,   magma_tally4_int_t lda,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb );

void
magma_tally4_ssetmatrix_1D_col_bcyclic(
    magma_tally4_int_t m, magma_tally4_int_t n,
    const float *hA,   magma_tally4_int_t lda,
    magma_tally4Float_ptr    dA[], magma_tally4_int_t ldda,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb );

void
magma_tally4_sgetmatrix_1D_row_bcyclic(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr const dA[], magma_tally4_int_t ldda,
    float                *hA,   magma_tally4_int_t lda,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb );

void
magma_tally4_ssetmatrix_1D_row_bcyclic(
    magma_tally4_int_t m, magma_tally4_int_t n,
    const float *hA,   magma_tally4_int_t lda,
    magma_tally4Float_ptr    dA[], magma_tally4_int_t ldda,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb );

// in src/ssytrd_mgpu.cpp
// TODO rename ssetmatrix_sy or similar
magma_tally4_int_t
magma_tally4_shtodhe(
    magma_tally4_int_t ngpu, magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    float     *A,   magma_tally4_int_t lda,
    magma_tally4Float_ptr dA[], magma_tally4_int_t ldda,
    magma_tally4_queue_t queues[][10],
    magma_tally4_int_t *info );

// in src/spotrf3_mgpu.cpp
// TODO same as magma_tally4_shtodhe?
magma_tally4_int_t
magma_tally4_shtodpo(
    magma_tally4_int_t ngpu, magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t off_i, magma_tally4_int_t off_j, magma_tally4_int_t nb,
    float     *A,   magma_tally4_int_t lda,
    magma_tally4Float_ptr dA[], magma_tally4_int_t ldda,
    magma_tally4_queue_t queues[][3],
    magma_tally4_int_t *info );

// in src/spotrf3_mgpu.cpp
// TODO rename sgetmatrix_sy or similar
magma_tally4_int_t
magma_tally4_sdtohpo(
    magma_tally4_int_t ngpu, magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t off_i, magma_tally4_int_t off_j, magma_tally4_int_t nb, magma_tally4_int_t NB,
    float     *A,   magma_tally4_int_t lda,
    magma_tally4Float_ptr dA[], magma_tally4_int_t ldda,
    magma_tally4_queue_t queues[][3],
    magma_tally4_int_t *info );


  /*
   * Multi-GPU BLAS functions (alphabetical order)
   */
void
magma_tally4blas_ssymm_mgpu(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_ptr dA[],    magma_tally4_int_t ldda,  magma_tally4_int_t offset,
    magma_tally4Float_ptr dB[],    magma_tally4_int_t lddb,
    float beta,
    magma_tally4Float_ptr dC[],    magma_tally4_int_t lddc,
    magma_tally4Float_ptr dwork[], magma_tally4_int_t dworksiz,
    float    *C,       magma_tally4_int_t ldc,
    float    *work[],  magma_tally4_int_t worksiz,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_event_t redevents[][20], magma_tally4_int_t nbevents );

void
magma_tally4blas_ssymm_mgpu_com(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_ptr dA[],    magma_tally4_int_t ldda,  magma_tally4_int_t offset,
    magma_tally4Float_ptr dB[],    magma_tally4_int_t lddb,
    float beta,
    magma_tally4Float_ptr dC[],    magma_tally4_int_t lddc,
    magma_tally4Float_ptr dwork[], magma_tally4_int_t dworksiz,
    float    *C,       magma_tally4_int_t ldc,
    float    *work[],  magma_tally4_int_t worksiz,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_event_t redevents[][Magma_tally4MaxGPUs*Magma_tally4MaxGPUs+10], magma_tally4_int_t nbevents,
    magma_tally4_int_t gnode[Magma_tally4MaxGPUs][Magma_tally4MaxGPUs+2], magma_tally4_int_t nbcmplx );

void
magma_tally4blas_ssymm_mgpu_spec(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_ptr dA[],    magma_tally4_int_t ldda,  magma_tally4_int_t offset,
    magma_tally4Float_ptr dB[],    magma_tally4_int_t lddb,
    float beta,
    magma_tally4Float_ptr dC[],    magma_tally4_int_t lddc,
    magma_tally4Float_ptr dwork[], magma_tally4_int_t dworksiz,
    float    *C,       magma_tally4_int_t ldc,
    float    *work[],  magma_tally4_int_t worksiz,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_event_t redevents[][Magma_tally4MaxGPUs*Magma_tally4MaxGPUs+10], magma_tally4_int_t nbevents,
    magma_tally4_int_t gnode[Magma_tally4MaxGPUs][Magma_tally4MaxGPUs+2], magma_tally4_int_t nbcmplx );

void
magma_tally4blas_ssymm_mgpu_spec33(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_ptr dA[],    magma_tally4_int_t ldda,  magma_tally4_int_t offset,
    magma_tally4Float_ptr dB[],    magma_tally4_int_t lddb,
    float beta,
    magma_tally4Float_ptr dC[],    magma_tally4_int_t lddc,
    magma_tally4Float_ptr dVIN[],  magma_tally4_int_t lddv, magma_tally4_int_t v_offset,
    magma_tally4Float_ptr dwork[], magma_tally4_int_t dworksiz,
    float *C,       magma_tally4_int_t ldc,
    float *work[],  magma_tally4_int_t worksiz,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_event_t redevents[][Magma_tally4MaxGPUs*Magma_tally4MaxGPUs+10], magma_tally4_int_t nbevents,
    magma_tally4_int_t gnode[Magma_tally4MaxGPUs][Magma_tally4MaxGPUs+2], magma_tally4_int_t nbcmplx );

magma_tally4_int_t
magma_tally4blas_ssymv_mgpu(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr const d_lA[], magma_tally4_int_t ldda, magma_tally4_int_t offset,
    float const *x,         magma_tally4_int_t incx,
    float beta,             
    float       *y,         magma_tally4_int_t incy,
    float       *hwork,     magma_tally4_int_t lhwork,
    magma_tally4Float_ptr    dwork[],   magma_tally4_int_t ldwork,
    magma_tally4_int_t ngpu,
    magma_tally4_int_t nb,
    magma_tally4_queue_t queues[] );

magma_tally4_int_t
magma_tally4blas_ssymv_mgpu_sync(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr const d_lA[], magma_tally4_int_t ldda, magma_tally4_int_t offset,
    float const *x,         magma_tally4_int_t incx,
    float beta,             
    float       *y,         magma_tally4_int_t incy,
    float       *hwork,     magma_tally4_int_t lhwork,
    magma_tally4Float_ptr    dwork[],   magma_tally4_int_t ldwork,
    magma_tally4_int_t ngpu,
    magma_tally4_int_t nb,
    magma_tally4_queue_t queues[] );

// Ichi's version, in src/ssytrd_mgpu.cpp
void
magma_tally4_ssyr2k_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t nb, magma_tally4_int_t n, magma_tally4_int_t k,
    float alpha,
    magma_tally4Float_ptr dB[], magma_tally4_int_t lddb, magma_tally4_int_t b_offset,
    float beta,
    magma_tally4Float_ptr dC[], magma_tally4_int_t lddc, magma_tally4_int_t c_offset,
    magma_tally4_int_t nqueue, magma_tally4_queue_t queues[][10] );

void
magma_tally4blas_ssyr2k_mgpu_spec(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t k,
    float alpha,
    magma_tally4Float_ptr dA[], magma_tally4_int_t ldda, magma_tally4_int_t a_offset,
    magma_tally4Float_ptr dB[], magma_tally4_int_t lddb, magma_tally4_int_t b_offset,
    float beta,
    magma_tally4Float_ptr dC[], magma_tally4_int_t lddc, magma_tally4_int_t c_offset,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb, magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue );

void
magma_tally4blas_ssyr2k_mgpu_spec324(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_ptr dVIN[], magma_tally4_int_t lddv, magma_tally4_int_t v_offset,
    magma_tally4Float_ptr dWIN[], magma_tally4_int_t lddw, magma_tally4_int_t w_offset,
    float beta,
    magma_tally4Float_ptr dC[],   magma_tally4_int_t lddc, magma_tally4_int_t c_offset,
    magma_tally4Float_ptr dwork[], magma_tally4_int_t lndwork,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_event_t redevents[][Magma_tally4MaxGPUs*Magma_tally4MaxGPUs+10], magma_tally4_int_t nbevents );

void
magma_tally4blas_ssyr2k_mgpu_spec325(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_ptr dVIN[], magma_tally4_int_t lddv, magma_tally4_int_t v_offset,
    magma_tally4Float_ptr dWIN[], magma_tally4_int_t lddw, magma_tally4_int_t w_offset,
    float beta,
    magma_tally4Float_ptr dC[],   magma_tally4_int_t lddc, magma_tally4_int_t c_offset,
    magma_tally4Float_ptr dwork[], magma_tally4_int_t lndwork,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb,
    float *harray[],
    magma_tally4Float_ptr *darray[],
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_event_t redevents[][Magma_tally4MaxGPUs*Magma_tally4MaxGPUs+10], magma_tally4_int_t nbevents );

void
magma_tally4blas_ssyr2k_mgpu2(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t k,
    float alpha,
    magma_tally4Float_ptr dA[], magma_tally4_int_t ldda, magma_tally4_int_t a_offset,
    magma_tally4Float_ptr dB[], magma_tally4_int_t lddb, magma_tally4_int_t b_offset,
    float beta,
    magma_tally4Float_ptr dC[], magma_tally4_int_t lddc, magma_tally4_int_t c_offset,
    magma_tally4_int_t ngpu, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue );

// in src/spotrf_mgpu_right.cpp
void
magma_tally4_ssyrk_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t nb, magma_tally4_int_t n, magma_tally4_int_t k,
    float alpha,
    magma_tally4Float_ptr dB[], magma_tally4_int_t lddb, magma_tally4_int_t b_offset,
    float beta,
    magma_tally4Float_ptr dC[], magma_tally4_int_t lddc, magma_tally4_int_t c_offset,
    magma_tally4_int_t nqueue, magma_tally4_queue_t queues[][10]);

// in src/spotrf_mgpu_right.cpp
void
magma_tally4_ssyrk_mgpu2(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t nb, magma_tally4_int_t n, magma_tally4_int_t k,
    float alpha,
    magma_tally4Float_ptr dB[], magma_tally4_int_t lddb, magma_tally4_int_t b_offset,
    float beta,
    magma_tally4Float_ptr dC[], magma_tally4_int_t lddc, magma_tally4_int_t c_offset,
    magma_tally4_int_t nqueue, magma_tally4_queue_t queues[][10]);


  /*
   * LAPACK auxiliary functions (alphabetical order)
   */
void
magma_tally4blas_sgeadd(
    magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr       dB, magma_tally4_int_t lddb );

void
magma_tally4blas_slacpy(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr       dB, magma_tally4_int_t lddb );

void
magma_tally4blas_slacpy_cnjg(
    magma_tally4_int_t n, float *dA1, magma_tally4_int_t lda1,
    float *dA2, magma_tally4_int_t lda2);

float
magma_tally4blas_slange(
    magma_tally4_norm_t norm,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dwork );

float
magma_tally4blas_slansy(
    magma_tally4_norm_t norm, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dwork );

float
magma_tally4blas_slansy(
    magma_tally4_norm_t norm, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dwork );

void
magma_tally4blas_slarfg(
    magma_tally4_int_t n,
    magma_tally4Float_ptr dalpha, magma_tally4Float_ptr dx, magma_tally4_int_t incx,
    magma_tally4Float_ptr dtau );

void
magma_tally4blas_slarfg_work(
    magma_tally4_int_t n,
    magma_tally4Float_ptr dalpha, magma_tally4Float_ptr dx, magma_tally4_int_t incx,
    magma_tally4Float_ptr dtau, magma_tally4Float_ptr dwork );

void
magma_tally4blas_slascl(
    magma_tally4_type_t type, magma_tally4_int_t kl, magma_tally4_int_t ku,
    float cfrom, float cto,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info );

void
magma_tally4blas_slascl_2x2(
    magma_tally4_type_t type, magma_tally4_int_t m,
    float *dW, magma_tally4_int_t lddw,
    float *dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info );

void
magma_tally4blas_slascl2(
    magma_tally4_type_t type,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr dD,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info );

void
magma_tally4blas_slascl_diag(
    magma_tally4_type_t type, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr dD, magma_tally4_int_t lddd,
          magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info );

void
magma_tally4blas_slaset(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    float offdiag, float diag,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda );

void
magma_tally4blas_slaset_band(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    float offdiag, float diag,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda);

void
magma_tally4blas_slaswp(
    magma_tally4_int_t n,
    magma_tally4Float_ptr dAT, magma_tally4_int_t ldda,
    magma_tally4_int_t k1, magma_tally4_int_t k2,
    const magma_tally4_int_t *ipiv, magma_tally4_int_t inci );

void
magma_tally4blas_slaswp2(
    magma_tally4_int_t n,
    magma_tally4Float_ptr dAT, magma_tally4_int_t ldda,
    magma_tally4_int_t k1, magma_tally4_int_t k2,
    magma_tally4Int_const_ptr d_ipiv, magma_tally4_int_t inci );

void
magma_tally4blas_slaswpx(
    magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldx, magma_tally4_int_t ldy,
    magma_tally4_int_t k1, magma_tally4_int_t k2,
    const magma_tally4_int_t *ipiv, magma_tally4_int_t inci );

void
magma_tally4blas_ssymmetrize(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda );

void
magma_tally4blas_ssymmetrize_tiles(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t ntile, magma_tally4_int_t mstride, magma_tally4_int_t nstride );

void
magma_tally4blas_strtri_diag(
    magma_tally4_uplo_t uplo, magma_tally4_diag_t diag, magma_tally4_int_t n,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr d_dinvA);

  /*
   * to cleanup (alphabetical order)
   */
void
magma_tally4blas_snrm2_adjust(
    magma_tally4_int_t k,
    magma_tally4Float_ptr dxnorm,
    magma_tally4Float_ptr dc);

void
magma_tally4blas_snrm2_check(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dxnorm,
    magma_tally4Float_ptr dlsticc);

void
magma_tally4blas_snrm2_cols(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dxnorm);

void
magma_tally4blas_snrm2_row_check_adjust(
    magma_tally4_int_t k, float tol,
    magma_tally4Float_ptr dxnorm,
    magma_tally4Float_ptr dxnorm2,
    magma_tally4Float_ptr dC, magma_tally4_int_t lddc,
    magma_tally4Float_ptr dlsticc);

void
magma_tally4_slarfbx_gpu(
    magma_tally4_int_t m, magma_tally4_int_t k,
    magma_tally4Float_ptr V,  magma_tally4_int_t ldv,
    magma_tally4Float_ptr dT, magma_tally4_int_t ldt,
    magma_tally4Float_ptr c,
    magma_tally4Float_ptr dwork);

void
magma_tally4_slarfg_gpu(
    magma_tally4_int_t n,
    magma_tally4Float_ptr dx0,
    magma_tally4Float_ptr dx,
    magma_tally4Float_ptr dtau,
    magma_tally4Float_ptr        dxnorm,
    magma_tally4Float_ptr dAkk );

void
magma_tally4_slarfgtx_gpu(
    magma_tally4_int_t n,
    magma_tally4Float_ptr dx0,
    magma_tally4Float_ptr dx,
    magma_tally4Float_ptr dtau,
    magma_tally4Float_ptr        dxnorm,
    magma_tally4Float_ptr dA, magma_tally4_int_t iter,
    magma_tally4Float_ptr V,  magma_tally4_int_t ldv,
    magma_tally4Float_ptr T,  magma_tally4_int_t ldt,
    magma_tally4Float_ptr dwork);

void
magma_tally4_slarfgx_gpu(
    magma_tally4_int_t n,
    magma_tally4Float_ptr dx0,
    magma_tally4Float_ptr dx,
    magma_tally4Float_ptr dtau,
    magma_tally4Float_ptr        dxnorm,
    magma_tally4Float_ptr dA, magma_tally4_int_t iter);

void
magma_tally4_slarfx_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr v,
    magma_tally4Float_ptr tau,
    magma_tally4Float_ptr C,  magma_tally4_int_t ldc,
    magma_tally4Float_ptr        xnorm,
    magma_tally4Float_ptr dT, magma_tally4_int_t iter,
    magma_tally4Float_ptr work);


  /*
   * Level 1 BLAS (alphabetical order)
   */
void
magma_tally4blas_sswap(
    magma_tally4_int_t n,
    magma_tally4Float_ptr dx, magma_tally4_int_t incx,
    magma_tally4Float_ptr dy, magma_tally4_int_t incy );

void
magma_tally4blas_sswapblk(
    magma_tally4_order_t order,
    magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t i1, magma_tally4_int_t i2,
    const magma_tally4_int_t *ipiv, magma_tally4_int_t inci,
    magma_tally4_int_t offset );

void
magma_tally4blas_sswapdblk(
    magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t inca,
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb, magma_tally4_int_t incb );

  /*
   * Level 2 BLAS (alphabetical order)
   */
void
magma_tally4blas_sgemv(
    magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx,
    float beta,
    magma_tally4Float_ptr       dy, magma_tally4_int_t incy );


void
magma_tally4blas_sgemv_conjv(
    magma_tally4_int_t m, magma_tally4_int_t n, float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx,
    float beta,
    magma_tally4Float_ptr dy, magma_tally4_int_t incy);

magma_tally4_int_t
magma_tally4blas_ssymv(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx,
    float beta,
    magma_tally4Float_ptr       dy, magma_tally4_int_t incy );

magma_tally4_int_t
magma_tally4blas_ssymv(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx,
    float beta,
    magma_tally4Float_ptr       dy, magma_tally4_int_t incy );

magma_tally4_int_t
magma_tally4blas_ssymv_work(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx,
    float beta,
    magma_tally4Float_ptr       dy, magma_tally4_int_t incy,
    magma_tally4Float_ptr       dwork, magma_tally4_int_t lwork,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4blas_ssymv_work(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx,
    float beta,
    magma_tally4Float_ptr       dy, magma_tally4_int_t incy,
    magma_tally4Float_ptr       dwork, magma_tally4_int_t lwork,
    magma_tally4_queue_t queue );

  /*
   * Level 3 BLAS (alphabetical order)
   */
void
magma_tally4blas_sgemm(
    magma_tally4_trans_t transA, magma_tally4_trans_t transB,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dB, magma_tally4_int_t lddb,
    float beta,
    magma_tally4Float_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4blas_sgemm_reduce(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dB, magma_tally4_int_t lddb,
    float beta,
    magma_tally4Float_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4blas_ssymm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo,
    magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dB, magma_tally4_int_t lddb,
    float beta,
    magma_tally4Float_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4blas_ssymm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo,
    magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dB, magma_tally4_int_t lddb,
    float beta,
    magma_tally4Float_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4blas_ssyr2k(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dB, magma_tally4_int_t lddb,
    float beta,
    magma_tally4Float_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4blas_ssyr2k(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dB, magma_tally4_int_t lddb,
    float  beta,
    magma_tally4Float_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4blas_ssyrk(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    float beta,
    magma_tally4Float_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4blas_ssyrk(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    float  alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    float  beta,
    magma_tally4Float_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4blas_strsm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t transA, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr       dB, magma_tally4_int_t lddb );

void
magma_tally4blas_strsm_outofplace(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t transA, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr       dB, magma_tally4_int_t lddb,
    magma_tally4_int_t flag, magma_tally4Float_ptr d_dinvA, magma_tally4Float_ptr dX );

void
magma_tally4blas_strsm_work(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t transA, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr       dB, magma_tally4_int_t lddb,
    magma_tally4_int_t flag, magma_tally4Float_ptr d_dinvA, magma_tally4Float_ptr dX );


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

#define magma_tally4_ssetvector(           n, hx_src, incx, dy_dst, incy ) \
        magma_tally4_ssetvector_internal(  n, hx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally4_sgetvector(           n, dx_src, incx, hy_dst, incy ) \
        magma_tally4_sgetvector_internal(  n, dx_src, incx, hy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally4_scopyvector(          n, dx_src, incx, dy_dst, incy ) \
        magma_tally4_scopyvector_internal( n, dx_src, incx, dy_dst, incy, __func__, __FILE__, __LINE__ )

#define magma_tally4_ssetvector_async(           n, hx_src, incx, dy_dst, incy, queue ) \
        magma_tally4_ssetvector_async_internal(  n, hx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_tally4_sgetvector_async(           n, dx_src, incx, hy_dst, incy, queue ) \
        magma_tally4_sgetvector_async_internal(  n, dx_src, incx, hy_dst, incy, queue, __func__, __FILE__, __LINE__ )

#define magma_tally4_scopyvector_async(          n, dx_src, incx, dy_dst, incy, queue ) \
        magma_tally4_scopyvector_async_internal( n, dx_src, incx, dy_dst, incy, queue, __func__, __FILE__, __LINE__ )

void
magma_tally4_ssetvector_internal(
    magma_tally4_int_t n,
    float const    *hx_src, magma_tally4_int_t incx,
    magma_tally4Float_ptr       dy_dst, magma_tally4_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally4_sgetvector_internal(
    magma_tally4_int_t n,
    magma_tally4Float_const_ptr dx_src, magma_tally4_int_t incx,
    float          *hy_dst, magma_tally4_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally4_scopyvector_internal(
    magma_tally4_int_t n,
    magma_tally4Float_const_ptr dx_src, magma_tally4_int_t incx,
    magma_tally4Float_ptr       dy_dst, magma_tally4_int_t incy,
    const char* func, const char* file, int line );

void
magma_tally4_ssetvector_async_internal(
    magma_tally4_int_t n,
    float const    *hx_src, magma_tally4_int_t incx,
    magma_tally4Float_ptr       dy_dst, magma_tally4_int_t incy,
    magma_tally4_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally4_sgetvector_async_internal(
    magma_tally4_int_t n,
    magma_tally4Float_const_ptr dx_src, magma_tally4_int_t incx,
    float          *hy_dst, magma_tally4_int_t incy,
    magma_tally4_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally4_scopyvector_async_internal(
    magma_tally4_int_t n,
    magma_tally4Float_const_ptr dx_src, magma_tally4_int_t incx,
    magma_tally4Float_ptr       dy_dst, magma_tally4_int_t incy,
    magma_tally4_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// copying sub-matrices (contiguous columns)
// set  copies host to device
// get  copies device to host
// copy copies device to device
// (with CUDA unified addressing, copy can be between same or different devices)
// Add the function, file, and line for error-reporting purposes.

#define magma_tally4_ssetmatrix(           m, n, hA_src, lda, dB_dst, lddb ) \
        magma_tally4_ssetmatrix_internal(  m, n, hA_src, lda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_tally4_sgetmatrix(           m, n, dA_src, ldda, hB_dst, ldb ) \
        magma_tally4_sgetmatrix_internal(  m, n, dA_src, ldda, hB_dst, ldb, __func__, __FILE__, __LINE__ )

#define magma_tally4_scopymatrix(          m, n, dA_src, ldda, dB_dst, lddb ) \
        magma_tally4_scopymatrix_internal( m, n, dA_src, ldda, dB_dst, lddb, __func__, __FILE__, __LINE__ )

#define magma_tally4_ssetmatrix_async(           m, n, hA_src, lda, dB_dst, lddb, queue ) \
        magma_tally4_ssetmatrix_async_internal(  m, n, hA_src, lda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

#define magma_tally4_sgetmatrix_async(           m, n, dA_src, ldda, hB_dst, ldb, queue ) \
        magma_tally4_sgetmatrix_async_internal(  m, n, dA_src, ldda, hB_dst, ldb, queue, __func__, __FILE__, __LINE__ )

#define magma_tally4_scopymatrix_async(          m, n, dA_src, ldda, dB_dst, lddb, queue ) \
        magma_tally4_scopymatrix_async_internal( m, n, dA_src, ldda, dB_dst, lddb, queue, __func__, __FILE__, __LINE__ )

void
magma_tally4_ssetmatrix_internal(
    magma_tally4_int_t m, magma_tally4_int_t n,
    float const    *hA_src, magma_tally4_int_t ldha,
    magma_tally4Float_ptr       dB_dst, magma_tally4_int_t lddb,
    const char* func, const char* file, int line );

void
magma_tally4_sgetmatrix_internal(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr dA_src, magma_tally4_int_t ldda,
    float          *hB_dst, magma_tally4_int_t ldhb,
    const char* func, const char* file, int line );

void
magma_tally4_scopymatrix_internal(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr dA_src, magma_tally4_int_t ldda,
    magma_tally4Float_ptr       dB_dst, magma_tally4_int_t lddb,
    const char* func, const char* file, int line );

void
magma_tally4_ssetmatrix_async_internal(
    magma_tally4_int_t m, magma_tally4_int_t n,
    float const    *hA_src, magma_tally4_int_t ldha,
    magma_tally4Float_ptr       dB_dst, magma_tally4_int_t lddb,
    magma_tally4_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally4_sgetmatrix_async_internal(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr dA_src, magma_tally4_int_t ldda,
    float          *hB_dst, magma_tally4_int_t ldhb,
    magma_tally4_queue_t queue,
    const char* func, const char* file, int line );

void
magma_tally4_scopymatrix_async_internal(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr dA_src, magma_tally4_int_t ldda,
    magma_tally4Float_ptr       dB_dst, magma_tally4_int_t lddb,
    magma_tally4_queue_t queue,
    const char* func, const char* file, int line );


// ========================================
// Level 1 BLAS (alphabetical order)

// in cublas_v2, result returned through output argument
magma_tally4_int_t
magma_tally4_isamax(
    magma_tally4_int_t n,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx );

// in cublas_v2, result returned through output argument
magma_tally4_int_t
magma_tally4_isamin(
    magma_tally4_int_t n,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx );

// in cublas_v2, result returned through output argument
float
magma_tally4_sasum(
    magma_tally4_int_t n,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx );

void
magma_tally4_saxpy(
    magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4Float_ptr       dy, magma_tally4_int_t incy );

void
magma_tally4_scopy(
    magma_tally4_int_t n,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4Float_ptr       dy, magma_tally4_int_t incy );

// in cublas_v2, result returned through output argument
float
magma_tally4_sdot(
    magma_tally4_int_t n,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4Float_const_ptr dy, magma_tally4_int_t incy );

// in cublas_v2, result returned through output argument
float
magma_tally4_sdot(
    magma_tally4_int_t n,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4Float_const_ptr dy, magma_tally4_int_t incy );

// in cublas_v2, result returned through output argument
float
magma_tally4_snrm2(
    magma_tally4_int_t n,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx );

void
magma_tally4_srot(
    magma_tally4_int_t n,
    magma_tally4Float_ptr dx, magma_tally4_int_t incx,
    magma_tally4Float_ptr dy, magma_tally4_int_t incy,
    float dc, float ds );

void
magma_tally4_srot(
    magma_tally4_int_t n,
    magma_tally4Float_ptr dx, magma_tally4_int_t incx,
    magma_tally4Float_ptr dy, magma_tally4_int_t incy,
    float dc, float ds );

#ifdef REAL
void
magma_tally4_srotm(
    magma_tally4_int_t n,
    magma_tally4Float_ptr dx, magma_tally4_int_t incx,
    magma_tally4Float_ptr dy, magma_tally4_int_t incy,
    magma_tally4Float_const_ptr param );

void
magma_tally4_srotmg(
    magma_tally4Float_ptr d1, magma_tally4Float_ptr       d2,
    magma_tally4Float_ptr x1, magma_tally4Float_const_ptr y1,
    magma_tally4Float_ptr param );
#endif

void
magma_tally4_sscal(
    magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_ptr dx, magma_tally4_int_t incx );

void
magma_tally4_sscal(
    magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_ptr dx, magma_tally4_int_t incx );

void
magma_tally4_sswap(
    magma_tally4_int_t n,
    magma_tally4Float_ptr dx, magma_tally4_int_t incx,
    magma_tally4Float_ptr dy, magma_tally4_int_t incy );

// ========================================
// Level 2 BLAS (alphabetical order)

void
magma_tally4_sgemv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx,
    float beta,
    magma_tally4Float_ptr       dy, magma_tally4_int_t incy );

void
magma_tally4_sger(
    magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4Float_const_ptr dy, magma_tally4_int_t incy,
    magma_tally4Float_ptr       dA, magma_tally4_int_t ldda );

void
magma_tally4_sger(
    magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4Float_const_ptr dy, magma_tally4_int_t incy,
    magma_tally4Float_ptr       dA, magma_tally4_int_t ldda );

void
magma_tally4_ssymv(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx,
    float beta,
    magma_tally4Float_ptr       dy, magma_tally4_int_t incy );

void
magma_tally4_ssyr(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4Float_ptr       dA, magma_tally4_int_t ldda );

void
magma_tally4_ssyr2(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4Float_const_ptr dy, magma_tally4_int_t incy,
    magma_tally4Float_ptr       dA, magma_tally4_int_t ldda );

void
magma_tally4_strmv(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_diag_t diag,
    magma_tally4_int_t n,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr       dx, magma_tally4_int_t incx );

void
magma_tally4_strsv(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_diag_t diag,
    magma_tally4_int_t n,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr       dx, magma_tally4_int_t incx );

// ========================================
// Level 3 BLAS (alphabetical order)

void
magma_tally4_sgemm(
    magma_tally4_trans_t transA, magma_tally4_trans_t transB,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dB, magma_tally4_int_t lddb,
    float beta,
    magma_tally4Float_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4_ssymm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo,
    magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dB, magma_tally4_int_t lddb,
    float beta,
    magma_tally4Float_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4_ssymm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo,
    magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dB, magma_tally4_int_t lddb,
    float beta,
    magma_tally4Float_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4_ssyr2k(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dB, magma_tally4_int_t lddb,
    float beta,
    magma_tally4Float_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4_ssyr2k(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_const_ptr dB, magma_tally4_int_t lddb,
    float beta,
    magma_tally4Float_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4_ssyrk(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    float beta,
    magma_tally4Float_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4_ssyrk(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    float beta,
    magma_tally4Float_ptr       dC, magma_tally4_int_t lddc );

void
magma_tally4_strmm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr       dB, magma_tally4_int_t lddb );

void
magma_tally4_strsm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr       dB, magma_tally4_int_t lddb );


#ifdef __cplusplus
}
#endif

#undef REAL

#endif  /* MAGMA_tally4BLAS_S_H */
