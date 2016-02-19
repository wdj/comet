/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally3blas_z_q.h normal z -> s, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_tally3BLAS_S_Q_H
#define MAGMA_tally3BLAS_S_Q_H
                    
#include "magma_tally3_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_tally3blas_stranspose_inplace_q(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_queue_t queue );

void
magma_tally3blas_stranspose_q(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3Float_ptr       dAT, magma_tally3_int_t lddat,
    magma_tally3_queue_t queue );

void
magma_tally3blas_sgetmatrix_transpose_q(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dAT,   magma_tally3_int_t ldda,
    float          *hA,    magma_tally3_int_t lda,
    magma_tally3Float_ptr       dwork, magma_tally3_int_t lddwork, magma_tally3_int_t nb,
    magma_tally3_queue_t queues[2] );

void
magma_tally3blas_ssetmatrix_transpose_q(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const float *hA,    magma_tally3_int_t lda,
    magma_tally3Float_ptr    dAT,   magma_tally3_int_t ldda,
    magma_tally3Float_ptr    dwork, magma_tally3_int_t lddwork, magma_tally3_int_t nb,
    magma_tally3_queue_t queues[2] );

  /*
   * Multi-GPU functions
   */



  /*
   * LAPACK auxiliary functions
   */
void
magma_tally3blas_sgeadd_q(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr       dB, magma_tally3_int_t lddb,
    magma_tally3_queue_t queue );

void
magma_tally3blas_slacpy_q(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr       dB, magma_tally3_int_t lddb,
    magma_tally3_queue_t queue );

float
magma_tally3blas_slange_q(
    magma_tally3_norm_t norm,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dwork,
    magma_tally3_queue_t queue );

float
magma_tally3blas_slansy_q(
    magma_tally3_norm_t norm, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dwork,
    magma_tally3_queue_t queue );

float
magma_tally3blas_slansy_q(
    magma_tally3_norm_t norm, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dwork,
    magma_tally3_queue_t queue );

void
magma_tally3blas_slarfg_q(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dalpha,
    magma_tally3Float_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_ptr dtau,
    magma_tally3_queue_t queue );

void
magma_tally3blas_slascl_q(
    magma_tally3_type_t type, magma_tally3_int_t kl, magma_tally3_int_t ku,
    float cfrom, float cto,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_queue_t queue,
    magma_tally3_int_t *info );

void
magma_tally3blas_slascl2_q(
    magma_tally3_type_t type,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dD,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_queue_t queue,
    magma_tally3_int_t *info );

void
magma_tally3blas_slaset_q(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    float offdiag, float diag,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_queue_t queue );

void
magma_tally3blas_slaset_band_q(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    float offdiag, float diag,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_queue_t queue);

void
magma_tally3blas_slaswp_q(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dAT, magma_tally3_int_t ldda,
    magma_tally3_int_t k1, magma_tally3_int_t k2,
    const magma_tally3_int_t *ipiv, magma_tally3_int_t inci,
    magma_tally3_queue_t queue );

void
magma_tally3blas_slaswpx_q(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldx, magma_tally3_int_t ldy,
    magma_tally3_int_t k1, magma_tally3_int_t k2,
    const magma_tally3_int_t *ipiv, magma_tally3_int_t inci,
    magma_tally3_queue_t queue );

void
magma_tally3blas_slaswp2_q(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dAT, magma_tally3_int_t ldda,
    magma_tally3_int_t k1, magma_tally3_int_t k2,
    magma_tally3Int_const_ptr d_ipiv, magma_tally3_int_t inci,
    magma_tally3_queue_t queue );

void
magma_tally3blas_ssymmetrize_q(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_queue_t queue );

void
magma_tally3blas_ssymmetrize_tiles_q(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t ntile, magma_tally3_int_t mstride, magma_tally3_int_t nstride,
    magma_tally3_queue_t queue );

void
magma_tally3blas_strtri_diag_q(
    magma_tally3_uplo_t uplo, magma_tally3_diag_t diag, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr d_dinvA,
    magma_tally3_queue_t queue );

  /*
   * Level 1 BLAS
   */
void
magma_tally3blas_sswap_q(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dx, magma_tally3_int_t incx,
    magma_tally3Float_ptr dy, magma_tally3_int_t incy,
    magma_tally3_queue_t queue );

void
magma_tally3blas_sswapblk_q(
    magma_tally3_order_t order,
    magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t i1, magma_tally3_int_t i2,
    const magma_tally3_int_t *ipiv, magma_tally3_int_t inci,
    magma_tally3_int_t offset,
    magma_tally3_queue_t queue );

void
magma_tally3blas_sswapdblk_q(
    magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t inca,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb, magma_tally3_int_t incb,
    magma_tally3_queue_t queue );

  /*
   * Level 2 BLAS
   */



  /*
   * Level 3 BLAS
   */



#ifdef __cplusplus
}
#endif

#undef REAL

#endif  /* MAGMA_tally3BLAS_S_H */
