/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally4blas_z_q.h normal z -> c, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_tally4BLAS_C_Q_H
#define MAGMA_tally4BLAS_C_Q_H
                    
#include "magma_tally4_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_tally4blas_ctranspose_inplace_q(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_queue_t queue );

void
magma_tally4blas_ctranspose_q(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr       dAT, magma_tally4_int_t lddat,
    magma_tally4_queue_t queue );

void
magma_tally4blas_cgetmatrix_transpose_q(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dAT,   magma_tally4_int_t ldda,
    magma_tally4FloatComplex          *hA,    magma_tally4_int_t lda,
    magma_tally4FloatComplex_ptr       dwork, magma_tally4_int_t lddwork, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[2] );

void
magma_tally4blas_csetmatrix_transpose_q(
    magma_tally4_int_t m, magma_tally4_int_t n,
    const magma_tally4FloatComplex *hA,    magma_tally4_int_t lda,
    magma_tally4FloatComplex_ptr    dAT,   magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr    dwork, magma_tally4_int_t lddwork, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[2] );

  /*
   * Multi-GPU functions
   */



  /*
   * LAPACK auxiliary functions
   */
void
magma_tally4blas_cgeadd_q(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr       dB, magma_tally4_int_t lddb,
    magma_tally4_queue_t queue );

void
magma_tally4blas_clacpy_q(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr       dB, magma_tally4_int_t lddb,
    magma_tally4_queue_t queue );

float
magma_tally4blas_clange_q(
    magma_tally4_norm_t norm,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dwork,
    magma_tally4_queue_t queue );

float
magma_tally4blas_clanhe_q(
    magma_tally4_norm_t norm, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dwork,
    magma_tally4_queue_t queue );

float
magma_tally4blas_clansy_q(
    magma_tally4_norm_t norm, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dwork,
    magma_tally4_queue_t queue );

void
magma_tally4blas_clarfg_q(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dalpha,
    magma_tally4FloatComplex_ptr dx, magma_tally4_int_t incx,
    magma_tally4FloatComplex_ptr dtau,
    magma_tally4_queue_t queue );

void
magma_tally4blas_clascl_q(
    magma_tally4_type_t type, magma_tally4_int_t kl, magma_tally4_int_t ku,
    float cfrom, float cto,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_queue_t queue,
    magma_tally4_int_t *info );

void
magma_tally4blas_clascl2_q(
    magma_tally4_type_t type,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr dD,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_queue_t queue,
    magma_tally4_int_t *info );

void
magma_tally4blas_claset_q(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex offdiag, magma_tally4FloatComplex diag,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_queue_t queue );

void
magma_tally4blas_claset_band_q(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex offdiag, magma_tally4FloatComplex diag,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_queue_t queue);

void
magma_tally4blas_claswp_q(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dAT, magma_tally4_int_t ldda,
    magma_tally4_int_t k1, magma_tally4_int_t k2,
    const magma_tally4_int_t *ipiv, magma_tally4_int_t inci,
    magma_tally4_queue_t queue );

void
magma_tally4blas_claswpx_q(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldx, magma_tally4_int_t ldy,
    magma_tally4_int_t k1, magma_tally4_int_t k2,
    const magma_tally4_int_t *ipiv, magma_tally4_int_t inci,
    magma_tally4_queue_t queue );

void
magma_tally4blas_claswp2_q(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dAT, magma_tally4_int_t ldda,
    magma_tally4_int_t k1, magma_tally4_int_t k2,
    magma_tally4Int_const_ptr d_ipiv, magma_tally4_int_t inci,
    magma_tally4_queue_t queue );

void
magma_tally4blas_csymmetrize_q(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_queue_t queue );

void
magma_tally4blas_csymmetrize_tiles_q(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t ntile, magma_tally4_int_t mstride, magma_tally4_int_t nstride,
    magma_tally4_queue_t queue );

void
magma_tally4blas_ctrtri_diag_q(
    magma_tally4_uplo_t uplo, magma_tally4_diag_t diag, magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr d_dinvA,
    magma_tally4_queue_t queue );

  /*
   * Level 1 BLAS
   */
void
magma_tally4blas_cswap_q(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dx, magma_tally4_int_t incx,
    magma_tally4FloatComplex_ptr dy, magma_tally4_int_t incy,
    magma_tally4_queue_t queue );

void
magma_tally4blas_cswapblk_q(
    magma_tally4_order_t order,
    magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t i1, magma_tally4_int_t i2,
    const magma_tally4_int_t *ipiv, magma_tally4_int_t inci,
    magma_tally4_int_t offset,
    magma_tally4_queue_t queue );

void
magma_tally4blas_cswapdblk_q(
    magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t inca,
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb, magma_tally4_int_t incb,
    magma_tally4_queue_t queue );

  /*
   * Level 2 BLAS
   */



  /*
   * Level 3 BLAS
   */



#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif  /* MAGMA_tally4BLAS_C_H */
