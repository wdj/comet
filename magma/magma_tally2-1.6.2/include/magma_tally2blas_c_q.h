/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally2blas_z_q.h normal z -> c, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_tally2BLAS_C_Q_H
#define MAGMA_tally2BLAS_C_Q_H
                    
#include "magma_tally2_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_tally2blas_ctranspose_inplace_q(
    magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue );

void
magma_tally2blas_ctranspose_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr       dAT, magma_tally2_int_t lddat,
    magma_tally2_queue_t queue );

void
magma_tally2blas_cgetmatrix_transpose_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr dAT,   magma_tally2_int_t ldda,
    magma_tally2FloatComplex          *hA,    magma_tally2_int_t lda,
    magma_tally2FloatComplex_ptr       dwork, magma_tally2_int_t lddwork, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[2] );

void
magma_tally2blas_csetmatrix_transpose_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const magma_tally2FloatComplex *hA,    magma_tally2_int_t lda,
    magma_tally2FloatComplex_ptr    dAT,   magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr    dwork, magma_tally2_int_t lddwork, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[2] );

  /*
   * Multi-GPU functions
   */



  /*
   * LAPACK auxiliary functions
   */
void
magma_tally2blas_cgeadd_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr       dB, magma_tally2_int_t lddb,
    magma_tally2_queue_t queue );

void
magma_tally2blas_clacpy_q(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr       dB, magma_tally2_int_t lddb,
    magma_tally2_queue_t queue );

float
magma_tally2blas_clange_q(
    magma_tally2_norm_t norm,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dwork,
    magma_tally2_queue_t queue );

float
magma_tally2blas_clanhe_q(
    magma_tally2_norm_t norm, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dwork,
    magma_tally2_queue_t queue );

float
magma_tally2blas_clansy_q(
    magma_tally2_norm_t norm, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dwork,
    magma_tally2_queue_t queue );

void
magma_tally2blas_clarfg_q(
    magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dalpha,
    magma_tally2FloatComplex_ptr dx, magma_tally2_int_t incx,
    magma_tally2FloatComplex_ptr dtau,
    magma_tally2_queue_t queue );

void
magma_tally2blas_clascl_q(
    magma_tally2_type_t type, magma_tally2_int_t kl, magma_tally2_int_t ku,
    float cfrom, float cto,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue,
    magma_tally2_int_t *info );

void
magma_tally2blas_clascl2_q(
    magma_tally2_type_t type,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_const_ptr dD,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue,
    magma_tally2_int_t *info );

void
magma_tally2blas_claset_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex offdiag, magma_tally2FloatComplex diag,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue );

void
magma_tally2blas_claset_band_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex offdiag, magma_tally2FloatComplex diag,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue);

void
magma_tally2blas_claswp_q(
    magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dAT, magma_tally2_int_t ldda,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci,
    magma_tally2_queue_t queue );

void
magma_tally2blas_claswpx_q(
    magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldx, magma_tally2_int_t ldy,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci,
    magma_tally2_queue_t queue );

void
magma_tally2blas_claswp2_q(
    magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dAT, magma_tally2_int_t ldda,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    magma_tally2Int_const_ptr d_ipiv, magma_tally2_int_t inci,
    magma_tally2_queue_t queue );

void
magma_tally2blas_csymmetrize_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue );

void
magma_tally2blas_csymmetrize_tiles_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t ntile, magma_tally2_int_t mstride, magma_tally2_int_t nstride,
    magma_tally2_queue_t queue );

void
magma_tally2blas_ctrtri_diag_q(
    magma_tally2_uplo_t uplo, magma_tally2_diag_t diag, magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr d_dinvA,
    magma_tally2_queue_t queue );

  /*
   * Level 1 BLAS
   */
void
magma_tally2blas_cswap_q(
    magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dx, magma_tally2_int_t incx,
    magma_tally2FloatComplex_ptr dy, magma_tally2_int_t incy,
    magma_tally2_queue_t queue );

void
magma_tally2blas_cswapblk_q(
    magma_tally2_order_t order,
    magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t i1, magma_tally2_int_t i2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci,
    magma_tally2_int_t offset,
    magma_tally2_queue_t queue );

void
magma_tally2blas_cswapdblk_q(
    magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t inca,
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb, magma_tally2_int_t incb,
    magma_tally2_queue_t queue );

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

#endif  /* MAGMA_tally2BLAS_C_H */
