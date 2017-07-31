/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
*/

#ifndef MAGMA_tally2BLAS_Z_Q_H
#define MAGMA_tally2BLAS_Z_Q_H
                    
#include "magma_tally2_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_tally2blas_ztranspose_inplace_q(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue );

void
magma_tally2blas_ztranspose_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr       dAT, magma_tally2_int_t lddat,
    magma_tally2_queue_t queue );

void
magma_tally2blas_zgetmatrix_transpose_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dAT,   magma_tally2_int_t ldda,
    magma_tally2DoubleComplex          *hA,    magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr       dwork, magma_tally2_int_t lddwork, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[2] );

void
magma_tally2blas_zsetmatrix_transpose_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const magma_tally2DoubleComplex *hA,    magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr    dAT,   magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr    dwork, magma_tally2_int_t lddwork, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[2] );

  /*
   * Multi-GPU functions
   */



  /*
   * LAPACK auxiliary functions
   */
void
magma_tally2blas_zgeadd_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr       dB, magma_tally2_int_t lddb,
    magma_tally2_queue_t queue );

void
magma_tally2blas_zlacpy_q(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr       dB, magma_tally2_int_t lddb,
    magma_tally2_queue_t queue );

double
magma_tally2blas_zlange_q(
    magma_tally2_norm_t norm,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dwork,
    magma_tally2_queue_t queue );

double
magma_tally2blas_zlanhe_q(
    magma_tally2_norm_t norm, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dwork,
    magma_tally2_queue_t queue );

double
magma_tally2blas_zlansy_q(
    magma_tally2_norm_t norm, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dwork,
    magma_tally2_queue_t queue );

void
magma_tally2blas_zlarfg_q(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dalpha,
    magma_tally2DoubleComplex_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr dtau,
    magma_tally2_queue_t queue );

void
magma_tally2blas_zlascl_q(
    magma_tally2_type_t type, magma_tally2_int_t kl, magma_tally2_int_t ku,
    double cfrom, double cto,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue,
    magma_tally2_int_t *info );

void
magma_tally2blas_zlascl2_q(
    magma_tally2_type_t type,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dD,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue,
    magma_tally2_int_t *info );

void
magma_tally2blas_zlaset_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex offdiag, magma_tally2DoubleComplex diag,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue );

void
magma_tally2blas_zlaset_band_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex offdiag, magma_tally2DoubleComplex diag,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue);

void
magma_tally2blas_zlaswp_q(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dAT, magma_tally2_int_t ldda,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci,
    magma_tally2_queue_t queue );

void
magma_tally2blas_zlaswpx_q(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldx, magma_tally2_int_t ldy,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci,
    magma_tally2_queue_t queue );

void
magma_tally2blas_zlaswp2_q(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dAT, magma_tally2_int_t ldda,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    magma_tally2Int_const_ptr d_ipiv, magma_tally2_int_t inci,
    magma_tally2_queue_t queue );

void
magma_tally2blas_zsymmetrize_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue );

void
magma_tally2blas_zsymmetrize_tiles_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t ntile, magma_tally2_int_t mstride, magma_tally2_int_t nstride,
    magma_tally2_queue_t queue );

void
magma_tally2blas_ztrtri_diag_q(
    magma_tally2_uplo_t uplo, magma_tally2_diag_t diag, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr d_dinvA,
    magma_tally2_queue_t queue );

  /*
   * Level 1 BLAS
   */
void
magma_tally2blas_zswap_q(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dx, magma_tally2_int_t incx,
    magma_tally2DoubleComplex_ptr dy, magma_tally2_int_t incy,
    magma_tally2_queue_t queue );

void
magma_tally2blas_zswapblk_q(
    magma_tally2_order_t order,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t i1, magma_tally2_int_t i2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci,
    magma_tally2_int_t offset,
    magma_tally2_queue_t queue );

void
magma_tally2blas_zswapdblk_q(
    magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t inca,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb, magma_tally2_int_t incb,
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

#endif  /* MAGMA_tally2BLAS_Z_H */
