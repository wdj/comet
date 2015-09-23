/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
*/

#ifndef MAGMA_minproductBLAS_Z_Q_H
#define MAGMA_minproductBLAS_Z_Q_H
                    
#include "magma_minproduct_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_minproductblas_ztranspose_inplace_q(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_ztranspose_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr       dAT, magma_minproduct_int_t lddat,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_zgetmatrix_transpose_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dAT,   magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex          *hA,    magma_minproduct_int_t lda,
    magma_minproductDoubleComplex_ptr       dwork, magma_minproduct_int_t lddwork, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[2] );

void
magma_minproductblas_zsetmatrix_transpose_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const magma_minproductDoubleComplex *hA,    magma_minproduct_int_t lda,
    magma_minproductDoubleComplex_ptr    dAT,   magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr    dwork, magma_minproduct_int_t lddwork, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[2] );

  /*
   * Multi-GPU functions
   */



  /*
   * LAPACK auxiliary functions
   */
void
magma_minproductblas_zgeadd_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr       dB, magma_minproduct_int_t lddb,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_zlacpy_q(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr       dB, magma_minproduct_int_t lddb,
    magma_minproduct_queue_t queue );

double
magma_minproductblas_zlange_q(
    magma_minproduct_norm_t norm,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_queue_t queue );

double
magma_minproductblas_zlanhe_q(
    magma_minproduct_norm_t norm, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_queue_t queue );

double
magma_minproductblas_zlansy_q(
    magma_minproduct_norm_t norm, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_zlarfg_q(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dalpha,
    magma_minproductDoubleComplex_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_ptr dtau,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_zlascl_q(
    magma_minproduct_type_t type, magma_minproduct_int_t kl, magma_minproduct_int_t ku,
    double cfrom, double cto,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue,
    magma_minproduct_int_t *info );

void
magma_minproductblas_zlascl2_q(
    magma_minproduct_type_t type,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dD,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue,
    magma_minproduct_int_t *info );

void
magma_minproductblas_zlaset_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex offdiag, magma_minproductDoubleComplex diag,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_zlaset_band_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex offdiag, magma_minproductDoubleComplex diag,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue);

void
magma_minproductblas_zlaswp_q(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dAT, magma_minproduct_int_t ldda,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_zlaswpx_q(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldx, magma_minproduct_int_t ldy,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_zlaswp2_q(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dAT, magma_minproduct_int_t ldda,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    magma_minproductInt_const_ptr d_ipiv, magma_minproduct_int_t inci,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_zsymmetrize_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_zsymmetrize_tiles_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t ntile, magma_minproduct_int_t mstride, magma_minproduct_int_t nstride,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_ztrtri_diag_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr d_dinvA,
    magma_minproduct_queue_t queue );

  /*
   * Level 1 BLAS
   */
void
magma_minproductblas_zswap_q(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex_ptr dy, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_zswapblk_q(
    magma_minproduct_order_t order,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t i1, magma_minproduct_int_t i2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci,
    magma_minproduct_int_t offset,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_zswapdblk_q(
    magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t inca,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb, magma_minproduct_int_t incb,
    magma_minproduct_queue_t queue );

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

#endif  /* MAGMA_minproductBLAS_Z_H */
