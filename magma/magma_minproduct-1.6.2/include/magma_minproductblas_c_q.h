/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproductblas_z_q.h normal z -> c, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_minproductBLAS_C_Q_H
#define MAGMA_minproductBLAS_C_Q_H
                    
#include "magma_minproduct_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_minproductblas_ctranspose_inplace_q(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_ctranspose_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr       dAT, magma_minproduct_int_t lddat,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_cgetmatrix_transpose_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dAT,   magma_minproduct_int_t ldda,
    magma_minproductFloatComplex          *hA,    magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr       dwork, magma_minproduct_int_t lddwork, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[2] );

void
magma_minproductblas_csetmatrix_transpose_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const magma_minproductFloatComplex *hA,    magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr    dAT,   magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr    dwork, magma_minproduct_int_t lddwork, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[2] );

  /*
   * Multi-GPU functions
   */



  /*
   * LAPACK auxiliary functions
   */
void
magma_minproductblas_cgeadd_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr       dB, magma_minproduct_int_t lddb,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_clacpy_q(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr       dB, magma_minproduct_int_t lddb,
    magma_minproduct_queue_t queue );

float
magma_minproductblas_clange_q(
    magma_minproduct_norm_t norm,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_queue_t queue );

float
magma_minproductblas_clanhe_q(
    magma_minproduct_norm_t norm, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_queue_t queue );

float
magma_minproductblas_clansy_q(
    magma_minproduct_norm_t norm, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_clarfg_q(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dalpha,
    magma_minproductFloatComplex_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr dtau,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_clascl_q(
    magma_minproduct_type_t type, magma_minproduct_int_t kl, magma_minproduct_int_t ku,
    float cfrom, float cto,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue,
    magma_minproduct_int_t *info );

void
magma_minproductblas_clascl2_q(
    magma_minproduct_type_t type,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dD,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue,
    magma_minproduct_int_t *info );

void
magma_minproductblas_claset_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex offdiag, magma_minproductFloatComplex diag,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_claset_band_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex offdiag, magma_minproductFloatComplex diag,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue);

void
magma_minproductblas_claswp_q(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dAT, magma_minproduct_int_t ldda,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_claswpx_q(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldx, magma_minproduct_int_t ldy,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_claswp2_q(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dAT, magma_minproduct_int_t ldda,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    magma_minproductInt_const_ptr d_ipiv, magma_minproduct_int_t inci,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_csymmetrize_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_csymmetrize_tiles_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t ntile, magma_minproduct_int_t mstride, magma_minproduct_int_t nstride,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_ctrtri_diag_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr d_dinvA,
    magma_minproduct_queue_t queue );

  /*
   * Level 1 BLAS
   */
void
magma_minproductblas_cswap_q(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloatComplex_ptr dy, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_cswapblk_q(
    magma_minproduct_order_t order,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t i1, magma_minproduct_int_t i2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci,
    magma_minproduct_int_t offset,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_cswapdblk_q(
    magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t inca,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb, magma_minproduct_int_t incb,
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

#endif  /* MAGMA_minproductBLAS_C_H */
