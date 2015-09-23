/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproductblas_z_q.h normal z -> d, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_minproductBLAS_D_Q_H
#define MAGMA_minproductBLAS_D_Q_H
                    
#include "magma_minproduct_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_minproductblas_dtranspose_inplace_q(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_dtranspose_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dAT, magma_minproduct_int_t lddat,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_dgetmatrix_transpose_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dAT,   magma_minproduct_int_t ldda,
    double          *hA,    magma_minproduct_int_t lda,
    magma_minproductDouble_ptr       dwork, magma_minproduct_int_t lddwork, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[2] );

void
magma_minproductblas_dsetmatrix_transpose_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const double *hA,    magma_minproduct_int_t lda,
    magma_minproductDouble_ptr    dAT,   magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr    dwork, magma_minproduct_int_t lddwork, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[2] );

  /*
   * Multi-GPU functions
   */



  /*
   * LAPACK auxiliary functions
   */
void
magma_minproductblas_dgeadd_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dB, magma_minproduct_int_t lddb,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_dlacpy_q(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dB, magma_minproduct_int_t lddb,
    magma_minproduct_queue_t queue );

double
magma_minproductblas_dlange_q(
    magma_minproduct_norm_t norm,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_queue_t queue );

double
magma_minproductblas_dlansy_q(
    magma_minproduct_norm_t norm, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_queue_t queue );

double
magma_minproductblas_dlansy_q(
    magma_minproduct_norm_t norm, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_dlarfg_q(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dalpha,
    magma_minproductDouble_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr dtau,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_dlascl_q(
    magma_minproduct_type_t type, magma_minproduct_int_t kl, magma_minproduct_int_t ku,
    double cfrom, double cto,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue,
    magma_minproduct_int_t *info );

void
magma_minproductblas_dlascl2_q(
    magma_minproduct_type_t type,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dD,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue,
    magma_minproduct_int_t *info );

void
magma_minproductblas_dlaset_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double offdiag, double diag,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_dlaset_band_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double offdiag, double diag,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue);

void
magma_minproductblas_dlaswp_q(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dAT, magma_minproduct_int_t ldda,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_dlaswpx_q(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldx, magma_minproduct_int_t ldy,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_dlaswp2_q(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dAT, magma_minproduct_int_t ldda,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    magma_minproductInt_const_ptr d_ipiv, magma_minproduct_int_t inci,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_dsymmetrize_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_dsymmetrize_tiles_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t ntile, magma_minproduct_int_t mstride, magma_minproduct_int_t nstride,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_dtrtri_diag_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr d_dinvA,
    magma_minproduct_queue_t queue );

  /*
   * Level 1 BLAS
   */
void
magma_minproductblas_dswap_q(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr dy, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_dswapblk_q(
    magma_minproduct_order_t order,
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t i1, magma_minproduct_int_t i2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci,
    magma_minproduct_int_t offset,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_dswapdblk_q(
    magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t inca,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb, magma_minproduct_int_t incb,
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

#undef REAL

#endif  /* MAGMA_minproductBLAS_D_H */
