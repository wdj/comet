/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproductblas_z_q.h normal z -> s, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_minproductBLAS_S_Q_H
#define MAGMA_minproductBLAS_S_Q_H
                    
#include "magma_minproduct_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_minproductblas_stranspose_inplace_q(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_stranspose_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr       dAT, magma_minproduct_int_t lddat,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_sgetmatrix_transpose_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dAT,   magma_minproduct_int_t ldda,
    float          *hA,    magma_minproduct_int_t lda,
    magma_minproductFloat_ptr       dwork, magma_minproduct_int_t lddwork, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[2] );

void
magma_minproductblas_ssetmatrix_transpose_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const float *hA,    magma_minproduct_int_t lda,
    magma_minproductFloat_ptr    dAT,   magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr    dwork, magma_minproduct_int_t lddwork, magma_minproduct_int_t nb,
    magma_minproduct_queue_t queues[2] );

  /*
   * Multi-GPU functions
   */



  /*
   * LAPACK auxiliary functions
   */
void
magma_minproductblas_sgeadd_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr       dB, magma_minproduct_int_t lddb,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_slacpy_q(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr       dB, magma_minproduct_int_t lddb,
    magma_minproduct_queue_t queue );

float
magma_minproductblas_slange_q(
    magma_minproduct_norm_t norm,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_queue_t queue );

float
magma_minproductblas_slansy_q(
    magma_minproduct_norm_t norm, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_queue_t queue );

float
magma_minproductblas_slansy_q(
    magma_minproduct_norm_t norm, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_slarfg_q(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dalpha,
    magma_minproductFloat_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr dtau,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_slascl_q(
    magma_minproduct_type_t type, magma_minproduct_int_t kl, magma_minproduct_int_t ku,
    float cfrom, float cto,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue,
    magma_minproduct_int_t *info );

void
magma_minproductblas_slascl2_q(
    magma_minproduct_type_t type,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dD,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue,
    magma_minproduct_int_t *info );

void
magma_minproductblas_slaset_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    float offdiag, float diag,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_slaset_band_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float offdiag, float diag,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue);

void
magma_minproductblas_slaswp_q(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dAT, magma_minproduct_int_t ldda,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_slaswpx_q(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldx, magma_minproduct_int_t ldy,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_slaswp2_q(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dAT, magma_minproduct_int_t ldda,
    magma_minproduct_int_t k1, magma_minproduct_int_t k2,
    magma_minproductInt_const_ptr d_ipiv, magma_minproduct_int_t inci,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_ssymmetrize_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_ssymmetrize_tiles_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t ntile, magma_minproduct_int_t mstride, magma_minproduct_int_t nstride,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_strtri_diag_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr d_dinvA,
    magma_minproduct_queue_t queue );

  /*
   * Level 1 BLAS
   */
void
magma_minproductblas_sswap_q(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dx, magma_minproduct_int_t incx,
    magma_minproductFloat_ptr dy, magma_minproduct_int_t incy,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_sswapblk_q(
    magma_minproduct_order_t order,
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t i1, magma_minproduct_int_t i2,
    const magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci,
    magma_minproduct_int_t offset,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_sswapdblk_q(
    magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t inca,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb, magma_minproduct_int_t incb,
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

#endif  /* MAGMA_minproductBLAS_S_H */
