/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally2blas_z_q.h normal z -> d, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_tally2BLAS_D_Q_H
#define MAGMA_tally2BLAS_D_Q_H
                    
#include "magma_tally2_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_tally2blas_dtranspose_inplace_q(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue );

void
magma_tally2blas_dtranspose_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2Double_ptr       dAT, magma_tally2_int_t lddat,
    magma_tally2_queue_t queue );

void
magma_tally2blas_dgetmatrix_transpose_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dAT,   magma_tally2_int_t ldda,
    double          *hA,    magma_tally2_int_t lda,
    magma_tally2Double_ptr       dwork, magma_tally2_int_t lddwork, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[2] );

void
magma_tally2blas_dsetmatrix_transpose_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const double *hA,    magma_tally2_int_t lda,
    magma_tally2Double_ptr    dAT,   magma_tally2_int_t ldda,
    magma_tally2Double_ptr    dwork, magma_tally2_int_t lddwork, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[2] );

  /*
   * Multi-GPU functions
   */



  /*
   * LAPACK auxiliary functions
   */
void
magma_tally2blas_dgeadd_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr       dB, magma_tally2_int_t lddb,
    magma_tally2_queue_t queue );

void
magma_tally2blas_dlacpy_q(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr       dB, magma_tally2_int_t lddb,
    magma_tally2_queue_t queue );

double
magma_tally2blas_dlange_q(
    magma_tally2_norm_t norm,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dwork,
    magma_tally2_queue_t queue );

double
magma_tally2blas_dlansy_q(
    magma_tally2_norm_t norm, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dwork,
    magma_tally2_queue_t queue );

double
magma_tally2blas_dlansy_q(
    magma_tally2_norm_t norm, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dwork,
    magma_tally2_queue_t queue );

void
magma_tally2blas_dlarfg_q(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dalpha,
    magma_tally2Double_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_ptr dtau,
    magma_tally2_queue_t queue );

void
magma_tally2blas_dlascl_q(
    magma_tally2_type_t type, magma_tally2_int_t kl, magma_tally2_int_t ku,
    double cfrom, double cto,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue,
    magma_tally2_int_t *info );

void
magma_tally2blas_dlascl2_q(
    magma_tally2_type_t type,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dD,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue,
    magma_tally2_int_t *info );

void
magma_tally2blas_dlaset_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    double offdiag, double diag,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue );

void
magma_tally2blas_dlaset_band_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    double offdiag, double diag,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue);

void
magma_tally2blas_dlaswp_q(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dAT, magma_tally2_int_t ldda,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci,
    magma_tally2_queue_t queue );

void
magma_tally2blas_dlaswpx_q(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldx, magma_tally2_int_t ldy,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci,
    magma_tally2_queue_t queue );

void
magma_tally2blas_dlaswp2_q(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dAT, magma_tally2_int_t ldda,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    magma_tally2Int_const_ptr d_ipiv, magma_tally2_int_t inci,
    magma_tally2_queue_t queue );

void
magma_tally2blas_dsymmetrize_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue );

void
magma_tally2blas_dsymmetrize_tiles_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t ntile, magma_tally2_int_t mstride, magma_tally2_int_t nstride,
    magma_tally2_queue_t queue );

void
magma_tally2blas_dtrtri_diag_q(
    magma_tally2_uplo_t uplo, magma_tally2_diag_t diag, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr d_dinvA,
    magma_tally2_queue_t queue );

  /*
   * Level 1 BLAS
   */
void
magma_tally2blas_dswap_q(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dx, magma_tally2_int_t incx,
    magma_tally2Double_ptr dy, magma_tally2_int_t incy,
    magma_tally2_queue_t queue );

void
magma_tally2blas_dswapblk_q(
    magma_tally2_order_t order,
    magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t i1, magma_tally2_int_t i2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci,
    magma_tally2_int_t offset,
    magma_tally2_queue_t queue );

void
magma_tally2blas_dswapdblk_q(
    magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t inca,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb, magma_tally2_int_t incb,
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

#undef REAL

#endif  /* MAGMA_tally2BLAS_D_H */
