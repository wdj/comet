/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally4blas_z_q.h normal z -> d, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_tally4BLAS_D_Q_H
#define MAGMA_tally4BLAS_D_Q_H
                    
#include "magma_tally4_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Transpose functions
   */
void
magma_tally4blas_dtranspose_inplace_q(
    magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_queue_t queue );

void
magma_tally4blas_dtranspose_q(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_const_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4Double_ptr       dAT, magma_tally4_int_t lddat,
    magma_tally4_queue_t queue );

void
magma_tally4blas_dgetmatrix_transpose_q(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_const_ptr dAT,   magma_tally4_int_t ldda,
    double          *hA,    magma_tally4_int_t lda,
    magma_tally4Double_ptr       dwork, magma_tally4_int_t lddwork, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[2] );

void
magma_tally4blas_dsetmatrix_transpose_q(
    magma_tally4_int_t m, magma_tally4_int_t n,
    const double *hA,    magma_tally4_int_t lda,
    magma_tally4Double_ptr    dAT,   magma_tally4_int_t ldda,
    magma_tally4Double_ptr    dwork, magma_tally4_int_t lddwork, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[2] );

  /*
   * Multi-GPU functions
   */



  /*
   * LAPACK auxiliary functions
   */
void
magma_tally4blas_dgeadd_q(
    magma_tally4_int_t m, magma_tally4_int_t n,
    double alpha,
    magma_tally4Double_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr       dB, magma_tally4_int_t lddb,
    magma_tally4_queue_t queue );

void
magma_tally4blas_dlacpy_q(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr       dB, magma_tally4_int_t lddb,
    magma_tally4_queue_t queue );

double
magma_tally4blas_dlange_q(
    magma_tally4_norm_t norm,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dwork,
    magma_tally4_queue_t queue );

double
magma_tally4blas_dlansy_q(
    magma_tally4_norm_t norm, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4Double_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dwork,
    magma_tally4_queue_t queue );

double
magma_tally4blas_dlansy_q(
    magma_tally4_norm_t norm, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4Double_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dwork,
    magma_tally4_queue_t queue );

void
magma_tally4blas_dlarfg_q(
    magma_tally4_int_t n,
    magma_tally4Double_ptr dalpha,
    magma_tally4Double_ptr dx, magma_tally4_int_t incx,
    magma_tally4Double_ptr dtau,
    magma_tally4_queue_t queue );

void
magma_tally4blas_dlascl_q(
    magma_tally4_type_t type, magma_tally4_int_t kl, magma_tally4_int_t ku,
    double cfrom, double cto,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_queue_t queue,
    magma_tally4_int_t *info );

void
magma_tally4blas_dlascl2_q(
    magma_tally4_type_t type,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_const_ptr dD,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_queue_t queue,
    magma_tally4_int_t *info );

void
magma_tally4blas_dlaset_q(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    double offdiag, double diag,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_queue_t queue );

void
magma_tally4blas_dlaset_band_q(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    double offdiag, double diag,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_queue_t queue);

void
magma_tally4blas_dlaswp_q(
    magma_tally4_int_t n,
    magma_tally4Double_ptr dAT, magma_tally4_int_t ldda,
    magma_tally4_int_t k1, magma_tally4_int_t k2,
    const magma_tally4_int_t *ipiv, magma_tally4_int_t inci,
    magma_tally4_queue_t queue );

void
magma_tally4blas_dlaswpx_q(
    magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldx, magma_tally4_int_t ldy,
    magma_tally4_int_t k1, magma_tally4_int_t k2,
    const magma_tally4_int_t *ipiv, magma_tally4_int_t inci,
    magma_tally4_queue_t queue );

void
magma_tally4blas_dlaswp2_q(
    magma_tally4_int_t n,
    magma_tally4Double_ptr dAT, magma_tally4_int_t ldda,
    magma_tally4_int_t k1, magma_tally4_int_t k2,
    magma_tally4Int_const_ptr d_ipiv, magma_tally4_int_t inci,
    magma_tally4_queue_t queue );

void
magma_tally4blas_dsymmetrize_q(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_queue_t queue );

void
magma_tally4blas_dsymmetrize_tiles_q(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t ntile, magma_tally4_int_t mstride, magma_tally4_int_t nstride,
    magma_tally4_queue_t queue );

void
magma_tally4blas_dtrtri_diag_q(
    magma_tally4_uplo_t uplo, magma_tally4_diag_t diag, magma_tally4_int_t n,
    magma_tally4Double_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr d_dinvA,
    magma_tally4_queue_t queue );

  /*
   * Level 1 BLAS
   */
void
magma_tally4blas_dswap_q(
    magma_tally4_int_t n,
    magma_tally4Double_ptr dx, magma_tally4_int_t incx,
    magma_tally4Double_ptr dy, magma_tally4_int_t incy,
    magma_tally4_queue_t queue );

void
magma_tally4blas_dswapblk_q(
    magma_tally4_order_t order,
    magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t i1, magma_tally4_int_t i2,
    const magma_tally4_int_t *ipiv, magma_tally4_int_t inci,
    magma_tally4_int_t offset,
    magma_tally4_queue_t queue );

void
magma_tally4blas_dswapdblk_q(
    magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t inca,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb, magma_tally4_int_t incb,
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

#undef REAL

#endif  /* MAGMA_tally4BLAS_D_H */
