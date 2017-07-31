/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally2blas_zc_q.h mixed zc -> ds, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_tally2BLAS_DS_Q_H
#define MAGMA_tally2BLAS_DS_Q_H

#include "magma_tally2_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Mixed precision */
void magma_tally2blas_dsaxpycp_q(
    magma_tally2_int_t m,
    magma_tally2Float_ptr  r,
    magma_tally2Double_ptr x,
    magma_tally2Double_const_ptr b,
    magma_tally2Double_ptr w,
    magma_tally2_queue_t queue );

void magma_tally2blas_daxpycp_q(
    magma_tally2_int_t m,
    magma_tally2Double_ptr r,
    magma_tally2Double_ptr x,
    magma_tally2Double_const_ptr b,
    magma_tally2_queue_t queue  );

void magma_tally2blas_dslaswp_q(
    magma_tally2_int_t n,
    magma_tally2Double_ptr A, magma_tally2_int_t lda,
    magma_tally2Float_ptr SA, magma_tally2_int_t m,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t incx,
    magma_tally2_queue_t queue );

void magma_tally2blas_dlag2s_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr A,  magma_tally2_int_t lda,
    magma_tally2Float_ptr       SA, magma_tally2_int_t ldsa,
    magma_tally2_int_t *info,
    magma_tally2_queue_t queue );

void magma_tally2blas_slag2d_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_const_ptr SA, magma_tally2_int_t ldsa,
    magma_tally2Double_ptr       A,  magma_tally2_int_t lda,
    magma_tally2_int_t *info,
    magma_tally2_queue_t queue );

void magma_tally2blas_dlat2s_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Double_const_ptr A,  magma_tally2_int_t lda,
    magma_tally2Float_ptr       SA, magma_tally2_int_t ldsa,
    magma_tally2_int_t *info,
    magma_tally2_queue_t queue );

void magma_tally2blas_slat2d_q(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Float_const_ptr SA, magma_tally2_int_t ldsa,
    magma_tally2Double_ptr       A,  magma_tally2_int_t lda,
    magma_tally2_int_t *info,
    magma_tally2_queue_t queue );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally2BLAS_DS_H */
