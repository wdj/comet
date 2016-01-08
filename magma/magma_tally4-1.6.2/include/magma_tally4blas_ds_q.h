/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally4blas_zc_q.h mixed zc -> ds, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_tally4BLAS_DS_Q_H
#define MAGMA_tally4BLAS_DS_Q_H

#include "magma_tally4_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Mixed precision */
void magma_tally4blas_dsaxpycp_q(
    magma_tally4_int_t m,
    magma_tally4Float_ptr  r,
    magma_tally4Double_ptr x,
    magma_tally4Double_const_ptr b,
    magma_tally4Double_ptr w,
    magma_tally4_queue_t queue );

void magma_tally4blas_daxpycp_q(
    magma_tally4_int_t m,
    magma_tally4Double_ptr r,
    magma_tally4Double_ptr x,
    magma_tally4Double_const_ptr b,
    magma_tally4_queue_t queue  );

void magma_tally4blas_dslaswp_q(
    magma_tally4_int_t n,
    magma_tally4Double_ptr A, magma_tally4_int_t lda,
    magma_tally4Float_ptr SA, magma_tally4_int_t m,
    const magma_tally4_int_t *ipiv, magma_tally4_int_t incx,
    magma_tally4_queue_t queue );

void magma_tally4blas_dlag2s_q(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_const_ptr A,  magma_tally4_int_t lda,
    magma_tally4Float_ptr       SA, magma_tally4_int_t ldsa,
    magma_tally4_int_t *info,
    magma_tally4_queue_t queue );

void magma_tally4blas_slag2d_q(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr SA, magma_tally4_int_t ldsa,
    magma_tally4Double_ptr       A,  magma_tally4_int_t lda,
    magma_tally4_int_t *info,
    magma_tally4_queue_t queue );

void magma_tally4blas_dlat2s_q(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Double_const_ptr A,  magma_tally4_int_t lda,
    magma_tally4Float_ptr       SA, magma_tally4_int_t ldsa,
    magma_tally4_int_t *info,
    magma_tally4_queue_t queue );

void magma_tally4blas_slat2d_q(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Float_const_ptr SA, magma_tally4_int_t ldsa,
    magma_tally4Double_ptr       A,  magma_tally4_int_t lda,
    magma_tally4_int_t *info,
    magma_tally4_queue_t queue );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally4BLAS_DS_H */
