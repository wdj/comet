/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally3blas_zc_q.h mixed zc -> ds, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_tally3BLAS_DS_Q_H
#define MAGMA_tally3BLAS_DS_Q_H

#include "magma_tally3_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Mixed precision */
void magma_tally3blas_dsaxpycp_q(
    magma_tally3_int_t m,
    magma_tally3Float_ptr  r,
    magma_tally3Double_ptr x,
    magma_tally3Double_const_ptr b,
    magma_tally3Double_ptr w,
    magma_tally3_queue_t queue );

void magma_tally3blas_daxpycp_q(
    magma_tally3_int_t m,
    magma_tally3Double_ptr r,
    magma_tally3Double_ptr x,
    magma_tally3Double_const_ptr b,
    magma_tally3_queue_t queue  );

void magma_tally3blas_dslaswp_q(
    magma_tally3_int_t n,
    magma_tally3Double_ptr A, magma_tally3_int_t lda,
    magma_tally3Float_ptr SA, magma_tally3_int_t m,
    const magma_tally3_int_t *ipiv, magma_tally3_int_t incx,
    magma_tally3_queue_t queue );

void magma_tally3blas_dlag2s_q(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_const_ptr A,  magma_tally3_int_t lda,
    magma_tally3Float_ptr       SA, magma_tally3_int_t ldsa,
    magma_tally3_int_t *info,
    magma_tally3_queue_t queue );

void magma_tally3blas_slag2d_q(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr SA, magma_tally3_int_t ldsa,
    magma_tally3Double_ptr       A,  magma_tally3_int_t lda,
    magma_tally3_int_t *info,
    magma_tally3_queue_t queue );

void magma_tally3blas_dlat2s_q(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Double_const_ptr A,  magma_tally3_int_t lda,
    magma_tally3Float_ptr       SA, magma_tally3_int_t ldsa,
    magma_tally3_int_t *info,
    magma_tally3_queue_t queue );

void magma_tally3blas_slat2d_q(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Float_const_ptr SA, magma_tally3_int_t ldsa,
    magma_tally3Double_ptr       A,  magma_tally3_int_t lda,
    magma_tally3_int_t *info,
    magma_tally3_queue_t queue );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally3BLAS_DS_H */
