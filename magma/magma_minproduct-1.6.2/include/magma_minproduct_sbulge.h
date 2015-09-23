/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproduct_zbulge.h normal z -> s, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_minproduct_SBULGE_H
#define MAGMA_minproduct_SBULGE_H

#include "magma_minproduct_types.h"

#ifdef __cplusplus
extern "C" {
#endif

magma_minproduct_int_t
magma_minproduct_sbulge_applyQ_v2(
    magma_minproduct_side_t side, 
    magma_minproduct_int_t NE, magma_minproduct_int_t N, 
    magma_minproduct_int_t NB, magma_minproduct_int_t Vblksiz, 
    magma_minproductFloat_ptr dE, magma_minproduct_int_t ldde, 
    float *V, magma_minproduct_int_t ldv, 
    float *T, magma_minproduct_int_t ldt, 
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sbulge_applyQ_v2_m(
    magma_minproduct_int_t ngpu, magma_minproduct_side_t side, 
    magma_minproduct_int_t NE, magma_minproduct_int_t N, 
    magma_minproduct_int_t NB, magma_minproduct_int_t Vblksiz, 
    float *E, magma_minproduct_int_t lde, 
    float *V, magma_minproduct_int_t ldv, 
    float *T, magma_minproduct_int_t ldt, 
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sbulge_back(
    magma_minproduct_uplo_t uplo, 
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    magma_minproduct_int_t ne, magma_minproduct_int_t Vblksiz,
    float *Z, magma_minproduct_int_t ldz,
    magma_minproductFloat_ptr dZ, magma_minproduct_int_t lddz,
    float *V, magma_minproduct_int_t ldv,
    float *TAU,
    float *T, magma_minproduct_int_t ldt,
    magma_minproduct_int_t* info);

magma_minproduct_int_t
magma_minproduct_sbulge_back_m(
    magma_minproduct_int_t ngpu, magma_minproduct_uplo_t uplo, 
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    magma_minproduct_int_t ne, magma_minproduct_int_t Vblksiz,
    float *Z, magma_minproduct_int_t ldz,
    float *V, magma_minproduct_int_t ldv, 
    float *TAU, 
    float *T, magma_minproduct_int_t ldt, 
    magma_minproduct_int_t* info);

void
magma_minproduct_strdtype1cbHLsym_withQ_v2(
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    float *A, magma_minproduct_int_t lda, 
    float *V, magma_minproduct_int_t ldv, 
    float *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, 
    magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz, 
    float *work);

void
magma_minproduct_strdtype2cbHLsym_withQ_v2(
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    float *A, magma_minproduct_int_t lda, 
    float *V, magma_minproduct_int_t ldv, 
    float *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, 
    magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz, 
    float *work);

void
magma_minproduct_strdtype3cbHLsym_withQ_v2(
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    float *A, magma_minproduct_int_t lda, 
    float *V, magma_minproduct_int_t ldv, 
    float *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, 
    magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz, 
    float *work);

magma_minproduct_int_t
magma_minproduct_sormqr_gpu_2stages(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dC, magma_minproduct_int_t lddc,
    magma_minproductFloat_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

// used only for old version and internal
magma_minproduct_int_t
magma_minproduct_ssytrd_bsy2trc_v5(
    magma_minproduct_int_t threads, magma_minproduct_int_t wantz, magma_minproduct_uplo_t uplo, 
    magma_minproduct_int_t ne, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    float *A, magma_minproduct_int_t lda, 
    float *D, float *E,
    magma_minproductFloat_ptr dT1, magma_minproduct_int_t ldt1);

magma_minproduct_int_t
magma_minproduct_sorgqr_2stage_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *tau,
    magma_minproductFloat_ptr dT,
    magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sbulge_get_lq2(
    magma_minproduct_int_t n, magma_minproduct_int_t threads);

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_minproduct_SBULGE_H */
