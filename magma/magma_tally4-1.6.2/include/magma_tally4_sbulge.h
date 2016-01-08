/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally4_zbulge.h normal z -> s, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_tally4_SBULGE_H
#define MAGMA_tally4_SBULGE_H

#include "magma_tally4_types.h"

#ifdef __cplusplus
extern "C" {
#endif

magma_tally4_int_t
magma_tally4_sbulge_applyQ_v2(
    magma_tally4_side_t side, 
    magma_tally4_int_t NE, magma_tally4_int_t N, 
    magma_tally4_int_t NB, magma_tally4_int_t Vblksiz, 
    magma_tally4Float_ptr dE, magma_tally4_int_t ldde, 
    float *V, magma_tally4_int_t ldv, 
    float *T, magma_tally4_int_t ldt, 
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sbulge_applyQ_v2_m(
    magma_tally4_int_t ngpu, magma_tally4_side_t side, 
    magma_tally4_int_t NE, magma_tally4_int_t N, 
    magma_tally4_int_t NB, magma_tally4_int_t Vblksiz, 
    float *E, magma_tally4_int_t lde, 
    float *V, magma_tally4_int_t ldv, 
    float *T, magma_tally4_int_t ldt, 
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sbulge_back(
    magma_tally4_uplo_t uplo, 
    magma_tally4_int_t n, magma_tally4_int_t nb, 
    magma_tally4_int_t ne, magma_tally4_int_t Vblksiz,
    float *Z, magma_tally4_int_t ldz,
    magma_tally4Float_ptr dZ, magma_tally4_int_t lddz,
    float *V, magma_tally4_int_t ldv,
    float *TAU,
    float *T, magma_tally4_int_t ldt,
    magma_tally4_int_t* info);

magma_tally4_int_t
magma_tally4_sbulge_back_m(
    magma_tally4_int_t ngpu, magma_tally4_uplo_t uplo, 
    magma_tally4_int_t n, magma_tally4_int_t nb, 
    magma_tally4_int_t ne, magma_tally4_int_t Vblksiz,
    float *Z, magma_tally4_int_t ldz,
    float *V, magma_tally4_int_t ldv, 
    float *TAU, 
    float *T, magma_tally4_int_t ldt, 
    magma_tally4_int_t* info);

void
magma_tally4_strdtype1cbHLsym_withQ_v2(
    magma_tally4_int_t n, magma_tally4_int_t nb, 
    float *A, magma_tally4_int_t lda, 
    float *V, magma_tally4_int_t ldv, 
    float *TAU,
    magma_tally4_int_t st, magma_tally4_int_t ed, 
    magma_tally4_int_t sweep, magma_tally4_int_t Vblksiz, 
    float *work);

void
magma_tally4_strdtype2cbHLsym_withQ_v2(
    magma_tally4_int_t n, magma_tally4_int_t nb, 
    float *A, magma_tally4_int_t lda, 
    float *V, magma_tally4_int_t ldv, 
    float *TAU,
    magma_tally4_int_t st, magma_tally4_int_t ed, 
    magma_tally4_int_t sweep, magma_tally4_int_t Vblksiz, 
    float *work);

void
magma_tally4_strdtype3cbHLsym_withQ_v2(
    magma_tally4_int_t n, magma_tally4_int_t nb, 
    float *A, magma_tally4_int_t lda, 
    float *V, magma_tally4_int_t ldv, 
    float *TAU,
    magma_tally4_int_t st, magma_tally4_int_t ed, 
    magma_tally4_int_t sweep, magma_tally4_int_t Vblksiz, 
    float *work);

magma_tally4_int_t
magma_tally4_sormqr_gpu_2stages(
    magma_tally4_side_t side, magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dC, magma_tally4_int_t lddc,
    magma_tally4Float_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

// used only for old version and internal
magma_tally4_int_t
magma_tally4_ssytrd_bsy2trc_v5(
    magma_tally4_int_t threads, magma_tally4_int_t wantz, magma_tally4_uplo_t uplo, 
    magma_tally4_int_t ne, magma_tally4_int_t n, magma_tally4_int_t nb,
    float *A, magma_tally4_int_t lda, 
    float *D, float *E,
    magma_tally4Float_ptr dT1, magma_tally4_int_t ldt1);

magma_tally4_int_t
magma_tally4_sorgqr_2stage_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *tau,
    magma_tally4Float_ptr dT,
    magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sbulge_get_lq2(
    magma_tally4_int_t n, magma_tally4_int_t threads);

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally4_SBULGE_H */
