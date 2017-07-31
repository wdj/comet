/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
*/

#ifndef MAGMA_tally2_ZBULGE_H
#define MAGMA_tally2_ZBULGE_H

#include "magma_tally2_types.h"

#ifdef __cplusplus
extern "C" {
#endif

magma_tally2_int_t
magma_tally2_zbulge_applyQ_v2(
    magma_tally2_side_t side, 
    magma_tally2_int_t NE, magma_tally2_int_t N, 
    magma_tally2_int_t NB, magma_tally2_int_t Vblksiz, 
    magma_tally2DoubleComplex_ptr dE, magma_tally2_int_t ldde, 
    magma_tally2DoubleComplex *V, magma_tally2_int_t ldv, 
    magma_tally2DoubleComplex *T, magma_tally2_int_t ldt, 
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zbulge_applyQ_v2_m(
    magma_tally2_int_t ngpu, magma_tally2_side_t side, 
    magma_tally2_int_t NE, magma_tally2_int_t N, 
    magma_tally2_int_t NB, magma_tally2_int_t Vblksiz, 
    magma_tally2DoubleComplex *E, magma_tally2_int_t lde, 
    magma_tally2DoubleComplex *V, magma_tally2_int_t ldv, 
    magma_tally2DoubleComplex *T, magma_tally2_int_t ldt, 
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zbulge_back(
    magma_tally2_uplo_t uplo, 
    magma_tally2_int_t n, magma_tally2_int_t nb, 
    magma_tally2_int_t ne, magma_tally2_int_t Vblksiz,
    magma_tally2DoubleComplex *Z, magma_tally2_int_t ldz,
    magma_tally2DoubleComplex_ptr dZ, magma_tally2_int_t lddz,
    magma_tally2DoubleComplex *V, magma_tally2_int_t ldv,
    magma_tally2DoubleComplex *TAU,
    magma_tally2DoubleComplex *T, magma_tally2_int_t ldt,
    magma_tally2_int_t* info);

magma_tally2_int_t
magma_tally2_zbulge_back_m(
    magma_tally2_int_t ngpu, magma_tally2_uplo_t uplo, 
    magma_tally2_int_t n, magma_tally2_int_t nb, 
    magma_tally2_int_t ne, magma_tally2_int_t Vblksiz,
    magma_tally2DoubleComplex *Z, magma_tally2_int_t ldz,
    magma_tally2DoubleComplex *V, magma_tally2_int_t ldv, 
    magma_tally2DoubleComplex *TAU, 
    magma_tally2DoubleComplex *T, magma_tally2_int_t ldt, 
    magma_tally2_int_t* info);

void
magma_tally2_ztrdtype1cbHLsym_withQ_v2(
    magma_tally2_int_t n, magma_tally2_int_t nb, 
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda, 
    magma_tally2DoubleComplex *V, magma_tally2_int_t ldv, 
    magma_tally2DoubleComplex *TAU,
    magma_tally2_int_t st, magma_tally2_int_t ed, 
    magma_tally2_int_t sweep, magma_tally2_int_t Vblksiz, 
    magma_tally2DoubleComplex *work);

void
magma_tally2_ztrdtype2cbHLsym_withQ_v2(
    magma_tally2_int_t n, magma_tally2_int_t nb, 
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda, 
    magma_tally2DoubleComplex *V, magma_tally2_int_t ldv, 
    magma_tally2DoubleComplex *TAU,
    magma_tally2_int_t st, magma_tally2_int_t ed, 
    magma_tally2_int_t sweep, magma_tally2_int_t Vblksiz, 
    magma_tally2DoubleComplex *work);

void
magma_tally2_ztrdtype3cbHLsym_withQ_v2(
    magma_tally2_int_t n, magma_tally2_int_t nb, 
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda, 
    magma_tally2DoubleComplex *V, magma_tally2_int_t ldv, 
    magma_tally2DoubleComplex *TAU,
    magma_tally2_int_t st, magma_tally2_int_t ed, 
    magma_tally2_int_t sweep, magma_tally2_int_t Vblksiz, 
    magma_tally2DoubleComplex *work);

magma_tally2_int_t
magma_tally2_zunmqr_gpu_2stages(
    magma_tally2_side_t side, magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dC, magma_tally2_int_t lddc,
    magma_tally2DoubleComplex_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

// used only for old version and internal
magma_tally2_int_t
magma_tally2_zhetrd_bhe2trc_v5(
    magma_tally2_int_t threads, magma_tally2_int_t wantz, magma_tally2_uplo_t uplo, 
    magma_tally2_int_t ne, magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda, 
    double *D, double *E,
    magma_tally2DoubleComplex_ptr dT1, magma_tally2_int_t ldt1);

magma_tally2_int_t
magma_tally2_zungqr_2stage_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex_ptr dT,
    magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zbulge_get_lq2(
    magma_tally2_int_t n, magma_tally2_int_t threads);

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally2_ZBULGE_H */
