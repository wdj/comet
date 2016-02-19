/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
*/

#ifndef MAGMA_tally3_ZBULGE_H
#define MAGMA_tally3_ZBULGE_H

#include "magma_tally3_types.h"

#ifdef __cplusplus
extern "C" {
#endif

magma_tally3_int_t
magma_tally3_zbulge_applyQ_v2(
    magma_tally3_side_t side, 
    magma_tally3_int_t NE, magma_tally3_int_t N, 
    magma_tally3_int_t NB, magma_tally3_int_t Vblksiz, 
    magma_tally3DoubleComplex_ptr dE, magma_tally3_int_t ldde, 
    magma_tally3DoubleComplex *V, magma_tally3_int_t ldv, 
    magma_tally3DoubleComplex *T, magma_tally3_int_t ldt, 
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zbulge_applyQ_v2_m(
    magma_tally3_int_t ngpu, magma_tally3_side_t side, 
    magma_tally3_int_t NE, magma_tally3_int_t N, 
    magma_tally3_int_t NB, magma_tally3_int_t Vblksiz, 
    magma_tally3DoubleComplex *E, magma_tally3_int_t lde, 
    magma_tally3DoubleComplex *V, magma_tally3_int_t ldv, 
    magma_tally3DoubleComplex *T, magma_tally3_int_t ldt, 
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zbulge_back(
    magma_tally3_uplo_t uplo, 
    magma_tally3_int_t n, magma_tally3_int_t nb, 
    magma_tally3_int_t ne, magma_tally3_int_t Vblksiz,
    magma_tally3DoubleComplex *Z, magma_tally3_int_t ldz,
    magma_tally3DoubleComplex_ptr dZ, magma_tally3_int_t lddz,
    magma_tally3DoubleComplex *V, magma_tally3_int_t ldv,
    magma_tally3DoubleComplex *TAU,
    magma_tally3DoubleComplex *T, magma_tally3_int_t ldt,
    magma_tally3_int_t* info);

magma_tally3_int_t
magma_tally3_zbulge_back_m(
    magma_tally3_int_t ngpu, magma_tally3_uplo_t uplo, 
    magma_tally3_int_t n, magma_tally3_int_t nb, 
    magma_tally3_int_t ne, magma_tally3_int_t Vblksiz,
    magma_tally3DoubleComplex *Z, magma_tally3_int_t ldz,
    magma_tally3DoubleComplex *V, magma_tally3_int_t ldv, 
    magma_tally3DoubleComplex *TAU, 
    magma_tally3DoubleComplex *T, magma_tally3_int_t ldt, 
    magma_tally3_int_t* info);

void
magma_tally3_ztrdtype1cbHLsym_withQ_v2(
    magma_tally3_int_t n, magma_tally3_int_t nb, 
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda, 
    magma_tally3DoubleComplex *V, magma_tally3_int_t ldv, 
    magma_tally3DoubleComplex *TAU,
    magma_tally3_int_t st, magma_tally3_int_t ed, 
    magma_tally3_int_t sweep, magma_tally3_int_t Vblksiz, 
    magma_tally3DoubleComplex *work);

void
magma_tally3_ztrdtype2cbHLsym_withQ_v2(
    magma_tally3_int_t n, magma_tally3_int_t nb, 
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda, 
    magma_tally3DoubleComplex *V, magma_tally3_int_t ldv, 
    magma_tally3DoubleComplex *TAU,
    magma_tally3_int_t st, magma_tally3_int_t ed, 
    magma_tally3_int_t sweep, magma_tally3_int_t Vblksiz, 
    magma_tally3DoubleComplex *work);

void
magma_tally3_ztrdtype3cbHLsym_withQ_v2(
    magma_tally3_int_t n, magma_tally3_int_t nb, 
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda, 
    magma_tally3DoubleComplex *V, magma_tally3_int_t ldv, 
    magma_tally3DoubleComplex *TAU,
    magma_tally3_int_t st, magma_tally3_int_t ed, 
    magma_tally3_int_t sweep, magma_tally3_int_t Vblksiz, 
    magma_tally3DoubleComplex *work);

magma_tally3_int_t
magma_tally3_zunmqr_gpu_2stages(
    magma_tally3_side_t side, magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dC, magma_tally3_int_t lddc,
    magma_tally3DoubleComplex_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

// used only for old version and internal
magma_tally3_int_t
magma_tally3_zhetrd_bhe2trc_v5(
    magma_tally3_int_t threads, magma_tally3_int_t wantz, magma_tally3_uplo_t uplo, 
    magma_tally3_int_t ne, magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda, 
    double *D, double *E,
    magma_tally3DoubleComplex_ptr dT1, magma_tally3_int_t ldt1);

magma_tally3_int_t
magma_tally3_zungqr_2stage_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex_ptr dT,
    magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zbulge_get_lq2(
    magma_tally3_int_t n, magma_tally3_int_t threads);

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally3_ZBULGE_H */
