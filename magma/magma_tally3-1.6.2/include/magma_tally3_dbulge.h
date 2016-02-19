/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally3_zbulge.h normal z -> d, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_tally3_DBULGE_H
#define MAGMA_tally3_DBULGE_H

#include "magma_tally3_types.h"

#ifdef __cplusplus
extern "C" {
#endif

magma_tally3_int_t
magma_tally3_dbulge_applyQ_v2(
    magma_tally3_side_t side, 
    magma_tally3_int_t NE, magma_tally3_int_t N, 
    magma_tally3_int_t NB, magma_tally3_int_t Vblksiz, 
    magma_tally3Double_ptr dE, magma_tally3_int_t ldde, 
    double *V, magma_tally3_int_t ldv, 
    double *T, magma_tally3_int_t ldt, 
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dbulge_applyQ_v2_m(
    magma_tally3_int_t ngpu, magma_tally3_side_t side, 
    magma_tally3_int_t NE, magma_tally3_int_t N, 
    magma_tally3_int_t NB, magma_tally3_int_t Vblksiz, 
    double *E, magma_tally3_int_t lde, 
    double *V, magma_tally3_int_t ldv, 
    double *T, magma_tally3_int_t ldt, 
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dbulge_back(
    magma_tally3_uplo_t uplo, 
    magma_tally3_int_t n, magma_tally3_int_t nb, 
    magma_tally3_int_t ne, magma_tally3_int_t Vblksiz,
    double *Z, magma_tally3_int_t ldz,
    magma_tally3Double_ptr dZ, magma_tally3_int_t lddz,
    double *V, magma_tally3_int_t ldv,
    double *TAU,
    double *T, magma_tally3_int_t ldt,
    magma_tally3_int_t* info);

magma_tally3_int_t
magma_tally3_dbulge_back_m(
    magma_tally3_int_t ngpu, magma_tally3_uplo_t uplo, 
    magma_tally3_int_t n, magma_tally3_int_t nb, 
    magma_tally3_int_t ne, magma_tally3_int_t Vblksiz,
    double *Z, magma_tally3_int_t ldz,
    double *V, magma_tally3_int_t ldv, 
    double *TAU, 
    double *T, magma_tally3_int_t ldt, 
    magma_tally3_int_t* info);

void
magma_tally3_dtrdtype1cbHLsym_withQ_v2(
    magma_tally3_int_t n, magma_tally3_int_t nb, 
    double *A, magma_tally3_int_t lda, 
    double *V, magma_tally3_int_t ldv, 
    double *TAU,
    magma_tally3_int_t st, magma_tally3_int_t ed, 
    magma_tally3_int_t sweep, magma_tally3_int_t Vblksiz, 
    double *work);

void
magma_tally3_dtrdtype2cbHLsym_withQ_v2(
    magma_tally3_int_t n, magma_tally3_int_t nb, 
    double *A, magma_tally3_int_t lda, 
    double *V, magma_tally3_int_t ldv, 
    double *TAU,
    magma_tally3_int_t st, magma_tally3_int_t ed, 
    magma_tally3_int_t sweep, magma_tally3_int_t Vblksiz, 
    double *work);

void
magma_tally3_dtrdtype3cbHLsym_withQ_v2(
    magma_tally3_int_t n, magma_tally3_int_t nb, 
    double *A, magma_tally3_int_t lda, 
    double *V, magma_tally3_int_t ldv, 
    double *TAU,
    magma_tally3_int_t st, magma_tally3_int_t ed, 
    magma_tally3_int_t sweep, magma_tally3_int_t Vblksiz, 
    double *work);

magma_tally3_int_t
magma_tally3_dormqr_gpu_2stages(
    magma_tally3_side_t side, magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dC, magma_tally3_int_t lddc,
    magma_tally3Double_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

// used only for old version and internal
magma_tally3_int_t
magma_tally3_dsytrd_bsy2trc_v5(
    magma_tally3_int_t threads, magma_tally3_int_t wantz, magma_tally3_uplo_t uplo, 
    magma_tally3_int_t ne, magma_tally3_int_t n, magma_tally3_int_t nb,
    double *A, magma_tally3_int_t lda, 
    double *D, double *E,
    magma_tally3Double_ptr dT1, magma_tally3_int_t ldt1);

magma_tally3_int_t
magma_tally3_dorgqr_2stage_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *tau,
    magma_tally3Double_ptr dT,
    magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dbulge_get_lq2(
    magma_tally3_int_t n, magma_tally3_int_t threads);

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally3_DBULGE_H */
