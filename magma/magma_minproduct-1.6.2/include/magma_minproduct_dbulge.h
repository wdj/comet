/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproduct_zbulge.h normal z -> d, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_minproduct_DBULGE_H
#define MAGMA_minproduct_DBULGE_H

#include "magma_minproduct_types.h"

#ifdef __cplusplus
extern "C" {
#endif

magma_minproduct_int_t
magma_minproduct_dbulge_applyQ_v2(
    magma_minproduct_side_t side, 
    magma_minproduct_int_t NE, magma_minproduct_int_t N, 
    magma_minproduct_int_t NB, magma_minproduct_int_t Vblksiz, 
    magma_minproductDouble_ptr dE, magma_minproduct_int_t ldde, 
    double *V, magma_minproduct_int_t ldv, 
    double *T, magma_minproduct_int_t ldt, 
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dbulge_applyQ_v2_m(
    magma_minproduct_int_t ngpu, magma_minproduct_side_t side, 
    magma_minproduct_int_t NE, magma_minproduct_int_t N, 
    magma_minproduct_int_t NB, magma_minproduct_int_t Vblksiz, 
    double *E, magma_minproduct_int_t lde, 
    double *V, magma_minproduct_int_t ldv, 
    double *T, magma_minproduct_int_t ldt, 
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dbulge_back(
    magma_minproduct_uplo_t uplo, 
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    magma_minproduct_int_t ne, magma_minproduct_int_t Vblksiz,
    double *Z, magma_minproduct_int_t ldz,
    magma_minproductDouble_ptr dZ, magma_minproduct_int_t lddz,
    double *V, magma_minproduct_int_t ldv,
    double *TAU,
    double *T, magma_minproduct_int_t ldt,
    magma_minproduct_int_t* info);

magma_minproduct_int_t
magma_minproduct_dbulge_back_m(
    magma_minproduct_int_t ngpu, magma_minproduct_uplo_t uplo, 
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    magma_minproduct_int_t ne, magma_minproduct_int_t Vblksiz,
    double *Z, magma_minproduct_int_t ldz,
    double *V, magma_minproduct_int_t ldv, 
    double *TAU, 
    double *T, magma_minproduct_int_t ldt, 
    magma_minproduct_int_t* info);

void
magma_minproduct_dtrdtype1cbHLsym_withQ_v2(
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    double *A, magma_minproduct_int_t lda, 
    double *V, magma_minproduct_int_t ldv, 
    double *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, 
    magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz, 
    double *work);

void
magma_minproduct_dtrdtype2cbHLsym_withQ_v2(
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    double *A, magma_minproduct_int_t lda, 
    double *V, magma_minproduct_int_t ldv, 
    double *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, 
    magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz, 
    double *work);

void
magma_minproduct_dtrdtype3cbHLsym_withQ_v2(
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    double *A, magma_minproduct_int_t lda, 
    double *V, magma_minproduct_int_t ldv, 
    double *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, 
    magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz, 
    double *work);

magma_minproduct_int_t
magma_minproduct_dormqr_gpu_2stages(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dC, magma_minproduct_int_t lddc,
    magma_minproductDouble_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

// used only for old version and internal
magma_minproduct_int_t
magma_minproduct_dsytrd_bsy2trc_v5(
    magma_minproduct_int_t threads, magma_minproduct_int_t wantz, magma_minproduct_uplo_t uplo, 
    magma_minproduct_int_t ne, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    double *A, magma_minproduct_int_t lda, 
    double *D, double *E,
    magma_minproductDouble_ptr dT1, magma_minproduct_int_t ldt1);

magma_minproduct_int_t
magma_minproduct_dorgqr_2stage_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *tau,
    magma_minproductDouble_ptr dT,
    magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dbulge_get_lq2(
    magma_minproduct_int_t n, magma_minproduct_int_t threads);

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_minproduct_DBULGE_H */
