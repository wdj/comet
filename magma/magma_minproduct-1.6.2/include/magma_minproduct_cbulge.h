/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproduct_zbulge.h normal z -> c, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_minproduct_CBULGE_H
#define MAGMA_minproduct_CBULGE_H

#include "magma_minproduct_types.h"

#ifdef __cplusplus
extern "C" {
#endif

magma_minproduct_int_t
magma_minproduct_cbulge_applyQ_v2(
    magma_minproduct_side_t side, 
    magma_minproduct_int_t NE, magma_minproduct_int_t N, 
    magma_minproduct_int_t NB, magma_minproduct_int_t Vblksiz, 
    magma_minproductFloatComplex_ptr dE, magma_minproduct_int_t ldde, 
    magma_minproductFloatComplex *V, magma_minproduct_int_t ldv, 
    magma_minproductFloatComplex *T, magma_minproduct_int_t ldt, 
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cbulge_applyQ_v2_m(
    magma_minproduct_int_t ngpu, magma_minproduct_side_t side, 
    magma_minproduct_int_t NE, magma_minproduct_int_t N, 
    magma_minproduct_int_t NB, magma_minproduct_int_t Vblksiz, 
    magma_minproductFloatComplex *E, magma_minproduct_int_t lde, 
    magma_minproductFloatComplex *V, magma_minproduct_int_t ldv, 
    magma_minproductFloatComplex *T, magma_minproduct_int_t ldt, 
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cbulge_back(
    magma_minproduct_uplo_t uplo, 
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    magma_minproduct_int_t ne, magma_minproduct_int_t Vblksiz,
    magma_minproductFloatComplex *Z, magma_minproduct_int_t ldz,
    magma_minproductFloatComplex_ptr dZ, magma_minproduct_int_t lddz,
    magma_minproductFloatComplex *V, magma_minproduct_int_t ldv,
    magma_minproductFloatComplex *TAU,
    magma_minproductFloatComplex *T, magma_minproduct_int_t ldt,
    magma_minproduct_int_t* info);

magma_minproduct_int_t
magma_minproduct_cbulge_back_m(
    magma_minproduct_int_t ngpu, magma_minproduct_uplo_t uplo, 
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    magma_minproduct_int_t ne, magma_minproduct_int_t Vblksiz,
    magma_minproductFloatComplex *Z, magma_minproduct_int_t ldz,
    magma_minproductFloatComplex *V, magma_minproduct_int_t ldv, 
    magma_minproductFloatComplex *TAU, 
    magma_minproductFloatComplex *T, magma_minproduct_int_t ldt, 
    magma_minproduct_int_t* info);

void
magma_minproduct_ctrdtype1cbHLsym_withQ_v2(
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda, 
    magma_minproductFloatComplex *V, magma_minproduct_int_t ldv, 
    magma_minproductFloatComplex *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, 
    magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz, 
    magma_minproductFloatComplex *work);

void
magma_minproduct_ctrdtype2cbHLsym_withQ_v2(
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda, 
    magma_minproductFloatComplex *V, magma_minproduct_int_t ldv, 
    magma_minproductFloatComplex *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, 
    magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz, 
    magma_minproductFloatComplex *work);

void
magma_minproduct_ctrdtype3cbHLsym_withQ_v2(
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda, 
    magma_minproductFloatComplex *V, magma_minproduct_int_t ldv, 
    magma_minproductFloatComplex *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, 
    magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz, 
    magma_minproductFloatComplex *work);

magma_minproduct_int_t
magma_minproduct_cunmqr_gpu_2stages(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dC, magma_minproduct_int_t lddc,
    magma_minproductFloatComplex_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

// used only for old version and internal
magma_minproduct_int_t
magma_minproduct_chetrd_bhe2trc_v5(
    magma_minproduct_int_t threads, magma_minproduct_int_t wantz, magma_minproduct_uplo_t uplo, 
    magma_minproduct_int_t ne, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda, 
    float *D, float *E,
    magma_minproductFloatComplex_ptr dT1, magma_minproduct_int_t ldt1);

magma_minproduct_int_t
magma_minproduct_cungqr_2stage_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex_ptr dT,
    magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cbulge_get_lq2(
    magma_minproduct_int_t n, magma_minproduct_int_t threads);

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_minproduct_CBULGE_H */
