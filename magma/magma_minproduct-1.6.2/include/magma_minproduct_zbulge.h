/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
*/

#ifndef MAGMA_minproduct_ZBULGE_H
#define MAGMA_minproduct_ZBULGE_H

#include "magma_minproduct_types.h"

#ifdef __cplusplus
extern "C" {
#endif

magma_minproduct_int_t
magma_minproduct_zbulge_applyQ_v2(
    magma_minproduct_side_t side, 
    magma_minproduct_int_t NE, magma_minproduct_int_t N, 
    magma_minproduct_int_t NB, magma_minproduct_int_t Vblksiz, 
    magma_minproductDoubleComplex_ptr dE, magma_minproduct_int_t ldde, 
    magma_minproductDoubleComplex *V, magma_minproduct_int_t ldv, 
    magma_minproductDoubleComplex *T, magma_minproduct_int_t ldt, 
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zbulge_applyQ_v2_m(
    magma_minproduct_int_t ngpu, magma_minproduct_side_t side, 
    magma_minproduct_int_t NE, magma_minproduct_int_t N, 
    magma_minproduct_int_t NB, magma_minproduct_int_t Vblksiz, 
    magma_minproductDoubleComplex *E, magma_minproduct_int_t lde, 
    magma_minproductDoubleComplex *V, magma_minproduct_int_t ldv, 
    magma_minproductDoubleComplex *T, magma_minproduct_int_t ldt, 
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zbulge_back(
    magma_minproduct_uplo_t uplo, 
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    magma_minproduct_int_t ne, magma_minproduct_int_t Vblksiz,
    magma_minproductDoubleComplex *Z, magma_minproduct_int_t ldz,
    magma_minproductDoubleComplex_ptr dZ, magma_minproduct_int_t lddz,
    magma_minproductDoubleComplex *V, magma_minproduct_int_t ldv,
    magma_minproductDoubleComplex *TAU,
    magma_minproductDoubleComplex *T, magma_minproduct_int_t ldt,
    magma_minproduct_int_t* info);

magma_minproduct_int_t
magma_minproduct_zbulge_back_m(
    magma_minproduct_int_t ngpu, magma_minproduct_uplo_t uplo, 
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    magma_minproduct_int_t ne, magma_minproduct_int_t Vblksiz,
    magma_minproductDoubleComplex *Z, magma_minproduct_int_t ldz,
    magma_minproductDoubleComplex *V, magma_minproduct_int_t ldv, 
    magma_minproductDoubleComplex *TAU, 
    magma_minproductDoubleComplex *T, magma_minproduct_int_t ldt, 
    magma_minproduct_int_t* info);

void
magma_minproduct_ztrdtype1cbHLsym_withQ_v2(
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda, 
    magma_minproductDoubleComplex *V, magma_minproduct_int_t ldv, 
    magma_minproductDoubleComplex *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, 
    magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz, 
    magma_minproductDoubleComplex *work);

void
magma_minproduct_ztrdtype2cbHLsym_withQ_v2(
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda, 
    magma_minproductDoubleComplex *V, magma_minproduct_int_t ldv, 
    magma_minproductDoubleComplex *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, 
    magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz, 
    magma_minproductDoubleComplex *work);

void
magma_minproduct_ztrdtype3cbHLsym_withQ_v2(
    magma_minproduct_int_t n, magma_minproduct_int_t nb, 
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda, 
    magma_minproductDoubleComplex *V, magma_minproduct_int_t ldv, 
    magma_minproductDoubleComplex *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, 
    magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz, 
    magma_minproductDoubleComplex *work);

magma_minproduct_int_t
magma_minproduct_zunmqr_gpu_2stages(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dC, magma_minproduct_int_t lddc,
    magma_minproductDoubleComplex_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

// used only for old version and internal
magma_minproduct_int_t
magma_minproduct_zhetrd_bhe2trc_v5(
    magma_minproduct_int_t threads, magma_minproduct_int_t wantz, magma_minproduct_uplo_t uplo, 
    magma_minproduct_int_t ne, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda, 
    double *D, double *E,
    magma_minproductDoubleComplex_ptr dT1, magma_minproduct_int_t ldt1);

magma_minproduct_int_t
magma_minproduct_zungqr_2stage_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex_ptr dT,
    magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zbulge_get_lq2(
    magma_minproduct_int_t n, magma_minproduct_int_t threads);

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_minproduct_ZBULGE_H */
