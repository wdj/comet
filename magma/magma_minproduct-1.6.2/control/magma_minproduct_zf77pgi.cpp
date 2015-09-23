/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Mathieu Faverge
*/
#include "common_magma_minproduct.h"

#undef REAL
#define COMPLEX

#ifdef ADD_

    #define MAGMA_minproduct_ZGEBRD magma_minproduct_zgebrd_
    #define MAGMA_minproduct_ZGEHRD magma_minproduct_zgehrd_
    #define MAGMA_minproduct_ZGELQF magma_minproduct_zgelqf_
    #define MAGMA_minproduct_ZGEQLF magma_minproduct_zgeqlf_
    #define MAGMA_minproduct_ZGEQRF magma_minproduct_zgeqrf_
    #define MAGMA_minproduct_ZGETRF magma_minproduct_zgetrf_
    #define MAGMA_minproduct_ZLABRD magma_minproduct_zlabrd_
    #define MAGMA_minproduct_ZLAHR2 magma_minproduct_zlahr2_
    #define MAGMA_minproduct_ZLAHRU magma_minproduct_zlahru_
    #define MAGMA_minproduct_ZPOTRF magma_minproduct_zpotrf_
    #define MAGMA_minproduct_ZHETRD magma_minproduct_zhetrd_
    
    #define MAGMA_minproduct_ZUNMQR_GPU magma_minproduct_zunmqr_gpu_
    #define MAGMA_minproduct_ZGEQRF_GPU  magma_minproduct_zgeqrf_gpu_
    #define MAGMA_minproduct_ZGEQRF2_GPU magma_minproduct_zgeqrf2_gpu_
    #define MAGMA_minproduct_ZGEQRS_GPU magma_minproduct_zgeqrs_gpu_
    #define MAGMA_minproduct_ZGETRF_GPU magma_minproduct_zgetrf_gpu_
    #define MAGMA_minproduct_ZGETRS_GPU magma_minproduct_zgetrs_gpu_
    #define MAGMA_minproduct_ZLARFB_GPU magma_minproduct_zlarfb_gpu_
    #define MAGMA_minproduct_ZPOTRF_GPU magma_minproduct_zpotrf_gpu_
    #define MAGMA_minproduct_ZPOTRS_GPU magma_minproduct_zpotrs_gpu_

#elif defined (NOCHANGE)

    #define MAGMA_minproduct_ZGEBRD magma_minproduct_zgebrd
    #define MAGMA_minproduct_ZGEHRD magma_minproduct_zgehrd
    #define MAGMA_minproduct_ZGELQF magma_minproduct_zgelqf
    #define MAGMA_minproduct_ZGEQLF magma_minproduct_zgeqlf
    #define MAGMA_minproduct_ZGEQRF magma_minproduct_zgeqrf
    #define MAGMA_minproduct_ZGETRF magma_minproduct_zgetrf
    #define MAGMA_minproduct_ZLABRD magma_minproduct_zlabrd
    #define MAGMA_minproduct_ZLAHR2 magma_minproduct_zlahr2
    #define MAGMA_minproduct_ZLAHRU magma_minproduct_zlahru
    #define MAGMA_minproduct_ZPOTRF magma_minproduct_zpotrf
    #define MAGMA_minproduct_ZHETRD magma_minproduct_zhetrd
    
    #define MAGMA_minproduct_ZUNMQR_GPU magma_minproduct_zunmqr_gpu
    #define MAGMA_minproduct_ZGEQRF_GPU  magma_minproduct_zgeqrf_gpu
    #define MAGMA_minproduct_ZGEQRF2_GPU magma_minproduct_zgeqrf2_gpu
    #define MAGMA_minproduct_ZGEQRS_GPU magma_minproduct_zgeqrs_gpu
    #define MAGMA_minproduct_ZGETRF_GPU magma_minproduct_zgetrf_gpu
    #define MAGMA_minproduct_ZGETRS_GPU magma_minproduct_zgetrs_gpu
    #define MAGMA_minproduct_ZLARFB_GPU magma_minproduct_zlarfb_gpu
    #define MAGMA_minproduct_ZPOTRF_GPU magma_minproduct_zpotrf_gpu
    #define MAGMA_minproduct_ZPOTRS_GPU magma_minproduct_zpotrs_gpu

#endif

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (simple interface)
 **/

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct function definitions / Data on CPU
*/
void MAGMA_minproduct_ZGEBRD( magma_minproduct_int_t *m, magma_minproduct_int_t *n, double2 *A, magma_minproduct_int_t *lda, double *d, double *e, double2 *tauq,  double2 *taup, double2 *work, magma_minproduct_int_t *lwork, double2 *da, magma_minproduct_int_t *info)
{ magma_minproduct_zgebrd( *m, *n, A, *lda, d, e, tauq,  taup, work, *lwork, da, info); }

void MAGMA_minproduct_ZGEHRD( magma_minproduct_int_t *n, magma_minproduct_int_t *ilo, magma_minproduct_int_t *ihi, double2 *A, magma_minproduct_int_t *lda, double2 *tau, double2 *work, magma_minproduct_int_t *lwork, double2 *da, magma_minproduct_int_t *info)
{ magma_minproduct_zgehrd( *n, *ilo, *ihi, A, *lda, tau, work, lwork, da, info); }

void MAGMA_minproduct_ZGELQF( magma_minproduct_int_t *m, magma_minproduct_int_t *n, double2 *A, magma_minproduct_int_t *lda, double2 *tau, double2 *work, magma_minproduct_int_t *lwork, magma_minproduct_int_t *info)
{ magma_minproduct_zgelqf( *m, *n, A, *lda, tau, work, *lwork, info); }

void MAGMA_minproduct_ZGEQLF( magma_minproduct_int_t *m, magma_minproduct_int_t *n, double2 *A, magma_minproduct_int_t *lda, double2 *tau, double2 *work, magma_minproduct_int_t *lwork, magma_minproduct_int_t *info)
{ magma_minproduct_zgeqlf( *m, *n, A, *lda, tau, work, *lwork, info); }

void MAGMA_minproduct_ZGEQRF( magma_minproduct_int_t *m, magma_minproduct_int_t *n, double2 *A, magma_minproduct_int_t *lda, double2 *tau, double2 *work, magma_minproduct_int_t *lwork, magma_minproduct_int_t *info)
{ magma_minproduct_zgeqrf( *m, *n, A, *lda, tau, work, *lwork, info); }

void MAGMA_minproduct_ZGETRF( magma_minproduct_int_t *m, magma_minproduct_int_t *n, double2 *A, magma_minproduct_int_t *lda, magma_minproduct_int_t *ipiv, magma_minproduct_int_t *info)
{ magma_minproduct_zgetrf( *m, *n, A, *lda, ipiv, info); }

void MAGMA_minproduct_ZLABRD( magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *nb, double2 *a, magma_minproduct_int_t *lda, double *d, double *e, double2 *tauq, double2 *taup, double2 *x, magma_minproduct_int_t *ldx, double2 *y, magma_minproduct_int_t *ldy, double2 *da, magma_minproduct_int_t *ldda, double2 *dx, magma_minproduct_int_t *lddx, double2 *dy, magma_minproduct_int_t *lddy)
{ magma_minproduct_zlabrd( *m, *n, *nb, a, *lda, d, e, tauq, taup, x, *ldx, y, *ldy, da, *ldda, dx, *lddx, dy, *lddy); }

void MAGMA_minproduct_ZLAHR2( magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *nb, double2 *da, double2 *dv, double2 *a, magma_minproduct_int_t *lda, double2 *tau, double2 *t, magma_minproduct_int_t *ldt, double2 *y, magma_minproduct_int_t *ldy)
{ magma_minproduct_zlahr2( *m, *n, *nb, da, dv, a, *lda, tau, t, *ldt, y, *ldy); }

void MAGMA_minproduct_ZLAHRU( magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *nb, double2 *a, magma_minproduct_int_t *lda, double2 *da, double2 *y, double2 *v, double2 *t, double2 *dwork)
{ magma_minproduct_zlahru( *m, *n, *nb, a, *lda, da, y, v, t, dwork); }

void MAGMA_minproduct_ZPOTRF( char *uplo, magma_minproduct_int_t *n, double2 *A, magma_minproduct_int_t *lda, magma_minproduct_int_t *info)
{ magma_minproduct_zpotrf( uplo[0], *n, A, *lda, info); }

void MAGMA_minproduct_ZHETRD( char *uplo, magma_minproduct_int_t *n, double2 *A, magma_minproduct_int_t *lda, double *d, double *e, double2 *tau, double2 *work, magma_minproduct_int_t *lwork, double2 *da, magma_minproduct_int_t *info)
{ magma_minproduct_zhetrd( uplo[0], *n, A, *lda, d, e, tau, work, *lwork, da, info); }


/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA_minproduct function definitions / Data on GPU
*/
void MAGMA_minproduct_ZUNMQR_GPU(char *side, char *trans, magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *k, double2 *a, magma_minproduct_int_t *lda, double2 *tau, double2 *c, magma_minproduct_int_t *ldc, double2 *work, magma_minproduct_int_t *lwork, double2 *td, magma_minproduct_int_t *nb, magma_minproduct_int_t *info)
{ magma_minproduct_zunmqr_gpu(side[0], trans[0], *m, *n, *k, a, *lda, tau, c, *ldc, work, *lwork, td, *nb, info); }

void MAGMA_minproduct_ZGEQRF2_GPU( magma_minproduct_int_t *m, magma_minproduct_int_t *n, double2 *A, magma_minproduct_int_t *lda, double2 *tau, magma_minproduct_int_t *info)
{ magma_minproduct_zgeqrf2_gpu( *m, *n, A, *lda, tau, info); }

void MAGMA_minproduct_ZGEQRF_GPU(magma_minproduct_int_t *m, magma_minproduct_int_t *n, double2 *A, magma_minproduct_int_t *lda, double2 *tau, double2 *dwork, magma_minproduct_int_t *info)
{ magma_minproduct_zgeqrf_gpu(*m, *n, A, *lda, tau, dwork, info); }

void MAGMA_minproduct_ZGEQRS_GPU( magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *nrhs, double2 *A, magma_minproduct_int_t *lda, double2 *tau, double2 *td, double2 *c, magma_minproduct_int_t *ldc, double2 *work, magma_minproduct_int_t *lwork, magma_minproduct_int_t *info)
{ magma_minproduct_zgeqrs_gpu( *m, *n, *nrhs, A, *lda, tau, td, c, *ldc, work, *lwork, info); }

void MAGMA_minproduct_ZGETRF_GPU( magma_minproduct_int_t *m, magma_minproduct_int_t *n, double2 *A, magma_minproduct_int_t *lda, magma_minproduct_int_t *ipiv, magma_minproduct_int_t *info)
{ magma_minproduct_zgetrf_gpu( *m, *n, A, *lda, ipiv, info); }

void MAGMA_minproduct_ZGETRS_GPU( char *trans, magma_minproduct_int_t *n, magma_minproduct_int_t *nrhs, double2 *A, magma_minproduct_int_t *lda, magma_minproduct_int_t *ipiv, double2 *b, magma_minproduct_int_t *ldb, magma_minproduct_int_t *info)
{ magma_minproduct_zgetrs_gpu( trans[0], *n, *nrhs, A, *lda, ipiv, b, *ldb, info ); }

void MAGMA_minproduct_ZLARFB_GPU( char *side, char *trans, char *direct, char *storev, magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *k, double2 *dv, magma_minproduct_int_t *ldv, double2 *dt, magma_minproduct_int_t *ldt, double2 *dc, magma_minproduct_int_t *ldc, double2 *dowrk, magma_minproduct_int_t *ldwork)
{ magma_minproduct_zlarfb_gpu( side[0], trans[0], direct[0], storev[0], *m, *n, *k, dv, *ldv, dt, *ldt, dc, *ldc, dowrk, *ldwork); }

void MAGMA_minproduct_ZPOTRF_GPU( char *uplo,  magma_minproduct_int_t *n, double2 *A, magma_minproduct_int_t *lda, magma_minproduct_int_t *info)
{ magma_minproduct_zpotrf_gpu( uplo[0], *n, A, *lda, info); }

void MAGMA_minproduct_ZPOTRS_GPU( char *uplo,  magma_minproduct_int_t *n, magma_minproduct_int_t *nrhs, double2 *A, magma_minproduct_int_t *lda, double2 *b, magma_minproduct_int_t *ldb, magma_minproduct_int_t *info)
{ magma_minproduct_zpotrs_gpu( uplo[0], *n, *nrhs, A, *lda, b, *ldb, info); }

#ifdef __cplusplus
}
#endif
