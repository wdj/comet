/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Mathieu Faverge
*/
#include "common_magma_tally3.h"

#undef REAL
#define COMPLEX

#ifdef ADD_

    #define MAGMA_tally3_ZGEBRD magma_tally3_zgebrd_
    #define MAGMA_tally3_ZGEHRD magma_tally3_zgehrd_
    #define MAGMA_tally3_ZGELQF magma_tally3_zgelqf_
    #define MAGMA_tally3_ZGEQLF magma_tally3_zgeqlf_
    #define MAGMA_tally3_ZGEQRF magma_tally3_zgeqrf_
    #define MAGMA_tally3_ZGETRF magma_tally3_zgetrf_
    #define MAGMA_tally3_ZLABRD magma_tally3_zlabrd_
    #define MAGMA_tally3_ZLAHR2 magma_tally3_zlahr2_
    #define MAGMA_tally3_ZLAHRU magma_tally3_zlahru_
    #define MAGMA_tally3_ZPOTRF magma_tally3_zpotrf_
    #define MAGMA_tally3_ZHETRD magma_tally3_zhetrd_
    
    #define MAGMA_tally3_ZUNMQR_GPU magma_tally3_zunmqr_gpu_
    #define MAGMA_tally3_ZGEQRF_GPU  magma_tally3_zgeqrf_gpu_
    #define MAGMA_tally3_ZGEQRF2_GPU magma_tally3_zgeqrf2_gpu_
    #define MAGMA_tally3_ZGEQRS_GPU magma_tally3_zgeqrs_gpu_
    #define MAGMA_tally3_ZGETRF_GPU magma_tally3_zgetrf_gpu_
    #define MAGMA_tally3_ZGETRS_GPU magma_tally3_zgetrs_gpu_
    #define MAGMA_tally3_ZLARFB_GPU magma_tally3_zlarfb_gpu_
    #define MAGMA_tally3_ZPOTRF_GPU magma_tally3_zpotrf_gpu_
    #define MAGMA_tally3_ZPOTRS_GPU magma_tally3_zpotrs_gpu_

#elif defined (NOCHANGE)

    #define MAGMA_tally3_ZGEBRD magma_tally3_zgebrd
    #define MAGMA_tally3_ZGEHRD magma_tally3_zgehrd
    #define MAGMA_tally3_ZGELQF magma_tally3_zgelqf
    #define MAGMA_tally3_ZGEQLF magma_tally3_zgeqlf
    #define MAGMA_tally3_ZGEQRF magma_tally3_zgeqrf
    #define MAGMA_tally3_ZGETRF magma_tally3_zgetrf
    #define MAGMA_tally3_ZLABRD magma_tally3_zlabrd
    #define MAGMA_tally3_ZLAHR2 magma_tally3_zlahr2
    #define MAGMA_tally3_ZLAHRU magma_tally3_zlahru
    #define MAGMA_tally3_ZPOTRF magma_tally3_zpotrf
    #define MAGMA_tally3_ZHETRD magma_tally3_zhetrd
    
    #define MAGMA_tally3_ZUNMQR_GPU magma_tally3_zunmqr_gpu
    #define MAGMA_tally3_ZGEQRF_GPU  magma_tally3_zgeqrf_gpu
    #define MAGMA_tally3_ZGEQRF2_GPU magma_tally3_zgeqrf2_gpu
    #define MAGMA_tally3_ZGEQRS_GPU magma_tally3_zgeqrs_gpu
    #define MAGMA_tally3_ZGETRF_GPU magma_tally3_zgetrf_gpu
    #define MAGMA_tally3_ZGETRS_GPU magma_tally3_zgetrs_gpu
    #define MAGMA_tally3_ZLARFB_GPU magma_tally3_zlarfb_gpu
    #define MAGMA_tally3_ZPOTRF_GPU magma_tally3_zpotrf_gpu
    #define MAGMA_tally3_ZPOTRS_GPU magma_tally3_zpotrs_gpu

#endif

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (simple interface)
 **/

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 function definitions / Data on CPU
*/
void MAGMA_tally3_ZGEBRD( magma_tally3_int_t *m, magma_tally3_int_t *n, double2 *A, magma_tally3_int_t *lda, double *d, double *e, double2 *tauq,  double2 *taup, double2 *work, magma_tally3_int_t *lwork, double2 *da, magma_tally3_int_t *info)
{ magma_tally3_zgebrd( *m, *n, A, *lda, d, e, tauq,  taup, work, *lwork, da, info); }

void MAGMA_tally3_ZGEHRD( magma_tally3_int_t *n, magma_tally3_int_t *ilo, magma_tally3_int_t *ihi, double2 *A, magma_tally3_int_t *lda, double2 *tau, double2 *work, magma_tally3_int_t *lwork, double2 *da, magma_tally3_int_t *info)
{ magma_tally3_zgehrd( *n, *ilo, *ihi, A, *lda, tau, work, lwork, da, info); }

void MAGMA_tally3_ZGELQF( magma_tally3_int_t *m, magma_tally3_int_t *n, double2 *A, magma_tally3_int_t *lda, double2 *tau, double2 *work, magma_tally3_int_t *lwork, magma_tally3_int_t *info)
{ magma_tally3_zgelqf( *m, *n, A, *lda, tau, work, *lwork, info); }

void MAGMA_tally3_ZGEQLF( magma_tally3_int_t *m, magma_tally3_int_t *n, double2 *A, magma_tally3_int_t *lda, double2 *tau, double2 *work, magma_tally3_int_t *lwork, magma_tally3_int_t *info)
{ magma_tally3_zgeqlf( *m, *n, A, *lda, tau, work, *lwork, info); }

void MAGMA_tally3_ZGEQRF( magma_tally3_int_t *m, magma_tally3_int_t *n, double2 *A, magma_tally3_int_t *lda, double2 *tau, double2 *work, magma_tally3_int_t *lwork, magma_tally3_int_t *info)
{ magma_tally3_zgeqrf( *m, *n, A, *lda, tau, work, *lwork, info); }

void MAGMA_tally3_ZGETRF( magma_tally3_int_t *m, magma_tally3_int_t *n, double2 *A, magma_tally3_int_t *lda, magma_tally3_int_t *ipiv, magma_tally3_int_t *info)
{ magma_tally3_zgetrf( *m, *n, A, *lda, ipiv, info); }

void MAGMA_tally3_ZLABRD( magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *nb, double2 *a, magma_tally3_int_t *lda, double *d, double *e, double2 *tauq, double2 *taup, double2 *x, magma_tally3_int_t *ldx, double2 *y, magma_tally3_int_t *ldy, double2 *da, magma_tally3_int_t *ldda, double2 *dx, magma_tally3_int_t *lddx, double2 *dy, magma_tally3_int_t *lddy)
{ magma_tally3_zlabrd( *m, *n, *nb, a, *lda, d, e, tauq, taup, x, *ldx, y, *ldy, da, *ldda, dx, *lddx, dy, *lddy); }

void MAGMA_tally3_ZLAHR2( magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *nb, double2 *da, double2 *dv, double2 *a, magma_tally3_int_t *lda, double2 *tau, double2 *t, magma_tally3_int_t *ldt, double2 *y, magma_tally3_int_t *ldy)
{ magma_tally3_zlahr2( *m, *n, *nb, da, dv, a, *lda, tau, t, *ldt, y, *ldy); }

void MAGMA_tally3_ZLAHRU( magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *nb, double2 *a, magma_tally3_int_t *lda, double2 *da, double2 *y, double2 *v, double2 *t, double2 *dwork)
{ magma_tally3_zlahru( *m, *n, *nb, a, *lda, da, y, v, t, dwork); }

void MAGMA_tally3_ZPOTRF( char *uplo, magma_tally3_int_t *n, double2 *A, magma_tally3_int_t *lda, magma_tally3_int_t *info)
{ magma_tally3_zpotrf( uplo[0], *n, A, *lda, info); }

void MAGMA_tally3_ZHETRD( char *uplo, magma_tally3_int_t *n, double2 *A, magma_tally3_int_t *lda, double *d, double *e, double2 *tau, double2 *work, magma_tally3_int_t *lwork, double2 *da, magma_tally3_int_t *info)
{ magma_tally3_zhetrd( uplo[0], *n, A, *lda, d, e, tau, work, *lwork, da, info); }


/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA_tally3 function definitions / Data on GPU
*/
void MAGMA_tally3_ZUNMQR_GPU(char *side, char *trans, magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *k, double2 *a, magma_tally3_int_t *lda, double2 *tau, double2 *c, magma_tally3_int_t *ldc, double2 *work, magma_tally3_int_t *lwork, double2 *td, magma_tally3_int_t *nb, magma_tally3_int_t *info)
{ magma_tally3_zunmqr_gpu(side[0], trans[0], *m, *n, *k, a, *lda, tau, c, *ldc, work, *lwork, td, *nb, info); }

void MAGMA_tally3_ZGEQRF2_GPU( magma_tally3_int_t *m, magma_tally3_int_t *n, double2 *A, magma_tally3_int_t *lda, double2 *tau, magma_tally3_int_t *info)
{ magma_tally3_zgeqrf2_gpu( *m, *n, A, *lda, tau, info); }

void MAGMA_tally3_ZGEQRF_GPU(magma_tally3_int_t *m, magma_tally3_int_t *n, double2 *A, magma_tally3_int_t *lda, double2 *tau, double2 *dwork, magma_tally3_int_t *info)
{ magma_tally3_zgeqrf_gpu(*m, *n, A, *lda, tau, dwork, info); }

void MAGMA_tally3_ZGEQRS_GPU( magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *nrhs, double2 *A, magma_tally3_int_t *lda, double2 *tau, double2 *td, double2 *c, magma_tally3_int_t *ldc, double2 *work, magma_tally3_int_t *lwork, magma_tally3_int_t *info)
{ magma_tally3_zgeqrs_gpu( *m, *n, *nrhs, A, *lda, tau, td, c, *ldc, work, *lwork, info); }

void MAGMA_tally3_ZGETRF_GPU( magma_tally3_int_t *m, magma_tally3_int_t *n, double2 *A, magma_tally3_int_t *lda, magma_tally3_int_t *ipiv, magma_tally3_int_t *info)
{ magma_tally3_zgetrf_gpu( *m, *n, A, *lda, ipiv, info); }

void MAGMA_tally3_ZGETRS_GPU( char *trans, magma_tally3_int_t *n, magma_tally3_int_t *nrhs, double2 *A, magma_tally3_int_t *lda, magma_tally3_int_t *ipiv, double2 *b, magma_tally3_int_t *ldb, magma_tally3_int_t *info)
{ magma_tally3_zgetrs_gpu( trans[0], *n, *nrhs, A, *lda, ipiv, b, *ldb, info ); }

void MAGMA_tally3_ZLARFB_GPU( char *side, char *trans, char *direct, char *storev, magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *k, double2 *dv, magma_tally3_int_t *ldv, double2 *dt, magma_tally3_int_t *ldt, double2 *dc, magma_tally3_int_t *ldc, double2 *dowrk, magma_tally3_int_t *ldwork)
{ magma_tally3_zlarfb_gpu( side[0], trans[0], direct[0], storev[0], *m, *n, *k, dv, *ldv, dt, *ldt, dc, *ldc, dowrk, *ldwork); }

void MAGMA_tally3_ZPOTRF_GPU( char *uplo,  magma_tally3_int_t *n, double2 *A, magma_tally3_int_t *lda, magma_tally3_int_t *info)
{ magma_tally3_zpotrf_gpu( uplo[0], *n, A, *lda, info); }

void MAGMA_tally3_ZPOTRS_GPU( char *uplo,  magma_tally3_int_t *n, magma_tally3_int_t *nrhs, double2 *A, magma_tally3_int_t *lda, double2 *b, magma_tally3_int_t *ldb, magma_tally3_int_t *info)
{ magma_tally3_zpotrs_gpu( uplo[0], *n, *nrhs, A, *lda, b, *ldb, info); }

#ifdef __cplusplus
}
#endif
