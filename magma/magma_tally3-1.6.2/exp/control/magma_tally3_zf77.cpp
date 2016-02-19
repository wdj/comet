/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#if defined(__GNUC__)
#include <stdint.h>
#endif /* __GNUC__ */
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include "magma_tally3.h"

/* 
 * typedef comming from fortran.h file provided in $CUDADIR/src directory
 * it will probably change with future release of cublas when they will use 64bits address
 */
typedef size_t devptr_t;

#define PRECISION_z

#ifdef PGI_FORTRAN
#define DEVPTR(__ptr) ((cuDoubleComplex*)(__ptr))
#else
#define DEVPTR(__ptr) ((cuDoubleComplex*)(uintptr_t)(*(__ptr)))
#endif


#ifndef MAGMA_tally3_FORTRAN_NAME
#if defined(ADD_)
#define MAGMA_tally3_FORTRAN_NAME(lcname, UCNAME)  magma_tally3f_##lcname##_
#elif defined(NOCHANGE)
#define MAGMA_tally3_FORTRAN_NAME(lcname, UCNAME)  magma_tally3f_##lcname
#elif defined(UPCASE)
#define MAGMA_tally3_FORTRAN_NAME(lcname, UCNAME)  MAGMA_tally3F_##UCNAME
#endif
#endif

#ifndef MAGMA_tally3_GPU_FORTRAN_NAME
#if defined(ADD_)
#define MAGMA_tally3_GPU_FORTRAN_NAME(lcname, UCNAME)  magma_tally3f_##lcname##_gpu_
#elif defined(NOCHANGE)
#define MAGMA_tally3_GPU_FORTRAN_NAME(lcname, UCNAME)  magma_tally3f_##lcname##_gpu
#elif defined(UPCASE)
#define MAGMA_tally3_GPU_FORTRAN_NAME(lcname, UCNAME)  MAGMA_tally3F_##UCNAME##_GPU
#endif
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 function definitions / Data on CPU
*/
#define MAGMA_tally3F_ZGEBRD  MAGMA_tally3_FORTRAN_NAME(zgebrd,  ZGEBRD ) 
#define MAGMA_tally3F_ZGEHRD2 MAGMA_tally3_FORTRAN_NAME(zgehrd2, ZGEHRD2)
#define MAGMA_tally3F_ZGEHRD  MAGMA_tally3_FORTRAN_NAME(zgehrd,  ZGEHRD )
#define MAGMA_tally3F_ZGELQF  MAGMA_tally3_FORTRAN_NAME(zgelqf,  ZGELQF )
#define MAGMA_tally3F_ZGEQLF  MAGMA_tally3_FORTRAN_NAME(zgeqlf,  ZGEQLF )
#define MAGMA_tally3F_ZGEQRF  MAGMA_tally3_FORTRAN_NAME(zgeqrf,  ZGEQRF )
#define MAGMA_tally3F_ZGESV   MAGMA_tally3_FORTRAN_NAME(zgesv,   ZGESV  )
#define MAGMA_tally3F_ZGETRF  MAGMA_tally3_FORTRAN_NAME(zgetrf,  ZGETRF )
#define MAGMA_tally3F_ZLATRD  MAGMA_tally3_FORTRAN_NAME(zlatrd,  ZLATRD )
#define MAGMA_tally3F_ZLAHR2  MAGMA_tally3_FORTRAN_NAME(zlahr2,  ZLAHR2 )
#define MAGMA_tally3F_ZLAHRU  MAGMA_tally3_FORTRAN_NAME(zlahru,  ZLAHRU )
#define MAGMA_tally3F_ZPOSV   MAGMA_tally3_FORTRAN_NAME(zposv,   ZPOSV  )
#define MAGMA_tally3F_ZPOTRF  MAGMA_tally3_FORTRAN_NAME(zpotrf,  ZPOTRF )
#define MAGMA_tally3F_ZHETRD  MAGMA_tally3_FORTRAN_NAME(zhetrd,  ZHETRD )
#define MAGMA_tally3F_ZUNGQR  MAGMA_tally3_FORTRAN_NAME(zungqr,  ZUNGQR )
#define MAGMA_tally3F_ZUNMQR  MAGMA_tally3_FORTRAN_NAME(zunmqr,  ZUNMQR )
#define MAGMA_tally3F_ZUNMTR  MAGMA_tally3_FORTRAN_NAME(zunmtr,  ZUNMTR )
#define MAGMA_tally3F_ZUNGHR  MAGMA_tally3_FORTRAN_NAME(zunghr,  ZUNGHR )
#define MAGMA_tally3F_ZGEEV   MAGMA_tally3_FORTRAN_NAME(zgeev,   ZGEEV  )
#define MAGMA_tally3F_ZGESVD  MAGMA_tally3_FORTRAN_NAME(zgesvd,  ZGESVD )
#define MAGMA_tally3F_ZHEEVD  MAGMA_tally3_FORTRAN_NAME(zheevd,  ZHEEVD )
#define MAGMA_tally3F_ZHEGVD  MAGMA_tally3_FORTRAN_NAME(zhegvd,  ZHEGVD )

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA_tally3 function definitions / Data on GPU
*/
#define MAGMA_tally3F_ZGELS_GPU   MAGMA_tally3_GPU_FORTRAN_NAME(zgels,   ZGELS  )
#define MAGMA_tally3F_ZGEQRF_GPU  MAGMA_tally3_GPU_FORTRAN_NAME(zgeqrf,  ZGEQRF ) 
#define MAGMA_tally3F_ZGEQRF2_GPU MAGMA_tally3_GPU_FORTRAN_NAME(zgeqrf2, ZGEQRF2)
#define MAGMA_tally3F_ZGEQRF3_GPU MAGMA_tally3_GPU_FORTRAN_NAME(zgeqrf3, ZGEQRF3)
#define MAGMA_tally3F_ZGEQRS_GPU  MAGMA_tally3_GPU_FORTRAN_NAME(zgeqrs,  ZGEQRS ) 
#define MAGMA_tally3F_ZGEQRS3_GPU MAGMA_tally3_GPU_FORTRAN_NAME(zgeqrs3, ZGEQRS3) 
#define MAGMA_tally3F_ZGESSM_GPU  MAGMA_tally3_GPU_FORTRAN_NAME(zgessm,  ZGESSM ) 
#define MAGMA_tally3F_ZGESV_GPU   MAGMA_tally3_GPU_FORTRAN_NAME(zgesv,   ZGESV  )  
#define MAGMA_tally3F_ZGETRF_INCPIV_GPU  MAGMA_tally3_GPU_FORTRAN_NAME(zgetrf_incpiv,  ZGETRF_INCPIV ) 
#define MAGMA_tally3F_ZGETRF_GPU  MAGMA_tally3_GPU_FORTRAN_NAME(zgetrf,  ZGETRF ) 
#define MAGMA_tally3F_ZGETRS_GPU  MAGMA_tally3_GPU_FORTRAN_NAME(zgetrs,  ZGETRS ) 
#define MAGMA_tally3F_ZLABRD_GPU  MAGMA_tally3_GPU_FORTRAN_NAME(zlabrd,  ZLABRD ) 
#define MAGMA_tally3F_ZLARFB_GPU  MAGMA_tally3_GPU_FORTRAN_NAME(zlarfb,  ZLARFB ) 
#define MAGMA_tally3F_ZPOSV_GPU   MAGMA_tally3_GPU_FORTRAN_NAME(zposv,   ZPOSV  )  
#define MAGMA_tally3F_ZPOTRF_GPU  MAGMA_tally3_GPU_FORTRAN_NAME(zpotrf,  ZPOTRF ) 
#define MAGMA_tally3F_ZPOTRS_GPU  MAGMA_tally3_GPU_FORTRAN_NAME(zpotrs,  ZPOTRS ) 
#define MAGMA_tally3F_ZSSSSM_GPU  MAGMA_tally3_GPU_FORTRAN_NAME(zssssm,  ZSSSSM ) 
#define MAGMA_tally3F_ZTSTRF_GPU  MAGMA_tally3_GPU_FORTRAN_NAME(ztstrf,  ZTSTRF ) 
#define MAGMA_tally3F_ZUNGQR_GPU  MAGMA_tally3_GPU_FORTRAN_NAME(zungqr,  ZUNGQR ) 
#define MAGMA_tally3F_ZUNMQR_GPU  MAGMA_tally3_GPU_FORTRAN_NAME(zunmqr,  ZUNMQR ) 

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (simple interface)
 **/

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 function definitions / Data on CPU
*/
void MAGMA_tally3F_ZGEBRD( magma_tally3_int_t *m, magma_tally3_int_t *n, cuDoubleComplex *A, 
                    magma_tally3_int_t *lda, double *d, double *e,
                    cuDoubleComplex *tauq, cuDoubleComplex *taup, 
                    cuDoubleComplex *work, magma_tally3_int_t *lwork, magma_tally3_int_t *info)
{
    magma_tally3_zgebrd( *m, *n, A, 
                  *lda, d, e,
                  tauq, taup, 
                  work, *lwork, info);
}
    
void MAGMA_tally3F_ZGEHRD2(magma_tally3_int_t *n, magma_tally3_int_t *ilo, magma_tally3_int_t *ihi,
                    cuDoubleComplex *A, magma_tally3_int_t *lda, cuDoubleComplex *tau, 
                    cuDoubleComplex *work, magma_tally3_int_t *lwork, magma_tally3_int_t *info)
{
    magma_tally3_zgehrd2(*n, *ilo, *ihi,
                  A, *lda, tau, 
                  work, lwork, info);
}
    
void MAGMA_tally3F_ZGEHRD( magma_tally3_int_t *n, magma_tally3_int_t *ilo, magma_tally3_int_t *ihi,
            cuDoubleComplex *A, magma_tally3_int_t *lda, cuDoubleComplex *tau,
            cuDoubleComplex *work, magma_tally3_int_t *lwork,
            cuDoubleComplex *d_T, magma_tally3_int_t *info)
{
  magma_tally3_zgehrd( *n, *ilo, *ihi,
        A, *lda, tau,
        work, *lwork,
        d_T, info);
}

void MAGMA_tally3F_ZGELQF( magma_tally3_int_t *m, magma_tally3_int_t *n, 
                    cuDoubleComplex *A,    magma_tally3_int_t *lda,   cuDoubleComplex *tau, 
                    cuDoubleComplex *work, magma_tally3_int_t *lwork, magma_tally3_int_t *info)
{
    magma_tally3_zgelqf( *m, *n, 
                  A,    *lda,   tau, 
                  work, *lwork, info);
}

void MAGMA_tally3F_ZGEQLF( magma_tally3_int_t *m, magma_tally3_int_t *n, 
                    cuDoubleComplex *A,    magma_tally3_int_t *lda,   cuDoubleComplex *tau, 
                    cuDoubleComplex *work, magma_tally3_int_t *lwork, magma_tally3_int_t *info)
{
    magma_tally3_zgeqlf( *m, *n, 
                  A,    *lda,   tau, 
                  work, *lwork, info);
}

void MAGMA_tally3F_ZGEQRF( magma_tally3_int_t *m, magma_tally3_int_t *n, cuDoubleComplex *A, 
                    magma_tally3_int_t *lda, cuDoubleComplex *tau, cuDoubleComplex *work, 
                    magma_tally3_int_t *lwork, magma_tally3_int_t *info)
{
    magma_tally3_zgeqrf( *m, *n, A, 
                  *lda, tau, work, 
                  *lwork, info);
}

void MAGMA_tally3F_ZGESV ( magma_tally3_int_t *n, magma_tally3_int_t *nrhs,
                    cuDoubleComplex *A, magma_tally3_int_t *lda, magma_tally3_int_t *ipiv,
                    cuDoubleComplex *B, magma_tally3_int_t *ldb, magma_tally3_int_t *info)
{
    magma_tally3_zgesv(  *n, *nrhs,
                  A, *lda, ipiv,
                  B, *ldb,
                  info);
}
    
void MAGMA_tally3F_ZGETRF( magma_tally3_int_t *m, magma_tally3_int_t *n, cuDoubleComplex *A, 
                    magma_tally3_int_t *lda, magma_tally3_int_t *ipiv, 
                    magma_tally3_int_t *info)
{
    magma_tally3_zgetrf( *m, *n, A, 
                  *lda, ipiv, 
                  info);
}

// void MAGMA_tally3F_ZLATRD( char *uplo, magma_tally3_int_t *n, magma_tally3_int_t *nb, cuDoubleComplex *a, 
//                     magma_tally3_int_t *lda, double *e, cuDoubleComplex *tau, 
//                     cuDoubleComplex *w, magma_tally3_int_t *ldw,
//                     cuDoubleComplex *da, magma_tally3_int_t *ldda, 
//                     cuDoubleComplex *dw, magma_tally3_int_t *lddw)
// {
//     magma_tally3_zlatrd( uplo[0], *n, *nb, a, 
//                   *lda, e, tau, 
//                   w, *ldw,
//                   da, *ldda, 
//                   dw, *lddw);
// }

  /* This has nothing to do here, it should be a GPU function */
// void MAGMA_tally3F_ZLAHR2( magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *nb, 
//                     cuDoubleComplex *da, cuDoubleComplex *dv, cuDoubleComplex *a, 
//                     magma_tally3_int_t *lda, cuDoubleComplex *tau, cuDoubleComplex *t, 
//                     magma_tally3_int_t *ldt, cuDoubleComplex *y, magma_tally3_int_t *ldy)
// {
//     magma_tally3_zlahr2( *m, *n, *nb, 
//                   da, dv, a, 
//                   *lda, tau, t, 
//                   *ldt, y, *ldy);
// }

// void MAGMA_tally3F_ZLAHRU( magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *nb, 
//                     cuDoubleComplex *a, magma_tally3_int_t *lda, 
//                     cuDoubleComplex *da, cuDoubleComplex *y, 
//                     cuDoubleComplex *v, cuDoubleComplex *t, 
//                     cuDoubleComplex *dwork)
// {
//     magma_tally3_zlahru( *m, *n, *nb, 
//                   a, *lda, 
//                   da, y, 
//                   v, t, 
//                   dwork);
// }

void MAGMA_tally3F_ZPOSV(  char *uplo, magma_tally3_int_t *n, magma_tally3_int_t *nrhs,
                    cuDoubleComplex *A, magma_tally3_int_t *lda,
                    cuDoubleComplex *B, magma_tally3_int_t *ldb, magma_tally3_int_t *info)
{
    magma_tally3_zposv(  uplo[0], *n, *nrhs,
                  A, *lda,
                  B, *ldb, info);
}

void MAGMA_tally3F_ZPOTRF( char *uplo, magma_tally3_int_t *n, cuDoubleComplex *A, 
                    magma_tally3_int_t *lda, magma_tally3_int_t *info)
{
    magma_tally3_zpotrf( uplo[0], *n, A, 
                  *lda, info);
}

void MAGMA_tally3F_ZHETRD( char *uplo, magma_tally3_int_t *n, cuDoubleComplex *A, 
                    magma_tally3_int_t *lda, double *d, double *e, 
                    cuDoubleComplex *tau, cuDoubleComplex *work, magma_tally3_int_t *lwork, 
                    magma_tally3_int_t *info)
{
    magma_tally3_zhetrd( uplo[0], *n, A, 
                  *lda, d, e, 
                  tau, work, *lwork, 
                  info);
}

// void MAGMA_tally3F_ZUNGQR( magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *k,
//                     cuDoubleComplex *a, magma_tally3_int_t *lda,
//                     cuDoubleComplex *tau, cuDoubleComplex *dwork,
//                     magma_tally3_int_t *nb, magma_tally3_int_t *info )
// {
//     magma_tally3_zungqr( *m, *n, *k,
//                   a, *lda,
//                   tau, dwork,
//                   *nb, info );
// }

void MAGMA_tally3F_ZUNMQR( char *side, char *trans, 
                    magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *k, 
                    cuDoubleComplex *a, magma_tally3_int_t *lda, cuDoubleComplex *tau, 
                    cuDoubleComplex *c, magma_tally3_int_t *ldc, 
                    cuDoubleComplex *work, magma_tally3_int_t *lwork, magma_tally3_int_t *info)
{
    magma_tally3_zunmqr( side[0], trans[0], 
                  *m, *n, *k, 
                  a, *lda, tau, 
                  c, *ldc, 
                  work, *lwork, info);
}

void MAGMA_tally3F_ZUNMTR( char *side, char *uplo, char *trans,
                    magma_tally3_int_t *m, magma_tally3_int_t *n,
                    cuDoubleComplex *a,    magma_tally3_int_t *lda,
                    cuDoubleComplex *tau,
                    cuDoubleComplex *c,    magma_tally3_int_t *ldc,
                    cuDoubleComplex *work, magma_tally3_int_t *lwork,
                    magma_tally3_int_t *info)
{
    magma_tally3_zunmtr( side[0], uplo[0], trans[0],
                  *m, *n,
                  a,    *lda,
                  tau,
                  c,    *ldc,
                  work, *lwork,
                  info);
}

// void MAGMA_tally3F_ZUNGHR( magma_tally3_int_t *n, magma_tally3_int_t *ilo, magma_tally3_int_t *ihi,
//                     cuDoubleComplex *a, magma_tally3_int_t *lda,
//                     cuDoubleComplex *tau,
//                     cuDoubleComplex *dT, magma_tally3_int_t *nb,
//                     magma_tally3_int_t *info)
// {
//     magma_tally3_zunghr( *n, *ilo, *ihi,
//                   a, *lda,
//                   tau,
//                   dT, *nb,
//                   info);
// }

#if defined(PRECISION_z) || defined(PRECISION_c)
void MAGMA_tally3F_ZGEEV( char *jobvl, char *jobvr, magma_tally3_int_t *n,
                   cuDoubleComplex *a, magma_tally3_int_t *lda,
                   cuDoubleComplex *w,
                   cuDoubleComplex *vl, magma_tally3_int_t *ldvl,
                   cuDoubleComplex *vr, magma_tally3_int_t *ldvr,
                   cuDoubleComplex *work, magma_tally3_int_t *lwork,
                   double *rwork, magma_tally3_int_t *info)
{
    magma_tally3_zgeev( jobvl[0], jobvr[0], *n,
                 a, *lda,
                 w,
                 vl, *ldvl,
                 vr, *ldvr,
                 work, *lwork,
                 rwork, info);
}

void MAGMA_tally3F_ZGESVD( char *jobu, char *jobvt, magma_tally3_int_t *m, magma_tally3_int_t *n,
                    cuDoubleComplex *a,    magma_tally3_int_t *lda, double *s, 
                    cuDoubleComplex *u,    magma_tally3_int_t *ldu, 
                    cuDoubleComplex *vt,   magma_tally3_int_t *ldvt,
                    cuDoubleComplex *work, magma_tally3_int_t *lwork,
                    double *rwork, magma_tally3_int_t *info )
{
    magma_tally3_zgesvd( jobu[0], jobvt[0], *m, *n,
                  a,    *lda, s, 
                  u,    *ldu, 
                  vt,   *ldvt,
                  work, *lwork,
                  rwork, info );
}
    
void MAGMA_tally3F_ZHEEVD( char *jobz, char *uplo, magma_tally3_int_t *n,
                    cuDoubleComplex *a,     magma_tally3_int_t *lda, double *w,
                    cuDoubleComplex *work,  magma_tally3_int_t *lwork,
                    double          *rwork, magma_tally3_int_t *lrwork,
                    magma_tally3_int_t *iwork, magma_tally3_int_t *liwork, magma_tally3_int_t *info)
{
    magma_tally3_zheevd( jobz[0], uplo[0], *n,
                  a, *lda, w,
                  work, *lwork,
                  rwork, *lrwork,
                  iwork, *liwork, info);
}

void MAGMA_tally3F_ZHEGVD(magma_tally3_int_t *itype, char *jobz, char *uplo, magma_tally3_int_t *n,
           cuDoubleComplex *a, magma_tally3_int_t *lda, 
           cuDoubleComplex *b, magma_tally3_int_t *ldb,
           double *w, cuDoubleComplex *work, magma_tally3_int_t *lwork,
           double *rwork, magma_tally3_int_t *lrwork,
           magma_tally3_int_t *iwork, magma_tally3_int_t *liwork, magma_tally3_int_t *info)
{
  magma_tally3_zhegvd( *itype, jobz[0], uplo[0], *n,
        a, *lda, b, *ldb,
        w, work, *lwork,
        rwork, *lrwork,
        iwork, *liwork, info);
}
    
#else
void MAGMA_tally3F_ZGEEV( char *jobvl, char *jobvr, magma_tally3_int_t *n,
                   cuDoubleComplex *a,    magma_tally3_int_t *lda,
                   cuDoubleComplex *wr, cuDoubleComplex *wi,
                   cuDoubleComplex *vl,   magma_tally3_int_t *ldvl,
                   cuDoubleComplex *vr,   magma_tally3_int_t *ldvr,
                   cuDoubleComplex *work, magma_tally3_int_t *lwork,
                   magma_tally3_int_t *info)
{
    magma_tally3_zgeev( jobvl[0], jobvr[0], *n,
                 a,    *lda,
                 wr, wi,
                 vl,   *ldvl,
                 vr,   *ldvr,
                 work, *lwork,
                 info);
}

void MAGMA_tally3F_ZGESVD( char *jobu, char *jobvt, magma_tally3_int_t *m, magma_tally3_int_t *n,
                    cuDoubleComplex *a,    magma_tally3_int_t *lda, double *s,
                    cuDoubleComplex *u,    magma_tally3_int_t *ldu, 
                    cuDoubleComplex *vt,   magma_tally3_int_t *ldvt,
                    cuDoubleComplex *work, magma_tally3_int_t *lwork,
                    magma_tally3_int_t *info )
{
    magma_tally3_zgesvd( jobu[0], jobvt[0], *m, *n,
                  a,    *lda, s, 
                  u,    *ldu, 
                  vt,   *ldvt,
                  work, *lwork,
                  info );
}

void MAGMA_tally3F_ZHEEVD( char *jobz, char *uplo, magma_tally3_int_t *n,
                    cuDoubleComplex *a, magma_tally3_int_t *lda, double *w,
                    cuDoubleComplex *work, magma_tally3_int_t *lwork,
                    magma_tally3_int_t *iwork, magma_tally3_int_t *liwork, magma_tally3_int_t *info)
{
    magma_tally3_zheevd( jobz[0], uplo[0], *n,
                  a, *lda, w,
                  work, *lwork,
                  iwork, *liwork, info);
}

void MAGMA_tally3F_ZHEGVD(magma_tally3_int_t *itype, char *jobz, char *uplo, magma_tally3_int_t *n,
                   cuDoubleComplex *a, magma_tally3_int_t *lda,
                   cuDoubleComplex *b, magma_tally3_int_t *ldb,
                   double *w, cuDoubleComplex *work, magma_tally3_int_t *lwork,
                   magma_tally3_int_t *iwork, magma_tally3_int_t *liwork, magma_tally3_int_t *info)
{
  magma_tally3_zhegvd( *itype, jobz[0], uplo[0], *n,
                a, *lda, b, *ldb,
        w, work, *lwork,
                iwork, *liwork, info);
}


#endif

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA_tally3 function definitions / Data on GPU
*/
void MAGMA_tally3F_ZGELS_GPU(  char *trans, magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *nrhs,
                        devptr_t *dA,    magma_tally3_int_t *ldda, 
                        devptr_t *dB,    magma_tally3_int_t *lddb, 
                        cuDoubleComplex *hwork, magma_tally3_int_t *lwork, 
                        magma_tally3_int_t *info)
{
    magma_tally3_zgels_gpu(  trans[0], *m, *n, *nrhs, 
                      DEVPTR(dA),    *ldda,  
                      DEVPTR(dB),    *lddb,  
                      hwork, *lwork,  info);
}

void MAGMA_tally3F_ZGEQRF_GPU( magma_tally3_int_t *m, magma_tally3_int_t *n, 
                        devptr_t *dA,  magma_tally3_int_t *ldda, 
                        cuDoubleComplex *tau, devptr_t *dT, 
                        magma_tally3_int_t *info)
{
    magma_tally3_zgeqrf_gpu( *m, *n,  
                      DEVPTR(dA),  *ldda,  
                      tau, 
                      DEVPTR(dT),  info);
}

void MAGMA_tally3F_ZGEQRF2_GPU(magma_tally3_int_t *m, magma_tally3_int_t *n, 
                        devptr_t *dA,  magma_tally3_int_t *ldda, 
                        cuDoubleComplex *tau, magma_tally3_int_t *info)
{
    magma_tally3_zgeqrf2_gpu(*m, *n,  
                      DEVPTR(dA),  *ldda,  
                      tau, info); 
}

void MAGMA_tally3F_ZGEQRF3_GPU(magma_tally3_int_t *m, magma_tally3_int_t *n, 
                        devptr_t *dA,  magma_tally3_int_t *ldda, 
                        cuDoubleComplex *tau, devptr_t *dT,
                        magma_tally3_int_t *info)
{
    magma_tally3_zgeqrf3_gpu(*m, *n,  
                      DEVPTR(dA),  *ldda,  
                      tau, DEVPTR(dT), info); 
}

void MAGMA_tally3F_ZGEQRS_GPU( magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *nrhs, 
                        devptr_t *dA,     magma_tally3_int_t *ldda, 
                        cuDoubleComplex *tau,   devptr_t *dT,
                        devptr_t *dB,    magma_tally3_int_t *lddb,
                        cuDoubleComplex *hwork, magma_tally3_int_t *lhwork, 
                        magma_tally3_int_t *info)
{
    magma_tally3_zgeqrs_gpu( *m, *n, *nrhs,  
                      DEVPTR(dA),     *ldda,  
                      tau,
                      DEVPTR(dT), 
                      DEVPTR(dB),    *lddb, 
                      hwork, *lhwork,  info);
}

void MAGMA_tally3F_ZGEQRS3_GPU(magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *nrhs, 
                        devptr_t *dA,     magma_tally3_int_t *ldda, 
                        cuDoubleComplex *tau,   devptr_t *dT,
                        devptr_t *dB,    magma_tally3_int_t *lddb,
                        cuDoubleComplex *hwork, magma_tally3_int_t *lhwork, 
                        magma_tally3_int_t *info)
{
    magma_tally3_zgeqrs3_gpu(*m, *n, *nrhs,  
                      DEVPTR(dA),     *ldda,  
                      tau,
                      DEVPTR(dT), 
                      DEVPTR(dB),    *lddb, 
                      hwork, *lhwork,  info);
}

void MAGMA_tally3F_ZGESSM_GPU( char *storev, magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *k, magma_tally3_int_t *ib, 
                        magma_tally3_int_t *ipiv, 
                        devptr_t *dL1, magma_tally3_int_t *lddl1, 
                        devptr_t *dL,  magma_tally3_int_t *lddl, 
                        devptr_t *dA,  magma_tally3_int_t *ldda, 
                        magma_tally3_int_t *info)
{
    magma_tally3_zgessm_gpu( storev[0], *m, *n, *k, *ib, ipiv,  
                      DEVPTR(dL1), *lddl1,  
                      DEVPTR(dL),  *lddl,  
                      DEVPTR(dA),  *ldda,  info);
}

void MAGMA_tally3F_ZGESV_GPU(  magma_tally3_int_t *n, magma_tally3_int_t *nrhs, 
                        devptr_t *dA, magma_tally3_int_t *ldda, magma_tally3_int_t *ipiv, 
                        devptr_t *dB, magma_tally3_int_t *lddb, magma_tally3_int_t *info)
{
    magma_tally3_zgesv_gpu(  *n, *nrhs,  
                      DEVPTR(dA), *ldda, ipiv,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMA_tally3F_ZGETRF_GPU( magma_tally3_int_t *m, magma_tally3_int_t *n, 
                        devptr_t *dA, magma_tally3_int_t *ldda, 
                        magma_tally3_int_t *ipiv, magma_tally3_int_t *info)
{
    magma_tally3_zgetrf_gpu( *m, *n,  
                      DEVPTR(dA), *ldda, ipiv, info);
}

void MAGMA_tally3F_ZGETRS_GPU( char *trans, magma_tally3_int_t *n, magma_tally3_int_t *nrhs, 
                        devptr_t *dA, magma_tally3_int_t *ldda, magma_tally3_int_t *ipiv, 
                        devptr_t *dB, magma_tally3_int_t *lddb, magma_tally3_int_t *info)
{
    magma_tally3_zgetrs_gpu( trans[0], *n, *nrhs,  
                      DEVPTR(dA), *ldda, ipiv,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMA_tally3F_ZLABRD_GPU( magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *nb, 
                        cuDoubleComplex *a, magma_tally3_int_t *lda, devptr_t *da, magma_tally3_int_t *ldda,
                        double *d, double *e, cuDoubleComplex *tauq, cuDoubleComplex *taup,  
                        cuDoubleComplex *x, magma_tally3_int_t *ldx, devptr_t *dx, magma_tally3_int_t *lddx, 
                        cuDoubleComplex *y, magma_tally3_int_t *ldy, devptr_t *dy, magma_tally3_int_t *lddy)
{
    magma_tally3_zlabrd_gpu( *m, *n, *nb,  
                      a, *lda, DEVPTR(da), *ldda, 
                      d, e, tauq, taup,   
                      x, *ldx, DEVPTR(dx), *lddx,  
                      y, *ldy, DEVPTR(dy), *lddy);
}

void MAGMA_tally3F_ZLARFB_GPU( char *side, char *trans, char *direct, char *storev, 
                        magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *k,
                        devptr_t *dv, magma_tally3_int_t *ldv, devptr_t *dt,    magma_tally3_int_t *ldt, 
                        devptr_t *dc, magma_tally3_int_t *ldc, devptr_t *dowrk, magma_tally3_int_t *ldwork )
{
    magma_tally3_zlarfb_gpu( side[0], trans[0], direct[0], storev[0],  *m, *n, *k, 
                      DEVPTR(dv), *ldv, DEVPTR(dt),    *ldt,  
                      DEVPTR(dc), *ldc, DEVPTR(dowrk), *ldwork);
}

void MAGMA_tally3F_ZPOSV_GPU(  char *uplo, magma_tally3_int_t *n, magma_tally3_int_t *nrhs, 
                        devptr_t *dA, magma_tally3_int_t *ldda, 
                        devptr_t *dB, magma_tally3_int_t *lddb, magma_tally3_int_t *info)
{
    magma_tally3_zposv_gpu(  uplo[0], *n, *nrhs,  
                      DEVPTR(dA), *ldda,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMA_tally3F_ZPOTRF_GPU( char *uplo,  magma_tally3_int_t *n, 
                        devptr_t *dA, magma_tally3_int_t *ldda, magma_tally3_int_t *info)
{
    magma_tally3_zpotrf_gpu( uplo[0],  *n,  
                      DEVPTR(dA), *ldda, info); }

void MAGMA_tally3F_ZPOTRS_GPU( char *uplo,  magma_tally3_int_t *n, magma_tally3_int_t *nrhs, 
                        devptr_t *dA, magma_tally3_int_t *ldda, 
                        devptr_t *dB, magma_tally3_int_t *lddb, magma_tally3_int_t *info)
{
    magma_tally3_zpotrs_gpu( uplo[0],  *n, *nrhs,  
                      DEVPTR(dA), *ldda,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMA_tally3F_ZSSSSM_GPU( char *storev, magma_tally3_int_t *m1, magma_tally3_int_t *n1, 
                        magma_tally3_int_t *m2, magma_tally3_int_t *n2, magma_tally3_int_t *k, magma_tally3_int_t *ib, 
                        devptr_t *dA1, magma_tally3_int_t *ldda1, 
                        devptr_t *dA2, magma_tally3_int_t *ldda2, 
                        devptr_t *dL1, magma_tally3_int_t *lddl1, 
                        devptr_t *dL2, magma_tally3_int_t *lddl2,
                        magma_tally3_int_t *IPIV, magma_tally3_int_t *info)
{
    magma_tally3_zssssm_gpu( storev[0], *m1, *n1,  *m2, *n2, *k, *ib,  
                      DEVPTR(dA1), *ldda1,  
                      DEVPTR(dA2), *ldda2,  
                      DEVPTR(dL1), *lddl1,  
                      DEVPTR(dL2), *lddl2,
                      IPIV, info);
}

void MAGMA_tally3F_ZUNGQR_GPU( magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *k, 
                        devptr_t *da, magma_tally3_int_t *ldda, 
                        cuDoubleComplex *tau, devptr_t *dwork, 
                        magma_tally3_int_t *nb, magma_tally3_int_t *info )
{
    magma_tally3_zungqr_gpu( *m, *n, *k,  
                      DEVPTR(da), *ldda, tau, 
                      DEVPTR(dwork), *nb, info );
}

void MAGMA_tally3F_ZUNMQR_GPU( char *side, char *trans, 
                        magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *k,
                        devptr_t *a,    magma_tally3_int_t *lda, cuDoubleComplex *tau, 
                        devptr_t *c,    magma_tally3_int_t *ldc,
                        devptr_t *work, magma_tally3_int_t *lwork, 
                        devptr_t *td,   magma_tally3_int_t *nb, magma_tally3_int_t *info)
{
    magma_tally3_zunmqr_gpu( side[0], trans[0], *m, *n, *k, 
                      DEVPTR(a),    *lda, tau,  
                      DEVPTR(c),    *ldc, 
                      DEVPTR(work), *lwork,  
                      DEVPTR(td),   *nb, info);
}

#ifdef __cplusplus
}
#endif
