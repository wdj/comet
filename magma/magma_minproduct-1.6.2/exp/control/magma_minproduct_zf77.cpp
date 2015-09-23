/*
    -- MAGMA_minproduct (version 1.6.1) --
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
#include "magma_minproduct.h"

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


#ifndef MAGMA_minproduct_FORTRAN_NAME
#if defined(ADD_)
#define MAGMA_minproduct_FORTRAN_NAME(lcname, UCNAME)  magma_minproductf_##lcname##_
#elif defined(NOCHANGE)
#define MAGMA_minproduct_FORTRAN_NAME(lcname, UCNAME)  magma_minproductf_##lcname
#elif defined(UPCASE)
#define MAGMA_minproduct_FORTRAN_NAME(lcname, UCNAME)  MAGMA_minproductF_##UCNAME
#endif
#endif

#ifndef MAGMA_minproduct_GPU_FORTRAN_NAME
#if defined(ADD_)
#define MAGMA_minproduct_GPU_FORTRAN_NAME(lcname, UCNAME)  magma_minproductf_##lcname##_gpu_
#elif defined(NOCHANGE)
#define MAGMA_minproduct_GPU_FORTRAN_NAME(lcname, UCNAME)  magma_minproductf_##lcname##_gpu
#elif defined(UPCASE)
#define MAGMA_minproduct_GPU_FORTRAN_NAME(lcname, UCNAME)  MAGMA_minproductF_##UCNAME##_GPU
#endif
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct function definitions / Data on CPU
*/
#define MAGMA_minproductF_ZGEBRD  MAGMA_minproduct_FORTRAN_NAME(zgebrd,  ZGEBRD ) 
#define MAGMA_minproductF_ZGEHRD2 MAGMA_minproduct_FORTRAN_NAME(zgehrd2, ZGEHRD2)
#define MAGMA_minproductF_ZGEHRD  MAGMA_minproduct_FORTRAN_NAME(zgehrd,  ZGEHRD )
#define MAGMA_minproductF_ZGELQF  MAGMA_minproduct_FORTRAN_NAME(zgelqf,  ZGELQF )
#define MAGMA_minproductF_ZGEQLF  MAGMA_minproduct_FORTRAN_NAME(zgeqlf,  ZGEQLF )
#define MAGMA_minproductF_ZGEQRF  MAGMA_minproduct_FORTRAN_NAME(zgeqrf,  ZGEQRF )
#define MAGMA_minproductF_ZGESV   MAGMA_minproduct_FORTRAN_NAME(zgesv,   ZGESV  )
#define MAGMA_minproductF_ZGETRF  MAGMA_minproduct_FORTRAN_NAME(zgetrf,  ZGETRF )
#define MAGMA_minproductF_ZLATRD  MAGMA_minproduct_FORTRAN_NAME(zlatrd,  ZLATRD )
#define MAGMA_minproductF_ZLAHR2  MAGMA_minproduct_FORTRAN_NAME(zlahr2,  ZLAHR2 )
#define MAGMA_minproductF_ZLAHRU  MAGMA_minproduct_FORTRAN_NAME(zlahru,  ZLAHRU )
#define MAGMA_minproductF_ZPOSV   MAGMA_minproduct_FORTRAN_NAME(zposv,   ZPOSV  )
#define MAGMA_minproductF_ZPOTRF  MAGMA_minproduct_FORTRAN_NAME(zpotrf,  ZPOTRF )
#define MAGMA_minproductF_ZHETRD  MAGMA_minproduct_FORTRAN_NAME(zhetrd,  ZHETRD )
#define MAGMA_minproductF_ZUNGQR  MAGMA_minproduct_FORTRAN_NAME(zungqr,  ZUNGQR )
#define MAGMA_minproductF_ZUNMQR  MAGMA_minproduct_FORTRAN_NAME(zunmqr,  ZUNMQR )
#define MAGMA_minproductF_ZUNMTR  MAGMA_minproduct_FORTRAN_NAME(zunmtr,  ZUNMTR )
#define MAGMA_minproductF_ZUNGHR  MAGMA_minproduct_FORTRAN_NAME(zunghr,  ZUNGHR )
#define MAGMA_minproductF_ZGEEV   MAGMA_minproduct_FORTRAN_NAME(zgeev,   ZGEEV  )
#define MAGMA_minproductF_ZGESVD  MAGMA_minproduct_FORTRAN_NAME(zgesvd,  ZGESVD )
#define MAGMA_minproductF_ZHEEVD  MAGMA_minproduct_FORTRAN_NAME(zheevd,  ZHEEVD )
#define MAGMA_minproductF_ZHEGVD  MAGMA_minproduct_FORTRAN_NAME(zhegvd,  ZHEGVD )

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA_minproduct function definitions / Data on GPU
*/
#define MAGMA_minproductF_ZGELS_GPU   MAGMA_minproduct_GPU_FORTRAN_NAME(zgels,   ZGELS  )
#define MAGMA_minproductF_ZGEQRF_GPU  MAGMA_minproduct_GPU_FORTRAN_NAME(zgeqrf,  ZGEQRF ) 
#define MAGMA_minproductF_ZGEQRF2_GPU MAGMA_minproduct_GPU_FORTRAN_NAME(zgeqrf2, ZGEQRF2)
#define MAGMA_minproductF_ZGEQRF3_GPU MAGMA_minproduct_GPU_FORTRAN_NAME(zgeqrf3, ZGEQRF3)
#define MAGMA_minproductF_ZGEQRS_GPU  MAGMA_minproduct_GPU_FORTRAN_NAME(zgeqrs,  ZGEQRS ) 
#define MAGMA_minproductF_ZGEQRS3_GPU MAGMA_minproduct_GPU_FORTRAN_NAME(zgeqrs3, ZGEQRS3) 
#define MAGMA_minproductF_ZGESSM_GPU  MAGMA_minproduct_GPU_FORTRAN_NAME(zgessm,  ZGESSM ) 
#define MAGMA_minproductF_ZGESV_GPU   MAGMA_minproduct_GPU_FORTRAN_NAME(zgesv,   ZGESV  )  
#define MAGMA_minproductF_ZGETRF_INCPIV_GPU  MAGMA_minproduct_GPU_FORTRAN_NAME(zgetrf_incpiv,  ZGETRF_INCPIV ) 
#define MAGMA_minproductF_ZGETRF_GPU  MAGMA_minproduct_GPU_FORTRAN_NAME(zgetrf,  ZGETRF ) 
#define MAGMA_minproductF_ZGETRS_GPU  MAGMA_minproduct_GPU_FORTRAN_NAME(zgetrs,  ZGETRS ) 
#define MAGMA_minproductF_ZLABRD_GPU  MAGMA_minproduct_GPU_FORTRAN_NAME(zlabrd,  ZLABRD ) 
#define MAGMA_minproductF_ZLARFB_GPU  MAGMA_minproduct_GPU_FORTRAN_NAME(zlarfb,  ZLARFB ) 
#define MAGMA_minproductF_ZPOSV_GPU   MAGMA_minproduct_GPU_FORTRAN_NAME(zposv,   ZPOSV  )  
#define MAGMA_minproductF_ZPOTRF_GPU  MAGMA_minproduct_GPU_FORTRAN_NAME(zpotrf,  ZPOTRF ) 
#define MAGMA_minproductF_ZPOTRS_GPU  MAGMA_minproduct_GPU_FORTRAN_NAME(zpotrs,  ZPOTRS ) 
#define MAGMA_minproductF_ZSSSSM_GPU  MAGMA_minproduct_GPU_FORTRAN_NAME(zssssm,  ZSSSSM ) 
#define MAGMA_minproductF_ZTSTRF_GPU  MAGMA_minproduct_GPU_FORTRAN_NAME(ztstrf,  ZTSTRF ) 
#define MAGMA_minproductF_ZUNGQR_GPU  MAGMA_minproduct_GPU_FORTRAN_NAME(zungqr,  ZUNGQR ) 
#define MAGMA_minproductF_ZUNMQR_GPU  MAGMA_minproduct_GPU_FORTRAN_NAME(zunmqr,  ZUNMQR ) 

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (simple interface)
 **/

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct function definitions / Data on CPU
*/
void MAGMA_minproductF_ZGEBRD( magma_minproduct_int_t *m, magma_minproduct_int_t *n, cuDoubleComplex *A, 
                    magma_minproduct_int_t *lda, double *d, double *e,
                    cuDoubleComplex *tauq, cuDoubleComplex *taup, 
                    cuDoubleComplex *work, magma_minproduct_int_t *lwork, magma_minproduct_int_t *info)
{
    magma_minproduct_zgebrd( *m, *n, A, 
                  *lda, d, e,
                  tauq, taup, 
                  work, *lwork, info);
}
    
void MAGMA_minproductF_ZGEHRD2(magma_minproduct_int_t *n, magma_minproduct_int_t *ilo, magma_minproduct_int_t *ihi,
                    cuDoubleComplex *A, magma_minproduct_int_t *lda, cuDoubleComplex *tau, 
                    cuDoubleComplex *work, magma_minproduct_int_t *lwork, magma_minproduct_int_t *info)
{
    magma_minproduct_zgehrd2(*n, *ilo, *ihi,
                  A, *lda, tau, 
                  work, lwork, info);
}
    
void MAGMA_minproductF_ZGEHRD( magma_minproduct_int_t *n, magma_minproduct_int_t *ilo, magma_minproduct_int_t *ihi,
            cuDoubleComplex *A, magma_minproduct_int_t *lda, cuDoubleComplex *tau,
            cuDoubleComplex *work, magma_minproduct_int_t *lwork,
            cuDoubleComplex *d_T, magma_minproduct_int_t *info)
{
  magma_minproduct_zgehrd( *n, *ilo, *ihi,
        A, *lda, tau,
        work, *lwork,
        d_T, info);
}

void MAGMA_minproductF_ZGELQF( magma_minproduct_int_t *m, magma_minproduct_int_t *n, 
                    cuDoubleComplex *A,    magma_minproduct_int_t *lda,   cuDoubleComplex *tau, 
                    cuDoubleComplex *work, magma_minproduct_int_t *lwork, magma_minproduct_int_t *info)
{
    magma_minproduct_zgelqf( *m, *n, 
                  A,    *lda,   tau, 
                  work, *lwork, info);
}

void MAGMA_minproductF_ZGEQLF( magma_minproduct_int_t *m, magma_minproduct_int_t *n, 
                    cuDoubleComplex *A,    magma_minproduct_int_t *lda,   cuDoubleComplex *tau, 
                    cuDoubleComplex *work, magma_minproduct_int_t *lwork, magma_minproduct_int_t *info)
{
    magma_minproduct_zgeqlf( *m, *n, 
                  A,    *lda,   tau, 
                  work, *lwork, info);
}

void MAGMA_minproductF_ZGEQRF( magma_minproduct_int_t *m, magma_minproduct_int_t *n, cuDoubleComplex *A, 
                    magma_minproduct_int_t *lda, cuDoubleComplex *tau, cuDoubleComplex *work, 
                    magma_minproduct_int_t *lwork, magma_minproduct_int_t *info)
{
    magma_minproduct_zgeqrf( *m, *n, A, 
                  *lda, tau, work, 
                  *lwork, info);
}

void MAGMA_minproductF_ZGESV ( magma_minproduct_int_t *n, magma_minproduct_int_t *nrhs,
                    cuDoubleComplex *A, magma_minproduct_int_t *lda, magma_minproduct_int_t *ipiv,
                    cuDoubleComplex *B, magma_minproduct_int_t *ldb, magma_minproduct_int_t *info)
{
    magma_minproduct_zgesv(  *n, *nrhs,
                  A, *lda, ipiv,
                  B, *ldb,
                  info);
}
    
void MAGMA_minproductF_ZGETRF( magma_minproduct_int_t *m, magma_minproduct_int_t *n, cuDoubleComplex *A, 
                    magma_minproduct_int_t *lda, magma_minproduct_int_t *ipiv, 
                    magma_minproduct_int_t *info)
{
    magma_minproduct_zgetrf( *m, *n, A, 
                  *lda, ipiv, 
                  info);
}

// void MAGMA_minproductF_ZLATRD( char *uplo, magma_minproduct_int_t *n, magma_minproduct_int_t *nb, cuDoubleComplex *a, 
//                     magma_minproduct_int_t *lda, double *e, cuDoubleComplex *tau, 
//                     cuDoubleComplex *w, magma_minproduct_int_t *ldw,
//                     cuDoubleComplex *da, magma_minproduct_int_t *ldda, 
//                     cuDoubleComplex *dw, magma_minproduct_int_t *lddw)
// {
//     magma_minproduct_zlatrd( uplo[0], *n, *nb, a, 
//                   *lda, e, tau, 
//                   w, *ldw,
//                   da, *ldda, 
//                   dw, *lddw);
// }

  /* This has nothing to do here, it should be a GPU function */
// void MAGMA_minproductF_ZLAHR2( magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *nb, 
//                     cuDoubleComplex *da, cuDoubleComplex *dv, cuDoubleComplex *a, 
//                     magma_minproduct_int_t *lda, cuDoubleComplex *tau, cuDoubleComplex *t, 
//                     magma_minproduct_int_t *ldt, cuDoubleComplex *y, magma_minproduct_int_t *ldy)
// {
//     magma_minproduct_zlahr2( *m, *n, *nb, 
//                   da, dv, a, 
//                   *lda, tau, t, 
//                   *ldt, y, *ldy);
// }

// void MAGMA_minproductF_ZLAHRU( magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *nb, 
//                     cuDoubleComplex *a, magma_minproduct_int_t *lda, 
//                     cuDoubleComplex *da, cuDoubleComplex *y, 
//                     cuDoubleComplex *v, cuDoubleComplex *t, 
//                     cuDoubleComplex *dwork)
// {
//     magma_minproduct_zlahru( *m, *n, *nb, 
//                   a, *lda, 
//                   da, y, 
//                   v, t, 
//                   dwork);
// }

void MAGMA_minproductF_ZPOSV(  char *uplo, magma_minproduct_int_t *n, magma_minproduct_int_t *nrhs,
                    cuDoubleComplex *A, magma_minproduct_int_t *lda,
                    cuDoubleComplex *B, magma_minproduct_int_t *ldb, magma_minproduct_int_t *info)
{
    magma_minproduct_zposv(  uplo[0], *n, *nrhs,
                  A, *lda,
                  B, *ldb, info);
}

void MAGMA_minproductF_ZPOTRF( char *uplo, magma_minproduct_int_t *n, cuDoubleComplex *A, 
                    magma_minproduct_int_t *lda, magma_minproduct_int_t *info)
{
    magma_minproduct_zpotrf( uplo[0], *n, A, 
                  *lda, info);
}

void MAGMA_minproductF_ZHETRD( char *uplo, magma_minproduct_int_t *n, cuDoubleComplex *A, 
                    magma_minproduct_int_t *lda, double *d, double *e, 
                    cuDoubleComplex *tau, cuDoubleComplex *work, magma_minproduct_int_t *lwork, 
                    magma_minproduct_int_t *info)
{
    magma_minproduct_zhetrd( uplo[0], *n, A, 
                  *lda, d, e, 
                  tau, work, *lwork, 
                  info);
}

// void MAGMA_minproductF_ZUNGQR( magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *k,
//                     cuDoubleComplex *a, magma_minproduct_int_t *lda,
//                     cuDoubleComplex *tau, cuDoubleComplex *dwork,
//                     magma_minproduct_int_t *nb, magma_minproduct_int_t *info )
// {
//     magma_minproduct_zungqr( *m, *n, *k,
//                   a, *lda,
//                   tau, dwork,
//                   *nb, info );
// }

void MAGMA_minproductF_ZUNMQR( char *side, char *trans, 
                    magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *k, 
                    cuDoubleComplex *a, magma_minproduct_int_t *lda, cuDoubleComplex *tau, 
                    cuDoubleComplex *c, magma_minproduct_int_t *ldc, 
                    cuDoubleComplex *work, magma_minproduct_int_t *lwork, magma_minproduct_int_t *info)
{
    magma_minproduct_zunmqr( side[0], trans[0], 
                  *m, *n, *k, 
                  a, *lda, tau, 
                  c, *ldc, 
                  work, *lwork, info);
}

void MAGMA_minproductF_ZUNMTR( char *side, char *uplo, char *trans,
                    magma_minproduct_int_t *m, magma_minproduct_int_t *n,
                    cuDoubleComplex *a,    magma_minproduct_int_t *lda,
                    cuDoubleComplex *tau,
                    cuDoubleComplex *c,    magma_minproduct_int_t *ldc,
                    cuDoubleComplex *work, magma_minproduct_int_t *lwork,
                    magma_minproduct_int_t *info)
{
    magma_minproduct_zunmtr( side[0], uplo[0], trans[0],
                  *m, *n,
                  a,    *lda,
                  tau,
                  c,    *ldc,
                  work, *lwork,
                  info);
}

// void MAGMA_minproductF_ZUNGHR( magma_minproduct_int_t *n, magma_minproduct_int_t *ilo, magma_minproduct_int_t *ihi,
//                     cuDoubleComplex *a, magma_minproduct_int_t *lda,
//                     cuDoubleComplex *tau,
//                     cuDoubleComplex *dT, magma_minproduct_int_t *nb,
//                     magma_minproduct_int_t *info)
// {
//     magma_minproduct_zunghr( *n, *ilo, *ihi,
//                   a, *lda,
//                   tau,
//                   dT, *nb,
//                   info);
// }

#if defined(PRECISION_z) || defined(PRECISION_c)
void MAGMA_minproductF_ZGEEV( char *jobvl, char *jobvr, magma_minproduct_int_t *n,
                   cuDoubleComplex *a, magma_minproduct_int_t *lda,
                   cuDoubleComplex *w,
                   cuDoubleComplex *vl, magma_minproduct_int_t *ldvl,
                   cuDoubleComplex *vr, magma_minproduct_int_t *ldvr,
                   cuDoubleComplex *work, magma_minproduct_int_t *lwork,
                   double *rwork, magma_minproduct_int_t *info)
{
    magma_minproduct_zgeev( jobvl[0], jobvr[0], *n,
                 a, *lda,
                 w,
                 vl, *ldvl,
                 vr, *ldvr,
                 work, *lwork,
                 rwork, info);
}

void MAGMA_minproductF_ZGESVD( char *jobu, char *jobvt, magma_minproduct_int_t *m, magma_minproduct_int_t *n,
                    cuDoubleComplex *a,    magma_minproduct_int_t *lda, double *s, 
                    cuDoubleComplex *u,    magma_minproduct_int_t *ldu, 
                    cuDoubleComplex *vt,   magma_minproduct_int_t *ldvt,
                    cuDoubleComplex *work, magma_minproduct_int_t *lwork,
                    double *rwork, magma_minproduct_int_t *info )
{
    magma_minproduct_zgesvd( jobu[0], jobvt[0], *m, *n,
                  a,    *lda, s, 
                  u,    *ldu, 
                  vt,   *ldvt,
                  work, *lwork,
                  rwork, info );
}
    
void MAGMA_minproductF_ZHEEVD( char *jobz, char *uplo, magma_minproduct_int_t *n,
                    cuDoubleComplex *a,     magma_minproduct_int_t *lda, double *w,
                    cuDoubleComplex *work,  magma_minproduct_int_t *lwork,
                    double          *rwork, magma_minproduct_int_t *lrwork,
                    magma_minproduct_int_t *iwork, magma_minproduct_int_t *liwork, magma_minproduct_int_t *info)
{
    magma_minproduct_zheevd( jobz[0], uplo[0], *n,
                  a, *lda, w,
                  work, *lwork,
                  rwork, *lrwork,
                  iwork, *liwork, info);
}

void MAGMA_minproductF_ZHEGVD(magma_minproduct_int_t *itype, char *jobz, char *uplo, magma_minproduct_int_t *n,
           cuDoubleComplex *a, magma_minproduct_int_t *lda, 
           cuDoubleComplex *b, magma_minproduct_int_t *ldb,
           double *w, cuDoubleComplex *work, magma_minproduct_int_t *lwork,
           double *rwork, magma_minproduct_int_t *lrwork,
           magma_minproduct_int_t *iwork, magma_minproduct_int_t *liwork, magma_minproduct_int_t *info)
{
  magma_minproduct_zhegvd( *itype, jobz[0], uplo[0], *n,
        a, *lda, b, *ldb,
        w, work, *lwork,
        rwork, *lrwork,
        iwork, *liwork, info);
}
    
#else
void MAGMA_minproductF_ZGEEV( char *jobvl, char *jobvr, magma_minproduct_int_t *n,
                   cuDoubleComplex *a,    magma_minproduct_int_t *lda,
                   cuDoubleComplex *wr, cuDoubleComplex *wi,
                   cuDoubleComplex *vl,   magma_minproduct_int_t *ldvl,
                   cuDoubleComplex *vr,   magma_minproduct_int_t *ldvr,
                   cuDoubleComplex *work, magma_minproduct_int_t *lwork,
                   magma_minproduct_int_t *info)
{
    magma_minproduct_zgeev( jobvl[0], jobvr[0], *n,
                 a,    *lda,
                 wr, wi,
                 vl,   *ldvl,
                 vr,   *ldvr,
                 work, *lwork,
                 info);
}

void MAGMA_minproductF_ZGESVD( char *jobu, char *jobvt, magma_minproduct_int_t *m, magma_minproduct_int_t *n,
                    cuDoubleComplex *a,    magma_minproduct_int_t *lda, double *s,
                    cuDoubleComplex *u,    magma_minproduct_int_t *ldu, 
                    cuDoubleComplex *vt,   magma_minproduct_int_t *ldvt,
                    cuDoubleComplex *work, magma_minproduct_int_t *lwork,
                    magma_minproduct_int_t *info )
{
    magma_minproduct_zgesvd( jobu[0], jobvt[0], *m, *n,
                  a,    *lda, s, 
                  u,    *ldu, 
                  vt,   *ldvt,
                  work, *lwork,
                  info );
}

void MAGMA_minproductF_ZHEEVD( char *jobz, char *uplo, magma_minproduct_int_t *n,
                    cuDoubleComplex *a, magma_minproduct_int_t *lda, double *w,
                    cuDoubleComplex *work, magma_minproduct_int_t *lwork,
                    magma_minproduct_int_t *iwork, magma_minproduct_int_t *liwork, magma_minproduct_int_t *info)
{
    magma_minproduct_zheevd( jobz[0], uplo[0], *n,
                  a, *lda, w,
                  work, *lwork,
                  iwork, *liwork, info);
}

void MAGMA_minproductF_ZHEGVD(magma_minproduct_int_t *itype, char *jobz, char *uplo, magma_minproduct_int_t *n,
                   cuDoubleComplex *a, magma_minproduct_int_t *lda,
                   cuDoubleComplex *b, magma_minproduct_int_t *ldb,
                   double *w, cuDoubleComplex *work, magma_minproduct_int_t *lwork,
                   magma_minproduct_int_t *iwork, magma_minproduct_int_t *liwork, magma_minproduct_int_t *info)
{
  magma_minproduct_zhegvd( *itype, jobz[0], uplo[0], *n,
                a, *lda, b, *ldb,
        w, work, *lwork,
                iwork, *liwork, info);
}


#endif

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA_minproduct function definitions / Data on GPU
*/
void MAGMA_minproductF_ZGELS_GPU(  char *trans, magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *nrhs,
                        devptr_t *dA,    magma_minproduct_int_t *ldda, 
                        devptr_t *dB,    magma_minproduct_int_t *lddb, 
                        cuDoubleComplex *hwork, magma_minproduct_int_t *lwork, 
                        magma_minproduct_int_t *info)
{
    magma_minproduct_zgels_gpu(  trans[0], *m, *n, *nrhs, 
                      DEVPTR(dA),    *ldda,  
                      DEVPTR(dB),    *lddb,  
                      hwork, *lwork,  info);
}

void MAGMA_minproductF_ZGEQRF_GPU( magma_minproduct_int_t *m, magma_minproduct_int_t *n, 
                        devptr_t *dA,  magma_minproduct_int_t *ldda, 
                        cuDoubleComplex *tau, devptr_t *dT, 
                        magma_minproduct_int_t *info)
{
    magma_minproduct_zgeqrf_gpu( *m, *n,  
                      DEVPTR(dA),  *ldda,  
                      tau, 
                      DEVPTR(dT),  info);
}

void MAGMA_minproductF_ZGEQRF2_GPU(magma_minproduct_int_t *m, magma_minproduct_int_t *n, 
                        devptr_t *dA,  magma_minproduct_int_t *ldda, 
                        cuDoubleComplex *tau, magma_minproduct_int_t *info)
{
    magma_minproduct_zgeqrf2_gpu(*m, *n,  
                      DEVPTR(dA),  *ldda,  
                      tau, info); 
}

void MAGMA_minproductF_ZGEQRF3_GPU(magma_minproduct_int_t *m, magma_minproduct_int_t *n, 
                        devptr_t *dA,  magma_minproduct_int_t *ldda, 
                        cuDoubleComplex *tau, devptr_t *dT,
                        magma_minproduct_int_t *info)
{
    magma_minproduct_zgeqrf3_gpu(*m, *n,  
                      DEVPTR(dA),  *ldda,  
                      tau, DEVPTR(dT), info); 
}

void MAGMA_minproductF_ZGEQRS_GPU( magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *nrhs, 
                        devptr_t *dA,     magma_minproduct_int_t *ldda, 
                        cuDoubleComplex *tau,   devptr_t *dT,
                        devptr_t *dB,    magma_minproduct_int_t *lddb,
                        cuDoubleComplex *hwork, magma_minproduct_int_t *lhwork, 
                        magma_minproduct_int_t *info)
{
    magma_minproduct_zgeqrs_gpu( *m, *n, *nrhs,  
                      DEVPTR(dA),     *ldda,  
                      tau,
                      DEVPTR(dT), 
                      DEVPTR(dB),    *lddb, 
                      hwork, *lhwork,  info);
}

void MAGMA_minproductF_ZGEQRS3_GPU(magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *nrhs, 
                        devptr_t *dA,     magma_minproduct_int_t *ldda, 
                        cuDoubleComplex *tau,   devptr_t *dT,
                        devptr_t *dB,    magma_minproduct_int_t *lddb,
                        cuDoubleComplex *hwork, magma_minproduct_int_t *lhwork, 
                        magma_minproduct_int_t *info)
{
    magma_minproduct_zgeqrs3_gpu(*m, *n, *nrhs,  
                      DEVPTR(dA),     *ldda,  
                      tau,
                      DEVPTR(dT), 
                      DEVPTR(dB),    *lddb, 
                      hwork, *lhwork,  info);
}

void MAGMA_minproductF_ZGESSM_GPU( char *storev, magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *k, magma_minproduct_int_t *ib, 
                        magma_minproduct_int_t *ipiv, 
                        devptr_t *dL1, magma_minproduct_int_t *lddl1, 
                        devptr_t *dL,  magma_minproduct_int_t *lddl, 
                        devptr_t *dA,  magma_minproduct_int_t *ldda, 
                        magma_minproduct_int_t *info)
{
    magma_minproduct_zgessm_gpu( storev[0], *m, *n, *k, *ib, ipiv,  
                      DEVPTR(dL1), *lddl1,  
                      DEVPTR(dL),  *lddl,  
                      DEVPTR(dA),  *ldda,  info);
}

void MAGMA_minproductF_ZGESV_GPU(  magma_minproduct_int_t *n, magma_minproduct_int_t *nrhs, 
                        devptr_t *dA, magma_minproduct_int_t *ldda, magma_minproduct_int_t *ipiv, 
                        devptr_t *dB, magma_minproduct_int_t *lddb, magma_minproduct_int_t *info)
{
    magma_minproduct_zgesv_gpu(  *n, *nrhs,  
                      DEVPTR(dA), *ldda, ipiv,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMA_minproductF_ZGETRF_GPU( magma_minproduct_int_t *m, magma_minproduct_int_t *n, 
                        devptr_t *dA, magma_minproduct_int_t *ldda, 
                        magma_minproduct_int_t *ipiv, magma_minproduct_int_t *info)
{
    magma_minproduct_zgetrf_gpu( *m, *n,  
                      DEVPTR(dA), *ldda, ipiv, info);
}

void MAGMA_minproductF_ZGETRS_GPU( char *trans, magma_minproduct_int_t *n, magma_minproduct_int_t *nrhs, 
                        devptr_t *dA, magma_minproduct_int_t *ldda, magma_minproduct_int_t *ipiv, 
                        devptr_t *dB, magma_minproduct_int_t *lddb, magma_minproduct_int_t *info)
{
    magma_minproduct_zgetrs_gpu( trans[0], *n, *nrhs,  
                      DEVPTR(dA), *ldda, ipiv,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMA_minproductF_ZLABRD_GPU( magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *nb, 
                        cuDoubleComplex *a, magma_minproduct_int_t *lda, devptr_t *da, magma_minproduct_int_t *ldda,
                        double *d, double *e, cuDoubleComplex *tauq, cuDoubleComplex *taup,  
                        cuDoubleComplex *x, magma_minproduct_int_t *ldx, devptr_t *dx, magma_minproduct_int_t *lddx, 
                        cuDoubleComplex *y, magma_minproduct_int_t *ldy, devptr_t *dy, magma_minproduct_int_t *lddy)
{
    magma_minproduct_zlabrd_gpu( *m, *n, *nb,  
                      a, *lda, DEVPTR(da), *ldda, 
                      d, e, tauq, taup,   
                      x, *ldx, DEVPTR(dx), *lddx,  
                      y, *ldy, DEVPTR(dy), *lddy);
}

void MAGMA_minproductF_ZLARFB_GPU( char *side, char *trans, char *direct, char *storev, 
                        magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *k,
                        devptr_t *dv, magma_minproduct_int_t *ldv, devptr_t *dt,    magma_minproduct_int_t *ldt, 
                        devptr_t *dc, magma_minproduct_int_t *ldc, devptr_t *dowrk, magma_minproduct_int_t *ldwork )
{
    magma_minproduct_zlarfb_gpu( side[0], trans[0], direct[0], storev[0],  *m, *n, *k, 
                      DEVPTR(dv), *ldv, DEVPTR(dt),    *ldt,  
                      DEVPTR(dc), *ldc, DEVPTR(dowrk), *ldwork);
}

void MAGMA_minproductF_ZPOSV_GPU(  char *uplo, magma_minproduct_int_t *n, magma_minproduct_int_t *nrhs, 
                        devptr_t *dA, magma_minproduct_int_t *ldda, 
                        devptr_t *dB, magma_minproduct_int_t *lddb, magma_minproduct_int_t *info)
{
    magma_minproduct_zposv_gpu(  uplo[0], *n, *nrhs,  
                      DEVPTR(dA), *ldda,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMA_minproductF_ZPOTRF_GPU( char *uplo,  magma_minproduct_int_t *n, 
                        devptr_t *dA, magma_minproduct_int_t *ldda, magma_minproduct_int_t *info)
{
    magma_minproduct_zpotrf_gpu( uplo[0],  *n,  
                      DEVPTR(dA), *ldda, info); }

void MAGMA_minproductF_ZPOTRS_GPU( char *uplo,  magma_minproduct_int_t *n, magma_minproduct_int_t *nrhs, 
                        devptr_t *dA, magma_minproduct_int_t *ldda, 
                        devptr_t *dB, magma_minproduct_int_t *lddb, magma_minproduct_int_t *info)
{
    magma_minproduct_zpotrs_gpu( uplo[0],  *n, *nrhs,  
                      DEVPTR(dA), *ldda,  
                      DEVPTR(dB), *lddb, info);
}

void MAGMA_minproductF_ZSSSSM_GPU( char *storev, magma_minproduct_int_t *m1, magma_minproduct_int_t *n1, 
                        magma_minproduct_int_t *m2, magma_minproduct_int_t *n2, magma_minproduct_int_t *k, magma_minproduct_int_t *ib, 
                        devptr_t *dA1, magma_minproduct_int_t *ldda1, 
                        devptr_t *dA2, magma_minproduct_int_t *ldda2, 
                        devptr_t *dL1, magma_minproduct_int_t *lddl1, 
                        devptr_t *dL2, magma_minproduct_int_t *lddl2,
                        magma_minproduct_int_t *IPIV, magma_minproduct_int_t *info)
{
    magma_minproduct_zssssm_gpu( storev[0], *m1, *n1,  *m2, *n2, *k, *ib,  
                      DEVPTR(dA1), *ldda1,  
                      DEVPTR(dA2), *ldda2,  
                      DEVPTR(dL1), *lddl1,  
                      DEVPTR(dL2), *lddl2,
                      IPIV, info);
}

void MAGMA_minproductF_ZUNGQR_GPU( magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *k, 
                        devptr_t *da, magma_minproduct_int_t *ldda, 
                        cuDoubleComplex *tau, devptr_t *dwork, 
                        magma_minproduct_int_t *nb, magma_minproduct_int_t *info )
{
    magma_minproduct_zungqr_gpu( *m, *n, *k,  
                      DEVPTR(da), *ldda, tau, 
                      DEVPTR(dwork), *nb, info );
}

void MAGMA_minproductF_ZUNMQR_GPU( char *side, char *trans, 
                        magma_minproduct_int_t *m, magma_minproduct_int_t *n, magma_minproduct_int_t *k,
                        devptr_t *a,    magma_minproduct_int_t *lda, cuDoubleComplex *tau, 
                        devptr_t *c,    magma_minproduct_int_t *ldc,
                        devptr_t *work, magma_minproduct_int_t *lwork, 
                        devptr_t *td,   magma_minproduct_int_t *nb, magma_minproduct_int_t *info)
{
    magma_minproduct_zunmqr_gpu( side[0], trans[0], *m, *n, *k, 
                      DEVPTR(a),    *lda, tau,  
                      DEVPTR(c),    *ldc, 
                      DEVPTR(work), *lwork,  
                      DEVPTR(td),   *nb, info);
}

#ifdef __cplusplus
}
#endif
