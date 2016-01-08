/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
*/

#ifndef MAGMA_tally4_ZLAPACK_H
#define MAGMA_tally4_ZLAPACK_H

#include "magma_tally4_types.h"
#include "magma_tally4_mangling.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- BLAS and LAPACK functions (alphabetical order)
*/
#define blasf77_izamax     FORTRAN_NAME( izamax, IZAMAX )
#define blasf77_zaxpy      FORTRAN_NAME( zaxpy,  ZAXPY  )
#define blasf77_zcopy      FORTRAN_NAME( zcopy,  ZCOPY  )
#define blasf77_zgemm      FORTRAN_NAME( zgemm,  ZGEMM  )
#define blasf77_zgemv      FORTRAN_NAME( zgemv,  ZGEMV  )
#define blasf77_zgerc      FORTRAN_NAME( zgerc,  ZGERC  )
#define blasf77_zgeru      FORTRAN_NAME( zgeru,  ZGERU  )
#define blasf77_zhemm      FORTRAN_NAME( zhemm,  ZHEMM  )
#define blasf77_zhemv      FORTRAN_NAME( zhemv,  ZHEMV  )
#define blasf77_zher       FORTRAN_NAME( zher,   ZHER   )
#define blasf77_zher2      FORTRAN_NAME( zher2,  ZHER2  )
#define blasf77_zher2k     FORTRAN_NAME( zher2k, ZHER2K )
#define blasf77_zherk      FORTRAN_NAME( zherk,  ZHERK  )
#define blasf77_zscal      FORTRAN_NAME( zscal,  ZSCAL  )
#define blasf77_zdscal     FORTRAN_NAME( zdscal, ZDSCAL )
#define blasf77_zswap      FORTRAN_NAME( zswap,  ZSWAP  )
#define blasf77_zsymm      FORTRAN_NAME( zsymm,  ZSYMM  )
#define blasf77_zsyr2k     FORTRAN_NAME( zsyr2k, ZSYR2K )
#define blasf77_zsyrk      FORTRAN_NAME( zsyrk,  ZSYRK  )
#define blasf77_zrotg      FORTRAN_NAME( zrotg,  ZROTG  )
#define blasf77_zrot       FORTRAN_NAME( zrot,   ZROT   )
#define blasf77_zdrot      FORTRAN_NAME( zdrot,  ZDROT  )
#define blasf77_ztrmm      FORTRAN_NAME( ztrmm,  ZTRMM  )
#define blasf77_ztrmv      FORTRAN_NAME( ztrmv,  ZTRMV  )
#define blasf77_ztrsm      FORTRAN_NAME( ztrsm,  ZTRSM  )
#define blasf77_ztrsv      FORTRAN_NAME( ztrsv,  ZTRSV  )

#define lapackf77_dlaed2   FORTRAN_NAME( dlaed2, DLAED2 )
#define lapackf77_dlaed4   FORTRAN_NAME( dlaed4, DLAED4 )
#define lapackf77_dlaln2   FORTRAN_NAME( dlaln2, DLALN2 )
#define lapackf77_dlamc3   FORTRAN_NAME( dlamc3, DLAMC3 )
#define lapackf77_dlamrg   FORTRAN_NAME( dlamrg, DLAMRG )
#define lapackf77_dlasrt   FORTRAN_NAME( dlasrt, DLASRT )
#define lapackf77_dstebz   FORTRAN_NAME( dstebz, DSTEBZ )

#define lapackf77_dbdsdc   FORTRAN_NAME( dbdsdc, DBDSDC )
#define lapackf77_zbdsqr   FORTRAN_NAME( zbdsqr, ZBDSQR )
#define lapackf77_zgebak   FORTRAN_NAME( zgebak, ZGEBAK )
#define lapackf77_zgebal   FORTRAN_NAME( zgebal, ZGEBAL )
#define lapackf77_zgebd2   FORTRAN_NAME( zgebd2, ZGEBD2 )
#define lapackf77_zgebrd   FORTRAN_NAME( zgebrd, ZGEBRD )
#define lapackf77_zgbbrd   FORTRAN_NAME( zgbbrd, ZGBBRD )
#define lapackf77_zgeev    FORTRAN_NAME( zgeev,  ZGEEV  )
#define lapackf77_zgehd2   FORTRAN_NAME( zgehd2, ZGEHD2 )
#define lapackf77_zgehrd   FORTRAN_NAME( zgehrd, ZGEHRD )
#define lapackf77_zgelqf   FORTRAN_NAME( zgelqf, ZGELQF )
#define lapackf77_zgels    FORTRAN_NAME( zgels,  ZGELS  )
#define lapackf77_zgeqlf   FORTRAN_NAME( zgeqlf, ZGEQLF )
#define lapackf77_zgeqp3   FORTRAN_NAME( zgeqp3, ZGEQP3 )
#define lapackf77_zgeqrf   FORTRAN_NAME( zgeqrf, ZGEQRF )
#define lapackf77_zgesdd   FORTRAN_NAME( zgesdd, ZGESDD )
#define lapackf77_zgesv    FORTRAN_NAME( zgesv,  ZGESV  )
#define lapackf77_zgesvd   FORTRAN_NAME( zgesvd, ZGESVD )
#define lapackf77_zgetrf   FORTRAN_NAME( zgetrf, ZGETRF )
#define lapackf77_zgetri   FORTRAN_NAME( zgetri, ZGETRI )
#define lapackf77_zgetrs   FORTRAN_NAME( zgetrs, ZGETRS )
#define lapackf77_zhetf2   FORTRAN_NAME( zhetf2, ZHETF2 )
#define lapackf77_zhetrs   FORTRAN_NAME( zhetrs, ZHETRS )
#define lapackf77_zhbtrd   FORTRAN_NAME( zhbtrd, ZHBTRD )
#define lapackf77_zheev    FORTRAN_NAME( zheev,  ZHEEV  )
#define lapackf77_zheevd   FORTRAN_NAME( zheevd, ZHEEVD )
#define lapackf77_zheevr   FORTRAN_NAME( zheevr, ZHEEVR )
#define lapackf77_zheevx   FORTRAN_NAME( zheevx, ZHEEVX )
#define lapackf77_zhegs2   FORTRAN_NAME( zhegs2, ZHEGS2 )
#define lapackf77_zhegst   FORTRAN_NAME( zhegst, ZHEGST )
#define lapackf77_zhegvd   FORTRAN_NAME( zhegvd, ZHEGVD )
#define lapackf77_zhesv    FORTRAN_NAME( zhesv,  ZHESV  )
#define lapackf77_zhetd2   FORTRAN_NAME( zhetd2, ZHETD2 )
#define lapackf77_zhetrd   FORTRAN_NAME( zhetrd, ZHETRD )
#define lapackf77_zhetrf   FORTRAN_NAME( zhetrf, ZHETRF )
#define lapackf77_zhseqr   FORTRAN_NAME( zhseqr, ZHSEQR )
#define lapackf77_zlabrd   FORTRAN_NAME( zlabrd, ZLABRD )
#define lapackf77_zlacgv   FORTRAN_NAME( zlacgv, ZLACGV )
#define lapackf77_zlacp2   FORTRAN_NAME( zlacp2, ZLACP2 )
#define lapackf77_zlacpy   FORTRAN_NAME( zlacpy, ZLACPY )
#define lapackf77_zlacrm   FORTRAN_NAME( zlacrm, ZLACRM )
#define lapackf77_zladiv   FORTRAN_NAME( zladiv, ZLADIV )
#define lapackf77_zlahef   FORTRAN_NAME( zlahef, ZLAHEF )
#define lapackf77_zlange   FORTRAN_NAME( zlange, ZLANGE )
#define lapackf77_zlanhe   FORTRAN_NAME( zlanhe, ZLANHE )
#define lapackf77_zlanht   FORTRAN_NAME( zlanht, ZLANHT )
#define lapackf77_zlansy   FORTRAN_NAME( zlansy, ZLANSY )
#define lapackf77_zlantr   FORTRAN_NAME( zlantr, ZLANTR )
#define lapackf77_dlapy3   FORTRAN_NAME( dlapy3, DLAPY3 )
#define lapackf77_zlaqp2   FORTRAN_NAME( zlaqp2, ZLAQP2 )
#define lapackf77_zlarcm   FORTRAN_NAME( zlarcm, ZLARCM )
#define lapackf77_zlarf    FORTRAN_NAME( zlarf,  ZLARF  )
#define lapackf77_zlarfb   FORTRAN_NAME( zlarfb, ZLARFB )
#define lapackf77_zlarfg   FORTRAN_NAME( zlarfg, ZLARFG )
#define lapackf77_zlarft   FORTRAN_NAME( zlarft, ZLARFT )
#define lapackf77_zlarnv   FORTRAN_NAME( zlarnv, ZLARNV )
#define lapackf77_zlartg   FORTRAN_NAME( zlartg, ZLARTG )
#define lapackf77_zlascl   FORTRAN_NAME( zlascl, ZLASCL )
#define lapackf77_zlaset   FORTRAN_NAME( zlaset, ZLASET )
#define lapackf77_zlaswp   FORTRAN_NAME( zlaswp, ZLASWP )
#define lapackf77_zlatrd   FORTRAN_NAME( zlatrd, ZLATRD )
#define lapackf77_zlatrs   FORTRAN_NAME( zlatrs, ZLATRS )
#define lapackf77_zlauum   FORTRAN_NAME( zlauum, ZLAUUM )
#define lapackf77_zlavhe   FORTRAN_NAME( zlavhe, ZLAVHE )
#define lapackf77_zposv    FORTRAN_NAME( zposv,  ZPOSV  )
#define lapackf77_zpotrf   FORTRAN_NAME( zpotrf, ZPOTRF )
#define lapackf77_zpotri   FORTRAN_NAME( zpotri, ZPOTRI )
#define lapackf77_zpotrs   FORTRAN_NAME( zpotrs, ZPOTRS )
#define lapackf77_zstedc   FORTRAN_NAME( zstedc, ZSTEDC )
#define lapackf77_zstein   FORTRAN_NAME( zstein, ZSTEIN )
#define lapackf77_zstemr   FORTRAN_NAME( zstemr, ZSTEMR )
#define lapackf77_zsteqr   FORTRAN_NAME( zsteqr, ZSTEQR )
#define lapackf77_zsymv    FORTRAN_NAME( zsymv,  ZSYMV  )
#define lapackf77_ztrevc   FORTRAN_NAME( ztrevc, ZTREVC )
#define lapackf77_ztrevc3  FORTRAN_NAME( ztrevc3, ZTREVC3 )
#define lapackf77_ztrtri   FORTRAN_NAME( ztrtri, ZTRTRI )
#define lapackf77_zung2r   FORTRAN_NAME( zung2r, ZUNG2R )
#define lapackf77_zungbr   FORTRAN_NAME( zungbr, ZUNGBR )
#define lapackf77_zunghr   FORTRAN_NAME( zunghr, ZUNGHR )
#define lapackf77_zunglq   FORTRAN_NAME( zunglq, ZUNGLQ )
#define lapackf77_zungql   FORTRAN_NAME( zungql, ZUNGQL )
#define lapackf77_zungqr   FORTRAN_NAME( zungqr, ZUNGQR )
#define lapackf77_zungtr   FORTRAN_NAME( zungtr, ZUNGTR )
#define lapackf77_zunm2r   FORTRAN_NAME( zunm2r, ZUNM2R )
#define lapackf77_zunmbr   FORTRAN_NAME( zunmbr, ZUNMBR )
#define lapackf77_zunmlq   FORTRAN_NAME( zunmlq, ZUNMLQ )
#define lapackf77_zunmql   FORTRAN_NAME( zunmql, ZUNMQL )
#define lapackf77_zunmqr   FORTRAN_NAME( zunmqr, ZUNMQR )
#define lapackf77_zunmtr   FORTRAN_NAME( zunmtr, ZUNMTR )

/* testing functions (alphabetical order) */
#define lapackf77_zbdt01   FORTRAN_NAME( zbdt01, ZBDT01 )
#define lapackf77_zget22   FORTRAN_NAME( zget22, ZGET22 )
#define lapackf77_zhet21   FORTRAN_NAME( zhet21, ZHET21 )
#define lapackf77_zhst01   FORTRAN_NAME( zhst01, ZHST01 )
#define lapackf77_zlarfx   FORTRAN_NAME( zlarfx, ZLARFX )
#define lapackf77_zlarfy   FORTRAN_NAME( zlarfy, ZLARFY )
#define lapackf77_zlatms   FORTRAN_NAME( zlatms, ZLATMS )
#define lapackf77_zqpt01   FORTRAN_NAME( zqpt01, ZQPT01 )
#define lapackf77_zqrt02   FORTRAN_NAME( zqrt02, ZQRT02 )
#define lapackf77_zstt21   FORTRAN_NAME( zstt21, ZSTT21 )
#define lapackf77_zunt01   FORTRAN_NAME( zunt01, ZUNT01 )

/*
 * BLAS functions (alphabetical order)
 */
magma_tally4_int_t blasf77_izamax(
                     const magma_tally4_int_t *n,
                     const magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx );

void blasf77_zaxpy(  const magma_tally4_int_t *n,
                     const magma_tally4DoubleComplex *alpha,
                     const magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx,
                           magma_tally4DoubleComplex *y, const magma_tally4_int_t *incy );

void blasf77_zcopy(  const magma_tally4_int_t *n,
                     const magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx,
                           magma_tally4DoubleComplex *y, const magma_tally4_int_t *incy );

void blasf77_zgemm(  const char *transa, const char *transb,
                     const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                     const magma_tally4DoubleComplex *alpha,
                     const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                     const magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb,
                     const magma_tally4DoubleComplex *beta,
                           magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc );

void blasf77_zgemv(  const char *transa,
                     const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                     const magma_tally4DoubleComplex *alpha,
                     const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                     const magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx,
                     const magma_tally4DoubleComplex *beta,
                           magma_tally4DoubleComplex *y, const magma_tally4_int_t *incy );

void blasf77_zgerc(  const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                     const magma_tally4DoubleComplex *alpha,
                     const magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx,
                     const magma_tally4DoubleComplex *y, const magma_tally4_int_t *incy,
                           magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda );

void blasf77_zgeru(  const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                     const magma_tally4DoubleComplex *alpha,
                     const magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx,
                     const magma_tally4DoubleComplex *y, const magma_tally4_int_t *incy,
                           magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda );

void blasf77_zhemm(  const char *side, const char *uplo,
                     const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                     const magma_tally4DoubleComplex *alpha,
                     const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                     const magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb,
                     const magma_tally4DoubleComplex *beta,
                           magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc );

void blasf77_zhemv(  const char *uplo,
                     const magma_tally4_int_t *n,
                     const magma_tally4DoubleComplex *alpha,
                     const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                     const magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx,
                     const magma_tally4DoubleComplex *beta,
                           magma_tally4DoubleComplex *y, const magma_tally4_int_t *incy );

void blasf77_zher(   const char *uplo,
                     const magma_tally4_int_t *n,
                     const double *alpha,
                     const magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx,
                           magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda );

void blasf77_zher2(  const char *uplo,
                     const magma_tally4_int_t *n,
                     const magma_tally4DoubleComplex *alpha,
                     const magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx,
                     const magma_tally4DoubleComplex *y, const magma_tally4_int_t *incy,
                           magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda );

void blasf77_zher2k( const char *uplo, const char *trans,
                     const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                     const magma_tally4DoubleComplex *alpha,
                     const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                     const magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb,
                     const double *beta,
                           magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc );

void blasf77_zherk(  const char *uplo, const char *trans,
                     const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                     const double *alpha,
                     const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                     const double *beta,
                           magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc );

void blasf77_zscal(  const magma_tally4_int_t *n,
                     const magma_tally4DoubleComplex *alpha,
                           magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx );

void blasf77_zdscal( const magma_tally4_int_t *n,
                     const double *alpha,
                           magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx );

void blasf77_zswap(  const magma_tally4_int_t *n,
                     magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx,
                     magma_tally4DoubleComplex *y, const magma_tally4_int_t *incy );

/* complex-symmetric (non-Hermitian) routines */
void blasf77_zsymm(  const char *side, const char *uplo,
                     const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                     const magma_tally4DoubleComplex *alpha,
                     const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                     const magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb,
                     const magma_tally4DoubleComplex *beta,
                           magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc );

void blasf77_zsyr2k( const char *uplo, const char *trans,
                     const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                     const magma_tally4DoubleComplex *alpha,
                     const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                     const magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb,
                     const magma_tally4DoubleComplex *beta,
                           magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc );

void blasf77_zsyrk(  const char *uplo, const char *trans,
                     const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                     const magma_tally4DoubleComplex *alpha,
                     const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                     const magma_tally4DoubleComplex *beta,
                           magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc );

void blasf77_zrotg(  magma_tally4DoubleComplex *ca, const magma_tally4DoubleComplex *cb,
                     double *c, magma_tally4DoubleComplex *s );
                     
void blasf77_zrot(   const magma_tally4_int_t *n,
                     magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx,
                     magma_tally4DoubleComplex *y, const magma_tally4_int_t *incy,
                     const double *c, const magma_tally4DoubleComplex *s );
                     
void blasf77_zdrot(  const magma_tally4_int_t *n,
                     magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx,
                     magma_tally4DoubleComplex *y, const magma_tally4_int_t *incy,
                     const double *c, const double *s );

void blasf77_ztrmm(  const char *side, const char *uplo, const char *transa, const char *diag,
                     const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                     const magma_tally4DoubleComplex *alpha,
                     const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                           magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb );

void blasf77_ztrmv(  const char *uplo, const char *transa, const char *diag,
                     const magma_tally4_int_t *n,
                     const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                           magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx );

void blasf77_ztrsm(  const char *side, const char *uplo, const char *transa, const char *diag,
                     const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                     const magma_tally4DoubleComplex *alpha,
                     const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                           magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb );

void blasf77_ztrsv(  const char *uplo, const char *transa, const char *diag,
                     const magma_tally4_int_t *n,
                     const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                           magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4 wrappers around BLAS functions (alphabetical order)
    The Fortran interface for these is not portable, so we
    provide a C interface identical to the Fortran interface.
*/

double magma_tally4_cblas_dzasum(
    magma_tally4_int_t n,
    const magma_tally4DoubleComplex *x, magma_tally4_int_t incx );

double magma_tally4_cblas_dznrm2(
    magma_tally4_int_t n,
    const magma_tally4DoubleComplex *x, magma_tally4_int_t incx );

magma_tally4DoubleComplex magma_tally4_cblas_zdotc(
    magma_tally4_int_t n,
    const magma_tally4DoubleComplex *x, magma_tally4_int_t incx,
    const magma_tally4DoubleComplex *y, magma_tally4_int_t incy );

magma_tally4DoubleComplex magma_tally4_cblas_zdotu(
    magma_tally4_int_t n,
    const magma_tally4DoubleComplex *x, magma_tally4_int_t incx,
    const magma_tally4DoubleComplex *y, magma_tally4_int_t incy );


/*
 * LAPACK functions (alphabetical order)
 */
#ifdef REAL
void   lapackf77_dbdsdc( const char *uplo, const char *compq,
                         const magma_tally4_int_t *n,
                         double *d, double *e,
                         double *U,  const magma_tally4_int_t *ldu,
                         double *VT, const magma_tally4_int_t *ldvt,
                         double *Q, magma_tally4_int_t *IQ,
                         double *work, magma_tally4_int_t *iwork,
                         magma_tally4_int_t *info );
#endif

void   lapackf77_zbdsqr( const char *uplo,
                         const magma_tally4_int_t *n, const magma_tally4_int_t *ncvt, const magma_tally4_int_t *nru,  const magma_tally4_int_t *ncc,
                         double *d, double *e,
                         magma_tally4DoubleComplex *Vt, const magma_tally4_int_t *ldvt,
                         magma_tally4DoubleComplex *U, const magma_tally4_int_t *ldu,
                         magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc,
                         double *work,
                         magma_tally4_int_t *info );

void   lapackf77_zgebak( const char *job, const char *side,
                         const magma_tally4_int_t *n,
                         const magma_tally4_int_t *ilo, const magma_tally4_int_t *ihi,
                         const double *scale, const magma_tally4_int_t *m,
                         magma_tally4DoubleComplex *V, const magma_tally4_int_t *ldv,
                         magma_tally4_int_t *info );

void   lapackf77_zgebal( const char *job,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4_int_t *ilo, magma_tally4_int_t *ihi,
                         double *scale,
                         magma_tally4_int_t *info );

void   lapackf77_zgebd2( const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *d, double *e,
                         magma_tally4DoubleComplex *tauq,
                         magma_tally4DoubleComplex *taup,
                         magma_tally4DoubleComplex *work,
                         magma_tally4_int_t *info );

void   lapackf77_zgebrd( const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *d, double *e,
                         magma_tally4DoubleComplex *tauq,
                         magma_tally4DoubleComplex *taup,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zgbbrd( const char *vect, const magma_tally4_int_t *m,
                         const magma_tally4_int_t *n, const magma_tally4_int_t *ncc,
                         const magma_tally4_int_t *kl, const magma_tally4_int_t *ku,
                         magma_tally4DoubleComplex *Ab, const magma_tally4_int_t *ldab,
                         double *d, double *e,
                         magma_tally4DoubleComplex *Q, const magma_tally4_int_t *ldq,
                         magma_tally4DoubleComplex *PT, const magma_tally4_int_t *ldpt,
                         magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc,
                         magma_tally4DoubleComplex *work,
                         #ifdef COMPLEX
                         double *rwork,
                         #endif
                         magma_tally4_int_t *info );

void   lapackf77_zgeev(  const char *jobvl, const char *jobvr,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A,    const magma_tally4_int_t *lda,
                         #ifdef COMPLEX
                         magma_tally4DoubleComplex *w,
                         #else
                         double *wr, double *wi,
                         #endif
                         magma_tally4DoubleComplex *Vl,   const magma_tally4_int_t *ldvl,
                         magma_tally4DoubleComplex *Vr,   const magma_tally4_int_t *ldvr,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         #ifdef COMPLEX
                         double *rwork,
                         #endif
                         magma_tally4_int_t *info );

void   lapackf77_zgehd2( const magma_tally4_int_t *n,
                         const magma_tally4_int_t *ilo, const magma_tally4_int_t *ihi,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work,
                         magma_tally4_int_t *info );

void   lapackf77_zgehrd( const magma_tally4_int_t *n,
                         const magma_tally4_int_t *ilo, const magma_tally4_int_t *ihi,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zgelqf( const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zgels(  const char *trans,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *nrhs,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zgeqlf( const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zgeqp3( const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4_int_t *jpvt,
                         magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         #ifdef COMPLEX
                         double *rwork,
                         #endif
                         magma_tally4_int_t *info );

void   lapackf77_zgeqrf( const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zgesdd( const char *jobz,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *s,
                         magma_tally4DoubleComplex *U,  const magma_tally4_int_t *ldu,
                         magma_tally4DoubleComplex *Vt, const magma_tally4_int_t *ldvt,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         #ifdef COMPLEX
                         double *rwork,
                         #endif
                         magma_tally4_int_t *iwork, magma_tally4_int_t *info );

void   lapackf77_zgesv(  const magma_tally4_int_t *n, const magma_tally4_int_t *nrhs,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4_int_t *ipiv,
                         magma_tally4DoubleComplex *B,  const magma_tally4_int_t *ldb,
                         magma_tally4_int_t *info );

void   lapackf77_zgesvd( const char *jobu, const char *jobvt,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *s,
                         magma_tally4DoubleComplex *U,  const magma_tally4_int_t *ldu,
                         magma_tally4DoubleComplex *Vt, const magma_tally4_int_t *ldvt,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         #ifdef COMPLEX
                         double *rwork,
                         #endif
                         magma_tally4_int_t *info );

void   lapackf77_zgetrf( const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4_int_t *ipiv,
                         magma_tally4_int_t *info );

void   lapackf77_zgetri( const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4_int_t *ipiv,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zgetrs( const char *trans,
                         const magma_tally4_int_t *n, const magma_tally4_int_t *nrhs,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4_int_t *ipiv,
                         magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb,
                         magma_tally4_int_t *info );

void   lapackf77_zhetf2( const char*, magma_tally4_int_t*, 
                         magma_tally4DoubleComplex*, magma_tally4_int_t*, magma_tally4_int_t*, magma_tally4_int_t* );

void   lapackf77_zhetrs( const char *uplo,
                         const magma_tally4_int_t *n, const magma_tally4_int_t *nrhs,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4_int_t *ipiv,
                         magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb,
                         magma_tally4_int_t *info );

void   lapackf77_zhbtrd( const char *vect, const char *uplo,
                         const magma_tally4_int_t *n, const magma_tally4_int_t *kd,
                         magma_tally4DoubleComplex *Ab, const magma_tally4_int_t *ldab,
                         double *d, double *e,
                         magma_tally4DoubleComplex *Q, const magma_tally4_int_t *ldq,
                         magma_tally4DoubleComplex *work,
                         magma_tally4_int_t *info );

void   lapackf77_zheev(  const char *jobz, const char *uplo,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *w,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         #ifdef COMPLEX
                         double *rwork,
                         #endif
                         magma_tally4_int_t *info );

void   lapackf77_zheevd( const char *jobz, const char *uplo,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *w,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         #ifdef COMPLEX
                         double *rwork, const magma_tally4_int_t *lrwork,
                         #endif
                         magma_tally4_int_t *iwork, const magma_tally4_int_t *liwork,
                         magma_tally4_int_t *info );

void   lapackf77_zheevr( const char *jobz, const char *range, const char *uplo,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *vl, double *vu, 
                         magma_tally4_int_t *il, magma_tally4_int_t *iu, double *abstol, 
                         magma_tally4_int_t *m, double *w, magma_tally4DoubleComplex *z__, 
                         magma_tally4_int_t *ldz, magma_tally4_int_t *isuppz, 
                         magma_tally4DoubleComplex *work, magma_tally4_int_t *lwork, 
                         double *rwork, magma_tally4_int_t *lrwork, 
                         magma_tally4_int_t *iwork, magma_tally4_int_t *liwork, magma_tally4_int_t *info);

void   lapackf77_zheevx( const char *jobz, const char *range, const char *uplo,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *vl, double *vu,
                         magma_tally4_int_t *il, magma_tally4_int_t *iu, double *abstol,
                         magma_tally4_int_t *m, double *w, magma_tally4DoubleComplex *z__,
                         magma_tally4_int_t *ldz, magma_tally4DoubleComplex *work, magma_tally4_int_t *lwork,
                         double *rwork, magma_tally4_int_t *iwork, magma_tally4_int_t *ifail, magma_tally4_int_t *info);

void   lapackf77_zhegs2( const magma_tally4_int_t *itype, const char *uplo,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb,
                         magma_tally4_int_t *info );

void   lapackf77_zhegst( const magma_tally4_int_t *itype, const char *uplo,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb,
                         magma_tally4_int_t *info );

void   lapackf77_zhegvd( const magma_tally4_int_t *itype, const char *jobz, const char *uplo,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb,
                         double *w,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         #ifdef COMPLEX
                         double *rwork, const magma_tally4_int_t *lrwork,
                         #endif
                         magma_tally4_int_t *iwork, const magma_tally4_int_t *liwork,
                         magma_tally4_int_t *info );

void   lapackf77_zhesv( const char *uplo, 
                        const magma_tally4_int_t *n, const magma_tally4_int_t *nrhs,
                        magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda, magma_tally4_int_t *ipiv,
                        magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb,
                        magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                        magma_tally4_int_t *info );

void   lapackf77_zhetd2( const char *uplo,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *d, double *e,
                         magma_tally4DoubleComplex *tau,
                         magma_tally4_int_t *info );

void   lapackf77_zhetrd( const char *uplo,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *d, double *e,
                         magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zhetrf( const char *uplo,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4_int_t *ipiv,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zhseqr( const char *job, const char *compz,
                         const magma_tally4_int_t *n,
                         const magma_tally4_int_t *ilo, const magma_tally4_int_t *ihi,
                         magma_tally4DoubleComplex *H, const magma_tally4_int_t *ldh,
                         #ifdef COMPLEX
                         magma_tally4DoubleComplex *w,
                         #else
                         double *wr, double *wi,
                         #endif
                         magma_tally4DoubleComplex *Z, const magma_tally4_int_t *ldz,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zlabrd( const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *nb,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *d, double *e,
                         magma_tally4DoubleComplex *tauq,
                         magma_tally4DoubleComplex *taup,
                         magma_tally4DoubleComplex *X, const magma_tally4_int_t *ldx,
                         magma_tally4DoubleComplex *Y, const magma_tally4_int_t *ldy );

#ifdef COMPLEX
void   lapackf77_zlacgv( const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx );
#endif

#ifdef COMPLEX
void   lapackf77_zlacp2( const char *uplo,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         const double *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb );
#endif

void   lapackf77_zlacpy( const char *uplo,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb );

#ifdef COMPLEX
void   lapackf77_zlacrm( const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const double             *B, const magma_tally4_int_t *ldb,
                         magma_tally4DoubleComplex       *C, const magma_tally4_int_t *ldc,
                         double *rwork );
#endif

#ifdef COMPLEX
void   lapackf77_zladiv( magma_tally4DoubleComplex *ret_val,
                         magma_tally4DoubleComplex *x,
                         magma_tally4DoubleComplex *y );
#else // REAL
void   lapackf77_zladiv( const double *a, const double *b,
                         const double *c, const double *d,
                         double *p, double *q );
#endif

void   lapackf77_zlahef( const char *uplo,
                         const magma_tally4_int_t *n, const magma_tally4_int_t *kn,
                         magma_tally4_int_t *kb,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t lda,
                         magma_tally4_int_t *ipiv,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *ldwork,
                         magma_tally4_int_t *info );

double lapackf77_zlange( const char *norm,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *work );

double lapackf77_zlanhe( const char *norm, const char *uplo,
                         const magma_tally4_int_t *n,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *work );

double lapackf77_zlanht( const char *norm, const magma_tally4_int_t *n,
                         const double *d, const magma_tally4DoubleComplex *e );

double lapackf77_zlansy( const char *norm, const char *uplo,
                         const magma_tally4_int_t *n,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *work );

double lapackf77_zlantr( const char *norm, const char *uplo, const char *diag,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *work );

void   lapackf77_zlaqp2( magma_tally4_int_t *m, magma_tally4_int_t *n, magma_tally4_int_t *offset,
                         magma_tally4DoubleComplex *a, magma_tally4_int_t *lda, magma_tally4_int_t *jpvt,
                         magma_tally4DoubleComplex *tau,
                         double *vn1, double *vn2,
                         magma_tally4DoubleComplex *work );

#ifdef COMPLEX
void   lapackf77_zlarcm( const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         const double             *A, const magma_tally4_int_t *lda,
                         const magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb,
                         magma_tally4DoubleComplex       *C, const magma_tally4_int_t *ldc,
                         double *rwork );
#endif

void   lapackf77_zlarf(  const char *side, const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *v, const magma_tally4_int_t *incv,
                         magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc,
                         magma_tally4DoubleComplex *work );

void   lapackf77_zlarfb( const char *side, const char *trans, const char *direct, const char *storev,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                         const magma_tally4DoubleComplex *V, const magma_tally4_int_t *ldv,
                         const magma_tally4DoubleComplex *T, const magma_tally4_int_t *ldt,
                         magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *ldwork );

void   lapackf77_zlarfg( const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *alpha,
                         magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx,
                         magma_tally4DoubleComplex *tau );

void   lapackf77_zlarft( const char *direct, const char *storev,
                         const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                         magma_tally4DoubleComplex *V, const magma_tally4_int_t *ldv,
                         const magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *T, const magma_tally4_int_t *ldt );

void   lapackf77_zlarnv( const magma_tally4_int_t *idist, magma_tally4_int_t *iseed, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *x );

void   lapackf77_zlartg( magma_tally4DoubleComplex *F,
                         magma_tally4DoubleComplex *G,
                         double *cs,
                         magma_tally4DoubleComplex *SN,
                         magma_tally4DoubleComplex *R );

void   lapackf77_zlascl( const char *type,
                         const magma_tally4_int_t *kl, const magma_tally4_int_t *ku,
                         double *cfrom,
                         double *cto,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4_int_t *info );

void   lapackf77_zlaset( const char *uplo,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         const magma_tally4DoubleComplex *alpha,
                         const magma_tally4DoubleComplex *beta,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda );

void   lapackf77_zlaswp( const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4_int_t *k1, const magma_tally4_int_t *k2,
                         magma_tally4_int_t *ipiv,
                         const magma_tally4_int_t *incx );

void   lapackf77_zlatrd( const char *uplo,
                         const magma_tally4_int_t *n, const magma_tally4_int_t *nb,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *e,
                         magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *ldwork );

void   lapackf77_zlatrs( const char *uplo, const char *trans, const char *diag,
                         const char *normin,
                         const magma_tally4_int_t *n,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *x, double *scale,
                         double *cnorm, magma_tally4_int_t *info );

void   lapackf77_zlauum( const char *uplo,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4_int_t *info );

void   lapackf77_zlavhe( const char *uplo, const char *trans, const char *diag,
                         magma_tally4_int_t *n, magma_tally4_int_t *nrhs,
                         magma_tally4DoubleComplex *A, magma_tally4_int_t *lda,
                         magma_tally4_int_t *ipiv,
                         magma_tally4DoubleComplex *B, magma_tally4_int_t *ldb,
                         magma_tally4_int_t *info );

void   lapackf77_zposv(  const char *uplo,
                         const magma_tally4_int_t *n, const magma_tally4_int_t *nrhs,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *B,  const magma_tally4_int_t *ldb,
                         magma_tally4_int_t *info );

void   lapackf77_zpotrf( const char *uplo,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4_int_t *info );

void   lapackf77_zpotri( const char *uplo,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4_int_t *info );

void   lapackf77_zpotrs( const char *uplo,
                         const magma_tally4_int_t *n, const magma_tally4_int_t *nrhs,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *B, const magma_tally4_int_t *ldb,
                         magma_tally4_int_t *info );

void   lapackf77_zstedc( const char *compz,
                         const magma_tally4_int_t *n,
                         double *d, double *e,
                         magma_tally4DoubleComplex *Z, const magma_tally4_int_t *ldz,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         #ifdef COMPLEX
                         double *rwork, const magma_tally4_int_t *lrwork,
                         #endif
                         magma_tally4_int_t *iwork, const magma_tally4_int_t *liwork,
                         magma_tally4_int_t *info );

void   lapackf77_zstein( const magma_tally4_int_t *n,
                         const double *d, const double *e,
                         const magma_tally4_int_t *m,
                         const double *w,
                         const magma_tally4_int_t *iblock,
                         const magma_tally4_int_t *isplit,
                         magma_tally4DoubleComplex *Z, const magma_tally4_int_t *ldz,
                         double *work, magma_tally4_int_t *iwork, magma_tally4_int_t *ifailv,
                         magma_tally4_int_t *info );

void   lapackf77_zstemr( const char *jobz, const char *range,
                         const magma_tally4_int_t *n,
                         double *d, double *e,
                         const double *vl, const double *vu,
                         const magma_tally4_int_t *il, const magma_tally4_int_t *iu,
                         magma_tally4_int_t *m,
                         double *w,
                         magma_tally4DoubleComplex *Z, const magma_tally4_int_t *ldz,
                         const magma_tally4_int_t *nzc, magma_tally4_int_t *isuppz, magma_tally4_int_t *tryrac,
                         double *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *iwork, const magma_tally4_int_t *liwork,
                         magma_tally4_int_t *info );

void   lapackf77_zsteqr( const char *compz,
                         const magma_tally4_int_t *n,
                         double *d, double *e,
                         magma_tally4DoubleComplex *Z, const magma_tally4_int_t *ldz,
                         double *work,
                         magma_tally4_int_t *info );

#ifdef COMPLEX
void   lapackf77_zsymv(  const char *uplo,
                         const magma_tally4_int_t *n,
                         const magma_tally4DoubleComplex *alpha,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4DoubleComplex *x, const magma_tally4_int_t *incx,
                         const magma_tally4DoubleComplex *beta,
                               magma_tally4DoubleComplex *y, const magma_tally4_int_t *incy );
#endif

void   lapackf77_ztrevc( const char *side, const char *howmny,
                         magma_tally4_int_t *select, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *T,  const magma_tally4_int_t *ldt,
                         magma_tally4DoubleComplex *Vl, const magma_tally4_int_t *ldvl,
                         magma_tally4DoubleComplex *Vr, const magma_tally4_int_t *ldvr,
                         const magma_tally4_int_t *mm, magma_tally4_int_t *m,
                         magma_tally4DoubleComplex *work,
                         #ifdef COMPLEX
                         double *rwork,
                         #endif
                         magma_tally4_int_t *info );

void   lapackf77_ztrevc3( const char *side, const char *howmny,
                          magma_tally4_int_t *select, const magma_tally4_int_t *n,
                          magma_tally4DoubleComplex *T,  const magma_tally4_int_t *ldt,
                          magma_tally4DoubleComplex *VL, const magma_tally4_int_t *ldvl, 
                          magma_tally4DoubleComplex *VR, const magma_tally4_int_t *ldvr,
                          const magma_tally4_int_t *mm,
                          const magma_tally4_int_t *mout,
                          magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                          #ifdef COMPLEX
                          double *rwork,
                          #endif
                          magma_tally4_int_t *info );

void   lapackf77_ztrtri( const char *uplo, const char *diag,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4_int_t *info );

void   lapackf77_zung2r( const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work,
                         magma_tally4_int_t *info );

void   lapackf77_zungbr( const char *vect,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zunghr( const magma_tally4_int_t *n,
                         const magma_tally4_int_t *ilo, const magma_tally4_int_t *ihi,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zunglq( const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zungql( const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zungqr( const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zungtr( const char *uplo,
                         const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zunm2r( const char *side, const char *trans,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc,
                         magma_tally4DoubleComplex *work,
                         magma_tally4_int_t *info );

void   lapackf77_zunmbr( const char *vect, const char *side, const char *trans,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zunmlq( const char *side, const char *trans,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zunmql( const char *side, const char *trans,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zunmqr( const char *side, const char *trans,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

void   lapackf77_zunmtr( const char *side, const char *uplo, const char *trans,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         const magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         const magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         magma_tally4_int_t *info );

/*
 * Real precision extras
 */
void   lapackf77_dstebz( const char *range, const char *order,
                         const magma_tally4_int_t *n,
                         double *vl, double *vu,
                         magma_tally4_int_t *il, magma_tally4_int_t *iu,
                         double *abstol,
                         double *d, double *e,
                         const magma_tally4_int_t *m, const magma_tally4_int_t *nsplit,
                         double *w,
                         magma_tally4_int_t *iblock, magma_tally4_int_t *isplit,
                         double *work,
                         magma_tally4_int_t *iwork,
                         magma_tally4_int_t *info );

void   lapackf77_dlaln2( const magma_tally4_int_t *ltrans,
                         const magma_tally4_int_t *na, const magma_tally4_int_t *nw,
                         const double *smin, const double *ca,
                         const double *a,  const magma_tally4_int_t *lda,
                         const double *d1, const double *d2,
                         const double *b,  const magma_tally4_int_t *ldb,
                         const double *wr, const double *wi,
                         double *x, const magma_tally4_int_t *ldx,
                         double *scale, double *xnorm, magma_tally4_int_t *info );

double lapackf77_dlamc3( double *a, double *b );

void   lapackf77_dlamrg( magma_tally4_int_t *n1, magma_tally4_int_t *n2,
                         double *a,
                         magma_tally4_int_t *dtrd1, magma_tally4_int_t *dtrd2, magma_tally4_int_t *index );

double lapackf77_dlapy3( double *x, double *y, double *z );

void   lapackf77_dlaed2( magma_tally4_int_t *k, magma_tally4_int_t *n, magma_tally4_int_t *cutpnt,
                         double *d, double *q, magma_tally4_int_t *ldq, magma_tally4_int_t *indxq,
                         double *rho, double *z,
                         double *dlmda, double *w, double *q2,
                         magma_tally4_int_t *indx, magma_tally4_int_t *indxc, magma_tally4_int_t *indxp,
                         magma_tally4_int_t *coltyp, magma_tally4_int_t *info);

void   lapackf77_dlaed4( magma_tally4_int_t *n, magma_tally4_int_t *i,
                         double *d,
                         double *z,
                         double *delta,
                         double *rho,
                         double *dlam, magma_tally4_int_t *info );

void   lapackf77_dlasrt( const char *id, const magma_tally4_int_t *n, double *d, magma_tally4_int_t *info );

/*
 * Testing functions
 */
#ifdef COMPLEX
void   lapackf77_zbdt01( const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *kd,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *Q, const magma_tally4_int_t *ldq,
                         double *d, double *e,
                         magma_tally4DoubleComplex *Pt, const magma_tally4_int_t *ldpt,
                         magma_tally4DoubleComplex *work,
                         double *rwork,
                         double *resid );

void   lapackf77_zget22( const char *transa, const char *transe, const char *transw, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *E, const magma_tally4_int_t *lde,
                         magma_tally4DoubleComplex *w,
                         magma_tally4DoubleComplex *work,
                         double *rwork,
                         double *result );

void   lapackf77_zhet21( const magma_tally4_int_t *itype, const char *uplo,
                         const magma_tally4_int_t *n, const magma_tally4_int_t *kband,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *d, double *e,
                         magma_tally4DoubleComplex *U, const magma_tally4_int_t *ldu,
                         magma_tally4DoubleComplex *V, const magma_tally4_int_t *ldv,
                         magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work,
                         double *rwork,
                         double *result );

void   lapackf77_zhst01( const magma_tally4_int_t *n, const magma_tally4_int_t *ilo, const magma_tally4_int_t *ihi,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *H, const magma_tally4_int_t *ldh,
                         magma_tally4DoubleComplex *Q, const magma_tally4_int_t *ldq,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         double *rwork,
                         double *result );

void   lapackf77_zstt21( const magma_tally4_int_t *n, const magma_tally4_int_t *kband,
                         double *AD,
                         double *AE,
                         double *SD,
                         double *SE,
                         magma_tally4DoubleComplex *U, const magma_tally4_int_t *ldu,
                         magma_tally4DoubleComplex *work,
                         double *rwork,
                         double *result );

void   lapackf77_zunt01( const char *rowcol, const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *U, const magma_tally4_int_t *ldu,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         double *rwork,
                         double *resid );
#else
void   lapackf77_zbdt01( const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *kd,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *Q, const magma_tally4_int_t *ldq,
                         double *d, double *e,
                         magma_tally4DoubleComplex *Pt, const magma_tally4_int_t *ldpt,
                         magma_tally4DoubleComplex *work,
                         double *resid );

void   lapackf77_zget22( const char *transa, const char *transe, const char *transw, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *E, const magma_tally4_int_t *lde,
                         magma_tally4DoubleComplex *wr,
                         magma_tally4DoubleComplex *wi,
                         double *work,
                         double *result );

void   lapackf77_zhet21( magma_tally4_int_t *itype, const char *uplo, const magma_tally4_int_t *n, const magma_tally4_int_t *kband,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         double *d, double *e,
                         magma_tally4DoubleComplex *U, const magma_tally4_int_t *ldu,
                         magma_tally4DoubleComplex *V, const magma_tally4_int_t *ldv,
                         magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work,
                         double *result );

void   lapackf77_zhst01( const magma_tally4_int_t *n, const magma_tally4_int_t *ilo, const magma_tally4_int_t *ihi,
                         magma_tally4DoubleComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *H, const magma_tally4_int_t *ldh,
                         magma_tally4DoubleComplex *Q, const magma_tally4_int_t *ldq,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         double *result );

void   lapackf77_zstt21( const magma_tally4_int_t *n, const magma_tally4_int_t *kband,
                         double *AD,
                         double *AE,
                         double *SD,
                         double *SE,
                         magma_tally4DoubleComplex *U, const magma_tally4_int_t *ldu,
                         magma_tally4DoubleComplex *work,
                         double *result );

void   lapackf77_zunt01( const char *rowcol, const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *U, const magma_tally4_int_t *ldu,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         double *resid );
#endif

void   lapackf77_zlarfy( const char *uplo, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *V, const magma_tally4_int_t *incv,
                         magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc,
                         magma_tally4DoubleComplex *work );

void   lapackf77_zlarfx( const char *side, const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4DoubleComplex *V,
                         magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *C, const magma_tally4_int_t *ldc,
                         magma_tally4DoubleComplex *work );

double lapackf77_zqpt01( const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                         magma_tally4DoubleComplex *A,
                         magma_tally4DoubleComplex *Af, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *tau, magma_tally4_int_t *jpvt,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork );

void   lapackf77_zqrt02( const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                         magma_tally4DoubleComplex *A,
                         magma_tally4DoubleComplex *AF,
                         magma_tally4DoubleComplex *Q,
                         magma_tally4DoubleComplex *R, const magma_tally4_int_t *lda,
                         magma_tally4DoubleComplex *tau,
                         magma_tally4DoubleComplex *work, const magma_tally4_int_t *lwork,
                         double *rwork,
                         double *result );

void   lapackf77_zlatms( magma_tally4_int_t *m, magma_tally4_int_t *n,
                         const char *dist, magma_tally4_int_t *iseed, const char *sym, double *d,
                         magma_tally4_int_t *mode, const double *cond, const double *dmax,
                         magma_tally4_int_t *kl, magma_tally4_int_t *ku, const char *pack,
                         magma_tally4DoubleComplex *a, magma_tally4_int_t *lda, magma_tally4DoubleComplex *work, magma_tally4_int_t *info );

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif /* MAGMA_tally4_ZLAPACK_H */
