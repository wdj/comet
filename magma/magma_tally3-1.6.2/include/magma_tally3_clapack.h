/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally3_zlapack.h normal z -> c, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally3_CLAPACK_H
#define MAGMA_tally3_CLAPACK_H

#include "magma_tally3_types.h"
#include "magma_tally3_mangling.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- BLAS and LAPACK functions (alphabetical order)
*/
#define blasf77_icamax     FORTRAN_NAME( icamax, ICAMAX )
#define blasf77_caxpy      FORTRAN_NAME( caxpy,  CAXPY  )
#define blasf77_ccopy      FORTRAN_NAME( ccopy,  CCOPY  )
#define blasf77_cgemm      FORTRAN_NAME( cgemm,  CGEMM  )
#define blasf77_cgemv      FORTRAN_NAME( cgemv,  CGEMV  )
#define blasf77_cgerc      FORTRAN_NAME( cgerc,  CGERC  )
#define blasf77_cgeru      FORTRAN_NAME( cgeru,  CGERU  )
#define blasf77_chemm      FORTRAN_NAME( chemm,  CHEMM  )
#define blasf77_chemv      FORTRAN_NAME( chemv,  CHEMV  )
#define blasf77_cher       FORTRAN_NAME( cher,   CHER   )
#define blasf77_cher2      FORTRAN_NAME( cher2,  CHER2  )
#define blasf77_cher2k     FORTRAN_NAME( cher2k, CHER2K )
#define blasf77_cherk      FORTRAN_NAME( cherk,  CHERK  )
#define blasf77_cscal      FORTRAN_NAME( cscal,  CSCAL  )
#define blasf77_csscal     FORTRAN_NAME( csscal, CSSCAL )
#define blasf77_cswap      FORTRAN_NAME( cswap,  CSWAP  )
#define blasf77_csymm      FORTRAN_NAME( csymm,  CSYMM  )
#define blasf77_csyr2k     FORTRAN_NAME( csyr2k, CSYR2K )
#define blasf77_csyrk      FORTRAN_NAME( csyrk,  CSYRK  )
#define blasf77_crotg      FORTRAN_NAME( crotg,  CROTG  )
#define blasf77_crot       FORTRAN_NAME( crot,   CROT   )
#define blasf77_csrot      FORTRAN_NAME( csrot,  CSROT  )
#define blasf77_ctrmm      FORTRAN_NAME( ctrmm,  CTRMM  )
#define blasf77_ctrmv      FORTRAN_NAME( ctrmv,  CTRMV  )
#define blasf77_ctrsm      FORTRAN_NAME( ctrsm,  CTRSM  )
#define blasf77_ctrsv      FORTRAN_NAME( ctrsv,  CTRSV  )

#define lapackf77_slaed2   FORTRAN_NAME( slaed2, SLAED2 )
#define lapackf77_slaed4   FORTRAN_NAME( slaed4, SLAED4 )
#define lapackf77_slaln2   FORTRAN_NAME( slaln2, SLALN2 )
#define lapackf77_slamc3   FORTRAN_NAME( slamc3, SLAMC3 )
#define lapackf77_slamrg   FORTRAN_NAME( slamrg, SLAMRG )
#define lapackf77_slasrt   FORTRAN_NAME( slasrt, SLASRT )
#define lapackf77_sstebz   FORTRAN_NAME( sstebz, SSTEBZ )

#define lapackf77_sbdsdc   FORTRAN_NAME( sbdsdc, SBDSDC )
#define lapackf77_cbdsqr   FORTRAN_NAME( cbdsqr, CBDSQR )
#define lapackf77_cgebak   FORTRAN_NAME( cgebak, CGEBAK )
#define lapackf77_cgebal   FORTRAN_NAME( cgebal, CGEBAL )
#define lapackf77_cgebd2   FORTRAN_NAME( cgebd2, CGEBD2 )
#define lapackf77_cgebrd   FORTRAN_NAME( cgebrd, CGEBRD )
#define lapackf77_cgbbrd   FORTRAN_NAME( cgbbrd, CGBBRD )
#define lapackf77_cgeev    FORTRAN_NAME( cgeev,  CGEEV  )
#define lapackf77_cgehd2   FORTRAN_NAME( cgehd2, CGEHD2 )
#define lapackf77_cgehrd   FORTRAN_NAME( cgehrd, CGEHRD )
#define lapackf77_cgelqf   FORTRAN_NAME( cgelqf, CGELQF )
#define lapackf77_cgels    FORTRAN_NAME( cgels,  CGELS  )
#define lapackf77_cgeqlf   FORTRAN_NAME( cgeqlf, CGEQLF )
#define lapackf77_cgeqp3   FORTRAN_NAME( cgeqp3, CGEQP3 )
#define lapackf77_cgeqrf   FORTRAN_NAME( cgeqrf, CGEQRF )
#define lapackf77_cgesdd   FORTRAN_NAME( cgesdd, CGESDD )
#define lapackf77_cgesv    FORTRAN_NAME( cgesv,  CGESV  )
#define lapackf77_cgesvd   FORTRAN_NAME( cgesvd, CGESVD )
#define lapackf77_cgetrf   FORTRAN_NAME( cgetrf, CGETRF )
#define lapackf77_cgetri   FORTRAN_NAME( cgetri, CGETRI )
#define lapackf77_cgetrs   FORTRAN_NAME( cgetrs, CGETRS )
#define lapackf77_chetf2   FORTRAN_NAME( chetf2, CHETF2 )
#define lapackf77_chetrs   FORTRAN_NAME( chetrs, CHETRS )
#define lapackf77_chbtrd   FORTRAN_NAME( chbtrd, CHBTRD )
#define lapackf77_cheev    FORTRAN_NAME( cheev,  CHEEV  )
#define lapackf77_cheevd   FORTRAN_NAME( cheevd, CHEEVD )
#define lapackf77_cheevr   FORTRAN_NAME( cheevr, CHEEVR )
#define lapackf77_cheevx   FORTRAN_NAME( cheevx, CHEEVX )
#define lapackf77_chegs2   FORTRAN_NAME( chegs2, CHEGS2 )
#define lapackf77_chegst   FORTRAN_NAME( chegst, CHEGST )
#define lapackf77_chegvd   FORTRAN_NAME( chegvd, CHEGVD )
#define lapackf77_chesv    FORTRAN_NAME( chesv,  CHESV  )
#define lapackf77_chetd2   FORTRAN_NAME( chetd2, CHETD2 )
#define lapackf77_chetrd   FORTRAN_NAME( chetrd, CHETRD )
#define lapackf77_chetrf   FORTRAN_NAME( chetrf, CHETRF )
#define lapackf77_chseqr   FORTRAN_NAME( chseqr, CHSEQR )
#define lapackf77_clabrd   FORTRAN_NAME( clabrd, CLABRD )
#define lapackf77_clacgv   FORTRAN_NAME( clacgv, CLACGV )
#define lapackf77_clacp2   FORTRAN_NAME( clacp2, CLACP2 )
#define lapackf77_clacpy   FORTRAN_NAME( clacpy, CLACPY )
#define lapackf77_clacrm   FORTRAN_NAME( clacrm, CLACRM )
#define lapackf77_cladiv   FORTRAN_NAME( cladiv, CLADIV )
#define lapackf77_clahef   FORTRAN_NAME( clahef, CLAHEF )
#define lapackf77_clange   FORTRAN_NAME( clange, CLANGE )
#define lapackf77_clanhe   FORTRAN_NAME( clanhe, CLANHE )
#define lapackf77_clanht   FORTRAN_NAME( clanht, CLANHT )
#define lapackf77_clansy   FORTRAN_NAME( clansy, CLANSY )
#define lapackf77_clantr   FORTRAN_NAME( clantr, CLANTR )
#define lapackf77_slapy3   FORTRAN_NAME( slapy3, SLAPY3 )
#define lapackf77_claqp2   FORTRAN_NAME( claqp2, CLAQP2 )
#define lapackf77_clarcm   FORTRAN_NAME( clarcm, CLARCM )
#define lapackf77_clarf    FORTRAN_NAME( clarf,  CLARF  )
#define lapackf77_clarfb   FORTRAN_NAME( clarfb, CLARFB )
#define lapackf77_clarfg   FORTRAN_NAME( clarfg, CLARFG )
#define lapackf77_clarft   FORTRAN_NAME( clarft, CLARFT )
#define lapackf77_clarnv   FORTRAN_NAME( clarnv, CLARNV )
#define lapackf77_clartg   FORTRAN_NAME( clartg, CLARTG )
#define lapackf77_clascl   FORTRAN_NAME( clascl, CLASCL )
#define lapackf77_claset   FORTRAN_NAME( claset, CLASET )
#define lapackf77_claswp   FORTRAN_NAME( claswp, CLASWP )
#define lapackf77_clatrd   FORTRAN_NAME( clatrd, CLATRD )
#define lapackf77_clatrs   FORTRAN_NAME( clatrs, CLATRS )
#define lapackf77_clauum   FORTRAN_NAME( clauum, CLAUUM )
#define lapackf77_clavhe   FORTRAN_NAME( clavhe, CLAVHE )
#define lapackf77_cposv    FORTRAN_NAME( cposv,  CPOSV  )
#define lapackf77_cpotrf   FORTRAN_NAME( cpotrf, CPOTRF )
#define lapackf77_cpotri   FORTRAN_NAME( cpotri, CPOTRI )
#define lapackf77_cpotrs   FORTRAN_NAME( cpotrs, CPOTRS )
#define lapackf77_cstedc   FORTRAN_NAME( cstedc, CSTEDC )
#define lapackf77_cstein   FORTRAN_NAME( cstein, CSTEIN )
#define lapackf77_cstemr   FORTRAN_NAME( cstemr, CSTEMR )
#define lapackf77_csteqr   FORTRAN_NAME( csteqr, CSTEQR )
#define lapackf77_csymv    FORTRAN_NAME( csymv,  CSYMV  )
#define lapackf77_ctrevc   FORTRAN_NAME( ctrevc, CTREVC )
#define lapackf77_ctrevc3  FORTRAN_NAME( ctrevc3, CTREVC3 )
#define lapackf77_ctrtri   FORTRAN_NAME( ctrtri, CTRTRI )
#define lapackf77_cung2r   FORTRAN_NAME( cung2r, CUNG2R )
#define lapackf77_cungbr   FORTRAN_NAME( cungbr, CUNGBR )
#define lapackf77_cunghr   FORTRAN_NAME( cunghr, CUNGHR )
#define lapackf77_cunglq   FORTRAN_NAME( cunglq, CUNGLQ )
#define lapackf77_cungql   FORTRAN_NAME( cungql, CUNGQL )
#define lapackf77_cungqr   FORTRAN_NAME( cungqr, CUNGQR )
#define lapackf77_cungtr   FORTRAN_NAME( cungtr, CUNGTR )
#define lapackf77_cunm2r   FORTRAN_NAME( cunm2r, CUNM2R )
#define lapackf77_cunmbr   FORTRAN_NAME( cunmbr, CUNMBR )
#define lapackf77_cunmlq   FORTRAN_NAME( cunmlq, CUNMLQ )
#define lapackf77_cunmql   FORTRAN_NAME( cunmql, CUNMQL )
#define lapackf77_cunmqr   FORTRAN_NAME( cunmqr, CUNMQR )
#define lapackf77_cunmtr   FORTRAN_NAME( cunmtr, CUNMTR )

/* testing functions (alphabetical order) */
#define lapackf77_cbdt01   FORTRAN_NAME( cbdt01, CBDT01 )
#define lapackf77_cget22   FORTRAN_NAME( cget22, CGET22 )
#define lapackf77_chet21   FORTRAN_NAME( chet21, CHET21 )
#define lapackf77_chst01   FORTRAN_NAME( chst01, CHST01 )
#define lapackf77_clarfx   FORTRAN_NAME( clarfx, CLARFX )
#define lapackf77_clarfy   FORTRAN_NAME( clarfy, CLARFY )
#define lapackf77_clatms   FORTRAN_NAME( clatms, CLATMS )
#define lapackf77_cqpt01   FORTRAN_NAME( cqpt01, CQPT01 )
#define lapackf77_cqrt02   FORTRAN_NAME( cqrt02, CQRT02 )
#define lapackf77_cstt21   FORTRAN_NAME( cstt21, CSTT21 )
#define lapackf77_cunt01   FORTRAN_NAME( cunt01, CUNT01 )

/*
 * BLAS functions (alphabetical order)
 */
magma_tally3_int_t blasf77_icamax(
                     const magma_tally3_int_t *n,
                     const magma_tally3FloatComplex *x, const magma_tally3_int_t *incx );

void blasf77_caxpy(  const magma_tally3_int_t *n,
                     const magma_tally3FloatComplex *alpha,
                     const magma_tally3FloatComplex *x, const magma_tally3_int_t *incx,
                           magma_tally3FloatComplex *y, const magma_tally3_int_t *incy );

void blasf77_ccopy(  const magma_tally3_int_t *n,
                     const magma_tally3FloatComplex *x, const magma_tally3_int_t *incx,
                           magma_tally3FloatComplex *y, const magma_tally3_int_t *incy );

void blasf77_cgemm(  const char *transa, const char *transb,
                     const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                     const magma_tally3FloatComplex *alpha,
                     const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                     const magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb,
                     const magma_tally3FloatComplex *beta,
                           magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc );

void blasf77_cgemv(  const char *transa,
                     const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                     const magma_tally3FloatComplex *alpha,
                     const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                     const magma_tally3FloatComplex *x, const magma_tally3_int_t *incx,
                     const magma_tally3FloatComplex *beta,
                           magma_tally3FloatComplex *y, const magma_tally3_int_t *incy );

void blasf77_cgerc(  const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                     const magma_tally3FloatComplex *alpha,
                     const magma_tally3FloatComplex *x, const magma_tally3_int_t *incx,
                     const magma_tally3FloatComplex *y, const magma_tally3_int_t *incy,
                           magma_tally3FloatComplex *A, const magma_tally3_int_t *lda );

void blasf77_cgeru(  const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                     const magma_tally3FloatComplex *alpha,
                     const magma_tally3FloatComplex *x, const magma_tally3_int_t *incx,
                     const magma_tally3FloatComplex *y, const magma_tally3_int_t *incy,
                           magma_tally3FloatComplex *A, const magma_tally3_int_t *lda );

void blasf77_chemm(  const char *side, const char *uplo,
                     const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                     const magma_tally3FloatComplex *alpha,
                     const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                     const magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb,
                     const magma_tally3FloatComplex *beta,
                           magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc );

void blasf77_chemv(  const char *uplo,
                     const magma_tally3_int_t *n,
                     const magma_tally3FloatComplex *alpha,
                     const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                     const magma_tally3FloatComplex *x, const magma_tally3_int_t *incx,
                     const magma_tally3FloatComplex *beta,
                           magma_tally3FloatComplex *y, const magma_tally3_int_t *incy );

void blasf77_cher(   const char *uplo,
                     const magma_tally3_int_t *n,
                     const float *alpha,
                     const magma_tally3FloatComplex *x, const magma_tally3_int_t *incx,
                           magma_tally3FloatComplex *A, const magma_tally3_int_t *lda );

void blasf77_cher2(  const char *uplo,
                     const magma_tally3_int_t *n,
                     const magma_tally3FloatComplex *alpha,
                     const magma_tally3FloatComplex *x, const magma_tally3_int_t *incx,
                     const magma_tally3FloatComplex *y, const magma_tally3_int_t *incy,
                           magma_tally3FloatComplex *A, const magma_tally3_int_t *lda );

void blasf77_cher2k( const char *uplo, const char *trans,
                     const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                     const magma_tally3FloatComplex *alpha,
                     const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                     const magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb,
                     const float *beta,
                           magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc );

void blasf77_cherk(  const char *uplo, const char *trans,
                     const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                     const float *alpha,
                     const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                     const float *beta,
                           magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc );

void blasf77_cscal(  const magma_tally3_int_t *n,
                     const magma_tally3FloatComplex *alpha,
                           magma_tally3FloatComplex *x, const magma_tally3_int_t *incx );

void blasf77_csscal( const magma_tally3_int_t *n,
                     const float *alpha,
                           magma_tally3FloatComplex *x, const magma_tally3_int_t *incx );

void blasf77_cswap(  const magma_tally3_int_t *n,
                     magma_tally3FloatComplex *x, const magma_tally3_int_t *incx,
                     magma_tally3FloatComplex *y, const magma_tally3_int_t *incy );

/* complex-symmetric (non-Hermitian) routines */
void blasf77_csymm(  const char *side, const char *uplo,
                     const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                     const magma_tally3FloatComplex *alpha,
                     const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                     const magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb,
                     const magma_tally3FloatComplex *beta,
                           magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc );

void blasf77_csyr2k( const char *uplo, const char *trans,
                     const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                     const magma_tally3FloatComplex *alpha,
                     const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                     const magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb,
                     const magma_tally3FloatComplex *beta,
                           magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc );

void blasf77_csyrk(  const char *uplo, const char *trans,
                     const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                     const magma_tally3FloatComplex *alpha,
                     const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                     const magma_tally3FloatComplex *beta,
                           magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc );

void blasf77_crotg(  magma_tally3FloatComplex *ca, const magma_tally3FloatComplex *cb,
                     float *c, magma_tally3FloatComplex *s );
                     
void blasf77_crot(   const magma_tally3_int_t *n,
                     magma_tally3FloatComplex *x, const magma_tally3_int_t *incx,
                     magma_tally3FloatComplex *y, const magma_tally3_int_t *incy,
                     const float *c, const magma_tally3FloatComplex *s );
                     
void blasf77_csrot(  const magma_tally3_int_t *n,
                     magma_tally3FloatComplex *x, const magma_tally3_int_t *incx,
                     magma_tally3FloatComplex *y, const magma_tally3_int_t *incy,
                     const float *c, const float *s );

void blasf77_ctrmm(  const char *side, const char *uplo, const char *transa, const char *diag,
                     const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                     const magma_tally3FloatComplex *alpha,
                     const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                           magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb );

void blasf77_ctrmv(  const char *uplo, const char *transa, const char *diag,
                     const magma_tally3_int_t *n,
                     const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                           magma_tally3FloatComplex *x, const magma_tally3_int_t *incx );

void blasf77_ctrsm(  const char *side, const char *uplo, const char *transa, const char *diag,
                     const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                     const magma_tally3FloatComplex *alpha,
                     const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                           magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb );

void blasf77_ctrsv(  const char *uplo, const char *transa, const char *diag,
                     const magma_tally3_int_t *n,
                     const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                           magma_tally3FloatComplex *x, const magma_tally3_int_t *incx );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3 wrappers around BLAS functions (alphabetical order)
    The Fortran interface for these is not portable, so we
    provide a C interface identical to the Fortran interface.
*/

float magma_tally3_cblas_scasum(
    magma_tally3_int_t n,
    const magma_tally3FloatComplex *x, magma_tally3_int_t incx );

float magma_tally3_cblas_scnrm2(
    magma_tally3_int_t n,
    const magma_tally3FloatComplex *x, magma_tally3_int_t incx );

magma_tally3FloatComplex magma_tally3_cblas_cdotc(
    magma_tally3_int_t n,
    const magma_tally3FloatComplex *x, magma_tally3_int_t incx,
    const magma_tally3FloatComplex *y, magma_tally3_int_t incy );

magma_tally3FloatComplex magma_tally3_cblas_cdotu(
    magma_tally3_int_t n,
    const magma_tally3FloatComplex *x, magma_tally3_int_t incx,
    const magma_tally3FloatComplex *y, magma_tally3_int_t incy );


/*
 * LAPACK functions (alphabetical order)
 */
#ifdef REAL
void   lapackf77_sbdsdc( const char *uplo, const char *compq,
                         const magma_tally3_int_t *n,
                         float *d, float *e,
                         float *U,  const magma_tally3_int_t *ldu,
                         float *VT, const magma_tally3_int_t *ldvt,
                         float *Q, magma_tally3_int_t *IQ,
                         float *work, magma_tally3_int_t *iwork,
                         magma_tally3_int_t *info );
#endif

void   lapackf77_cbdsqr( const char *uplo,
                         const magma_tally3_int_t *n, const magma_tally3_int_t *ncvt, const magma_tally3_int_t *nru,  const magma_tally3_int_t *ncc,
                         float *d, float *e,
                         magma_tally3FloatComplex *Vt, const magma_tally3_int_t *ldvt,
                         magma_tally3FloatComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc,
                         float *work,
                         magma_tally3_int_t *info );

void   lapackf77_cgebak( const char *job, const char *side,
                         const magma_tally3_int_t *n,
                         const magma_tally3_int_t *ilo, const magma_tally3_int_t *ihi,
                         const float *scale, const magma_tally3_int_t *m,
                         magma_tally3FloatComplex *V, const magma_tally3_int_t *ldv,
                         magma_tally3_int_t *info );

void   lapackf77_cgebal( const char *job,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3_int_t *ilo, magma_tally3_int_t *ihi,
                         float *scale,
                         magma_tally3_int_t *info );

void   lapackf77_cgebd2( const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *d, float *e,
                         magma_tally3FloatComplex *tauq,
                         magma_tally3FloatComplex *taup,
                         magma_tally3FloatComplex *work,
                         magma_tally3_int_t *info );

void   lapackf77_cgebrd( const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *d, float *e,
                         magma_tally3FloatComplex *tauq,
                         magma_tally3FloatComplex *taup,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cgbbrd( const char *vect, const magma_tally3_int_t *m,
                         const magma_tally3_int_t *n, const magma_tally3_int_t *ncc,
                         const magma_tally3_int_t *kl, const magma_tally3_int_t *ku,
                         magma_tally3FloatComplex *Ab, const magma_tally3_int_t *ldab,
                         float *d, float *e,
                         magma_tally3FloatComplex *Q, const magma_tally3_int_t *ldq,
                         magma_tally3FloatComplex *PT, const magma_tally3_int_t *ldpt,
                         magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc,
                         magma_tally3FloatComplex *work,
                         #ifdef COMPLEX
                         float *rwork,
                         #endif
                         magma_tally3_int_t *info );

void   lapackf77_cgeev(  const char *jobvl, const char *jobvr,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A,    const magma_tally3_int_t *lda,
                         #ifdef COMPLEX
                         magma_tally3FloatComplex *w,
                         #else
                         float *wr, float *wi,
                         #endif
                         magma_tally3FloatComplex *Vl,   const magma_tally3_int_t *ldvl,
                         magma_tally3FloatComplex *Vr,   const magma_tally3_int_t *ldvr,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         #ifdef COMPLEX
                         float *rwork,
                         #endif
                         magma_tally3_int_t *info );

void   lapackf77_cgehd2( const magma_tally3_int_t *n,
                         const magma_tally3_int_t *ilo, const magma_tally3_int_t *ihi,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work,
                         magma_tally3_int_t *info );

void   lapackf77_cgehrd( const magma_tally3_int_t *n,
                         const magma_tally3_int_t *ilo, const magma_tally3_int_t *ihi,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cgelqf( const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cgels(  const char *trans,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *nrhs,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cgeqlf( const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cgeqp3( const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3_int_t *jpvt,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         #ifdef COMPLEX
                         float *rwork,
                         #endif
                         magma_tally3_int_t *info );

void   lapackf77_cgeqrf( const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cgesdd( const char *jobz,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *s,
                         magma_tally3FloatComplex *U,  const magma_tally3_int_t *ldu,
                         magma_tally3FloatComplex *Vt, const magma_tally3_int_t *ldvt,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         #ifdef COMPLEX
                         float *rwork,
                         #endif
                         magma_tally3_int_t *iwork, magma_tally3_int_t *info );

void   lapackf77_cgesv(  const magma_tally3_int_t *n, const magma_tally3_int_t *nrhs,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3_int_t *ipiv,
                         magma_tally3FloatComplex *B,  const magma_tally3_int_t *ldb,
                         magma_tally3_int_t *info );

void   lapackf77_cgesvd( const char *jobu, const char *jobvt,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *s,
                         magma_tally3FloatComplex *U,  const magma_tally3_int_t *ldu,
                         magma_tally3FloatComplex *Vt, const magma_tally3_int_t *ldvt,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         #ifdef COMPLEX
                         float *rwork,
                         #endif
                         magma_tally3_int_t *info );

void   lapackf77_cgetrf( const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3_int_t *ipiv,
                         magma_tally3_int_t *info );

void   lapackf77_cgetri( const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3_int_t *ipiv,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cgetrs( const char *trans,
                         const magma_tally3_int_t *n, const magma_tally3_int_t *nrhs,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3_int_t *ipiv,
                         magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb,
                         magma_tally3_int_t *info );

void   lapackf77_chetf2( const char*, magma_tally3_int_t*, 
                         magma_tally3FloatComplex*, magma_tally3_int_t*, magma_tally3_int_t*, magma_tally3_int_t* );

void   lapackf77_chetrs( const char *uplo,
                         const magma_tally3_int_t *n, const magma_tally3_int_t *nrhs,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3_int_t *ipiv,
                         magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb,
                         magma_tally3_int_t *info );

void   lapackf77_chbtrd( const char *vect, const char *uplo,
                         const magma_tally3_int_t *n, const magma_tally3_int_t *kd,
                         magma_tally3FloatComplex *Ab, const magma_tally3_int_t *ldab,
                         float *d, float *e,
                         magma_tally3FloatComplex *Q, const magma_tally3_int_t *ldq,
                         magma_tally3FloatComplex *work,
                         magma_tally3_int_t *info );

void   lapackf77_cheev(  const char *jobz, const char *uplo,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *w,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         #ifdef COMPLEX
                         float *rwork,
                         #endif
                         magma_tally3_int_t *info );

void   lapackf77_cheevd( const char *jobz, const char *uplo,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *w,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         #ifdef COMPLEX
                         float *rwork, const magma_tally3_int_t *lrwork,
                         #endif
                         magma_tally3_int_t *iwork, const magma_tally3_int_t *liwork,
                         magma_tally3_int_t *info );

void   lapackf77_cheevr( const char *jobz, const char *range, const char *uplo,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *vl, float *vu, 
                         magma_tally3_int_t *il, magma_tally3_int_t *iu, float *abstol, 
                         magma_tally3_int_t *m, float *w, magma_tally3FloatComplex *z__, 
                         magma_tally3_int_t *ldz, magma_tally3_int_t *isuppz, 
                         magma_tally3FloatComplex *work, magma_tally3_int_t *lwork, 
                         float *rwork, magma_tally3_int_t *lrwork, 
                         magma_tally3_int_t *iwork, magma_tally3_int_t *liwork, magma_tally3_int_t *info);

void   lapackf77_cheevx( const char *jobz, const char *range, const char *uplo,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *vl, float *vu,
                         magma_tally3_int_t *il, magma_tally3_int_t *iu, float *abstol,
                         magma_tally3_int_t *m, float *w, magma_tally3FloatComplex *z__,
                         magma_tally3_int_t *ldz, magma_tally3FloatComplex *work, magma_tally3_int_t *lwork,
                         float *rwork, magma_tally3_int_t *iwork, magma_tally3_int_t *ifail, magma_tally3_int_t *info);

void   lapackf77_chegs2( const magma_tally3_int_t *itype, const char *uplo,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb,
                         magma_tally3_int_t *info );

void   lapackf77_chegst( const magma_tally3_int_t *itype, const char *uplo,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb,
                         magma_tally3_int_t *info );

void   lapackf77_chegvd( const magma_tally3_int_t *itype, const char *jobz, const char *uplo,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb,
                         float *w,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         #ifdef COMPLEX
                         float *rwork, const magma_tally3_int_t *lrwork,
                         #endif
                         magma_tally3_int_t *iwork, const magma_tally3_int_t *liwork,
                         magma_tally3_int_t *info );

void   lapackf77_chesv( const char *uplo, 
                        const magma_tally3_int_t *n, const magma_tally3_int_t *nrhs,
                        magma_tally3FloatComplex *A, const magma_tally3_int_t *lda, magma_tally3_int_t *ipiv,
                        magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb,
                        magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                        magma_tally3_int_t *info );

void   lapackf77_chetd2( const char *uplo,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *d, float *e,
                         magma_tally3FloatComplex *tau,
                         magma_tally3_int_t *info );

void   lapackf77_chetrd( const char *uplo,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *d, float *e,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_chetrf( const char *uplo,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3_int_t *ipiv,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_chseqr( const char *job, const char *compz,
                         const magma_tally3_int_t *n,
                         const magma_tally3_int_t *ilo, const magma_tally3_int_t *ihi,
                         magma_tally3FloatComplex *H, const magma_tally3_int_t *ldh,
                         #ifdef COMPLEX
                         magma_tally3FloatComplex *w,
                         #else
                         float *wr, float *wi,
                         #endif
                         magma_tally3FloatComplex *Z, const magma_tally3_int_t *ldz,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_clabrd( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *nb,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *d, float *e,
                         magma_tally3FloatComplex *tauq,
                         magma_tally3FloatComplex *taup,
                         magma_tally3FloatComplex *X, const magma_tally3_int_t *ldx,
                         magma_tally3FloatComplex *Y, const magma_tally3_int_t *ldy );

#ifdef COMPLEX
void   lapackf77_clacgv( const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *x, const magma_tally3_int_t *incx );
#endif

#ifdef COMPLEX
void   lapackf77_clacp2( const char *uplo,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         const float *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb );
#endif

void   lapackf77_clacpy( const char *uplo,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb );

#ifdef COMPLEX
void   lapackf77_clacrm( const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const float             *B, const magma_tally3_int_t *ldb,
                         magma_tally3FloatComplex       *C, const magma_tally3_int_t *ldc,
                         float *rwork );
#endif

#ifdef COMPLEX
void   lapackf77_cladiv( magma_tally3FloatComplex *ret_val,
                         magma_tally3FloatComplex *x,
                         magma_tally3FloatComplex *y );
#else // REAL
void   lapackf77_cladiv( const float *a, const float *b,
                         const float *c, const float *d,
                         float *p, float *q );
#endif

void   lapackf77_clahef( const char *uplo,
                         const magma_tally3_int_t *n, const magma_tally3_int_t *kn,
                         magma_tally3_int_t *kb,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t lda,
                         magma_tally3_int_t *ipiv,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *ldwork,
                         magma_tally3_int_t *info );

float lapackf77_clange( const char *norm,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *work );

float lapackf77_clanhe( const char *norm, const char *uplo,
                         const magma_tally3_int_t *n,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *work );

float lapackf77_clanht( const char *norm, const magma_tally3_int_t *n,
                         const float *d, const magma_tally3FloatComplex *e );

float lapackf77_clansy( const char *norm, const char *uplo,
                         const magma_tally3_int_t *n,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *work );

float lapackf77_clantr( const char *norm, const char *uplo, const char *diag,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *work );

void   lapackf77_claqp2( magma_tally3_int_t *m, magma_tally3_int_t *n, magma_tally3_int_t *offset,
                         magma_tally3FloatComplex *a, magma_tally3_int_t *lda, magma_tally3_int_t *jpvt,
                         magma_tally3FloatComplex *tau,
                         float *vn1, float *vn2,
                         magma_tally3FloatComplex *work );

#ifdef COMPLEX
void   lapackf77_clarcm( const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         const float             *A, const magma_tally3_int_t *lda,
                         const magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb,
                         magma_tally3FloatComplex       *C, const magma_tally3_int_t *ldc,
                         float *rwork );
#endif

void   lapackf77_clarf(  const char *side, const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *v, const magma_tally3_int_t *incv,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc,
                         magma_tally3FloatComplex *work );

void   lapackf77_clarfb( const char *side, const char *trans, const char *direct, const char *storev,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         const magma_tally3FloatComplex *V, const magma_tally3_int_t *ldv,
                         const magma_tally3FloatComplex *T, const magma_tally3_int_t *ldt,
                         magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *ldwork );

void   lapackf77_clarfg( const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *alpha,
                         magma_tally3FloatComplex *x, const magma_tally3_int_t *incx,
                         magma_tally3FloatComplex *tau );

void   lapackf77_clarft( const char *direct, const char *storev,
                         const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         magma_tally3FloatComplex *V, const magma_tally3_int_t *ldv,
                         const magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *T, const magma_tally3_int_t *ldt );

void   lapackf77_clarnv( const magma_tally3_int_t *idist, magma_tally3_int_t *iseed, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *x );

void   lapackf77_clartg( magma_tally3FloatComplex *F,
                         magma_tally3FloatComplex *G,
                         float *cs,
                         magma_tally3FloatComplex *SN,
                         magma_tally3FloatComplex *R );

void   lapackf77_clascl( const char *type,
                         const magma_tally3_int_t *kl, const magma_tally3_int_t *ku,
                         float *cfrom,
                         float *cto,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3_int_t *info );

void   lapackf77_claset( const char *uplo,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         const magma_tally3FloatComplex *alpha,
                         const magma_tally3FloatComplex *beta,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda );

void   lapackf77_claswp( const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3_int_t *k1, const magma_tally3_int_t *k2,
                         magma_tally3_int_t *ipiv,
                         const magma_tally3_int_t *incx );

void   lapackf77_clatrd( const char *uplo,
                         const magma_tally3_int_t *n, const magma_tally3_int_t *nb,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *e,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *ldwork );

void   lapackf77_clatrs( const char *uplo, const char *trans, const char *diag,
                         const char *normin,
                         const magma_tally3_int_t *n,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *x, float *scale,
                         float *cnorm, magma_tally3_int_t *info );

void   lapackf77_clauum( const char *uplo,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3_int_t *info );

void   lapackf77_clavhe( const char *uplo, const char *trans, const char *diag,
                         magma_tally3_int_t *n, magma_tally3_int_t *nrhs,
                         magma_tally3FloatComplex *A, magma_tally3_int_t *lda,
                         magma_tally3_int_t *ipiv,
                         magma_tally3FloatComplex *B, magma_tally3_int_t *ldb,
                         magma_tally3_int_t *info );

void   lapackf77_cposv(  const char *uplo,
                         const magma_tally3_int_t *n, const magma_tally3_int_t *nrhs,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *B,  const magma_tally3_int_t *ldb,
                         magma_tally3_int_t *info );

void   lapackf77_cpotrf( const char *uplo,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3_int_t *info );

void   lapackf77_cpotri( const char *uplo,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3_int_t *info );

void   lapackf77_cpotrs( const char *uplo,
                         const magma_tally3_int_t *n, const magma_tally3_int_t *nrhs,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *B, const magma_tally3_int_t *ldb,
                         magma_tally3_int_t *info );

void   lapackf77_cstedc( const char *compz,
                         const magma_tally3_int_t *n,
                         float *d, float *e,
                         magma_tally3FloatComplex *Z, const magma_tally3_int_t *ldz,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         #ifdef COMPLEX
                         float *rwork, const magma_tally3_int_t *lrwork,
                         #endif
                         magma_tally3_int_t *iwork, const magma_tally3_int_t *liwork,
                         magma_tally3_int_t *info );

void   lapackf77_cstein( const magma_tally3_int_t *n,
                         const float *d, const float *e,
                         const magma_tally3_int_t *m,
                         const float *w,
                         const magma_tally3_int_t *iblock,
                         const magma_tally3_int_t *isplit,
                         magma_tally3FloatComplex *Z, const magma_tally3_int_t *ldz,
                         float *work, magma_tally3_int_t *iwork, magma_tally3_int_t *ifailv,
                         magma_tally3_int_t *info );

void   lapackf77_cstemr( const char *jobz, const char *range,
                         const magma_tally3_int_t *n,
                         float *d, float *e,
                         const float *vl, const float *vu,
                         const magma_tally3_int_t *il, const magma_tally3_int_t *iu,
                         magma_tally3_int_t *m,
                         float *w,
                         magma_tally3FloatComplex *Z, const magma_tally3_int_t *ldz,
                         const magma_tally3_int_t *nzc, magma_tally3_int_t *isuppz, magma_tally3_int_t *tryrac,
                         float *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *iwork, const magma_tally3_int_t *liwork,
                         magma_tally3_int_t *info );

void   lapackf77_csteqr( const char *compz,
                         const magma_tally3_int_t *n,
                         float *d, float *e,
                         magma_tally3FloatComplex *Z, const magma_tally3_int_t *ldz,
                         float *work,
                         magma_tally3_int_t *info );

#ifdef COMPLEX
void   lapackf77_csymv(  const char *uplo,
                         const magma_tally3_int_t *n,
                         const magma_tally3FloatComplex *alpha,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3FloatComplex *x, const magma_tally3_int_t *incx,
                         const magma_tally3FloatComplex *beta,
                               magma_tally3FloatComplex *y, const magma_tally3_int_t *incy );
#endif

void   lapackf77_ctrevc( const char *side, const char *howmny,
                         magma_tally3_int_t *select, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *T,  const magma_tally3_int_t *ldt,
                         magma_tally3FloatComplex *Vl, const magma_tally3_int_t *ldvl,
                         magma_tally3FloatComplex *Vr, const magma_tally3_int_t *ldvr,
                         const magma_tally3_int_t *mm, magma_tally3_int_t *m,
                         magma_tally3FloatComplex *work,
                         #ifdef COMPLEX
                         float *rwork,
                         #endif
                         magma_tally3_int_t *info );

void   lapackf77_ctrevc3( const char *side, const char *howmny,
                          magma_tally3_int_t *select, const magma_tally3_int_t *n,
                          magma_tally3FloatComplex *T,  const magma_tally3_int_t *ldt,
                          magma_tally3FloatComplex *VL, const magma_tally3_int_t *ldvl, 
                          magma_tally3FloatComplex *VR, const magma_tally3_int_t *ldvr,
                          const magma_tally3_int_t *mm,
                          const magma_tally3_int_t *mout,
                          magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                          #ifdef COMPLEX
                          float *rwork,
                          #endif
                          magma_tally3_int_t *info );

void   lapackf77_ctrtri( const char *uplo, const char *diag,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3_int_t *info );

void   lapackf77_cung2r( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work,
                         magma_tally3_int_t *info );

void   lapackf77_cungbr( const char *vect,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cunghr( const magma_tally3_int_t *n,
                         const magma_tally3_int_t *ilo, const magma_tally3_int_t *ihi,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cunglq( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cungql( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cungqr( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cungtr( const char *uplo,
                         const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cunm2r( const char *side, const char *trans,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc,
                         magma_tally3FloatComplex *work,
                         magma_tally3_int_t *info );

void   lapackf77_cunmbr( const char *vect, const char *side, const char *trans,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cunmlq( const char *side, const char *trans,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cunmql( const char *side, const char *trans,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cunmqr( const char *side, const char *trans,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

void   lapackf77_cunmtr( const char *side, const char *uplo, const char *trans,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         const magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         const magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         magma_tally3_int_t *info );

/*
 * Real precision extras
 */
void   lapackf77_sstebz( const char *range, const char *order,
                         const magma_tally3_int_t *n,
                         float *vl, float *vu,
                         magma_tally3_int_t *il, magma_tally3_int_t *iu,
                         float *abstol,
                         float *d, float *e,
                         const magma_tally3_int_t *m, const magma_tally3_int_t *nsplit,
                         float *w,
                         magma_tally3_int_t *iblock, magma_tally3_int_t *isplit,
                         float *work,
                         magma_tally3_int_t *iwork,
                         magma_tally3_int_t *info );

void   lapackf77_slaln2( const magma_tally3_int_t *ltrans,
                         const magma_tally3_int_t *na, const magma_tally3_int_t *nw,
                         const float *smin, const float *ca,
                         const float *a,  const magma_tally3_int_t *lda,
                         const float *d1, const float *d2,
                         const float *b,  const magma_tally3_int_t *ldb,
                         const float *wr, const float *wi,
                         float *x, const magma_tally3_int_t *ldx,
                         float *scale, float *xnorm, magma_tally3_int_t *info );

float lapackf77_slamc3( float *a, float *b );

void   lapackf77_slamrg( magma_tally3_int_t *n1, magma_tally3_int_t *n2,
                         float *a,
                         magma_tally3_int_t *dtrd1, magma_tally3_int_t *dtrd2, magma_tally3_int_t *index );

float lapackf77_slapy3( float *x, float *y, float *z );

void   lapackf77_slaed2( magma_tally3_int_t *k, magma_tally3_int_t *n, magma_tally3_int_t *cutpnt,
                         float *d, float *q, magma_tally3_int_t *ldq, magma_tally3_int_t *indxq,
                         float *rho, float *z,
                         float *dlmda, float *w, float *q2,
                         magma_tally3_int_t *indx, magma_tally3_int_t *indxc, magma_tally3_int_t *indxp,
                         magma_tally3_int_t *coltyp, magma_tally3_int_t *info);

void   lapackf77_slaed4( magma_tally3_int_t *n, magma_tally3_int_t *i,
                         float *d,
                         float *z,
                         float *delta,
                         float *rho,
                         float *dlam, magma_tally3_int_t *info );

void   lapackf77_slasrt( const char *id, const magma_tally3_int_t *n, float *d, magma_tally3_int_t *info );

/*
 * Testing functions
 */
#ifdef COMPLEX
void   lapackf77_cbdt01( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *kd,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *Q, const magma_tally3_int_t *ldq,
                         float *d, float *e,
                         magma_tally3FloatComplex *Pt, const magma_tally3_int_t *ldpt,
                         magma_tally3FloatComplex *work,
                         float *rwork,
                         float *resid );

void   lapackf77_cget22( const char *transa, const char *transe, const char *transw, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *E, const magma_tally3_int_t *lde,
                         magma_tally3FloatComplex *w,
                         magma_tally3FloatComplex *work,
                         float *rwork,
                         float *result );

void   lapackf77_chet21( const magma_tally3_int_t *itype, const char *uplo,
                         const magma_tally3_int_t *n, const magma_tally3_int_t *kband,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *d, float *e,
                         magma_tally3FloatComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3FloatComplex *V, const magma_tally3_int_t *ldv,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work,
                         float *rwork,
                         float *result );

void   lapackf77_chst01( const magma_tally3_int_t *n, const magma_tally3_int_t *ilo, const magma_tally3_int_t *ihi,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *H, const magma_tally3_int_t *ldh,
                         magma_tally3FloatComplex *Q, const magma_tally3_int_t *ldq,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         float *rwork,
                         float *result );

void   lapackf77_cstt21( const magma_tally3_int_t *n, const magma_tally3_int_t *kband,
                         float *AD,
                         float *AE,
                         float *SD,
                         float *SE,
                         magma_tally3FloatComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3FloatComplex *work,
                         float *rwork,
                         float *result );

void   lapackf77_cunt01( const char *rowcol, const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         float *rwork,
                         float *resid );
#else
void   lapackf77_cbdt01( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *kd,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *Q, const magma_tally3_int_t *ldq,
                         float *d, float *e,
                         magma_tally3FloatComplex *Pt, const magma_tally3_int_t *ldpt,
                         magma_tally3FloatComplex *work,
                         float *resid );

void   lapackf77_cget22( const char *transa, const char *transe, const char *transw, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *E, const magma_tally3_int_t *lde,
                         magma_tally3FloatComplex *wr,
                         magma_tally3FloatComplex *wi,
                         float *work,
                         float *result );

void   lapackf77_chet21( magma_tally3_int_t *itype, const char *uplo, const magma_tally3_int_t *n, const magma_tally3_int_t *kband,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *d, float *e,
                         magma_tally3FloatComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3FloatComplex *V, const magma_tally3_int_t *ldv,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work,
                         float *result );

void   lapackf77_chst01( const magma_tally3_int_t *n, const magma_tally3_int_t *ilo, const magma_tally3_int_t *ihi,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *H, const magma_tally3_int_t *ldh,
                         magma_tally3FloatComplex *Q, const magma_tally3_int_t *ldq,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         float *result );

void   lapackf77_cstt21( const magma_tally3_int_t *n, const magma_tally3_int_t *kband,
                         float *AD,
                         float *AE,
                         float *SD,
                         float *SE,
                         magma_tally3FloatComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3FloatComplex *work,
                         float *result );

void   lapackf77_cunt01( const char *rowcol, const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         float *resid );
#endif

void   lapackf77_clarfy( const char *uplo, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *V, const magma_tally3_int_t *incv,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc,
                         magma_tally3FloatComplex *work );

void   lapackf77_clarfx( const char *side, const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *V,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc,
                         magma_tally3FloatComplex *work );

float lapackf77_cqpt01( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         magma_tally3FloatComplex *A,
                         magma_tally3FloatComplex *Af, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *tau, magma_tally3_int_t *jpvt,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork );

void   lapackf77_cqrt02( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         magma_tally3FloatComplex *A,
                         magma_tally3FloatComplex *AF,
                         magma_tally3FloatComplex *Q,
                         magma_tally3FloatComplex *R, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         float *rwork,
                         float *result );

void   lapackf77_clatms( magma_tally3_int_t *m, magma_tally3_int_t *n,
                         const char *dist, magma_tally3_int_t *iseed, const char *sym, float *d,
                         magma_tally3_int_t *mode, const float *cond, const float *dmax,
                         magma_tally3_int_t *kl, magma_tally3_int_t *ku, const char *pack,
                         magma_tally3FloatComplex *a, magma_tally3_int_t *lda, magma_tally3FloatComplex *work, magma_tally3_int_t *info );

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif /* MAGMA_tally3_CLAPACK_H */
