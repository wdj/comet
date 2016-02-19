/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @precisions normal z -> s d c
       
       This is simply a copy of part of magma_tally3_zlapack.h,
       with the { printf(...); } function body added to each function.
*/
#include <stdio.h>

#include "magma_tally3.h"
#include "magma_tally3_lapack.h"

#ifdef __cplusplus
extern "C" {
#endif

#define COMPLEX

static const char* format = "Cannot check results: %s unavailable, since there was no Fortran compiler.\n";

/*
 * Testing functions
 */
#ifdef COMPLEX
void   lapackf77_zbdt01( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *kd,
                         magma_tally3DoubleComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3DoubleComplex *Q, const magma_tally3_int_t *ldq,
                         double *d, double *e,
                         magma_tally3DoubleComplex *Pt, const magma_tally3_int_t *ldpt,
                         magma_tally3DoubleComplex *work,
                         double *rwork,
                         double *resid )
                         { printf( format, __func__ ); }

void   lapackf77_zget22( const char *transa, const char *transe, const char *transw, const magma_tally3_int_t *n,
                         magma_tally3DoubleComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3DoubleComplex *E, const magma_tally3_int_t *lde,
                         magma_tally3DoubleComplex *w,
                         magma_tally3DoubleComplex *work,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zhet21( const magma_tally3_int_t *itype, const char *uplo,
                         const magma_tally3_int_t *n, const magma_tally3_int_t *kband,
                         magma_tally3DoubleComplex *A, const magma_tally3_int_t *lda,
                         double *d, double *e,
                         magma_tally3DoubleComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3DoubleComplex *V, const magma_tally3_int_t *ldv,
                         magma_tally3DoubleComplex *tau,
                         magma_tally3DoubleComplex *work,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zhst01( const magma_tally3_int_t *n, const magma_tally3_int_t *ilo, const magma_tally3_int_t *ihi,
                         magma_tally3DoubleComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3DoubleComplex *H, const magma_tally3_int_t *ldh,
                         magma_tally3DoubleComplex *Q, const magma_tally3_int_t *ldq,
                         magma_tally3DoubleComplex *work, const magma_tally3_int_t *lwork,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zstt21( const magma_tally3_int_t *n, const magma_tally3_int_t *kband,
                         double *AD,
                         double *AE,
                         double *SD,
                         double *SE,
                         magma_tally3DoubleComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3DoubleComplex *work,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zunt01( const char *rowcol, const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3DoubleComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3DoubleComplex *work, const magma_tally3_int_t *lwork,
                         double *rwork,
                         double *resid )
                         { printf( format, __func__ ); }
#else
void   lapackf77_zbdt01( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *kd,
                         magma_tally3DoubleComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3DoubleComplex *Q, const magma_tally3_int_t *ldq,
                         double *d, double *e,
                         magma_tally3DoubleComplex *Pt, const magma_tally3_int_t *ldpt,
                         magma_tally3DoubleComplex *work,
                         double *resid )
                         { printf( format, __func__ ); }

void   lapackf77_zget22( const char *transa, const char *transe, const char *transw, const magma_tally3_int_t *n,
                         magma_tally3DoubleComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3DoubleComplex *E, const magma_tally3_int_t *lde,
                         magma_tally3DoubleComplex *wr,
                         magma_tally3DoubleComplex *wi,
                         double *work,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zhet21( magma_tally3_int_t *itype, const char *uplo, const magma_tally3_int_t *n, const magma_tally3_int_t *kband,
                         magma_tally3DoubleComplex *A, const magma_tally3_int_t *lda,
                         double *d, double *e,
                         magma_tally3DoubleComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3DoubleComplex *V, const magma_tally3_int_t *ldv,
                         magma_tally3DoubleComplex *tau,
                         magma_tally3DoubleComplex *work,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zhst01( const magma_tally3_int_t *n, const magma_tally3_int_t *ilo, const magma_tally3_int_t *ihi,
                         magma_tally3DoubleComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3DoubleComplex *H, const magma_tally3_int_t *ldh,
                         magma_tally3DoubleComplex *Q, const magma_tally3_int_t *ldq,
                         magma_tally3DoubleComplex *work, const magma_tally3_int_t *lwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zstt21( const magma_tally3_int_t *n, const magma_tally3_int_t *kband,
                         double *AD,
                         double *AE,
                         double *SD,
                         double *SE,
                         magma_tally3DoubleComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3DoubleComplex *work,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zunt01( const char *rowcol, const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3DoubleComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3DoubleComplex *work, const magma_tally3_int_t *lwork,
                         double *resid )
                         { printf( format, __func__ ); }
#endif

void   lapackf77_zlarfy( const char *uplo, const magma_tally3_int_t *n,
                         magma_tally3DoubleComplex *V, const magma_tally3_int_t *incv,
                         magma_tally3DoubleComplex *tau,
                         magma_tally3DoubleComplex *C, const magma_tally3_int_t *ldc,
                         magma_tally3DoubleComplex *work )
                         { printf( format, __func__ ); }

void   lapackf77_zlarfx( const char *side, const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3DoubleComplex *V,
                         magma_tally3DoubleComplex *tau,
                         magma_tally3DoubleComplex *C, const magma_tally3_int_t *ldc,
                         magma_tally3DoubleComplex *work )
                         { printf( format, __func__ ); }

double lapackf77_zqpt01( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         magma_tally3DoubleComplex *A,
                         magma_tally3DoubleComplex *Af, const magma_tally3_int_t *lda,
                         magma_tally3DoubleComplex *tau, magma_tally3_int_t *jpvt,
                         magma_tally3DoubleComplex *work, const magma_tally3_int_t *lwork )
                         { printf( format, __func__ ); return -1; }

void   lapackf77_zqrt02( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         magma_tally3DoubleComplex *A,
                         magma_tally3DoubleComplex *AF,
                         magma_tally3DoubleComplex *Q,
                         magma_tally3DoubleComplex *R, const magma_tally3_int_t *lda,
                         magma_tally3DoubleComplex *tau,
                         magma_tally3DoubleComplex *work, const magma_tally3_int_t *lwork,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zlatms( magma_tally3_int_t *m, magma_tally3_int_t *n,
                         const char *dist, magma_tally3_int_t *iseed, const char *sym, double *d,
                         magma_tally3_int_t *mode, const double *cond, const double *dmax,
                         magma_tally3_int_t *kl, magma_tally3_int_t *ku, const char *pack,
                         magma_tally3DoubleComplex *a, magma_tally3_int_t *lda, magma_tally3DoubleComplex *work, magma_tally3_int_t *info )
                         { printf( format, __func__ ); }

#ifdef __cplusplus
}
#endif
