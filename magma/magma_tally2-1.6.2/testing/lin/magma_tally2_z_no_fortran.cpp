/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @precisions normal z -> s d c
       
       This is simply a copy of part of magma_tally2_zlapack.h,
       with the { printf(...); } function body added to each function.
*/
#include <stdio.h>

#include "magma_tally2.h"
#include "magma_tally2_lapack.h"

#ifdef __cplusplus
extern "C" {
#endif

#define COMPLEX

static const char* format = "Cannot check results: %s unavailable, since there was no Fortran compiler.\n";

/*
 * Testing functions
 */
#ifdef COMPLEX
void   lapackf77_zbdt01( const magma_tally2_int_t *m, const magma_tally2_int_t *n, const magma_tally2_int_t *kd,
                         magma_tally2DoubleComplex *A, const magma_tally2_int_t *lda,
                         magma_tally2DoubleComplex *Q, const magma_tally2_int_t *ldq,
                         double *d, double *e,
                         magma_tally2DoubleComplex *Pt, const magma_tally2_int_t *ldpt,
                         magma_tally2DoubleComplex *work,
                         double *rwork,
                         double *resid )
                         { printf( format, __func__ ); }

void   lapackf77_zget22( const char *transa, const char *transe, const char *transw, const magma_tally2_int_t *n,
                         magma_tally2DoubleComplex *A, const magma_tally2_int_t *lda,
                         magma_tally2DoubleComplex *E, const magma_tally2_int_t *lde,
                         magma_tally2DoubleComplex *w,
                         magma_tally2DoubleComplex *work,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zhet21( const magma_tally2_int_t *itype, const char *uplo,
                         const magma_tally2_int_t *n, const magma_tally2_int_t *kband,
                         magma_tally2DoubleComplex *A, const magma_tally2_int_t *lda,
                         double *d, double *e,
                         magma_tally2DoubleComplex *U, const magma_tally2_int_t *ldu,
                         magma_tally2DoubleComplex *V, const magma_tally2_int_t *ldv,
                         magma_tally2DoubleComplex *tau,
                         magma_tally2DoubleComplex *work,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zhst01( const magma_tally2_int_t *n, const magma_tally2_int_t *ilo, const magma_tally2_int_t *ihi,
                         magma_tally2DoubleComplex *A, const magma_tally2_int_t *lda,
                         magma_tally2DoubleComplex *H, const magma_tally2_int_t *ldh,
                         magma_tally2DoubleComplex *Q, const magma_tally2_int_t *ldq,
                         magma_tally2DoubleComplex *work, const magma_tally2_int_t *lwork,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zstt21( const magma_tally2_int_t *n, const magma_tally2_int_t *kband,
                         double *AD,
                         double *AE,
                         double *SD,
                         double *SE,
                         magma_tally2DoubleComplex *U, const magma_tally2_int_t *ldu,
                         magma_tally2DoubleComplex *work,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zunt01( const char *rowcol, const magma_tally2_int_t *m, const magma_tally2_int_t *n,
                         magma_tally2DoubleComplex *U, const magma_tally2_int_t *ldu,
                         magma_tally2DoubleComplex *work, const magma_tally2_int_t *lwork,
                         double *rwork,
                         double *resid )
                         { printf( format, __func__ ); }
#else
void   lapackf77_zbdt01( const magma_tally2_int_t *m, const magma_tally2_int_t *n, const magma_tally2_int_t *kd,
                         magma_tally2DoubleComplex *A, const magma_tally2_int_t *lda,
                         magma_tally2DoubleComplex *Q, const magma_tally2_int_t *ldq,
                         double *d, double *e,
                         magma_tally2DoubleComplex *Pt, const magma_tally2_int_t *ldpt,
                         magma_tally2DoubleComplex *work,
                         double *resid )
                         { printf( format, __func__ ); }

void   lapackf77_zget22( const char *transa, const char *transe, const char *transw, const magma_tally2_int_t *n,
                         magma_tally2DoubleComplex *A, const magma_tally2_int_t *lda,
                         magma_tally2DoubleComplex *E, const magma_tally2_int_t *lde,
                         magma_tally2DoubleComplex *wr,
                         magma_tally2DoubleComplex *wi,
                         double *work,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zhet21( magma_tally2_int_t *itype, const char *uplo, const magma_tally2_int_t *n, const magma_tally2_int_t *kband,
                         magma_tally2DoubleComplex *A, const magma_tally2_int_t *lda,
                         double *d, double *e,
                         magma_tally2DoubleComplex *U, const magma_tally2_int_t *ldu,
                         magma_tally2DoubleComplex *V, const magma_tally2_int_t *ldv,
                         magma_tally2DoubleComplex *tau,
                         magma_tally2DoubleComplex *work,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zhst01( const magma_tally2_int_t *n, const magma_tally2_int_t *ilo, const magma_tally2_int_t *ihi,
                         magma_tally2DoubleComplex *A, const magma_tally2_int_t *lda,
                         magma_tally2DoubleComplex *H, const magma_tally2_int_t *ldh,
                         magma_tally2DoubleComplex *Q, const magma_tally2_int_t *ldq,
                         magma_tally2DoubleComplex *work, const magma_tally2_int_t *lwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zstt21( const magma_tally2_int_t *n, const magma_tally2_int_t *kband,
                         double *AD,
                         double *AE,
                         double *SD,
                         double *SE,
                         magma_tally2DoubleComplex *U, const magma_tally2_int_t *ldu,
                         magma_tally2DoubleComplex *work,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zunt01( const char *rowcol, const magma_tally2_int_t *m, const magma_tally2_int_t *n,
                         magma_tally2DoubleComplex *U, const magma_tally2_int_t *ldu,
                         magma_tally2DoubleComplex *work, const magma_tally2_int_t *lwork,
                         double *resid )
                         { printf( format, __func__ ); }
#endif

void   lapackf77_zlarfy( const char *uplo, const magma_tally2_int_t *n,
                         magma_tally2DoubleComplex *V, const magma_tally2_int_t *incv,
                         magma_tally2DoubleComplex *tau,
                         magma_tally2DoubleComplex *C, const magma_tally2_int_t *ldc,
                         magma_tally2DoubleComplex *work )
                         { printf( format, __func__ ); }

void   lapackf77_zlarfx( const char *side, const magma_tally2_int_t *m, const magma_tally2_int_t *n,
                         magma_tally2DoubleComplex *V,
                         magma_tally2DoubleComplex *tau,
                         magma_tally2DoubleComplex *C, const magma_tally2_int_t *ldc,
                         magma_tally2DoubleComplex *work )
                         { printf( format, __func__ ); }

double lapackf77_zqpt01( const magma_tally2_int_t *m, const magma_tally2_int_t *n, const magma_tally2_int_t *k,
                         magma_tally2DoubleComplex *A,
                         magma_tally2DoubleComplex *Af, const magma_tally2_int_t *lda,
                         magma_tally2DoubleComplex *tau, magma_tally2_int_t *jpvt,
                         magma_tally2DoubleComplex *work, const magma_tally2_int_t *lwork )
                         { printf( format, __func__ ); return -1; }

void   lapackf77_zqrt02( const magma_tally2_int_t *m, const magma_tally2_int_t *n, const magma_tally2_int_t *k,
                         magma_tally2DoubleComplex *A,
                         magma_tally2DoubleComplex *AF,
                         magma_tally2DoubleComplex *Q,
                         magma_tally2DoubleComplex *R, const magma_tally2_int_t *lda,
                         magma_tally2DoubleComplex *tau,
                         magma_tally2DoubleComplex *work, const magma_tally2_int_t *lwork,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zlatms( magma_tally2_int_t *m, magma_tally2_int_t *n,
                         const char *dist, magma_tally2_int_t *iseed, const char *sym, double *d,
                         magma_tally2_int_t *mode, const double *cond, const double *dmax,
                         magma_tally2_int_t *kl, magma_tally2_int_t *ku, const char *pack,
                         magma_tally2DoubleComplex *a, magma_tally2_int_t *lda, magma_tally2DoubleComplex *work, magma_tally2_int_t *info )
                         { printf( format, __func__ ); }

#ifdef __cplusplus
}
#endif
