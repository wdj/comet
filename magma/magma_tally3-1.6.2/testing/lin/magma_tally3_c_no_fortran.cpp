/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @generated from magma_tally3_z_no_fortran.cpp normal z -> c, Fri Jan 30 19:00:27 2015
       
       This is simply a copy of part of magma_tally3_clapack.h,
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
void   lapackf77_cbdt01( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *kd,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *Q, const magma_tally3_int_t *ldq,
                         float *d, float *e,
                         magma_tally3FloatComplex *Pt, const magma_tally3_int_t *ldpt,
                         magma_tally3FloatComplex *work,
                         float *rwork,
                         float *resid )
                         { printf( format, __func__ ); }

void   lapackf77_cget22( const char *transa, const char *transe, const char *transw, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *E, const magma_tally3_int_t *lde,
                         magma_tally3FloatComplex *w,
                         magma_tally3FloatComplex *work,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_chet21( const magma_tally3_int_t *itype, const char *uplo,
                         const magma_tally3_int_t *n, const magma_tally3_int_t *kband,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *d, float *e,
                         magma_tally3FloatComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3FloatComplex *V, const magma_tally3_int_t *ldv,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_chst01( const magma_tally3_int_t *n, const magma_tally3_int_t *ilo, const magma_tally3_int_t *ihi,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *H, const magma_tally3_int_t *ldh,
                         magma_tally3FloatComplex *Q, const magma_tally3_int_t *ldq,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_cstt21( const magma_tally3_int_t *n, const magma_tally3_int_t *kband,
                         float *AD,
                         float *AE,
                         float *SD,
                         float *SE,
                         magma_tally3FloatComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3FloatComplex *work,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_cunt01( const char *rowcol, const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         float *rwork,
                         float *resid )
                         { printf( format, __func__ ); }
#else
void   lapackf77_cbdt01( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *kd,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *Q, const magma_tally3_int_t *ldq,
                         float *d, float *e,
                         magma_tally3FloatComplex *Pt, const magma_tally3_int_t *ldpt,
                         magma_tally3FloatComplex *work,
                         float *resid )
                         { printf( format, __func__ ); }

void   lapackf77_cget22( const char *transa, const char *transe, const char *transw, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *E, const magma_tally3_int_t *lde,
                         magma_tally3FloatComplex *wr,
                         magma_tally3FloatComplex *wi,
                         float *work,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_chet21( magma_tally3_int_t *itype, const char *uplo, const magma_tally3_int_t *n, const magma_tally3_int_t *kband,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         float *d, float *e,
                         magma_tally3FloatComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3FloatComplex *V, const magma_tally3_int_t *ldv,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_chst01( const magma_tally3_int_t *n, const magma_tally3_int_t *ilo, const magma_tally3_int_t *ihi,
                         magma_tally3FloatComplex *A, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *H, const magma_tally3_int_t *ldh,
                         magma_tally3FloatComplex *Q, const magma_tally3_int_t *ldq,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_cstt21( const magma_tally3_int_t *n, const magma_tally3_int_t *kband,
                         float *AD,
                         float *AE,
                         float *SD,
                         float *SE,
                         magma_tally3FloatComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3FloatComplex *work,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_cunt01( const char *rowcol, const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *U, const magma_tally3_int_t *ldu,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         float *resid )
                         { printf( format, __func__ ); }
#endif

void   lapackf77_clarfy( const char *uplo, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *V, const magma_tally3_int_t *incv,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc,
                         magma_tally3FloatComplex *work )
                         { printf( format, __func__ ); }

void   lapackf77_clarfx( const char *side, const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         magma_tally3FloatComplex *V,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *C, const magma_tally3_int_t *ldc,
                         magma_tally3FloatComplex *work )
                         { printf( format, __func__ ); }

float lapackf77_cqpt01( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         magma_tally3FloatComplex *A,
                         magma_tally3FloatComplex *Af, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *tau, magma_tally3_int_t *jpvt,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork )
                         { printf( format, __func__ ); return -1; }

void   lapackf77_cqrt02( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         magma_tally3FloatComplex *A,
                         magma_tally3FloatComplex *AF,
                         magma_tally3FloatComplex *Q,
                         magma_tally3FloatComplex *R, const magma_tally3_int_t *lda,
                         magma_tally3FloatComplex *tau,
                         magma_tally3FloatComplex *work, const magma_tally3_int_t *lwork,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_clatms( magma_tally3_int_t *m, magma_tally3_int_t *n,
                         const char *dist, magma_tally3_int_t *iseed, const char *sym, float *d,
                         magma_tally3_int_t *mode, const float *cond, const float *dmax,
                         magma_tally3_int_t *kl, magma_tally3_int_t *ku, const char *pack,
                         magma_tally3FloatComplex *a, magma_tally3_int_t *lda, magma_tally3FloatComplex *work, magma_tally3_int_t *info )
                         { printf( format, __func__ ); }

#ifdef __cplusplus
}
#endif
