/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @generated from magma_tally2_z_no_fortran.cpp normal z -> c, Fri Jan 30 19:00:27 2015
       
       This is simply a copy of part of magma_tally2_clapack.h,
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
void   lapackf77_cbdt01( const magma_tally2_int_t *m, const magma_tally2_int_t *n, const magma_tally2_int_t *kd,
                         magma_tally2FloatComplex *A, const magma_tally2_int_t *lda,
                         magma_tally2FloatComplex *Q, const magma_tally2_int_t *ldq,
                         float *d, float *e,
                         magma_tally2FloatComplex *Pt, const magma_tally2_int_t *ldpt,
                         magma_tally2FloatComplex *work,
                         float *rwork,
                         float *resid )
                         { printf( format, __func__ ); }

void   lapackf77_cget22( const char *transa, const char *transe, const char *transw, const magma_tally2_int_t *n,
                         magma_tally2FloatComplex *A, const magma_tally2_int_t *lda,
                         magma_tally2FloatComplex *E, const magma_tally2_int_t *lde,
                         magma_tally2FloatComplex *w,
                         magma_tally2FloatComplex *work,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_chet21( const magma_tally2_int_t *itype, const char *uplo,
                         const magma_tally2_int_t *n, const magma_tally2_int_t *kband,
                         magma_tally2FloatComplex *A, const magma_tally2_int_t *lda,
                         float *d, float *e,
                         magma_tally2FloatComplex *U, const magma_tally2_int_t *ldu,
                         magma_tally2FloatComplex *V, const magma_tally2_int_t *ldv,
                         magma_tally2FloatComplex *tau,
                         magma_tally2FloatComplex *work,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_chst01( const magma_tally2_int_t *n, const magma_tally2_int_t *ilo, const magma_tally2_int_t *ihi,
                         magma_tally2FloatComplex *A, const magma_tally2_int_t *lda,
                         magma_tally2FloatComplex *H, const magma_tally2_int_t *ldh,
                         magma_tally2FloatComplex *Q, const magma_tally2_int_t *ldq,
                         magma_tally2FloatComplex *work, const magma_tally2_int_t *lwork,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_cstt21( const magma_tally2_int_t *n, const magma_tally2_int_t *kband,
                         float *AD,
                         float *AE,
                         float *SD,
                         float *SE,
                         magma_tally2FloatComplex *U, const magma_tally2_int_t *ldu,
                         magma_tally2FloatComplex *work,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_cunt01( const char *rowcol, const magma_tally2_int_t *m, const magma_tally2_int_t *n,
                         magma_tally2FloatComplex *U, const magma_tally2_int_t *ldu,
                         magma_tally2FloatComplex *work, const magma_tally2_int_t *lwork,
                         float *rwork,
                         float *resid )
                         { printf( format, __func__ ); }
#else
void   lapackf77_cbdt01( const magma_tally2_int_t *m, const magma_tally2_int_t *n, const magma_tally2_int_t *kd,
                         magma_tally2FloatComplex *A, const magma_tally2_int_t *lda,
                         magma_tally2FloatComplex *Q, const magma_tally2_int_t *ldq,
                         float *d, float *e,
                         magma_tally2FloatComplex *Pt, const magma_tally2_int_t *ldpt,
                         magma_tally2FloatComplex *work,
                         float *resid )
                         { printf( format, __func__ ); }

void   lapackf77_cget22( const char *transa, const char *transe, const char *transw, const magma_tally2_int_t *n,
                         magma_tally2FloatComplex *A, const magma_tally2_int_t *lda,
                         magma_tally2FloatComplex *E, const magma_tally2_int_t *lde,
                         magma_tally2FloatComplex *wr,
                         magma_tally2FloatComplex *wi,
                         float *work,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_chet21( magma_tally2_int_t *itype, const char *uplo, const magma_tally2_int_t *n, const magma_tally2_int_t *kband,
                         magma_tally2FloatComplex *A, const magma_tally2_int_t *lda,
                         float *d, float *e,
                         magma_tally2FloatComplex *U, const magma_tally2_int_t *ldu,
                         magma_tally2FloatComplex *V, const magma_tally2_int_t *ldv,
                         magma_tally2FloatComplex *tau,
                         magma_tally2FloatComplex *work,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_chst01( const magma_tally2_int_t *n, const magma_tally2_int_t *ilo, const magma_tally2_int_t *ihi,
                         magma_tally2FloatComplex *A, const magma_tally2_int_t *lda,
                         magma_tally2FloatComplex *H, const magma_tally2_int_t *ldh,
                         magma_tally2FloatComplex *Q, const magma_tally2_int_t *ldq,
                         magma_tally2FloatComplex *work, const magma_tally2_int_t *lwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_cstt21( const magma_tally2_int_t *n, const magma_tally2_int_t *kband,
                         float *AD,
                         float *AE,
                         float *SD,
                         float *SE,
                         magma_tally2FloatComplex *U, const magma_tally2_int_t *ldu,
                         magma_tally2FloatComplex *work,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_cunt01( const char *rowcol, const magma_tally2_int_t *m, const magma_tally2_int_t *n,
                         magma_tally2FloatComplex *U, const magma_tally2_int_t *ldu,
                         magma_tally2FloatComplex *work, const magma_tally2_int_t *lwork,
                         float *resid )
                         { printf( format, __func__ ); }
#endif

void   lapackf77_clarfy( const char *uplo, const magma_tally2_int_t *n,
                         magma_tally2FloatComplex *V, const magma_tally2_int_t *incv,
                         magma_tally2FloatComplex *tau,
                         magma_tally2FloatComplex *C, const magma_tally2_int_t *ldc,
                         magma_tally2FloatComplex *work )
                         { printf( format, __func__ ); }

void   lapackf77_clarfx( const char *side, const magma_tally2_int_t *m, const magma_tally2_int_t *n,
                         magma_tally2FloatComplex *V,
                         magma_tally2FloatComplex *tau,
                         magma_tally2FloatComplex *C, const magma_tally2_int_t *ldc,
                         magma_tally2FloatComplex *work )
                         { printf( format, __func__ ); }

float lapackf77_cqpt01( const magma_tally2_int_t *m, const magma_tally2_int_t *n, const magma_tally2_int_t *k,
                         magma_tally2FloatComplex *A,
                         magma_tally2FloatComplex *Af, const magma_tally2_int_t *lda,
                         magma_tally2FloatComplex *tau, magma_tally2_int_t *jpvt,
                         magma_tally2FloatComplex *work, const magma_tally2_int_t *lwork )
                         { printf( format, __func__ ); return -1; }

void   lapackf77_cqrt02( const magma_tally2_int_t *m, const magma_tally2_int_t *n, const magma_tally2_int_t *k,
                         magma_tally2FloatComplex *A,
                         magma_tally2FloatComplex *AF,
                         magma_tally2FloatComplex *Q,
                         magma_tally2FloatComplex *R, const magma_tally2_int_t *lda,
                         magma_tally2FloatComplex *tau,
                         magma_tally2FloatComplex *work, const magma_tally2_int_t *lwork,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_clatms( magma_tally2_int_t *m, magma_tally2_int_t *n,
                         const char *dist, magma_tally2_int_t *iseed, const char *sym, float *d,
                         magma_tally2_int_t *mode, const float *cond, const float *dmax,
                         magma_tally2_int_t *kl, magma_tally2_int_t *ku, const char *pack,
                         magma_tally2FloatComplex *a, magma_tally2_int_t *lda, magma_tally2FloatComplex *work, magma_tally2_int_t *info )
                         { printf( format, __func__ ); }

#ifdef __cplusplus
}
#endif
