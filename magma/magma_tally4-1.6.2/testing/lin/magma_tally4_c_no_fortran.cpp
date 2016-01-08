/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @generated from magma_tally4_z_no_fortran.cpp normal z -> c, Fri Jan 30 19:00:27 2015
       
       This is simply a copy of part of magma_tally4_clapack.h,
       with the { printf(...); } function body added to each function.
*/
#include <stdio.h>

#include "magma_tally4.h"
#include "magma_tally4_lapack.h"

#ifdef __cplusplus
extern "C" {
#endif

#define COMPLEX

static const char* format = "Cannot check results: %s unavailable, since there was no Fortran compiler.\n";

/*
 * Testing functions
 */
#ifdef COMPLEX
void   lapackf77_cbdt01( const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *kd,
                         magma_tally4FloatComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4FloatComplex *Q, const magma_tally4_int_t *ldq,
                         float *d, float *e,
                         magma_tally4FloatComplex *Pt, const magma_tally4_int_t *ldpt,
                         magma_tally4FloatComplex *work,
                         float *rwork,
                         float *resid )
                         { printf( format, __func__ ); }

void   lapackf77_cget22( const char *transa, const char *transe, const char *transw, const magma_tally4_int_t *n,
                         magma_tally4FloatComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4FloatComplex *E, const magma_tally4_int_t *lde,
                         magma_tally4FloatComplex *w,
                         magma_tally4FloatComplex *work,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_chet21( const magma_tally4_int_t *itype, const char *uplo,
                         const magma_tally4_int_t *n, const magma_tally4_int_t *kband,
                         magma_tally4FloatComplex *A, const magma_tally4_int_t *lda,
                         float *d, float *e,
                         magma_tally4FloatComplex *U, const magma_tally4_int_t *ldu,
                         magma_tally4FloatComplex *V, const magma_tally4_int_t *ldv,
                         magma_tally4FloatComplex *tau,
                         magma_tally4FloatComplex *work,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_chst01( const magma_tally4_int_t *n, const magma_tally4_int_t *ilo, const magma_tally4_int_t *ihi,
                         magma_tally4FloatComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4FloatComplex *H, const magma_tally4_int_t *ldh,
                         magma_tally4FloatComplex *Q, const magma_tally4_int_t *ldq,
                         magma_tally4FloatComplex *work, const magma_tally4_int_t *lwork,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_cstt21( const magma_tally4_int_t *n, const magma_tally4_int_t *kband,
                         float *AD,
                         float *AE,
                         float *SD,
                         float *SE,
                         magma_tally4FloatComplex *U, const magma_tally4_int_t *ldu,
                         magma_tally4FloatComplex *work,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_cunt01( const char *rowcol, const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4FloatComplex *U, const magma_tally4_int_t *ldu,
                         magma_tally4FloatComplex *work, const magma_tally4_int_t *lwork,
                         float *rwork,
                         float *resid )
                         { printf( format, __func__ ); }
#else
void   lapackf77_cbdt01( const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *kd,
                         magma_tally4FloatComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4FloatComplex *Q, const magma_tally4_int_t *ldq,
                         float *d, float *e,
                         magma_tally4FloatComplex *Pt, const magma_tally4_int_t *ldpt,
                         magma_tally4FloatComplex *work,
                         float *resid )
                         { printf( format, __func__ ); }

void   lapackf77_cget22( const char *transa, const char *transe, const char *transw, const magma_tally4_int_t *n,
                         magma_tally4FloatComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4FloatComplex *E, const magma_tally4_int_t *lde,
                         magma_tally4FloatComplex *wr,
                         magma_tally4FloatComplex *wi,
                         float *work,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_chet21( magma_tally4_int_t *itype, const char *uplo, const magma_tally4_int_t *n, const magma_tally4_int_t *kband,
                         magma_tally4FloatComplex *A, const magma_tally4_int_t *lda,
                         float *d, float *e,
                         magma_tally4FloatComplex *U, const magma_tally4_int_t *ldu,
                         magma_tally4FloatComplex *V, const magma_tally4_int_t *ldv,
                         magma_tally4FloatComplex *tau,
                         magma_tally4FloatComplex *work,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_chst01( const magma_tally4_int_t *n, const magma_tally4_int_t *ilo, const magma_tally4_int_t *ihi,
                         magma_tally4FloatComplex *A, const magma_tally4_int_t *lda,
                         magma_tally4FloatComplex *H, const magma_tally4_int_t *ldh,
                         magma_tally4FloatComplex *Q, const magma_tally4_int_t *ldq,
                         magma_tally4FloatComplex *work, const magma_tally4_int_t *lwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_cstt21( const magma_tally4_int_t *n, const magma_tally4_int_t *kband,
                         float *AD,
                         float *AE,
                         float *SD,
                         float *SE,
                         magma_tally4FloatComplex *U, const magma_tally4_int_t *ldu,
                         magma_tally4FloatComplex *work,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_cunt01( const char *rowcol, const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4FloatComplex *U, const magma_tally4_int_t *ldu,
                         magma_tally4FloatComplex *work, const magma_tally4_int_t *lwork,
                         float *resid )
                         { printf( format, __func__ ); }
#endif

void   lapackf77_clarfy( const char *uplo, const magma_tally4_int_t *n,
                         magma_tally4FloatComplex *V, const magma_tally4_int_t *incv,
                         magma_tally4FloatComplex *tau,
                         magma_tally4FloatComplex *C, const magma_tally4_int_t *ldc,
                         magma_tally4FloatComplex *work )
                         { printf( format, __func__ ); }

void   lapackf77_clarfx( const char *side, const magma_tally4_int_t *m, const magma_tally4_int_t *n,
                         magma_tally4FloatComplex *V,
                         magma_tally4FloatComplex *tau,
                         magma_tally4FloatComplex *C, const magma_tally4_int_t *ldc,
                         magma_tally4FloatComplex *work )
                         { printf( format, __func__ ); }

float lapackf77_cqpt01( const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                         magma_tally4FloatComplex *A,
                         magma_tally4FloatComplex *Af, const magma_tally4_int_t *lda,
                         magma_tally4FloatComplex *tau, magma_tally4_int_t *jpvt,
                         magma_tally4FloatComplex *work, const magma_tally4_int_t *lwork )
                         { printf( format, __func__ ); return -1; }

void   lapackf77_cqrt02( const magma_tally4_int_t *m, const magma_tally4_int_t *n, const magma_tally4_int_t *k,
                         magma_tally4FloatComplex *A,
                         magma_tally4FloatComplex *AF,
                         magma_tally4FloatComplex *Q,
                         magma_tally4FloatComplex *R, const magma_tally4_int_t *lda,
                         magma_tally4FloatComplex *tau,
                         magma_tally4FloatComplex *work, const magma_tally4_int_t *lwork,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_clatms( magma_tally4_int_t *m, magma_tally4_int_t *n,
                         const char *dist, magma_tally4_int_t *iseed, const char *sym, float *d,
                         magma_tally4_int_t *mode, const float *cond, const float *dmax,
                         magma_tally4_int_t *kl, magma_tally4_int_t *ku, const char *pack,
                         magma_tally4FloatComplex *a, magma_tally4_int_t *lda, magma_tally4FloatComplex *work, magma_tally4_int_t *info )
                         { printf( format, __func__ ); }

#ifdef __cplusplus
}
#endif
