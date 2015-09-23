/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @generated from magma_minproduct_z_no_fortran.cpp normal z -> c, Fri Jan 30 19:00:27 2015
       
       This is simply a copy of part of magma_minproduct_clapack.h,
       with the { printf(...); } function body added to each function.
*/
#include <stdio.h>

#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"

#ifdef __cplusplus
extern "C" {
#endif

#define COMPLEX

static const char* format = "Cannot check results: %s unavailable, since there was no Fortran compiler.\n";

/*
 * Testing functions
 */
#ifdef COMPLEX
void   lapackf77_cbdt01( const magma_minproduct_int_t *m, const magma_minproduct_int_t *n, const magma_minproduct_int_t *kd,
                         magma_minproductFloatComplex *A, const magma_minproduct_int_t *lda,
                         magma_minproductFloatComplex *Q, const magma_minproduct_int_t *ldq,
                         float *d, float *e,
                         magma_minproductFloatComplex *Pt, const magma_minproduct_int_t *ldpt,
                         magma_minproductFloatComplex *work,
                         float *rwork,
                         float *resid )
                         { printf( format, __func__ ); }

void   lapackf77_cget22( const char *transa, const char *transe, const char *transw, const magma_minproduct_int_t *n,
                         magma_minproductFloatComplex *A, const magma_minproduct_int_t *lda,
                         magma_minproductFloatComplex *E, const magma_minproduct_int_t *lde,
                         magma_minproductFloatComplex *w,
                         magma_minproductFloatComplex *work,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_chet21( const magma_minproduct_int_t *itype, const char *uplo,
                         const magma_minproduct_int_t *n, const magma_minproduct_int_t *kband,
                         magma_minproductFloatComplex *A, const magma_minproduct_int_t *lda,
                         float *d, float *e,
                         magma_minproductFloatComplex *U, const magma_minproduct_int_t *ldu,
                         magma_minproductFloatComplex *V, const magma_minproduct_int_t *ldv,
                         magma_minproductFloatComplex *tau,
                         magma_minproductFloatComplex *work,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_chst01( const magma_minproduct_int_t *n, const magma_minproduct_int_t *ilo, const magma_minproduct_int_t *ihi,
                         magma_minproductFloatComplex *A, const magma_minproduct_int_t *lda,
                         magma_minproductFloatComplex *H, const magma_minproduct_int_t *ldh,
                         magma_minproductFloatComplex *Q, const magma_minproduct_int_t *ldq,
                         magma_minproductFloatComplex *work, const magma_minproduct_int_t *lwork,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_cstt21( const magma_minproduct_int_t *n, const magma_minproduct_int_t *kband,
                         float *AD,
                         float *AE,
                         float *SD,
                         float *SE,
                         magma_minproductFloatComplex *U, const magma_minproduct_int_t *ldu,
                         magma_minproductFloatComplex *work,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_cunt01( const char *rowcol, const magma_minproduct_int_t *m, const magma_minproduct_int_t *n,
                         magma_minproductFloatComplex *U, const magma_minproduct_int_t *ldu,
                         magma_minproductFloatComplex *work, const magma_minproduct_int_t *lwork,
                         float *rwork,
                         float *resid )
                         { printf( format, __func__ ); }
#else
void   lapackf77_cbdt01( const magma_minproduct_int_t *m, const magma_minproduct_int_t *n, const magma_minproduct_int_t *kd,
                         magma_minproductFloatComplex *A, const magma_minproduct_int_t *lda,
                         magma_minproductFloatComplex *Q, const magma_minproduct_int_t *ldq,
                         float *d, float *e,
                         magma_minproductFloatComplex *Pt, const magma_minproduct_int_t *ldpt,
                         magma_minproductFloatComplex *work,
                         float *resid )
                         { printf( format, __func__ ); }

void   lapackf77_cget22( const char *transa, const char *transe, const char *transw, const magma_minproduct_int_t *n,
                         magma_minproductFloatComplex *A, const magma_minproduct_int_t *lda,
                         magma_minproductFloatComplex *E, const magma_minproduct_int_t *lde,
                         magma_minproductFloatComplex *wr,
                         magma_minproductFloatComplex *wi,
                         float *work,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_chet21( magma_minproduct_int_t *itype, const char *uplo, const magma_minproduct_int_t *n, const magma_minproduct_int_t *kband,
                         magma_minproductFloatComplex *A, const magma_minproduct_int_t *lda,
                         float *d, float *e,
                         magma_minproductFloatComplex *U, const magma_minproduct_int_t *ldu,
                         magma_minproductFloatComplex *V, const magma_minproduct_int_t *ldv,
                         magma_minproductFloatComplex *tau,
                         magma_minproductFloatComplex *work,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_chst01( const magma_minproduct_int_t *n, const magma_minproduct_int_t *ilo, const magma_minproduct_int_t *ihi,
                         magma_minproductFloatComplex *A, const magma_minproduct_int_t *lda,
                         magma_minproductFloatComplex *H, const magma_minproduct_int_t *ldh,
                         magma_minproductFloatComplex *Q, const magma_minproduct_int_t *ldq,
                         magma_minproductFloatComplex *work, const magma_minproduct_int_t *lwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_cstt21( const magma_minproduct_int_t *n, const magma_minproduct_int_t *kband,
                         float *AD,
                         float *AE,
                         float *SD,
                         float *SE,
                         magma_minproductFloatComplex *U, const magma_minproduct_int_t *ldu,
                         magma_minproductFloatComplex *work,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_cunt01( const char *rowcol, const magma_minproduct_int_t *m, const magma_minproduct_int_t *n,
                         magma_minproductFloatComplex *U, const magma_minproduct_int_t *ldu,
                         magma_minproductFloatComplex *work, const magma_minproduct_int_t *lwork,
                         float *resid )
                         { printf( format, __func__ ); }
#endif

void   lapackf77_clarfy( const char *uplo, const magma_minproduct_int_t *n,
                         magma_minproductFloatComplex *V, const magma_minproduct_int_t *incv,
                         magma_minproductFloatComplex *tau,
                         magma_minproductFloatComplex *C, const magma_minproduct_int_t *ldc,
                         magma_minproductFloatComplex *work )
                         { printf( format, __func__ ); }

void   lapackf77_clarfx( const char *side, const magma_minproduct_int_t *m, const magma_minproduct_int_t *n,
                         magma_minproductFloatComplex *V,
                         magma_minproductFloatComplex *tau,
                         magma_minproductFloatComplex *C, const magma_minproduct_int_t *ldc,
                         magma_minproductFloatComplex *work )
                         { printf( format, __func__ ); }

float lapackf77_cqpt01( const magma_minproduct_int_t *m, const magma_minproduct_int_t *n, const magma_minproduct_int_t *k,
                         magma_minproductFloatComplex *A,
                         magma_minproductFloatComplex *Af, const magma_minproduct_int_t *lda,
                         magma_minproductFloatComplex *tau, magma_minproduct_int_t *jpvt,
                         magma_minproductFloatComplex *work, const magma_minproduct_int_t *lwork )
                         { printf( format, __func__ ); return -1; }

void   lapackf77_cqrt02( const magma_minproduct_int_t *m, const magma_minproduct_int_t *n, const magma_minproduct_int_t *k,
                         magma_minproductFloatComplex *A,
                         magma_minproductFloatComplex *AF,
                         magma_minproductFloatComplex *Q,
                         magma_minproductFloatComplex *R, const magma_minproduct_int_t *lda,
                         magma_minproductFloatComplex *tau,
                         magma_minproductFloatComplex *work, const magma_minproduct_int_t *lwork,
                         float *rwork,
                         float *result )
                         { printf( format, __func__ ); }

void   lapackf77_clatms( magma_minproduct_int_t *m, magma_minproduct_int_t *n,
                         const char *dist, magma_minproduct_int_t *iseed, const char *sym, float *d,
                         magma_minproduct_int_t *mode, const float *cond, const float *dmax,
                         magma_minproduct_int_t *kl, magma_minproduct_int_t *ku, const char *pack,
                         magma_minproductFloatComplex *a, magma_minproduct_int_t *lda, magma_minproductFloatComplex *work, magma_minproduct_int_t *info )
                         { printf( format, __func__ ); }

#ifdef __cplusplus
}
#endif
