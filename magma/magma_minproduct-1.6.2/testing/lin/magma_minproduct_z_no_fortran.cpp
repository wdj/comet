/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @precisions normal z -> s d c
       
       This is simply a copy of part of magma_minproduct_zlapack.h,
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
void   lapackf77_zbdt01( const magma_minproduct_int_t *m, const magma_minproduct_int_t *n, const magma_minproduct_int_t *kd,
                         magma_minproductDoubleComplex *A, const magma_minproduct_int_t *lda,
                         magma_minproductDoubleComplex *Q, const magma_minproduct_int_t *ldq,
                         double *d, double *e,
                         magma_minproductDoubleComplex *Pt, const magma_minproduct_int_t *ldpt,
                         magma_minproductDoubleComplex *work,
                         double *rwork,
                         double *resid )
                         { printf( format, __func__ ); }

void   lapackf77_zget22( const char *transa, const char *transe, const char *transw, const magma_minproduct_int_t *n,
                         magma_minproductDoubleComplex *A, const magma_minproduct_int_t *lda,
                         magma_minproductDoubleComplex *E, const magma_minproduct_int_t *lde,
                         magma_minproductDoubleComplex *w,
                         magma_minproductDoubleComplex *work,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zhet21( const magma_minproduct_int_t *itype, const char *uplo,
                         const magma_minproduct_int_t *n, const magma_minproduct_int_t *kband,
                         magma_minproductDoubleComplex *A, const magma_minproduct_int_t *lda,
                         double *d, double *e,
                         magma_minproductDoubleComplex *U, const magma_minproduct_int_t *ldu,
                         magma_minproductDoubleComplex *V, const magma_minproduct_int_t *ldv,
                         magma_minproductDoubleComplex *tau,
                         magma_minproductDoubleComplex *work,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zhst01( const magma_minproduct_int_t *n, const magma_minproduct_int_t *ilo, const magma_minproduct_int_t *ihi,
                         magma_minproductDoubleComplex *A, const magma_minproduct_int_t *lda,
                         magma_minproductDoubleComplex *H, const magma_minproduct_int_t *ldh,
                         magma_minproductDoubleComplex *Q, const magma_minproduct_int_t *ldq,
                         magma_minproductDoubleComplex *work, const magma_minproduct_int_t *lwork,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zstt21( const magma_minproduct_int_t *n, const magma_minproduct_int_t *kband,
                         double *AD,
                         double *AE,
                         double *SD,
                         double *SE,
                         magma_minproductDoubleComplex *U, const magma_minproduct_int_t *ldu,
                         magma_minproductDoubleComplex *work,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zunt01( const char *rowcol, const magma_minproduct_int_t *m, const magma_minproduct_int_t *n,
                         magma_minproductDoubleComplex *U, const magma_minproduct_int_t *ldu,
                         magma_minproductDoubleComplex *work, const magma_minproduct_int_t *lwork,
                         double *rwork,
                         double *resid )
                         { printf( format, __func__ ); }
#else
void   lapackf77_zbdt01( const magma_minproduct_int_t *m, const magma_minproduct_int_t *n, const magma_minproduct_int_t *kd,
                         magma_minproductDoubleComplex *A, const magma_minproduct_int_t *lda,
                         magma_minproductDoubleComplex *Q, const magma_minproduct_int_t *ldq,
                         double *d, double *e,
                         magma_minproductDoubleComplex *Pt, const magma_minproduct_int_t *ldpt,
                         magma_minproductDoubleComplex *work,
                         double *resid )
                         { printf( format, __func__ ); }

void   lapackf77_zget22( const char *transa, const char *transe, const char *transw, const magma_minproduct_int_t *n,
                         magma_minproductDoubleComplex *A, const magma_minproduct_int_t *lda,
                         magma_minproductDoubleComplex *E, const magma_minproduct_int_t *lde,
                         magma_minproductDoubleComplex *wr,
                         magma_minproductDoubleComplex *wi,
                         double *work,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zhet21( magma_minproduct_int_t *itype, const char *uplo, const magma_minproduct_int_t *n, const magma_minproduct_int_t *kband,
                         magma_minproductDoubleComplex *A, const magma_minproduct_int_t *lda,
                         double *d, double *e,
                         magma_minproductDoubleComplex *U, const magma_minproduct_int_t *ldu,
                         magma_minproductDoubleComplex *V, const magma_minproduct_int_t *ldv,
                         magma_minproductDoubleComplex *tau,
                         magma_minproductDoubleComplex *work,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zhst01( const magma_minproduct_int_t *n, const magma_minproduct_int_t *ilo, const magma_minproduct_int_t *ihi,
                         magma_minproductDoubleComplex *A, const magma_minproduct_int_t *lda,
                         magma_minproductDoubleComplex *H, const magma_minproduct_int_t *ldh,
                         magma_minproductDoubleComplex *Q, const magma_minproduct_int_t *ldq,
                         magma_minproductDoubleComplex *work, const magma_minproduct_int_t *lwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zstt21( const magma_minproduct_int_t *n, const magma_minproduct_int_t *kband,
                         double *AD,
                         double *AE,
                         double *SD,
                         double *SE,
                         magma_minproductDoubleComplex *U, const magma_minproduct_int_t *ldu,
                         magma_minproductDoubleComplex *work,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zunt01( const char *rowcol, const magma_minproduct_int_t *m, const magma_minproduct_int_t *n,
                         magma_minproductDoubleComplex *U, const magma_minproduct_int_t *ldu,
                         magma_minproductDoubleComplex *work, const magma_minproduct_int_t *lwork,
                         double *resid )
                         { printf( format, __func__ ); }
#endif

void   lapackf77_zlarfy( const char *uplo, const magma_minproduct_int_t *n,
                         magma_minproductDoubleComplex *V, const magma_minproduct_int_t *incv,
                         magma_minproductDoubleComplex *tau,
                         magma_minproductDoubleComplex *C, const magma_minproduct_int_t *ldc,
                         magma_minproductDoubleComplex *work )
                         { printf( format, __func__ ); }

void   lapackf77_zlarfx( const char *side, const magma_minproduct_int_t *m, const magma_minproduct_int_t *n,
                         magma_minproductDoubleComplex *V,
                         magma_minproductDoubleComplex *tau,
                         magma_minproductDoubleComplex *C, const magma_minproduct_int_t *ldc,
                         magma_minproductDoubleComplex *work )
                         { printf( format, __func__ ); }

double lapackf77_zqpt01( const magma_minproduct_int_t *m, const magma_minproduct_int_t *n, const magma_minproduct_int_t *k,
                         magma_minproductDoubleComplex *A,
                         magma_minproductDoubleComplex *Af, const magma_minproduct_int_t *lda,
                         magma_minproductDoubleComplex *tau, magma_minproduct_int_t *jpvt,
                         magma_minproductDoubleComplex *work, const magma_minproduct_int_t *lwork )
                         { printf( format, __func__ ); return -1; }

void   lapackf77_zqrt02( const magma_minproduct_int_t *m, const magma_minproduct_int_t *n, const magma_minproduct_int_t *k,
                         magma_minproductDoubleComplex *A,
                         magma_minproductDoubleComplex *AF,
                         magma_minproductDoubleComplex *Q,
                         magma_minproductDoubleComplex *R, const magma_minproduct_int_t *lda,
                         magma_minproductDoubleComplex *tau,
                         magma_minproductDoubleComplex *work, const magma_minproduct_int_t *lwork,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_zlatms( magma_minproduct_int_t *m, magma_minproduct_int_t *n,
                         const char *dist, magma_minproduct_int_t *iseed, const char *sym, double *d,
                         magma_minproduct_int_t *mode, const double *cond, const double *dmax,
                         magma_minproduct_int_t *kl, magma_minproduct_int_t *ku, const char *pack,
                         magma_minproductDoubleComplex *a, magma_minproduct_int_t *lda, magma_minproductDoubleComplex *work, magma_minproduct_int_t *info )
                         { printf( format, __func__ ); }

#ifdef __cplusplus
}
#endif
