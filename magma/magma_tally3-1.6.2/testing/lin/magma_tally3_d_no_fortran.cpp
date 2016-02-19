/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @generated from magma_tally3_z_no_fortran.cpp normal z -> d, Fri Jan 30 19:00:27 2015
       
       This is simply a copy of part of magma_tally3_dlapack.h,
       with the { printf(...); } function body added to each function.
*/
#include <stdio.h>

#include "magma_tally3.h"
#include "magma_tally3_lapack.h"

#ifdef __cplusplus
extern "C" {
#endif

#define REAL

static const char* format = "Cannot check results: %s unavailable, since there was no Fortran compiler.\n";

/*
 * Testing functions
 */
#ifdef COMPLEX
void   lapackf77_dbdt01( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *kd,
                         double *A, const magma_tally3_int_t *lda,
                         double *Q, const magma_tally3_int_t *ldq,
                         double *d, double *e,
                         double *Pt, const magma_tally3_int_t *ldpt,
                         double *work,
                         double *rwork,
                         double *resid )
                         { printf( format, __func__ ); }

void   lapackf77_dget22( const char *transa, const char *transe, const char *transw, const magma_tally3_int_t *n,
                         double *A, const magma_tally3_int_t *lda,
                         double *E, const magma_tally3_int_t *lde,
                         double *w,
                         double *work,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_dsyt21( const magma_tally3_int_t *itype, const char *uplo,
                         const magma_tally3_int_t *n, const magma_tally3_int_t *kband,
                         double *A, const magma_tally3_int_t *lda,
                         double *d, double *e,
                         double *U, const magma_tally3_int_t *ldu,
                         double *V, const magma_tally3_int_t *ldv,
                         double *tau,
                         double *work,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_dhst01( const magma_tally3_int_t *n, const magma_tally3_int_t *ilo, const magma_tally3_int_t *ihi,
                         double *A, const magma_tally3_int_t *lda,
                         double *H, const magma_tally3_int_t *ldh,
                         double *Q, const magma_tally3_int_t *ldq,
                         double *work, const magma_tally3_int_t *lwork,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_dstt21( const magma_tally3_int_t *n, const magma_tally3_int_t *kband,
                         double *AD,
                         double *AE,
                         double *SD,
                         double *SE,
                         double *U, const magma_tally3_int_t *ldu,
                         double *work,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_dort01( const char *rowcol, const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         double *U, const magma_tally3_int_t *ldu,
                         double *work, const magma_tally3_int_t *lwork,
                         double *rwork,
                         double *resid )
                         { printf( format, __func__ ); }
#else
void   lapackf77_dbdt01( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *kd,
                         double *A, const magma_tally3_int_t *lda,
                         double *Q, const magma_tally3_int_t *ldq,
                         double *d, double *e,
                         double *Pt, const magma_tally3_int_t *ldpt,
                         double *work,
                         double *resid )
                         { printf( format, __func__ ); }

void   lapackf77_dget22( const char *transa, const char *transe, const char *transw, const magma_tally3_int_t *n,
                         double *A, const magma_tally3_int_t *lda,
                         double *E, const magma_tally3_int_t *lde,
                         double *wr,
                         double *wi,
                         double *work,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_dsyt21( magma_tally3_int_t *itype, const char *uplo, const magma_tally3_int_t *n, const magma_tally3_int_t *kband,
                         double *A, const magma_tally3_int_t *lda,
                         double *d, double *e,
                         double *U, const magma_tally3_int_t *ldu,
                         double *V, const magma_tally3_int_t *ldv,
                         double *tau,
                         double *work,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_dhst01( const magma_tally3_int_t *n, const magma_tally3_int_t *ilo, const magma_tally3_int_t *ihi,
                         double *A, const magma_tally3_int_t *lda,
                         double *H, const magma_tally3_int_t *ldh,
                         double *Q, const magma_tally3_int_t *ldq,
                         double *work, const magma_tally3_int_t *lwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_dstt21( const magma_tally3_int_t *n, const magma_tally3_int_t *kband,
                         double *AD,
                         double *AE,
                         double *SD,
                         double *SE,
                         double *U, const magma_tally3_int_t *ldu,
                         double *work,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_dort01( const char *rowcol, const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         double *U, const magma_tally3_int_t *ldu,
                         double *work, const magma_tally3_int_t *lwork,
                         double *resid )
                         { printf( format, __func__ ); }
#endif

void   lapackf77_dlarfy( const char *uplo, const magma_tally3_int_t *n,
                         double *V, const magma_tally3_int_t *incv,
                         double *tau,
                         double *C, const magma_tally3_int_t *ldc,
                         double *work )
                         { printf( format, __func__ ); }

void   lapackf77_dlarfx( const char *side, const magma_tally3_int_t *m, const magma_tally3_int_t *n,
                         double *V,
                         double *tau,
                         double *C, const magma_tally3_int_t *ldc,
                         double *work )
                         { printf( format, __func__ ); }

double lapackf77_dqpt01( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         double *A,
                         double *Af, const magma_tally3_int_t *lda,
                         double *tau, magma_tally3_int_t *jpvt,
                         double *work, const magma_tally3_int_t *lwork )
                         { printf( format, __func__ ); return -1; }

void   lapackf77_dqrt02( const magma_tally3_int_t *m, const magma_tally3_int_t *n, const magma_tally3_int_t *k,
                         double *A,
                         double *AF,
                         double *Q,
                         double *R, const magma_tally3_int_t *lda,
                         double *tau,
                         double *work, const magma_tally3_int_t *lwork,
                         double *rwork,
                         double *result )
                         { printf( format, __func__ ); }

void   lapackf77_dlatms( magma_tally3_int_t *m, magma_tally3_int_t *n,
                         const char *dist, magma_tally3_int_t *iseed, const char *sym, double *d,
                         magma_tally3_int_t *mode, const double *cond, const double *dmax,
                         magma_tally3_int_t *kl, magma_tally3_int_t *ku, const char *pack,
                         double *a, magma_tally3_int_t *lda, double *work, magma_tally3_int_t *info )
                         { printf( format, __func__ ); }

#ifdef __cplusplus
}
#endif
