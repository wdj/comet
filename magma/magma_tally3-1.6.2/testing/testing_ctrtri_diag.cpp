/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_ztrtri_diag.cpp normal z -> c, Fri Jan 30 19:00:23 2015
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "testings.h"

#define     h_A(i_, j_) (h_A     + (i_) + (j_)*lda)
#define h_dinvA(i_, j_) (h_dinvA + (i_) + (j_)*nb)


/* ////////////////////////////////////////////////////////////////////////////
   -- like axpy for matrices: B += alpha*A.
*/
void cgeadd(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    const magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex       *B, magma_tally3_int_t ldb )
{
    #define A(i_, j_) (A + (i_) + (j_)*lda)
    #define B(i_, j_) (B + (i_) + (j_)*ldb)
    
    const magma_tally3_int_t ione = 1;
    
    for( int j=0; j < n; ++j ) {
        blasf77_caxpy( &m, &alpha, A(0,j), &ione, B(0,j), &ione );
    }
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing ctrtri
*/
int main( int argc, char** argv )
{
    TESTING_INIT();

    real_Double_t   gflops, magma_tally3_perf, magma_tally3_time=0;  //, cpu_perf=0, cpu_time=0;
    float          magma_tally3_error, norm_invA, work[1];
    magma_tally3_int_t N, lda, ldda, info;
    magma_tally3_int_t jb, nb, nblock, sizeA, size_inv;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t *ipiv;

    magma_tally3FloatComplex *h_A, *h_dinvA;
    magma_tally3FloatComplex_ptr d_A, d_dinvA;
    magma_tally3FloatComplex c_neg_one = MAGMA_tally3_C_NEG_ONE;
    magma_tally3_int_t status = 0;
    
    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    float tol = opts.tolerance * lapackf77_slamch("E");
    const char *uplo_ = lapack_uplo_const_tally3(opts.uplo);

    // this is the NB hard coded into ctrtri_diag.
    nb = 128;
    
    printf("uplo = %s, diag = %s\n",
           lapack_uplo_const_tally3(opts.uplo), lapack_diag_const_tally3(opts.diag) );
    printf("    N  MAGMA_tally3 Gflop/s (ms)   MAGMA_tally3 error\n");
    printf("=======================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda = N;
            ldda = ((lda+31)/32)*32;
            nblock = (N+nb-1)/nb;
            gflops = nblock * FLOPS_CTRTRI( nb ) / 1e9;
            
            TESTING_MALLOC_CPU( h_A,    magma_tally3FloatComplex, lda*N );
            TESTING_MALLOC_CPU( ipiv,   magma_tally3_int_t,        N     );
            
            size_inv = nblock*nb*nb;
            TESTING_MALLOC_DEV( d_A,    magma_tally3FloatComplex, ldda*N );
            TESTING_MALLOC_DEV( d_dinvA, magma_tally3FloatComplex, size_inv );
            TESTING_MALLOC_CPU( h_dinvA, magma_tally3FloatComplex, size_inv );
            
            /* Initialize the matrices */
            /* Factor A into LU to get well-conditioned triangular matrix.
             * Copy L to U, since L seems okay when used with non-unit diagonal
             * (i.e., from U), while U fails when used with unit diagonal. */
            sizeA = lda*N;            
            lapackf77_clarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_cgetrf( &N, &N, h_A, &lda, ipiv, &info );
            for( int j = 0; j < N; ++j ) {
                for( int i = 0; i < j; ++i ) {
                    *h_A(i,j) = *h_A(j,i);
                }
            }
            
            /* =====================================================================
               Performs operation using MAGMA_tally3BLAS
               =================================================================== */
            magma_tally3_csetmatrix( N, N, h_A, lda, d_A, ldda );
            
            magma_tally3_time = magma_tally3_sync_wtime( NULL );
            magma_tally3blas_ctrtri_diag( opts.uplo, opts.diag, N, d_A, ldda, d_dinvA );
            magma_tally3_time = magma_tally3_sync_wtime( NULL ) - magma_tally3_time;
            magma_tally3_perf = gflops / magma_tally3_time;
            
            magma_tally3_cgetvector( size_inv, d_dinvA, 1, h_dinvA, 1 );
            
            if ( opts.verbose ) {
                printf( "A%d=", (int) N );
                magma_tally3_cprint( N, N, h_A, lda );
                printf( "d_dinvA%d=", (int) N );
                magma_tally3_cprint( min(N+4, nb), min(N+4, nblock*nb), h_dinvA, nb );
            }
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                //cpu_time = magma_tally3_wtime();
                lapackf77_ctrtri(
                    lapack_uplo_const_tally3(opts.uplo), lapack_diag_const_tally3(opts.diag),
                    &N, h_A, &lda, &info );
                //cpu_time = magma_tally3_wtime() - cpu_time;
                //cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.check ) {
                // |invA - invA_magma_tally3| / |invA|, accumulated over all diagonal blocks
                magma_tally3_error = 0;
                norm_invA   = 0;
                for( int i=0; i < N; i += nb ) {
                    jb = min( nb, N-i );
                    cgeadd( jb, jb, c_neg_one, h_A(i, i), lda, h_dinvA(0, i), nb );
                    magma_tally3_error = max( magma_tally3_error, lapackf77_clantr( "M", uplo_, Magma_tally3NonUnitStr, &jb, &jb, h_dinvA(0, i), &nb,  work ));
                    norm_invA   = max( norm_invA,   lapackf77_clantr( "M", uplo_, Magma_tally3NonUnitStr, &jb, &jb, h_A(i, i),     &lda, work ));
                }
                magma_tally3_error /= norm_invA;
                
                // CPU is doing N-by-N inverse, while GPU is doing (N/NB) NB-by-NB inverses.
                // So don't compare performance.
                printf("%5d   %7.2f (%7.2f)   %8.2e   %s\n",
                        (int) N,
                        magma_tally3_perf,  1000.*magma_tally3_time,
                        //cpu_perf,    1000.*cpu_time,
                        magma_tally3_error,
                        (magma_tally3_error < tol ? "ok" : "failed"));
                status += ! (magma_tally3_error < tol);
            }
            else {
                printf("%5d   %7.2f (%7.2f)      ---\n",
                        (int) N,
                        magma_tally3_perf,  1000.*magma_tally3_time );
            }
            
            TESTING_FREE_CPU( h_A     );
            TESTING_FREE_CPU( ipiv    );
            
            TESTING_FREE_DEV( d_A     );
            TESTING_FREE_DEV( d_dinvA );
            TESTING_FREE_CPU( h_dinvA );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
