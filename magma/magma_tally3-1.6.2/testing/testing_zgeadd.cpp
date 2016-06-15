/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> c d s
       @author Mark Gates
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "testings.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zgeadd
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double          error, work[1];
    magma_tally3DoubleComplex *h_A, *h_B, *d_A, *d_B;
    magma_tally3DoubleComplex alpha = MAGMA_tally3_Z_MAKE( 3.1415, 2.718 );
    magma_tally3DoubleComplex c_neg_one = MAGMA_tally3_Z_NEG_ONE;
    
    magma_tally3_int_t M, N, size, lda, ldda;
    magma_tally3_int_t ione = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t status = 0;

    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );

    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    /* Uncomment these lines to check parameters.
     * magma_tally3_xerbla calls lapack's xerbla to print out error. */
    //magma_tally3blas_zgeadd( -1,  N, alpha, d_A, ldda, d_B, ldda );
    //magma_tally3blas_zgeadd(  M, -1, alpha, d_A, ldda, d_B, ldda );
    //magma_tally3blas_zgeadd(  M,  N, alpha, d_A, M-1,  d_B, ldda );
    //magma_tally3blas_zgeadd(  M,  N, alpha, d_A, ldda, d_B, N-1  );

    printf("    M     N   CPU GFlop/s (ms)    GPU GFlop/s (ms)    |Bl-Bm|/|Bl|\n");
    printf("=========================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            lda    = M;
            ldda   = ((M+31)/32)*32;
            size   = lda*N;
            gflops = 2.*M*N / 1e9;
            
            TESTING_MALLOC_CPU( h_A, magma_tally3DoubleComplex, lda *N );
            TESTING_MALLOC_CPU( h_B, magma_tally3DoubleComplex, lda *N );
            
            TESTING_MALLOC_DEV( d_A, magma_tally3DoubleComplex, ldda*N );
            TESTING_MALLOC_DEV( d_B, magma_tally3DoubleComplex, ldda*N );
            
            lapackf77_zlarnv( &ione, ISEED, &size, h_A );
            lapackf77_zlarnv( &ione, ISEED, &size, h_B );
            
            /* ====================================================================
               Performs operation using MAGMA_tally3
               =================================================================== */
            magma_tally3_zsetmatrix( M, N, h_A, lda, d_A, ldda );
            magma_tally3_zsetmatrix( M, N, h_B, lda, d_B, ldda );
            
            gpu_time = magma_tally3_sync_wtime( NULL );
            magma_tally3blas_zgeadd( M, N, alpha, d_A, ldda, d_B, ldda );
            gpu_time = magma_tally3_sync_wtime( NULL ) - gpu_time;
            gpu_perf = gflops / gpu_time;
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_tally3_wtime();
            for( int j = 0; j < N; ++j ) {
                blasf77_zaxpy( &M, &alpha, &h_A[j*lda], &ione, &h_B[j*lda], &ione );
            }
            cpu_time = magma_tally3_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            
            /* =====================================================================
               Check result
               =================================================================== */
            magma_tally3_zgetmatrix( M, N, d_B, ldda, h_A, lda );
            
            error = lapackf77_zlange( "F", &M, &N, h_B, &lda, work );
            blasf77_zaxpy( &size, &c_neg_one, h_A, &ione, h_B, &ione );
            error = lapackf77_zlange( "F", &M, &N, h_B, &lda, work ) / error;
            
            printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                   (int) M, (int) N,
                   cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000.,
                   error, (error < tol ? "ok" : "failed"));
            status += ! (error < tol);
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_B );
            
            TESTING_FREE_DEV( d_A );
            TESTING_FREE_DEV( d_B );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
