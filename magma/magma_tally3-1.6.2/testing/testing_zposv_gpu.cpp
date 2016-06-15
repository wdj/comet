/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> c d s
*/
// includes, system
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zposv_gpu
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, cpu_perf, cpu_time, gpu_perf, gpu_time;
    double          error, Rnorm, Anorm, Xnorm, *work;
    magma_tally3DoubleComplex c_one     = MAGMA_tally3_Z_ONE;
    magma_tally3DoubleComplex c_neg_one = MAGMA_tally3_Z_NEG_ONE;
    magma_tally3DoubleComplex *h_A, *h_B, *h_X;
    magma_tally3DoubleComplex_ptr d_A, d_B;
    magma_tally3_int_t N, lda, ldb, ldda, lddb, info, sizeA, sizeB;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t status = 0;
    
    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("uplo = %s\n", lapack_uplo_const_tally3(opts.uplo) );
    printf("    N  NRHS   CPU Gflop/s (sec)   GPU GFlop/s (sec)   ||B - AX|| / N*||A||*||X||\n");
    printf("================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N   = opts.nsize[itest];
            lda = ldb = N;
            ldda = ((N+31)/32)*32;
            lddb = ldda;
            gflops = ( FLOPS_ZPOTRF( N ) + FLOPS_ZPOTRS( N, opts.nrhs ) ) / 1e9;
            
            TESTING_MALLOC_CPU( h_A, magma_tally3DoubleComplex, lda*N         );
            TESTING_MALLOC_CPU( h_B, magma_tally3DoubleComplex, ldb*opts.nrhs );
            TESTING_MALLOC_CPU( h_X, magma_tally3DoubleComplex, ldb*opts.nrhs );
            TESTING_MALLOC_CPU( work, double, N );
            
            TESTING_MALLOC_DEV( d_A, magma_tally3DoubleComplex, ldda*N         );
            TESTING_MALLOC_DEV( d_B, magma_tally3DoubleComplex, lddb*opts.nrhs );
            
            /* ====================================================================
               Initialize the matrix
               =================================================================== */
            sizeA = lda*N;
            sizeB = ldb*opts.nrhs;
            lapackf77_zlarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_zlarnv( &ione, ISEED, &sizeB, h_B );
            magma_tally3_zmake_hpd( N, h_A, lda );
            
            magma_tally3_zsetmatrix( N, N,         h_A, N, d_A, ldda );
            magma_tally3_zsetmatrix( N, opts.nrhs, h_B, N, d_B, lddb );
            
            /* ====================================================================
               Performs operation using MAGMA_tally3
               =================================================================== */
            gpu_time = magma_tally3_wtime();
            magma_tally3_zposv_gpu( opts.uplo, N, opts.nrhs, d_A, ldda, d_B, lddb, &info );
            gpu_time = magma_tally3_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally3_zpotrf_gpu returned error %d: %s.\n",
                       (int) info, magma_tally3_strerror( info ));

            /* =====================================================================
               Residual
               =================================================================== */
            magma_tally3_zgetmatrix( N, opts.nrhs, d_B, lddb, h_X, ldb );
            
            Anorm = lapackf77_zlange("I", &N, &N,         h_A, &lda, work);
            Xnorm = lapackf77_zlange("I", &N, &opts.nrhs, h_X, &ldb, work);
            
            blasf77_zgemm( Magma_tally3NoTransStr, Magma_tally3NoTransStr, &N, &opts.nrhs, &N,
                           &c_one,     h_A, &lda,
                                       h_X, &ldb,
                           &c_neg_one, h_B, &ldb );
            
            Rnorm = lapackf77_zlange("I", &N, &opts.nrhs, h_B, &ldb, work);
            error = Rnorm/(N*Anorm*Xnorm);
            status += ! (error < tol);
            
            /* ====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally3_wtime();
                lapackf77_zposv( lapack_uplo_const_tally3(opts.uplo), &N, &opts.nrhs, h_A, &lda, h_B, &ldb, &info );
                cpu_time = magma_tally3_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_zposv returned error %d: %s.\n",
                           (int) info, magma_tally3_strerror( info ));
                
                printf( "%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                        (int) N, (int) opts.nrhs, cpu_perf, cpu_time, gpu_perf, gpu_time,
                        error, (error < tol ? "ok" : "failed"));
            }
            else {
                printf( "%5d %5d     ---   (  ---  )   %7.2f (%7.2f)   %8.2e   %s\n",
                        (int) N, (int) opts.nrhs, gpu_perf, gpu_time,
                        error, (error < tol ? "ok" : "failed"));
            }
            
            TESTING_FREE_CPU( h_A  );
            TESTING_FREE_CPU( h_B  );
            TESTING_FREE_CPU( h_X  );
            TESTING_FREE_CPU( work );
            
            TESTING_FREE_DEV( d_A  );
            TESTING_FREE_DEV( d_B  );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
