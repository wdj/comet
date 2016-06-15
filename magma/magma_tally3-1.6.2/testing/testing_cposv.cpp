/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zposv.cpp normal z -> c, Fri Jan 30 19:00:24 2015
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
   -- Testing cposv
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, cpu_perf, cpu_time, gpu_perf, gpu_time;
    float          error, Rnorm, Anorm, Xnorm, *work;
    magma_tally3FloatComplex c_one     = MAGMA_tally3_C_ONE;
    magma_tally3FloatComplex c_neg_one = MAGMA_tally3_C_NEG_ONE;
    magma_tally3FloatComplex *h_A, *h_R, *h_B, *h_X;
    magma_tally3_int_t N, lda, ldb, info, sizeA, sizeB;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t status = 0;
    
    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");
    
    printf("ngpu = %d, uplo = %s\n", (int) opts.ngpu, lapack_uplo_const_tally3(opts.uplo) );
    printf("    N  NRHS   CPU Gflop/s (sec)   GPU GFlop/s (sec)   ||B - AX|| / N*||A||*||X||\n");
    printf("================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N   = opts.nsize[itest];
            lda = ldb = N;
            gflops = ( FLOPS_CPOTRF( N ) + FLOPS_CPOTRS( N, opts.nrhs ) ) / 1e9;
            
            TESTING_MALLOC_CPU( h_A, magma_tally3FloatComplex, lda*N         );
            TESTING_MALLOC_CPU( h_R, magma_tally3FloatComplex, lda*N         );
            TESTING_MALLOC_CPU( h_B, magma_tally3FloatComplex, ldb*opts.nrhs );
            TESTING_MALLOC_CPU( h_X, magma_tally3FloatComplex, ldb*opts.nrhs );
            TESTING_MALLOC_CPU( work, float, N );
            
            /* ====================================================================
               Initialize the matrix
               =================================================================== */
            sizeA = lda*N;
            sizeB = ldb*opts.nrhs;
            lapackf77_clarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_clarnv( &ione, ISEED, &sizeB, h_B );
            magma_tally3_cmake_hpd( N, h_A, lda );
            
            // copy A to R and B to X; save A and B for residual
            lapackf77_clacpy( Magma_tally3UpperLowerStr, &N, &N,         h_A, &lda, h_R, &lda );
            lapackf77_clacpy( Magma_tally3UpperLowerStr, &N, &opts.nrhs, h_B, &ldb, h_X, &ldb );
            
            /* ====================================================================
               Performs operation using MAGMA_tally3
               =================================================================== */
            gpu_time = magma_tally3_wtime();
            magma_tally3_cposv( opts.uplo, N, opts.nrhs, h_R, lda, h_X, ldb, &info );
            gpu_time = magma_tally3_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally3_cpotrf returned error %d: %s.\n",
                       (int) info, magma_tally3_strerror( info ));
            
            /* =====================================================================
               Residual
               =================================================================== */
            Anorm = lapackf77_clange("I", &N, &N,         h_A, &lda, work);
            Xnorm = lapackf77_clange("I", &N, &opts.nrhs, h_X, &ldb, work);
            
            blasf77_cgemm( Magma_tally3NoTransStr, Magma_tally3NoTransStr, &N, &opts.nrhs, &N,
                           &c_one,     h_A, &lda,
                                       h_X, &ldb,
                           &c_neg_one, h_B, &ldb );
            
            Rnorm = lapackf77_clange("I", &N, &opts.nrhs, h_B, &ldb, work);
            error = Rnorm/(N*Anorm*Xnorm);
            status += ! (error < tol);
            
            /* ====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally3_wtime();
                lapackf77_cposv( lapack_uplo_const_tally3(opts.uplo), &N, &opts.nrhs, h_A, &lda, h_B, &ldb, &info );
                cpu_time = magma_tally3_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_cposv returned error %d: %s.\n",
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
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_R );
            TESTING_FREE_CPU( h_B );
            TESTING_FREE_CPU( h_X );
            TESTING_FREE_CPU( work );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
