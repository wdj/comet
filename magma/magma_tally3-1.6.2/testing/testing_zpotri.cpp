/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
  
       @precisions normal z -> c d s
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

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zpotri
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    magma_tally3DoubleComplex *h_A, *h_R;
    magma_tally3DoubleComplex c_neg_one = MAGMA_tally3_Z_NEG_ONE;
    magma_tally3_int_t N, n2, lda, info;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    double      work[1], error;
    magma_tally3_int_t status = 0;

    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("uplo = %s\n", lapack_uplo_const_tally3(opts.uplo) );
    printf("    N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R||_F / ||A||_F\n");
    printf("=================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            n2     = lda*N;
            gflops = FLOPS_ZPOTRI( N ) / 1e9;
            
            TESTING_MALLOC_CPU( h_A, magma_tally3DoubleComplex, n2 );
            TESTING_MALLOC_PIN( h_R, magma_tally3DoubleComplex, n2 );
            
            /* ====================================================================
               Initialize the matrix
               =================================================================== */
            lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
            magma_tally3_zmake_hpd( N, h_A, lda );
            lapackf77_zlacpy( Magma_tally3UpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
            
            /* ====================================================================
               Performs operation using MAGMA_tally3
               =================================================================== */
            if ( opts.warmup ) {
                magma_tally3_zpotrf( opts.uplo, N, h_R, lda, &info );
                magma_tally3_zpotri( opts.uplo, N, h_R, lda, &info );
                lapackf77_zlacpy( Magma_tally3UpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
            }
            
            /* factorize matrix */
            magma_tally3_zpotrf( opts.uplo, N, h_R, lda, &info );
            
            // check for exact singularity
            //h_R[ 10 + 10*lda ] = MAGMA_tally3_Z_MAKE( 0.0, 0.0 );
            
            gpu_time = magma_tally3_wtime();
            magma_tally3_zpotri( opts.uplo, N, h_R, lda, &info );
            gpu_time = magma_tally3_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally3_zpotri returned error %d: %s.\n",
                       (int) info, magma_tally3_strerror( info ));
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                lapackf77_zpotrf( lapack_uplo_const_tally3(opts.uplo), &N, h_A, &lda, &info );
                
                cpu_time = magma_tally3_wtime();
                lapackf77_zpotri( lapack_uplo_const_tally3(opts.uplo), &N, h_A, &lda, &info );
                cpu_time = magma_tally3_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_zpotri returned error %d: %s.\n",
                           (int) info, magma_tally3_strerror( info ));
                
                /* =====================================================================
                   Check the result compared to LAPACK
                   =================================================================== */
                error = lapackf77_zlange("f", &N, &N, h_A, &N, work);
                blasf77_zaxpy(&n2, &c_neg_one, h_A, &ione, h_R, &ione);
                error = lapackf77_zlange("f", &N, &N, h_R, &N, work) / error;
                printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                       (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time,
                       error, (error < tol ? "ok" : "failed") );
                status += ! (error < tol);
            }
            else {
                printf("%5d     ---   (  ---  )   %7.2f (%7.2f)     ---\n",
                       (int) N, gpu_perf, gpu_time );
            }
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_PIN( h_R );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
