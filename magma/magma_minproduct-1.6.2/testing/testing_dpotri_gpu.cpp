/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
  
       @generated from testing_zpotri_gpu.cpp normal z -> d, Fri Jan 30 19:00:24 2015
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dpotrf
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double *h_A, *h_R;
    magma_minproductDouble_ptr d_A;
    double c_neg_one = MAGMA_minproduct_D_NEG_ONE;
    magma_minproduct_int_t N, n2, lda, ldda, info;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    double      work[1], error;
    magma_minproduct_int_t status = 0;

    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("uplo = %s\n", lapack_uplo_const(opts.uplo) );
    printf("    N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R||_F / ||A||_F\n");
    printf("=================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            n2     = lda*N;
            ldda   = ((N+31)/32)*32;
            gflops = FLOPS_DPOTRI( N ) / 1e9;
            
            TESTING_MALLOC_CPU( h_A, double, n2 );
            TESTING_MALLOC_PIN( h_R, double, n2 );
            TESTING_MALLOC_DEV( d_A, double, ldda*N );
            
            /* Initialize the matrix */
            lapackf77_dlarnv( &ione, ISEED, &n2, h_A );
            magma_minproduct_dmake_hpd( N, h_A, lda );
            lapackf77_dlacpy( Magma_minproductUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
            
            /* ====================================================================
               Performs operation using MAGMA_minproduct
               =================================================================== */
            /* factorize matrix */
            magma_minproduct_dsetmatrix( N, N, h_A, lda, d_A, ldda );
            magma_minproduct_dpotrf_gpu( opts.uplo, N, d_A, ldda, &info );
            
            // check for exact singularity
            //magma_minproduct_dgetmatrix( N, N, d_A, ldda, h_R, lda );
            //h_R[ 10 + 10*lda ] = MAGMA_minproduct_D_MAKE( 0.0, 0.0 );
            //magma_minproduct_dsetmatrix( N, N, h_R, lda, d_A, ldda );
            
            gpu_time = magma_minproduct_wtime();
            magma_minproduct_dpotri_gpu( opts.uplo, N, d_A, ldda, &info );
            gpu_time = magma_minproduct_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_minproduct_dpotri_gpu returned error %d: %s.\n",
                       (int) info, magma_minproduct_strerror( info ));
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                lapackf77_dpotrf( lapack_uplo_const(opts.uplo), &N, h_A, &lda, &info );
                
                cpu_time = magma_minproduct_wtime();
                lapackf77_dpotri( lapack_uplo_const(opts.uplo), &N, h_A, &lda, &info );
                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_dpotri returned error %d: %s.\n",
                           (int) info, magma_minproduct_strerror( info ));
                
                /* =====================================================================
                   Check the result compared to LAPACK
                   =================================================================== */
                magma_minproduct_dgetmatrix( N, N, d_A, ldda, h_R, lda );
                error = lapackf77_dlange("f", &N, &N, h_A, &lda, work);
                blasf77_daxpy(&n2, &c_neg_one, h_A, &ione, h_R, &ione);
                error = lapackf77_dlange("f", &N, &N, h_R, &lda, work) / error;
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
            TESTING_FREE_DEV( d_A );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
