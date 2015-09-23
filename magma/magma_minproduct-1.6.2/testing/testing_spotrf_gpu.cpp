/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zpotrf_gpu.cpp normal z -> s, Fri Jan 30 19:00:24 2015
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
   -- Testing spotrf
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float *h_A, *h_R;
    magma_minproductFloat_ptr d_A;
    magma_minproduct_int_t N, n2, lda, ldda, info;
    float c_neg_one = MAGMA_minproduct_S_NEG_ONE;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    float      work[1], error;
    magma_minproduct_int_t     status = 0;

    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    float tol = opts.tolerance * lapackf77_slamch("E");
    
    printf("uplo = %s\n", lapack_uplo_const(opts.uplo) );
    printf("  N     CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R_magma_minproduct - R_lapack||_F / ||R_lapack||_F\n");
    printf("========================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N   = opts.nsize[itest];
            lda = N;
            n2  = lda*N;
            ldda = ((N+31)/32)*32;
            gflops = FLOPS_SPOTRF( N ) / 1e9;
            
            TESTING_MALLOC_CPU( h_A, float, n2     );
            TESTING_MALLOC_PIN( h_R, float, n2     );
            TESTING_MALLOC_DEV( d_A, float, ldda*N );
            
            /* Initialize the matrix */
            lapackf77_slarnv( &ione, ISEED, &n2, h_A );
            magma_minproduct_smake_hpd( N, h_A, lda );
            lapackf77_slacpy( Magma_minproductUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
            magma_minproduct_ssetmatrix( N, N, h_A, lda, d_A, ldda );
            
            /* ====================================================================
               Performs operation using MAGMA_minproduct
               =================================================================== */
            gpu_time = magma_minproduct_wtime();
            magma_minproduct_spotrf_gpu( opts.uplo, N, d_A, ldda, &info );
            gpu_time = magma_minproduct_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_minproduct_spotrf_gpu returned error %d: %s.\n",
                       (int) info, magma_minproduct_strerror( info ));
            
            if ( opts.lapack ) {
                /* =====================================================================
                   Performs operation using LAPACK
                   =================================================================== */
                cpu_time = magma_minproduct_wtime();
                lapackf77_spotrf( lapack_uplo_const(opts.uplo), &N, h_A, &lda, &info );
                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_spotrf returned error %d: %s.\n",
                           (int) info, magma_minproduct_strerror( info ));
                
                /* =====================================================================
                   Check the result compared to LAPACK
                   =================================================================== */
                magma_minproduct_sgetmatrix( N, N, d_A, ldda, h_R, lda );
                error = lapackf77_slange("f", &N, &N, h_A, &lda, work);
                blasf77_saxpy(&n2, &c_neg_one, h_A, &ione, h_R, &ione);
                error = lapackf77_slange("f", &N, &N, h_R, &lda, work) / error;
                
                printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                       (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time,
                       error, (error < tol ? "ok" : "failed") );
                status += ! (error < tol);
            }
            else {
                printf("%5d     ---   (  ---  )   %7.2f (%7.2f)     ---  \n",
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
