/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zgeqr2_gpu.cpp normal z -> c, Fri Jan 30 19:00:25 2015
       @author Stan Tomov

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
   -- Testing cgeqrf
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float           error, work[1];
    magma_minproductFloatComplex  c_neg_one = MAGMA_minproduct_C_NEG_ONE;
    magma_minproductFloatComplex *h_A, *h_R, *tau, *dtau, *h_work, tmp[1];
    magma_minproductFloatComplex_ptr d_A;
    magma_minproductFloat_ptr dwork;
    magma_minproduct_int_t M, N, n2, lda, ldda, lwork, info, min_mn;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t status = 0;

    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );

    float tol = opts.tolerance * lapackf77_slamch("E");
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    printf("  M     N     CPU GFlop/s (ms)    GPU GFlop/s (ms)    ||R||_F / ||A||_F\n");
    printf("=======================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N;
            ldda   = ((M+31)/32)*32;
            gflops = FLOPS_CGEQRF( M, N ) / 1e9;
            
            // query for workspace size
            lwork = -1;
            lapackf77_cgeqrf(&M, &N, NULL, &M, NULL, tmp, &lwork, &info);
            lwork = (magma_minproduct_int_t)MAGMA_minproduct_C_REAL( tmp[0] );
            
            TESTING_MALLOC_CPU( tau,    magma_minproductFloatComplex, min_mn );
            TESTING_MALLOC_CPU( h_A,    magma_minproductFloatComplex, n2     );
            TESTING_MALLOC_CPU( h_work, magma_minproductFloatComplex, lwork  );
            
            TESTING_MALLOC_PIN( h_R,    magma_minproductFloatComplex, n2     );
            
            TESTING_MALLOC_DEV( d_A,    magma_minproductFloatComplex, ldda*N );
            TESTING_MALLOC_DEV( dtau,   magma_minproductFloatComplex, min_mn );
            TESTING_MALLOC_DEV( dwork,  float, min_mn );
            
            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &n2, h_A );
            lapackf77_clacpy( Magma_minproductUpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );
            magma_minproduct_csetmatrix( M, N, h_R, lda, d_A, ldda );
            
            // warmup
            if ( opts.warmup ) {
                magma_minproduct_cgeqr2_gpu( M, N, d_A, ldda, dtau, dwork, &info );
                magma_minproduct_csetmatrix( M, N, h_R, lda, d_A, ldda );
            }
            
            /* ====================================================================
               Performs operation using MAGMA_minproduct
               =================================================================== */
            gpu_time = magma_minproduct_sync_wtime( 0 );

            magma_minproduct_cgeqr2_gpu( M, N, d_A, ldda, dtau, dwork, &info );

            gpu_time = magma_minproduct_sync_wtime( 0 ) - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_minproduct_cgeqr2_gpu returned error %d: %s.\n",
                       (int) info, magma_minproduct_strerror( info ));
            
            if ( opts.lapack ) {
                /* =====================================================================
                   Performs operation using LAPACK
                   =================================================================== */
                cpu_time = magma_minproduct_wtime();
                lapackf77_cgeqrf(&M, &N, h_A, &lda, tau, h_work, &lwork, &info);
                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_cgeqrf returned error %d: %s.\n",
                           (int) info, magma_minproduct_strerror( info ));
                
                /* =====================================================================
                   Check the result compared to LAPACK
                   =================================================================== */
                magma_minproduct_cgetmatrix( M, N, d_A, ldda, h_R, M );
                error = lapackf77_clange("f", &M, &N, h_A, &lda, work);
                blasf77_caxpy(&n2, &c_neg_one, h_A, &ione, h_R, &ione);
                error = lapackf77_clange("f", &M, &N, h_R, &lda, work) / error;
                
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                       (int) M, (int) N, cpu_perf, 1000.*cpu_time, gpu_perf, 1000.*gpu_time,
                       error, (error < tol ? "ok" : "failed"));
                status += ! (error < tol);
            }
            else {
                printf("%5d %5d     ---   (  ---  )   %7.2f (%7.2f)     ---  \n",
                       (int) M, (int) N, gpu_perf, 1000.*gpu_time );
            }
            
            TESTING_FREE_CPU( tau    );
            TESTING_FREE_CPU( h_A    );
            TESTING_FREE_CPU( h_work );
            
            TESTING_FREE_PIN( h_R   );
            
            TESTING_FREE_DEV( d_A   );
            TESTING_FREE_DEV( dtau  );
            TESTING_FREE_DEV( dwork );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    
    TESTING_FINALIZE();
    return status;
}
