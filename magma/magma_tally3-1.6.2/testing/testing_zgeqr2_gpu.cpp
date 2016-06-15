/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
       @author Stan Tomov

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
   -- Testing zgeqrf
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double           error, work[1];
    magma_tally3DoubleComplex  c_neg_one = MAGMA_tally3_Z_NEG_ONE;
    magma_tally3DoubleComplex *h_A, *h_R, *tau, *dtau, *h_work, tmp[1];
    magma_tally3DoubleComplex_ptr d_A;
    magma_tally3Double_ptr dwork;
    magma_tally3_int_t M, N, n2, lda, ldda, lwork, info, min_mn;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t status = 0;

    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );

    double tol = opts.tolerance * lapackf77_dlamch("E");
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
            gflops = FLOPS_ZGEQRF( M, N ) / 1e9;
            
            // query for workspace size
            lwork = -1;
            lapackf77_zgeqrf(&M, &N, NULL, &M, NULL, tmp, &lwork, &info);
            lwork = (magma_tally3_int_t)MAGMA_tally3_Z_REAL( tmp[0] );
            
            TESTING_MALLOC_CPU( tau,    magma_tally3DoubleComplex, min_mn );
            TESTING_MALLOC_CPU( h_A,    magma_tally3DoubleComplex, n2     );
            TESTING_MALLOC_CPU( h_work, magma_tally3DoubleComplex, lwork  );
            
            TESTING_MALLOC_PIN( h_R,    magma_tally3DoubleComplex, n2     );
            
            TESTING_MALLOC_DEV( d_A,    magma_tally3DoubleComplex, ldda*N );
            TESTING_MALLOC_DEV( dtau,   magma_tally3DoubleComplex, min_mn );
            TESTING_MALLOC_DEV( dwork,  double, min_mn );
            
            /* Initialize the matrix */
            lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
            lapackf77_zlacpy( Magma_tally3UpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );
            magma_tally3_zsetmatrix( M, N, h_R, lda, d_A, ldda );
            
            // warmup
            if ( opts.warmup ) {
                magma_tally3_zgeqr2_gpu( M, N, d_A, ldda, dtau, dwork, &info );
                magma_tally3_zsetmatrix( M, N, h_R, lda, d_A, ldda );
            }
            
            /* ====================================================================
               Performs operation using MAGMA_tally3
               =================================================================== */
            gpu_time = magma_tally3_sync_wtime( 0 );

            magma_tally3_zgeqr2_gpu( M, N, d_A, ldda, dtau, dwork, &info );

            gpu_time = magma_tally3_sync_wtime( 0 ) - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally3_zgeqr2_gpu returned error %d: %s.\n",
                       (int) info, magma_tally3_strerror( info ));
            
            if ( opts.lapack ) {
                /* =====================================================================
                   Performs operation using LAPACK
                   =================================================================== */
                cpu_time = magma_tally3_wtime();
                lapackf77_zgeqrf(&M, &N, h_A, &lda, tau, h_work, &lwork, &info);
                cpu_time = magma_tally3_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_zgeqrf returned error %d: %s.\n",
                           (int) info, magma_tally3_strerror( info ));
                
                /* =====================================================================
                   Check the result compared to LAPACK
                   =================================================================== */
                magma_tally3_zgetmatrix( M, N, d_A, ldda, h_R, M );
                error = lapackf77_zlange("f", &M, &N, h_A, &lda, work);
                blasf77_zaxpy(&n2, &c_neg_one, h_A, &ione, h_R, &ione);
                error = lapackf77_zlange("f", &M, &N, h_R, &lda, work) / error;
                
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
