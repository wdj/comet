/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zgeqp3_gpu.cpp normal z -> c, Fri Jan 30 19:00:25 2015

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"
#include "testings.h"

#define PRECISION_c

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cgeqp3_gpu
*/
int main( int argc, char** argv)
{
    TESTING_INIT();
    
    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf=0, cpu_time=0;
    magma_tally4FloatComplex *h_A, *h_R, *tau, *h_work;
    magma_tally4FloatComplex_ptr d_A, dtau, d_work;
    magma_tally4_int_t *jpvt;
    magma_tally4_int_t M, N, K, n2, lda, lwork, j, info, min_mn, nb;
    magma_tally4_int_t ione     = 1;
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    magma_tally4_int_t status = 0;
    
    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );

    float tol = opts.tolerance * lapackf77_slamch("E");
    
    printf("  M     N     CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||A*P - Q*R||_F\n");
    printf("=====================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            K = opts.ksize[itest];
            if ( M < K || N < K || K <= 0 ) {
                printf( "%5d %5d %5d   skipping because cgeqp3 requires M >= K, N >= K, K(the rank) >= 0\n",
                        (int) M, (int) N, (int) K );
                continue;
            }           
 
            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N;
            nb     = magma_tally4_get_cgeqp3_nb( min_mn );
            gflops = FLOPS_CGEQRF( M, N ) / 1e9;
            
            lwork = ( N+1 )*nb;
            #if defined(PRECISION_d) || defined(PRECISION_s)
            lwork += 2*N;
            #endif
            if ( opts.check )
                lwork = max( lwork, M*N + N );
            
            #if defined(PRECISION_z) || defined(PRECISION_c)
            float *rwork;
            magma_tally4Float_ptr drwork;
            TESTING_MALLOC_DEV( drwork, float, 2*N + (N+1)*nb);
            TESTING_MALLOC_CPU( rwork,  float, 2*N );
            #endif
            TESTING_MALLOC_CPU( jpvt,   magma_tally4_int_t,        N      );
            TESTING_MALLOC_CPU( tau,    magma_tally4FloatComplex, min_mn );
            TESTING_MALLOC_CPU( h_A,    magma_tally4FloatComplex, n2     );
            
            TESTING_MALLOC_PIN( h_R,    magma_tally4FloatComplex, n2     );
            TESTING_MALLOC_PIN( h_work, magma_tally4FloatComplex, lwork  );
            
            TESTING_MALLOC_DEV( dtau,   magma_tally4FloatComplex, min_mn );
            TESTING_MALLOC_DEV( d_A,    magma_tally4FloatComplex, lda*N  );
            TESTING_MALLOC_DEV( d_work, magma_tally4FloatComplex, lwork  );
            
            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &n2, h_R );

            /* Make h_A of rank K */
            magma_tally4FloatComplex alpha = MAGMA_tally4_C_MAKE(  1., 0. );
            magma_tally4FloatComplex beta  = MAGMA_tally4_C_MAKE(  0., 0. );
            blasf77_cgemm("N", "N", &M, &N, &K, &alpha, h_R, &lda, h_R, &lda, 
                          &beta, h_A, &lda); 

            lapackf77_clacpy( Magma_tally4UpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                for( j = 0; j < N; j++)
                    jpvt[j] = 0;
                
                cpu_time = magma_tally4_wtime();
                #if defined(PRECISION_z) || defined(PRECISION_c)
                lapackf77_cgeqp3(&M, &N, h_R, &lda, jpvt, tau, h_work, &lwork, rwork, &info);
                #else
                lapackf77_cgeqp3(&M, &N, h_R, &lda, jpvt, tau, h_work, &lwork, &info);
                #endif
                cpu_time = magma_tally4_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapack_cgeqp3 returned error %d: %s.\n",
                           (int) info, magma_tally4_strerror( info ));
            }
            
            /* ====================================================================
               Performs operation using MAGMA_tally4
               =================================================================== */
            lapackf77_clacpy( Magma_tally4UpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );
            for( j = 0; j < N; j++)
                jpvt[j] = 0;
            
            /* copy A to gpu */
            magma_tally4_csetmatrix( M, N, h_R, lda, d_A, lda );

            /* call gpu-interface */
            gpu_time = magma_tally4_wtime();
            #if defined(PRECISION_z) || defined(PRECISION_c)
            magma_tally4_cgeqp3_gpu(M, N, d_A, lda, jpvt, dtau, d_work, lwork, drwork, &info);
            #else
            magma_tally4_cgeqp3_gpu(M, N, d_A, lda, jpvt, dtau, d_work, lwork, &info);
            #endif
            gpu_time = magma_tally4_wtime() - gpu_time;
            
            /* copy outputs to cpu */
            magma_tally4_cgetmatrix( M, N, d_A, lda, h_R, lda );
            magma_tally4_cgetvector( min_mn, dtau, 1, tau, 1 );
            
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally4_cgeqp3 returned error %d: %s.\n",
                       (int) info, magma_tally4_strerror( info ));
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.lapack ) {
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)",
                       (int) M, (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time );
            }
            else {
                printf("%5d %5d     ---   (  ---  )   %7.2f (%7.2f)",
                       (int) M, (int) N, gpu_perf, gpu_time );
            }
            if ( opts.check ) {
                float error, ulp;
                ulp = lapackf77_slamch( "P" );
                
                // Compute norm( A*P - Q*R )
                error = lapackf77_cqpt01( &M, &N, &min_mn, h_A, h_R, &lda,
                                          tau, jpvt, h_work, &lwork );
                error *= ulp;
                printf("   %8.2e   %s\n", error, (error < tol ? "ok" : "failed"));
                status += ! (error < tol);
            }
            else {
                printf("     ---  \n");
            }
            
            #if defined(PRECISION_z) || defined(PRECISION_c)
            TESTING_FREE_CPU( rwork  );
            TESTING_FREE_DEV( drwork );
            #endif
            TESTING_FREE_CPU( jpvt   );
            TESTING_FREE_CPU( tau    );
            TESTING_FREE_CPU( h_A    );
            
            TESTING_FREE_PIN( h_R    );
            TESTING_FREE_PIN( h_work );
            
            TESTING_FREE_DEV( dtau   );
            TESTING_FREE_DEV( d_A    );
            TESTING_FREE_DEV( d_work );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    
    TESTING_FINALIZE();
    return status;
}