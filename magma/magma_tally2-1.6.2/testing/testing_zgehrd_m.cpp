/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_tally2.h"
#include "magma_tally2_lapack.h"
#include "testings.h"

#define PRECISION_z

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zgehrd_m
*/
int main( int argc, char** argv)
{
    TESTING_INIT();
    magma_tally2_setdevice( 0 );  // without this, T -> dT copy fails
    
    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    magma_tally2DoubleComplex *h_A, *h_R, *h_Q, *h_work, *tau, *twork, *T, *dT;
    #if defined(PRECISION_z) || defined(PRECISION_c)
    double      *rwork;
    #endif
    double      eps, result[2];
    magma_tally2_int_t N, n2, lda, nb, lwork, ltwork, info;
    magma_tally2_int_t ione     = 1;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    magma_tally2_int_t status = 0;
    
    eps   = lapackf77_dlamch( "E" );
    
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );
    
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("    N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |A-QHQ'|/N|A|   |I-QQ'|/N\n");
    printf("=========================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            n2     = lda*N;
            nb     = magma_tally2_get_zgehrd_nb(N);
            // magma_tally2 needs larger workspace than lapack, esp. multi-gpu verison
            lwork  = N*(nb + nb*Magma_tally2MaxGPUs);
            gflops = FLOPS_ZGEHRD( N ) / 1e9;
            
            TESTING_MALLOC_CPU( h_A,    magma_tally2DoubleComplex, n2    );
            TESTING_MALLOC_CPU( tau,    magma_tally2DoubleComplex, N     );
            TESTING_MALLOC_CPU( T,      magma_tally2DoubleComplex, nb*N  );
            
            TESTING_MALLOC_PIN( h_R,    magma_tally2DoubleComplex, n2    );
            TESTING_MALLOC_PIN( h_work, magma_tally2DoubleComplex, lwork );
            
            TESTING_MALLOC_DEV( dT,     magma_tally2DoubleComplex, nb*N  );
            
            /* Initialize the matrices */
            lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
            lapackf77_zlacpy( Magma_tally2UpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
            
            /* ====================================================================
               Performs operation using MAGMA_tally2
               =================================================================== */
            gpu_time = magma_tally2_wtime();
            magma_tally2_zgehrd_m( N, ione, N, h_R, lda, tau, h_work, lwork, T, &info);
            gpu_time = magma_tally2_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally2_zgehrd_m returned error %d: %s.\n",
                       (int) info, magma_tally2_strerror( info ));
            
            /* =====================================================================
               Check the factorization
               =================================================================== */
            if ( opts.check ) {
                ltwork = 2*(N*N);
                TESTING_MALLOC_PIN( h_Q,   magma_tally2DoubleComplex, lda*N  );
                TESTING_MALLOC_CPU( twork, magma_tally2DoubleComplex, ltwork );
                #if defined(PRECISION_z) || defined(PRECISION_c)
                TESTING_MALLOC_CPU( rwork, double, N );
                #endif
                
                lapackf77_zlacpy(Magma_tally2UpperLowerStr, &N, &N, h_R, &lda, h_Q, &lda);
                for( int j = 0; j < N-1; ++j )
                    for( int i = j+2; i < N; ++i )
                        h_R[i+j*lda] = MAGMA_tally2_Z_ZERO;
                
                magma_tally2_zsetmatrix( nb, N, T, nb, dT, nb );
                
                magma_tally2_zunghr(N, ione, N, h_Q, lda, tau, dT, nb, &info);
                if ( info != 0 ) {
                    printf("magma_tally2_zunghr returned error %d: %s.\n",
                           (int) info, magma_tally2_strerror( info ));
                    exit(1);
                }
                #if defined(PRECISION_z) || defined(PRECISION_c)
                lapackf77_zhst01(&N, &ione, &N,
                                 h_A, &lda, h_R, &lda,
                                 h_Q, &lda, twork, &ltwork, rwork, result);
                #else
                lapackf77_zhst01(&N, &ione, &N,
                                 h_A, &lda, h_R, &lda,
                                 h_Q, &lda, twork, &ltwork, result);
                #endif
                
                TESTING_FREE_PIN( h_Q   );
                TESTING_FREE_CPU( twork );
                #if defined(PRECISION_z) || defined(PRECISION_c)
                TESTING_FREE_CPU( rwork );
                #endif
            }
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally2_wtime();
                lapackf77_zgehrd(&N, &ione, &N, h_R, &lda, tau, h_work, &lwork, &info);
                cpu_time = magma_tally2_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_zgehrd returned error %d: %s.\n",
                           (int) info, magma_tally2_strerror( info ));
            }
            
            /* =====================================================================
               Print performance and error.
               =================================================================== */
            if ( opts.lapack ) {
                printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)",
                       (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time );
            }
            else {
                printf("%5d     ---   (  ---  )   %7.2f (%7.2f)",
                       (int) N, gpu_perf, gpu_time );
            }
            if ( opts.check ) {
                printf("   %8.2e        %8.2e   %s\n",
                       result[0]*eps, result[1]*eps,
                       ( ( (result[0]*eps < tol) && (result[1]*eps < tol) ) ? "ok" : "failed")  );
                status += ! (result[0]*eps < tol);
                status += ! (result[1]*eps < tol);
            }
            else {
                printf("     ---             ---\n");
            }
            
            TESTING_FREE_CPU( h_A    );
            TESTING_FREE_CPU( tau    );
            TESTING_FREE_CPU( T      );
            
            TESTING_FREE_PIN( h_R    );
            TESTING_FREE_PIN( h_work );
            
            TESTING_FREE_DEV( dT     );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    
    TESTING_FINALIZE();
    return status;
}
