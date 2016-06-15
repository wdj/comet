/*
    -- MAGMA_tally4 (version 1.6.1) --
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
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zgeqrf
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    const double             d_neg_one = MAGMA_tally4_D_NEG_ONE;
    const double             d_one     = MAGMA_tally4_D_ONE;
    const magma_tally4DoubleComplex c_neg_one = MAGMA_tally4_Z_NEG_ONE;
    const magma_tally4DoubleComplex c_one     = MAGMA_tally4_Z_ONE;
    const magma_tally4DoubleComplex c_zero    = MAGMA_tally4_Z_ZERO;
    const magma_tally4_int_t        ione      = 1;
    
    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf=0, cpu_time=0;
    double           Anorm, error=0, error2=0;
    magma_tally4DoubleComplex *h_A, *h_R, *tau, *h_work, tmp[1];
    magma_tally4_int_t M, N, n2, lda, lwork, info, min_mn, nb;
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    
    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );

    magma_tally4_int_t status = 0;
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("ngpu %d\n", (int) opts.ngpu );
    printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |R - Q^H*A|   |I - Q^H*Q|\n");
    printf("===============================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N;
            nb     = magma_tally4_get_zgeqrf_nb(M);
            gflops = FLOPS_ZGEQRF( M, N ) / 1e9;
            
            // query for workspace size
            lwork = -1;
            lapackf77_zgeqrf(&M, &N, NULL, &M, NULL, tmp, &lwork, &info);
            lwork = (magma_tally4_int_t)MAGMA_tally4_Z_REAL( tmp[0] );
            lwork = max( lwork, max( N*nb, 2*nb*nb ));
            
            TESTING_MALLOC_CPU( tau,    magma_tally4DoubleComplex, min_mn );
            TESTING_MALLOC_CPU( h_A,    magma_tally4DoubleComplex, n2     );
            TESTING_MALLOC_CPU( h_work, magma_tally4DoubleComplex, lwork  );
            
            TESTING_MALLOC_PIN( h_R,    magma_tally4DoubleComplex, n2     );
            
            /* Initialize the matrix */
            lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
            lapackf77_zlacpy( Magma_tally4UpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );
            
            if ( opts.warmup ) {
                magma_tally4_zgeqrf(M, N, h_R, lda, tau, h_work, lwork, &info);
                lapackf77_zlacpy( Magma_tally4UpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );
            }

            /* ====================================================================
               Performs operation using MAGMA_tally4
               =================================================================== */
            gpu_time = magma_tally4_wtime();
            magma_tally4_zgeqrf(M, N, h_R, lda, tau, h_work, lwork, &info);
            gpu_time = magma_tally4_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally4_zgeqrf returned error %d: %s.\n",
                       (int) info, magma_tally4_strerror( info ));
            
            /* =====================================================================
               Check the result, following zqrt01 except using the reduced Q.
               This works for any M,N (square, tall, wide).
               =================================================================== */
            if ( opts.check ) {
                magma_tally4_int_t ldq = M;
                magma_tally4_int_t ldr = min_mn;
                magma_tally4DoubleComplex *Q, *R;
                double *work;
                TESTING_MALLOC_CPU( Q,    magma_tally4DoubleComplex, ldq*min_mn );  // M by K
                TESTING_MALLOC_CPU( R,    magma_tally4DoubleComplex, ldr*N );       // K by N
                TESTING_MALLOC_CPU( work, double,             min_mn );
                
                // generate M by K matrix Q, where K = min(M,N)
                lapackf77_zlacpy( "Lower", &M, &min_mn, h_R, &lda, Q, &ldq );
                lapackf77_zungqr( &M, &min_mn, &min_mn, Q, &ldq, tau, h_work, &lwork, &info );
                assert( info == 0 );
                
                // copy K by N matrix R
                lapackf77_zlaset( "Lower", &min_mn, &N, &c_zero, &c_zero, R, &ldr );
                lapackf77_zlacpy( "Upper", &min_mn, &N, h_R, &lda,        R, &ldr );
                
                // error = || R - Q^H*A || / (N * ||A||)
                blasf77_zgemm( "Conj", "NoTrans", &min_mn, &N, &M,
                               &c_neg_one, Q, &ldq, h_A, &lda, &c_one, R, &ldr );
                Anorm = lapackf77_zlange( "1", &M,      &N, h_A, &lda, work );
                error = lapackf77_zlange( "1", &min_mn, &N, R,   &ldr, work );
                if ( N > 0 && Anorm > 0 )
                    error /= (N*Anorm);
                
                // set R = I (K by K identity), then R = I - Q^H*Q
                // error = || I - Q^H*Q || / N
                lapackf77_zlaset( "Upper", &min_mn, &min_mn, &c_zero, &c_one, R, &ldr );
                blasf77_zherk( "Upper", "Conj", &min_mn, &M, &d_neg_one, Q, &ldq, &d_one, R, &ldr );
                error2 = lapackf77_zlanhe( "1", "Upper", &min_mn, R, &ldr, work );
                if ( N > 0 )
                    error2 /= N;
                
                TESTING_FREE_CPU( Q    );  Q    = NULL;
                TESTING_FREE_CPU( R    );  R    = NULL;
                TESTING_FREE_CPU( work );  work = NULL;
            }
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally4_wtime();
                lapackf77_zgeqrf( &M, &N, h_A, &lda, tau, h_work, &lwork, &info );
                cpu_time = magma_tally4_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_zgeqrf returned error %d: %s.\n",
                           (int) info, magma_tally4_strerror( info ));
            }
            
            /* =====================================================================
               Print performance and error.
               =================================================================== */
            printf("%5d %5d   ", (int) M, (int) N );
            if ( opts.lapack ) {
                printf( "%7.2f (%7.2f)", cpu_perf, cpu_time );
            }
            else {
                printf("  ---   (  ---  )" );
            }
            printf( "   %7.2f (%7.2f)   ", gpu_perf, gpu_time );
            if ( opts.check ) {
                bool okay = (error < tol && error2 < tol);
                status += ! okay;
                printf( "%11.2e   %11.2e   %s\n", error, error2, (okay ? "ok" : "failed") );
            }
            else {
                printf( "    ---\n" );
            }
            
            TESTING_FREE_CPU( tau    );
            TESTING_FREE_CPU( h_A    );
            TESTING_FREE_CPU( h_work );
            
            TESTING_FREE_PIN( h_R    );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
