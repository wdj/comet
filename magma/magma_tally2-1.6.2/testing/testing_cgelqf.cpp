/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zgelqf.cpp normal z -> c, Fri Jan 30 19:00:25 2015

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

#define COMPLEX

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cgelqf
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    const float             d_neg_one = MAGMA_tally2_D_NEG_ONE;
    const float             d_one     = MAGMA_tally2_D_ONE;
    const magma_tally2FloatComplex c_neg_one = MAGMA_tally2_C_NEG_ONE;
    const magma_tally2FloatComplex c_one     = MAGMA_tally2_C_ONE;
    const magma_tally2FloatComplex c_zero    = MAGMA_tally2_C_ZERO;
    const magma_tally2_int_t        ione      = 1;
    
    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf=0, cpu_time=0;
    float           Anorm, error=0, error2=0;
    magma_tally2FloatComplex *h_A, *h_R, *tau, *h_work, tmp[1];
    magma_tally2_int_t M, N, n2, lda, lwork, info, min_mn, nb;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    magma_tally2_int_t status = 0;

    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );

    float tol = opts.tolerance * lapackf77_slamch("E");
    
    printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |L - A*Q^H|   |I - Q*Q^H|\n");
    printf("===============================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N;
            nb     = magma_tally2_get_cgeqrf_nb(M);
            gflops = FLOPS_CGELQF( M, N ) / 1e9;

            // query for workspace size
            lwork = -1;
            lapackf77_cgelqf(&M, &N, NULL, &M, NULL, tmp, &lwork, &info);
            lwork = (magma_tally2_int_t)MAGMA_tally2_C_REAL( tmp[0] );
            lwork = max( lwork, M*nb );

            TESTING_MALLOC_CPU( tau,    magma_tally2FloatComplex, min_mn );
            TESTING_MALLOC_CPU( h_A,    magma_tally2FloatComplex, n2     );
            
            TESTING_MALLOC_PIN( h_R,    magma_tally2FloatComplex, n2     );
            TESTING_MALLOC_PIN( h_work, magma_tally2FloatComplex, lwork  );
            
            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &n2, h_A );
            lapackf77_clacpy( Magma_tally2FullStr, &M, &N, h_A, &lda, h_R,  &lda );
            
            /* ====================================================================
               Performs operation using MAGMA_tally2
               =================================================================== */
            gpu_time = magma_tally2_wtime();
            magma_tally2_cgelqf( M, N, h_R, lda, tau, h_work, lwork, &info);
            gpu_time = magma_tally2_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally2_cgelqf returned error %d: %s.\n",
                       (int) info, magma_tally2_strerror( info ));

            /* =====================================================================
               Check the result, following zlqt01 except using the reduced Q.
               This works for any M,N (square, tall, wide).
               =================================================================== */
            if ( opts.check ) {
                magma_tally2_int_t ldq = min_mn;
                magma_tally2_int_t ldl = M;
                magma_tally2FloatComplex *Q, *L;
                float *work;
                TESTING_MALLOC_CPU( Q,    magma_tally2FloatComplex, ldq*N );       // K by N
                TESTING_MALLOC_CPU( L,    magma_tally2FloatComplex, ldl*min_mn );  // M by K
                TESTING_MALLOC_CPU( work, float,             min_mn );
                
                // generate K by N matrix Q, where K = min(M,N)
                lapackf77_clacpy( "Upper", &min_mn, &N, h_R, &lda, Q, &ldq );
                lapackf77_cunglq( &min_mn, &N, &min_mn, Q, &ldq, tau, h_work, &lwork, &info );
                assert( info == 0 );
                
                // copy N by K matrix L
                lapackf77_claset( "Upper", &M, &min_mn, &c_zero, &c_zero, L, &ldl );
                lapackf77_clacpy( "Lower", &M, &min_mn, h_R, &lda,        L, &ldl );
                
                // error = || L - A*Q^H || / (N * ||A||)
                blasf77_cgemm( "NoTrans", "Conj", &M, &min_mn, &N,
                               &c_neg_one, h_A, &lda, Q, &ldq, &c_one, L, &ldl );
                Anorm = lapackf77_clange( "1", &M, &N,      h_A, &lda, work );
                error = lapackf77_clange( "1", &M, &min_mn, L,   &ldl, work );
                if ( N > 0 && Anorm > 0 )
                    error /= (N*Anorm);
                
                // set L = I (K by K), then L = I - Q*Q^H
                // error = || I - Q*Q^H || / N
                lapackf77_claset( "Upper", &min_mn, &min_mn, &c_zero, &c_one, L, &ldl );
                blasf77_cherk( "Upper", "NoTrans", &min_mn, &N, &d_neg_one, Q, &ldq, &d_one, L, &ldl );
                error2 = lapackf77_clanhe( "1", "Upper", &min_mn, L, &ldl, work );
                if ( N > 0 )
                    error2 /= N;
                
                TESTING_FREE_CPU( Q    );  Q    = NULL;
                TESTING_FREE_CPU( L    );  L    = NULL;
                TESTING_FREE_CPU( work );  work = NULL;
            }
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally2_wtime();
                lapackf77_cgelqf( &M, &N, h_A, &lda, tau, h_work, &lwork, &info );
                cpu_time = magma_tally2_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapack_cgelqf returned error %d: %s.\n",
                           (int) info, magma_tally2_strerror( info ));
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
            
            TESTING_FREE_CPU( tau );
            TESTING_FREE_CPU( h_A );
            
            TESTING_FREE_PIN( h_R    );
            TESTING_FREE_PIN( h_work );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
