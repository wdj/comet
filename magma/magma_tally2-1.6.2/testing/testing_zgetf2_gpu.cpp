/*
    -- MAGMA_tally2 (version 1.6.1) --
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
#include "magma_tally2.h"
#include "magma_tally2_lapack.h"
#include "testings.h"


double get_LU_error(magma_tally2_int_t M, magma_tally2_int_t N,
                    magma_tally2DoubleComplex *A,  magma_tally2_int_t lda,
                    magma_tally2DoubleComplex *LU, magma_tally2_int_t *IPIV)
{
    magma_tally2_int_t min_mn = min(M,N);
    magma_tally2_int_t ione   = 1;
    magma_tally2_int_t i, j;
    magma_tally2DoubleComplex alpha = MAGMA_tally2_Z_ONE;
    magma_tally2DoubleComplex beta  = MAGMA_tally2_Z_ZERO;
    magma_tally2DoubleComplex *L, *U;
    double work[1], matnorm, residual;
    
    TESTING_MALLOC_CPU( L, magma_tally2DoubleComplex, M*min_mn);
    TESTING_MALLOC_CPU( U, magma_tally2DoubleComplex, min_mn*N);
    memset( L, 0, M*min_mn*sizeof(magma_tally2DoubleComplex) );
    memset( U, 0, min_mn*N*sizeof(magma_tally2DoubleComplex) );

    lapackf77_zlaswp( &N, A, &lda, &ione, &min_mn, IPIV, &ione);
    lapackf77_zlacpy( Magma_tally2LowerStr, &M, &min_mn, LU, &lda, L, &M      );
    lapackf77_zlacpy( Magma_tally2UpperStr, &min_mn, &N, LU, &lda, U, &min_mn );

    for(j=0; j<min_mn; j++)
        L[j+j*M] = MAGMA_tally2_Z_MAKE( 1., 0. );
    
    matnorm = lapackf77_zlange("f", &M, &N, A, &lda, work);

    blasf77_zgemm("N", "N", &M, &N, &min_mn,
                  &alpha, L, &M, U, &min_mn, &beta, LU, &lda);

    for( j = 0; j < N; j++ ) {
        for( i = 0; i < M; i++ ) {
            LU[i+j*lda] = MAGMA_tally2_Z_SUB( LU[i+j*lda], A[i+j*lda] );
        }
    }
    residual = lapackf77_zlange("f", &M, &N, LU, &lda, work);

    TESTING_FREE_CPU( L );
    TESTING_FREE_CPU( U );

    return residual / (matnorm * N);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zgetrf
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf=0, cpu_time=0;
    double          error;
    magma_tally2DoubleComplex *h_A, *h_R;
    magma_tally2DoubleComplex_ptr d_A;
    magma_tally2_int_t     *ipiv;
    magma_tally2_int_t M, N, n2, lda, ldda, info, min_mn;
    magma_tally2_int_t ione     = 1;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    magma_tally2_int_t status = 0;

    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );

    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("    M     N   CPU GFlop/s (ms)    GPU GFlop/s (ms)  Copy time (ms)  ||PA-LU||/(||A||*N)\n");
    printf("=======================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N;
            ldda   = ((M+31)/32)*32;
            gflops = FLOPS_ZGETRF( M, N ) / 1e9;
            
            if ( N > 512 ) {
                printf( "%5d %5d   skipping because zgetf2 does not support N > 512\n", (int) M, (int) N );
                continue;
            }
            
            TESTING_MALLOC_CPU( ipiv, magma_tally2_int_t,        min_mn );
            TESTING_MALLOC_CPU( h_A,  magma_tally2DoubleComplex, n2     );
            TESTING_MALLOC_PIN( h_R,  magma_tally2DoubleComplex, n2     );
            TESTING_MALLOC_DEV( d_A,  magma_tally2DoubleComplex, ldda*N );
            
            /* Initialize the matrix */
            lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
            lapackf77_zlacpy( Magma_tally2UpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );

            real_Double_t set_time = magma_tally2_wtime();
            magma_tally2_zsetmatrix( M, N, h_R, lda, d_A, ldda );
            set_time =  magma_tally2_wtime() - set_time;

            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally2_wtime();
                lapackf77_zgetrf(&M, &N, h_A, &lda, ipiv, &info);
                cpu_time = magma_tally2_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_zgetrf returned error %d: %s.\n",
                           (int) info, magma_tally2_strerror( info ));
            }
            
            /* ====================================================================
               Performs operation using MAGMA_tally2
               =================================================================== */
            gpu_time = magma_tally2_wtime();
            magma_tally2_zgetf2_gpu( M, N, d_A, ldda, ipiv, &info);
            gpu_time = magma_tally2_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally2_zgetf2_gpu returned error %d: %s.\n",
                       (int) info, magma_tally2_strerror( info ));
            
            real_Double_t get_time = magma_tally2_wtime();
            magma_tally2_zgetmatrix( M, N, d_A, ldda, h_A, lda );
            get_time =  magma_tally2_wtime() - get_time;

            /* =====================================================================
               Check the factorization
               =================================================================== */
            if ( opts.lapack ) {
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %7.2f",
                       (int) M, (int) N, cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000.,
                       set_time*1000.+get_time*1000.);
            }
            else {
                printf("%5d %5d     ---   (  ---  )   %7.2f (%7.2f)   %7.2f",
                       (int) M, (int) N, gpu_perf, gpu_time*1000., set_time*1000.+get_time*1000. );
            }
            if ( opts.check ) {
                magma_tally2_zgetmatrix( M, N, d_A, ldda, h_A, lda );
                error = get_LU_error( M, N, h_R, lda, h_A, ipiv );
                printf("   %8.2e   %s\n", error, (error < tol ? "ok" : "failed") );
                status += ! (error < tol);
            }
            else {
                printf("     ---  \n");
            }
            
            TESTING_FREE_CPU( ipiv );
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
