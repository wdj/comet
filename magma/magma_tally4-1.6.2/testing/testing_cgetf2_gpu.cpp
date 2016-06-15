/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zgetf2_gpu.cpp normal z -> c, Fri Jan 30 19:00:24 2015
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


float get_LU_error(magma_tally4_int_t M, magma_tally4_int_t N,
                    magma_tally4FloatComplex *A,  magma_tally4_int_t lda,
                    magma_tally4FloatComplex *LU, magma_tally4_int_t *IPIV)
{
    magma_tally4_int_t min_mn = min(M,N);
    magma_tally4_int_t ione   = 1;
    magma_tally4_int_t i, j;
    magma_tally4FloatComplex alpha = MAGMA_tally4_C_ONE;
    magma_tally4FloatComplex beta  = MAGMA_tally4_C_ZERO;
    magma_tally4FloatComplex *L, *U;
    float work[1], matnorm, residual;
    
    TESTING_MALLOC_CPU( L, magma_tally4FloatComplex, M*min_mn);
    TESTING_MALLOC_CPU( U, magma_tally4FloatComplex, min_mn*N);
    memset( L, 0, M*min_mn*sizeof(magma_tally4FloatComplex) );
    memset( U, 0, min_mn*N*sizeof(magma_tally4FloatComplex) );

    lapackf77_claswp( &N, A, &lda, &ione, &min_mn, IPIV, &ione);
    lapackf77_clacpy( Magma_tally4LowerStr, &M, &min_mn, LU, &lda, L, &M      );
    lapackf77_clacpy( Magma_tally4UpperStr, &min_mn, &N, LU, &lda, U, &min_mn );

    for(j=0; j<min_mn; j++)
        L[j+j*M] = MAGMA_tally4_C_MAKE( 1., 0. );
    
    matnorm = lapackf77_clange("f", &M, &N, A, &lda, work);

    blasf77_cgemm("N", "N", &M, &N, &min_mn,
                  &alpha, L, &M, U, &min_mn, &beta, LU, &lda);

    for( j = 0; j < N; j++ ) {
        for( i = 0; i < M; i++ ) {
            LU[i+j*lda] = MAGMA_tally4_C_SUB( LU[i+j*lda], A[i+j*lda] );
        }
    }
    residual = lapackf77_clange("f", &M, &N, LU, &lda, work);

    TESTING_FREE_CPU( L );
    TESTING_FREE_CPU( U );

    return residual / (matnorm * N);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cgetrf
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf=0, cpu_time=0;
    float          error;
    magma_tally4FloatComplex *h_A, *h_R;
    magma_tally4FloatComplex_ptr d_A;
    magma_tally4_int_t     *ipiv;
    magma_tally4_int_t M, N, n2, lda, ldda, info, min_mn;
    magma_tally4_int_t ione     = 1;
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    magma_tally4_int_t status = 0;

    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );

    float tol = opts.tolerance * lapackf77_slamch("E");
    
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
            gflops = FLOPS_CGETRF( M, N ) / 1e9;
            
            if ( N > 512 ) {
                printf( "%5d %5d   skipping because cgetf2 does not support N > 512\n", (int) M, (int) N );
                continue;
            }
            
            TESTING_MALLOC_CPU( ipiv, magma_tally4_int_t,        min_mn );
            TESTING_MALLOC_CPU( h_A,  magma_tally4FloatComplex, n2     );
            TESTING_MALLOC_PIN( h_R,  magma_tally4FloatComplex, n2     );
            TESTING_MALLOC_DEV( d_A,  magma_tally4FloatComplex, ldda*N );
            
            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &n2, h_A );
            lapackf77_clacpy( Magma_tally4UpperLowerStr, &M, &N, h_A, &lda, h_R, &lda );

            real_Double_t set_time = magma_tally4_wtime();
            magma_tally4_csetmatrix( M, N, h_R, lda, d_A, ldda );
            set_time =  magma_tally4_wtime() - set_time;

            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally4_wtime();
                lapackf77_cgetrf(&M, &N, h_A, &lda, ipiv, &info);
                cpu_time = magma_tally4_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_cgetrf returned error %d: %s.\n",
                           (int) info, magma_tally4_strerror( info ));
            }
            
            /* ====================================================================
               Performs operation using MAGMA_tally4
               =================================================================== */
            gpu_time = magma_tally4_wtime();
            magma_tally4_cgetf2_gpu( M, N, d_A, ldda, ipiv, &info);
            gpu_time = magma_tally4_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally4_cgetf2_gpu returned error %d: %s.\n",
                       (int) info, magma_tally4_strerror( info ));
            
            real_Double_t get_time = magma_tally4_wtime();
            magma_tally4_cgetmatrix( M, N, d_A, ldda, h_A, lda );
            get_time =  magma_tally4_wtime() - get_time;

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
                magma_tally4_cgetmatrix( M, N, d_A, ldda, h_A, lda );
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