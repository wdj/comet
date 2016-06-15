/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> c d s
       @author Mark Gates

*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "testings.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zgeadd_batched
   Code is very similar to testing_zlacpy_batched.cpp
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double           error, work[1];
    magma_tally3DoubleComplex  c_neg_one = MAGMA_tally3_Z_NEG_ONE;
    magma_tally3DoubleComplex *h_A, *h_B;
    magma_tally3DoubleComplex_ptr d_A, d_B;
    magma_tally3DoubleComplex **hAarray, **hBarray, **dAarray, **dBarray;
    magma_tally3DoubleComplex alpha = MAGMA_tally3_Z_MAKE( 3.1415, 2.718 );
    magma_tally3_int_t M, N, mb, nb, size, lda, ldda, mstride, nstride, ntile;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t status = 0;
    
    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );

    double tol = opts.tolerance * lapackf77_dlamch("E");
    mb = (opts.nb == 0 ? 32 : opts.nb);
    nb = (opts.nb == 0 ? 64 : opts.nb);
    mstride = 2*mb;
    nstride = 3*nb;
    
    printf("mb=%d, nb=%d, mstride=%d, nstride=%d\n", (int) mb, (int) nb, (int) mstride, (int) nstride );
    printf("    M     N ntile   CPU GFlop/s (ms)    GPU GFlop/s (ms)    error   \n");
    printf("====================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            lda    = M;
            ldda   = ((M+31)/32)*32;
            size   = lda*N;
            
            if ( N < nb || M < nb ) {
                ntile = 0;
            } else {
                ntile = min( (M - nb)/mstride + 1,
                             (N - nb)/nstride + 1 );
            }
            gflops = 2.*mb*nb*ntile / 1e9;
            
            TESTING_MALLOC_CPU( h_A, magma_tally3DoubleComplex, lda *N );
            TESTING_MALLOC_CPU( h_B, magma_tally3DoubleComplex, lda *N );
            TESTING_MALLOC_DEV( d_A, magma_tally3DoubleComplex, ldda*N );
            TESTING_MALLOC_DEV( d_B, magma_tally3DoubleComplex, ldda*N );
            
            TESTING_MALLOC_CPU( hAarray, magma_tally3DoubleComplex*, ntile );
            TESTING_MALLOC_CPU( hBarray, magma_tally3DoubleComplex*, ntile );
            TESTING_MALLOC_DEV( dAarray, magma_tally3DoubleComplex*, ntile );
            TESTING_MALLOC_DEV( dBarray, magma_tally3DoubleComplex*, ntile );
            
            lapackf77_zlarnv( &ione, ISEED, &size, h_A );
            lapackf77_zlarnv( &ione, ISEED, &size, h_B );

            /* ====================================================================
               Performs operation using MAGMA_tally3
               =================================================================== */
            magma_tally3_zsetmatrix( M, N, h_A, lda, d_A, ldda );
            magma_tally3_zsetmatrix( M, N, h_B, lda, d_B, ldda );
            
            // setup pointers
            for( int tile = 0; tile < ntile; ++tile ) {
                int offset = tile*mstride + tile*nstride*ldda;
                hAarray[tile] = &d_A[offset];
                hBarray[tile] = &d_B[offset];
            }
            magma_tally3_setvector( ntile, sizeof(magma_tally3DoubleComplex*), hAarray, 1, dAarray, 1 );
            magma_tally3_setvector( ntile, sizeof(magma_tally3DoubleComplex*), hBarray, 1, dBarray, 1 );
            
            gpu_time = magma_tally3_sync_wtime( 0 );
            magma_tally3blas_zgeadd_batched( mb, nb, alpha, dAarray, ldda, dBarray, ldda, ntile );
            gpu_time = magma_tally3_sync_wtime( 0 ) - gpu_time;
            gpu_perf = gflops / gpu_time;
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_tally3_wtime();
            for( int tile = 0; tile < ntile; ++tile ) {
                int offset = tile*mstride + tile*nstride*lda;
                for( int j = 0; j < nb; ++j ) {
                    blasf77_zaxpy( &mb, &alpha,
                                   &h_A[offset + j*lda], &ione,
                                   &h_B[offset + j*lda], &ione );
                }
            }
            cpu_time = magma_tally3_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            
            /* =====================================================================
               Check the result
               =================================================================== */
            magma_tally3_zgetmatrix( M, N, d_B, ldda, h_A, lda );
            
            error = lapackf77_zlange( "F", &M, &N, h_B, &lda, work );
            blasf77_zaxpy(&size, &c_neg_one, h_A, &ione, h_B, &ione);
            error = lapackf77_zlange("f", &M, &N, h_B, &lda, work) / error;

            printf("%5d %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                   (int) M, (int) N, (int) ntile,
                   cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000.,
                   error, (error < tol ? "ok" : "failed"));
            status += ! (error < tol);
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_B );
            TESTING_FREE_DEV( d_A );
            TESTING_FREE_DEV( d_B );
            
            TESTING_FREE_CPU( hAarray );
            TESTING_FREE_CPU( hBarray );
            TESTING_FREE_DEV( dAarray );
            TESTING_FREE_DEV( dBarray );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
