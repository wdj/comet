/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zlacpy_batched.cpp normal z -> c, Fri Jan 30 19:00:26 2015
       @author Mark Gates

*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "testings.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing clacpy_batched
   Code is very similar to testing_cgeadd_batched.cpp
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gbytes, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float           error, work[1];
    magma_minproductFloatComplex  c_neg_one = MAGMA_minproduct_C_NEG_ONE;
    magma_minproductFloatComplex *h_A, *h_B;
    magma_minproductFloatComplex_ptr d_A, d_B;
    magma_minproductFloatComplex **hAarray, **hBarray, **dAarray, **dBarray;
    magma_minproduct_int_t M, N, mb, nb, size, lda, ldda, mstride, nstride, ntile;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t status = 0;
    
    magma_minproduct_queue_t queue = magma_minproduct_stream;
    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );

    mb = (opts.nb == 0 ? 32 : opts.nb);
    nb = (opts.nb == 0 ? 64 : opts.nb);
    mstride = 2*mb;
    nstride = 3*nb;
    
    printf("mb=%d, nb=%d, mstride=%d, nstride=%d\n", (int) mb, (int) nb, (int) mstride, (int) nstride );
    printf("    M    N ntile   CPU GFlop/s (ms)    GPU GFlop/s (ms)    check\n");
    printf("=================================================================\n");
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
            gbytes = 2.*mb*nb*ntile / 1e9;
            
            TESTING_MALLOC_CPU( h_A, magma_minproductFloatComplex, lda *N );
            TESTING_MALLOC_CPU( h_B, magma_minproductFloatComplex, lda *N );
            TESTING_MALLOC_DEV( d_A, magma_minproductFloatComplex, ldda*N );
            TESTING_MALLOC_DEV( d_B, magma_minproductFloatComplex, ldda*N );
            
            TESTING_MALLOC_CPU( hAarray, magma_minproductFloatComplex*, ntile );
            TESTING_MALLOC_CPU( hBarray, magma_minproductFloatComplex*, ntile );
            TESTING_MALLOC_DEV( dAarray, magma_minproductFloatComplex*, ntile );
            TESTING_MALLOC_DEV( dBarray, magma_minproductFloatComplex*, ntile );
            
            lapackf77_clarnv( &ione, ISEED, &size, h_A );
            lapackf77_clarnv( &ione, ISEED, &size, h_B );

            /* ====================================================================
               Performs operation using MAGMA_minproduct
               =================================================================== */
            magma_minproduct_csetmatrix( M, N, h_A, lda, d_A, ldda );
            magma_minproduct_csetmatrix( M, N, h_B, lda, d_B, ldda );
            
            // setup pointers
            for( int tile = 0; tile < ntile; ++tile ) {
                int offset = tile*mstride + tile*nstride*ldda;
                hAarray[tile] = &d_A[offset];
                hBarray[tile] = &d_B[offset];
            }
            magma_minproduct_setvector( ntile, sizeof(magma_minproductFloatComplex*), hAarray, 1, dAarray, 1 );
            magma_minproduct_setvector( ntile, sizeof(magma_minproductFloatComplex*), hBarray, 1, dBarray, 1 );
            
            gpu_time = magma_minproduct_sync_wtime( 0 );
            magma_minproductblas_clacpy_batched( Magma_minproductUpperLower, mb, nb, dAarray, ldda, dBarray, ldda, ntile, queue );
            gpu_time = magma_minproduct_sync_wtime( 0 ) - gpu_time;
            gpu_perf = gbytes / gpu_time;
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_minproduct_wtime();
            for( int tile = 0; tile < ntile; ++tile ) {
                int offset = tile*mstride + tile*nstride*lda;
                lapackf77_clacpy( Magma_minproductUpperLowerStr, &mb, &nb,
                                  &h_A[offset], &lda,
                                  &h_B[offset], &lda );
            }
            cpu_time = magma_minproduct_wtime() - cpu_time;
            cpu_perf = gbytes / cpu_time;
            
            /* =====================================================================
               Check the result
               =================================================================== */
            magma_minproduct_cgetmatrix( M, N, d_B, ldda, h_A, lda );
            
            blasf77_caxpy(&size, &c_neg_one, h_A, &ione, h_B, &ione);
            error = lapackf77_clange("f", &M, &N, h_B, &lda, work);

            printf("%5d %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %s\n",
                   (int) M, (int) N, (int) ntile,
                   cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000.,
                   (error == 0. ? "ok" : "failed") );
            status += ! (error == 0.);
            
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
