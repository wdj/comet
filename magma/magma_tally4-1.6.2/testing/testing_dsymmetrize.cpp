/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zsymmetrize.cpp normal z -> d, Fri Jan 30 19:00:24 2015
       @author Mark Gates

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dsymmetrize
   Code is very similar to testing_dtranspose.cpp
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gbytes, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double           error, work[1];
    double  c_neg_one = MAGMA_tally4_D_NEG_ONE;
    double *h_A, *h_R;
    magma_tally4Double_ptr d_A;
    magma_tally4_int_t N, size, lda, ldda;
    magma_tally4_int_t ione     = 1;
    magma_tally4_int_t status = 0;
    
    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );

    printf("uplo = %s\n", lapack_uplo_const_tally4(opts.uplo) );
    printf("    N   CPU GByte/s (ms)    GPU GByte/s (ms)    check\n");
    printf("=====================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            ldda   = ((N+31)/32)*32;
            size   = lda*N;
            // load strictly lower triangle, save strictly upper triangle
            gbytes = sizeof(double) * 1.*N*(N-1) / 1e9;
    
            TESTING_MALLOC_CPU( h_A, double, size   );
            TESTING_MALLOC_CPU( h_R, double, size   );
            
            TESTING_MALLOC_DEV( d_A, double, ldda*N );
            
            /* Initialize the matrix */
            for( int j = 0; j < N; ++j ) {
                for( int i = 0; i < N; ++i ) {
                    h_A[i + j*lda] = MAGMA_tally4_D_MAKE( i + j/10000., j );
                }
            }
            
            /* ====================================================================
               Performs operation using MAGMA_tally4
               =================================================================== */
            magma_tally4_dsetmatrix( N, N, h_A, lda, d_A, ldda );
            
            gpu_time = magma_tally4_sync_wtime( 0 );
            //magma_tally4blas_dsymmetrize( opts.uplo, N-2, d_A+1+ldda, ldda );  // inset by 1 row & col
            magma_tally4blas_dsymmetrize( opts.uplo, N, d_A, ldda );
            gpu_time = magma_tally4_sync_wtime( 0 ) - gpu_time;
            gpu_perf = gbytes / gpu_time;
            
            /* =====================================================================
               Performs operation using naive in-place algorithm
               (LAPACK doesn't implement symmetrize)
               =================================================================== */
            cpu_time = magma_tally4_wtime();
            //for( int j = 1; j < N-1; ++j ) {    // inset by 1 row & col
            //    for( int i = 1; i < j; ++i ) {
            for( int j = 0; j < N; ++j ) {
                for( int i = 0; i < j; ++i ) {
                    if ( opts.uplo == Magma_tally4Lower ) {
                        h_A[i + j*lda] = MAGMA_tally4_D_CNJG( h_A[j + i*lda] );
                    }
                    else {
                        h_A[j + i*lda] = MAGMA_tally4_D_CNJG( h_A[i + j*lda] );
                    }
                }
            }
            cpu_time = magma_tally4_wtime() - cpu_time;
            cpu_perf = gbytes / cpu_time;
            
            /* =====================================================================
               Check the result
               =================================================================== */
            magma_tally4_dgetmatrix( N, N, d_A, ldda, h_R, lda );
            
            blasf77_daxpy(&size, &c_neg_one, h_A, &ione, h_R, &ione);
            error = lapackf77_dlange("f", &N, &N, h_R, &lda, work);

            printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %s\n",
                   (int) N, cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000.,
                   (error == 0. ? "ok" : "failed") );
            status += ! (error == 0.);
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_R );
            
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
