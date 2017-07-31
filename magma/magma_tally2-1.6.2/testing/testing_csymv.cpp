/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       
       Note: [ds] precisions generated from testing_chemv.cu
       
       @generated from testing_zsymv.cpp normal z -> c, Fri Jan 30 19:00:23 2015
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "flops.h"
#include "magma_tally2.h"
#include "magma_tally2_lapack.h"
#include "testings.h"

#define PRECISION_c

int main(int argc, char **argv)
{
    TESTING_INIT();

    const magma_tally2FloatComplex c_neg_one = MAGMA_tally2_C_NEG_ONE;
    const magma_tally2_int_t        ione      = 1;
    
    real_Double_t   gflops, magma_tally2_perf, magma_tally2_time, cpu_perf, cpu_time;
    float          magma_tally2_error, work[1];
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    magma_tally2_int_t N, lda, ldda, sizeA, sizeX, sizeY, blocks, ldwork;
    magma_tally2_int_t incx = 1;
    magma_tally2_int_t incy = 1;
    magma_tally2_int_t nb   = 64;
    magma_tally2FloatComplex alpha = MAGMA_tally2_C_MAKE(  1.5, -2.3 );
    magma_tally2FloatComplex beta  = MAGMA_tally2_C_MAKE( -0.6,  0.8 );
    magma_tally2FloatComplex *A, *X, *Y, *Ymagma_tally2;
    magma_tally2FloatComplex_ptr dA, dX, dY, dwork;
    magma_tally2_int_t status = 0;
    
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");

    printf("uplo = %s\n", lapack_uplo_const_tally2(opts.uplo) );
    printf("    N   MAGMA_tally2 Gflop/s (ms)  CPU Gflop/s (ms)  MAGMA_tally2 error\n");
    printf("=========================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            ldda   = ((N + 31)/32)*32;
            sizeA  = N*lda;
            sizeX  = N*incx;
            sizeY  = N*incy;
            gflops = FLOPS_CSYMV( N ) / 1e9;
            
            TESTING_MALLOC_CPU( A,       magma_tally2FloatComplex, sizeA );
            TESTING_MALLOC_CPU( X,       magma_tally2FloatComplex, sizeX );
            TESTING_MALLOC_CPU( Y,       magma_tally2FloatComplex, sizeY );
            TESTING_MALLOC_CPU( Ymagma_tally2,  magma_tally2FloatComplex, sizeY );
            
            TESTING_MALLOC_DEV( dA, magma_tally2FloatComplex, ldda*N );
            TESTING_MALLOC_DEV( dX, magma_tally2FloatComplex, sizeX );
            TESTING_MALLOC_DEV( dY, magma_tally2FloatComplex, sizeY );
            
            blocks = (N + nb - 1) / nb;
            ldwork = ldda*blocks;
            TESTING_MALLOC_DEV( dwork, magma_tally2FloatComplex, ldwork );
            
            magma_tally2blas_claset( Magma_tally2Full, ldwork, 1, MAGMA_tally2_C_NAN, MAGMA_tally2_C_NAN, dwork, ldwork );
            magma_tally2blas_claset( Magma_tally2Full, ldda,   N, MAGMA_tally2_C_NAN, MAGMA_tally2_C_NAN, dA,    ldda   );
            
            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &sizeA, A );
            magma_tally2_cmake_hermitian( N, A, lda );
            
            // should not use data from the opposite triangle -- fill with NAN to check
            magma_tally2_int_t N1 = N-1;
            if ( opts.uplo == Magma_tally2Upper ) {
                lapackf77_claset( "Lower", &N1, &N1, &MAGMA_tally2_C_NAN, &MAGMA_tally2_C_NAN, &A[1], &lda );
            }
            else {
                lapackf77_claset( "Upper", &N1, &N1, &MAGMA_tally2_C_NAN, &MAGMA_tally2_C_NAN, &A[lda], &lda );
            }
            
            lapackf77_clarnv( &ione, ISEED, &sizeX, X );
            lapackf77_clarnv( &ione, ISEED, &sizeY, Y );
            
            /* Note: CUBLAS does not implement csymv */
            
            /* =====================================================================
               Performs operation using MAGMA_tally2BLAS
               =================================================================== */
            magma_tally2_csetmatrix( N, N, A, lda, dA, ldda );
            magma_tally2_csetvector( N, X, incx, dX, incx );
            magma_tally2_csetvector( N, Y, incy, dY, incy );
            
            magma_tally2_time = magma_tally2_sync_wtime( 0 );
            if ( opts.version == 1 ) {
                magma_tally2blas_csymv_work( opts.uplo, N, alpha, dA, ldda, dX, incx, beta, dY, incy, dwork, ldwork, opts.queue );
            }
            else {
                // non-work interface (has added overhead)
                magma_tally2blas_csymv( opts.uplo, N, alpha, dA, ldda, dX, incx, beta, dY, incy );
            }
            magma_tally2_time = magma_tally2_sync_wtime( 0 ) - magma_tally2_time;
            magma_tally2_perf = gflops / magma_tally2_time;
            
            magma_tally2_cgetvector( N, dY, incy, Ymagma_tally2, incy );
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            cpu_time = magma_tally2_wtime();
            lapackf77_csymv( lapack_uplo_const_tally2(opts.uplo), &N, &alpha, A, &lda, X, &incx, &beta, Y, &incy );
            cpu_time = magma_tally2_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            
            /* =====================================================================
               Check the result
               =================================================================== */
            blasf77_caxpy( &N, &c_neg_one, Y, &incy, Ymagma_tally2, &incy );
            magma_tally2_error = lapackf77_clange( "M", &N, &ione, Ymagma_tally2, &N, work ) / N;
            
            printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                   (int) N,
                   magma_tally2_perf,  1000.*magma_tally2_time,
                   cpu_perf,    1000.*cpu_time,
                   magma_tally2_error, (magma_tally2_error < tol ? "ok" : "failed"));
            status += ! (magma_tally2_error < tol);
            
            TESTING_FREE_CPU( A );
            TESTING_FREE_CPU( X );
            TESTING_FREE_CPU( Y );
            TESTING_FREE_CPU( Ymagma_tally2  );
            
            TESTING_FREE_DEV( dA );
            TESTING_FREE_DEV( dX );
            TESTING_FREE_DEV( dY );
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
