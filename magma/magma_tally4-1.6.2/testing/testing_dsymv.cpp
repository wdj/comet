/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zhemv.cpp normal z -> d, Fri Jan 30 19:00:23 2015
       
       @author Mark Gates
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "testings.h"  // before magma_tally4.h, to include cublas_v2
#include "flops.h"
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"

#define PRECISION_d

int main(int argc, char **argv)
{
    TESTING_INIT();

    const double c_neg_one = MAGMA_tally4_D_NEG_ONE;
    const magma_tally4_int_t        ione      = 1;
    
    real_Double_t   atomics_perf, atomics_time;
    real_Double_t   gflops, magma_tally4_perf, magma_tally4_time, cublas_perf, cublas_time, cpu_perf, cpu_time;
    double          magma_tally4_error, atomics_error, cublas_error, work[1];
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    magma_tally4_int_t N, lda, ldda, sizeA, sizeX, sizeY, blocks, ldwork;
    magma_tally4_int_t incx = 1;
    magma_tally4_int_t incy = 1;
    magma_tally4_int_t nb   = 64;
    double alpha = MAGMA_tally4_D_MAKE(  1.5, -2.3 );
    double beta  = MAGMA_tally4_D_MAKE( -0.6,  0.8 );
    double *A, *X, *Y, *Yatomics, *Ycublas, *Ymagma_tally4;
    magma_tally4Double_ptr dA, dX, dY, dwork;
    magma_tally4_int_t status = 0;
    
    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );
    
    double tol = opts.tolerance * lapackf77_dlamch("E");

    printf("uplo = %s\n", lapack_uplo_const_tally4(opts.uplo) );
    printf("    N   MAGMA_tally4 Gflop/s (ms)    Atomics Gflop/s      CUBLAS Gflop/s       CPU Gflop/s   MAGMA_tally4 error  Atomics    CUBLAS\n");
    printf("======================================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            ldda   = ((N + 31)/32)*32;
            sizeA  = N*lda;
            sizeX  = N*incx;
            sizeY  = N*incy;
            gflops = FLOPS_DSYMV( N ) / 1e9;
            
            TESTING_MALLOC_CPU( A,        double, sizeA );
            TESTING_MALLOC_CPU( X,        double, sizeX );
            TESTING_MALLOC_CPU( Y,        double, sizeY );
            TESTING_MALLOC_CPU( Yatomics, double, sizeY );
            TESTING_MALLOC_CPU( Ycublas,  double, sizeY );
            TESTING_MALLOC_CPU( Ymagma_tally4,   double, sizeY );
            
            TESTING_MALLOC_DEV( dA, double, ldda*N );
            TESTING_MALLOC_DEV( dX, double, sizeX );
            TESTING_MALLOC_DEV( dY, double, sizeY );
            
            blocks = (N + nb - 1) / nb;
            ldwork = ldda*blocks;
            TESTING_MALLOC_DEV( dwork, double, ldwork );
            
            magma_tally4blas_dlaset( Magma_tally4Full, ldwork, 1, MAGMA_tally4_D_NAN, MAGMA_tally4_D_NAN, dwork, ldwork );
            magma_tally4blas_dlaset( Magma_tally4Full, ldda,   N, MAGMA_tally4_D_NAN, MAGMA_tally4_D_NAN, dA,    ldda   );
            
            /* Initialize the matrix */
            lapackf77_dlarnv( &ione, ISEED, &sizeA, A );
            magma_tally4_dmake_symmetric( N, A, lda );
            
            // should not use data from the opposite triangle -- fill with NAN to check
            magma_tally4_int_t N1 = N-1;
            if ( opts.uplo == Magma_tally4Upper ) {
                lapackf77_dlaset( "Lower", &N1, &N1, &MAGMA_tally4_D_NAN, &MAGMA_tally4_D_NAN, &A[1], &lda );
            }
            else {
                lapackf77_dlaset( "Upper", &N1, &N1, &MAGMA_tally4_D_NAN, &MAGMA_tally4_D_NAN, &A[lda], &lda );
            }
            
            lapackf77_dlarnv( &ione, ISEED, &sizeX, X );
            lapackf77_dlarnv( &ione, ISEED, &sizeY, Y );
            
            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            magma_tally4_dsetmatrix( N, N, A, lda, dA, ldda );
            magma_tally4_dsetvector( N, X, incx, dX, incx );
            magma_tally4_dsetvector( N, Y, incy, dY, incy );
            
            cublas_time = magma_tally4_sync_wtime( 0 );
            cublasDsymv( opts.handle, cublas_uplo_const_tally4(opts.uplo),
                         N, &alpha, dA, ldda, dX, incx, &beta, dY, incy );
            cublas_time = magma_tally4_sync_wtime( 0 ) - cublas_time;
            cublas_perf = gflops / cublas_time;
            
            magma_tally4_dgetvector( N, dY, incy, Ycublas, incy );
            
            /* =====================================================================
               Performs operation using CUBLAS - using atomics
               =================================================================== */
            cublasSetAtomicsMode( opts.handle, CUBLAS_ATOMICS_ALLOWED );
            magma_tally4_dsetvector( N, Y, incy, dY, incy );
            
            atomics_time = magma_tally4_sync_wtime( 0 );
            cublasDsymv( opts.handle, cublas_uplo_const_tally4(opts.uplo),
                         N, &alpha, dA, ldda, dX, incx, &beta, dY, incy );
            atomics_time = magma_tally4_sync_wtime( 0 ) - atomics_time;
            atomics_perf = gflops / atomics_time;
            
            magma_tally4_dgetvector( N, dY, incy, Yatomics, incy );
            cublasSetAtomicsMode( opts.handle, CUBLAS_ATOMICS_NOT_ALLOWED );
            
            /* =====================================================================
               Performs operation using MAGMA_tally4BLAS
               =================================================================== */
            magma_tally4_dsetvector( N, Y, incy, dY, incy );
            
            magma_tally4_time = magma_tally4_sync_wtime( 0 );
            if ( opts.version == 1 ) {
                magma_tally4blas_dsymv_work( opts.uplo, N, alpha, dA, ldda, dX, incx, beta, dY, incy, dwork, ldwork, opts.queue );
            }
            else {
                // non-work interface (has added overhead)
                magma_tally4blas_dsymv( opts.uplo, N, alpha, dA, ldda, dX, incx, beta, dY, incy );
            }
            magma_tally4_time = magma_tally4_sync_wtime( 0 ) - magma_tally4_time;
            magma_tally4_perf = gflops / magma_tally4_time;
            
            magma_tally4_dgetvector( N, dY, incy, Ymagma_tally4, incy );
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            cpu_time = magma_tally4_wtime();
            blasf77_dsymv( lapack_uplo_const_tally4(opts.uplo), &N, &alpha, A, &lda, X, &incx, &beta, Y, &incy );
            cpu_time = magma_tally4_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            
            /* =====================================================================
               Check the result
               =================================================================== */
            blasf77_daxpy( &N, &c_neg_one, Y, &incy, Ymagma_tally4, &incy );
            magma_tally4_error = lapackf77_dlange( "M", &N, &ione, Ymagma_tally4, &N, work ) / N;
            
            blasf77_daxpy( &N, &c_neg_one, Y, &incy, Ycublas, &incy );
            cublas_error = lapackf77_dlange( "M", &N, &ione, Ycublas, &N, work ) / N;
            
            blasf77_daxpy( &N, &c_neg_one, Y, &incy, Yatomics, &incy );
            atomics_error = lapackf77_dlange( "M", &N, &ione, Yatomics, &N, work ) / N;
            
            bool ok = (magma_tally4_error < tol && cublas_error < tol && atomics_error < tol);
            status += ! ok;
            printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %8.2e   %8.2e   %s\n",
                   (int) N,
                   magma_tally4_perf,   1000.*magma_tally4_time,
                   atomics_perf, 1000.*atomics_time,
                   cublas_perf,  1000.*cublas_time,
                   cpu_perf,     1000.*cpu_time,
                   magma_tally4_error, cublas_error, atomics_error,
                   (ok ? "ok" : "failed"));
            
            TESTING_FREE_CPU( A );
            TESTING_FREE_CPU( X );
            TESTING_FREE_CPU( Y );
            TESTING_FREE_CPU( Ycublas  );
            TESTING_FREE_CPU( Yatomics );
            TESTING_FREE_CPU( Ymagma_tally4   );
            
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
