/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zhemv.cpp normal z -> s, Fri Jan 30 19:00:23 2015
       
       @author Mark Gates
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "testings.h"  // before magma_minproduct.h, to include cublas_v2
#include "flops.h"
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"

#define PRECISION_s

int main(int argc, char **argv)
{
    TESTING_INIT();

    const float c_neg_one = MAGMA_minproduct_S_NEG_ONE;
    const magma_minproduct_int_t        ione      = 1;
    
    real_Double_t   atomics_perf, atomics_time;
    real_Double_t   gflops, magma_minproduct_perf, magma_minproduct_time, cublas_perf, cublas_time, cpu_perf, cpu_time;
    float          magma_minproduct_error, atomics_error, cublas_error, work[1];
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t N, lda, ldda, sizeA, sizeX, sizeY, blocks, ldwork;
    magma_minproduct_int_t incx = 1;
    magma_minproduct_int_t incy = 1;
    magma_minproduct_int_t nb   = 64;
    float alpha = MAGMA_minproduct_S_MAKE(  1.5, -2.3 );
    float beta  = MAGMA_minproduct_S_MAKE( -0.6,  0.8 );
    float *A, *X, *Y, *Yatomics, *Ycublas, *Ymagma_minproduct;
    magma_minproductFloat_ptr dA, dX, dY, dwork;
    magma_minproduct_int_t status = 0;
    
    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");

    printf("uplo = %s\n", lapack_uplo_const(opts.uplo) );
    printf("    N   MAGMA_minproduct Gflop/s (ms)    Atomics Gflop/s      CUBLAS Gflop/s       CPU Gflop/s   MAGMA_minproduct error  Atomics    CUBLAS\n");
    printf("======================================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            ldda   = ((N + 31)/32)*32;
            sizeA  = N*lda;
            sizeX  = N*incx;
            sizeY  = N*incy;
            gflops = FLOPS_SSYMV( N ) / 1e9;
            
            TESTING_MALLOC_CPU( A,        float, sizeA );
            TESTING_MALLOC_CPU( X,        float, sizeX );
            TESTING_MALLOC_CPU( Y,        float, sizeY );
            TESTING_MALLOC_CPU( Yatomics, float, sizeY );
            TESTING_MALLOC_CPU( Ycublas,  float, sizeY );
            TESTING_MALLOC_CPU( Ymagma_minproduct,   float, sizeY );
            
            TESTING_MALLOC_DEV( dA, float, ldda*N );
            TESTING_MALLOC_DEV( dX, float, sizeX );
            TESTING_MALLOC_DEV( dY, float, sizeY );
            
            blocks = (N + nb - 1) / nb;
            ldwork = ldda*blocks;
            TESTING_MALLOC_DEV( dwork, float, ldwork );
            
            magma_minproductblas_slaset( Magma_minproductFull, ldwork, 1, MAGMA_minproduct_S_NAN, MAGMA_minproduct_S_NAN, dwork, ldwork );
            magma_minproductblas_slaset( Magma_minproductFull, ldda,   N, MAGMA_minproduct_S_NAN, MAGMA_minproduct_S_NAN, dA,    ldda   );
            
            /* Initialize the matrix */
            lapackf77_slarnv( &ione, ISEED, &sizeA, A );
            magma_minproduct_smake_symmetric( N, A, lda );
            
            // should not use data from the opposite triangle -- fill with NAN to check
            magma_minproduct_int_t N1 = N-1;
            if ( opts.uplo == Magma_minproductUpper ) {
                lapackf77_slaset( "Lower", &N1, &N1, &MAGMA_minproduct_S_NAN, &MAGMA_minproduct_S_NAN, &A[1], &lda );
            }
            else {
                lapackf77_slaset( "Upper", &N1, &N1, &MAGMA_minproduct_S_NAN, &MAGMA_minproduct_S_NAN, &A[lda], &lda );
            }
            
            lapackf77_slarnv( &ione, ISEED, &sizeX, X );
            lapackf77_slarnv( &ione, ISEED, &sizeY, Y );
            
            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            magma_minproduct_ssetmatrix( N, N, A, lda, dA, ldda );
            magma_minproduct_ssetvector( N, X, incx, dX, incx );
            magma_minproduct_ssetvector( N, Y, incy, dY, incy );
            
            cublas_time = magma_minproduct_sync_wtime( 0 );
            cublasSsymv( opts.handle, cublas_uplo_const(opts.uplo),
                         N, &alpha, dA, ldda, dX, incx, &beta, dY, incy );
            cublas_time = magma_minproduct_sync_wtime( 0 ) - cublas_time;
            cublas_perf = gflops / cublas_time;
            
            magma_minproduct_sgetvector( N, dY, incy, Ycublas, incy );
            
            /* =====================================================================
               Performs operation using CUBLAS - using atomics
               =================================================================== */
            cublasSetAtomicsMode( opts.handle, CUBLAS_ATOMICS_ALLOWED );
            magma_minproduct_ssetvector( N, Y, incy, dY, incy );
            
            atomics_time = magma_minproduct_sync_wtime( 0 );
            cublasSsymv( opts.handle, cublas_uplo_const(opts.uplo),
                         N, &alpha, dA, ldda, dX, incx, &beta, dY, incy );
            atomics_time = magma_minproduct_sync_wtime( 0 ) - atomics_time;
            atomics_perf = gflops / atomics_time;
            
            magma_minproduct_sgetvector( N, dY, incy, Yatomics, incy );
            cublasSetAtomicsMode( opts.handle, CUBLAS_ATOMICS_NOT_ALLOWED );
            
            /* =====================================================================
               Performs operation using MAGMA_minproductBLAS
               =================================================================== */
            magma_minproduct_ssetvector( N, Y, incy, dY, incy );
            
            magma_minproduct_time = magma_minproduct_sync_wtime( 0 );
            if ( opts.version == 1 ) {
                magma_minproductblas_ssymv_work( opts.uplo, N, alpha, dA, ldda, dX, incx, beta, dY, incy, dwork, ldwork, opts.queue );
            }
            else {
                // non-work interface (has added overhead)
                magma_minproductblas_ssymv( opts.uplo, N, alpha, dA, ldda, dX, incx, beta, dY, incy );
            }
            magma_minproduct_time = magma_minproduct_sync_wtime( 0 ) - magma_minproduct_time;
            magma_minproduct_perf = gflops / magma_minproduct_time;
            
            magma_minproduct_sgetvector( N, dY, incy, Ymagma_minproduct, incy );
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            cpu_time = magma_minproduct_wtime();
            blasf77_ssymv( lapack_uplo_const(opts.uplo), &N, &alpha, A, &lda, X, &incx, &beta, Y, &incy );
            cpu_time = magma_minproduct_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            
            /* =====================================================================
               Check the result
               =================================================================== */
            blasf77_saxpy( &N, &c_neg_one, Y, &incy, Ymagma_minproduct, &incy );
            magma_minproduct_error = lapackf77_slange( "M", &N, &ione, Ymagma_minproduct, &N, work ) / N;
            
            blasf77_saxpy( &N, &c_neg_one, Y, &incy, Ycublas, &incy );
            cublas_error = lapackf77_slange( "M", &N, &ione, Ycublas, &N, work ) / N;
            
            blasf77_saxpy( &N, &c_neg_one, Y, &incy, Yatomics, &incy );
            atomics_error = lapackf77_slange( "M", &N, &ione, Yatomics, &N, work ) / N;
            
            bool ok = (magma_minproduct_error < tol && cublas_error < tol && atomics_error < tol);
            status += ! ok;
            printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %8.2e   %8.2e   %s\n",
                   (int) N,
                   magma_minproduct_perf,   1000.*magma_minproduct_time,
                   atomics_perf, 1000.*atomics_time,
                   cublas_perf,  1000.*cublas_time,
                   cpu_perf,     1000.*cpu_time,
                   magma_minproduct_error, cublas_error, atomics_error,
                   (ok ? "ok" : "failed"));
            
            TESTING_FREE_CPU( A );
            TESTING_FREE_CPU( X );
            TESTING_FREE_CPU( Y );
            TESTING_FREE_CPU( Ycublas  );
            TESTING_FREE_CPU( Yatomics );
            TESTING_FREE_CPU( Ymagma_minproduct   );
            
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
