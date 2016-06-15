/*
    -- MAGMA_tally3 (version 1.6.1) --
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
#include "testings.h"  // before magma_tally3.h, to include cublas_v2
#include "flops.h"
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"

#define PRECISION_z

int main(int argc, char **argv)
{
    TESTING_INIT();

    real_Double_t   gflops, magma_tally3_perf, magma_tally3_time, dev_perf, dev_time, cpu_perf, cpu_time;
    double          magma_tally3_error, dev_error, work[1];
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t M, N, Xm, Ym, lda, sizeA, sizeX, sizeY;
    magma_tally3_int_t incx = 1;
    magma_tally3_int_t incy = 1;
    magma_tally3DoubleComplex c_neg_one = MAGMA_tally3_Z_NEG_ONE;
    magma_tally3DoubleComplex alpha = MAGMA_tally3_Z_MAKE(  1.5, -2.3 );
    magma_tally3DoubleComplex beta  = MAGMA_tally3_Z_MAKE( -0.6,  0.8 );
    magma_tally3DoubleComplex *A, *X, *Y, *Ydev, *Ymagma_tally3;
    magma_tally3DoubleComplex_ptr dA, dX, dY;
    magma_tally3_int_t status = 0;
    
    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    
    double tol = opts.tolerance * lapackf77_dlamch("E");

    printf("trans = %s\n", lapack_trans_const_tally3(opts.transA) );
    #ifdef HAVE_CUBLAS
        printf("    M     N   MAGMA_tally3 Gflop/s (ms)  %s Gflop/s (ms)   CPU Gflop/s (ms)  MAGMA_tally3 error  %s error\n",
                g_platform_str, g_platform_str );
    #else
        printf("    M     N   %s Gflop/s (ms)   CPU Gflop/s (ms)  %s error\n",
                g_platform_str, g_platform_str );
    #endif
    printf("===================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            lda    = ((M+31)/32)*32;
            gflops = FLOPS_ZGEMV( M, N ) / 1e9;

            if ( opts.transA == Magma_tally3NoTrans ) {
                Xm = N;
                Ym = M;
            } else {
                Xm = M;
                Ym = N;
            }

            sizeA = lda*N;
            sizeX = incx*Xm;
            sizeY = incy*Ym;
            
            TESTING_MALLOC_CPU( A,       magma_tally3DoubleComplex, sizeA );
            TESTING_MALLOC_CPU( X,       magma_tally3DoubleComplex, sizeX );
            TESTING_MALLOC_CPU( Y,       magma_tally3DoubleComplex, sizeY );
            TESTING_MALLOC_CPU( Ydev,    magma_tally3DoubleComplex, sizeY );
            TESTING_MALLOC_CPU( Ymagma_tally3,  magma_tally3DoubleComplex, sizeY );
            
            TESTING_MALLOC_DEV( dA, magma_tally3DoubleComplex, sizeA );
            TESTING_MALLOC_DEV( dX, magma_tally3DoubleComplex, sizeX );
            TESTING_MALLOC_DEV( dY, magma_tally3DoubleComplex, sizeY );
            
            /* Initialize the matrix */
            lapackf77_zlarnv( &ione, ISEED, &sizeA, A );
            lapackf77_zlarnv( &ione, ISEED, &sizeX, X );
            lapackf77_zlarnv( &ione, ISEED, &sizeY, Y );
            
            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            magma_tally3_zsetmatrix( M, N, A, lda, dA, lda );
            magma_tally3_zsetvector( Xm, X, incx, dX, incx );
            magma_tally3_zsetvector( Ym, Y, incy, dY, incy );
            
            #ifdef HAVE_CUBLAS
                dev_time = magma_tally3_sync_wtime( 0 );
                cublasZgemv( opts.handle, cublas_trans_const_tally3(opts.transA),
                             M, N, &alpha, dA, lda, dX, incx, &beta, dY, incy );
                dev_time = magma_tally3_sync_wtime( 0 ) - dev_time;
            #else
                dev_time = magma_tally3_sync_wtime( opts.queue );
                magma_tally3_zgemv( opts.transA, M, N,
                             &alpha, dA, lda,
                                     dX, incx,
                             &beta,  dY, incy );
                dev_time = magma_tally3_sync_wtime( opts.queue ) - dev_time;
            #endif
            dev_perf = gflops / dev_time;
            
            magma_tally3_zgetvector( Ym, dY, incy, Ydev, incy );
            
            /* =====================================================================
               Performs operation using MAGMA_tally3BLAS (currently only with CUDA)
               =================================================================== */
            #ifdef HAVE_CUBLAS
                magma_tally3_zsetvector( Ym, Y, incy, dY, incy );
                
                magma_tally3_time = magma_tally3_sync_wtime( 0 );
                magma_tally3blas_zgemv( opts.transA, M, N, alpha, dA, lda, dX, incx, beta, dY, incy );
                magma_tally3_time = magma_tally3_sync_wtime( 0 ) - magma_tally3_time;
                magma_tally3_perf = gflops / magma_tally3_time;
                
                magma_tally3_zgetvector( Ym, dY, incy, Ymagma_tally3, incy );
            #endif
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            cpu_time = magma_tally3_wtime();
            blasf77_zgemv( lapack_trans_const_tally3(opts.transA), &M, &N,
                           &alpha, A, &lda,
                                   X, &incx,
                           &beta,  Y, &incy );
            cpu_time = magma_tally3_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            
            /* =====================================================================
               Check the result
               =================================================================== */
            double Anorm = lapackf77_zlange( "F", &M, &N, A, &lda, work );
            double Xnorm = lapackf77_zlange( "F", &Xm, &ione, X, &Xm, work );
            
            blasf77_zaxpy( &Ym, &c_neg_one, Y, &incy, Ydev, &incy );
            dev_error = lapackf77_zlange( "F", &Ym, &ione, Ydev, &Ym, work ) / (Anorm * Xnorm);
            
            #ifdef HAVE_CUBLAS
                blasf77_zaxpy( &Ym, &c_neg_one, Y, &incy, Ymagma_tally3, &incy );
                magma_tally3_error = lapackf77_zlange( "F", &Ym, &ione, Ymagma_tally3, &Ym, work ) / (Anorm * Xnorm);
                
                printf("%5d %5d   %7.2f (%7.2f)    %7.2f (%7.2f)   %7.2f (%7.2f)    %8.2e     %8.2e   %s\n",
                       (int) M, (int) N,
                       magma_tally3_perf,  1000.*magma_tally3_time,
                       dev_perf,    1000.*dev_time,
                       cpu_perf,    1000.*cpu_time,
                       magma_tally3_error, dev_error,
                       (magma_tally3_error < tol && dev_error < tol ? "ok" : "failed"));
                status += ! (magma_tally3_error < tol && dev_error < tol);
            #else
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)    %8.2e   %s\n",
                       (int) M, (int) N,
                       dev_perf,    1000.*dev_time,
                       cpu_perf,    1000.*cpu_time,
                       dev_error,
                       (dev_error < tol ? "ok" : "failed"));
                status += ! (dev_error < tol);
            #endif
            
            TESTING_FREE_CPU( A );
            TESTING_FREE_CPU( X );
            TESTING_FREE_CPU( Y );
            TESTING_FREE_CPU( Ydev    );
            TESTING_FREE_CPU( Ymagma_tally3  );
            
            TESTING_FREE_DEV( dA );
            TESTING_FREE_DEV( dX );
            TESTING_FREE_DEV( dY );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    
    TESTING_FINALIZE();
    return status;
}
