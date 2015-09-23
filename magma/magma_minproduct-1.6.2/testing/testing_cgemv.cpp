/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zgemv.cpp normal z -> c, Fri Jan 30 19:00:24 2015
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

#define PRECISION_c

int main(int argc, char **argv)
{
    TESTING_INIT();

    real_Double_t   gflops, magma_minproduct_perf, magma_minproduct_time, dev_perf, dev_time, cpu_perf, cpu_time;
    float          magma_minproduct_error, dev_error, work[1];
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t M, N, Xm, Ym, lda, sizeA, sizeX, sizeY;
    magma_minproduct_int_t incx = 1;
    magma_minproduct_int_t incy = 1;
    magma_minproductFloatComplex c_neg_one = MAGMA_minproduct_C_NEG_ONE;
    magma_minproductFloatComplex alpha = MAGMA_minproduct_C_MAKE(  1.5, -2.3 );
    magma_minproductFloatComplex beta  = MAGMA_minproduct_C_MAKE( -0.6,  0.8 );
    magma_minproductFloatComplex *A, *X, *Y, *Ydev, *Ymagma_minproduct;
    magma_minproductFloatComplex_ptr dA, dX, dY;
    magma_minproduct_int_t status = 0;
    
    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");

    printf("trans = %s\n", lapack_trans_const(opts.transA) );
    #ifdef HAVE_CUBLAS
        printf("    M     N   MAGMA_minproduct Gflop/s (ms)  %s Gflop/s (ms)   CPU Gflop/s (ms)  MAGMA_minproduct error  %s error\n",
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
            gflops = FLOPS_CGEMV( M, N ) / 1e9;

            if ( opts.transA == Magma_minproductNoTrans ) {
                Xm = N;
                Ym = M;
            } else {
                Xm = M;
                Ym = N;
            }

            sizeA = lda*N;
            sizeX = incx*Xm;
            sizeY = incy*Ym;
            
            TESTING_MALLOC_CPU( A,       magma_minproductFloatComplex, sizeA );
            TESTING_MALLOC_CPU( X,       magma_minproductFloatComplex, sizeX );
            TESTING_MALLOC_CPU( Y,       magma_minproductFloatComplex, sizeY );
            TESTING_MALLOC_CPU( Ydev,    magma_minproductFloatComplex, sizeY );
            TESTING_MALLOC_CPU( Ymagma_minproduct,  magma_minproductFloatComplex, sizeY );
            
            TESTING_MALLOC_DEV( dA, magma_minproductFloatComplex, sizeA );
            TESTING_MALLOC_DEV( dX, magma_minproductFloatComplex, sizeX );
            TESTING_MALLOC_DEV( dY, magma_minproductFloatComplex, sizeY );
            
            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &sizeA, A );
            lapackf77_clarnv( &ione, ISEED, &sizeX, X );
            lapackf77_clarnv( &ione, ISEED, &sizeY, Y );
            
            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            magma_minproduct_csetmatrix( M, N, A, lda, dA, lda );
            magma_minproduct_csetvector( Xm, X, incx, dX, incx );
            magma_minproduct_csetvector( Ym, Y, incy, dY, incy );
            
            #ifdef HAVE_CUBLAS
                dev_time = magma_minproduct_sync_wtime( 0 );
                cublasCgemv( opts.handle, cublas_trans_const(opts.transA),
                             M, N, &alpha, dA, lda, dX, incx, &beta, dY, incy );
                dev_time = magma_minproduct_sync_wtime( 0 ) - dev_time;
            #else
                dev_time = magma_minproduct_sync_wtime( opts.queue );
                magma_minproduct_cgemv( opts.transA, M, N,
                             &alpha, dA, lda,
                                     dX, incx,
                             &beta,  dY, incy );
                dev_time = magma_minproduct_sync_wtime( opts.queue ) - dev_time;
            #endif
            dev_perf = gflops / dev_time;
            
            magma_minproduct_cgetvector( Ym, dY, incy, Ydev, incy );
            
            /* =====================================================================
               Performs operation using MAGMA_minproductBLAS (currently only with CUDA)
               =================================================================== */
            #ifdef HAVE_CUBLAS
                magma_minproduct_csetvector( Ym, Y, incy, dY, incy );
                
                magma_minproduct_time = magma_minproduct_sync_wtime( 0 );
                magma_minproductblas_cgemv( opts.transA, M, N, alpha, dA, lda, dX, incx, beta, dY, incy );
                magma_minproduct_time = magma_minproduct_sync_wtime( 0 ) - magma_minproduct_time;
                magma_minproduct_perf = gflops / magma_minproduct_time;
                
                magma_minproduct_cgetvector( Ym, dY, incy, Ymagma_minproduct, incy );
            #endif
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            cpu_time = magma_minproduct_wtime();
            blasf77_cgemv( lapack_trans_const(opts.transA), &M, &N,
                           &alpha, A, &lda,
                                   X, &incx,
                           &beta,  Y, &incy );
            cpu_time = magma_minproduct_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            
            /* =====================================================================
               Check the result
               =================================================================== */
            float Anorm = lapackf77_clange( "F", &M, &N, A, &lda, work );
            float Xnorm = lapackf77_clange( "F", &Xm, &ione, X, &Xm, work );
            
            blasf77_caxpy( &Ym, &c_neg_one, Y, &incy, Ydev, &incy );
            dev_error = lapackf77_clange( "F", &Ym, &ione, Ydev, &Ym, work ) / (Anorm * Xnorm);
            
            #ifdef HAVE_CUBLAS
                blasf77_caxpy( &Ym, &c_neg_one, Y, &incy, Ymagma_minproduct, &incy );
                magma_minproduct_error = lapackf77_clange( "F", &Ym, &ione, Ymagma_minproduct, &Ym, work ) / (Anorm * Xnorm);
                
                printf("%5d %5d   %7.2f (%7.2f)    %7.2f (%7.2f)   %7.2f (%7.2f)    %8.2e     %8.2e   %s\n",
                       (int) M, (int) N,
                       magma_minproduct_perf,  1000.*magma_minproduct_time,
                       dev_perf,    1000.*dev_time,
                       cpu_perf,    1000.*cpu_time,
                       magma_minproduct_error, dev_error,
                       (magma_minproduct_error < tol && dev_error < tol ? "ok" : "failed"));
                status += ! (magma_minproduct_error < tol && dev_error < tol);
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
            TESTING_FREE_CPU( Ymagma_minproduct  );
            
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
