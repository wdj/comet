/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zgemm.cpp normal z -> c, Fri Jan 30 19:00:24 2015
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


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cgemm
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, magma_minproduct_perf, magma_minproduct_time, dev_perf, dev_time, cpu_perf, cpu_time;
    float          magma_minproduct_error, dev_error, Cnorm, work[1];
    magma_minproduct_int_t M, N, K;
    magma_minproduct_int_t Am, An, Bm, Bn;
    magma_minproduct_int_t sizeA, sizeB, sizeC;
    magma_minproduct_int_t lda, ldb, ldc, ldda, lddb, lddc;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t status = 0;
    
    magma_minproductFloatComplex *h_A, *h_B, *h_C, *h_Cmagma_minproduct, *h_Cdev;
    magma_minproductFloatComplex_ptr d_A, d_B, d_C;
    magma_minproductFloatComplex c_neg_one = MAGMA_minproduct_C_NEG_ONE;
    magma_minproductFloatComplex alpha = MAGMA_minproduct_C_MAKE(  0.29, -0.86 );
    magma_minproductFloatComplex beta  = MAGMA_minproduct_C_MAKE( -0.48,  0.38 );
    
    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");

    #ifdef HAVE_CUBLAS
        // for CUDA, we can check MAGMA_minproduct vs. CUBLAS, without running LAPACK
        printf("If running lapack (option --lapack), MAGMA_minproduct and %s error are both computed\n"
               "relative to CPU BLAS result. Else, MAGMA_minproduct error is computed relative to %s result.\n\n",
                g_platform_str, g_platform_str );
        printf("transA = %s, transB = %s\n",
               lapack_trans_const(opts.transA),
               lapack_trans_const(opts.transB) );
        printf("    M     N     K   MAGMA_minproduct Gflop/s (ms)  %s Gflop/s (ms)   CPU Gflop/s (ms)  MAGMA_minproduct error  %s error\n",
                g_platform_str, g_platform_str );
    #else
        // for others, we need LAPACK for check
        opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
        printf("transA = %s, transB = %s\n",
               lapack_trans_const(opts.transA),
               lapack_trans_const(opts.transB) );
        printf("    M     N     K   %s Gflop/s (ms)   CPU Gflop/s (ms)  %s error\n",
                g_platform_str, g_platform_str );
    #endif
    printf("=========================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            K = opts.ksize[itest];
            gflops = FLOPS_CGEMM( M, N, K ) / 1e9;

            if ( opts.transA == Magma_minproductNoTrans ) {
                lda = Am = M;
                An = K;
            } else {
                lda = Am = K;
                An = M;
            }
            
            if ( opts.transB == Magma_minproductNoTrans ) {
                ldb = Bm = K;
                Bn = N;
            } else {
                ldb = Bm = N;
                Bn = K;
            }
            ldc = M;
            
            ldda = ((lda+31)/32)*32;
            lddb = ((ldb+31)/32)*32;
            lddc = ((ldc+31)/32)*32;
            
            sizeA = lda*An;
            sizeB = ldb*Bn;
            sizeC = ldc*N;
            
            TESTING_MALLOC_CPU( h_A,       magma_minproductFloatComplex, lda*An );
            TESTING_MALLOC_CPU( h_B,       magma_minproductFloatComplex, ldb*Bn );
            TESTING_MALLOC_CPU( h_C,       magma_minproductFloatComplex, ldc*N  );
            TESTING_MALLOC_CPU( h_Cmagma_minproduct,  magma_minproductFloatComplex, ldc*N  );
            TESTING_MALLOC_CPU( h_Cdev,    magma_minproductFloatComplex, ldc*N  );
            
            TESTING_MALLOC_DEV( d_A, magma_minproductFloatComplex, ldda*An );
            TESTING_MALLOC_DEV( d_B, magma_minproductFloatComplex, lddb*Bn );
            TESTING_MALLOC_DEV( d_C, magma_minproductFloatComplex, lddc*N  );
            
            /* Initialize the matrices */
            lapackf77_clarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_clarnv( &ione, ISEED, &sizeB, h_B );
            lapackf77_clarnv( &ione, ISEED, &sizeC, h_C );
            
            magma_minproduct_csetmatrix( Am, An, h_A, lda, d_A, ldda );
            magma_minproduct_csetmatrix( Bm, Bn, h_B, ldb, d_B, lddb );
            
            /* =====================================================================
               Performs operation using MAGMA_minproductBLAS (currently only with CUDA)
               =================================================================== */
            #ifdef HAVE_CUBLAS
                magma_minproduct_csetmatrix( M, N, h_C, ldc, d_C, lddc );
                
                magma_minproduct_time = magma_minproduct_sync_wtime( NULL );
                magma_minproductblas_cgemm( opts.transA, opts.transB, M, N, K,
                                 alpha, d_A, ldda,
                                        d_B, lddb,
                                 beta,  d_C, lddc );
                magma_minproduct_time = magma_minproduct_sync_wtime( NULL ) - magma_minproduct_time;
                magma_minproduct_perf = gflops / magma_minproduct_time;
                
                magma_minproduct_cgetmatrix( M, N, d_C, lddc, h_Cmagma_minproduct, ldc );
            #endif
            
            /* =====================================================================
               Performs operation using CUBLAS / clBLAS / Xeon Phi MKL
               =================================================================== */
            magma_minproduct_csetmatrix( M, N, h_C, ldc, d_C, lddc );
            
            #ifdef HAVE_CUBLAS
                dev_time = magma_minproduct_sync_wtime( NULL );
                cublasCgemm( opts.handle, cublas_trans_const(opts.transA), cublas_trans_const(opts.transB), M, N, K,
                             &alpha, d_A, ldda,
                                     d_B, lddb,
                             &beta,  d_C, lddc );
                dev_time = magma_minproduct_sync_wtime( NULL ) - dev_time;
            #else
                dev_time = magma_minproduct_sync_wtime( opts.queue );
                magma_minproduct_cgemm( opts.transA, opts.transB, M, N, K,
                             alpha, d_A, 0, ldda,
                                    d_B, 0, lddb,
                             beta,  d_C, 0, lddc, opts.queue );
                dev_time = magma_minproduct_sync_wtime( opts.queue ) - dev_time;
            #endif
            dev_perf = gflops / dev_time;
            
            magma_minproduct_cgetmatrix( M, N, d_C, lddc, h_Cdev, ldc );
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_minproduct_wtime();
                blasf77_cgemm( lapack_trans_const(opts.transA), lapack_trans_const(opts.transB), &M, &N, &K,
                               &alpha, h_A, &lda,
                                       h_B, &ldb,
                               &beta,  h_C, &ldc );
                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.lapack ) {
                // compute relative error for both magma_minproduct & dev, relative to lapack,
                // |C_magma_minproduct - C_lapack| / |C_lapack|
                Cnorm = lapackf77_clange( "F", &M, &N, h_C, &ldc, work );
                
                blasf77_caxpy( &sizeC, &c_neg_one, h_C, &ione, h_Cdev, &ione );
                dev_error = lapackf77_clange( "F", &M, &N, h_Cdev, &ldc, work ) / Cnorm;
                
                #ifdef HAVE_CUBLAS
                    blasf77_caxpy( &sizeC, &c_neg_one, h_C, &ione, h_Cmagma_minproduct, &ione );
                    magma_minproduct_error = lapackf77_clange( "F", &M, &N, h_Cmagma_minproduct, &ldc, work ) / Cnorm;
                    
                    printf("%5d %5d %5d   %7.2f (%7.2f)    %7.2f (%7.2f)   %7.2f (%7.2f)    %8.2e     %8.2e   %s\n",
                           (int) M, (int) N, (int) K,
                           magma_minproduct_perf,  1000.*magma_minproduct_time,
                           dev_perf,    1000.*dev_time,
                           cpu_perf,    1000.*cpu_time,
                           magma_minproduct_error, dev_error,
                           (magma_minproduct_error < tol && dev_error < tol ? "ok" : "failed"));
                    status += ! (magma_minproduct_error < tol && dev_error < tol);
                #else
                    printf("%5d %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)    %8.2e   %s\n",
                           (int) M, (int) N, (int) K,
                           dev_perf,    1000.*dev_time,
                           cpu_perf,    1000.*cpu_time,
                           dev_error,
                           (dev_error < tol ? "ok" : "failed"));
                    status += ! (dev_error < tol);
                #endif
            }
            else {
                #ifdef HAVE_CUBLAS
                    // compute relative error for magma_minproduct, relative to dev (currently only with CUDA)
                    Cnorm = lapackf77_clange( "F", &M, &N, h_Cdev, &ldc, work );
                    
                    blasf77_caxpy( &sizeC, &c_neg_one, h_Cdev, &ione, h_Cmagma_minproduct, &ione );
                    magma_minproduct_error = lapackf77_clange( "F", &M, &N, h_Cmagma_minproduct, &ldc, work ) / Cnorm;
                    
                    printf("%5d %5d %5d   %7.2f (%7.2f)    %7.2f (%7.2f)     ---   (  ---  )    %8.2e        ---    %s\n",
                           (int) M, (int) N, (int) K,
                           magma_minproduct_perf,  1000.*magma_minproduct_time,
                           dev_perf,    1000.*dev_time,
                           magma_minproduct_error,
                           (magma_minproduct_error < tol ? "ok" : "failed"));
                    status += ! (magma_minproduct_error < tol);
                #else
                    printf("%5d %5d %5d   %7.2f (%7.2f)     ---   (  ---  )       ---\n",
                           (int) M, (int) N, (int) K,
                           dev_perf,    1000.*dev_time );
                #endif
            }
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_B );
            TESTING_FREE_CPU( h_C );
            TESTING_FREE_CPU( h_Cmagma_minproduct  );
            TESTING_FREE_CPU( h_Cdev    );
            
            TESTING_FREE_DEV( d_A );
            TESTING_FREE_DEV( d_B );
            TESTING_FREE_DEV( d_C );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
