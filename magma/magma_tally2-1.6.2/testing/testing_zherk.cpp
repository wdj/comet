/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> c d s
       @author Chongxiao Cao
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "testings.h"  // before magma_tally2.h, to include cublas_v2
#include "flops.h"
#include "magma_tally2.h"
#include "magma_tally2_lapack.h"

#include "magma_tally2_threadsetting.h"  // to work around MKL bug


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zherk
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, cublas_perf, cublas_time, cpu_perf, cpu_time;
    double          cublas_error, Cnorm, work[1];
    magma_tally2_int_t N, K;
    magma_tally2_int_t Ak, An;
    magma_tally2_int_t sizeA, sizeC;
    magma_tally2_int_t lda, ldc, ldda, lddc;
    magma_tally2_int_t ione     = 1;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    
    magma_tally2DoubleComplex *h_A, *h_C, *h_Ccublas;
    magma_tally2DoubleComplex_ptr d_A, d_C;
    magma_tally2DoubleComplex c_neg_one = MAGMA_tally2_Z_NEG_ONE;
    double alpha = MAGMA_tally2_D_MAKE(  0.29, -0.86 );
    double beta  = MAGMA_tally2_D_MAKE( -0.48,  0.38 );
    magma_tally2_int_t status = 0;
    
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("If running lapack (option --lapack), CUBLAS error is computed\n"
           "relative to CPU BLAS result.\n\n");
    printf("uplo = %s, transA = %s\n",
           lapack_uplo_const_tally2(opts.uplo), lapack_trans_const_tally2(opts.transA) );
    printf("    N     K   CUBLAS Gflop/s (ms)   CPU Gflop/s (ms)  CUBLAS error\n");
    printf("==================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            K = opts.ksize[itest];
            gflops = FLOPS_ZHERK(K, N) / 1e9;

            if ( opts.transA == Magma_tally2NoTrans ) {
                lda = An = N;
                Ak = K;
            } else {
                lda = An = K;
                Ak = N;
            }
            
            ldc = N;
            
            ldda = ((lda+31)/32)*32;
            lddc = ((ldc+31)/32)*32;
            
            sizeA = lda*Ak;
            sizeC = ldc*N;
            
            TESTING_MALLOC_CPU( h_A,       magma_tally2DoubleComplex, lda*Ak );
            TESTING_MALLOC_CPU( h_C,       magma_tally2DoubleComplex, ldc*N  );
            TESTING_MALLOC_CPU( h_Ccublas, magma_tally2DoubleComplex, ldc*N  );
            
            TESTING_MALLOC_DEV( d_A, magma_tally2DoubleComplex, ldda*Ak );
            TESTING_MALLOC_DEV( d_C, magma_tally2DoubleComplex, lddc*N  );
            
            /* Initialize the matrices */
            lapackf77_zlarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_zlarnv( &ione, ISEED, &sizeC, h_C );
            
            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            magma_tally2_zsetmatrix( An, Ak, h_A, lda, d_A, ldda );
            magma_tally2_zsetmatrix( N, N, h_C, ldc, d_C, lddc );

            cublas_time = magma_tally2_sync_wtime( NULL );
            cublasZherk( opts.handle, cublas_uplo_const_tally2(opts.uplo), cublas_trans_const_tally2(opts.transA), N, K,
                         &alpha, d_A, ldda,
                         &beta,  d_C, lddc );
            cublas_time = magma_tally2_sync_wtime( NULL ) - cublas_time;
            cublas_perf = gflops / cublas_time;
            
            magma_tally2_zgetmatrix( N, N, d_C, lddc, h_Ccublas, ldc );
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally2_wtime();
                blasf77_zherk( lapack_uplo_const_tally2(opts.uplo), lapack_trans_const_tally2(opts.transA), &N, &K,
                               &alpha, h_A, &lda,
                               &beta,  h_C, &ldc );
                cpu_time = magma_tally2_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.lapack ) {
                #ifdef MAGMA_tally2_WITH_MKL
                // MKL (11.1.2) has bug in multi-threaded zlanhe; use single thread to work around
                int threads = magma_tally2_get_lapack_numthreads();
                magma_tally2_set_lapack_numthreads( 1 );
                #endif
                
                // compute relative error for both magma_tally2 & cublas, relative to lapack,
                // |C_magma_tally2 - C_lapack| / |C_lapack|
                Cnorm = lapackf77_zlanhe("fro", lapack_uplo_const_tally2(opts.uplo), &N, h_C, &ldc, work);

                blasf77_zaxpy( &sizeC, &c_neg_one, h_C, &ione, h_Ccublas, &ione );
                cublas_error = lapackf77_zlanhe( "fro", lapack_uplo_const_tally2(opts.uplo), &N, h_Ccublas, &ldc, work ) / Cnorm;
                
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)    %8.2e   %s\n",
                       (int) N, (int) K,
                       cublas_perf, 1000.*cublas_time,
                       cpu_perf,    1000.*cpu_time,
                       cublas_error, (cublas_error < tol ? "ok" : "failed"));
                status += ! (cublas_error < tol);
                
                #ifdef MAGMA_tally2_WITH_MKL
                // end single thread to work around MKL bug
                magma_tally2_set_lapack_numthreads( threads );
                #endif
            }
            else {
                printf("%5d %5d   %7.2f (%7.2f)    ---   (  ---  )    ---     ---\n",
                       (int) N, (int) K,
                       cublas_perf, 1000.*cublas_time);
            }
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_C );
            TESTING_FREE_CPU( h_Ccublas );
            
            TESTING_FREE_DEV( d_A );
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
