/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zher2k.cpp normal z -> d, Fri Jan 30 19:00:23 2015
       @author Chongxiao Cao
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


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dsyr2k
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, cublas_perf, cublas_time, cpu_perf, cpu_time;
    double          cublas_error, Cnorm, work[1];
    magma_tally3_int_t N, K;
    magma_tally3_int_t Ak, An, Bk, Bn;
    magma_tally3_int_t sizeA, sizeB, sizeC;
    magma_tally3_int_t lda, ldb, ldc, ldda, lddb, lddc;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    
    double *h_A, *h_B, *h_C, *h_Ccublas;
    magma_tally3Double_ptr d_A, d_B, d_C;
    double c_neg_one = MAGMA_tally3_D_NEG_ONE;
    double alpha = MAGMA_tally3_D_MAKE(  0.29, -0.86 );
    double beta  = MAGMA_tally3_D_MAKE( -0.48,  0.38 );
    magma_tally3_int_t status = 0;
    
    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("If running lapack (option --lapack), CUBLAS error is computed\n"
           "relative to CPU BLAS result.\n\n");
    printf("uplo = %s, transA = %s\n",
           lapack_uplo_const_tally3(opts.uplo), lapack_trans_const_tally3(opts.transA) );
    printf("    N     K   CUBLAS Gflop/s (ms)   CPU Gflop/s (ms)  CUBLAS error\n");
    printf("==================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.msize[itest];
            K = opts.ksize[itest];
            gflops = FLOPS_DSYR2K(K, N) / 1e9;

            if ( opts.transA == Magma_tally3NoTrans ) {
                lda = An = N;
                Ak = K;
                ldb = Bn = N;
                Bk = K;
            } else {
                lda = An = K;
                Ak = N;
                ldb = Bn = K;
                Bk = N;
            }
            
            ldc = N;
            
            ldda = ((lda+31)/32)*32;
            lddb = ((ldb+31)/32)*32;
            lddc = ((ldc+31)/32)*32;
            
            sizeA = lda*Ak;
            sizeB = ldb*Ak;
            sizeC = ldc*N;
            
            TESTING_MALLOC_CPU( h_A,       double, lda*Ak );
            TESTING_MALLOC_CPU( h_B,       double, ldb*Bk );
            TESTING_MALLOC_CPU( h_C,       double, ldc*N  );
            TESTING_MALLOC_CPU( h_Ccublas, double, ldc*N  );
            
            TESTING_MALLOC_DEV( d_A, double, ldda*Ak );
            TESTING_MALLOC_DEV( d_B, double, lddb*Bk );
            TESTING_MALLOC_DEV( d_C, double, lddc*N  );
            
            /* Initialize the matrices */
            lapackf77_dlarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_dlarnv( &ione, ISEED, &sizeB, h_B );
            lapackf77_dlarnv( &ione, ISEED, &sizeC, h_C );
            
            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            magma_tally3_dsetmatrix( An, Ak, h_A, lda, d_A, ldda );
            magma_tally3_dsetmatrix( Bn, Bk, h_B, ldb, d_B, lddb );
            magma_tally3_dsetmatrix( N, N, h_C, ldc, d_C, lddc );
            
            cublas_time = magma_tally3_sync_wtime( NULL );
            cublasDsyr2k( opts.handle, cublas_uplo_const_tally3(opts.uplo), cublas_trans_const_tally3(opts.transA), N, K,
                          &alpha, d_A, ldda,
                                  d_B, lddb,
                          &beta,  d_C, lddc );
            cublas_time = magma_tally3_sync_wtime( NULL ) - cublas_time;
            cublas_perf = gflops / cublas_time;
            
            magma_tally3_dgetmatrix( N, N, d_C, lddc, h_Ccublas, ldc );
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally3_wtime();
                blasf77_dsyr2k( lapack_uplo_const_tally3(opts.uplo), lapack_trans_const_tally3(opts.transA), &N, &K,
                               &alpha, h_A, &lda,
                                       h_B, &ldb,
                               &beta,  h_C, &ldc );
                cpu_time = magma_tally3_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.lapack ) {
                // compute relative error for both magma_tally3 & cublas, relative to lapack,
                // |C_magma_tally3 - C_lapack| / |C_lapack|
                Cnorm = lapackf77_dlange( "M", &N, &N, h_C, &ldc, work );
                
                blasf77_daxpy( &sizeC, &c_neg_one, h_C, &ione, h_Ccublas, &ione );
                cublas_error = lapackf77_dlange( "M", &N, &N, h_Ccublas, &ldc, work ) / Cnorm;
                
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)    %8.2e   %s\n",
                       (int) N, (int) K,
                       cublas_perf, 1000.*cublas_time,
                       cpu_perf,    1000.*cpu_time,
                       cublas_error, (cublas_error < tol ? "ok" : "failed"));
                status += ! (cublas_error < tol);
            }
            else {
                printf("%5d %5d   %7.2f (%7.2f)    ---   (  ---  )    ---     ---\n",
                       (int) N, (int) K,
                       cublas_perf, 1000.*cublas_time);
            }
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_B );
            TESTING_FREE_CPU( h_C );
            TESTING_FREE_CPU( h_Ccublas );
            
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
