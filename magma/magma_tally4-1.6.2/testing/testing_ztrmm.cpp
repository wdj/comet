/*
    -- MAGMA_tally4 (version 1.6.1) --
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
#include "testings.h"  // before magma_tally4.h, to include cublas_v2
#include "flops.h"
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing ztrmm
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, cublas_perf, cublas_time, cpu_perf, cpu_time;
    double          cublas_error, Cnorm, work[1];
    magma_tally4_int_t M, N;
    magma_tally4_int_t Ak;
    magma_tally4_int_t sizeA, sizeB;
    magma_tally4_int_t lda, ldb, ldda, lddb;
    magma_tally4_int_t ione     = 1;
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    
    magma_tally4DoubleComplex *h_A, *h_B, *h_Bcublas;
    magma_tally4DoubleComplex_ptr d_A, d_B;
    magma_tally4DoubleComplex c_neg_one = MAGMA_tally4_Z_NEG_ONE;
    magma_tally4DoubleComplex alpha = MAGMA_tally4_Z_MAKE(  0.29, -0.86 );
    magma_tally4_int_t status = 0;
    
    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("If running lapack (option --lapack), CUBLAS error is computed\n"
           "relative to CPU BLAS result.\n\n");
    printf("side = %s, uplo = %s, transA = %s, diag = %s \n",
           lapack_side_const_tally4(opts.side), lapack_uplo_const_tally4(opts.uplo),
           lapack_trans_const_tally4(opts.transA), lapack_diag_const_tally4(opts.diag) );
    printf("    M     N   CUBLAS Gflop/s (ms)   CPU Gflop/s (ms)  CUBLAS error\n");
    printf("==================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            gflops = FLOPS_ZTRMM(opts.side, M, N) / 1e9;

            if ( opts.side == Magma_tally4Left ) {
                lda = M;
                Ak = M;
            } else {
                lda = N;
                Ak = N;
            }
            
            ldb = M;
            
            ldda = ((lda+31)/32)*32;
            lddb = ((ldb+31)/32)*32;
            
            sizeA = lda*Ak;
            sizeB = ldb*N;
            
            TESTING_MALLOC_CPU( h_A,       magma_tally4DoubleComplex, lda*Ak );
            TESTING_MALLOC_CPU( h_B,       magma_tally4DoubleComplex, ldb*N  );
            TESTING_MALLOC_CPU( h_Bcublas, magma_tally4DoubleComplex, ldb*N  );
            
            TESTING_MALLOC_DEV( d_A, magma_tally4DoubleComplex, ldda*Ak );
            TESTING_MALLOC_DEV( d_B, magma_tally4DoubleComplex, lddb*N  );
            
            /* Initialize the matrices */
            lapackf77_zlarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_zlarnv( &ione, ISEED, &sizeB, h_B );
            
            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            magma_tally4_zsetmatrix( Ak, Ak, h_A, lda, d_A, ldda );
            magma_tally4_zsetmatrix( M, N, h_B, ldb, d_B, lddb );
            
            // note cublas does trmm out-of-place (i.e., adds output matrix C),
            // but allows C=B to do in-place.
            cublas_time = magma_tally4_sync_wtime( NULL );
            cublasZtrmm( opts.handle, cublas_side_const_tally4(opts.side), cublas_uplo_const_tally4(opts.uplo),
                         cublas_trans_const_tally4(opts.transA), cublas_diag_const_tally4(opts.diag),
                         M, N, 
                         &alpha, d_A, ldda,
                                 d_B, lddb,
                                 d_B, lddb );
            cublas_time = magma_tally4_sync_wtime( NULL ) - cublas_time;
            cublas_perf = gflops / cublas_time;
            
            magma_tally4_zgetmatrix( M, N, d_B, lddb, h_Bcublas, ldb );
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally4_wtime();
                blasf77_ztrmm( lapack_side_const_tally4(opts.side), lapack_uplo_const_tally4(opts.uplo),
                               lapack_trans_const_tally4(opts.transA), lapack_diag_const_tally4(opts.diag), 
                               &M, &N,
                               &alpha, h_A, &lda,
                                       h_B, &ldb );
                cpu_time = magma_tally4_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.lapack ) {
                // compute relative error for both magma_tally4 & cublas, relative to lapack,
                // |C_magma_tally4 - C_lapack| / |C_lapack|
                Cnorm = lapackf77_zlange( "M", &M, &N, h_B, &ldb, work );
                
                blasf77_zaxpy( &sizeB, &c_neg_one, h_B, &ione, h_Bcublas, &ione );
                cublas_error = lapackf77_zlange( "M", &M, &N, h_Bcublas, &ldb, work ) / Cnorm;
                
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)    %8.2e   %s\n",
                       (int) M, (int) N,
                       cublas_perf, 1000.*cublas_time,
                       cpu_perf,    1000.*cpu_time,
                       cublas_error, (cublas_error < tol ? "ok" : "failed"));
                status += ! (cublas_error < tol);
            }
            else {
                printf("%5d %5d   %7.2f (%7.2f)    ---   (  ---  )    ---     ---\n",
                       (int) M, (int) N,
                       cublas_perf, 1000.*cublas_time);
            }
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_B );
            TESTING_FREE_CPU( h_Bcublas );
            
            TESTING_FREE_DEV( d_A );
            TESTING_FREE_DEV( d_B );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
