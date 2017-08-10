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

#define h_A(i,j) (h_A + (i) + (j)*lda)

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing ztrsm
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, magma_tally2_perf, magma_tally2_time=0, cublas_perf, cublas_time, cpu_perf=0, cpu_time=0;
    double          magma_tally2_error, cublas_error, lapack_error, work[1];
    magma_tally2_int_t M, N, info;
    magma_tally2_int_t Ak;
    magma_tally2_int_t sizeA, sizeB;
    magma_tally2_int_t lda, ldb, ldda, lddb;
    magma_tally2_int_t ione     = 1;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    magma_tally2_int_t *ipiv;

    magma_tally2DoubleComplex *h_A, *h_B, *h_Bcublas, *h_Bmagma_tally2, *h_Blapack, *h_X;
    magma_tally2DoubleComplex_ptr d_A, d_B;
    magma_tally2DoubleComplex c_neg_one = MAGMA_tally2_Z_NEG_ONE;
    magma_tally2DoubleComplex c_one = MAGMA_tally2_Z_ONE;
    magma_tally2DoubleComplex alpha = MAGMA_tally2_Z_MAKE(  0.29, -0.86 );
    magma_tally2_int_t status = 0;
    
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );
    
    double tol = opts.tolerance * lapackf77_dlamch("E");

    printf("side = %s, uplo = %s, transA = %s, diag = %s \n",
           lapack_side_const_tally2(opts.side), lapack_uplo_const_tally2(opts.uplo),
           lapack_trans_const_tally2(opts.transA), lapack_diag_const_tally2(opts.diag) );
    printf("    M     N  MAGMA_tally2 Gflop/s (ms)  CUBLAS Gflop/s (ms)   CPU Gflop/s (ms)      MAGMA_tally2     CUBLAS   LAPACK error\n");
    printf("============================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            gflops = FLOPS_ZTRSM(opts.side, M, N) / 1e9;

            if ( opts.side == Magma_tally2Left ) {
                lda = M;
                Ak  = M;
            } else {
                lda = N;
                Ak  = N;
            }
            
            ldb = M;
            
            ldda = ((lda+31)/32)*32;
            lddb = ((ldb+31)/32)*32;
            
            sizeA = lda*Ak;
            sizeB = ldb*N;
            
            TESTING_MALLOC_CPU( h_A,       magma_tally2DoubleComplex, lda*Ak  );
            TESTING_MALLOC_CPU( h_B,       magma_tally2DoubleComplex, ldb*N   );
            TESTING_MALLOC_CPU( h_X,       magma_tally2DoubleComplex, ldb*N   );
            TESTING_MALLOC_CPU( h_Blapack, magma_tally2DoubleComplex, ldb*N   );
            TESTING_MALLOC_CPU( h_Bcublas, magma_tally2DoubleComplex, ldb*N   );
            TESTING_MALLOC_CPU( h_Bmagma_tally2,  magma_tally2DoubleComplex, ldb*N   );
            TESTING_MALLOC_CPU( ipiv,      magma_tally2_int_t,        Ak      );
            
            TESTING_MALLOC_DEV( d_A,       magma_tally2DoubleComplex, ldda*Ak );
            TESTING_MALLOC_DEV( d_B,       magma_tally2DoubleComplex, lddb*N  );
            
            /* Initialize the matrices */
            /* Factor A into LU to get well-conditioned triangular matrix.
             * Copy L to U, since L seems okay when used with non-unit diagonal
             * (i.e., from U), while U fails when used with unit diagonal. */
            lapackf77_zlarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_zgetrf( &Ak, &Ak, h_A, &lda, ipiv, &info );
            for( int j = 0; j < Ak; ++j ) {
                for( int i = 0; i < j; ++i ) {
                    *h_A(i,j) = *h_A(j,i);
                }
            }
            
            lapackf77_zlarnv( &ione, ISEED, &sizeB, h_B );
            memcpy( h_Blapack, h_B, sizeB*sizeof(magma_tally2DoubleComplex) );
            
            /* =====================================================================
               Performs operation using MAGMA_tally2BLAS
               =================================================================== */
            magma_tally2_zsetmatrix( Ak, Ak, h_A, lda, d_A, ldda );
            magma_tally2_zsetmatrix( M, N, h_B, ldb, d_B, lddb );
            
            magma_tally2_time = magma_tally2_sync_wtime( NULL );
            magma_tally2blas_ztrsm( opts.side, opts.uplo, opts.transA, opts.diag, 
                             M, N,
                             alpha, d_A, ldda,
                                    d_B, lddb );
            magma_tally2_time = magma_tally2_sync_wtime( NULL ) - magma_tally2_time;
            magma_tally2_perf = gflops / magma_tally2_time;
            
            magma_tally2_zgetmatrix( M, N, d_B, lddb, h_Bmagma_tally2, ldb );
            
            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            magma_tally2_zsetmatrix( M, N, h_B, ldb, d_B, lddb );
            
            cublas_time = magma_tally2_sync_wtime( NULL );
            cublasZtrsm( opts.handle, cublas_side_const_tally2(opts.side), cublas_uplo_const_tally2(opts.uplo),
                         cublas_trans_const_tally2(opts.transA), cublas_diag_const_tally2(opts.diag),
                         M, N, 
                         &alpha, d_A, ldda,
                                 d_B, lddb );
            cublas_time = magma_tally2_sync_wtime( NULL ) - cublas_time;
            cublas_perf = gflops / cublas_time;
            
            magma_tally2_zgetmatrix( M, N, d_B, lddb, h_Bcublas, ldb );
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally2_wtime();
                blasf77_ztrsm( lapack_side_const_tally2(opts.side), lapack_uplo_const_tally2(opts.uplo),
                               lapack_trans_const_tally2(opts.transA), lapack_diag_const_tally2(opts.diag), 
                               &M, &N,
                               &alpha, h_A, &lda,
                                       h_Blapack, &ldb );
                cpu_time = magma_tally2_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            // ||b - 1/alpha*A*x|| / (||A||*||x||)
            magma_tally2DoubleComplex alpha2 = MAGMA_tally2_Z_DIV( c_one, alpha );
            double normR, normX, normA;
            normA = lapackf77_zlange( "M", &Ak, &Ak, h_A, &lda, work );
            
            // check magma_tally2
            memcpy( h_X, h_Bmagma_tally2, sizeB*sizeof(magma_tally2DoubleComplex) );
            blasf77_ztrmm( lapack_side_const_tally2(opts.side), lapack_uplo_const_tally2(opts.uplo),
                           lapack_trans_const_tally2(opts.transA), lapack_diag_const_tally2(opts.diag), 
                           &M, &N,
                           &alpha2, h_A, &lda,
                                    h_X, &ldb );

            blasf77_zaxpy( &sizeB, &c_neg_one, h_B, &ione, h_X, &ione );
            normR = lapackf77_zlange( "M", &M, &N, h_X,      &ldb, work );
            normX = lapackf77_zlange( "M", &M, &N, h_Bmagma_tally2, &ldb, work );
            magma_tally2_error = normR/(normX*normA);

            // check cublas
            memcpy( h_X, h_Bcublas, sizeB*sizeof(magma_tally2DoubleComplex) );
            blasf77_ztrmm( lapack_side_const_tally2(opts.side), lapack_uplo_const_tally2(opts.uplo),
                           lapack_trans_const_tally2(opts.transA), lapack_diag_const_tally2(opts.diag), 
                           &M, &N,
                           &alpha2, h_A, &lda,
                                    h_X, &ldb );

            blasf77_zaxpy( &sizeB, &c_neg_one, h_B, &ione, h_X, &ione );
            normR = lapackf77_zlange( "M", &M, &N, h_X,       &ldb, work );
            normX = lapackf77_zlange( "M", &M, &N, h_Bcublas, &ldb, work );            
            cublas_error = normR/(normX*normA);

            if ( opts.lapack ) {
                // check lapack
                // this verifies that the matrix wasn't so bad that it couldn't be solved accurately.
                memcpy( h_X, h_Blapack, sizeB*sizeof(magma_tally2DoubleComplex) );
                blasf77_ztrmm( lapack_side_const_tally2(opts.side), lapack_uplo_const_tally2(opts.uplo),
                               lapack_trans_const_tally2(opts.transA), lapack_diag_const_tally2(opts.diag), 
                               &M, &N,
                               &alpha2, h_A, &lda,
                                        h_X, &ldb );
    
                blasf77_zaxpy( &sizeB, &c_neg_one, h_B, &ione, h_X, &ione );
                normR = lapackf77_zlange( "M", &M, &N, h_X,       &ldb, work );
                normX = lapackf77_zlange( "M", &M, &N, h_Blapack, &ldb, work );
                lapack_error = normR/(normX*normA);
                
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %8.2e   %8.2e   %s\n",
                        (int) M, (int) N,
                        magma_tally2_perf,  1000.*magma_tally2_time,
                        cublas_perf, 1000.*cublas_time,
                        cpu_perf,    1000.*cpu_time,
                        magma_tally2_error, cublas_error, lapack_error,
                        (magma_tally2_error < tol && cublas_error < tol? "ok" : "failed"));
                status += ! (magma_tally2_error < tol && cublas_error < tol);
            }
            else {
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)     ---   (  ---  )   %8.2e   %8.2e     ---      %s\n",
                        (int) M, (int) N,
                        magma_tally2_perf,  1000.*magma_tally2_time,
                        cublas_perf, 1000.*cublas_time,
                        magma_tally2_error, cublas_error,
                        (magma_tally2_error < tol && cublas_error < tol ? "ok" : "failed"));
                status += ! (magma_tally2_error < tol && cublas_error < tol);
            }
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_B );
            TESTING_FREE_CPU( h_X );
            TESTING_FREE_CPU( h_Blapack );
            TESTING_FREE_CPU( h_Bcublas );
            TESTING_FREE_CPU( h_Bmagma_tally2  );
            TESTING_FREE_CPU( ipiv );
            
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