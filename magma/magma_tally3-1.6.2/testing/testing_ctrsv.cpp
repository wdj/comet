/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_ztrsv.cpp normal z -> c, Fri Jan 30 19:00:23 2015
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

#define h_A(i,j) (h_A + (i) + (j)*lda)

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing ctrsm
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, cublas_perf, cublas_time, cpu_perf=0, cpu_time=0;
    float          cublas_error, normA, normx, normr, work[1];
    magma_tally3_int_t N, info;
    magma_tally3_int_t sizeA;
    magma_tally3_int_t lda, ldda;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t *ipiv;

    magma_tally3FloatComplex *h_A, *h_b, *h_x, *h_xcublas;
    magma_tally3FloatComplex_ptr d_A, d_x;
    magma_tally3FloatComplex c_neg_one = MAGMA_tally3_C_NEG_ONE;
    magma_tally3_int_t status = 0;
    
    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");
    
    printf("uplo = %s, transA = %s, diag = %s\n",
           lapack_uplo_const_tally3(opts.uplo), lapack_trans_const_tally3(opts.transA), lapack_diag_const_tally3(opts.diag) );
    printf("    N  CUBLAS Gflop/s (ms)   CPU Gflop/s (ms)   CUBLAS error\n");
    printf("============================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            gflops = FLOPS_CTRSM(opts.side, N, 1) / 1e9;
            lda    = N;
            ldda   = ((lda+31)/32)*32;
            sizeA  = lda*N;
            
            TESTING_MALLOC_CPU( ipiv,      magma_tally3_int_t,        N     );
            TESTING_MALLOC_CPU( h_A,       magma_tally3FloatComplex, lda*N );
            TESTING_MALLOC_CPU( h_b,       magma_tally3FloatComplex, N     );
            TESTING_MALLOC_CPU( h_x,       magma_tally3FloatComplex, N     );
            TESTING_MALLOC_CPU( h_xcublas, magma_tally3FloatComplex, N     );
            
            TESTING_MALLOC_DEV( d_A, magma_tally3FloatComplex, ldda*N );
            TESTING_MALLOC_DEV( d_x, magma_tally3FloatComplex, N      );
            
            /* Initialize the matrices */
            /* Factor A into LU to get well-conditioned triangular matrix.
             * Copy L to U, since L seems okay when used with non-unit diagonal
             * (i.e., from U), while U fails when used with unit diagonal. */
            lapackf77_clarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_cgetrf( &N, &N, h_A, &lda, ipiv, &info );
            for( int j = 0; j < N; ++j ) {
                for( int i = 0; i < j; ++i ) {
                    *h_A(i,j) = *h_A(j,i);
                }
            }
            
            lapackf77_clarnv( &ione, ISEED, &N, h_b );
            blasf77_ccopy( &N, h_b, &ione, h_x, &ione );
            
            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            magma_tally3_csetmatrix( N, N, h_A, lda, d_A, ldda );
            magma_tally3_csetvector( N, h_x, 1, d_x, 1 );
            
            cublas_time = magma_tally3_sync_wtime( NULL );
            cublasCtrsv( opts.handle, cublas_uplo_const_tally3(opts.uplo),
                         cublas_trans_const_tally3(opts.transA), cublas_diag_const_tally3(opts.diag),
                         N,
                         d_A, ldda,
                         d_x, 1 );
            cublas_time = magma_tally3_sync_wtime( NULL ) - cublas_time;
            cublas_perf = gflops / cublas_time;
            
            magma_tally3_cgetvector( N, d_x, 1, h_xcublas, 1 );
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally3_wtime();
                blasf77_ctrsv( lapack_uplo_const_tally3(opts.uplo), lapack_trans_const_tally3(opts.transA), lapack_diag_const_tally3(opts.diag),
                               &N,
                               h_A, &lda,
                               h_x, &ione );
                cpu_time = magma_tally3_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            // ||b - Ax|| / (||A||*||x||)
            // error for CUBLAS
            normA = lapackf77_clange( "F", &N, &N, h_A, &lda, work );
            
            normx = lapackf77_clange( "F", &N, &ione, h_xcublas, &ione, work );
            blasf77_ctrmv( lapack_uplo_const_tally3(opts.uplo), lapack_trans_const_tally3(opts.transA), lapack_diag_const_tally3(opts.diag),
                           &N,
                           h_A, &lda,
                           h_xcublas, &ione );
            blasf77_caxpy( &N, &c_neg_one, h_b, &ione, h_xcublas, &ione );
            normr = lapackf77_clange( "F", &N, &ione, h_xcublas, &N, work );
            cublas_error = normr / (normA*normx);

            if ( opts.lapack ) {
                printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                        (int) N,
                        cublas_perf, 1000.*cublas_time,
                        cpu_perf,    1000.*cpu_time,
                        cublas_error, (cublas_error < tol ? "ok" : "failed"));
                status += ! (cublas_error < tol);
            }
            else {
                printf("%5d   %7.2f (%7.2f)     ---  (  ---  )   %8.2e   %s\n",
                        (int) N,
                        cublas_perf, 1000.*cublas_time,
                        cublas_error, (cublas_error < tol ? "ok" : "failed"));
                status += ! (cublas_error < tol);
            }
            
            TESTING_FREE_CPU( ipiv );
            TESTING_FREE_CPU( h_A  );
            TESTING_FREE_CPU( h_b  );
            TESTING_FREE_CPU( h_x  );
            TESTING_FREE_CPU( h_xcublas );
            
            TESTING_FREE_DEV( d_A );
            TESTING_FREE_DEV( d_x );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
