/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zhesv_nopiv_gpu.cpp normal z -> c, Fri Jan 30 19:00:24 2015
       
       @author Mark Gates
       @author Adrien Remy
*/
// includes, system
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing chesv_nopiv_gpu
*/
int main(int argc, char **argv)
{
    TESTING_INIT();

    real_Double_t   gflops, cpu_perf, cpu_time, gpu_perf, gpu_time;
    float          error, Rnorm, Anorm, Xnorm, *work;
    magma_tally3FloatComplex c_one     = MAGMA_tally3_C_ONE;
    magma_tally3FloatComplex c_neg_one = MAGMA_tally3_C_NEG_ONE;
    magma_tally3FloatComplex *h_A, *h_B, *h_X, temp, *hwork;
    magma_tally3FloatComplex_ptr d_A, d_B;
    magma_tally3_int_t *ipiv;
    magma_tally3_int_t N, nrhs, lda, ldb, ldda, lddb, info, sizeA, sizeB, lwork;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t status = 0;
    
    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");
    
    nrhs = opts.nrhs;
    
    printf("    N  NRHS   CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||B - AX|| / N*||A||*||X||\n");
    printf("================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            ldb    = lda;
            ldda   = ((N+31)/32)*32;
            lddb   = ldda;
            gflops = ( FLOPS_CGETRF( N, N ) + FLOPS_CGETRS( N, nrhs ) ) / 1e9;
            
            TESTING_MALLOC_CPU( h_A, magma_tally3FloatComplex, lda*N    );
            TESTING_MALLOC_CPU( h_B, magma_tally3FloatComplex, ldb*nrhs );
            TESTING_MALLOC_CPU( h_X, magma_tally3FloatComplex, ldb*nrhs );
            TESTING_MALLOC_CPU( work, float,      N );
            TESTING_MALLOC_CPU( ipiv, magma_tally3_int_t, N );
            
            TESTING_MALLOC_DEV( d_A, magma_tally3FloatComplex, ldda*N    );
            TESTING_MALLOC_DEV( d_B, magma_tally3FloatComplex, lddb*nrhs );
            
            /* Initialize the matrices */
            sizeA = lda*N;
            sizeB = ldb*nrhs;
            lapackf77_clarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_clarnv( &ione, ISEED, &sizeB, h_B );
            
            bool nopiv = true;
            if ( nopiv ) {
                magma_tally3_cmake_hpd( N, h_A, lda );  // SPD / HPD does not require pivoting
            }
            else {
                magma_tally3_cmake_hermitian( N, h_A, lda );  // symmetric/Hermitian generally requires pivoting
            }
            
            magma_tally3_csetmatrix( N, N,    h_A, lda, d_A, ldda );
            magma_tally3_csetmatrix( N, nrhs, h_B, ldb, d_B, lddb );
            
            /* ====================================================================
               Performs operation using MAGMA_tally3
               =================================================================== */
            gpu_time = magma_tally3_wtime();

            magma_tally3_chesv_nopiv_gpu( opts.uplo, N, nrhs, d_A, ldda, d_B, lddb, &info );
            gpu_time = magma_tally3_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally3_cgesv_gpu returned error %d: %s.\n",
                       (int) info, magma_tally3_strerror( info ));
            
            //=====================================================================
            // Residual
            //=====================================================================
            magma_tally3_cgetmatrix( N, nrhs, d_B, lddb, h_X, ldb );
            
            Anorm = lapackf77_clange("I", &N, &N,    h_A, &lda, work);
            Xnorm = lapackf77_clange("I", &N, &nrhs, h_X, &ldb, work);
            
            blasf77_cgemm( Magma_tally3NoTransStr, Magma_tally3NoTransStr, &N, &nrhs, &N,
                           &c_one,     h_A, &lda,
                                       h_X, &ldb,
                           &c_neg_one, h_B, &ldb);
            
            Rnorm = lapackf77_clange("I", &N, &nrhs, h_B, &ldb, work);
            error = Rnorm/(N*Anorm*Xnorm);
            status += ! (error < tol);
            
            /* ====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                lwork = -1;
                lapackf77_chesv( lapack_uplo_const_tally3(opts.uplo), &N,&nrhs,
                                 h_A, &lda, ipiv, h_B, &ldb, &temp, &lwork, &info );
                lwork = (magma_tally3_int_t) MAGMA_tally3_C_REAL( temp );
                TESTING_MALLOC_PIN( hwork, magma_tally3FloatComplex, lwork );

                cpu_time = magma_tally3_wtime();
                lapackf77_chesv( lapack_uplo_const_tally3(opts.uplo), &N, &nrhs,
                                 h_A, &lda, ipiv, h_B, &ldb, hwork, &lwork, &info );
                cpu_time = magma_tally3_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_chesv returned error %d: %s.\n",
                           (int) info, magma_tally3_strerror( info ));
                
                printf( "%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                        (int) N, (int) nrhs, cpu_perf, cpu_time, gpu_perf, gpu_time,
                        error, (error < tol ? "ok" : "failed"));
                TESTING_FREE_CPU( hwork );
            }
            else {
                printf( "%5d %5d     ---   (  ---  )   %7.2f (%7.2f)   %8.2e   %s\n",
                        (int) N, (int) nrhs, gpu_perf, gpu_time,
                        error, (error < tol ? "ok" : "failed"));
            }
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_B );
            TESTING_FREE_CPU( h_X );
            TESTING_FREE_CPU( work );
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
