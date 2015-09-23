/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zhesv_nopiv_gpu.cpp normal z -> d, Fri Jan 30 19:00:24 2015
       
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
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dsysv_nopiv_gpu
*/
int main(int argc, char **argv)
{
    TESTING_INIT();

    real_Double_t   gflops, cpu_perf, cpu_time, gpu_perf, gpu_time;
    double          error, Rnorm, Anorm, Xnorm, *work;
    double c_one     = MAGMA_minproduct_D_ONE;
    double c_neg_one = MAGMA_minproduct_D_NEG_ONE;
    double *h_A, *h_B, *h_X, temp, *hwork;
    magma_minproductDouble_ptr d_A, d_B;
    magma_minproduct_int_t *ipiv;
    magma_minproduct_int_t N, nrhs, lda, ldb, ldda, lddb, info, sizeA, sizeB, lwork;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t status = 0;
    
    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
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
            gflops = ( FLOPS_DGETRF( N, N ) + FLOPS_DGETRS( N, nrhs ) ) / 1e9;
            
            TESTING_MALLOC_CPU( h_A, double, lda*N    );
            TESTING_MALLOC_CPU( h_B, double, ldb*nrhs );
            TESTING_MALLOC_CPU( h_X, double, ldb*nrhs );
            TESTING_MALLOC_CPU( work, double,      N );
            TESTING_MALLOC_CPU( ipiv, magma_minproduct_int_t, N );
            
            TESTING_MALLOC_DEV( d_A, double, ldda*N    );
            TESTING_MALLOC_DEV( d_B, double, lddb*nrhs );
            
            /* Initialize the matrices */
            sizeA = lda*N;
            sizeB = ldb*nrhs;
            lapackf77_dlarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_dlarnv( &ione, ISEED, &sizeB, h_B );
            
            bool nopiv = true;
            if ( nopiv ) {
                magma_minproduct_dmake_hpd( N, h_A, lda );  // SPD / HPD does not require pivoting
            }
            else {
                magma_minproduct_dmake_symmetric( N, h_A, lda );  // symmetric/symmetric generally requires pivoting
            }
            
            magma_minproduct_dsetmatrix( N, N,    h_A, lda, d_A, ldda );
            magma_minproduct_dsetmatrix( N, nrhs, h_B, ldb, d_B, lddb );
            
            /* ====================================================================
               Performs operation using MAGMA_minproduct
               =================================================================== */
            gpu_time = magma_minproduct_wtime();

            magma_minproduct_dsysv_nopiv_gpu( opts.uplo, N, nrhs, d_A, ldda, d_B, lddb, &info );
            gpu_time = magma_minproduct_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_minproduct_dgesv_gpu returned error %d: %s.\n",
                       (int) info, magma_minproduct_strerror( info ));
            
            //=====================================================================
            // Residual
            //=====================================================================
            magma_minproduct_dgetmatrix( N, nrhs, d_B, lddb, h_X, ldb );
            
            Anorm = lapackf77_dlange("I", &N, &N,    h_A, &lda, work);
            Xnorm = lapackf77_dlange("I", &N, &nrhs, h_X, &ldb, work);
            
            blasf77_dgemm( Magma_minproductNoTransStr, Magma_minproductNoTransStr, &N, &nrhs, &N,
                           &c_one,     h_A, &lda,
                                       h_X, &ldb,
                           &c_neg_one, h_B, &ldb);
            
            Rnorm = lapackf77_dlange("I", &N, &nrhs, h_B, &ldb, work);
            error = Rnorm/(N*Anorm*Xnorm);
            status += ! (error < tol);
            
            /* ====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                lwork = -1;
                lapackf77_dsysv( lapack_uplo_const(opts.uplo), &N,&nrhs,
                                 h_A, &lda, ipiv, h_B, &ldb, &temp, &lwork, &info );
                lwork = (magma_minproduct_int_t) MAGMA_minproduct_D_REAL( temp );
                TESTING_MALLOC_PIN( hwork, double, lwork );

                cpu_time = magma_minproduct_wtime();
                lapackf77_dsysv( lapack_uplo_const(opts.uplo), &N, &nrhs,
                                 h_A, &lda, ipiv, h_B, &ldb, hwork, &lwork, &info );
                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_dsysv returned error %d: %s.\n",
                           (int) info, magma_minproduct_strerror( info ));
                
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
