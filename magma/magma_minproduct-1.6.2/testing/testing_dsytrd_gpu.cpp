/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Stan Tomov
       @author Azzam Haidar
       @author Mark Gates

       @generated from testing_zhetrd_gpu.cpp normal z -> d, Fri Jan 30 19:00:26 2015

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "testings.h"

#define REAL

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dsytrd_gpu
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gflops, gpu_perf, cpu_perf, gpu_time, cpu_time;
    double           eps;
    double *h_A, *h_R, *h_Q, *h_work, *work;
    magma_minproductDouble_ptr d_R, dwork;
    double *tau;
    double          *diag, *offdiag;
    double           result[2] = {0., 0.};
    magma_minproduct_int_t N, n2, lda, ldda, lwork, info, nb, ldwork;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t itwo     = 2;
    magma_minproduct_int_t ithree   = 3;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t status = 0;
    
    #ifdef COMPLEX
    double *rwork;
    #endif

    eps = lapackf77_dlamch( "E" );

    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );

    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("Running version %d; available are (specified through --version num):\n", 
           (int) opts.version);
    printf("1 - uses DSYMV from CUBLAS (default)\n");
    printf("2 - uses DSYMV from MAGMA_minproduct BLAS that requires extra space\n\n");

    printf("uplo = %s\n", lapack_uplo_const(opts.uplo) );
    printf("  N     CPU GFlop/s (sec)   GPU GFlop/s (sec)   |A-QHQ'|/N|A|   |I-QQ'|/N\n");
    printf("===========================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            ldda   = roundup( N, opts.roundup );  // by default, round to multiple of 32
            n2     = N*lda;
            nb     = magma_minproduct_get_dsytrd_nb(N);
            lwork  = N*nb;  /* We suppose the magma_minproduct nb is bigger than lapack nb */
            gflops = FLOPS_DSYTRD( N ) / 1e9;
            ldwork = ldda*ceildiv(N,64) + 2*ldda*nb;
            
            TESTING_MALLOC_CPU( h_A,     double, lda*N );
            TESTING_MALLOC_CPU( tau,     double, N     );
            TESTING_MALLOC_CPU( diag,    double, N   );
            TESTING_MALLOC_CPU( offdiag, double, N-1 );
            
            TESTING_MALLOC_PIN( h_R,     double, lda*N );
            TESTING_MALLOC_PIN( h_work,  double, lwork );
            
            TESTING_MALLOC_DEV( d_R,     double, ldda*N );
            TESTING_MALLOC_DEV( dwork,   double, ldwork );
            
            /* ====================================================================
               Initialize the matrix
               =================================================================== */
            lapackf77_dlarnv( &ione, ISEED, &n2, h_A );
            magma_minproduct_dmake_symmetric( N, h_A, lda );
            magma_minproduct_dsetmatrix( N, N, h_A, lda, d_R, ldda );
            
            /* ====================================================================
               Performs operation using MAGMA_minproduct
               =================================================================== */
            gpu_time = magma_minproduct_wtime();
            if (opts.version == 1) {
                magma_minproduct_dsytrd_gpu( opts.uplo, N, d_R, ldda, diag, offdiag,
                                  tau, h_R, lda, h_work, lwork, &info );
            }
            else {
                magma_minproduct_dsytrd2_gpu( opts.uplo, N, d_R, ldda, diag, offdiag,
                                   tau, h_R, lda, h_work, lwork, dwork, ldwork, &info );
            }
            gpu_time = magma_minproduct_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_minproduct_dsytrd_gpu returned error %d: %s.\n",
                       (int) info, magma_minproduct_strerror( info ));
            
            /* =====================================================================
               Check the factorization
               =================================================================== */
            if ( opts.check ) {
                TESTING_MALLOC_CPU( h_Q,  double, lda*N );
                TESTING_MALLOC_CPU( work, double, 2*N*N );
                #ifdef COMPLEX
                TESTING_MALLOC_CPU( rwork, double, N );
                #endif
                
                magma_minproduct_dgetmatrix( N, N, d_R, ldda, h_R, lda );
                magma_minproduct_dgetmatrix( N, N, d_R, ldda, h_Q, lda );
                lapackf77_dorgtr( lapack_uplo_const(opts.uplo), &N, h_Q, &lda, tau, h_work, &lwork, &info );
                
                lapackf77_dsyt21( &itwo, lapack_uplo_const(opts.uplo), &N, &ione,
                                  h_A, &lda, diag, offdiag,
                                  h_Q, &lda, h_R, &lda,
                                  tau, work,
                                  #ifdef COMPLEX
                                  rwork,
                                  #endif
                                  &result[0] );
                
                lapackf77_dsyt21( &ithree, lapack_uplo_const(opts.uplo), &N, &ione,
                                  h_A, &lda, diag, offdiag,
                                  h_Q, &lda, h_R, &lda,
                                  tau, work,
                                  #ifdef COMPLEX
                                  rwork,
                                  #endif
                                  &result[1] );
                result[0] *= eps;
                result[1] *= eps;
                
                TESTING_FREE_CPU( h_Q  );
                TESTING_FREE_CPU( work );
                #ifdef COMPLEX
                TESTING_FREE_CPU( rwork );
                #endif
            }
                        
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_minproduct_wtime();
                lapackf77_dsytrd( lapack_uplo_const(opts.uplo), &N, h_A, &lda, diag, offdiag, tau,
                                  h_work, &lwork, &info );
                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_dsytrd returned error %d: %s.\n",
                           (int) info, magma_minproduct_strerror( info ));
            }
            
            /* =====================================================================
               Print performance and error.
               =================================================================== */
            if ( opts.lapack ) {
                printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)",
                       (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time );
            } else {
                printf("%5d     ---   (  ---  )   %7.2f (%7.2f)",
                       (int) N, gpu_perf, gpu_time );
            }
            if ( opts.check ) {
                printf("   %8.2e        %8.2e   %s\n", result[0], result[1],
                        ((result[0] < tol && result[1] < tol) ? "ok" : "failed")  );
                status += ! (result[0] < tol && result[1] < tol);
            } else {
                printf("     ---             ---\n");
            }
            
            TESTING_FREE_CPU( h_A     );
            TESTING_FREE_CPU( tau     );
            TESTING_FREE_CPU( diag    );
            TESTING_FREE_CPU( offdiag );
            
            TESTING_FREE_PIN( h_R    );
            TESTING_FREE_PIN( h_work );
            
            TESTING_FREE_DEV( d_R   );
            TESTING_FREE_DEV( dwork );
            
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
