/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @author Mark Gates

       @generated from testing_zhetrd_mgpu.cpp normal z -> c, Fri Jan 30 19:00:26 2015

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "testings.h"

#define COMPLEX

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing chetrd
*/

int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gflops, gpu_perf, cpu_perf, gpu_time, cpu_time;
    magma_tally3FloatComplex *h_A, *h_R, *h_Q, *h_work, *work;
    magma_tally3FloatComplex *tau;
    float *diag, *offdiag;
    #ifdef COMPLEX
    float *rwork = NULL;
    #endif
    float result[2] = {0., 0.};
    magma_tally3_int_t N, n2, lda, lwork, info, nb;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t itwo     = 2;
    magma_tally3_int_t ithree   = 3;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t status = 0;
    magma_tally3_int_t k = 1;  // TODO: UNKNOWN, UNDOCUMENTED VARIABLE (number of streams?)

    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");
    float eps = lapackf77_slamch( "E" );

    /* To avoid uninitialized variable warning */
    h_Q   = NULL;
    work  = NULL;

    printf("uplo = %s, ngpu %d\n", lapack_uplo_const_tally3(opts.uplo), (int) opts.ngpu );
    printf("  N     CPU GFlop/s (sec)   GPU GFlop/s (sec)   |A-QHQ'|/N|A|   |I-QQ'|/N\n");
    printf("===========================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N      = opts.nsize[itest];
            lda    = N;
            n2     = N*lda;
            nb     = magma_tally3_get_chetrd_nb(N);
            /* We suppose the magma_tally3 nb is bigger than lapack nb */
            lwork  = N*nb;
            gflops = FLOPS_CHETRD( N ) / 1e9;
            
            /* Allocate host memory for the matrix */
            TESTING_MALLOC_PIN( h_R,     magma_tally3FloatComplex, lda*N );
            TESTING_MALLOC_PIN( h_work,  magma_tally3FloatComplex, lwork );
            
            TESTING_MALLOC_CPU( h_A,     magma_tally3FloatComplex, lda*N );
            TESTING_MALLOC_CPU( tau,     magma_tally3FloatComplex, N     );
            TESTING_MALLOC_CPU( diag,    float, N   );
            TESTING_MALLOC_CPU( offdiag, float, N-1 );
            
            /* ====================================================================
               Initialize the matrix
               =================================================================== */
            lapackf77_clarnv( &ione, ISEED, &n2, h_A );
            magma_tally3_cmake_hermitian( N, h_A, lda );
            
            lapackf77_clacpy( Magma_tally3UpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
            
            /* ====================================================================
               Performs operation using MAGMA_tally3
               =================================================================== */
            gpu_time = magma_tally3_wtime();
            if ( opts.ngpu == 0 ) {
                magma_tally3_chetrd(opts.uplo, N, h_R, lda, diag, offdiag,
                             tau, h_work, lwork, &info);
            } else {
                magma_tally3_chetrd_mgpu(opts.ngpu, k, opts.uplo, N, h_R, lda, diag, offdiag,
                                  tau, h_work, lwork, &info);
            }
            gpu_time = magma_tally3_wtime() - gpu_time;
            if ( info != 0 ) {
                printf("magma_tally3_chetrd_mgpu returned error %d: %s.\n",
                       (int) info, magma_tally3_strerror( info ));
            }
            
            gpu_perf = gflops / gpu_time;
            
            /* =====================================================================
               Check the factorization
               =================================================================== */
            if ( opts.check ) {
                TESTING_MALLOC_CPU( h_Q,  magma_tally3FloatComplex, lda*N );
                TESTING_MALLOC_CPU( work, magma_tally3FloatComplex, 2*N*N );
                #ifdef COMPLEX
                TESTING_MALLOC_CPU( rwork, float, N );
                #endif
                
                lapackf77_clacpy( lapack_uplo_const_tally3(opts.uplo), &N, &N, h_R, &lda, h_Q, &lda);
                lapackf77_cungtr( lapack_uplo_const_tally3(opts.uplo), &N, h_Q, &lda, tau, h_work, &lwork, &info);

                lapackf77_chet21( &itwo, lapack_uplo_const_tally3(opts.uplo), &N, &ione,
                                  h_A, &lda, diag, offdiag,
                                  h_Q, &lda, h_R, &lda,
                                  tau, work,
                                  #ifdef COMPLEX
                                  rwork,
                                  #endif
                                  &result[0] );
                
                lapackf77_chet21( &ithree, lapack_uplo_const_tally3(opts.uplo), &N, &ione,
                                  h_A, &lda, diag, offdiag,
                                  h_Q, &lda, h_R, &lda,
                                  tau, work,
                                  #ifdef COMPLEX
                                  rwork,
                                  #endif
                                  &result[1] );
                result[0] *= eps;
                result[1] *= eps;

                TESTING_FREE_CPU( h_Q );
                TESTING_FREE_CPU( work );
                #ifdef COMPLEX
                TESTING_FREE_CPU( rwork );
                #endif
            }
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally3_wtime();
                lapackf77_chetrd(lapack_uplo_const_tally3(opts.uplo), &N, h_A, &lda, diag, offdiag, tau,
                                 h_work, &lwork, &info);
                cpu_time = magma_tally3_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if ( info != 0 ) {
                    printf("lapackf77_chetrd returned error %d: %s.\n",
                           (int) info, magma_tally3_strerror( info ));
                }
            }
            
            /* =====================================================================
               Print performance and error.
               =================================================================== */
            if ( opts.lapack ) {
                printf("%5d   %7.2f (%7.2f)", (int) N, cpu_perf, cpu_time );
            }
            else {
                printf("%5d     ---   (  ---  )", (int) N );
            }
            printf("   %7.2f (%7.2f)", gpu_perf, gpu_time );
            if ( opts.check ) {
                printf("   %8.2e        %8.2e   %s\n",
                       result[0], result[1], (result[0] < tol && result[1] < tol ? "ok" : "failed") );
                status += ! (result[0] < tol && result[1] < tol);
            }
            else {
                printf("     ---  \n");
            }

            TESTING_FREE_PIN( h_R );
            TESTING_FREE_PIN( h_work );

            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( tau );
            TESTING_FREE_CPU( diag );
            TESTING_FREE_CPU( offdiag );
            
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}