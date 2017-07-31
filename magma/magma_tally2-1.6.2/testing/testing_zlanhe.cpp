/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> c d s
       @author Mark Gates
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "magma_tally2.h"
#include "magma_tally2_lapack.h"
#include "testings.h"

#include "magma_tally2_threadsetting.h"  // to work around MKL bug

#define PRECISION_z

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zlanhe
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gbytes, gpu_perf, gpu_time, cpu_perf, cpu_time;
    magma_tally2DoubleComplex *h_A;
    double *h_work;
    magma_tally2DoubleComplex_ptr d_A;
    magma_tally2Double_ptr d_work;
    magma_tally2_int_t N, n2, lda, ldda;
    magma_tally2_int_t idist    = 3;  // normal distribution (otherwise max norm is always ~ 1)
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    double      error, norm_magma_tally2, norm_lapack;
    magma_tally2_int_t status = 0;
    bool mkl_warning = false;

    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );
    
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    magma_tally2_uplo_t uplo[] = { Magma_tally2Lower, Magma_tally2Upper };
    magma_tally2_norm_t norm[] = { Magma_tally2InfNorm, Magma_tally2OneNorm, Magma_tally2MaxNorm };
    
    // Double-Complex inf-norm not supported on Tesla (CUDA arch 1.x)
#if defined(PRECISION_z)
    magma_tally2_int_t arch = magma_tally2_getdevice_arch();
    if ( arch < 200 ) {
        printf("!!!! NOTE: Double-Complex %s and %s norm are not supported\n"
               "!!!! on CUDA architecture %d; requires arch >= 200.\n"
               "!!!! It should report \"parameter number 1 had an illegal value\" below.\n\n",
               Magma_tally2InfNormStr, Magma_tally2OneNormStr, (int) arch );
        for( int inorm = 0; inorm < 2; ++inorm ) {
        for( int iuplo = 0; iuplo < 2; ++iuplo ) {
            printf( "Testing that magma_tally2blas_zlanhe( %s, %s, ... ) returns -1 error...\n",
                    lapack_norm_const_tally2( norm[inorm] ),
                    lapack_uplo_const_tally2( uplo[iuplo] ));
            norm_magma_tally2 = magma_tally2blas_zlanhe( norm[inorm], uplo[iuplo], 1, NULL, 1, NULL );
            if ( norm_magma_tally2 != -1 ) {
                printf( "expected magma_tally2blas_zlanhe to return -1 error, but got %f\n", norm_magma_tally2 );
                status = 1;
            }
        }}
        printf( "...return values %s\n\n", (status == 0 ? "ok" : "failed") );
    }
#endif

    #ifdef MAGMA_tally2_WITH_MKL
    printf( "\nNote: using single thread to work around MKL zlanhe bug.\n\n" );
    #endif
    
    printf("    N   norm   uplo   CPU GByte/s (ms)    GPU GByte/s (ms)    error   \n");
    printf("=======================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
      for( int inorm = 0; inorm < 3; ++inorm ) {
      for( int iuplo = 0; iuplo < 2; ++iuplo ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N   = opts.nsize[itest];
            lda = N;
            n2  = lda*N;
            ldda = roundup( N, opts.roundup );
            // read upper or lower triangle
            gbytes = 0.5*(N+1)*N*sizeof(magma_tally2DoubleComplex) / 1e9;
            
            TESTING_MALLOC_CPU( h_A,    magma_tally2DoubleComplex, n2 );
            TESTING_MALLOC_CPU( h_work, double, N );
            
            TESTING_MALLOC_DEV( d_A,    magma_tally2DoubleComplex, ldda*N );
            TESTING_MALLOC_DEV( d_work, double, N );
            
            /* Initialize the matrix */
            lapackf77_zlarnv( &idist, ISEED, &n2, h_A );
            
            magma_tally2_zsetmatrix( N, N, h_A, lda, d_A, ldda );
            
            /* ====================================================================
               Performs operation using MAGMA_tally2
               =================================================================== */
            gpu_time = magma_tally2_wtime();
            norm_magma_tally2 = magma_tally2blas_zlanhe( norm[inorm], uplo[iuplo], N, d_A, ldda, d_work );
            gpu_time = magma_tally2_wtime() - gpu_time;
            gpu_perf = gbytes / gpu_time;
            if (norm_magma_tally2 == -1) {
                printf( "%5d   %4c   skipped because it isn't supported on this GPU\n",
                        (int) N, lapacke_norm_const_tally2( norm[inorm] ));
                continue;
            }
            if (norm_magma_tally2 < 0)
                printf("magma_tally2blas_zlanhe returned error %f: %s.\n",
                       norm_magma_tally2, magma_tally2_strerror( (int) norm_magma_tally2 ));
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            #ifdef MAGMA_tally2_WITH_MKL
            // MKL (11.1.2) has bug in multi-threaded zlanhe; use single thread to work around
            int threads = magma_tally2_get_lapack_numthreads();
            magma_tally2_set_lapack_numthreads( 1 );
            #endif
            
            cpu_time = magma_tally2_wtime();
            norm_lapack = lapackf77_zlanhe(
                lapack_norm_const_tally2( norm[inorm] ),
                lapack_uplo_const_tally2( uplo[iuplo] ),
                &N, h_A, &lda, h_work );
            cpu_time = magma_tally2_wtime() - cpu_time;
            cpu_perf = gbytes / cpu_time;
            if (norm_lapack < 0)
                printf("lapackf77_zlanhe returned error %f: %s.\n",
                       norm_lapack, magma_tally2_strerror( (int) norm_lapack ));
            
            #ifdef MAGMA_tally2_WITH_MKL
            // end single thread to work around MKL bug
            magma_tally2_set_lapack_numthreads( threads );
            #endif
            
            /* =====================================================================
               Check the result compared to LAPACK
               Note: MKL (11.1.0) has bug for uplo=Lower with multiple threads.
               Try with $MKL_NUM_THREADS = 1.
               =================================================================== */
            error = fabs( norm_magma_tally2 - norm_lapack ) / norm_lapack;
            double tol2 = tol;
            if ( norm[inorm] == Magma_tally2MaxNorm ) {
                // max-norm depends on only one element, so for Real precisions,
                // MAGMA_tally2 and LAPACK should exactly agree (tol2 = 0),
                // while Complex precisions incur roundoff in cuCabs.
                #if defined(PRECISION_s) || defined(PRECISION_d)
                tol2 = 0;
                #endif
            }
            
            bool okay = (error <= tol2);
            printf("%5d   %4c   %4c   %7.2f (%7.2f)   %7.2f (%7.2f)   %#9.3g   %s\n",
                   (int) N,
                   lapacke_norm_const_tally2( norm[inorm] ),
                   lapacke_uplo_const_tally2( uplo[iuplo] ),
                   cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000.,
                   error, (okay ? "ok" : "failed") );
            status += ! okay;
            
            if ( ! okay ) {
                mkl_warning = true;
            }
            
            TESTING_FREE_CPU( h_A    );
            TESTING_FREE_CPU( h_work );
            
            TESTING_FREE_DEV( d_A    );
            TESTING_FREE_DEV( d_work );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
      }} // end iuplo, inorm, iter
      printf( "\n" );
    }
    
    if ( mkl_warning ) {
        printf("* MKL (e.g., 11.1.0) has a bug in zlanhe with multiple threads.\n"
               "  Try again with MKL_NUM_THREADS=1.\n" );
    }
    
    TESTING_FINALIZE();
    return status;
}
