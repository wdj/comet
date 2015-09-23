/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zlanhe.cpp normal z -> c, Fri Jan 30 19:00:24 2015
       @author Mark Gates
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "testings.h"

#include "magma_minproduct_threadsetting.h"  // to work around MKL bug

#define PRECISION_c

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing clanhe
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gbytes, gpu_perf, gpu_time, cpu_perf, cpu_time;
    magma_minproductFloatComplex *h_A;
    float *h_work;
    magma_minproductFloatComplex_ptr d_A;
    magma_minproductFloat_ptr d_work;
    magma_minproduct_int_t N, n2, lda, ldda;
    magma_minproduct_int_t idist    = 3;  // normal distribution (otherwise max norm is always ~ 1)
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    float      error, norm_magma_minproduct, norm_lapack;
    magma_minproduct_int_t status = 0;
    bool mkl_warning = false;

    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");
    
    magma_minproduct_uplo_t uplo[] = { Magma_minproductLower, Magma_minproductUpper };
    magma_minproduct_norm_t norm[] = { Magma_minproductInfNorm, Magma_minproductOneNorm, Magma_minproductMaxNorm };
    
    // Double-Complex inf-norm not supported on Tesla (CUDA arch 1.x)
#if defined(PRECISION_z)
    magma_minproduct_int_t arch = magma_minproduct_getdevice_arch();
    if ( arch < 200 ) {
        printf("!!!! NOTE: Double-Complex %s and %s norm are not supported\n"
               "!!!! on CUDA architecture %d; requires arch >= 200.\n"
               "!!!! It should report \"parameter number 1 had an illegal value\" below.\n\n",
               Magma_minproductInfNormStr, Magma_minproductOneNormStr, (int) arch );
        for( int inorm = 0; inorm < 2; ++inorm ) {
        for( int iuplo = 0; iuplo < 2; ++iuplo ) {
            printf( "Testing that magma_minproductblas_clanhe( %s, %s, ... ) returns -1 error...\n",
                    lapack_norm_const( norm[inorm] ),
                    lapack_uplo_const( uplo[iuplo] ));
            norm_magma_minproduct = magma_minproductblas_clanhe( norm[inorm], uplo[iuplo], 1, NULL, 1, NULL );
            if ( norm_magma_minproduct != -1 ) {
                printf( "expected magma_minproductblas_clanhe to return -1 error, but got %f\n", norm_magma_minproduct );
                status = 1;
            }
        }}
        printf( "...return values %s\n\n", (status == 0 ? "ok" : "failed") );
    }
#endif

    #ifdef MAGMA_minproduct_WITH_MKL
    printf( "\nNote: using single thread to work around MKL clanhe bug.\n\n" );
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
            gbytes = 0.5*(N+1)*N*sizeof(magma_minproductFloatComplex) / 1e9;
            
            TESTING_MALLOC_CPU( h_A,    magma_minproductFloatComplex, n2 );
            TESTING_MALLOC_CPU( h_work, float, N );
            
            TESTING_MALLOC_DEV( d_A,    magma_minproductFloatComplex, ldda*N );
            TESTING_MALLOC_DEV( d_work, float, N );
            
            /* Initialize the matrix */
            lapackf77_clarnv( &idist, ISEED, &n2, h_A );
            
            magma_minproduct_csetmatrix( N, N, h_A, lda, d_A, ldda );
            
            /* ====================================================================
               Performs operation using MAGMA_minproduct
               =================================================================== */
            gpu_time = magma_minproduct_wtime();
            norm_magma_minproduct = magma_minproductblas_clanhe( norm[inorm], uplo[iuplo], N, d_A, ldda, d_work );
            gpu_time = magma_minproduct_wtime() - gpu_time;
            gpu_perf = gbytes / gpu_time;
            if (norm_magma_minproduct == -1) {
                printf( "%5d   %4c   skipped because it isn't supported on this GPU\n",
                        (int) N, lapacke_norm_const( norm[inorm] ));
                continue;
            }
            if (norm_magma_minproduct < 0)
                printf("magma_minproductblas_clanhe returned error %f: %s.\n",
                       norm_magma_minproduct, magma_minproduct_strerror( (int) norm_magma_minproduct ));
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            #ifdef MAGMA_minproduct_WITH_MKL
            // MKL (11.1.2) has bug in multi-threaded clanhe; use single thread to work around
            int threads = magma_minproduct_get_lapack_numthreads();
            magma_minproduct_set_lapack_numthreads( 1 );
            #endif
            
            cpu_time = magma_minproduct_wtime();
            norm_lapack = lapackf77_clanhe(
                lapack_norm_const( norm[inorm] ),
                lapack_uplo_const( uplo[iuplo] ),
                &N, h_A, &lda, h_work );
            cpu_time = magma_minproduct_wtime() - cpu_time;
            cpu_perf = gbytes / cpu_time;
            if (norm_lapack < 0)
                printf("lapackf77_clanhe returned error %f: %s.\n",
                       norm_lapack, magma_minproduct_strerror( (int) norm_lapack ));
            
            #ifdef MAGMA_minproduct_WITH_MKL
            // end single thread to work around MKL bug
            magma_minproduct_set_lapack_numthreads( threads );
            #endif
            
            /* =====================================================================
               Check the result compared to LAPACK
               Note: MKL (11.1.0) has bug for uplo=Lower with multiple threads.
               Try with $MKL_NUM_THREADS = 1.
               =================================================================== */
            error = fabs( norm_magma_minproduct - norm_lapack ) / norm_lapack;
            float tol2 = tol;
            if ( norm[inorm] == Magma_minproductMaxNorm ) {
                // max-norm depends on only one element, so for Real precisions,
                // MAGMA_minproduct and LAPACK should exactly agree (tol2 = 0),
                // while Complex precisions incur roundoff in cuCabsf.
                #if defined(PRECISION_s) || defined(PRECISION_d)
                tol2 = 0;
                #endif
            }
            
            bool okay = (error <= tol2);
            printf("%5d   %4c   %4c   %7.2f (%7.2f)   %7.2f (%7.2f)   %#9.3g   %s\n",
                   (int) N,
                   lapacke_norm_const( norm[inorm] ),
                   lapacke_uplo_const( uplo[iuplo] ),
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
        printf("* MKL (e.g., 11.1.0) has a bug in clanhe with multiple threads.\n"
               "  Try again with MKL_NUM_THREADS=1.\n" );
    }
    
    TESTING_FINALIZE();
    return status;
}
