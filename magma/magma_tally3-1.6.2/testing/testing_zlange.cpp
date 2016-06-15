/*
    -- MAGMA_tally3 (version 1.6.1) --
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
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zlange
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gbytes, gpu_perf, gpu_time, cpu_perf, cpu_time;
    magma_tally3DoubleComplex *h_A;
    double *h_work;
    magma_tally3DoubleComplex_ptr d_A;
    double *d_work;
    magma_tally3_int_t M, N, n2, lda, ldda, lwork;
    magma_tally3_int_t idist    = 3;  // normal distribution (otherwise max norm is always ~ 1)
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    double      error, norm_magma_tally3, norm_lapack;
    magma_tally3_int_t status = 0;

    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    // Only one norm supported for now, but leave this here for future support
    // of different norms. See similar code in testing_zlanhe.cpp.
    magma_tally3_norm_t norm[] = { Magma_tally3MaxNorm, Magma_tally3OneNorm, Magma_tally3InfNorm };
    
    printf("    M     N   norm   CPU GByte/s (ms)    GPU GByte/s (ms)    error   \n");
    printf("=====================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int inorm = 0; inorm < 3; ++inorm ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M   = opts.msize[itest];
            N   = opts.nsize[itest];
            lda = M;
            n2  = lda*N;
            ldda = roundup( M, opts.roundup );
            if ( norm[inorm] == Magma_tally3OneNorm )
                lwork = N;
            else
                lwork = M;
            // read whole matrix
            gbytes = M*N*sizeof(magma_tally3DoubleComplex) / 1e9;
            
            TESTING_MALLOC_CPU( h_A,    magma_tally3DoubleComplex, n2 );
            TESTING_MALLOC_CPU( h_work, double, M );
            
            TESTING_MALLOC_DEV( d_A,    magma_tally3DoubleComplex, ldda*N );
            TESTING_MALLOC_DEV( d_work, double, lwork );
            
            /* Initialize the matrix */
            lapackf77_zlarnv( &idist, ISEED, &n2, h_A );
            
            // uncomment to test handling NaN
            // MAGMA_tally3 propogates NaN; LAPACK does not necesarily.
            //h_A[ 1 + 1*lda ] = MAGMA_tally3_Z_NAN;
            
            magma_tally3_zsetmatrix( M, N, h_A, lda, d_A, ldda );
            
            /* ====================================================================
               Performs operation using MAGMA_tally3
               =================================================================== */
            gpu_time = magma_tally3_wtime();
            norm_magma_tally3 = magma_tally3blas_zlange( norm[inorm], M, N, d_A, ldda, d_work );
            gpu_time = magma_tally3_wtime() - gpu_time;
            gpu_perf = gbytes / gpu_time;
            if (norm_magma_tally3 < 0)
                printf("magma_tally3blas_zlange returned error %f: %s.\n",
                       norm_magma_tally3, magma_tally3_strerror( (int) norm_magma_tally3 ));
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_tally3_wtime();
            norm_lapack = lapackf77_zlange( lapack_norm_const_tally3(norm[inorm]), &M, &N, h_A, &lda, h_work );
            cpu_time = magma_tally3_wtime() - cpu_time;
            cpu_perf = gbytes / cpu_time;
            if (norm_lapack < 0)
                printf("lapackf77_zlange returned error %f: %s.\n",
                       norm_lapack, magma_tally3_strerror( (int) norm_lapack ));
            
            //printf( "norm %12.8f, lapack %12.8f,  A[1,1] %12.8f+%12.8f  ",
            //         norm_magma_tally3, norm_lapack,
            //         MAGMA_tally3_Z_REAL( h_A[1+1*lda] ),
            //         MAGMA_tally3_Z_IMAG( h_A[1+1*lda] ) );
            
            /* =====================================================================
               Check the result compared to LAPACK
               Max norm should be identical; others should be within tolerance.
               =================================================================== */
            error = fabs( norm_magma_tally3 - norm_lapack ) / norm_lapack;
            double tol2 = tol;
            if ( norm[inorm] == Magma_tally3MaxNorm ) {
                // max-norm depends on only one element, so for Real precisions,
                // MAGMA_tally3 and LAPACK should exactly agree (tol2 = 0),
                // while Complex precisions incur roundoff in cuCabs.
                #if defined(PRECISION_s) || defined(PRECISION_d)
                tol2 = 0;
                #endif
            }
            
            bool okay = (error <= tol2);
            printf("%5d %5d   %4c   %7.2f (%7.2f)   %7.2f (%7.2f)   %#9.3g   %s\n",
                   (int) M, (int) N, lapacke_norm_const_tally3(norm[inorm]),
                   cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000.,
                   error, (okay ? "ok" : "failed") );
            status += ! okay;
            
            TESTING_FREE_CPU( h_A    );
            TESTING_FREE_CPU( h_work );
            
            TESTING_FREE_DEV( d_A    );
            TESTING_FREE_DEV( d_work );
            fflush( stdout );
        }} // end inorm, iter
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
