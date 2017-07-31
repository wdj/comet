/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zlange.cpp normal z -> s, Fri Jan 30 19:00:23 2015
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

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing slange
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gbytes, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float *h_A;
    float *h_work;
    magma_tally2Float_ptr d_A;
    float *d_work;
    magma_tally2_int_t M, N, n2, lda, ldda, lwork;
    magma_tally2_int_t idist    = 3;  // normal distribution (otherwise max norm is always ~ 1)
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    float      error, norm_magma_tally2, norm_lapack;
    magma_tally2_int_t status = 0;

    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");
    
    // Only one norm supported for now, but leave this here for future support
    // of different norms. See similar code in testing_slansy.cpp.
    magma_tally2_norm_t norm[] = { Magma_tally2MaxNorm, Magma_tally2OneNorm, Magma_tally2InfNorm };
    
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
            if ( norm[inorm] == Magma_tally2OneNorm )
                lwork = N;
            else
                lwork = M;
            // read whole matrix
            gbytes = M*N*sizeof(float) / 1e9;
            
            TESTING_MALLOC_CPU( h_A,    float, n2 );
            TESTING_MALLOC_CPU( h_work, float, M );
            
            TESTING_MALLOC_DEV( d_A,    float, ldda*N );
            TESTING_MALLOC_DEV( d_work, float, lwork );
            
            /* Initialize the matrix */
            lapackf77_slarnv( &idist, ISEED, &n2, h_A );
            
            // uncomment to test handling NaN
            // MAGMA_tally2 propogates NaN; LAPACK does not necesarily.
            //h_A[ 1 + 1*lda ] = MAGMA_tally2_S_NAN;
            
            magma_tally2_ssetmatrix( M, N, h_A, lda, d_A, ldda );
            
            /* ====================================================================
               Performs operation using MAGMA_tally2
               =================================================================== */
            gpu_time = magma_tally2_wtime();
            norm_magma_tally2 = magma_tally2blas_slange( norm[inorm], M, N, d_A, ldda, d_work );
            gpu_time = magma_tally2_wtime() - gpu_time;
            gpu_perf = gbytes / gpu_time;
            if (norm_magma_tally2 < 0)
                printf("magma_tally2blas_slange returned error %f: %s.\n",
                       norm_magma_tally2, magma_tally2_strerror( (int) norm_magma_tally2 ));
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_tally2_wtime();
            norm_lapack = lapackf77_slange( lapack_norm_const_tally2(norm[inorm]), &M, &N, h_A, &lda, h_work );
            cpu_time = magma_tally2_wtime() - cpu_time;
            cpu_perf = gbytes / cpu_time;
            if (norm_lapack < 0)
                printf("lapackf77_slange returned error %f: %s.\n",
                       norm_lapack, magma_tally2_strerror( (int) norm_lapack ));
            
            //printf( "norm %12.8f, lapack %12.8f,  A[1,1] %12.8f+%12.8f  ",
            //         norm_magma_tally2, norm_lapack,
            //         MAGMA_tally2_S_REAL( h_A[1+1*lda] ),
            //         MAGMA_tally2_S_IMAG( h_A[1+1*lda] ) );
            
            /* =====================================================================
               Check the result compared to LAPACK
               Max norm should be identical; others should be within tolerance.
               =================================================================== */
            error = fabs( norm_magma_tally2 - norm_lapack ) / norm_lapack;
            float tol2 = tol;
            if ( norm[inorm] == Magma_tally2MaxNorm ) {
                // max-norm depends on only one element, so for Real precisions,
                // MAGMA_tally2 and LAPACK should exactly agree (tol2 = 0),
                // while Complex precisions incur roundoff in fabsf.
                #if defined(PRECISION_s) || defined(PRECISION_d)
                tol2 = 0;
                #endif
            }
            
            bool okay = (error <= tol2);
            printf("%5d %5d   %4c   %7.2f (%7.2f)   %7.2f (%7.2f)   %#9.3g   %s\n",
                   (int) M, (int) N, lapacke_norm_const_tally2(norm[inorm]),
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
