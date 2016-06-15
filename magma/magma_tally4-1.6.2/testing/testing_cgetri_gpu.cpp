/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zgetri_gpu.cpp normal z -> c, Fri Jan 30 19:00:25 2015
       @author Mark Gates
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cgetrf
*/
int main( int argc, char** argv )
{
    TESTING_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    magma_tally4FloatComplex *h_A, *h_R, *work;
    magma_tally4FloatComplex_ptr d_A, dwork;
    magma_tally4FloatComplex c_neg_one = MAGMA_tally4_C_NEG_ONE;
    magma_tally4_int_t N, n2, lda, ldda, info, lwork, ldwork;
    magma_tally4_int_t ione     = 1;
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    magma_tally4FloatComplex tmp;
    float error, rwork[1];
    magma_tally4_int_t *ipiv;
    magma_tally4_int_t status = 0;
    
    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    // need looser bound (3000*eps instead of 30*eps) for tests
    // TODO: should compute ||I - A*A^{-1}|| / (n*||A||*||A^{-1}||)
    opts.tolerance = max( 3000., opts.tolerance );
    float tol = opts.tolerance * lapackf77_slamch("E");
    
    printf("    N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R||_F / (N*||A||_F)\n");
    printf("=================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            n2     = lda*N;
            ldda   = ((N+31)/32)*32;
            ldwork = N * magma_tally4_get_cgetri_nb( N );
            gflops = FLOPS_CGETRI( N ) / 1e9;
            
            // query for workspace size
            lwork = -1;
            lapackf77_cgetri( &N, NULL, &lda, NULL, &tmp, &lwork, &info );
            if (info != 0)
                printf("lapackf77_cgetri returned error %d: %s.\n",
                       (int) info, magma_tally4_strerror( info ));
            lwork = int( MAGMA_tally4_C_REAL( tmp ));
            
            TESTING_MALLOC_CPU( ipiv,  magma_tally4_int_t,        N      );
            TESTING_MALLOC_CPU( work,  magma_tally4FloatComplex, lwork  );
            TESTING_MALLOC_CPU( h_A,   magma_tally4FloatComplex, n2     );
            
            TESTING_MALLOC_PIN( h_R,   magma_tally4FloatComplex, n2     );
            
            TESTING_MALLOC_DEV( d_A,   magma_tally4FloatComplex, ldda*N );
            TESTING_MALLOC_DEV( dwork, magma_tally4FloatComplex, ldwork );
            
            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &n2, h_A );
            error = lapackf77_clange( "f", &N, &N, h_A, &lda, rwork );  // norm(A)
            
            /* Factor the matrix. Both MAGMA_tally4 and LAPACK will use this factor. */
            magma_tally4_csetmatrix( N, N, h_A, lda, d_A, ldda );
            magma_tally4_cgetrf_gpu( N, N, d_A, ldda, ipiv, &info );
            magma_tally4_cgetmatrix( N, N, d_A, ldda, h_A, lda );
            if ( info != 0 )
                printf("magma_tally4_cgetrf_gpu returned error %d: %s.\n",
                       (int) info, magma_tally4_strerror( info ));
            
            // check for exact singularity
            //h_A[ 10 + 10*lda ] = MAGMA_tally4_C_MAKE( 0.0, 0.0 );
            //magma_tally4_csetmatrix( N, N, h_A, lda, d_A, ldda );
            
            /* ====================================================================
               Performs operation using MAGMA_tally4
               =================================================================== */
            gpu_time = magma_tally4_wtime();
            magma_tally4_cgetri_gpu( N, d_A, ldda, ipiv, dwork, ldwork, &info );
            gpu_time = magma_tally4_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally4_cgetri_gpu returned error %d: %s.\n",
                       (int) info, magma_tally4_strerror( info ));
            
            magma_tally4_cgetmatrix( N, N, d_A, ldda, h_R, lda );
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally4_wtime();
                lapackf77_cgetri( &N, h_A, &lda, ipiv, work, &lwork, &info );
                cpu_time = magma_tally4_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_cgetri returned error %d: %s.\n",
                           (int) info, magma_tally4_strerror( info ));
                
                /* =====================================================================
                   Check the result compared to LAPACK
                   =================================================================== */
                blasf77_caxpy( &n2, &c_neg_one, h_A, &ione, h_R, &ione );
                error = lapackf77_clange( "f", &N, &N, h_R, &lda, rwork ) / (N*error);
                
                printf( "%5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                        (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time,
                        error, (error < tol ? "ok" : "failed"));
                status += ! (error < tol);
            }
            else {
                printf( "%5d     ---   (  ---  )   %7.2f (%7.2f)     ---\n",
                        (int) N, gpu_perf, gpu_time );
            }
            
            TESTING_FREE_CPU( ipiv  );
            TESTING_FREE_CPU( work  );
            TESTING_FREE_CPU( h_A   );
            
            TESTING_FREE_PIN( h_R   );
            
            TESTING_FREE_DEV( d_A   );
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
