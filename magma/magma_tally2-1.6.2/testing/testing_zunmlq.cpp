/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @precisions normal z -> c d s
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// includes, project
#include "flops.h"
#include "magma_tally2.h"
#include "magma_tally2_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zunmlq
*/
int main( int argc, char** argv )
{
    TESTING_INIT();
    
    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double error, work[1];
    magma_tally2DoubleComplex c_neg_one = MAGMA_tally2_Z_NEG_ONE;
    magma_tally2_int_t ione = 1;
    magma_tally2_int_t mm, m, n, k, size, info;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    magma_tally2_int_t nb, ldc, lda, lwork, lwork_max;
    magma_tally2DoubleComplex *C, *R, *A, *W, *tau;
    magma_tally2_int_t status = 0;
    
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );
    
    // need slightly looser bound (60*eps instead of 30*eps) for some tests
    opts.tolerance = max( 60., opts.tolerance );
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    // test all combinations of input parameters
    magma_tally2_side_t  side [] = { Magma_tally2Left,       Magma_tally2Right   };
    magma_tally2_trans_t trans[] = { Magma_tally2_ConjTrans, Magma_tally2NoTrans };

    printf("    M     N     K   side   trans   CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R||_F / ||QC||_F\n");
    printf("===============================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
      for( int iside = 0; iside < 2; ++iside ) {
      for( int itran = 0; itran < 2; ++itran ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            m = opts.msize[itest];
            n = opts.nsize[itest];
            k = opts.ksize[itest];
            nb  = magma_tally2_get_zgelqf_nb( min( m, n ));
            ldc = m;
            // A is k x m (left) or k x n (right)
            mm = (side[iside] == Magma_tally2Left ? m : n);
            lda = k;
            gflops = FLOPS_ZUNMLQ( m, n, k, side[iside] ) / 1e9;
            
            if ( side[iside] == Magma_tally2Left && m < k ) {
                printf( "%5d %5d %5d   %4c   %5c   skipping because side=left  and m < k\n",
                        (int) m, (int) n, (int) k,
                        lapacke_side_const_tally2( side[iside] ),
                        lapacke_trans_const_tally2( trans[itran] ) );
                continue;
            }
            if ( side[iside] == Magma_tally2Right && n < k ) {
                printf( "%5d %5d %5d   %4c   %5c   skipping because side=right and n < k\n",
                        (int) m, (int) n, (int) k,
                        lapacke_side_const_tally2( side[iside] ),
                        lapacke_trans_const_tally2( trans[itran] ) );
                continue;
            }
            
            // need at least 2*nb*nb for gelqf
            lwork_max = max( max( m*nb, n*nb ), 2*nb*nb );
            
            TESTING_MALLOC_CPU( C,   magma_tally2DoubleComplex, ldc*n );
            TESTING_MALLOC_CPU( R,   magma_tally2DoubleComplex, ldc*n );
            TESTING_MALLOC_CPU( A,   magma_tally2DoubleComplex, lda*mm );
            TESTING_MALLOC_CPU( W,   magma_tally2DoubleComplex, lwork_max );
            TESTING_MALLOC_CPU( tau, magma_tally2DoubleComplex, k );
            
            // C is full, m x n
            size = ldc*n;
            lapackf77_zlarnv( &ione, ISEED, &size, C );
            lapackf77_zlacpy( "Full", &m, &n, C, &ldc, R, &ldc );
            
            size = lda*mm;
            lapackf77_zlarnv( &ione, ISEED, &size, A );
            
            // compute LQ factorization to get Householder vectors in A, tau
            magma_tally2_zgelqf( k, mm, A, lda, tau, W, lwork_max, &info );
            if (info != 0)
                printf("magma_tally2_zgelqf returned error %d: %s.\n",
                       (int) info, magma_tally2_strerror( info ));
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_tally2_wtime();
            lapackf77_zunmlq( lapack_side_const_tally2( side[iside] ), lapack_trans_const_tally2( trans[itran] ),
                              &m, &n, &k,
                              A, &lda, tau, C, &ldc, W, &lwork_max, &info );
            cpu_time = magma_tally2_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            if (info != 0)
                printf("lapackf77_zunmlq returned error %d: %s.\n",
                       (int) info, magma_tally2_strerror( info ));
            
            /* ====================================================================
               Performs operation using MAGMA_tally2
               =================================================================== */
            // query for workspace size
            lwork = -1;
            magma_tally2_zunmlq( side[iside], trans[itran],
                          m, n, k,
                          A, lda, tau, R, ldc, W, lwork, &info );
            if (info != 0)
                printf("magma_tally2_zunmlq (lwork query) returned error %d: %s.\n",
                       (int) info, magma_tally2_strerror( info ));
            lwork = (magma_tally2_int_t) MAGMA_tally2_Z_REAL( W[0] );
            if ( lwork < 0 || lwork > lwork_max ) {
                printf("optimal lwork %d > lwork_max %d\n", (int) lwork, (int) lwork_max );
                lwork = lwork_max;
            }
            
            gpu_time = magma_tally2_wtime();
            magma_tally2_zunmlq( side[iside], trans[itran],
                          m, n, k,
                          A, lda, tau, R, ldc, W, lwork, &info );
            gpu_time = magma_tally2_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally2_zunmlq returned error %d: %s.\n",
                       (int) info, magma_tally2_strerror( info ));
                        
            /* =====================================================================
               compute relative error |QC_magma_tally2 - QC_lapack| / |QC_lapack|
               =================================================================== */
            error = lapackf77_zlange( "Fro", &m, &n, C, &ldc, work );
            size = ldc*n;
            blasf77_zaxpy( &size, &c_neg_one, C, &ione, R, &ione );
            error = lapackf77_zlange( "Fro", &m, &n, R, &ldc, work ) / error;
            
            printf( "%5d %5d %5d   %4c   %5c   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                    (int) m, (int) n, (int) k,
                    lapacke_side_const_tally2( side[iside] ),
                    lapacke_trans_const_tally2( trans[itran] ),
                    cpu_perf, cpu_time, gpu_perf, gpu_time,
                    error, (error < tol ? "ok" : "failed") );
            status += ! (error < tol);
            
            TESTING_FREE_CPU( C );
            TESTING_FREE_CPU( R );
            TESTING_FREE_CPU( A );
            TESTING_FREE_CPU( W );
            TESTING_FREE_CPU( tau );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
      }}  // end iside, itran
      printf( "\n" );
    }
    
    TESTING_FINALIZE();
    return status;
}
