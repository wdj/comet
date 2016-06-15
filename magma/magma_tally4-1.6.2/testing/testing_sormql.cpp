/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @generated from testing_zunmql.cpp normal z -> s, Fri Jan 30 19:00:26 2015
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// includes, project
#include "flops.h"
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sormql
*/
int main( int argc, char** argv )
{
    TESTING_INIT();
    
    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float error, work[1];
    float c_neg_one = MAGMA_tally4_S_NEG_ONE;
    magma_tally4_int_t ione = 1;
    magma_tally4_int_t mm, m, n, k, size, info;
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    magma_tally4_int_t nb, ldc, lda, lwork, lwork_max;
    float *C, *R, *A, *W, *tau;
    magma_tally4_int_t status = 0;
    
    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );
    
    // need slightly looser bound (60*eps instead of 30*eps) for some tests
    opts.tolerance = max( 60., opts.tolerance );
    float tol = opts.tolerance * lapackf77_slamch("E");
    
    // test all combinations of input parameters
    magma_tally4_side_t  side [] = { Magma_tally4Left,       Magma_tally4Right   };
    magma_tally4_trans_t trans[] = { Magma_tally4Trans, Magma_tally4NoTrans };

    printf("    M     N     K   side   trans   CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R||_F / ||QC||_F\n");
    printf("===============================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
      for( int iside = 0; iside < 2; ++iside ) {
      for( int itran = 0; itran < 2; ++itran ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            m = opts.msize[itest];
            n = opts.nsize[itest];
            k = opts.ksize[itest];
            nb  = magma_tally4_get_sgeqlf_nb( m );
            ldc = m;
            // A is m x k (left) or n x k (right)
            mm = (side[iside] == Magma_tally4Left ? m : n);
            lda = mm;
            gflops = FLOPS_SORMQL( m, n, k, side[iside] ) / 1e9;
            
            if ( side[iside] == Magma_tally4Left && m < k ) {
                printf( "%5d %5d %5d   %4c   %5c   skipping because side=left  and m < k\n",
                        (int) m, (int) n, (int) k,
                        lapacke_side_const_tally4( side[iside] ),
                        lapacke_trans_const_tally4( trans[itran] ) );
                continue;
            }
            if ( side[iside] == Magma_tally4Right && n < k ) {
                printf( "%5d %5d %5d  %4c   %5c    skipping because side=right and n < k\n",
                        (int) m, (int) n, (int) k,
                        lapacke_side_const_tally4( side[iside] ),
                        lapacke_trans_const_tally4( trans[itran] ) );
                continue;
            }
            
            // need at least 2*nb*nb for geqlf
            lwork_max = max( max( m*nb, n*nb ), 2*nb*nb );
            
            TESTING_MALLOC_CPU( C,   float, ldc*n );
            TESTING_MALLOC_CPU( R,   float, ldc*n );
            TESTING_MALLOC_CPU( A,   float, lda*k );
            TESTING_MALLOC_CPU( W,   float, lwork_max );
            TESTING_MALLOC_CPU( tau, float, k );
            
            // C is full, m x n
            size = ldc*n;
            lapackf77_slarnv( &ione, ISEED, &size, C );
            lapackf77_slacpy( "Full", &m, &n, C, &ldc, R, &ldc );
            
            size = lda*k;
            lapackf77_slarnv( &ione, ISEED, &size, A );
            
            // compute QL factorization to get Householder vectors in A, tau
            magma_tally4_sgeqlf( mm, k, A, lda, tau, W, lwork_max, &info );
            if (info != 0)
                printf("magma_tally4_sgeqlf returned error %d: %s.\n",
                       (int) info, magma_tally4_strerror( info ));
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_tally4_wtime();
            lapackf77_sormql( lapack_side_const_tally4( side[iside] ), lapack_trans_const_tally4( trans[itran] ),
                              &m, &n, &k,
                              A, &lda, tau, C, &ldc, W, &lwork_max, &info );
            cpu_time = magma_tally4_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            if (info != 0)
                printf("lapackf77_sormql returned error %d: %s.\n",
                       (int) info, magma_tally4_strerror( info ));
            
            /* ====================================================================
               Performs operation using MAGMA_tally4
               =================================================================== */
            // query for workspace size
            lwork = -1;
            magma_tally4_sormql( side[iside], trans[itran],
                          m, n, k,
                          A, lda, tau, R, ldc, W, lwork, &info );
            if (info != 0)
                printf("magma_tally4_sormql (lwork query) returned error %d: %s.\n",
                       (int) info, magma_tally4_strerror( info ));
            lwork = (magma_tally4_int_t) MAGMA_tally4_S_REAL( W[0] );
            if ( lwork < 0 || lwork > lwork_max ) {
                printf("optimal lwork %d > lwork_max %d\n", (int) lwork, (int) lwork_max );
                lwork = lwork_max;
            }
            
            gpu_time = magma_tally4_wtime();
            magma_tally4_sormql( side[iside], trans[itran],
                          m, n, k,
                          A, lda, tau, R, ldc, W, lwork, &info );
            gpu_time = magma_tally4_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally4_sormql returned error %d: %s.\n",
                       (int) info, magma_tally4_strerror( info ));
                        
            /* =====================================================================
               compute relative error |QC_magma_tally4 - QC_lapack| / |QC_lapack|
               =================================================================== */
            error = lapackf77_slange( "Fro", &m, &n, C, &ldc, work );
            size = ldc*n;
            blasf77_saxpy( &size, &c_neg_one, C, &ione, R, &ione );
            error = lapackf77_slange( "Fro", &m, &n, R, &ldc, work ) / error;
            
            printf( "%5d %5d %5d   %4c   %5c   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                    (int) m, (int) n, (int) k,
                    lapacke_side_const_tally4( side[iside] ),
                    lapacke_trans_const_tally4( trans[itran] ),
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
