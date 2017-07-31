/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @generated from testing_zunmqr_gpu.cpp normal z -> c, Fri Jan 30 19:00:25 2015
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
   -- Testing cunmqr_gpu
*/
int main( int argc, char** argv )
{
    TESTING_INIT();
    
    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float error, work[1];
    magma_tally2FloatComplex c_neg_one = MAGMA_tally2_C_NEG_ONE;
    magma_tally2_int_t ione = 1;
    magma_tally2_int_t m, n, k, size, info;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    magma_tally2_int_t nb, ldc, lda, lwork, lwork_max, dt_size;
    magma_tally2FloatComplex *C, *R, *A, *W, *tau;
    magma_tally2FloatComplex_ptr dC, dA, dT;
    magma_tally2_int_t status = 0;
    
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = 2. * opts.tolerance * lapackf77_slamch("E");
    
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
            nb  = magma_tally2_get_cgeqrf_nb( m );
            ldc = ((m + 31)/32)*32;
            lda = ((max(m,n) + 31)/32)*32;
            gflops = FLOPS_CUNMQR( m, n, k, side[iside] ) / 1e9;
            
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
            
            if ( side[iside] == Magma_tally2Left ) {
                // side = left
                lwork_max = (m - k + nb)*(n + nb) + n*nb;
                dt_size = ( 2*min(m,k) + ((max(m,n) + 31)/32)*32 )*nb;
            }
            else {
                // side = right
                lwork_max = (n - k + nb)*(m + nb) + m*nb;
                dt_size = ( 2*min(n,k) + ((max(m,n) + 31)/32)*32 )*nb;
            }
            
            TESTING_MALLOC_CPU( C,   magma_tally2FloatComplex, ldc*n );
            TESTING_MALLOC_CPU( R,   magma_tally2FloatComplex, ldc*n );
            TESTING_MALLOC_CPU( A,   magma_tally2FloatComplex, lda*k );
            TESTING_MALLOC_CPU( W,   magma_tally2FloatComplex, lwork_max );
            TESTING_MALLOC_CPU( tau, magma_tally2FloatComplex, k );
            
            TESTING_MALLOC_DEV( dC, magma_tally2FloatComplex, ldc*n );
            TESTING_MALLOC_DEV( dA, magma_tally2FloatComplex, lda*k );
            TESTING_MALLOC_DEV( dT, magma_tally2FloatComplex, dt_size );
            
            // C is full, m x n
            size = ldc*n;
            lapackf77_clarnv( &ione, ISEED, &size, C );
            magma_tally2_csetmatrix( m, n, C, ldc, dC, ldc );
            
            // A is m x k (left) or n x k (right)
            lda = (side[iside] == Magma_tally2Left ? m : n);
            size = lda*k;
            lapackf77_clarnv( &ione, ISEED, &size, A );
            
            // compute QR factorization to get Householder vectors in dA, tau, dT
            magma_tally2_csetmatrix( lda, k, A,  lda, dA, lda );
            magma_tally2_cgeqrf_gpu( lda, k, dA, lda, tau, dT, &info );
            magma_tally2_cgetmatrix( lda, k, dA, lda, A,  lda );
            if (info != 0)
                printf("magma_tally2_cgeqrf_gpu returned error %d: %s.\n",
                       (int) info, magma_tally2_strerror( info ));
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_tally2_wtime();
            lapackf77_cunmqr( lapack_side_const_tally2( side[iside] ), lapack_trans_const_tally2( trans[itran] ),
                              &m, &n, &k,
                              A, &lda, tau, C, &ldc, W, &lwork_max, &info );
            cpu_time = magma_tally2_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            if (info != 0)
                printf("lapackf77_cunmqr returned error %d: %s.\n",
                       (int) info, magma_tally2_strerror( info ));
            
            /* ====================================================================
               Performs operation using MAGMA_tally2
               =================================================================== */
            // query for workspace size
            lwork = -1;
            magma_tally2_cunmqr_gpu( side[iside], trans[itran],
                              m, n, k,
                              dA, lda, tau, dC, ldc, W, lwork, dT, nb, &info );
            if (info != 0)
                printf("magma_tally2_cunmqr_gpu (lwork query) returned error %d: %s.\n",
                       (int) info, magma_tally2_strerror( info ));
            lwork = (magma_tally2_int_t) MAGMA_tally2_C_REAL( W[0] );
            if ( lwork < 0 || lwork > lwork_max )
                printf("invalid lwork %d, lwork_max %d\n", (int) lwork, (int) lwork_max );
            
            gpu_time = magma_tally2_sync_wtime( 0 );  // sync needed for L,N and R,T cases
            magma_tally2_cunmqr_gpu( side[iside], trans[itran],
                              m, n, k,
                              dA, lda, tau, dC, ldc, W, lwork, dT, nb, &info );
            gpu_time = magma_tally2_sync_wtime( 0 ) - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally2_cunmqr_gpu returned error %d: %s.\n",
                       (int) info, magma_tally2_strerror( info ));
            
            magma_tally2_cgetmatrix( m, n, dC, ldc, R, ldc );
            
            /* =====================================================================
               compute relative error |QC_magma_tally2 - QC_lapack| / |QC_lapack|
               =================================================================== */
            error = lapackf77_clange( "Fro", &m, &n, C, &ldc, work );
            size = ldc*n;
            blasf77_caxpy( &size, &c_neg_one, C, &ione, R, &ione );
            error = lapackf77_clange( "Fro", &m, &n, R, &ldc, work ) / error;
            
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
            
            TESTING_FREE_DEV( dC );
            TESTING_FREE_DEV( dA );
            TESTING_FREE_DEV( dT );
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
