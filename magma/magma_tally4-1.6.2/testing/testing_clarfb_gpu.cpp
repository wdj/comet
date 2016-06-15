/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @generated from testing_zlarfb_gpu.cpp normal z -> c, Fri Jan 30 19:00:25 2015
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <algorithm>  // std::swap

// includes, project
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing clarfb_gpu
*/
int main( int argc, char** argv )
{
    TESTING_INIT();
    
    magma_tally4FloatComplex c_zero    = MAGMA_tally4_C_ZERO;
    magma_tally4FloatComplex c_one     = MAGMA_tally4_C_ONE;
    magma_tally4FloatComplex c_neg_one = MAGMA_tally4_C_NEG_ONE;
    magma_tally4_int_t M, N, K, size, ldc, ldv, ldt, ldw, nv;
    magma_tally4_int_t ione =  1;
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    float error, work[1];
    magma_tally4_int_t status = 0;
    
    // test all combinations of input parameters
    magma_tally4_side_t   side  [] = { Magma_tally4Left,       Magma_tally4Right    };
    magma_tally4_trans_t  trans [] = { Magma_tally4ConjTrans,  Magma_tally4NoTrans  };
    magma_tally4_direct_t direct[] = { Magma_tally4Forward,    Magma_tally4Backward };
    magma_tally4_storev_t storev[] = { Magma_tally4Columnwise, Magma_tally4Rowwise  };

    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");
    
    printf("    M     N     K   storev   side   direct   trans    ||R||_F / ||HC||_F\n");
    printf("========================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
      M = opts.msize[itest];
      N = opts.nsize[itest];
      K = opts.ksize[itest];
      if ( M < K || N < K || K <= 0 ) {
          printf( "%5d %5d %5d   skipping because clarfb requires M >= K, N >= K, K >= 0\n",
                  (int) M, (int) N, (int) K );
          continue;
      }
      for( int istor = 0; istor < 2; ++istor ) {
      for( int iside = 0; iside < 2; ++iside ) {
      for( int idir  = 0; idir  < 2; ++idir  ) {
      for( int itran = 0; itran < 2; ++itran ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {            
            ldc = ((M+31)/32)*32;
            ldt = ((K+31)/32)*32;
            ldw = (side[iside] == Magma_tally4Left ? N : M);
            // (ldv, nv) get swapped later if rowwise
            ldv = (side[iside] == Magma_tally4Left ? M : N);
            nv  = K;
            
            // Allocate memory for matrices
            magma_tally4FloatComplex *C, *R, *V, *T, *W;
            TESTING_MALLOC_CPU( C, magma_tally4FloatComplex, ldc*N );
            TESTING_MALLOC_CPU( R, magma_tally4FloatComplex, ldc*N );
            TESTING_MALLOC_CPU( V, magma_tally4FloatComplex, ldv*K );
            TESTING_MALLOC_CPU( T, magma_tally4FloatComplex, ldt*K );
            TESTING_MALLOC_CPU( W, magma_tally4FloatComplex, ldw*K );
            
            magma_tally4FloatComplex_ptr dC, dV, dT, dW;
            TESTING_MALLOC_DEV( dC, magma_tally4FloatComplex, ldc*N );
            TESTING_MALLOC_DEV( dV, magma_tally4FloatComplex, ldv*K );
            TESTING_MALLOC_DEV( dT, magma_tally4FloatComplex, ldt*K );
            TESTING_MALLOC_DEV( dW, magma_tally4FloatComplex, ldw*K );
            
            // C is M x N.
            size = ldc*N;
            lapackf77_clarnv( &ione, ISEED, &size, C );
            //printf( "C=" );  magma_tally4_cprint( M, N, C, ldc );
            
            // V is ldv x nv. See larfb docs for description.
            // if column-wise and left,  M x K
            // if column-wise and right, N x K
            // if row-wise and left,     K x M
            // if row-wise and right,    K x N
            size = ldv*nv;
            lapackf77_clarnv( &ione, ISEED, &size, V );
            if ( storev[istor] == Magma_tally4Columnwise ) {
                if ( direct[idir] == Magma_tally4Forward ) {
                    lapackf77_claset( Magma_tally4UpperStr, &K, &K, &c_zero, &c_one, V, &ldv );
                }
                else {
                    lapackf77_claset( Magma_tally4LowerStr, &K, &K, &c_zero, &c_one, &V[(ldv-K)], &ldv );
                }
            }
            else {
                // rowwise, swap V's dimensions
                std::swap( ldv, nv );
                if ( direct[idir] == Magma_tally4Forward ) {
                    lapackf77_claset( Magma_tally4LowerStr, &K, &K, &c_zero, &c_one, V, &ldv );
                }
                else {
                    lapackf77_claset( Magma_tally4UpperStr, &K, &K, &c_zero, &c_one, &V[(nv-K)*ldv], &ldv );
                }
            }
            //printf( "# ldv %d, nv %d\n", ldv, nv );
            //printf( "V=" );  magma_tally4_cprint( ldv, nv, V, ldv );
            
            // T is K x K, upper triangular for forward, and lower triangular for backward
            magma_tally4_int_t k1 = K-1;
            size = ldt*K;
            lapackf77_clarnv( &ione, ISEED, &size, T );
            if ( direct[idir] == Magma_tally4Forward ) {
                lapackf77_claset( Magma_tally4LowerStr, &k1, &k1, &c_zero, &c_zero, &T[1], &ldt );
            }
            else {
                lapackf77_claset( Magma_tally4UpperStr, &k1, &k1, &c_zero, &c_zero, &T[1*ldt], &ldt );
            }
            //printf( "T=" );  magma_tally4_cprint( K, K, T, ldt );
            
            magma_tally4_csetmatrix( M,   N,  C, ldc, dC, ldc );
            magma_tally4_csetmatrix( ldv, nv, V, ldv, dV, ldv );
            magma_tally4_csetmatrix( K,   K,  T, ldt, dT, ldt );
            
            lapackf77_clarfb( lapack_side_const_tally4( side[iside] ), lapack_trans_const_tally4( trans[itran] ),
                              lapack_direct_const_tally4( direct[idir] ), lapack_storev_const_tally4( storev[istor] ),
                              &M, &N, &K,
                              V, &ldv, T, &ldt, C, &ldc, W, &ldw );
            //printf( "HC=" );  magma_tally4_cprint( M, N, C, ldc );
            
            magma_tally4_clarfb_gpu( side[iside], trans[itran], direct[idir], storev[istor],
                              M, N, K,
                              dV, ldv, dT, ldt, dC, ldc, dW, ldw );
            magma_tally4_cgetmatrix( M, N, dC, ldc, R, ldc );
            //printf( "dHC=" );  magma_tally4_cprint( M, N, R, ldc );
            
            // compute relative error |HC_magma_tally4 - HC_lapack| / |HC_lapack|
            error = lapackf77_clange( "Fro", &M, &N, C, &ldc, work );
            size = ldc*N;
            blasf77_caxpy( &size, &c_neg_one, C, &ione, R, &ione );
            error = lapackf77_clange( "Fro", &M, &N, R, &ldc, work ) / error;
            printf( "%5d %5d %5d      %c       %c       %c       %c      %8.2e   %s\n",
                    (int) M, (int) N, (int) K,
                    lapacke_storev_const_tally4(storev[istor]), lapacke_side_const_tally4(side[iside]),
                    lapacke_direct_const_tally4(direct[idir]), lapacke_trans_const_tally4(trans[itran]),
                   error, (error < tol ? "ok" : "failed") );
            status += ! (error < tol);
            
            TESTING_FREE_CPU( C );
            TESTING_FREE_CPU( R );
            TESTING_FREE_CPU( V );
            TESTING_FREE_CPU( T );
            TESTING_FREE_CPU( W );
            
            TESTING_FREE_DEV( dC );
            TESTING_FREE_DEV( dV );
            TESTING_FREE_DEV( dT );
            TESTING_FREE_DEV( dW );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
      }}}}
      printf( "\n" );
    }
    
    TESTING_FINALIZE();
    return status;
}
