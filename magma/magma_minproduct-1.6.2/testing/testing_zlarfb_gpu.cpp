/*
    -- MAGMA_minproduct (version 1.6.1) --
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

#include <algorithm>  // std::swap

// includes, project
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zlarfb_gpu
*/
int main( int argc, char** argv )
{
    TESTING_INIT();
    
    magma_minproductDoubleComplex c_zero    = MAGMA_minproduct_Z_ZERO;
    magma_minproductDoubleComplex c_one     = MAGMA_minproduct_Z_ONE;
    magma_minproductDoubleComplex c_neg_one = MAGMA_minproduct_Z_NEG_ONE;
    magma_minproduct_int_t M, N, K, size, ldc, ldv, ldt, ldw, nv;
    magma_minproduct_int_t ione =  1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    double error, work[1];
    magma_minproduct_int_t status = 0;
    
    // test all combinations of input parameters
    magma_minproduct_side_t   side  [] = { Magma_minproductLeft,       Magma_minproductRight    };
    magma_minproduct_trans_t  trans [] = { Magma_minproductConjTrans,  Magma_minproductNoTrans  };
    magma_minproduct_direct_t direct[] = { Magma_minproductForward,    Magma_minproductBackward };
    magma_minproduct_storev_t storev[] = { Magma_minproductColumnwise, Magma_minproductRowwise  };

    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("    M     N     K   storev   side   direct   trans    ||R||_F / ||HC||_F\n");
    printf("========================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
      M = opts.msize[itest];
      N = opts.nsize[itest];
      K = opts.ksize[itest];
      if ( M < K || N < K || K <= 0 ) {
          printf( "%5d %5d %5d   skipping because zlarfb requires M >= K, N >= K, K >= 0\n",
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
            ldw = (side[iside] == Magma_minproductLeft ? N : M);
            // (ldv, nv) get swapped later if rowwise
            ldv = (side[iside] == Magma_minproductLeft ? M : N);
            nv  = K;
            
            // Allocate memory for matrices
            magma_minproductDoubleComplex *C, *R, *V, *T, *W;
            TESTING_MALLOC_CPU( C, magma_minproductDoubleComplex, ldc*N );
            TESTING_MALLOC_CPU( R, magma_minproductDoubleComplex, ldc*N );
            TESTING_MALLOC_CPU( V, magma_minproductDoubleComplex, ldv*K );
            TESTING_MALLOC_CPU( T, magma_minproductDoubleComplex, ldt*K );
            TESTING_MALLOC_CPU( W, magma_minproductDoubleComplex, ldw*K );
            
            magma_minproductDoubleComplex_ptr dC, dV, dT, dW;
            TESTING_MALLOC_DEV( dC, magma_minproductDoubleComplex, ldc*N );
            TESTING_MALLOC_DEV( dV, magma_minproductDoubleComplex, ldv*K );
            TESTING_MALLOC_DEV( dT, magma_minproductDoubleComplex, ldt*K );
            TESTING_MALLOC_DEV( dW, magma_minproductDoubleComplex, ldw*K );
            
            // C is M x N.
            size = ldc*N;
            lapackf77_zlarnv( &ione, ISEED, &size, C );
            //printf( "C=" );  magma_minproduct_zprint( M, N, C, ldc );
            
            // V is ldv x nv. See larfb docs for description.
            // if column-wise and left,  M x K
            // if column-wise and right, N x K
            // if row-wise and left,     K x M
            // if row-wise and right,    K x N
            size = ldv*nv;
            lapackf77_zlarnv( &ione, ISEED, &size, V );
            if ( storev[istor] == Magma_minproductColumnwise ) {
                if ( direct[idir] == Magma_minproductForward ) {
                    lapackf77_zlaset( Magma_minproductUpperStr, &K, &K, &c_zero, &c_one, V, &ldv );
                }
                else {
                    lapackf77_zlaset( Magma_minproductLowerStr, &K, &K, &c_zero, &c_one, &V[(ldv-K)], &ldv );
                }
            }
            else {
                // rowwise, swap V's dimensions
                std::swap( ldv, nv );
                if ( direct[idir] == Magma_minproductForward ) {
                    lapackf77_zlaset( Magma_minproductLowerStr, &K, &K, &c_zero, &c_one, V, &ldv );
                }
                else {
                    lapackf77_zlaset( Magma_minproductUpperStr, &K, &K, &c_zero, &c_one, &V[(nv-K)*ldv], &ldv );
                }
            }
            //printf( "# ldv %d, nv %d\n", ldv, nv );
            //printf( "V=" );  magma_minproduct_zprint( ldv, nv, V, ldv );
            
            // T is K x K, upper triangular for forward, and lower triangular for backward
            magma_minproduct_int_t k1 = K-1;
            size = ldt*K;
            lapackf77_zlarnv( &ione, ISEED, &size, T );
            if ( direct[idir] == Magma_minproductForward ) {
                lapackf77_zlaset( Magma_minproductLowerStr, &k1, &k1, &c_zero, &c_zero, &T[1], &ldt );
            }
            else {
                lapackf77_zlaset( Magma_minproductUpperStr, &k1, &k1, &c_zero, &c_zero, &T[1*ldt], &ldt );
            }
            //printf( "T=" );  magma_minproduct_zprint( K, K, T, ldt );
            
            magma_minproduct_zsetmatrix( M,   N,  C, ldc, dC, ldc );
            magma_minproduct_zsetmatrix( ldv, nv, V, ldv, dV, ldv );
            magma_minproduct_zsetmatrix( K,   K,  T, ldt, dT, ldt );
            
            lapackf77_zlarfb( lapack_side_const( side[iside] ), lapack_trans_const( trans[itran] ),
                              lapack_direct_const( direct[idir] ), lapack_storev_const( storev[istor] ),
                              &M, &N, &K,
                              V, &ldv, T, &ldt, C, &ldc, W, &ldw );
            //printf( "HC=" );  magma_minproduct_zprint( M, N, C, ldc );
            
            magma_minproduct_zlarfb_gpu( side[iside], trans[itran], direct[idir], storev[istor],
                              M, N, K,
                              dV, ldv, dT, ldt, dC, ldc, dW, ldw );
            magma_minproduct_zgetmatrix( M, N, dC, ldc, R, ldc );
            //printf( "dHC=" );  magma_minproduct_zprint( M, N, R, ldc );
            
            // compute relative error |HC_magma_minproduct - HC_lapack| / |HC_lapack|
            error = lapackf77_zlange( "Fro", &M, &N, C, &ldc, work );
            size = ldc*N;
            blasf77_zaxpy( &size, &c_neg_one, C, &ione, R, &ione );
            error = lapackf77_zlange( "Fro", &M, &N, R, &ldc, work ) / error;
            printf( "%5d %5d %5d      %c       %c       %c       %c      %8.2e   %s\n",
                    (int) M, (int) N, (int) K,
                    lapacke_storev_const(storev[istor]), lapacke_side_const(side[iside]),
                    lapacke_direct_const(direct[idir]), lapacke_trans_const(trans[itran]),
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
