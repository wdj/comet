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


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zprint
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    magma_tally2DoubleComplex *hA, *dA;
    //magma_tally2_int_t ione     = 1;
    //magma_tally2_int_t ISEED[4] = {0,0,0,1};
    magma_tally2_int_t M, N, lda, ldda;  //size
    magma_tally2_int_t status = 0;
    
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );

    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M     = opts.msize[itest];
            N     = opts.nsize[itest];
            lda   = M;
            ldda  = ((M + 31)/32)*32;
            //size  = lda*N;

            /* Allocate host memory for the matrix */
            TESTING_MALLOC_CPU( hA, magma_tally2DoubleComplex, lda *N );
            TESTING_MALLOC_DEV( dA, magma_tally2DoubleComplex, ldda*N );
        
            //lapackf77_zlarnv( &ione, ISEED, &size, hA );
            for( int j = 0; j < N; ++j ) {
                for( int i = 0; i < M; ++i ) {
                    hA[i + j*lda] = MAGMA_tally2_Z_MAKE( i + j*0.01, 0. );
                }
            }
            magma_tally2_zsetmatrix( M, N, hA, lda, dA, ldda );
            
            printf( "A=" );
            magma_tally2_zprint( M, N, hA, lda );
            
            printf( "dA=" );
            magma_tally2_zprint_gpu( M, N, dA, ldda );
            
            TESTING_FREE_CPU( hA );
            TESTING_FREE_DEV( dA );
        }
    }

    TESTING_FINALIZE();
    return status;
}
