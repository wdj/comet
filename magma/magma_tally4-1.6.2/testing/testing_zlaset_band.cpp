/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
       @author Mark Gates
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zlaset_band
   Code is very similar to testing_zlacpy.cpp
*/
int main( int argc, char** argv)
{
    TESTING_INIT();
    
    #define h_A(i_,j_) (h_A + (i_) + (j_)*lda)
    #define d_A(i_,j_) (d_A + (i_) + (j_)*ldda)

    real_Double_t    gbytes, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double           error, work[1];
    magma_tally4DoubleComplex  c_neg_one = MAGMA_tally4_Z_NEG_ONE;
    magma_tally4DoubleComplex *h_A, *h_R;
    magma_tally4DoubleComplex_ptr d_A;
    magma_tally4DoubleComplex offdiag = MAGMA_tally4_Z_MAKE( 1.2000, 6.7000 );
    magma_tally4DoubleComplex diag    = MAGMA_tally4_Z_MAKE( 3.1415, 2.7183 );
    magma_tally4_int_t M, N, nb, cnt, size, lda, ldda;
    magma_tally4_int_t ione     = 1;
    magma_tally4_int_t status = 0;
    
    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );
    
    nb = (opts.nb == 0 ? 32 : opts.nb);

    magma_tally4_uplo_t uplo[] = { Magma_tally4Lower, Magma_tally4Upper, Magma_tally4Full };
    
    printf("K = nb = %d\n", (int) nb );
    printf("uplo       M     N   CPU GByte/s (ms)    GPU GByte/s (ms)    check\n");
    printf("==================================================================\n");
    for( int iuplo = 0; iuplo < 2; ++iuplo ) {
      for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            int inset = 0;
            M = opts.msize[itest] + 2*inset;
            N = opts.nsize[itest] + 2*inset;
            lda    = M;
            ldda   = ((M+31)/32)*32;
            size   = lda*N;
            
            TESTING_MALLOC_CPU( h_A, magma_tally4DoubleComplex, size   );
            TESTING_MALLOC_CPU( h_R, magma_tally4DoubleComplex, size   );
            
            TESTING_MALLOC_DEV( d_A, magma_tally4DoubleComplex, ldda*N );
            
            /* Initialize the matrix */
            for( int j = 0; j < N; ++j ) {
                for( int i = 0; i < M; ++i ) {
                    h_A[i + j*lda] = MAGMA_tally4_Z_MAKE( i + j/10000., j );
                }
            }
            magma_tally4_zsetmatrix( M, N, h_A, lda, d_A, ldda );
            
            /* =====================================================================
               Performs operation on CPU
               Also count number of elements touched.
               =================================================================== */
            cpu_time = magma_tally4_wtime();
            
            cnt = 0;
            for( int j=inset; j < N-inset; ++j ) {
                for( int k=0; k < nb; ++k ) {  // set k-th sub- or super-diagonal
                    if ( k == 0 && j < M-inset ) {
                        *h_A(j,j)   = diag;
                        cnt += 1;
                    }
                    else if ( uplo[iuplo] == Magma_tally4Lower && j+k < M-inset ) {
                        *h_A(j+k,j) = offdiag;
                        cnt += 1;
                    }
                    else if ( uplo[iuplo] == Magma_tally4Upper && j-k >= inset && j-k < M-inset ) {
                        *h_A(j-k,j) = offdiag;
                        cnt += 1;
                    }
                }
            }
            
            gbytes = cnt / 1e9;
            
            cpu_time = magma_tally4_wtime() - cpu_time;
            cpu_perf = gbytes / cpu_time;
            
            /* ====================================================================
               Performs operation using MAGMA_tally4
               =================================================================== */
            gpu_time = magma_tally4_sync_wtime( 0 );
            
            int mm = M - 2*inset;
            int nn = N - 2*inset;
            magma_tally4blas_zlaset_band( uplo[iuplo], mm, nn, nb, offdiag, diag, d_A(inset,inset), ldda );
            
            gpu_time = magma_tally4_sync_wtime( 0 ) - gpu_time;
            gpu_perf = gbytes / gpu_time;
            
            /* =====================================================================
               Check the result
               =================================================================== */
            magma_tally4_zgetmatrix( M, N, d_A, ldda, h_R, lda );
                        
            //printf( "h_R=" );  magma_tally4_zprint( M, N, h_R, lda );
            //printf( "h_A=" );  magma_tally4_zprint( M, N, h_A, lda );

            blasf77_zaxpy(&size, &c_neg_one, h_A, &ione, h_R, &ione);
            error = lapackf77_zlange("f", &M, &N, h_R, &lda, work);
            
            printf("%4c   %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %s\n",
                   lapacke_uplo_const_tally4( uplo[iuplo] ), (int) M, (int) N,
                   cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000.,
                   (error == 0. ? "ok" : "failed") );
            status += ! (error == 0.);
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_R );
            
            TESTING_FREE_DEV( d_A );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
      }
      printf( "\n" );
    }

    TESTING_FINALIZE();
    return status;
}
