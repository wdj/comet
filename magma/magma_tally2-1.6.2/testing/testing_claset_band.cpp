/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zlaset_band.cpp normal z -> c, Fri Jan 30 19:00:24 2015
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
   -- Testing claset_band
   Code is very similar to testing_clacpy.cpp
*/
int main( int argc, char** argv)
{
    TESTING_INIT();
    
    #define h_A(i_,j_) (h_A + (i_) + (j_)*lda)
    #define d_A(i_,j_) (d_A + (i_) + (j_)*ldda)

    real_Double_t    gbytes, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float           error, work[1];
    magma_tally2FloatComplex  c_neg_one = MAGMA_tally2_C_NEG_ONE;
    magma_tally2FloatComplex *h_A, *h_R;
    magma_tally2FloatComplex_ptr d_A;
    magma_tally2FloatComplex offdiag = MAGMA_tally2_C_MAKE( 1.2000, 6.7000 );
    magma_tally2FloatComplex diag    = MAGMA_tally2_C_MAKE( 3.1415, 2.7183 );
    magma_tally2_int_t M, N, nb, cnt, size, lda, ldda;
    magma_tally2_int_t ione     = 1;
    magma_tally2_int_t status = 0;
    
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );
    
    nb = (opts.nb == 0 ? 32 : opts.nb);

    magma_tally2_uplo_t uplo[] = { Magma_tally2Lower, Magma_tally2Upper, Magma_tally2Full };
    
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
            
            TESTING_MALLOC_CPU( h_A, magma_tally2FloatComplex, size   );
            TESTING_MALLOC_CPU( h_R, magma_tally2FloatComplex, size   );
            
            TESTING_MALLOC_DEV( d_A, magma_tally2FloatComplex, ldda*N );
            
            /* Initialize the matrix */
            for( int j = 0; j < N; ++j ) {
                for( int i = 0; i < M; ++i ) {
                    h_A[i + j*lda] = MAGMA_tally2_C_MAKE( i + j/10000., j );
                }
            }
            magma_tally2_csetmatrix( M, N, h_A, lda, d_A, ldda );
            
            /* =====================================================================
               Performs operation on CPU
               Also count number of elements touched.
               =================================================================== */
            cpu_time = magma_tally2_wtime();
            
            cnt = 0;
            for( int j=inset; j < N-inset; ++j ) {
                for( int k=0; k < nb; ++k ) {  // set k-th sub- or super-diagonal
                    if ( k == 0 && j < M-inset ) {
                        *h_A(j,j)   = diag;
                        cnt += 1;
                    }
                    else if ( uplo[iuplo] == Magma_tally2Lower && j+k < M-inset ) {
                        *h_A(j+k,j) = offdiag;
                        cnt += 1;
                    }
                    else if ( uplo[iuplo] == Magma_tally2Upper && j-k >= inset && j-k < M-inset ) {
                        *h_A(j-k,j) = offdiag;
                        cnt += 1;
                    }
                }
            }
            
            gbytes = cnt / 1e9;
            
            cpu_time = magma_tally2_wtime() - cpu_time;
            cpu_perf = gbytes / cpu_time;
            
            /* ====================================================================
               Performs operation using MAGMA_tally2
               =================================================================== */
            gpu_time = magma_tally2_sync_wtime( 0 );
            
            int mm = M - 2*inset;
            int nn = N - 2*inset;
            magma_tally2blas_claset_band( uplo[iuplo], mm, nn, nb, offdiag, diag, d_A(inset,inset), ldda );
            
            gpu_time = magma_tally2_sync_wtime( 0 ) - gpu_time;
            gpu_perf = gbytes / gpu_time;
            
            /* =====================================================================
               Check the result
               =================================================================== */
            magma_tally2_cgetmatrix( M, N, d_A, ldda, h_R, lda );
                        
            //printf( "h_R=" );  magma_tally2_cprint( M, N, h_R, lda );
            //printf( "h_A=" );  magma_tally2_cprint( M, N, h_A, lda );

            blasf77_caxpy(&size, &c_neg_one, h_A, &ione, h_R, &ione);
            error = lapackf77_clange("f", &M, &N, h_R, &lda, work);
            
            printf("%4c   %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %s\n",
                   lapacke_uplo_const_tally2( uplo[iuplo] ), (int) M, (int) N,
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
