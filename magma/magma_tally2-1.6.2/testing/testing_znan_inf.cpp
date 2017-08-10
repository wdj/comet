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
   -- Testing znan_inf
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    #define hA(i,j) (hA + (i) + (j)*lda)
    
    magma_tally2DoubleComplex *hA, *dA;
    magma_tally2_int_t ione     = 1;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    magma_tally2_int_t M, N, lda, ldda, size;
    magma_tally2_int_t *ii, *jj;
    magma_tally2_int_t i, j, cnt, tmp;
    magma_tally2_int_t status = 0;
    
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );

    magma_tally2_uplo_t uplo[] = { Magma_tally2Lower, Magma_tally2Upper, Magma_tally2Full };
    
    printf("uplo     M     N      CPU nan + inf             GPU nan + inf          actual nan + inf        \n");
    printf("===============================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
      for( int iuplo = 0; iuplo < 3; ++iuplo ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M     = opts.msize[itest];
            N     = opts.nsize[itest];
            lda   = M;
            ldda  = ((M + 31)/32)*32;
            size  = lda*N;

            /* Allocate memory for the matrix */
            TESTING_MALLOC_CPU( hA, magma_tally2DoubleComplex, lda *N );
            TESTING_MALLOC_DEV( dA, magma_tally2DoubleComplex, ldda*N );
            
            /* Initialize the matrix */
            lapackf77_zlarnv( &ione, ISEED, &size, hA );
            
            // up to half of matrix is NAN, and
            // up to half of matrix is INF.
            magma_tally2_int_t cnt_nan = (magma_tally2_int_t)( (rand() / ((double)RAND_MAX)) * 0.5 * M*N );
            magma_tally2_int_t cnt_inf = (magma_tally2_int_t)( (rand() / ((double)RAND_MAX)) * 0.5 * M*N );
            magma_tally2_int_t total = cnt_nan + cnt_inf;
            assert( cnt_nan >= 0 );
            assert( cnt_inf >= 0 );
            assert( total <= M*N );
            
            // fill in indices
            TESTING_MALLOC_CPU( ii, magma_tally2_int_t, size );
            TESTING_MALLOC_CPU( jj, magma_tally2_int_t, size );
            for( cnt=0; cnt < size; ++cnt ) {
                ii[cnt] = cnt % M;
                jj[cnt] = cnt / M;
            }
            // shuffle indices
            for( cnt=0; cnt < total; ++cnt ) {
                i = int( rand() / ((double)RAND_MAX) * size );
                tmp=ii[cnt];  ii[cnt]=ii[i];  ii[i]=tmp;
                tmp=jj[cnt];  jj[cnt]=jj[i];  jj[i]=tmp;
            }
            // fill in NAN and INF
            // for uplo, count NAN and INF in triangular portion of A
            int c_nan=0;
            int c_inf=0;
            for( cnt=0; cnt < cnt_nan; ++cnt ) {
                i = ii[cnt];
                j = jj[cnt];
                *hA(i,j) = MAGMA_tally2_Z_NAN;
                if ( uplo[iuplo] == Magma_tally2Lower && i >= j ) { c_nan++; }
                if ( uplo[iuplo] == Magma_tally2Upper && i <= j ) { c_nan++; }
            }
            for( cnt=cnt_nan; cnt < cnt_nan + cnt_inf; ++cnt ) {
                i = ii[cnt];
                j = jj[cnt];
                *hA(i,j) = MAGMA_tally2_Z_INF;
                if ( uplo[iuplo] == Magma_tally2Lower && i >= j ) { c_inf++; }
                if ( uplo[iuplo] == Magma_tally2Upper && i <= j ) { c_inf++; }
            }
            if ( uplo[iuplo] == Magma_tally2Lower || uplo[iuplo] == Magma_tally2Upper ) {
                cnt_nan = c_nan;
                cnt_inf = c_inf;
                total = cnt_nan + cnt_inf;
            }
            
            //printf( "nan %g + %gi\n", MAGMA_tally2_Z_REAL( MAGMA_tally2_Z_NAN ), MAGMA_tally2_Z_REAL( MAGMA_tally2_Z_NAN ) );
            //printf( "inf %g + %gi\n", MAGMA_tally2_Z_REAL( MAGMA_tally2_Z_INF ), MAGMA_tally2_Z_REAL( MAGMA_tally2_Z_INF ) );
            //magma_tally2_zprint( M, N, hA, lda );
            
            magma_tally2_zsetmatrix( M, N, hA, lda, dA, ldda );
                        
            /* ====================================================================
               Performs operation using MAGMA_tally2
               =================================================================== */
            magma_tally2_int_t c_cpu_nan=-1, c_cpu_inf=-1;
            magma_tally2_int_t c_gpu_nan=-1, c_gpu_inf=-1;
            
            magma_tally2_int_t c_cpu = magma_tally2_znan_inf    ( uplo[iuplo], M, N, hA, lda,  &c_cpu_nan, &c_cpu_inf );
            magma_tally2_int_t c_gpu = magma_tally2_znan_inf_gpu( uplo[iuplo], M, N, dA, ldda, &c_gpu_nan, &c_gpu_inf );
            
            magma_tally2_int_t c_cpu2 = magma_tally2_znan_inf    ( uplo[iuplo], M, N, hA, lda,  NULL, NULL );
            magma_tally2_int_t c_gpu2 = magma_tally2_znan_inf_gpu( uplo[iuplo], M, N, dA, ldda, NULL, NULL );
            
            /* =====================================================================
               Check the result
               =================================================================== */
            bool ok = ( c_cpu == c_gpu )
                   && ( c_cpu == c_cpu2 )
                   && ( c_gpu == c_gpu2 )
                   && ( c_cpu == c_cpu_nan + c_cpu_inf )
                   && ( c_gpu == c_gpu_nan + c_gpu_inf )
                   && ( c_cpu_nan == cnt_nan )
                   && ( c_cpu_inf == cnt_inf )
                   && ( c_gpu_nan == cnt_nan )
                   && ( c_gpu_inf == cnt_inf );
            
            printf( "%4c %5d %5d   %10d + %-10d   %10d + %-10d   %10d + %-10d  %s\n",
                    lapacke_uplo_const_tally2( uplo[iuplo] ), (int) M, (int) N,
                    (int) c_cpu_nan, (int) c_cpu_inf,
                    (int) c_gpu_nan, (int) c_gpu_inf,
                    (int) cnt_nan,   (int) cnt_inf,
                    (ok ? "ok" : "failed"));
            status += ! ok;
            
            TESTING_FREE_CPU( hA );
            TESTING_FREE_DEV( dA );
            
            TESTING_FREE_CPU( ii );
            TESTING_FREE_CPU( jj );
        }
      }
      printf( "\n" );
    }

    TESTING_FINALIZE();
    return status;
}