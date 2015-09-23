/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
  
       @generated from testing_znan_inf.cpp normal z -> d, Fri Jan 30 19:00:24 2015
       @author Mark Gates
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "testings.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing znan_inf
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    #define hA(i,j) (hA + (i) + (j)*lda)
    
    double *hA, *dA;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t M, N, lda, ldda, size;
    magma_minproduct_int_t *ii, *jj;
    magma_minproduct_int_t i, j, cnt, tmp;
    magma_minproduct_int_t status = 0;
    
    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );

    magma_minproduct_uplo_t uplo[] = { Magma_minproductLower, Magma_minproductUpper, Magma_minproductFull };
    
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
            TESTING_MALLOC_CPU( hA, double, lda *N );
            TESTING_MALLOC_DEV( dA, double, ldda*N );
            
            /* Initialize the matrix */
            lapackf77_dlarnv( &ione, ISEED, &size, hA );
            
            // up to half of matrix is NAN, and
            // up to half of matrix is INF.
            magma_minproduct_int_t cnt_nan = (magma_minproduct_int_t)( (rand() / ((double)RAND_MAX)) * 0.5 * M*N );
            magma_minproduct_int_t cnt_inf = (magma_minproduct_int_t)( (rand() / ((double)RAND_MAX)) * 0.5 * M*N );
            magma_minproduct_int_t total = cnt_nan + cnt_inf;
            assert( cnt_nan >= 0 );
            assert( cnt_inf >= 0 );
            assert( total <= M*N );
            
            // fill in indices
            TESTING_MALLOC_CPU( ii, magma_minproduct_int_t, size );
            TESTING_MALLOC_CPU( jj, magma_minproduct_int_t, size );
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
                *hA(i,j) = MAGMA_minproduct_D_NAN;
                if ( uplo[iuplo] == Magma_minproductLower && i >= j ) { c_nan++; }
                if ( uplo[iuplo] == Magma_minproductUpper && i <= j ) { c_nan++; }
            }
            for( cnt=cnt_nan; cnt < cnt_nan + cnt_inf; ++cnt ) {
                i = ii[cnt];
                j = jj[cnt];
                *hA(i,j) = MAGMA_minproduct_D_INF;
                if ( uplo[iuplo] == Magma_minproductLower && i >= j ) { c_inf++; }
                if ( uplo[iuplo] == Magma_minproductUpper && i <= j ) { c_inf++; }
            }
            if ( uplo[iuplo] == Magma_minproductLower || uplo[iuplo] == Magma_minproductUpper ) {
                cnt_nan = c_nan;
                cnt_inf = c_inf;
                total = cnt_nan + cnt_inf;
            }
            
            //printf( "nan %g + %gi\n", MAGMA_minproduct_D_REAL( MAGMA_minproduct_D_NAN ), MAGMA_minproduct_D_REAL( MAGMA_minproduct_D_NAN ) );
            //printf( "inf %g + %gi\n", MAGMA_minproduct_D_REAL( MAGMA_minproduct_D_INF ), MAGMA_minproduct_D_REAL( MAGMA_minproduct_D_INF ) );
            //magma_minproduct_dprint( M, N, hA, lda );
            
            magma_minproduct_dsetmatrix( M, N, hA, lda, dA, ldda );
                        
            /* ====================================================================
               Performs operation using MAGMA_minproduct
               =================================================================== */
            magma_minproduct_int_t c_cpu_nan=-1, c_cpu_inf=-1;
            magma_minproduct_int_t c_gpu_nan=-1, c_gpu_inf=-1;
            
            magma_minproduct_int_t c_cpu = magma_minproduct_dnan_inf    ( uplo[iuplo], M, N, hA, lda,  &c_cpu_nan, &c_cpu_inf );
            magma_minproduct_int_t c_gpu = magma_minproduct_dnan_inf_gpu( uplo[iuplo], M, N, dA, ldda, &c_gpu_nan, &c_gpu_inf );
            
            magma_minproduct_int_t c_cpu2 = magma_minproduct_dnan_inf    ( uplo[iuplo], M, N, hA, lda,  NULL, NULL );
            magma_minproduct_int_t c_gpu2 = magma_minproduct_dnan_inf_gpu( uplo[iuplo], M, N, dA, ldda, NULL, NULL );
            
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
                    lapacke_uplo_const( uplo[iuplo] ), (int) M, (int) N,
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
