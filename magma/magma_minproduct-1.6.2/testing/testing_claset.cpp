/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zlaset.cpp normal z -> c, Fri Jan 30 19:00:24 2015
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
#include "magma_minproduct_operators.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing claset
   Code is very similar to testing_clacpy.cpp
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gbytes, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float           error, work[1];
    magma_minproductFloatComplex  c_neg_one = MAGMA_minproduct_C_NEG_ONE;
    magma_minproductFloatComplex *h_A, *h_R;
    magma_minproductFloatComplex_ptr d_A;
    magma_minproductFloatComplex offdiag, diag;
    magma_minproduct_int_t M, N, size, lda, ldda;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t status = 0;
    
    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );

    magma_minproduct_uplo_t uplo[] = { Magma_minproductLower, Magma_minproductUpper, Magma_minproductFull };

    printf("uplo      M     N  offdiag  diag    CPU GByte/s (ms)    GPU GByte/s (ms)   check\n");
    printf("================================================================================\n");
    for( int iuplo = 0; iuplo < 3; ++iuplo ) {
      for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
          for( int ival = 0; ival < 4; ++ival ) {
            // test combinations of zero & non-zero:
            // ival  offdiag  diag
            // 0     0        0
            // 1     0        3.14
            // 2     1.23     0
            // 3     1.23     3.14
            offdiag = MAGMA_minproduct_C_MAKE( 1.2345, 6.7890 ) * (ival / 2);
            diag    = MAGMA_minproduct_C_MAKE( 3.1415, 2.7183 ) * (ival % 2);
            
            M = opts.msize[itest];
            N = opts.nsize[itest];
            //M += 2;  // space for insets
            //N += 2;
            lda    = M;
            ldda   = roundup( M, opts.roundup );
            size   = lda*N;
            if ( uplo[iuplo] == Magma_minproductLower || uplo[iuplo] == Magma_minproductUpper ) {
                // save triangle (with diagonal)
                // TODO wrong for trapezoid
                gbytes = sizeof(magma_minproductFloatComplex) * 0.5*N*(N+1) / 1e9;
            }
            else {
                // save entire matrix
                gbytes = sizeof(magma_minproductFloatComplex) * 1.*M*N / 1e9;
            }
    
            TESTING_MALLOC_CPU( h_A, magma_minproductFloatComplex, size   );
            TESTING_MALLOC_CPU( h_R, magma_minproductFloatComplex, size   );
            
            TESTING_MALLOC_DEV( d_A, magma_minproductFloatComplex, ldda*N );
            
            /* Initialize the matrix */
            for( int j = 0; j < N; ++j ) {
                for( int i = 0; i < M; ++i ) {
                    h_A[i + j*lda] = MAGMA_minproduct_C_MAKE( i + j/10000., j );
                }
            }
            
            /* ====================================================================
               Performs operation using MAGMA_minproduct
               =================================================================== */
            magma_minproduct_csetmatrix( M, N, h_A, lda, d_A, ldda );
            
            gpu_time = magma_minproduct_sync_wtime( 0 );
            //magma_minproductblas_claset( uplo[iuplo], M-2, N-2, offdiag, diag, d_A+1+ldda, ldda );  // inset by 1 row & col
            magma_minproductblas_claset( uplo[iuplo], M, N, offdiag, diag, d_A, ldda );
            gpu_time = magma_minproduct_sync_wtime( 0 ) - gpu_time;
            gpu_perf = gbytes / gpu_time;
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_minproduct_wtime();
            //magma_minproduct_int_t M2 = M-2;  // inset by 1 row & col
            //magma_minproduct_int_t N2 = N-2;
            //lapackf77_claset( lapack_uplo_const( uplo[iuplo] ), &M2, &N2, &offdiag, &diag, h_A+1+lda, &lda );
            lapackf77_claset( lapack_uplo_const( uplo[iuplo] ), &M, &N, &offdiag, &diag, h_A, &lda );
            cpu_time = magma_minproduct_wtime() - cpu_time;
            cpu_perf = gbytes / cpu_time;
            
            if ( opts.verbose ) {
                printf( "A= " );  magma_minproduct_cprint(     M, N, h_A, lda );
                printf( "dA=" );  magma_minproduct_cprint_gpu( M, N, d_A, ldda );
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            magma_minproduct_cgetmatrix( M, N, d_A, ldda, h_R, lda );
            
            blasf77_caxpy(&size, &c_neg_one, h_A, &ione, h_R, &ione);
            error = lapackf77_clange("f", &M, &N, h_R, &lda, work);

            bool okay = (error == 0);
            status += ! okay;
            printf("%5s %5d %5d  %7.2f  %4.2f   %7.2f (%7.2f)   %7.2f (%7.2f)   %s\n",
                   lapack_uplo_const( uplo[iuplo] ), (int) M, (int) N,
                   real(offdiag), real(diag),
                   cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000.,
                   (okay ? "ok" : "failed") );
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_R );
            
            TESTING_FREE_DEV( d_A );
            fflush( stdout );
          }
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
