/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zlacpy.cpp normal z -> d, Fri Jan 30 19:00:23 2015
       @author Mark Gates
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dlacpy
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gbytes, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double           error, work[1];
    double  c_neg_one = MAGMA_tally3_D_NEG_ONE;
    double *h_A, *h_B, *h_R;
    magma_tally3Double_ptr d_A, d_B;
    magma_tally3_int_t M, N, size, lda, ldb, ldda, lddb;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t status = 0;
    
    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );

    magma_tally3_uplo_t uplo[] = { Magma_tally3Lower, Magma_tally3Upper, Magma_tally3Full };
    
    printf("uplo      M     N   CPU GByte/s (ms)    GPU GByte/s (ms)    check\n");
    printf("=================================================================\n");
    for( int iuplo = 0; iuplo < 3; ++iuplo ) {
      for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            lda    = M;
            ldb    = lda;
            ldda   = ((M+31)/32)*32;
            lddb   = ldda;
            size   = lda*N;
            if ( uplo[iuplo] == Magma_tally3Lower || uplo[iuplo] == Magma_tally3Upper ) {
                // load and save triangle (with diagonal)
                gbytes = sizeof(double) * 1.*N*(N+1) / 1e9;
            }
            else {
              // load entire matrix, save entire matrix
              gbytes = sizeof(double) * 2.*M*N / 1e9;
            }
    
            TESTING_MALLOC_CPU( h_A, double, size   );
            TESTING_MALLOC_CPU( h_B, double, size   );
            TESTING_MALLOC_CPU( h_R, double, size   );
            
            TESTING_MALLOC_DEV( d_A, double, ldda*N );
            TESTING_MALLOC_DEV( d_B, double, ldda*N );
            
            /* Initialize the matrix */
            for( int j = 0; j < N; ++j ) {
                for( int i = 0; i < M; ++i ) {
                    h_A[i + j*lda] = MAGMA_tally3_D_MAKE( i + j/10000., j );
                    h_B[i + j*ldb] = MAGMA_tally3_D_MAKE( i - j/10000. + 10000., j );
                }
            }
            
            /* ====================================================================
               Performs operation using MAGMA_tally3
               =================================================================== */
            magma_tally3_dsetmatrix( M, N, h_A, lda, d_A, ldda );
            magma_tally3_dsetmatrix( M, N, h_B, ldb, d_B, lddb );
            
            gpu_time = magma_tally3_sync_wtime( 0 );
            //magma_tally3blas_dlacpy( uplo[iuplo], M-2, N-2, d_A+1+ldda, ldda, d_B+1+lddb, lddb );  // inset by 1 row & col
            magma_tally3blas_dlacpy( uplo[iuplo], M, N, d_A, ldda, d_B, lddb );
            gpu_time = magma_tally3_sync_wtime( 0 ) - gpu_time;
            gpu_perf = gbytes / gpu_time;
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_tally3_wtime();
            //magma_tally3_int_t M2 = M-2;  // inset by 1 row & col
            //magma_tally3_int_t N2 = N-2;
            //lapackf77_dlacpy( lapack_uplo_const_tally3(uplo[iuplo]), &M2, &N2, h_A+1+lda, &lda, h_B+1+ldb, &ldb );
            lapackf77_dlacpy( lapack_uplo_const_tally3(uplo[iuplo]), &M, &N, h_A, &lda, h_B, &ldb );
            cpu_time = magma_tally3_wtime() - cpu_time;
            cpu_perf = gbytes / cpu_time;
            
            if ( opts.verbose ) {
                printf( "A= " );  magma_tally3_dprint(     M, N, h_A, lda );
                printf( "B= " );  magma_tally3_dprint(     M, N, h_B, lda );
                printf( "dA=" );  magma_tally3_dprint_gpu( M, N, d_A, ldda );
                printf( "dB=" );  magma_tally3_dprint_gpu( M, N, d_B, ldda );
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            magma_tally3_dgetmatrix( M, N, d_B, ldda, h_R, lda );
            
            blasf77_daxpy(&size, &c_neg_one, h_B, &ione, h_R, &ione);
            error = lapackf77_dlange("f", &M, &N, h_R, &lda, work);

            printf("%5s %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %s\n",
                   lapack_uplo_const_tally3(uplo[iuplo]), (int) M, (int) N,
                   cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000.,
                   (error == 0. ? "ok" : "failed") );
            status += ! (error == 0.);
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_B );
            TESTING_FREE_CPU( h_R );
            
            TESTING_FREE_DEV( d_A );
            TESTING_FREE_DEV( d_B );
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
