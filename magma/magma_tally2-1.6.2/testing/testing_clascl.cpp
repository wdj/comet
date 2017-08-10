/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zlascl.cpp normal z -> c, Fri Jan 30 19:00:24 2015
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
   -- Testing clascl
   Code is very similar to testing_clacpy.cpp
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gbytes, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float           error, work[1];
    magma_tally2FloatComplex  c_neg_one = MAGMA_tally2_C_NEG_ONE;
    magma_tally2FloatComplex *h_A, *h_R;
    magma_tally2FloatComplex_ptr d_A;
    float cto, cfrom;
    magma_tally2_int_t M, N, size, lda, ldda, info;
    magma_tally2_int_t ione     = 1;
    magma_tally2_int_t status = 0;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );

    magma_tally2_uplo_t uplo[] = { Magma_tally2Lower, Magma_tally2Upper, Magma_tally2Full };
    
    float sfmin = lapackf77_slamch("sfmin");
    float bignum = 1 / sfmin;
    
    printf("uplo      M     N    CPU GByte/s (ms)    GPU GByte/s (ms)   check\n");
    printf("====================================================================\n");
    for( int iuplo = 0; iuplo < 3; ++iuplo ) {
      for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            //M += 2;  // space for insets
            //N += 2;
            lda    = M;
            ldda   = ((M+31)/32)*32;
            size   = lda*N;
            if ( uplo[iuplo] == Magma_tally2Lower || uplo[iuplo] == Magma_tally2Upper ) {
                // read & write triangle (with diagonal)
                // TODO wrong for trapezoid
                gbytes = sizeof(magma_tally2FloatComplex) * 0.5*2.*N*(N+1) / 1e9;
            }
            else {
                // read & write entire matrix
                gbytes = sizeof(magma_tally2FloatComplex) * 2.*M*N / 1e9;
            }
    
            TESTING_MALLOC_CPU( h_A, magma_tally2FloatComplex, size   );
            TESTING_MALLOC_CPU( h_R, magma_tally2FloatComplex, size   );
            
            TESTING_MALLOC_DEV( d_A, magma_tally2FloatComplex, ldda*N );
            
            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &size, h_A );
            
            if ( 0 ) {
                // test that the over/underflow algorithm is working (but slower)
                // scale by 1e-4
                cto   = 1e-4;
                cfrom = 1;
                lapackf77_clascl( "G", &ione, &ione, &cfrom, &cto, &M, &N, h_A, &lda, &info );
                assert( info == 0 );
                
                // this (cto/cfrom) is inf,
                // but (1e-4 * cto/cfrom) is ~ 1e308 in float and ~ 1e37 in float
                cto   = 100.*sqrt( bignum );
                cfrom = 1/cto;
            }
            else {
                cto   = 1.2345;
                cfrom = 3.1415;
            }
                        
            /* ====================================================================
               Performs operation using MAGMA_tally2
               =================================================================== */
            magma_tally2_csetmatrix( M, N, h_A, lda, d_A, ldda );
            
            gpu_time = magma_tally2_sync_wtime( 0 );
            //magma_tally2blas_clascl( uplo[iuplo], 1, 1, cfrom, cto, M-2, N-2, d_A+1+ldda, ldda, &info );  // inset by 1 row & col
            magma_tally2blas_clascl( uplo[iuplo], 1, 1, cfrom, cto, M, N, d_A, ldda, &info );
            gpu_time = magma_tally2_sync_wtime( 0 ) - gpu_time;
            gpu_perf = gbytes / gpu_time;
            if (info != 0)
                printf("magma_tally2blas_clascl returned error %d: %s.\n",
                       (int) info, magma_tally2_strerror( info ));
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_tally2_wtime();
            //magma_tally2_int_t M2 = M-2;  // inset by 1 row & col
            //magma_tally2_int_t N2 = N-2;
            //lapackf77_clascl( lapack_uplo_const_tally2( uplo[iuplo] ), &ione, &ione, &cfrom, &cto, &M2, &N2, h_A+1+lda, &lda, &info );
            lapackf77_clascl( lapack_uplo_const_tally2( uplo[iuplo] ), &ione, &ione, &cfrom, &cto, &M, &N, h_A, &lda, &info );
            cpu_time = magma_tally2_wtime() - cpu_time;
            cpu_perf = gbytes / cpu_time;
            if (info != 0)
                printf("lapackf77_clascl returned error %d: %s.\n",
                       (int) info, magma_tally2_strerror( info ));
            
            /* =====================================================================
               Check the result
               =================================================================== */
            magma_tally2_cgetmatrix( M, N, d_A, ldda, h_R, lda );
            //magma_tally2_cprint( M, N, h_R, lda );
            
            blasf77_caxpy(&size, &c_neg_one, h_A, &ione, h_R, &ione);
            error = lapackf77_clange("f", &M, &N, h_R, &lda, work);

            printf("%5s %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                   lapack_uplo_const_tally2( uplo[iuplo] ), (int) M, (int) N,
                   cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000.,
                   error, (error == 0. ? "ok" : "failed") );
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