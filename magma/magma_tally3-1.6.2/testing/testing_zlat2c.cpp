/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @precisions mixed zc -> ds
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// includes, project
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "magma_tally3_operators.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zlat2c and clat2z
*/
int main( int argc, char** argv )
{
    #define  A(i_,j_) ( A + (i_) + (j_)*lda)
    #define SA(i_,j_) (SA + (i_) + (j_)*lda)
    
    TESTING_INIT();
    
    real_Double_t   gbytes, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double error, work[1];
    float serror, swork[1];
    magma_tally3DoubleComplex c_neg_one = MAGMA_tally3_Z_NEG_ONE;
    magma_tally3FloatComplex  s_neg_one = MAGMA_tally3_C_NEG_ONE;
    magma_tally3_int_t ione = 1;
    magma_tally3_int_t n, lda, ldda, size, info;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t status = 0;
    magma_tally3FloatComplex   *SA, *SR;
    magma_tally3DoubleComplex   *A,  *R;
    magma_tally3FloatComplex  *dSA;
    magma_tally3DoubleComplex_ptr dA;
    
    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    
    magma_tally3_uplo_t uplo[] = { Magma_tally3Lower, Magma_tally3Upper };
    
    printf("func    uplo     N     CPU GB/s (ms)       GPU GB/s (ms)     ||R||_F\n");
    printf("=====================================================================\n");
    for( int iuplo = 0; iuplo < 2; ++iuplo ) {
      for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            n = opts.nsize[itest];
            lda  = n;
            ldda = ((n+31)/32)*32;
            // 0.5*(n+1)*n double-complex loads and 0.5*(n+1)*n single-complex stores (and vice-versa for clat2z)
            gbytes = (real_Double_t) 0.5*(n+1)*n * (sizeof(magma_tally3DoubleComplex) + sizeof(magma_tally3FloatComplex)) / 1e9;
            size = ldda*n;  // ldda >= lda
            
            TESTING_MALLOC_CPU(  SA, magma_tally3FloatComplex,  size );
            TESTING_MALLOC_CPU(   A, magma_tally3DoubleComplex, size );
            TESTING_MALLOC_CPU(  SR, magma_tally3FloatComplex,  size );
            TESTING_MALLOC_CPU(   R, magma_tally3DoubleComplex, size );
            
            TESTING_MALLOC_DEV( dSA, magma_tally3FloatComplex,  size );
            TESTING_MALLOC_DEV(  dA, magma_tally3DoubleComplex, size );
            
            lapackf77_zlarnv( &ione, ISEED, &size,  A );
            lapackf77_clarnv( &ione, ISEED, &size, SA );
            
            magma_tally3_zsetmatrix( n, n, A,  lda, dA,  ldda );
            magma_tally3_csetmatrix( n, n, SA, lda, dSA, ldda );
            
            /* =====================================================================
               Performs operation using LAPACK zlat2c
               =================================================================== */
            info = 0;
            cpu_time = magma_tally3_wtime();
            lapackf77_zlat2c( lapack_uplo_const_tally3(uplo[iuplo]), &n, A, &lda, SA, &lda, &info );
            cpu_time = magma_tally3_wtime() - cpu_time;
            cpu_perf = gbytes / cpu_time;
            if (info != 0)
                printf("lapackf77_zlat2c returned error %d: %s.\n",
                       (int) info, magma_tally3_strerror( info ));
            
            /* ====================================================================
               Performs operation using MAGMA_tally3 zlat2c
               =================================================================== */
            gpu_time = magma_tally3_sync_wtime(0);
            magma_tally3blas_zlat2c( uplo[iuplo], n, dA, ldda, dSA, ldda, &info );
            gpu_time = magma_tally3_sync_wtime(0) - gpu_time;
            gpu_perf = gbytes / gpu_time;
            if (info != 0)
                printf("magma_tally3blas_zlat2c returned error %d: %s.\n",
                       (int) info, magma_tally3_strerror( info ));
            
            magma_tally3_cgetmatrix( n, n, dSA, ldda, SR, lda );
            
            if ( opts.verbose ) {
                printf( "A=  " );  magma_tally3_zprint( n, n, A,  lda );
                printf( "SA= " );  magma_tally3_cprint( n, n, SA, lda );
                printf( "dA= " );  magma_tally3_zprint_gpu( n, n, dA,  ldda );
                printf( "dSA=" );  magma_tally3_cprint_gpu( n, n, dSA, ldda );
            }
            
            /* =====================================================================
               compute error |SA_magma_tally3 - SA_lapack|
               should be zero if both are IEEE compliant
               =================================================================== */
            blasf77_caxpy( &size, &s_neg_one, SA, &ione, SR, &ione );
            serror = lapackf77_clange( "Fro", &n, &n, SR, &lda, swork );
            
            printf( "zlat2c %5s %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                    lapack_uplo_const_tally3(uplo[iuplo]), (int) n,
                    cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000.,
                    serror, (serror == 0 ? "ok" : "failed") );
            status += ! (serror == 0);
            
            
            /* =====================================================================
               Reset matrices
               =================================================================== */
            lapackf77_zlarnv( &ione, ISEED, &size,  A );
            lapackf77_clarnv( &ione, ISEED, &size, SA );
            
            magma_tally3_zsetmatrix( n, n, A,  lda, dA,  ldda );
            magma_tally3_csetmatrix( n, n, SA, lda, dSA, ldda );
            
            /* =====================================================================
               Performs operation using LAPACK clat2z
               LAPACK doesn't implement clat2z; use our own simple implementation.
               =================================================================== */
            cpu_time = magma_tally3_wtime();
            if ( uplo[iuplo] == Magma_tally3Lower ) {
                for( int j=0; j < n; ++j ) {
                    for( int i=j; i < n; ++i ) {
                        *A(i,j) = MAGMA_tally3_Z_MAKE( real(*SA(i,j)), imag(*SA(i,j)) );
                    }
                }
            }
            else { // upper
                for( int j=0; j < n; ++j ) {
                    for( int i=0; i <= j; ++i ) {
                        *A(i,j) = MAGMA_tally3_Z_MAKE( real(*SA(i,j)), imag(*SA(i,j)) );
                    }
                }
            }
            cpu_time = magma_tally3_wtime() - cpu_time;
            cpu_perf = gbytes / cpu_time;
            if (info != 0)
                printf("lapackf77_clat2z returned error %d: %s.\n",
                       (int) info, magma_tally3_strerror( info ));
            
            /* ====================================================================
               Performs operation using MAGMA_tally3 clat2z
               =================================================================== */
            magma_tally3_csetmatrix( n, n, SA, lda, dSA, ldda );
            
            gpu_time = magma_tally3_sync_wtime(0);
            magma_tally3blas_clat2z( uplo[iuplo], n, dSA, ldda, dA, ldda, &info );
            gpu_time = magma_tally3_sync_wtime(0) - gpu_time;
            gpu_perf = gbytes / gpu_time;
            if (info != 0)
                printf("magma_tally3blas_clat2z returned error %d: %s.\n",
                       (int) info, magma_tally3_strerror( info ));
            
            magma_tally3_zgetmatrix( n, n, dA, ldda, R, lda );
            
            if ( opts.verbose ) {
                printf( "A=  " );  magma_tally3_zprint( n, n, A,  lda );
                printf( "SA= " );  magma_tally3_cprint( n, n, SA, lda );
                printf( "dA= " );  magma_tally3_zprint_gpu( n, n, dA,  ldda );
                printf( "dSA=" );  magma_tally3_cprint_gpu( n, n, dSA, ldda );
            }
            
            /* =====================================================================
               compute error |A_magma_tally3 - A_lapack|
               should be zero if both are IEEE compliant
               =================================================================== */
            blasf77_zaxpy( &size, &c_neg_one, A, &ione, R, &ione );
            error = lapackf77_zlange( "Fro", &n, &n, R, &lda, work );
            
            printf( "clat2z %5s %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                    lapack_uplo_const_tally3(uplo[iuplo]), (int) n,
                    cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000.,
                    error, (error == 0 ? "ok" : "failed") );
            status += ! (error == 0);
            
            TESTING_FREE_CPU(  SA );
            TESTING_FREE_CPU(   A );
            TESTING_FREE_CPU(  SR );
            TESTING_FREE_CPU(   R );
                                 
            TESTING_FREE_DEV( dSA );
            TESTING_FREE_DEV(  dA );
            printf( "\n" );
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
