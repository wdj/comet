/*
   -- MAGMA_tally4 (version 1.6.1) --
   Univ. of Tennessee, Knoxville
   Univ. of California, Berkeley
   Univ. of Colorado, Denver
   @date January 2015

   @author Azzam Haidar
   @author Tingxing Dong

   @precisions normal z -> s d c
 */
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"
#include "testings.h"

double get_LU_error(magma_tally4_int_t M, magma_tally4_int_t N,
                    magma_tally4DoubleComplex *A,  magma_tally4_int_t lda,
                    magma_tally4DoubleComplex *LU, magma_tally4_int_t *IPIV)
{
    magma_tally4_int_t min_mn = min(M,N);
    magma_tally4_int_t ione   = 1;
    magma_tally4_int_t i, j;
    magma_tally4DoubleComplex alpha = MAGMA_tally4_Z_ONE;
    magma_tally4DoubleComplex beta  = MAGMA_tally4_Z_ZERO;
    magma_tally4DoubleComplex *L, *U;
    double work[1], matnorm, residual;
    
    TESTING_MALLOC_CPU( L, magma_tally4DoubleComplex, M*min_mn);
    TESTING_MALLOC_CPU( U, magma_tally4DoubleComplex, min_mn*N);
    memset( L, 0, M*min_mn*sizeof(magma_tally4DoubleComplex) );
    memset( U, 0, min_mn*N*sizeof(magma_tally4DoubleComplex) );

    lapackf77_zlaswp( &N, A, &lda, &ione, &min_mn, IPIV, &ione);
    lapackf77_zlacpy( Magma_tally4LowerStr, &M, &min_mn, LU, &lda, L, &M      );
    lapackf77_zlacpy( Magma_tally4UpperStr, &min_mn, &N, LU, &lda, U, &min_mn );

    for(j=0; j<min_mn; j++)
        L[j+j*M] = MAGMA_tally4_Z_MAKE( 1., 0. );
    
    matnorm = lapackf77_zlange("f", &M, &N, A, &lda, work);

    blasf77_zgemm("N", "N", &M, &N, &min_mn,
                  &alpha, L, &M, U, &min_mn, &beta, LU, &lda);

    for( j = 0; j < N; j++ ) {
        for( i = 0; i < M; i++ ) {
            LU[i+j*lda] = MAGMA_tally4_Z_SUB( LU[i+j*lda], A[i+j*lda] );
        }
    }
    residual = lapackf77_zlange("f", &M, &N, LU, &lda, work);

    TESTING_FREE_CPU(L);
    TESTING_FREE_CPU(U);

    return residual / (matnorm * N);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zgetrf_batched
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, magma_tally4_perf, magma_tally4_time, cublas_perf=0., cublas_time=0., cpu_perf=0, cpu_time=0;
    double          error; 
    magma_tally4_int_t cublas_enable = 0;
    magma_tally4DoubleComplex *h_A, *h_R;
    magma_tally4DoubleComplex *dA_magma_tally4;
    magma_tally4DoubleComplex **dA_array = NULL;

    magma_tally4_int_t     **dipiv_array = NULL;
    magma_tally4_int_t     *ipiv;
    magma_tally4_int_t     *dipiv_magma_tally4, *dinfo_magma_tally4;
    

    magma_tally4_int_t M, N, n2, lda, ldda, min_mn, info;
    magma_tally4_int_t ione     = 1; 
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    magma_tally4_int_t batchCount = 1;

    magma_tally4_queue_t queue = magma_tally4_stream;
    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );
    //opts.lapack |= opts.check; 

    batchCount = opts.batchcount ;
    magma_tally4_int_t columns;
    
    printf("BatchCount      M     N     CPU GFlop/s (ms)    MAGMA_tally4 GFlop/s (ms)  CUBLAS GFlop/s (ms)  ||PA-LU||/(||A||*N)\n");
    printf("=========================================================================\n");
    for( int i = 0; i < opts.ntest; ++i ) {
    
      for( int iter = 0; iter < opts.niter; ++iter ) {
            
            M = opts.msize[i];
            N = opts.nsize[i];
            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N * batchCount;
            ldda   = ((M+31)/32)*32;
            gflops = FLOPS_ZGETRF( M, N ) / 1e9 * batchCount;
            

            TESTING_MALLOC_CPU(    ipiv, magma_tally4_int_t,     min_mn * batchCount);
            TESTING_MALLOC_CPU(    h_A,  magma_tally4DoubleComplex, n2     );
            TESTING_MALLOC_PIN( h_R,  magma_tally4DoubleComplex, n2     );
            TESTING_MALLOC_DEV(  dA_magma_tally4,  magma_tally4DoubleComplex, ldda*N * batchCount);
            TESTING_MALLOC_DEV(  dipiv_magma_tally4,  magma_tally4_int_t, min_mn * batchCount);
            TESTING_MALLOC_DEV(  dinfo_magma_tally4,  magma_tally4_int_t, batchCount);

            magma_tally4_malloc((void**)&dA_array, batchCount * sizeof(*dA_array));
            magma_tally4_malloc((void**)&dipiv_array, batchCount * sizeof(*dipiv_array));

            /* Initialize the matrix */
            lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
            columns = N * batchCount;
            lapackf77_zlacpy( Magma_tally4UpperLowerStr, &M, &columns, h_A, &lda, h_R, &lda );
            magma_tally4_zsetmatrix( M, columns, h_R, lda, dA_magma_tally4, ldda );
            


            /* ====================================================================
               Performs operation using MAGMA_tally4
               =================================================================== */
            zset_pointer(dA_array, dA_magma_tally4, ldda, 0, 0, ldda*N, batchCount, queue);
            magma_tally4_time = magma_tally4_sync_wtime(0);
            info = magma_tally4_zgetrf_nopiv_batched( M, N, dA_array, ldda, dinfo_magma_tally4, batchCount, queue);
            magma_tally4_time = magma_tally4_sync_wtime(0) - magma_tally4_time;
            magma_tally4_perf = gflops / magma_tally4_time;
            // check correctness of results throught "dinfo_magma_tally4" and correctness of argument throught "info"
            magma_tally4_int_t *cpu_info = (magma_tally4_int_t*) malloc(batchCount*sizeof(magma_tally4_int_t));
            magma_tally4_getvector( batchCount, sizeof(magma_tally4_int_t), dinfo_magma_tally4, 1, cpu_info, 1);
            for(int i=0; i<batchCount; i++)
            {
                if(cpu_info[i] != 0 ){
                    printf("magma_tally4_zgetrf_batched matrix %d returned internal error %d\n",i, (int)cpu_info[i] );
                }
            }
            if (info != 0)
                printf("magma_tally4_zgetrf_batched returned argument error %d: %s.\n", (int) info, magma_tally4_strerror( info ));

            /* ====================================================================
               Performs operation using CUBLAS
               =================================================================== */

            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally4_wtime();
                //#pragma unroll
                for(int i=0; i<batchCount; i++)
                {
                  lapackf77_zgetrf(&M, &N, h_A + i*lda*N, &lda, ipiv + i * min_mn, &info);
                }
                cpu_time = magma_tally4_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_zgetrf returned error %d: %s.\n",
                           (int) info, magma_tally4_strerror( info ));
            }
            




            /* =====================================================================
               Check the factorization
               =================================================================== */
            if ( opts.lapack ) {
                printf("%10d   %5d  %5d     %7.2f (%7.2f)   %7.2f (%7.2f)    %7.2f (%7.2f)",
                       (int) batchCount, (int) M, (int) N, cpu_perf, cpu_time*1000., magma_tally4_perf, magma_tally4_time*1000., cublas_perf*cublas_enable, cublas_time*1000.*cublas_enable  );
            }
            else {
                printf("%10d   %5d  %5d     ---   (  ---  )   %7.2f (%7.2f)    %7.2f (%7.2f)",
                       (int) batchCount, (int) M, (int) N, magma_tally4_perf, magma_tally4_time*1000., cublas_perf*cublas_enable, cublas_time*1000.*cublas_enable );
            }

            double err = 0.0;
            if ( opts.check ) {
               
                // initialize ipiv to 1,2,3,4,5,6
                for(int i=0; i<batchCount; i++)
                {
                    for(int k=0;k<min_mn; k++){
                        ipiv[i*min_mn+k] = k+1;
                    }
                }

                magma_tally4_zgetmatrix( M, N*batchCount, dA_magma_tally4, ldda, h_A, lda );
                for(int i=0; i<batchCount; i++)
                {
                  error = get_LU_error( M, N, h_R + i * lda*N, lda, h_A + i * lda*N, ipiv + i * min_mn);                 
                  if ( isnan(error) || isinf(error) ) {
                      err = error;
                      break;
                  }
                  err = max(fabs(error),err);
                }
                 printf("   %8.2e\n", err );
            }
            else {
                printf("     ---  \n");
            }
            
            TESTING_FREE_CPU( ipiv );
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_PIN( h_R );
            TESTING_FREE_DEV( dA_magma_tally4 );
            TESTING_FREE_DEV( dinfo_magma_tally4 );
            TESTING_FREE_DEV( dipiv_magma_tally4 );
            TESTING_FREE_DEV( dipiv_array );
            TESTING_FREE_DEV( dA_array );
            free(cpu_info);
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    TESTING_FINALIZE();
    return 0;
}
