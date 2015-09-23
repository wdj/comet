/*
   -- MAGMA_minproduct (version 1.6.1) --
   Univ. of Tennessee, Knoxville
   Univ. of California, Berkeley
   Univ. of Colorado, Denver
   @date January 2015

   @author Azzam Haidar
   @author Tingxing Dong

   @generated from testing_zgetrf_nopiv_batched.cpp normal z -> s, Fri Jan 30 19:00:26 2015
 */
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "testings.h"

float get_LU_error(magma_minproduct_int_t M, magma_minproduct_int_t N,
                    float *A,  magma_minproduct_int_t lda,
                    float *LU, magma_minproduct_int_t *IPIV)
{
    magma_minproduct_int_t min_mn = min(M,N);
    magma_minproduct_int_t ione   = 1;
    magma_minproduct_int_t i, j;
    float alpha = MAGMA_minproduct_S_ONE;
    float beta  = MAGMA_minproduct_S_ZERO;
    float *L, *U;
    float work[1], matnorm, residual;
    
    TESTING_MALLOC_CPU( L, float, M*min_mn);
    TESTING_MALLOC_CPU( U, float, min_mn*N);
    memset( L, 0, M*min_mn*sizeof(float) );
    memset( U, 0, min_mn*N*sizeof(float) );

    lapackf77_slaswp( &N, A, &lda, &ione, &min_mn, IPIV, &ione);
    lapackf77_slacpy( Magma_minproductLowerStr, &M, &min_mn, LU, &lda, L, &M      );
    lapackf77_slacpy( Magma_minproductUpperStr, &min_mn, &N, LU, &lda, U, &min_mn );

    for(j=0; j<min_mn; j++)
        L[j+j*M] = MAGMA_minproduct_S_MAKE( 1., 0. );
    
    matnorm = lapackf77_slange("f", &M, &N, A, &lda, work);

    blasf77_sgemm("N", "N", &M, &N, &min_mn,
                  &alpha, L, &M, U, &min_mn, &beta, LU, &lda);

    for( j = 0; j < N; j++ ) {
        for( i = 0; i < M; i++ ) {
            LU[i+j*lda] = MAGMA_minproduct_S_SUB( LU[i+j*lda], A[i+j*lda] );
        }
    }
    residual = lapackf77_slange("f", &M, &N, LU, &lda, work);

    TESTING_FREE_CPU(L);
    TESTING_FREE_CPU(U);

    return residual / (matnorm * N);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sgetrf_batched
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, magma_minproduct_perf, magma_minproduct_time, cublas_perf=0., cublas_time=0., cpu_perf=0, cpu_time=0;
    float          error; 
    magma_minproduct_int_t cublas_enable = 0;
    float *h_A, *h_R;
    float *dA_magma_minproduct;
    float **dA_array = NULL;

    magma_minproduct_int_t     **dipiv_array = NULL;
    magma_minproduct_int_t     *ipiv;
    magma_minproduct_int_t     *dipiv_magma_minproduct, *dinfo_magma_minproduct;
    

    magma_minproduct_int_t M, N, n2, lda, ldda, min_mn, info;
    magma_minproduct_int_t ione     = 1; 
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t batchCount = 1;

    magma_minproduct_queue_t queue = magma_minproduct_stream;
    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    //opts.lapack |= opts.check; 

    batchCount = opts.batchcount ;
    magma_minproduct_int_t columns;
    
    printf("BatchCount      M     N     CPU GFlop/s (ms)    MAGMA_minproduct GFlop/s (ms)  CUBLAS GFlop/s (ms)  ||PA-LU||/(||A||*N)\n");
    printf("=========================================================================\n");
    for( int i = 0; i < opts.ntest; ++i ) {
    
      for( int iter = 0; iter < opts.niter; ++iter ) {
            
            M = opts.msize[i];
            N = opts.nsize[i];
            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N * batchCount;
            ldda   = ((M+31)/32)*32;
            gflops = FLOPS_SGETRF( M, N ) / 1e9 * batchCount;
            

            TESTING_MALLOC_CPU(    ipiv, magma_minproduct_int_t,     min_mn * batchCount);
            TESTING_MALLOC_CPU(    h_A,  float, n2     );
            TESTING_MALLOC_PIN( h_R,  float, n2     );
            TESTING_MALLOC_DEV(  dA_magma_minproduct,  float, ldda*N * batchCount);
            TESTING_MALLOC_DEV(  dipiv_magma_minproduct,  magma_minproduct_int_t, min_mn * batchCount);
            TESTING_MALLOC_DEV(  dinfo_magma_minproduct,  magma_minproduct_int_t, batchCount);

            magma_minproduct_malloc((void**)&dA_array, batchCount * sizeof(*dA_array));
            magma_minproduct_malloc((void**)&dipiv_array, batchCount * sizeof(*dipiv_array));

            /* Initialize the matrix */
            lapackf77_slarnv( &ione, ISEED, &n2, h_A );
            columns = N * batchCount;
            lapackf77_slacpy( Magma_minproductUpperLowerStr, &M, &columns, h_A, &lda, h_R, &lda );
            magma_minproduct_ssetmatrix( M, columns, h_R, lda, dA_magma_minproduct, ldda );
            


            /* ====================================================================
               Performs operation using MAGMA_minproduct
               =================================================================== */
            sset_pointer(dA_array, dA_magma_minproduct, ldda, 0, 0, ldda*N, batchCount, queue);
            magma_minproduct_time = magma_minproduct_sync_wtime(0);
            info = magma_minproduct_sgetrf_nopiv_batched( M, N, dA_array, ldda, dinfo_magma_minproduct, batchCount, queue);
            magma_minproduct_time = magma_minproduct_sync_wtime(0) - magma_minproduct_time;
            magma_minproduct_perf = gflops / magma_minproduct_time;
            // check correctness of results throught "dinfo_magma_minproduct" and correctness of argument throught "info"
            magma_minproduct_int_t *cpu_info = (magma_minproduct_int_t*) malloc(batchCount*sizeof(magma_minproduct_int_t));
            magma_minproduct_getvector( batchCount, sizeof(magma_minproduct_int_t), dinfo_magma_minproduct, 1, cpu_info, 1);
            for(int i=0; i<batchCount; i++)
            {
                if(cpu_info[i] != 0 ){
                    printf("magma_minproduct_sgetrf_batched matrix %d returned internal error %d\n",i, (int)cpu_info[i] );
                }
            }
            if (info != 0)
                printf("magma_minproduct_sgetrf_batched returned argument error %d: %s.\n", (int) info, magma_minproduct_strerror( info ));

            /* ====================================================================
               Performs operation using CUBLAS
               =================================================================== */

            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_minproduct_wtime();
                //#pragma unroll
                for(int i=0; i<batchCount; i++)
                {
                  lapackf77_sgetrf(&M, &N, h_A + i*lda*N, &lda, ipiv + i * min_mn, &info);
                }
                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_sgetrf returned error %d: %s.\n",
                           (int) info, magma_minproduct_strerror( info ));
            }
            




            /* =====================================================================
               Check the factorization
               =================================================================== */
            if ( opts.lapack ) {
                printf("%10d   %5d  %5d     %7.2f (%7.2f)   %7.2f (%7.2f)    %7.2f (%7.2f)",
                       (int) batchCount, (int) M, (int) N, cpu_perf, cpu_time*1000., magma_minproduct_perf, magma_minproduct_time*1000., cublas_perf*cublas_enable, cublas_time*1000.*cublas_enable  );
            }
            else {
                printf("%10d   %5d  %5d     ---   (  ---  )   %7.2f (%7.2f)    %7.2f (%7.2f)",
                       (int) batchCount, (int) M, (int) N, magma_minproduct_perf, magma_minproduct_time*1000., cublas_perf*cublas_enable, cublas_time*1000.*cublas_enable );
            }

            float err = 0.0;
            if ( opts.check ) {
               
                // initialize ipiv to 1,2,3,4,5,6
                for(int i=0; i<batchCount; i++)
                {
                    for(int k=0;k<min_mn; k++){
                        ipiv[i*min_mn+k] = k+1;
                    }
                }

                magma_minproduct_sgetmatrix( M, N*batchCount, dA_magma_minproduct, ldda, h_A, lda );
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
            TESTING_FREE_DEV( dA_magma_minproduct );
            TESTING_FREE_DEV( dinfo_magma_minproduct );
            TESTING_FREE_DEV( dipiv_magma_minproduct );
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
