/*
   -- MAGMA_tally2 (version 1.6.1) --
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
#include "magma_tally2.h"
#include "magma_tally2_lapack.h"
#include "testings.h"

double get_LU_error(magma_tally2_int_t M, magma_tally2_int_t N,
                    magma_tally2DoubleComplex *A,  magma_tally2_int_t lda,
                    magma_tally2DoubleComplex *LU, magma_tally2_int_t *IPIV)
{
    magma_tally2_int_t min_mn = min(M,N);
    magma_tally2_int_t ione   = 1;
    magma_tally2_int_t i, j;
    magma_tally2DoubleComplex alpha = MAGMA_tally2_Z_ONE;
    magma_tally2DoubleComplex beta  = MAGMA_tally2_Z_ZERO;
    magma_tally2DoubleComplex *L, *U;
    double work[1], matnorm, residual;
    
    TESTING_MALLOC_CPU( L, magma_tally2DoubleComplex, M*min_mn);
    TESTING_MALLOC_CPU( U, magma_tally2DoubleComplex, min_mn*N);
    memset( L, 0, M*min_mn*sizeof(magma_tally2DoubleComplex) );
    memset( U, 0, min_mn*N*sizeof(magma_tally2DoubleComplex) );

    lapackf77_zlaswp( &N, A, &lda, &ione, &min_mn, IPIV, &ione);
    lapackf77_zlacpy( Magma_tally2LowerStr, &M, &min_mn, LU, &lda, L, &M      );
    lapackf77_zlacpy( Magma_tally2UpperStr, &min_mn, &N, LU, &lda, U, &min_mn );

    for(j=0; j<min_mn; j++)
        L[j+j*M] = MAGMA_tally2_Z_MAKE( 1., 0. );
    
    matnorm = lapackf77_zlange("f", &M, &N, A, &lda, work);

    blasf77_zgemm("N", "N", &M, &N, &min_mn,
                  &alpha, L, &M, U, &min_mn, &beta, LU, &lda);

    for( j = 0; j < N; j++ ) {
        for( i = 0; i < M; i++ ) {
            LU[i+j*lda] = MAGMA_tally2_Z_SUB( LU[i+j*lda], A[i+j*lda] );
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

    real_Double_t   gflops, magma_tally2_perf, magma_tally2_time, cublas_perf, cublas_time, cpu_perf=0, cpu_time=0;
    double          error=0.0; 
    magma_tally2_int_t cublas_enable = 0;
    magma_tally2DoubleComplex *h_A, *h_R, *h_Amagma_tally2;
    magma_tally2DoubleComplex *dA;
    magma_tally2DoubleComplex **dA_array = NULL;

    magma_tally2_int_t     **dipiv_array = NULL;
    magma_tally2_int_t     *ipiv;
    magma_tally2_int_t     *dipiv_magma_tally2, *dinfo_magma_tally2;
    int             *dipiv_cublas, *dinfo_cublas;
    

    magma_tally2_int_t M, N, n2, lda, ldda, min_mn, info;
    magma_tally2_int_t ione     = 1; 
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    magma_tally2_int_t batchCount = 1;

    magma_tally2_queue_t queue = magma_tally2_stream;
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );
    //opts.lapack |= opts.check; 
    magma_tally2_int_t     status = 0;
    double tol = opts.tolerance * lapackf77_dlamch("E");

    batchCount = opts.batchcount ;
    magma_tally2_int_t columns;
    
    printf("BatchCount      M     N     CPU GFlop/s (ms)    MAGMA_tally2 GFlop/s (ms)  CUBLAS GFlop/s (ms)  ||PA-LU||/(||A||*N)\n");
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
            

            TESTING_MALLOC_CPU(    ipiv, magma_tally2_int_t,     min_mn * batchCount);
            TESTING_MALLOC_CPU(    h_A,  magma_tally2DoubleComplex, n2     );
            TESTING_MALLOC_CPU(    h_Amagma_tally2,  magma_tally2DoubleComplex, n2     );
            TESTING_MALLOC_PIN(    h_R,  magma_tally2DoubleComplex, n2     );
            
            TESTING_MALLOC_DEV(  dA,  magma_tally2DoubleComplex, ldda*N * batchCount);
            TESTING_MALLOC_DEV(  dipiv_magma_tally2,  magma_tally2_int_t, min_mn * batchCount);
            TESTING_MALLOC_DEV(  dinfo_magma_tally2,  magma_tally2_int_t, batchCount);
            TESTING_MALLOC_DEV(  dipiv_cublas,  magma_tally2_int_t, min_mn * batchCount);
            TESTING_MALLOC_DEV(  dinfo_cublas,  magma_tally2_int_t, batchCount);

            magma_tally2_malloc((void**)&dA_array, batchCount * sizeof(*dA_array));
            magma_tally2_malloc((void**)&dipiv_array, batchCount * sizeof(*dipiv_array));

            /* Initialize the matrix */
            lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
            columns = N * batchCount;
            lapackf77_zlacpy( Magma_tally2UpperLowerStr, &M, &columns, h_A, &lda, h_R, &lda );
            
            /* ====================================================================
               Performs operation using MAGMA_tally2
               =================================================================== */

            magma_tally2_zsetmatrix( M, columns, h_R, lda, dA, ldda );
            zset_pointer(dA_array, dA, ldda, 0, 0, ldda*N, batchCount, queue);
            set_ipointer(dipiv_array, dipiv_magma_tally2, 1, 0, 0, min_mn, batchCount, queue);
            
            magma_tally2_time = magma_tally2_sync_wtime(0);
            info = magma_tally2_zgetrf_batched( M, N, dA_array, ldda, dipiv_array,  dinfo_magma_tally2, batchCount, queue);
            magma_tally2_time = magma_tally2_sync_wtime(0) - magma_tally2_time;
            magma_tally2_perf = gflops / magma_tally2_time;
            
            magma_tally2_zgetmatrix( M, N*batchCount, dA, ldda, h_Amagma_tally2, lda );
            
            // check correctness of results throught "dinfo_magma_tally2" and correctness of argument throught "info"
            magma_tally2_int_t *cpu_info = (magma_tally2_int_t*) malloc(batchCount*sizeof(magma_tally2_int_t));
            magma_tally2_getvector( batchCount, sizeof(magma_tally2_int_t), dinfo_magma_tally2, 1, cpu_info, 1);
            
            for(int i=0; i<batchCount; i++)
            {
                if(cpu_info[i] != 0 ){
                    printf("magma_tally2_zgetrf_batched matrix %d returned internal error %d\n",i, (int)cpu_info[i] );
                }
            }
            
            if (info != 0)
                printf("magma_tally2_zgetrf_batched returned argument error %d: %s.\n", (int) info, magma_tally2_strerror( info ));

            /* ====================================================================
               Performs operation using CUBLAS
               =================================================================== */

            magma_tally2_zsetmatrix( M, columns, h_R, lda, dA,  ldda );
            zset_pointer(dA_array, dA, ldda, 0, 0, ldda * N, batchCount, queue);

            cublas_time = magma_tally2_sync_wtime(0);
            if(M == N )
            {
                cublasZgetrfBatched( opts.handle, N, dA_array, ldda, dipiv_cublas,  dinfo_cublas, batchCount);
                cublas_enable = 1;
            }
            else
            {
                printf("M != N, CUBLAS required M == N ; CUBLAS is disabled\n");
            }
            cublas_time = magma_tally2_sync_wtime(0) - cublas_time;
            cublas_perf = gflops / cublas_time;


            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally2_wtime();
                //#pragma unroll
                for(int i=0; i<batchCount; i++)
                {
                  lapackf77_zgetrf(&M, &N, h_A + i*lda*N, &lda, ipiv + i * min_mn, &info);
                }
                cpu_time = magma_tally2_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_zgetrf returned error %d: %s.\n",
                           (int) info, magma_tally2_strerror( info ));
            }
            
            /* =====================================================================
               Check the factorization
               =================================================================== */
            if ( opts.lapack ) {
                printf("%10d   %5d  %5d     %7.2f (%7.2f)   %7.2f (%7.2f)    %7.2f (%7.2f)",
                       (int) batchCount, (int) M, (int) N, cpu_perf, cpu_time*1000., magma_tally2_perf, magma_tally2_time*1000., cublas_perf*cublas_enable, cublas_time*1000.*cublas_enable  );
            }
            else {
                printf("%10d   %5d  %5d     ---   (  ---  )   %7.2f (%7.2f)    %7.2f (%7.2f)",
                       (int) batchCount, (int) M, (int) N, magma_tally2_perf, magma_tally2_time*1000., cublas_perf*cublas_enable, cublas_time*1000.*cublas_enable );
            }

            double err = 0.0;
            if ( opts.check ) {
                magma_tally2_getvector( min_mn * batchCount, sizeof(magma_tally2_int_t), dipiv_magma_tally2, 1, ipiv, 1 );
                int stop=0;
                for(int i=0; i<batchCount; i++)
                {
                    
                    for(int k=0;k<min_mn; k++){
                        if(ipiv[i*min_mn+k] < 1 || ipiv[i*min_mn+k] > M )
                        {
                                printf("error for matrix %d ipiv @ %d = %d\n",i,k,(int)ipiv[i*min_mn+k]);
                                stop=1;
                        }
                    }
                    if(stop==1){
                        err=-1.0;
                        break;
                    }
                    
                    error = get_LU_error( M, N, h_R + i * lda*N, lda, h_Amagma_tally2 + i * lda*N, ipiv + i * min_mn);                 
                    if ( isnan(error) || isinf(error) ) {
                        err = error;
                        break;
                    }
                    err = max(fabs(error),err);
                }
                printf("   %8.2e   %s\n", err, (err < tol ? "ok" : "failed") );
                status += ! (error < tol);
            }
            else {
                printf("     ---  \n");
            }
            
            TESTING_FREE_CPU( ipiv );
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_Amagma_tally2 );
            TESTING_FREE_PIN( h_R );

            TESTING_FREE_DEV( dA );
            TESTING_FREE_DEV( dinfo_magma_tally2 );
            TESTING_FREE_DEV( dipiv_magma_tally2 );
            TESTING_FREE_DEV( dipiv_cublas );
            TESTING_FREE_DEV( dinfo_cublas );
            TESTING_FREE_DEV( dipiv_array );
            TESTING_FREE_DEV( dA_array );
            free(cpu_info);
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    TESTING_FINALIZE();
    return status;
}
