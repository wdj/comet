/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zherk_batched.cpp normal z -> s, Fri Jan 30 19:00:26 2015
       @author Chongxiao Cao
       @author Tingxing Dong
*/


// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "testings.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing ssyrk_batched
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, magma_tally3_perf, magma_tally3_time, cpu_perf=0., cpu_time=0.;
    float          current_error, magma_tally3_error, Cnorm, work[1];
    magma_tally3_int_t N, K;
    magma_tally3_int_t Ak, An;
    magma_tally3_int_t sizeA, sizeC;
    magma_tally3_int_t lda, ldc, ldda, lddc;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t NN;
    magma_tally3_int_t batchCount;
 
    float *h_A, *h_C, *h_Cmagma_tally3;
    float *d_A, *d_C;
    float c_neg_one = MAGMA_tally3_S_NEG_ONE;
    float alpha = 0.29;
    float beta  = -0.48;
    float **A_array = NULL;
    float **C_array = NULL;
    magma_tally3_int_t status = 0;

    magma_tally3_queue_t queue = magma_tally3_stream;
    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    batchCount = opts.batchcount;

    float tol = opts.tolerance * lapackf77_slamch("E");
    
    printf("If running lapack (option --lapack), MAGMA_tally3 error is computed\n"
           "relative to CPU BLAS result.\n\n");
    printf("uplo = %s, transA = %s\n",
           lapack_uplo_const_tally3(opts.uplo), lapack_trans_const_tally3(opts.transA) );
    printf(" BatchCount    N     K   MAGMA_tally3 Gflop/s (ms)  CPU Gflop/s (ms)  MAGMA_tally3 error \n");
    printf("=========================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            K = opts.ksize[itest];
            gflops = FLOPS_SSYRK( K, N ) / 1e9 * batchCount;

            if ( opts.transA == Magma_tally3NoTrans ) {
                lda = An = N;
                Ak = K;
            } else {
                lda = An = K;
                Ak = N;
            }

            ldc = N;

            ldda = ((lda+31)/32)*32;
            lddc = ((ldc+31)/32)*32;
            
            NN = N * batchCount;

            
            sizeA = lda*Ak*batchCount;
            sizeC = ldc*N*batchCount;
            
            TESTING_MALLOC_CPU( h_A,  float, sizeA );
            TESTING_MALLOC_CPU( h_C,  float, sizeC );
            TESTING_MALLOC_CPU( h_Cmagma_tally3,  float, sizeC  );
            
            TESTING_MALLOC_DEV( d_A, float, ldda*An*batchCount );
            TESTING_MALLOC_DEV( d_C, float, lddc*N*batchCount );

            magma_tally3_malloc((void**)&A_array, batchCount * sizeof(*A_array));
            magma_tally3_malloc((void**)&C_array, batchCount * sizeof(*C_array));

            /* Initialize the matrices */
            lapackf77_slarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_slarnv( &ione, ISEED, &sizeC, h_C );
            
            /* =====================================================================
               Performs operation using MAGMA_tally3BLAS
               =================================================================== */
            magma_tally3_ssetmatrix( An, Ak*batchCount, h_A, lda, d_A, ldda );
            magma_tally3_ssetmatrix( N, N*batchCount, h_C, ldc, d_C, lddc );
            
            sset_pointer(A_array, d_A, lda, 0, 0, ldda*Ak, batchCount, queue);
            sset_pointer(C_array, d_C, ldc, 0, 0, lddc*N,  batchCount, queue);

            magma_tally3_time = magma_tally3_sync_wtime( NULL );
            magma_tally3blas_ssyrk_batched(opts.uplo, opts.transA, N, K,
                             alpha, A_array, ldda,
                             beta,  C_array, lddc, batchCount, queue);
                             
            magma_tally3_time = magma_tally3_sync_wtime( NULL ) - magma_tally3_time;
            magma_tally3_perf = gflops / magma_tally3_time;
            
            magma_tally3_sgetmatrix( N, NN, d_C, lddc, h_Cmagma_tally3, ldc );
            
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally3_wtime();
                for(int i=0; i<batchCount; i++)
                {
                   blasf77_ssyrk(
                               lapack_uplo_const_tally3(opts.uplo), lapack_trans_const_tally3(opts.transA),
                               &N, &K,
                               &alpha, h_A + i*lda*Ak, &lda,
                               &beta,  h_C + i*ldc*N, &ldc );
                }
                cpu_time = magma_tally3_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.check ) {
                // compute relative error for magma_tally3, relative to lapack,
                // |C_magma_tally3 - C_lapack| / |C_lapack|
                sizeC = ldc*N;
                magma_tally3_error = MAGMA_tally3_D_ZERO;
                for(int i=0; i<batchCount; i++)
                {
                    Cnorm = lapackf77_slansy("fro", lapack_uplo_const_tally3(opts.uplo), &N, h_C+i*ldc*N, &ldc, work);
                    blasf77_saxpy( &sizeC, &c_neg_one, h_C+i*ldc*N, &ione, h_Cmagma_tally3+i*ldc*N, &ione );
                    current_error = lapackf77_slansy( "fro", lapack_uplo_const_tally3(opts.uplo), &N, h_Cmagma_tally3+i*ldc*N, &ldc, work ) / Cnorm;
                    if ( isnan(current_error) || isinf(current_error) ) {
                        magma_tally3_error = current_error;
                        break;
                    }
                    magma_tally3_error = max(magma_tally3_error, current_error);
                }

                printf("%5d %5d %5d  %7.2f (%7.2f)   %7.2f (%7.2f)    %8.2e   %s\n",
                       (int) batchCount, (int) N, (int) K,
                       magma_tally3_perf, 1000.*magma_tally3_time,
                       cpu_perf,    1000.*cpu_time,
                       magma_tally3_error, (magma_tally3_error < tol ? "ok" : "failed"));

                status += ! (magma_tally3_error < tol);

            }
            else {
                printf("%5d %5d %5d   %7.2f (%7.2f)    ---   (  ---  )    ---     ---\n",
                       (int) batchCount, (int) N, (int) K,
                       magma_tally3_perf, 1000.*magma_tally3_time);

            }
            
            TESTING_FREE_CPU( h_A  );
            TESTING_FREE_CPU( h_C  );
            TESTING_FREE_CPU( h_Cmagma_tally3  );

            
            TESTING_FREE_DEV( d_A );
            TESTING_FREE_DEV( d_C );
            TESTING_FREE_DEV( A_array );
            TESTING_FREE_DEV( C_array );
            fflush( stdout);

        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
