/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> c d s
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
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"
#include "testings.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zherk_batched
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, magma_tally4_perf, magma_tally4_time, cpu_perf=0., cpu_time=0.;
    double          current_error, magma_tally4_error, Cnorm, work[1];
    magma_tally4_int_t N, K;
    magma_tally4_int_t Ak, An;
    magma_tally4_int_t sizeA, sizeC;
    magma_tally4_int_t lda, ldc, ldda, lddc;
    magma_tally4_int_t ione     = 1;
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    magma_tally4_int_t NN;
    magma_tally4_int_t batchCount;
 
    magma_tally4DoubleComplex *h_A, *h_C, *h_Cmagma_tally4;
    magma_tally4DoubleComplex *d_A, *d_C;
    magma_tally4DoubleComplex c_neg_one = MAGMA_tally4_Z_NEG_ONE;
    double alpha = 0.29;
    double beta  = -0.48;
    magma_tally4DoubleComplex **A_array = NULL;
    magma_tally4DoubleComplex **C_array = NULL;
    magma_tally4_int_t status = 0;

    magma_tally4_queue_t queue = magma_tally4_stream;
    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    batchCount = opts.batchcount;

    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("If running lapack (option --lapack), MAGMA_tally4 error is computed\n"
           "relative to CPU BLAS result.\n\n");
    printf("uplo = %s, transA = %s\n",
           lapack_uplo_const_tally4(opts.uplo), lapack_trans_const_tally4(opts.transA) );
    printf(" BatchCount    N     K   MAGMA_tally4 Gflop/s (ms)  CPU Gflop/s (ms)  MAGMA_tally4 error \n");
    printf("=========================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            K = opts.ksize[itest];
            gflops = FLOPS_ZHERK( K, N ) / 1e9 * batchCount;

            if ( opts.transA == Magma_tally4NoTrans ) {
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
            
            TESTING_MALLOC_CPU( h_A,  magma_tally4DoubleComplex, sizeA );
            TESTING_MALLOC_CPU( h_C,  magma_tally4DoubleComplex, sizeC );
            TESTING_MALLOC_CPU( h_Cmagma_tally4,  magma_tally4DoubleComplex, sizeC  );
            
            TESTING_MALLOC_DEV( d_A, magma_tally4DoubleComplex, ldda*An*batchCount );
            TESTING_MALLOC_DEV( d_C, magma_tally4DoubleComplex, lddc*N*batchCount );

            magma_tally4_malloc((void**)&A_array, batchCount * sizeof(*A_array));
            magma_tally4_malloc((void**)&C_array, batchCount * sizeof(*C_array));

            /* Initialize the matrices */
            lapackf77_zlarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_zlarnv( &ione, ISEED, &sizeC, h_C );
            
            /* =====================================================================
               Performs operation using MAGMA_tally4BLAS
               =================================================================== */
            magma_tally4_zsetmatrix( An, Ak*batchCount, h_A, lda, d_A, ldda );
            magma_tally4_zsetmatrix( N, N*batchCount, h_C, ldc, d_C, lddc );
            
            zset_pointer(A_array, d_A, lda, 0, 0, ldda*Ak, batchCount, queue);
            zset_pointer(C_array, d_C, ldc, 0, 0, lddc*N,  batchCount, queue);

            magma_tally4_time = magma_tally4_sync_wtime( NULL );
            magma_tally4blas_zherk_batched(opts.uplo, opts.transA, N, K,
                             alpha, A_array, ldda,
                             beta,  C_array, lddc, batchCount, queue);
                             
            magma_tally4_time = magma_tally4_sync_wtime( NULL ) - magma_tally4_time;
            magma_tally4_perf = gflops / magma_tally4_time;
            
            magma_tally4_zgetmatrix( N, NN, d_C, lddc, h_Cmagma_tally4, ldc );
            
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally4_wtime();
                for(int i=0; i<batchCount; i++)
                {
                   blasf77_zherk(
                               lapack_uplo_const_tally4(opts.uplo), lapack_trans_const_tally4(opts.transA),
                               &N, &K,
                               &alpha, h_A + i*lda*Ak, &lda,
                               &beta,  h_C + i*ldc*N, &ldc );
                }
                cpu_time = magma_tally4_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.check ) {
                // compute relative error for magma_tally4, relative to lapack,
                // |C_magma_tally4 - C_lapack| / |C_lapack|
                sizeC = ldc*N;
                magma_tally4_error = MAGMA_tally4_D_ZERO;
                for(int i=0; i<batchCount; i++)
                {
                    Cnorm = lapackf77_zlanhe("fro", lapack_uplo_const_tally4(opts.uplo), &N, h_C+i*ldc*N, &ldc, work);
                    blasf77_zaxpy( &sizeC, &c_neg_one, h_C+i*ldc*N, &ione, h_Cmagma_tally4+i*ldc*N, &ione );
                    current_error = lapackf77_zlanhe( "fro", lapack_uplo_const_tally4(opts.uplo), &N, h_Cmagma_tally4+i*ldc*N, &ldc, work ) / Cnorm;
                    if ( isnan(current_error) || isinf(current_error) ) {
                        magma_tally4_error = current_error;
                        break;
                    }
                    magma_tally4_error = max(magma_tally4_error, current_error);
                }

                printf("%5d %5d %5d  %7.2f (%7.2f)   %7.2f (%7.2f)    %8.2e   %s\n",
                       (int) batchCount, (int) N, (int) K,
                       magma_tally4_perf, 1000.*magma_tally4_time,
                       cpu_perf,    1000.*cpu_time,
                       magma_tally4_error, (magma_tally4_error < tol ? "ok" : "failed"));

                status += ! (magma_tally4_error < tol);

            }
            else {
                printf("%5d %5d %5d   %7.2f (%7.2f)    ---   (  ---  )    ---     ---\n",
                       (int) batchCount, (int) N, (int) K,
                       magma_tally4_perf, 1000.*magma_tally4_time);

            }
            
            TESTING_FREE_CPU( h_A  );
            TESTING_FREE_CPU( h_C  );
            TESTING_FREE_CPU( h_Cmagma_tally4  );

            
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
