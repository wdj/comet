/*
    -- MAGMA_minproduct (version 1.6.1) --
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
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "testings.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zherk_batched
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, magma_minproduct_perf, magma_minproduct_time, cpu_perf=0., cpu_time=0.;
    double          current_error, magma_minproduct_error, Cnorm, work[1];
    magma_minproduct_int_t N, K;
    magma_minproduct_int_t Ak, An;
    magma_minproduct_int_t sizeA, sizeC;
    magma_minproduct_int_t lda, ldc, ldda, lddc;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t NN;
    magma_minproduct_int_t batchCount;
 
    magma_minproductDoubleComplex *h_A, *h_C, *h_Cmagma_minproduct;
    magma_minproductDoubleComplex *d_A, *d_C;
    magma_minproductDoubleComplex c_neg_one = MAGMA_minproduct_Z_NEG_ONE;
    double alpha = 0.29;
    double beta  = -0.48;
    magma_minproductDoubleComplex **A_array = NULL;
    magma_minproductDoubleComplex **C_array = NULL;
    magma_minproduct_int_t status = 0;

    magma_minproduct_queue_t queue = magma_minproduct_stream;
    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    batchCount = opts.batchcount;

    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("If running lapack (option --lapack), MAGMA_minproduct error is computed\n"
           "relative to CPU BLAS result.\n\n");
    printf("uplo = %s, transA = %s\n",
           lapack_uplo_const(opts.uplo), lapack_trans_const(opts.transA) );
    printf(" BatchCount    N     K   MAGMA_minproduct Gflop/s (ms)  CPU Gflop/s (ms)  MAGMA_minproduct error \n");
    printf("=========================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            K = opts.ksize[itest];
            gflops = FLOPS_ZHERK( K, N ) / 1e9 * batchCount;

            if ( opts.transA == Magma_minproductNoTrans ) {
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
            
            TESTING_MALLOC_CPU( h_A,  magma_minproductDoubleComplex, sizeA );
            TESTING_MALLOC_CPU( h_C,  magma_minproductDoubleComplex, sizeC );
            TESTING_MALLOC_CPU( h_Cmagma_minproduct,  magma_minproductDoubleComplex, sizeC  );
            
            TESTING_MALLOC_DEV( d_A, magma_minproductDoubleComplex, ldda*An*batchCount );
            TESTING_MALLOC_DEV( d_C, magma_minproductDoubleComplex, lddc*N*batchCount );

            magma_minproduct_malloc((void**)&A_array, batchCount * sizeof(*A_array));
            magma_minproduct_malloc((void**)&C_array, batchCount * sizeof(*C_array));

            /* Initialize the matrices */
            lapackf77_zlarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_zlarnv( &ione, ISEED, &sizeC, h_C );
            
            /* =====================================================================
               Performs operation using MAGMA_minproductBLAS
               =================================================================== */
            magma_minproduct_zsetmatrix( An, Ak*batchCount, h_A, lda, d_A, ldda );
            magma_minproduct_zsetmatrix( N, N*batchCount, h_C, ldc, d_C, lddc );
            
            zset_pointer(A_array, d_A, lda, 0, 0, ldda*Ak, batchCount, queue);
            zset_pointer(C_array, d_C, ldc, 0, 0, lddc*N,  batchCount, queue);

            magma_minproduct_time = magma_minproduct_sync_wtime( NULL );
            magma_minproductblas_zherk_batched(opts.uplo, opts.transA, N, K,
                             alpha, A_array, ldda,
                             beta,  C_array, lddc, batchCount, queue);
                             
            magma_minproduct_time = magma_minproduct_sync_wtime( NULL ) - magma_minproduct_time;
            magma_minproduct_perf = gflops / magma_minproduct_time;
            
            magma_minproduct_zgetmatrix( N, NN, d_C, lddc, h_Cmagma_minproduct, ldc );
            
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_minproduct_wtime();
                for(int i=0; i<batchCount; i++)
                {
                   blasf77_zherk(
                               lapack_uplo_const(opts.uplo), lapack_trans_const(opts.transA),
                               &N, &K,
                               &alpha, h_A + i*lda*Ak, &lda,
                               &beta,  h_C + i*ldc*N, &ldc );
                }
                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.check ) {
                // compute relative error for magma_minproduct, relative to lapack,
                // |C_magma_minproduct - C_lapack| / |C_lapack|
                sizeC = ldc*N;
                magma_minproduct_error = MAGMA_minproduct_D_ZERO;
                for(int i=0; i<batchCount; i++)
                {
                    Cnorm = lapackf77_zlanhe("fro", lapack_uplo_const(opts.uplo), &N, h_C+i*ldc*N, &ldc, work);
                    blasf77_zaxpy( &sizeC, &c_neg_one, h_C+i*ldc*N, &ione, h_Cmagma_minproduct+i*ldc*N, &ione );
                    current_error = lapackf77_zlanhe( "fro", lapack_uplo_const(opts.uplo), &N, h_Cmagma_minproduct+i*ldc*N, &ldc, work ) / Cnorm;
                    if ( isnan(current_error) || isinf(current_error) ) {
                        magma_minproduct_error = current_error;
                        break;
                    }
                    magma_minproduct_error = max(magma_minproduct_error, current_error);
                }

                printf("%5d %5d %5d  %7.2f (%7.2f)   %7.2f (%7.2f)    %8.2e   %s\n",
                       (int) batchCount, (int) N, (int) K,
                       magma_minproduct_perf, 1000.*magma_minproduct_time,
                       cpu_perf,    1000.*cpu_time,
                       magma_minproduct_error, (magma_minproduct_error < tol ? "ok" : "failed"));

                status += ! (magma_minproduct_error < tol);

            }
            else {
                printf("%5d %5d %5d   %7.2f (%7.2f)    ---   (  ---  )    ---     ---\n",
                       (int) batchCount, (int) N, (int) K,
                       magma_minproduct_perf, 1000.*magma_minproduct_time);

            }
            
            TESTING_FREE_CPU( h_A  );
            TESTING_FREE_CPU( h_C  );
            TESTING_FREE_CPU( h_Cmagma_minproduct  );

            
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
