/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> c d s
       @author Mark Gates
       @author Azzam Haidar
       @author Tingxing Dong
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "testings.h"  // before magma_minproduct.h, to include cublas_v2
#include "flops.h"
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zgemm_batched
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, magma_minproduct_perf, magma_minproduct_time, cublas_perf, cublas_time, cpu_perf, cpu_time;
    double          magma_minproduct_error, cublas_error, magma_minproduct_err, cublas_err, Cnorm, work[1];
    magma_minproduct_int_t M, N, K;
    magma_minproduct_int_t Am, An, Bm, Bn;
    magma_minproduct_int_t sizeA, sizeB, sizeC;
    magma_minproduct_int_t lda, ldb, ldc, ldda, lddb, lddc;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t status = 0;
    magma_minproduct_int_t NN;
    magma_minproduct_int_t batchCount;

    magma_minproductDoubleComplex *h_A, *h_B, *h_C, *h_Cmagma_minproduct, *h_Ccublas;
    magma_minproductDoubleComplex *d_A, *d_B, *d_C;
    magma_minproductDoubleComplex c_neg_one = MAGMA_minproduct_Z_NEG_ONE;
    magma_minproductDoubleComplex alpha = MAGMA_minproduct_Z_MAKE(  0.29, -0.86 );
    magma_minproductDoubleComplex beta  = MAGMA_minproduct_Z_MAKE( -0.48,  0.38 );
    magma_minproductDoubleComplex **A_array = NULL;
    magma_minproductDoubleComplex **B_array = NULL;
    magma_minproductDoubleComplex **C_array = NULL;

    magma_minproduct_queue_t queue = magma_minproduct_stream;
    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    batchCount = opts.batchcount;
    cublasHandle_t handle = opts.handle;

    //double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("If running lapack (option --lapack), MAGMA_minproduct and CUBLAS error are both computed\n"
           "relative to CPU BLAS result. Else, MAGMA_minproduct error is computed relative to CUBLAS result.\n\n"
           "transA = %s, transB = %s\n", 
           lapack_trans_const(opts.transA),
           lapack_trans_const(opts.transB));
    printf("BatchCount    M     N     K   MAGMA_minproduct Gflop/s (ms)  CUBLAS Gflop/s (ms)  CPU Gflop/s (ms)  MAGMA_minproduct error  CUBLAS error\n");
    printf("=========================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            K = opts.ksize[itest];
            gflops = FLOPS_ZGEMM( M, N, K ) / 1e9 * batchCount;

            if ( opts.transA == Magma_minproductNoTrans ) {
                lda = Am = M;
                An = K;
            } else {
                lda = Am = K;
                An = M;
            }
            
            if ( opts.transB == Magma_minproductNoTrans ) {
                ldb = Bm = K;
                Bn = N;
            } else {
                ldb = Bm = N;
                Bn = K;
            }
            ldc = M;
            
            NN = N * batchCount;

            ldda = ((lda+31)/32)*32;
            lddb = ((ldb+31)/32)*32;
            lddc = ((ldc+31)/32)*32;

            sizeA = lda*An*batchCount;
            sizeB = ldb*Bn*batchCount;
            sizeC = ldc*N*batchCount;
            
            TESTING_MALLOC_CPU( h_A,  magma_minproductDoubleComplex, sizeA );
            TESTING_MALLOC_CPU( h_B,  magma_minproductDoubleComplex, sizeB );
            TESTING_MALLOC_CPU( h_C,  magma_minproductDoubleComplex, sizeC  );
            TESTING_MALLOC_CPU( h_Cmagma_minproduct,  magma_minproductDoubleComplex, sizeC  );
            TESTING_MALLOC_CPU( h_Ccublas, magma_minproductDoubleComplex, sizeC  );

            TESTING_MALLOC_DEV( d_A, magma_minproductDoubleComplex, ldda*An*batchCount );
            TESTING_MALLOC_DEV( d_B, magma_minproductDoubleComplex, lddb*Bn*batchCount );
            TESTING_MALLOC_DEV( d_C, magma_minproductDoubleComplex, lddc*N*batchCount  );

            magma_minproduct_malloc((void**)&A_array, batchCount * sizeof(*A_array));
            magma_minproduct_malloc((void**)&B_array, batchCount * sizeof(*B_array));
            magma_minproduct_malloc((void**)&C_array, batchCount * sizeof(*C_array));

            /* Initialize the matrices */
            lapackf77_zlarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_zlarnv( &ione, ISEED, &sizeB, h_B );
            lapackf77_zlarnv( &ione, ISEED, &sizeC, h_C );
            
            /* =====================================================================
               Performs operation using MAGMA_minproductBLAS
               =================================================================== */
            magma_minproduct_zsetmatrix( Am, An*batchCount, h_A, lda, d_A, ldda );
            magma_minproduct_zsetmatrix( Bm, Bn*batchCount, h_B, ldb, d_B, lddb );
            magma_minproduct_zsetmatrix( M, N*batchCount, h_C, ldc, d_C, lddc );
            
            zset_pointer(A_array, d_A, ldda, 0, 0, ldda*An, batchCount, queue);
            zset_pointer(B_array, d_B, lddb, 0, 0, lddb*Bn, batchCount, queue);
            zset_pointer(C_array, d_C, lddc, 0, 0, lddc*N,  batchCount, queue);

            magma_minproduct_time = magma_minproduct_sync_wtime( NULL );
            magma_minproductblas_zgemm_batched(opts.transA, opts.transB, M, N, K,
                             alpha, A_array, ldda,
                                    B_array, lddb,
                             beta,  C_array, lddc, batchCount, queue);
            magma_minproduct_time = magma_minproduct_sync_wtime( NULL ) - magma_minproduct_time;
            magma_minproduct_perf = gflops / magma_minproduct_time;            
            magma_minproduct_zgetmatrix( M, N*batchCount, d_C, lddc, h_Cmagma_minproduct, ldc );
            
            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */

            magma_minproduct_zsetmatrix( M, N*batchCount, h_C, ldc, d_C, lddc );
            
            cublas_time = magma_minproduct_sync_wtime( NULL );

            cublasZgemmBatched(handle, cublas_trans_const(opts.transA), cublas_trans_const(opts.transB), M, N, K,
                               &alpha, (const magma_minproductDoubleComplex**) A_array, ldda,
                               (const magma_minproductDoubleComplex**) B_array, lddb,
                               &beta,  C_array, lddc, batchCount );

            cublas_time = magma_minproduct_sync_wtime( NULL ) - cublas_time;
            cublas_perf = gflops / cublas_time;
            
            magma_minproduct_zgetmatrix( M, N*batchCount, d_C, lddc, h_Ccublas, ldc );
          
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_minproduct_wtime();
                for(int i=0; i<batchCount; i++)
                {
                   blasf77_zgemm(
                               lapack_trans_const(opts.transA), lapack_trans_const(opts.transB),
                               &M, &N, &K,
                               &alpha, h_A + i*lda*An, &lda,
                                       h_B + i*ldb*Bn, &ldb,
                               &beta,  h_C + i*ldc*N, &ldc );
                }
                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.lapack ) {
                // compute relative error for both magma_minproduct & cublas, relative to lapack,
                // |C_magma_minproduct - C_lapack| / |C_lapack|
                magma_minproduct_error = 0.0;
                cublas_error = 0.0;

                for(int s=0; s<batchCount; s++)
                {
                    magma_minproduct_int_t C_batchSize = ldc * N;
 
                    Cnorm = lapackf77_zlange( "M", &M, &N, h_C + s*C_batchSize, &ldc, work );

                    blasf77_zaxpy( &C_batchSize, &c_neg_one, h_C + s*C_batchSize, &ione, h_Cmagma_minproduct + s*C_batchSize, &ione );
                    magma_minproduct_err = lapackf77_zlange( "M", &M, &N, h_Cmagma_minproduct + s*C_batchSize, &ldc, work ) / Cnorm; 

                    if ( isnan(magma_minproduct_err) || isinf(magma_minproduct_err) ) {
                      magma_minproduct_error = magma_minproduct_err;
                      break;
                    }
                    magma_minproduct_error = max(fabs(magma_minproduct_err), magma_minproduct_error); 

                    blasf77_zaxpy( &C_batchSize, &c_neg_one, h_C + s*C_batchSize, &ione, h_Ccublas + s*C_batchSize, &ione );
                    cublas_err = lapackf77_zlange( "M", &M, &N, h_Ccublas + s*C_batchSize, &ldc, work ) / Cnorm; 
                    
                   if ( isnan(cublas_err) || isinf(cublas_err) ) {
                      cublas_error = cublas_err;
                      break;
                    }
                    cublas_error = max(fabs(cublas_err), cublas_error); 

                }

                    printf("%10d %5d %5d %5d  %7.2f (%7.2f)    %7.2f (%7.2f)   %7.2f (%7.2f)      %8.2e     %8.2e  \n",
                       (int) batchCount, (int) M, (int) N, (int) K, 
                       magma_minproduct_perf,  1000.*magma_minproduct_time,
                       cublas_perf, 1000.*cublas_time,
                       cpu_perf,    1000.*cpu_time,
                       magma_minproduct_error, cublas_error);
            }
            else {
                // compute relative error for magma_minproduct, relative to cublas

                    Cnorm = lapackf77_zlange( "M", &M, &NN, h_Ccublas, &ldc, work );
                    blasf77_zaxpy( &sizeC, &c_neg_one, h_Ccublas, &ione, h_Cmagma_minproduct, &ione );
                    magma_minproduct_error = lapackf77_zlange( "M", &M, &NN, h_Cmagma_minproduct, &ldc, work ) / Cnorm;

                    printf("%10d %5d %5d %5d  %7.2f (%7.2f)    %7.2f (%7.2f)   ---   (  ---  )    %8.2e     ---\n",
                       (int) batchCount, (int) M, (int) N, (int) K,
                       magma_minproduct_perf,  1000.*magma_minproduct_time,
                       cublas_perf, 1000.*cublas_time,
                       magma_minproduct_error );
            }
            
            TESTING_FREE_CPU( h_A  );
            TESTING_FREE_CPU( h_B  );
            TESTING_FREE_CPU( h_C  );
            TESTING_FREE_CPU( h_Cmagma_minproduct  );
            TESTING_FREE_CPU( h_Ccublas );

            TESTING_FREE_DEV( d_A );
            TESTING_FREE_DEV( d_B );
            TESTING_FREE_DEV( d_C );
            TESTING_FREE_DEV( A_array );
            TESTING_FREE_DEV( B_array );
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
