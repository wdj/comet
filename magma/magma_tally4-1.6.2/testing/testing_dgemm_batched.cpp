/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zgemm_batched.cpp normal z -> d, Fri Jan 30 19:00:26 2015
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
#include "testings.h"  // before magma_tally4.h, to include cublas_v2
#include "flops.h"
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dgemm_batched
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, magma_tally4_perf, magma_tally4_time, cublas_perf, cublas_time, cpu_perf, cpu_time;
    double          magma_tally4_error, cublas_error, magma_tally4_err, cublas_err, Cnorm, work[1];
    magma_tally4_int_t M, N, K;
    magma_tally4_int_t Am, An, Bm, Bn;
    magma_tally4_int_t sizeA, sizeB, sizeC;
    magma_tally4_int_t lda, ldb, ldc, ldda, lddb, lddc;
    magma_tally4_int_t ione     = 1;
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    magma_tally4_int_t status = 0;
    magma_tally4_int_t NN;
    magma_tally4_int_t batchCount;

    double *h_A, *h_B, *h_C, *h_Cmagma_tally4, *h_Ccublas;
    double *d_A, *d_B, *d_C;
    double c_neg_one = MAGMA_tally4_D_NEG_ONE;
    double alpha = MAGMA_tally4_D_MAKE(  0.29, -0.86 );
    double beta  = MAGMA_tally4_D_MAKE( -0.48,  0.38 );
    double **A_array = NULL;
    double **B_array = NULL;
    double **C_array = NULL;

    magma_tally4_queue_t queue = magma_tally4_stream;
    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );
    batchCount = opts.batchcount;
    cublasHandle_t handle = opts.handle;

    //double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("If running lapack (option --lapack), MAGMA_tally4 and CUBLAS error are both computed\n"
           "relative to CPU BLAS result. Else, MAGMA_tally4 error is computed relative to CUBLAS result.\n\n"
           "transA = %s, transB = %s\n", 
           lapack_trans_const_tally4(opts.transA),
           lapack_trans_const_tally4(opts.transB));
    printf("BatchCount    M     N     K   MAGMA_tally4 Gflop/s (ms)  CUBLAS Gflop/s (ms)  CPU Gflop/s (ms)  MAGMA_tally4 error  CUBLAS error\n");
    printf("=========================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            K = opts.ksize[itest];
            gflops = FLOPS_DGEMM( M, N, K ) / 1e9 * batchCount;

            if ( opts.transA == Magma_tally4NoTrans ) {
                lda = Am = M;
                An = K;
            } else {
                lda = Am = K;
                An = M;
            }
            
            if ( opts.transB == Magma_tally4NoTrans ) {
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
            
            TESTING_MALLOC_CPU( h_A,  double, sizeA );
            TESTING_MALLOC_CPU( h_B,  double, sizeB );
            TESTING_MALLOC_CPU( h_C,  double, sizeC  );
            TESTING_MALLOC_CPU( h_Cmagma_tally4,  double, sizeC  );
            TESTING_MALLOC_CPU( h_Ccublas, double, sizeC  );

            TESTING_MALLOC_DEV( d_A, double, ldda*An*batchCount );
            TESTING_MALLOC_DEV( d_B, double, lddb*Bn*batchCount );
            TESTING_MALLOC_DEV( d_C, double, lddc*N*batchCount  );

            magma_tally4_malloc((void**)&A_array, batchCount * sizeof(*A_array));
            magma_tally4_malloc((void**)&B_array, batchCount * sizeof(*B_array));
            magma_tally4_malloc((void**)&C_array, batchCount * sizeof(*C_array));

            /* Initialize the matrices */
            lapackf77_dlarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_dlarnv( &ione, ISEED, &sizeB, h_B );
            lapackf77_dlarnv( &ione, ISEED, &sizeC, h_C );
            
            /* =====================================================================
               Performs operation using MAGMA_tally4BLAS
               =================================================================== */
            magma_tally4_dsetmatrix( Am, An*batchCount, h_A, lda, d_A, ldda );
            magma_tally4_dsetmatrix( Bm, Bn*batchCount, h_B, ldb, d_B, lddb );
            magma_tally4_dsetmatrix( M, N*batchCount, h_C, ldc, d_C, lddc );
            
            dset_pointer(A_array, d_A, ldda, 0, 0, ldda*An, batchCount, queue);
            dset_pointer(B_array, d_B, lddb, 0, 0, lddb*Bn, batchCount, queue);
            dset_pointer(C_array, d_C, lddc, 0, 0, lddc*N,  batchCount, queue);

            magma_tally4_time = magma_tally4_sync_wtime( NULL );
            magma_tally4blas_dgemm_batched(opts.transA, opts.transB, M, N, K,
                             alpha, A_array, ldda,
                                    B_array, lddb,
                             beta,  C_array, lddc, batchCount, queue);
            magma_tally4_time = magma_tally4_sync_wtime( NULL ) - magma_tally4_time;
            magma_tally4_perf = gflops / magma_tally4_time;            
            magma_tally4_dgetmatrix( M, N*batchCount, d_C, lddc, h_Cmagma_tally4, ldc );
            
            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */

            magma_tally4_dsetmatrix( M, N*batchCount, h_C, ldc, d_C, lddc );
            
            cublas_time = magma_tally4_sync_wtime( NULL );

            cublasDgemmBatched(handle, cublas_trans_const_tally4(opts.transA), cublas_trans_const_tally4(opts.transB), M, N, K,
                               &alpha, (const double**) A_array, ldda,
                               (const double**) B_array, lddb,
                               &beta,  C_array, lddc, batchCount );

            cublas_time = magma_tally4_sync_wtime( NULL ) - cublas_time;
            cublas_perf = gflops / cublas_time;
            
            magma_tally4_dgetmatrix( M, N*batchCount, d_C, lddc, h_Ccublas, ldc );
          
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally4_wtime();
                for(int i=0; i<batchCount; i++)
                {
                   blasf77_dgemm(
                               lapack_trans_const_tally4(opts.transA), lapack_trans_const_tally4(opts.transB),
                               &M, &N, &K,
                               &alpha, h_A + i*lda*An, &lda,
                                       h_B + i*ldb*Bn, &ldb,
                               &beta,  h_C + i*ldc*N, &ldc );
                }
                cpu_time = magma_tally4_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.lapack ) {
                // compute relative error for both magma_tally4 & cublas, relative to lapack,
                // |C_magma_tally4 - C_lapack| / |C_lapack|
                magma_tally4_error = 0.0;
                cublas_error = 0.0;

                for(int s=0; s<batchCount; s++)
                {
                    magma_tally4_int_t C_batchSize = ldc * N;
 
                    Cnorm = lapackf77_dlange( "M", &M, &N, h_C + s*C_batchSize, &ldc, work );

                    blasf77_daxpy( &C_batchSize, &c_neg_one, h_C + s*C_batchSize, &ione, h_Cmagma_tally4 + s*C_batchSize, &ione );
                    magma_tally4_err = lapackf77_dlange( "M", &M, &N, h_Cmagma_tally4 + s*C_batchSize, &ldc, work ) / Cnorm; 

                    if ( isnan(magma_tally4_err) || isinf(magma_tally4_err) ) {
                      magma_tally4_error = magma_tally4_err;
                      break;
                    }
                    magma_tally4_error = max(fabs(magma_tally4_err), magma_tally4_error); 

                    blasf77_daxpy( &C_batchSize, &c_neg_one, h_C + s*C_batchSize, &ione, h_Ccublas + s*C_batchSize, &ione );
                    cublas_err = lapackf77_dlange( "M", &M, &N, h_Ccublas + s*C_batchSize, &ldc, work ) / Cnorm; 
                    
                   if ( isnan(cublas_err) || isinf(cublas_err) ) {
                      cublas_error = cublas_err;
                      break;
                    }
                    cublas_error = max(fabs(cublas_err), cublas_error); 

                }

                    printf("%10d %5d %5d %5d  %7.2f (%7.2f)    %7.2f (%7.2f)   %7.2f (%7.2f)      %8.2e     %8.2e  \n",
                       (int) batchCount, (int) M, (int) N, (int) K, 
                       magma_tally4_perf,  1000.*magma_tally4_time,
                       cublas_perf, 1000.*cublas_time,
                       cpu_perf,    1000.*cpu_time,
                       magma_tally4_error, cublas_error);
            }
            else {
                // compute relative error for magma_tally4, relative to cublas

                    Cnorm = lapackf77_dlange( "M", &M, &NN, h_Ccublas, &ldc, work );
                    blasf77_daxpy( &sizeC, &c_neg_one, h_Ccublas, &ione, h_Cmagma_tally4, &ione );
                    magma_tally4_error = lapackf77_dlange( "M", &M, &NN, h_Cmagma_tally4, &ldc, work ) / Cnorm;

                    printf("%10d %5d %5d %5d  %7.2f (%7.2f)    %7.2f (%7.2f)   ---   (  ---  )    %8.2e     ---\n",
                       (int) batchCount, (int) M, (int) N, (int) K,
                       magma_tally4_perf,  1000.*magma_tally4_time,
                       cublas_perf, 1000.*cublas_time,
                       magma_tally4_error );
            }
            
            TESTING_FREE_CPU( h_A  );
            TESTING_FREE_CPU( h_B  );
            TESTING_FREE_CPU( h_C  );
            TESTING_FREE_CPU( h_Cmagma_tally4  );
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
