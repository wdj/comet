/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zherk_batched.cpp normal z -> d, Fri Jan 30 19:00:26 2015
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
#include "magma_tally2.h"
#include "magma_tally2_lapack.h"
#include "testings.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dsyrk_batched
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, magma_tally2_perf, magma_tally2_time, cpu_perf=0., cpu_time=0.;
    double          current_error, magma_tally2_error, Cnorm, work[1];
    magma_tally2_int_t N, K;
    magma_tally2_int_t Ak, An;
    magma_tally2_int_t sizeA, sizeC;
    magma_tally2_int_t lda, ldc, ldda, lddc;
    magma_tally2_int_t ione     = 1;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    magma_tally2_int_t NN;
    magma_tally2_int_t batchCount;
 
    double *h_A, *h_C, *h_Cmagma_tally2;
    double *d_A, *d_C;
    double c_neg_one = MAGMA_tally2_D_NEG_ONE;
    double alpha = 0.29;
    double beta  = -0.48;
    double **A_array = NULL;
    double **C_array = NULL;
    magma_tally2_int_t status = 0;

    magma_tally2_queue_t queue = magma_tally2_stream;
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    batchCount = opts.batchcount;

    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("If running lapack (option --lapack), MAGMA_tally2 error is computed\n"
           "relative to CPU BLAS result.\n\n");
    printf("uplo = %s, transA = %s\n",
           lapack_uplo_const_tally2(opts.uplo), lapack_trans_const_tally2(opts.transA) );
    printf(" BatchCount    N     K   MAGMA_tally2 Gflop/s (ms)  CPU Gflop/s (ms)  MAGMA_tally2 error \n");
    printf("=========================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            K = opts.ksize[itest];
            gflops = FLOPS_DSYRK( K, N ) / 1e9 * batchCount;

            if ( opts.transA == Magma_tally2NoTrans ) {
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
            
            TESTING_MALLOC_CPU( h_A,  double, sizeA );
            TESTING_MALLOC_CPU( h_C,  double, sizeC );
            TESTING_MALLOC_CPU( h_Cmagma_tally2,  double, sizeC  );
            
            TESTING_MALLOC_DEV( d_A, double, ldda*An*batchCount );
            TESTING_MALLOC_DEV( d_C, double, lddc*N*batchCount );

            magma_tally2_malloc((void**)&A_array, batchCount * sizeof(*A_array));
            magma_tally2_malloc((void**)&C_array, batchCount * sizeof(*C_array));

            /* Initialize the matrices */
            lapackf77_dlarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_dlarnv( &ione, ISEED, &sizeC, h_C );
            
            /* =====================================================================
               Performs operation using MAGMA_tally2BLAS
               =================================================================== */
            magma_tally2_dsetmatrix( An, Ak*batchCount, h_A, lda, d_A, ldda );
            magma_tally2_dsetmatrix( N, N*batchCount, h_C, ldc, d_C, lddc );
            
            dset_pointer(A_array, d_A, lda, 0, 0, ldda*Ak, batchCount, queue);
            dset_pointer(C_array, d_C, ldc, 0, 0, lddc*N,  batchCount, queue);

            magma_tally2_time = magma_tally2_sync_wtime( NULL );
            magma_tally2blas_dsyrk_batched(opts.uplo, opts.transA, N, K,
                             alpha, A_array, ldda,
                             beta,  C_array, lddc, batchCount, queue);
                             
            magma_tally2_time = magma_tally2_sync_wtime( NULL ) - magma_tally2_time;
            magma_tally2_perf = gflops / magma_tally2_time;
            
            magma_tally2_dgetmatrix( N, NN, d_C, lddc, h_Cmagma_tally2, ldc );
            
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally2_wtime();
                for(int i=0; i<batchCount; i++)
                {
                   blasf77_dsyrk(
                               lapack_uplo_const_tally2(opts.uplo), lapack_trans_const_tally2(opts.transA),
                               &N, &K,
                               &alpha, h_A + i*lda*Ak, &lda,
                               &beta,  h_C + i*ldc*N, &ldc );
                }
                cpu_time = magma_tally2_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.check ) {
                // compute relative error for magma_tally2, relative to lapack,
                // |C_magma_tally2 - C_lapack| / |C_lapack|
                sizeC = ldc*N;
                magma_tally2_error = MAGMA_tally2_D_ZERO;
                for(int i=0; i<batchCount; i++)
                {
                    Cnorm = lapackf77_dlansy("fro", lapack_uplo_const_tally2(opts.uplo), &N, h_C+i*ldc*N, &ldc, work);
                    blasf77_daxpy( &sizeC, &c_neg_one, h_C+i*ldc*N, &ione, h_Cmagma_tally2+i*ldc*N, &ione );
                    current_error = lapackf77_dlansy( "fro", lapack_uplo_const_tally2(opts.uplo), &N, h_Cmagma_tally2+i*ldc*N, &ldc, work ) / Cnorm;
                    if ( isnan(current_error) || isinf(current_error) ) {
                        magma_tally2_error = current_error;
                        break;
                    }
                    magma_tally2_error = max(magma_tally2_error, current_error);
                }

                printf("%5d %5d %5d  %7.2f (%7.2f)   %7.2f (%7.2f)    %8.2e   %s\n",
                       (int) batchCount, (int) N, (int) K,
                       magma_tally2_perf, 1000.*magma_tally2_time,
                       cpu_perf,    1000.*cpu_time,
                       magma_tally2_error, (magma_tally2_error < tol ? "ok" : "failed"));

                status += ! (magma_tally2_error < tol);

            }
            else {
                printf("%5d %5d %5d   %7.2f (%7.2f)    ---   (  ---  )    ---     ---\n",
                       (int) batchCount, (int) N, (int) K,
                       magma_tally2_perf, 1000.*magma_tally2_time);

            }
            
            TESTING_FREE_CPU( h_A  );
            TESTING_FREE_CPU( h_C  );
            TESTING_FREE_CPU( h_Cmagma_tally2  );

            
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
