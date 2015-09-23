/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zgemv_batched.cpp normal z -> s, Fri Jan 30 19:00:26 2015
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
   -- Testing sgemm_batched
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, magma_minproduct_perf, magma_minproduct_time, cpu_perf, cpu_time;
    float          magma_minproduct_error, magma_minproduct_err, Ynorm, work[1];
    magma_minproduct_int_t M, N, Xm, Ym, lda, ldda;
    magma_minproduct_int_t sizeA, sizeX, sizeY;
    magma_minproduct_int_t incx = 1;
    magma_minproduct_int_t incy = 1;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t status = 0;
    magma_minproduct_int_t batchCount;

    float *h_A, *h_X, *h_Y, *h_Ymagma_minproduct;
    float *d_A, *d_X, *d_Y;
    float c_neg_one = MAGMA_minproduct_S_NEG_ONE;
    float alpha = MAGMA_minproduct_S_MAKE(  0.29, -0.86 );
    float beta  = MAGMA_minproduct_S_MAKE( -0.48,  0.38 );
    float **A_array = NULL;
    float **X_array = NULL;
    float **Y_array = NULL;


    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    batchCount = opts.batchcount;
    opts.lapack |= opts.check; 

    //float tol = opts.tolerance * lapackf77_slamch("E");

    printf("trans = %s\n", lapack_trans_const(opts.transA) );

    printf("BatchCount    M     N     MAGMA_minproduct Gflop/s (ms)  CPU Gflop/s (ms)  MAGMA_minproduct error\n");

    printf("===================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            lda    = ((M+31)/32)*32;
            gflops = FLOPS_SGEMV( M, N ) / 1e9 * batchCount;

            if ( opts.transA == Magma_minproductNoTrans ) {
                Xm = N;
                Ym = M;
            } else {
                Xm = M;
                Ym = N;
            }

            sizeA = lda*N*batchCount;
            sizeX = incx*Xm*batchCount;
            sizeY = incy*Ym*batchCount;

            ldda = ((lda+31)/32)*32;

            TESTING_MALLOC_CPU( h_A,  float, sizeA );
            TESTING_MALLOC_CPU( h_X,  float, sizeX );
            TESTING_MALLOC_CPU( h_Y,  float, sizeY  );
            TESTING_MALLOC_CPU( h_Ymagma_minproduct,  float, sizeY  );


            TESTING_MALLOC_DEV( d_A, float, ldda*N*batchCount );
            TESTING_MALLOC_DEV( d_X, float, sizeX );
            TESTING_MALLOC_DEV( d_Y, float, sizeY );

            magma_minproduct_malloc((void**)&A_array, batchCount * sizeof(*A_array));
            magma_minproduct_malloc((void**)&X_array, batchCount * sizeof(*X_array));
            magma_minproduct_malloc((void**)&Y_array, batchCount * sizeof(*Y_array));

            /* Initialize the matrices */
            lapackf77_slarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_slarnv( &ione, ISEED, &sizeX, h_X );
            lapackf77_slarnv( &ione, ISEED, &sizeY, h_Y );
            
            /* =====================================================================
               Performs operation using MAGMA_minproductBLAS
               =================================================================== */
            magma_minproduct_ssetmatrix( M, N*batchCount, h_A, lda, d_A, ldda );
            magma_minproduct_ssetvector( Xm*batchCount, h_X, incx, d_X, incx );
            magma_minproduct_ssetvector( Ym*batchCount, h_Y, incy, d_Y, incy );
            
            sset_pointer(A_array, d_A, ldda, 0, 0, ldda*N, batchCount, magma_minproduct_stream);
            sset_pointer(X_array, d_X, 1, 0, 0, incx*Xm, batchCount, magma_minproduct_stream);
            sset_pointer(Y_array, d_Y, 1, 0, 0, incy*Ym, batchCount, magma_minproduct_stream);

            magma_minproduct_time = magma_minproduct_sync_wtime( NULL );
            magma_minproductblas_sgemv_batched(opts.transA, M, N,
                             alpha, A_array, ldda,
                                    X_array, incx,
                             beta,  Y_array, incy, batchCount, magma_minproduct_stream);
            magma_minproduct_time = magma_minproduct_sync_wtime( NULL ) - magma_minproduct_time;
            magma_minproduct_perf = gflops / magma_minproduct_time;            
            magma_minproduct_sgetvector( Ym*batchCount, d_Y, incy, h_Ymagma_minproduct, incy );
                      
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_minproduct_wtime();
                for(int i=0; i<batchCount; i++)
                {
                   blasf77_sgemv(
                               lapack_trans_const(opts.transA), 
                               &M, &N,
                               &alpha, h_A + i*lda*N, &lda,
                                       h_X + i*Xm, &incx,
                               &beta,  h_Y + i*Ym, &incy );
                }
                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.lapack ) {
                // compute relative error for both magma_minproduct  relative to lapack,
                // |C_magma_minproduct - C_lapack| / |C_lapack|
                magma_minproduct_error = 0.0;

                for(int s=0; s<batchCount; s++)
                {

                    Ynorm = lapackf77_slange( "M", &M, &ione, h_Y + s*Ym, &incy, work );

                    blasf77_saxpy( &Ym, &c_neg_one, h_Y + s*Ym, &ione, h_Ymagma_minproduct + s*Ym, &ione );
                    magma_minproduct_err = lapackf77_slange( "M", &M, &ione, h_Ymagma_minproduct + s*Ym, &incy, work ) / Ynorm; 

                    if ( isnan(magma_minproduct_err) || isinf(magma_minproduct_err) ) {
                      magma_minproduct_error = magma_minproduct_err;
                      break;
                    }
                    magma_minproduct_error = max(fabs(magma_minproduct_err), magma_minproduct_error); 

                }

                    printf("%10d %5d %5d  %7.2f (%7.2f)    %7.2f (%7.2f)   %8.2e  \n",
                       (int) batchCount, (int) M, (int) N, 
                       magma_minproduct_perf,  1000.*magma_minproduct_time,
                       cpu_perf,    1000.*cpu_time,
                       magma_minproduct_error);
            }
            else {

                    printf("%10d %5d %5d  %7.2f (%7.2f)    ---   (  ---  )    ---\n",
                       (int) batchCount, (int) M, (int) N, 
                       magma_minproduct_perf,  1000.*magma_minproduct_time);
            }
            
            TESTING_FREE_CPU( h_A  );
            TESTING_FREE_CPU( h_X  );
            TESTING_FREE_CPU( h_Y  );
            TESTING_FREE_CPU( h_Ymagma_minproduct  );


            TESTING_FREE_DEV( d_A );
            TESTING_FREE_DEV( d_X );
            TESTING_FREE_DEV( d_Y );
            TESTING_FREE_DEV( A_array );
            TESTING_FREE_DEV( X_array );
            TESTING_FREE_DEV( Y_array );

            
            fflush( stdout);

        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
