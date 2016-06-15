/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zgemv_batched.cpp normal z -> c, Fri Jan 30 19:00:26 2015
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
#include "testings.h"  // before magma_tally3.h, to include cublas_v2
#include "flops.h"
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cgemm_batched
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, magma_tally3_perf, magma_tally3_time, cpu_perf, cpu_time;
    float          magma_tally3_error, magma_tally3_err, Ynorm, work[1];
    magma_tally3_int_t M, N, Xm, Ym, lda, ldda;
    magma_tally3_int_t sizeA, sizeX, sizeY;
    magma_tally3_int_t incx = 1;
    magma_tally3_int_t incy = 1;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t status = 0;
    magma_tally3_int_t batchCount;

    magma_tally3FloatComplex *h_A, *h_X, *h_Y, *h_Ymagma_tally3;
    magma_tally3FloatComplex *d_A, *d_X, *d_Y;
    magma_tally3FloatComplex c_neg_one = MAGMA_tally3_C_NEG_ONE;
    magma_tally3FloatComplex alpha = MAGMA_tally3_C_MAKE(  0.29, -0.86 );
    magma_tally3FloatComplex beta  = MAGMA_tally3_C_MAKE( -0.48,  0.38 );
    magma_tally3FloatComplex **A_array = NULL;
    magma_tally3FloatComplex **X_array = NULL;
    magma_tally3FloatComplex **Y_array = NULL;


    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    batchCount = opts.batchcount;
    opts.lapack |= opts.check; 

    //float tol = opts.tolerance * lapackf77_slamch("E");

    printf("trans = %s\n", lapack_trans_const_tally3(opts.transA) );

    printf("BatchCount    M     N     MAGMA_tally3 Gflop/s (ms)  CPU Gflop/s (ms)  MAGMA_tally3 error\n");

    printf("===================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            lda    = ((M+31)/32)*32;
            gflops = FLOPS_CGEMV( M, N ) / 1e9 * batchCount;

            if ( opts.transA == Magma_tally3NoTrans ) {
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

            TESTING_MALLOC_CPU( h_A,  magma_tally3FloatComplex, sizeA );
            TESTING_MALLOC_CPU( h_X,  magma_tally3FloatComplex, sizeX );
            TESTING_MALLOC_CPU( h_Y,  magma_tally3FloatComplex, sizeY  );
            TESTING_MALLOC_CPU( h_Ymagma_tally3,  magma_tally3FloatComplex, sizeY  );


            TESTING_MALLOC_DEV( d_A, magma_tally3FloatComplex, ldda*N*batchCount );
            TESTING_MALLOC_DEV( d_X, magma_tally3FloatComplex, sizeX );
            TESTING_MALLOC_DEV( d_Y, magma_tally3FloatComplex, sizeY );

            magma_tally3_malloc((void**)&A_array, batchCount * sizeof(*A_array));
            magma_tally3_malloc((void**)&X_array, batchCount * sizeof(*X_array));
            magma_tally3_malloc((void**)&Y_array, batchCount * sizeof(*Y_array));

            /* Initialize the matrices */
            lapackf77_clarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_clarnv( &ione, ISEED, &sizeX, h_X );
            lapackf77_clarnv( &ione, ISEED, &sizeY, h_Y );
            
            /* =====================================================================
               Performs operation using MAGMA_tally3BLAS
               =================================================================== */
            magma_tally3_csetmatrix( M, N*batchCount, h_A, lda, d_A, ldda );
            magma_tally3_csetvector( Xm*batchCount, h_X, incx, d_X, incx );
            magma_tally3_csetvector( Ym*batchCount, h_Y, incy, d_Y, incy );
            
            cset_pointer(A_array, d_A, ldda, 0, 0, ldda*N, batchCount, magma_tally3_stream);
            cset_pointer(X_array, d_X, 1, 0, 0, incx*Xm, batchCount, magma_tally3_stream);
            cset_pointer(Y_array, d_Y, 1, 0, 0, incy*Ym, batchCount, magma_tally3_stream);

            magma_tally3_time = magma_tally3_sync_wtime( NULL );
            magma_tally3blas_cgemv_batched(opts.transA, M, N,
                             alpha, A_array, ldda,
                                    X_array, incx,
                             beta,  Y_array, incy, batchCount, magma_tally3_stream);
            magma_tally3_time = magma_tally3_sync_wtime( NULL ) - magma_tally3_time;
            magma_tally3_perf = gflops / magma_tally3_time;            
            magma_tally3_cgetvector( Ym*batchCount, d_Y, incy, h_Ymagma_tally3, incy );
                      
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally3_wtime();
                for(int i=0; i<batchCount; i++)
                {
                   blasf77_cgemv(
                               lapack_trans_const_tally3(opts.transA), 
                               &M, &N,
                               &alpha, h_A + i*lda*N, &lda,
                                       h_X + i*Xm, &incx,
                               &beta,  h_Y + i*Ym, &incy );
                }
                cpu_time = magma_tally3_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.lapack ) {
                // compute relative error for both magma_tally3  relative to lapack,
                // |C_magma_tally3 - C_lapack| / |C_lapack|
                magma_tally3_error = 0.0;

                for(int s=0; s<batchCount; s++)
                {

                    Ynorm = lapackf77_clange( "M", &M, &ione, h_Y + s*Ym, &incy, work );

                    blasf77_caxpy( &Ym, &c_neg_one, h_Y + s*Ym, &ione, h_Ymagma_tally3 + s*Ym, &ione );
                    magma_tally3_err = lapackf77_clange( "M", &M, &ione, h_Ymagma_tally3 + s*Ym, &incy, work ) / Ynorm; 

                    if ( isnan(magma_tally3_err) || isinf(magma_tally3_err) ) {
                      magma_tally3_error = magma_tally3_err;
                      break;
                    }
                    magma_tally3_error = max(fabs(magma_tally3_err), magma_tally3_error); 

                }

                    printf("%10d %5d %5d  %7.2f (%7.2f)    %7.2f (%7.2f)   %8.2e  \n",
                       (int) batchCount, (int) M, (int) N, 
                       magma_tally3_perf,  1000.*magma_tally3_time,
                       cpu_perf,    1000.*cpu_time,
                       magma_tally3_error);
            }
            else {

                    printf("%10d %5d %5d  %7.2f (%7.2f)    ---   (  ---  )    ---\n",
                       (int) batchCount, (int) M, (int) N, 
                       magma_tally3_perf,  1000.*magma_tally3_time);
            }
            
            TESTING_FREE_CPU( h_A  );
            TESTING_FREE_CPU( h_X  );
            TESTING_FREE_CPU( h_Y  );
            TESTING_FREE_CPU( h_Ymagma_tally3  );


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
