/*
 *  -- MAGMA_minproduct (version 1.6.1) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     @date January 2015
 *
 * @precisions normal z -> c d s
 *
 **/
/* includes, system */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* includes, project */
#include "flops.h"
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "magma_minproduct_types.h"
#include "testings.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zpotrf_mgpu
*/
int main( int argc, char** argv) 
{
    /* Initialize */
    TESTING_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double      work[1], error;
    magma_minproductDoubleComplex *h_A, *h_R;
    magma_minproductDoubleComplex_ptr d_lA[4] = {NULL, NULL, NULL, NULL};
    magma_minproduct_int_t N, n2, lda, ldda, info;
    magma_minproduct_int_t j, k, ngpu0 = 1, ngpu;
    magma_minproductDoubleComplex c_neg_one = MAGMA_minproduct_Z_NEG_ONE;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t nb, nk, n_local, ldn_local;

    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    ngpu0 = opts.ngpu;

    printf("ngpu = %d, uplo = %s\n", (int) opts.ngpu, lapack_uplo_const(opts.uplo) );
    printf("  N     CPU GFlop/s (sec)   MAGMA_minproduct GFlop/s (sec)   ||R_magma_minproduct - R_lapack||_F / ||R_lapack||_F\n");
    printf("=============================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N   = opts.nsize[itest];
            lda = N;
            n2  = lda*N;
            ldda = ((N+31)/32)*32;
            gflops = FLOPS_ZPOTRF( N ) / 1e9;

            magma_minproduct_setdevice(0);
            TESTING_MALLOC(    h_A, magma_minproductDoubleComplex, n2 );
            TESTING_HOSTALLOC( h_R, magma_minproductDoubleComplex, n2 );

            //TESTING_MALLOC_DEV(  d_A, magma_minproductDoubleComplex, ldda*N );
            nb = magma_minproduct_get_zpotrf_nb(N);
            if( ngpu0 > N / nb ) {
                ngpu = N / nb;
                if( N % nb != 0 ) ngpu++;
                printf( " * too many gpus for the matrix size, using %d gpus\n", (int) ngpu );
            } else {
                ngpu = ngpu0;
            }

            for(j = 0; j < ngpu; j++){
                n_local = nb*(N /(nb * ngpu));
                if (j < (N / nb) % ngpu)
                    n_local += nb;
                else if (j == (N / nb) % ngpu)
                    n_local += N % nb;

                ldn_local = ((n_local + 31) / 32) * 32;
                ldn_local = (ldn_local % 256 == 0) ? ldn_local + 32 : ldn_local;

                magma_minproduct_setdevice(j);
                TESTING_DEVALLOC( d_lA[j], magma_minproductDoubleComplex, ldda * ldn_local );
            }

            /* Initialize the matrix */
            if(opts.check) {
                lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
                magma_minproduct_zmake_hpd( N, h_A, lda );
                lapackf77_zlacpy( Magma_minproductUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
            } else {
                lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
                magma_minproduct_zmake_hpd( N, h_A, lda );
            }

            /* ====================================================================
               Performs operation using MAGMA_minproduct 
               =================================================================== */
            /* distribute matrix to gpus */
            //magma_minproduct_zprint( N,N, h_A, lda );
            //if( opts.uplo == Magma_minproductUpper) {
                for(j = 0; j < N; j += nb){
                    k = (j / nb) % ngpu;
                    magma_minproduct_setdevice(k);
                    nk = min(nb, N - j);
                    magma_minproduct_zsetmatrix( N, nk,
                            h_A + j * lda,                       lda,
                            d_lA[k] + j / (nb * ngpu) * nb * ldda, ldda);
                }
            /*} else {
            }*/

            gpu_time = magma_minproduct_wtime();
            //magma_minproduct_zpotrf_mgpu(ngpu, opts.uplo, N, d_lA, ldda, &info);
            magma_minproduct_zpotrf_mgpu_right(ngpu, opts.uplo, N, d_lA, ldda, &info);
            gpu_time = magma_minproduct_wtime() - gpu_time;
            if (info != 0)
                printf("magma_minproduct_zpotrf_gpu returned error %d.\n", (int) info);
            gpu_perf = gflops / gpu_time;

            if ( opts.check && info == 0) {
                /* =====================================================================
                   Performs operation using LAPACK
                   =================================================================== */
                cpu_time = magma_minproduct_wtime();
                lapackf77_zpotrf( lapack_uplo_const(opts.uplo), &N, h_A, &lda, &info );
                //printf( " ==== LAPACK ====\n" );
                //magma_minproduct_zprint( N,N, h_A, lda );
                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_zpotrf returned error %d.\n", (int) info);

                /* =====================================================================
                   Check the result compared to LAPACK
                   =================================================================== */
                /* gather matrix from gpus */
                //if( opts.uplo == Magma_minproductUpper ) {
                    for(j = 0; j < N; j += nb){
                        k = (j / nb) % ngpu;
                        magma_minproduct_setdevice(k);
                        nk = min(nb, N - j);
                        magma_minproduct_zgetmatrix( N, nk,
                                d_lA[k] + j / (nb * ngpu) * nb * ldda, ldda,
                                h_R + j * lda,                        lda );
                    }
                /*} else {
                }*/
                magma_minproduct_setdevice(0);
                //printf( " ==== MAGMA_minproduct ====\n" );
                //magma_minproduct_zprint( N,N, h_R, lda );


                //error = lapackf77_zlange("f", &N, &N, h_A, &lda, work);
                error = lapackf77_zlanhe("f", "L", &N, h_A, &lda, work);
                blasf77_zaxpy(&n2, &c_neg_one, h_A, &ione, h_R, &ione);
                //error = lapackf77_zlange("f", &N, &N, h_R, &lda, work) / error;
                error = lapackf77_zlanhe("f", "L", &N, h_R, &lda, work) / error;

                printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e\n",
                        (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time, error );
            }
            else {
                printf("%5d     ---   (  ---  )   %7.2f (%7.2f)     ---  \n",
                        (int) N, gpu_perf, gpu_time );
            }

            for(j = 0; j < ngpu; j++) {
                magma_minproduct_setdevice(j);
                TESTING_DEVFREE( d_lA[j] );
            }
            magma_minproduct_setdevice(0);
            TESTING_FREE( h_A );
            TESTING_HOSTFREE( h_R );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    /* Shutdown */
    TESTING_FINALIZE();

    return 0;
}

