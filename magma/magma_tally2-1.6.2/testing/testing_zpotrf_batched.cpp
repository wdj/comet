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

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zpotrf_batched
*/

int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    magma_tally2DoubleComplex *h_A, *h_R;
    magma_tally2DoubleComplex *d_A;
    magma_tally2_int_t N, n2, lda, ldda, info;
    magma_tally2DoubleComplex c_neg_one = MAGMA_tally2_Z_NEG_ONE;
    magma_tally2_int_t ione     = 1;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    double      work[1], error;
    magma_tally2_int_t status = 0;
    magma_tally2DoubleComplex **d_A_array = NULL;
    magma_tally2_int_t *dinfo_magma_tally2;

    magma_tally2_int_t batchCount;

    magma_tally2_queue_t queue = magma_tally2_stream;
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    batchCount = opts.batchcount;
    double tol = opts.tolerance * lapackf77_dlamch("E");

    printf("BatchCount    N      CPU GFlop/s (ms)      GPU GFlop/s (ms)    ||R_magma_tally2 - R_lapack||_F / ||R_lapack||_F\n");
    printf("========================================================\n");
    for( int i = 0; i < opts.ntest; ++i ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N   = opts.nsize[i];
            ldda = lda = ((N+31)/32)*32;
            n2  = lda* N  * batchCount;

            gflops = batchCount * FLOPS_ZPOTRF( N ) / 1e9 ;

            TESTING_MALLOC_CPU( h_A, magma_tally2DoubleComplex, n2);
            TESTING_MALLOC_PIN( h_R, magma_tally2DoubleComplex, n2);
            TESTING_MALLOC_DEV(  d_A, magma_tally2DoubleComplex, ldda * N * batchCount);
            TESTING_MALLOC_DEV(  dinfo_magma_tally2,  magma_tally2_int_t, batchCount);
            
            magma_tally2_malloc((void**)&d_A_array, batchCount * sizeof(*d_A_array));

            /* Initialize the matrix */
            lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
            for(int i=0; i<batchCount; i++)
            {
               magma_tally2_zmake_hpd( N, h_A + i * lda * N, lda );// need modification
            }
            
            magma_tally2_int_t columns = N * batchCount;
            lapackf77_zlacpy( Magma_tally2UpperLowerStr, &N, &(columns), h_A, &lda, h_R, &lda );
            magma_tally2_zsetmatrix( N, columns, h_A, lda, d_A, ldda );


            /* ====================================================================
               Performs operation using MAGMA_tally2
               =================================================================== */
            zset_pointer(d_A_array, d_A, ldda, 0, 0, ldda * N, batchCount, queue);
            gpu_time = magma_tally2_sync_wtime(NULL);
            info = magma_tally2_zpotrf_batched( opts.uplo, N, d_A_array, ldda, dinfo_magma_tally2, batchCount, queue);
            gpu_time = magma_tally2_sync_wtime(NULL) - gpu_time;
            gpu_perf = gflops / gpu_time;
            magma_tally2_int_t *cpu_info = (magma_tally2_int_t*) malloc(batchCount*sizeof(magma_tally2_int_t));
            magma_tally2_getvector( batchCount, sizeof(magma_tally2_int_t), dinfo_magma_tally2, 1, cpu_info, 1);
            for(int i=0; i<batchCount; i++)
            {
                if(cpu_info[i] != 0 ){
                    printf("magma_tally2_zpotrf_batched matrix %d returned internal error %d\n",i, (int)cpu_info[i] );
                }
            }
            if (info != 0)
                printf("magma_tally2_zpotrf_batched returned argument error %d: %s.\n", (int) info, magma_tally2_strerror( info ));

            if ( opts.lapack ) {

                /* =====================================================================
                   Performs operation using LAPACK
                   =================================================================== */
                cpu_time = magma_tally2_wtime();
                for(int i=0; i<batchCount; i++)
                {
                   lapackf77_zpotrf( lapack_uplo_const_tally2(opts.uplo), &N, h_A + i * lda * N, &lda, &info );
                }
                cpu_time = magma_tally2_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_zpotrf returned error %d: %s.\n",
                           (int) info, magma_tally2_strerror( info ));

                /* =====================================================================
                   Check the result compared to LAPACK
                   =================================================================== */
                 magma_tally2_zgetmatrix( N, columns, d_A, ldda, h_R, lda );
                 magma_tally2_int_t NN = lda*N;
                 char const uplo = 'l'; // lapack_uplo_const_tally2(opts.uplo)
                 double err = 0.0;
                 for(int i=0; i<batchCount; i++)
                 { 
                     error = lapackf77_zlanhe("f", &uplo, &N, h_A + i * lda*N, &lda, work);                
                     blasf77_zaxpy(&NN, &c_neg_one, h_A + i * lda*N, &ione, h_R + i  * lda*N, &ione);
                     error = lapackf77_zlanhe("f", &uplo, &N, h_R + i * lda*N, &lda, work) / error;
                     if ( isnan(error) || isinf(error) ) {
                         err = error;
                         break;
                     }
                     err = max(fabs(error),err);
                 }
              

                printf("%5d      %5d    %7.2f (%7.2f)     %7.2f (%7.2f)     %8.2e   %s\n",
                       (int)batchCount, (int) N, cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000., err,  (error < tol ? "ok" : "failed"));
                status += ! (err < tol);
                
            }
            else {
                printf("%5d      %5d    ---   (  ---  )   %7.2f (%7.2f)     ---  \n",
                       (int)batchCount, (int) N, gpu_perf, gpu_time*1000. );
            }
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_PIN( h_R );
            TESTING_FREE_DEV( d_A );
            TESTING_FREE_DEV( d_A_array );
            TESTING_FREE_DEV( dinfo_magma_tally2 );
            free(cpu_info);
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;

}



