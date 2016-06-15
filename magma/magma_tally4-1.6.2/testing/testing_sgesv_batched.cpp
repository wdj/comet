/*
   -- MAGMA_tally4 (version 1.6.1) --
   Univ. of Tennessee, Knoxville
   Univ. of California, Berkeley
   Univ. of Colorado, Denver
   @date January 2015

   @author Mark gates
   @author Azzam Haidar
   @author Tingxing Dong

   @generated from testing_zgesv_batched.cpp normal z -> s, Fri Jan 30 19:00:26 2015
 */
// includes, system
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sgesv_batched
*/
int main(int argc, char **argv)
{
    TESTING_INIT();

    real_Double_t   gflops, cpu_perf, cpu_time, gpu_perf, gpu_time;
    float          err = 0.0, Rnorm, Anorm, Xnorm, *work;
    float c_one     = MAGMA_tally4_S_ONE;
    float c_neg_one = MAGMA_tally4_S_NEG_ONE;
    float *h_A, *h_B, *h_X;
    magma_tally4Float_ptr d_A, d_B;
    magma_tally4_int_t *ipiv, *dipiv, *dinfo_array;
    magma_tally4_int_t N, nrhs, lda, ldb, ldda, lddb, info, sizeA, sizeB;
    magma_tally4_int_t ione     = 1;
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    magma_tally4_int_t status = 0;
    magma_tally4_int_t batchCount = 1;

    float **dA_array = NULL;
    float **dB_array = NULL;
    magma_tally4_int_t     **dipiv_array = NULL;

    magma_tally4_queue_t queue = magma_tally4_stream;
    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");
    
    nrhs = opts.nrhs;
    batchCount = opts.batchcount ;

    printf("BatchCount    N  NRHS   CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||B - AX|| / N*||A||*||X||\n");
    printf("================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            ldb    = lda;
            ldda   = ((N+31)/32)*32;
            lddb   = ldda;
            gflops = ( FLOPS_SGETRF( N, N ) + FLOPS_SGETRS( N, nrhs ) ) / 1e9 * batchCount;
            
            sizeA = lda*N*batchCount;
            sizeB = ldb*nrhs*batchCount;

            TESTING_MALLOC_CPU( h_A, float, sizeA );
            TESTING_MALLOC_CPU( h_B, float, sizeB );
            TESTING_MALLOC_CPU( h_X, float, sizeB );
            TESTING_MALLOC_CPU( work, float,      N);
            TESTING_MALLOC_CPU( ipiv, magma_tally4_int_t, N);
            
            TESTING_MALLOC_DEV( d_A, float, ldda*N*batchCount    );
            TESTING_MALLOC_DEV( d_B, float, lddb*nrhs*batchCount );
            TESTING_MALLOC_DEV( dipiv, magma_tally4_int_t, N * batchCount );
            TESTING_MALLOC_DEV( dinfo_array, magma_tally4_int_t, batchCount );

            magma_tally4_malloc((void**)&dA_array, batchCount * sizeof(*dA_array));
            magma_tally4_malloc((void**)&dB_array, batchCount * sizeof(*dB_array));
            magma_tally4_malloc((void**)&dipiv_array, batchCount * sizeof(*dipiv_array));


            /* Initialize the matrices */
            lapackf77_slarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_slarnv( &ione, ISEED, &sizeB, h_B );
            
            magma_tally4_ssetmatrix( N, N*batchCount,    h_A, lda, d_A, ldda );
            magma_tally4_ssetmatrix( N, nrhs*batchCount, h_B, ldb, d_B, lddb );
            
            /* ====================================================================
               Performs operation using MAGMA_tally4
               =================================================================== */
            sset_pointer(dA_array, d_A, ldda, 0, 0, ldda*N, batchCount, queue);
            sset_pointer(dB_array, d_B, lddb, 0, 0, lddb*nrhs, batchCount, queue);
            set_ipointer(dipiv_array, dipiv, 1, 0, 0, N, batchCount, queue);

            gpu_time = magma_tally4_wtime();
            //magma_tally4_sgesv_gpu( N, nrhs, d_A, ldda, ipiv, d_B, lddb, &info );
            info = magma_tally4_sgesv_batched(N, nrhs, dA_array, ldda, dipiv_array, dB_array, lddb, dinfo_array, batchCount, queue); 
            gpu_time = magma_tally4_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            // check correctness of results throught "dinfo_magma_tally4" and correctness of argument throught "info"
            magma_tally4_int_t *cpu_info = (magma_tally4_int_t*) malloc(batchCount*sizeof(magma_tally4_int_t));
            magma_tally4_getvector( batchCount, sizeof(magma_tally4_int_t), dinfo_array, 1, cpu_info, 1);
            for(int i=0; i<batchCount; i++)
            {
                if(cpu_info[i] != 0 ){
                    printf("magma_tally4_sgesv_batched matrix %d returned internal error %d\n",i, (int)cpu_info[i] );
                }
            }
            if (info != 0)
                printf("magma_tally4_sgesv_batched returned argument error %d: %s.\n", (int) info, magma_tally4_strerror( info ));
            
            //=====================================================================
            // Residual
            //=====================================================================
            magma_tally4_sgetmatrix( N, nrhs*batchCount, d_B, lddb, h_X, ldb );

            for(magma_tally4_int_t s=0; s<batchCount; s++)
            {
                Anorm = lapackf77_slange("I", &N, &N,    h_A + s * lda * N, &lda, work);
                Xnorm = lapackf77_slange("I", &N, &nrhs, h_X + s * ldb * nrhs, &ldb, work);
            
                blasf77_sgemm( Magma_tally4NoTransStr, Magma_tally4NoTransStr, &N, &nrhs, &N,
                           &c_one,     h_A + s * lda * N, &lda,
                                       h_X + s * ldb * nrhs, &ldb,
                           &c_neg_one, h_B + s * ldb * nrhs, &ldb);
            
                Rnorm = lapackf77_slange("I", &N, &nrhs, h_B + s * ldb * nrhs, &ldb, work);
                float error = Rnorm/(N*Anorm*Xnorm);
                
                if ( isnan(error) || isinf(error) ) {
                    err = error;
                    break;
                }
                err = max(err, error);            
            }
            status += ! (err < tol);

            /* ====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally4_wtime();
                for(magma_tally4_int_t s=0; s<batchCount; s++)
                {
                    lapackf77_sgesv( &N, &nrhs, h_A + s * lda * N, &lda, ipiv, h_B + s * ldb * nrhs, &ldb, &info );
                }
                cpu_time = magma_tally4_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_sgesv returned err %d: %s.\n",
                           (int) info, magma_tally4_strerror( info ));
                
                printf( "%10d    %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                        (int)batchCount, (int) N, (int) nrhs, cpu_perf, cpu_time, gpu_perf, gpu_time,
                        err, (err < tol ? "ok" : "failed"));
            }
            else {
                printf( "%10d    %5d %5d     ---   (  ---  )   %7.2f (%7.2f)   %8.2e   %s\n",
                        (int)batchCount, (int) N, (int) nrhs, gpu_perf, gpu_time,
                        err, (err < tol ? "ok" : "failed"));
            }
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_B );
            TESTING_FREE_CPU( h_X );
            TESTING_FREE_CPU( work );
            TESTING_FREE_CPU( ipiv );
            
            TESTING_FREE_DEV( d_A );
            TESTING_FREE_DEV( d_B );

            TESTING_FREE_DEV( dipiv );
            TESTING_FREE_DEV( dinfo_array );

            magma_tally4_free(dA_array);
            magma_tally4_free(dB_array);
            magma_tally4_free(dipiv_array);
            free(cpu_info);
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
