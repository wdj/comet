/*
   -- MAGMA_minproduct (version 1.6.1) --
   Univ. of Tennessee, Knoxville
   Univ. of California, Berkeley
   Univ. of Colorado, Denver
   @date January 2015

   @author Mark gates
   @author Azzam Haidar
   @author Tingxing Dong

   @generated from testing_zgesv_batched.cpp normal z -> d, Fri Jan 30 19:00:26 2015
 */
// includes, system
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dgesv_batched
*/
int main(int argc, char **argv)
{
    TESTING_INIT();

    real_Double_t   gflops, cpu_perf, cpu_time, gpu_perf, gpu_time;
    double          err = 0.0, Rnorm, Anorm, Xnorm, *work;
    double c_one     = MAGMA_minproduct_D_ONE;
    double c_neg_one = MAGMA_minproduct_D_NEG_ONE;
    double *h_A, *h_B, *h_X;
    magma_minproductDouble_ptr d_A, d_B;
    magma_minproduct_int_t *ipiv, *dipiv, *dinfo_array;
    magma_minproduct_int_t N, nrhs, lda, ldb, ldda, lddb, info, sizeA, sizeB;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t status = 0;
    magma_minproduct_int_t batchCount = 1;

    double **dA_array = NULL;
    double **dB_array = NULL;
    magma_minproduct_int_t     **dipiv_array = NULL;

    magma_minproduct_queue_t queue = magma_minproduct_stream;
    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
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
            gflops = ( FLOPS_DGETRF( N, N ) + FLOPS_DGETRS( N, nrhs ) ) / 1e9 * batchCount;
            
            sizeA = lda*N*batchCount;
            sizeB = ldb*nrhs*batchCount;

            TESTING_MALLOC_CPU( h_A, double, sizeA );
            TESTING_MALLOC_CPU( h_B, double, sizeB );
            TESTING_MALLOC_CPU( h_X, double, sizeB );
            TESTING_MALLOC_CPU( work, double,      N);
            TESTING_MALLOC_CPU( ipiv, magma_minproduct_int_t, N);
            
            TESTING_MALLOC_DEV( d_A, double, ldda*N*batchCount    );
            TESTING_MALLOC_DEV( d_B, double, lddb*nrhs*batchCount );
            TESTING_MALLOC_DEV( dipiv, magma_minproduct_int_t, N * batchCount );
            TESTING_MALLOC_DEV( dinfo_array, magma_minproduct_int_t, batchCount );

            magma_minproduct_malloc((void**)&dA_array, batchCount * sizeof(*dA_array));
            magma_minproduct_malloc((void**)&dB_array, batchCount * sizeof(*dB_array));
            magma_minproduct_malloc((void**)&dipiv_array, batchCount * sizeof(*dipiv_array));


            /* Initialize the matrices */
            lapackf77_dlarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_dlarnv( &ione, ISEED, &sizeB, h_B );
            
            magma_minproduct_dsetmatrix( N, N*batchCount,    h_A, lda, d_A, ldda );
            magma_minproduct_dsetmatrix( N, nrhs*batchCount, h_B, ldb, d_B, lddb );
            
            /* ====================================================================
               Performs operation using MAGMA_minproduct
               =================================================================== */
            dset_pointer(dA_array, d_A, ldda, 0, 0, ldda*N, batchCount, queue);
            dset_pointer(dB_array, d_B, lddb, 0, 0, lddb*nrhs, batchCount, queue);
            set_ipointer(dipiv_array, dipiv, 1, 0, 0, N, batchCount, queue);

            gpu_time = magma_minproduct_wtime();
            //magma_minproduct_dgesv_gpu( N, nrhs, d_A, ldda, ipiv, d_B, lddb, &info );
            info = magma_minproduct_dgesv_batched(N, nrhs, dA_array, ldda, dipiv_array, dB_array, lddb, dinfo_array, batchCount, queue); 
            gpu_time = magma_minproduct_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            // check correctness of results throught "dinfo_magma_minproduct" and correctness of argument throught "info"
            magma_minproduct_int_t *cpu_info = (magma_minproduct_int_t*) malloc(batchCount*sizeof(magma_minproduct_int_t));
            magma_minproduct_getvector( batchCount, sizeof(magma_minproduct_int_t), dinfo_array, 1, cpu_info, 1);
            for(int i=0; i<batchCount; i++)
            {
                if(cpu_info[i] != 0 ){
                    printf("magma_minproduct_dgesv_batched matrix %d returned internal error %d\n",i, (int)cpu_info[i] );
                }
            }
            if (info != 0)
                printf("magma_minproduct_dgesv_batched returned argument error %d: %s.\n", (int) info, magma_minproduct_strerror( info ));
            
            //=====================================================================
            // Residual
            //=====================================================================
            magma_minproduct_dgetmatrix( N, nrhs*batchCount, d_B, lddb, h_X, ldb );

            for(magma_minproduct_int_t s=0; s<batchCount; s++)
            {
                Anorm = lapackf77_dlange("I", &N, &N,    h_A + s * lda * N, &lda, work);
                Xnorm = lapackf77_dlange("I", &N, &nrhs, h_X + s * ldb * nrhs, &ldb, work);
            
                blasf77_dgemm( Magma_minproductNoTransStr, Magma_minproductNoTransStr, &N, &nrhs, &N,
                           &c_one,     h_A + s * lda * N, &lda,
                                       h_X + s * ldb * nrhs, &ldb,
                           &c_neg_one, h_B + s * ldb * nrhs, &ldb);
            
                Rnorm = lapackf77_dlange("I", &N, &nrhs, h_B + s * ldb * nrhs, &ldb, work);
                double error = Rnorm/(N*Anorm*Xnorm);
                
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
                cpu_time = magma_minproduct_wtime();
                for(magma_minproduct_int_t s=0; s<batchCount; s++)
                {
                    lapackf77_dgesv( &N, &nrhs, h_A + s * lda * N, &lda, ipiv, h_B + s * ldb * nrhs, &ldb, &info );
                }
                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_dgesv returned err %d: %s.\n",
                           (int) info, magma_minproduct_strerror( info ));
                
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

            magma_minproduct_free(dA_array);
            magma_minproduct_free(dB_array);
            magma_minproduct_free(dipiv_array);
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
