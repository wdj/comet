/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> c d s
       @author Mark Gates
*/
// includes, system
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_tally2.h"
#include "magma_tally2_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zgesv_gpu
*/
int main(int argc, char **argv)
{
    TESTING_INIT();

    real_Double_t   gflops, cpu_perf, cpu_time, gpu_perf, gpu_time;
    double          error, Rnorm, Anorm, Xnorm, *work;
    magma_tally2DoubleComplex c_one     = MAGMA_tally2_Z_ONE;
    magma_tally2DoubleComplex c_neg_one = MAGMA_tally2_Z_NEG_ONE;
    magma_tally2DoubleComplex *h_A, *h_B, *h_X;
    magma_tally2DoubleComplex *d_A, *d_B;
    magma_tally2DoubleComplex **d_A_array, **d_B_array;
    magma_tally2_int_t *ipiv;
    magma_tally2_int_t N, n2, nrhs, lda, ldb, ldda, lddb, info, sizeA, sizeB;
    magma_tally2_int_t ione     = 1;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    magma_tally2_int_t status = 0;
    magma_tally2_int_t     *dinfo_magma_tally2;
    
    magma_tally2_int_t batchCount = 1;
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );
    
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    magma_tally2_queue_t queue = magma_tally2_stream;
    batchCount = opts.batchcount ;
    magma_tally2_int_t columns;
    nrhs = opts.nrhs;
    
    printf("  Batchcount  N  NRHS   CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||B - AX|| / N*||A||*||X||\n");
    printf("================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            ldb    = lda;
            ldda   = ((N+31)/32)*32;
            n2     = lda*N * batchCount;
            lddb   = ldda;
            gflops = ( FLOPS_ZGETRF( N, N ) + FLOPS_ZGETRS( N, nrhs ) ) / 1e9 * batchCount;
            
            TESTING_MALLOC_CPU( h_A, magma_tally2DoubleComplex, n2    );
            TESTING_MALLOC_CPU( h_B, magma_tally2DoubleComplex, ldb*nrhs*batchCount );
            TESTING_MALLOC_CPU( h_X, magma_tally2DoubleComplex, ldb*nrhs*batchCount );
            TESTING_MALLOC_CPU( work, double,      N );
            TESTING_MALLOC_CPU( ipiv, magma_tally2_int_t, N );
            
            TESTING_MALLOC_DEV(  dinfo_magma_tally2,  magma_tally2_int_t, batchCount);
            
            TESTING_MALLOC_DEV( d_A, magma_tally2DoubleComplex, ldda*N*batchCount    );
            TESTING_MALLOC_DEV( d_B, magma_tally2DoubleComplex, lddb*nrhs*batchCount );
            

            magma_tally2_malloc((void**)&d_A_array, batchCount * sizeof(*d_A_array));
            magma_tally2_malloc((void**)&d_B_array, batchCount * sizeof(*d_B_array));

            /* Initialize the matrices */
            sizeA = n2;
            sizeB = ldb*nrhs*batchCount;
            lapackf77_zlarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_zlarnv( &ione, ISEED, &sizeB, h_B );
            
            columns = N * batchCount;
            magma_tally2_zsetmatrix( N, columns,    h_A, lda, d_A, ldda );
            magma_tally2_zsetmatrix( N, nrhs*batchCount, h_B, ldb, d_B, lddb );
            
            /* ====================================================================
               Performs operation using MAGMA_tally2
               =================================================================== */
            
            
            
            zset_pointer(d_A_array, d_A, ldda, 0, 0, ldda*N, batchCount, queue);
            zset_pointer(d_B_array, d_B, lddb, 0, 0, lddb*nrhs, batchCount, queue);
            
            
            
            gpu_time = magma_tally2_wtime();
            magma_tally2_zgesv_nopiv_batched( N, nrhs, d_A_array, ldda, d_B_array, lddb, dinfo_magma_tally2, batchCount, queue );
            gpu_time = magma_tally2_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            // check correctness of results throught "dinfo_magma_tally2" and correctness of argument throught "info"
            magma_tally2_int_t *cpu_info = (magma_tally2_int_t*) malloc(batchCount*sizeof(magma_tally2_int_t));
            magma_tally2_getvector( batchCount, sizeof(magma_tally2_int_t), dinfo_magma_tally2, 1, cpu_info, 1);
            for(int i=0; i<batchCount; i++)
            {
                if(cpu_info[i] != 0 ){
                    printf("magma_tally2_zgetrf_batched matrix %d returned internal error %d\n",i, (int)cpu_info[i] );
                }
            }
            if (info != 0)
                printf("magma_tally2_zgetrf_batched returned argument error %d: %s.\n", (int) info, magma_tally2_strerror( info ));
            
            //=====================================================================
            // Residual
            //=====================================================================
            double err = 0;
            error = 0;
            magma_tally2_zgetmatrix( N, nrhs*batchCount, d_B, lddb, h_X, ldb );
            for(magma_tally2_int_t s = 0; s < batchCount; s++)
            {
                Anorm = lapackf77_zlange("I", &N, &N,    h_A+s*lda*N, &lda, work);
                Xnorm = lapackf77_zlange("I", &N, &nrhs, h_X+s*ldb*nrhs, &ldb, work);

                blasf77_zgemm( Magma_tally2NoTransStr, Magma_tally2NoTransStr, &N, &nrhs, &N,
                        &c_one,     h_A+s*lda*N, &lda,
                        h_X+s*ldb*nrhs, &ldb,
                        &c_neg_one, h_B+s*ldb*nrhs, &ldb);

                Rnorm = lapackf77_zlange("I", &N, &nrhs, h_B+s*ldb*nrhs, &ldb, work);
                err += Rnorm/(N*Anorm*Xnorm);
                if ( isnan(error) || isinf(error) ) {
                    err = error;
                    break;
                }
            }
            //printf("before error\n");
            error = err/(double)batchCount;
            //error = max(error, err);
            status += ! (error < tol);
            
            /* ====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally2_wtime();
                lapackf77_zgesv( &N, &nrhs, h_A, &lda, ipiv, h_B, &ldb, &info );
                cpu_time = magma_tally2_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_zgesv returned error %d: %s.\n",
                           (int) info, magma_tally2_strerror( info ));
                
                printf( "%5d %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                        (int)batchCount, (int) N, (int) nrhs, cpu_perf, cpu_time, gpu_perf, gpu_time,
                        error, (error < tol ? "ok" : "failed"));
            }
            else {
                printf( "%5d %5d %5d     ---   (  ---  )   %7.2f (%7.2f)   %8.2e   %s\n",
                        (int) batchCount, (int) N, (int) nrhs, gpu_perf, gpu_time,
                        error, (error < tol ? "ok" : "failed"));
            }
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_B );
            TESTING_FREE_CPU( h_X );
            TESTING_FREE_CPU( work );
            TESTING_FREE_CPU( ipiv );
            
            TESTING_FREE_DEV( d_A );
            TESTING_FREE_DEV( d_B );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
