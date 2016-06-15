/*
    -- MAGMA_tally3 (version 1.6) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @author Tingxing Dong
       @author Azzam Haidar

       @generated from testing_zgeqrf_batched.cpp normal z -> s, Fri May  1 21:31:40 2015

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <cuda.h>  // for CUDA_VERSION

// includes, project
#include "flops.h"
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "testings.h"
#include "batched_kernel_param.h"


void get_QR_error(magma_tally3_int_t M, magma_tally3_int_t N, magma_tally3_int_t min_mn,
                    float *h_R,  float *h_A, magma_tally3_int_t lda,
                    float *tau,
                    float *Q,  magma_tally3_int_t ldq,
                    float *R,  magma_tally3_int_t ldr,
                    float *h_work,  magma_tally3_int_t lwork,
                    float *work, float *error, float *error2)
{
    /* h_R:input the factorized matrix by lapack QR,
       h_A:input the original matrix copy
       tau: input
    */
    
    const float             d_neg_one = MAGMA_tally3_D_NEG_ONE;
    const float             d_one     = MAGMA_tally3_D_ONE;
    const float c_neg_one = MAGMA_tally3_S_NEG_ONE;
    const float c_one     = MAGMA_tally3_S_ONE;
    const float c_zero    = MAGMA_tally3_S_ZERO;
    float           Anorm;
    
    magma_tally3_int_t info;
    
    // generate M by K matrix Q, where K = min(M,N)
    lapackf77_slacpy( "Lower", &M, &min_mn, h_R, &lda, Q, &ldq );
    lapackf77_sorgqr( &M, &min_mn, &min_mn, Q, &ldq, tau, h_work, &lwork, &info );
    assert( info == 0 );
    
    // copy K by N matrix R
    lapackf77_slaset( "Lower", &min_mn, &N, &c_zero, &c_zero, R, &ldr );
    lapackf77_slacpy( "Upper", &min_mn, &N, h_R, &lda,        R, &ldr );
    
    // error = || R - Q^H*A || / (N * ||A||)
    blasf77_sgemm( "Conj", "NoTrans", &min_mn, &N, &M,
    &c_neg_one, Q, &ldq, h_A, &lda, &c_one, R, &ldr );
    
    Anorm = lapackf77_slange( "1", &M,      &N, h_A, &lda, work );
    *error = lapackf77_slange( "1", &min_mn, &N, R,   &ldr, work );
    
    if ( N > 0 && Anorm > 0 )
        *error /= (N*Anorm);
    
    // set R = I (K by K identity), then R = I - Q^H*Q
    // error = || I - Q^H*Q || / N
    lapackf77_slaset( "Upper", &min_mn, &min_mn, &c_zero, &c_one, R, &ldr );
    blasf77_ssyrk( "Upper", "Conj", &min_mn, &M, &d_neg_one, Q, &ldq, &d_one, R, &ldr );
    *error2 = lapackf77_slansy( "1", "Upper", &min_mn, R, &ldr, work );
    if ( N > 0 )
        *error2 /= N;
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing sgeqrf_batched
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gflops, magma_tally3_perf, magma_tally3_time, cublas_perf=0, cublas_time=0, cpu_perf, cpu_time;
    float           magma_tally3_error=0.0, cublas_error=0.0, magma_tally3_error2=0.0, cublas_error2=0.0, error, error2;

    float *h_A, *h_R, *h_Amagma_tally3, *tau, *h_work, tmp[1];
    float *d_A, *dtau_magma_tally3, *dtau_cublas;

    float **dA_array = NULL;
    float **dtau_array = NULL;

    magma_tally3_int_t   *dinfo_magma_tally3, *dinfo_cublas;

    magma_tally3_int_t M, N, lda, ldda, lwork, n2, info, min_mn;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t status = 0;

    magma_tally3_int_t batchCount = 1;
    magma_tally3_int_t column;

    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    batchCount = opts.batchcount;

    float tol = opts.tolerance * lapackf77_slamch("E");
    
    magma_tally3_queue_t stream[2];
    magma_tally3_queue_create( &stream[0] );
    magma_tally3_queue_create( &stream[1] );
    printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |R - Q^H*A|   |I - Q^H*Q|\n");
    printf("BatchCount     M    N     MAGMA_tally3 GFlop/s (ms)   CUBLAS GFlop/s (ms)   CPU GFlop/s (ms)   |R - Q^H*A|_mag   |I - Q^H*Q|_mag   |R - Q^H*A|_cub   |I - Q^H*Q|_cub \n");
    printf("============================================================================================================================================================= \n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M     = opts.msize[itest];
            N     = opts.nsize[itest];
            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N * batchCount;
            ldda = M;
            ldda   = magma_tally3_roundup( M, opts.roundup );  // multiple of 32 by default

            gflops = (FLOPS_SGEQRF( M, N ) + FLOPS_SGEQRT( M, N )) / 1e9 * batchCount;

            /* Allocate memory for the matrix */
            TESTING_MALLOC_CPU( tau,   float, min_mn * batchCount );
            TESTING_MALLOC_CPU( h_A,   float, n2     );
            TESTING_MALLOC_CPU( h_Amagma_tally3,   float, n2     );
            TESTING_MALLOC_PIN( h_R,   float, n2     );
        
            TESTING_MALLOC_DEV( d_A,   float, ldda*N * batchCount );

            TESTING_MALLOC_DEV( dtau_magma_tally3,  float, min_mn * batchCount);
            TESTING_MALLOC_DEV( dtau_cublas,  float, min_mn * batchCount);

            TESTING_MALLOC_DEV(  dinfo_magma_tally3,  magma_tally3_int_t, batchCount);
            TESTING_MALLOC_DEV(  dinfo_cublas,  magma_tally3_int_t, batchCount);

            magma_tally3_malloc((void**)&dA_array, batchCount * sizeof(*dA_array));
            magma_tally3_malloc((void**)&dtau_array, batchCount * sizeof(*dtau_array));
        
            // to determine the size of lwork
            lwork = -1;
            lapackf77_sgeqrf(&M, &N, NULL, &M, NULL, tmp, &lwork, &info);
            lwork = (magma_tally3_int_t)MAGMA_tally3_S_REAL( tmp[0] );
            lwork = max(lwork, N*N);
           
            TESTING_MALLOC_CPU( h_work, float, lwork * batchCount);

            column = N * batchCount;
            /* Initialize the matrix */
            lapackf77_slarnv( &ione, ISEED, &n2, h_A );
            for(int i=0; i<batchCount; i++)
			{
				magma_tally3_smake_hpd( N, h_A + i * lda * N, lda );// need modification
			}
			
			lapackf77_slacpy( Magma_tally3UpperLowerStr, &M, &column, h_A, &lda, h_R, &lda );
       
            /* ====================================================================
               Performs operation using MAGMA_tally3
               =================================================================== */
            magma_tally3_ssetmatrix( M, column, h_R, lda,  d_A, ldda );
            sset_pointer(dA_array, d_A, 1, 0, 0, ldda*N, batchCount, opts.queue);
            sset_pointer(dtau_array, dtau_magma_tally3, 1, 0, 0, min_mn, batchCount, opts.queue);
    
            magma_tally3_time = magma_tally3_sync_wtime(0);
    
            info = magma_tally3_sgeqrf_batched(M, N, dA_array, ldda, dtau_array, dinfo_magma_tally3, batchCount, opts.queue);

            magma_tally3_time = magma_tally3_sync_wtime(0) - magma_tally3_time;
            magma_tally3_perf = gflops / magma_tally3_time;

            magma_tally3_sgetmatrix( M, column, d_A, ldda, h_Amagma_tally3, lda);

            if (info != 0)
                printf("magma_tally3_sgeqrf_batched returned error %d: %s.\n",
                       (int) info, magma_tally3_strerror( info ));
            
            /* ====================================================================
               Performs operation using CUBLAS
               =================================================================== */

            /* cublasSgeqrfBatched is only available from CUBLAS v6.5 */
            #if CUDA_VERSION >= 6050
            magma_tally3_ssetmatrix( M, column, h_R, lda,  d_A, ldda );
            sset_pointer(dA_array, d_A, 1, 0, 0, ldda*N, batchCount, opts.queue);
            sset_pointer(dtau_array, dtau_cublas, 1, 0, 0, min_mn, batchCount, opts.queue);

            cublas_time = magma_tally3_sync_wtime(0);
    
            int cublas_info;  // int, not magma_tally3_int_t
            cublasSgeqrfBatched(opts.handle, M, N, dA_array, ldda, dtau_array, &cublas_info, batchCount);

            cublas_time = magma_tally3_sync_wtime(0) - cublas_time;
            cublas_perf = gflops / cublas_time;

            if (cublas_info != 0)
                printf("cublasSgeqrfBatched returned error %d: %s.\n",
                       (int) cublas_info, magma_tally3_strerror( cublas_info ));
            #endif

            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.check ) {
                cpu_time = magma_tally3_wtime();

                for (int i=0; i < batchCount; i++)
                {
                    lapackf77_sgeqrf(&M, &N, h_A + i*lda*N, &lda, tau + i*min_mn, h_work + i * lwork, &lwork, &info);
                  
                    /* lapackf77_slarft( Magma_tally3ForwardStr, Magma_tally3ColumnwiseStr,
                                     &M, &N, h_A + i*lda*N, &lda, tau + i*min_mn, h_work + i * lwork, &N); */
                }

                cpu_time = magma_tally3_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_sgeqrf returned error %d: %s.\n",
                           (int) info, magma_tally3_strerror( info ));
                
                /* =====================================================================
                   Check the MAGMA_tally3 CUBLAS result compared to LAPACK
                   =================================================================== */
                magma_tally3_int_t ldq = M;
                magma_tally3_int_t ldr = min_mn;
                float *Q, *R;
                float *work;

                TESTING_MALLOC_CPU( Q,    float, ldq*min_mn );  // M by K
                TESTING_MALLOC_CPU( R,    float, ldr*N );       // K by N
                TESTING_MALLOC_CPU( work, float,             min_mn );

                /* check magma_tally3 result */
                magma_tally3_sgetvector(min_mn*batchCount, dtau_magma_tally3, 1, tau, 1);
                for (int i=0; i < batchCount; i++)
                {
                    get_QR_error(M, N, min_mn,
                             h_Amagma_tally3 + i*lda*N, h_R + i*lda*N, lda, tau + i*min_mn,
                             Q, ldq, R, ldr, h_work, lwork,
                             work, &error, &error2);

                    if ( isnan(error) || isinf(error) ) {
                        magma_tally3_error = error;
                        break;
                    }
                    magma_tally3_error  = max( fabs(error),  magma_tally3_error  );
                    magma_tally3_error2 = max( fabs(error2), magma_tally3_error2 );
                }

                /* check cublas result */
                #if CUDA_VERSION >= 6050
                magma_tally3_sgetvector(min_mn*batchCount, dtau_magma_tally3, 1, tau, 1);
                magma_tally3_sgetmatrix( M, column, d_A, ldda, h_A, lda);
                for (int i=0; i < batchCount; i++)
                {
                    get_QR_error(M, N, min_mn,
                             h_A + i*lda*N, h_R + i*lda*N, lda, tau + i*min_mn,
                             Q, ldq, R, ldr, h_work, lwork,
                             work, &error, &error2);

                    if ( isnan(error) || isinf(error) ) {
                        cublas_error = error;
                        break;
                    }
                    cublas_error  = max( fabs(error),  cublas_error  );
                    cublas_error2 = max( fabs(error2), cublas_error2 );
                }
                #endif

                TESTING_FREE_CPU( Q    );  Q    = NULL;
                TESTING_FREE_CPU( R    );  R    = NULL;
                TESTING_FREE_CPU( work );  work = NULL;

                bool okay = (magma_tally3_error < tol && magma_tally3_error2 < tol);
                //bool okay_cublas = (cublas_error < tol && cublas_error2 < tol);
                status += ! okay;

                printf("%5d       %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %7.2f (%7.2f)   %15.2e   %15.2e   %15.2e   %15.2e   %s\n",
                       (int)batchCount, (int) M, (int) N,
                       magma_tally3_perf,  1000.*magma_tally3_time,
                       cublas_perf, 1000.*cublas_time,
                       cpu_perf,    1000.*cpu_time,
                       magma_tally3_error, magma_tally3_error2,
                       cublas_error, cublas_error2,
                       (okay ? "ok" : "failed") );
            }
            else {
                printf("%5d       %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)      ---   (  ---  )   ---  \n",
                       (int)batchCount, (int) M, (int) N,
                       magma_tally3_perf,  1000.*magma_tally3_time,
                       cublas_perf, 1000.*cublas_time );
            }
            
            TESTING_FREE_CPU( tau    );
            TESTING_FREE_CPU( h_A    );
            TESTING_FREE_CPU( h_Amagma_tally3);
            TESTING_FREE_CPU( h_work );
            TESTING_FREE_PIN( h_R    );
        
            TESTING_FREE_DEV( d_A   );
            TESTING_FREE_DEV( dtau_magma_tally3  );
            TESTING_FREE_DEV( dtau_cublas );

            TESTING_FREE_DEV( dinfo_magma_tally3 );
            TESTING_FREE_DEV( dinfo_cublas );

            TESTING_FREE_DEV( dA_array   );
            TESTING_FREE_DEV( dtau_array  );

            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    
    magma_tally3_queue_destroy( stream[0] );
    magma_tally3_queue_destroy( stream[1] );

    TESTING_FINALIZE();
    return status;
}
