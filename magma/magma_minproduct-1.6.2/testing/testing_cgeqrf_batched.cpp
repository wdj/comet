/*
    -- MAGMA_minproduct (version 1.6) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @author Tingxing Dong
       @author Azzam Haidar

       @generated from testing_zgeqrf_batched.cpp normal z -> c, Fri May  1 21:31:38 2015

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <cuda.h>  // for CUDA_VERSION

// includes, project
#include "flops.h"
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "testings.h"
#include "batched_kernel_param.h"


void get_QR_error(magma_minproduct_int_t M, magma_minproduct_int_t N, magma_minproduct_int_t min_mn,
                    magma_minproductFloatComplex *h_R,  magma_minproductFloatComplex *h_A, magma_minproduct_int_t lda,
                    magma_minproductFloatComplex *tau,
                    magma_minproductFloatComplex *Q,  magma_minproduct_int_t ldq,
                    magma_minproductFloatComplex *R,  magma_minproduct_int_t ldr,
                    magma_minproductFloatComplex *h_work,  magma_minproduct_int_t lwork,
                    float *work, float *error, float *error2)
{
    /* h_R:input the factorized matrix by lapack QR,
       h_A:input the original matrix copy
       tau: input
    */
    
    const float             d_neg_one = MAGMA_minproduct_D_NEG_ONE;
    const float             d_one     = MAGMA_minproduct_D_ONE;
    const magma_minproductFloatComplex c_neg_one = MAGMA_minproduct_C_NEG_ONE;
    const magma_minproductFloatComplex c_one     = MAGMA_minproduct_C_ONE;
    const magma_minproductFloatComplex c_zero    = MAGMA_minproduct_C_ZERO;
    float           Anorm;
    
    magma_minproduct_int_t info;
    
    // generate M by K matrix Q, where K = min(M,N)
    lapackf77_clacpy( "Lower", &M, &min_mn, h_R, &lda, Q, &ldq );
    lapackf77_cungqr( &M, &min_mn, &min_mn, Q, &ldq, tau, h_work, &lwork, &info );
    assert( info == 0 );
    
    // copy K by N matrix R
    lapackf77_claset( "Lower", &min_mn, &N, &c_zero, &c_zero, R, &ldr );
    lapackf77_clacpy( "Upper", &min_mn, &N, h_R, &lda,        R, &ldr );
    
    // error = || R - Q^H*A || / (N * ||A||)
    blasf77_cgemm( "Conj", "NoTrans", &min_mn, &N, &M,
    &c_neg_one, Q, &ldq, h_A, &lda, &c_one, R, &ldr );
    
    Anorm = lapackf77_clange( "1", &M,      &N, h_A, &lda, work );
    *error = lapackf77_clange( "1", &min_mn, &N, R,   &ldr, work );
    
    if ( N > 0 && Anorm > 0 )
        *error /= (N*Anorm);
    
    // set R = I (K by K identity), then R = I - Q^H*Q
    // error = || I - Q^H*Q || / N
    lapackf77_claset( "Upper", &min_mn, &min_mn, &c_zero, &c_one, R, &ldr );
    blasf77_cherk( "Upper", "Conj", &min_mn, &M, &d_neg_one, Q, &ldq, &d_one, R, &ldr );
    *error2 = lapackf77_clanhe( "1", "Upper", &min_mn, R, &ldr, work );
    if ( N > 0 )
        *error2 /= N;
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cgeqrf_batched
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gflops, magma_minproduct_perf, magma_minproduct_time, cublas_perf=0, cublas_time=0, cpu_perf, cpu_time;
    float           magma_minproduct_error=0.0, cublas_error=0.0, magma_minproduct_error2=0.0, cublas_error2=0.0, error, error2;

    magma_minproductFloatComplex *h_A, *h_R, *h_Amagma_minproduct, *tau, *h_work, tmp[1];
    magma_minproductFloatComplex *d_A, *dtau_magma_minproduct, *dtau_cublas;

    magma_minproductFloatComplex **dA_array = NULL;
    magma_minproductFloatComplex **dtau_array = NULL;

    magma_minproduct_int_t   *dinfo_magma_minproduct, *dinfo_cublas;

    magma_minproduct_int_t M, N, lda, ldda, lwork, n2, info, min_mn;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t status = 0;

    magma_minproduct_int_t batchCount = 1;
    magma_minproduct_int_t column;

    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    batchCount = opts.batchcount;

    float tol = opts.tolerance * lapackf77_slamch("E");
    
    magma_minproduct_queue_t stream[2];
    magma_minproduct_queue_create( &stream[0] );
    magma_minproduct_queue_create( &stream[1] );
    printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |R - Q^H*A|   |I - Q^H*Q|\n");
    printf("BatchCount     M    N     MAGMA_minproduct GFlop/s (ms)   CUBLAS GFlop/s (ms)   CPU GFlop/s (ms)   |R - Q^H*A|_mag   |I - Q^H*Q|_mag   |R - Q^H*A|_cub   |I - Q^H*Q|_cub \n");
    printf("============================================================================================================================================================= \n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M     = opts.msize[itest];
            N     = opts.nsize[itest];
            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N * batchCount;
            ldda = M;
            ldda   = magma_minproduct_roundup( M, opts.roundup );  // multiple of 32 by default

            gflops = (FLOPS_CGEQRF( M, N ) + FLOPS_CGEQRT( M, N )) / 1e9 * batchCount;

            /* Allocate memory for the matrix */
            TESTING_MALLOC_CPU( tau,   magma_minproductFloatComplex, min_mn * batchCount );
            TESTING_MALLOC_CPU( h_A,   magma_minproductFloatComplex, n2     );
            TESTING_MALLOC_CPU( h_Amagma_minproduct,   magma_minproductFloatComplex, n2     );
            TESTING_MALLOC_PIN( h_R,   magma_minproductFloatComplex, n2     );
        
            TESTING_MALLOC_DEV( d_A,   magma_minproductFloatComplex, ldda*N * batchCount );

            TESTING_MALLOC_DEV( dtau_magma_minproduct,  magma_minproductFloatComplex, min_mn * batchCount);
            TESTING_MALLOC_DEV( dtau_cublas,  magma_minproductFloatComplex, min_mn * batchCount);

            TESTING_MALLOC_DEV(  dinfo_magma_minproduct,  magma_minproduct_int_t, batchCount);
            TESTING_MALLOC_DEV(  dinfo_cublas,  magma_minproduct_int_t, batchCount);

            magma_minproduct_malloc((void**)&dA_array, batchCount * sizeof(*dA_array));
            magma_minproduct_malloc((void**)&dtau_array, batchCount * sizeof(*dtau_array));
        
            // to determine the size of lwork
            lwork = -1;
            lapackf77_cgeqrf(&M, &N, NULL, &M, NULL, tmp, &lwork, &info);
            lwork = (magma_minproduct_int_t)MAGMA_minproduct_C_REAL( tmp[0] );
            lwork = max(lwork, N*N);
           
            TESTING_MALLOC_CPU( h_work, magma_minproductFloatComplex, lwork * batchCount);

            column = N * batchCount;
            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &n2, h_A );
            for(int i=0; i<batchCount; i++)
			{
				magma_minproduct_cmake_hpd( N, h_A + i * lda * N, lda );// need modification
			}
			
			lapackf77_clacpy( Magma_minproductUpperLowerStr, &M, &column, h_A, &lda, h_R, &lda );
       
            /* ====================================================================
               Performs operation using MAGMA_minproduct
               =================================================================== */
            magma_minproduct_csetmatrix( M, column, h_R, lda,  d_A, ldda );
            cset_pointer(dA_array, d_A, 1, 0, 0, ldda*N, batchCount, opts.queue);
            cset_pointer(dtau_array, dtau_magma_minproduct, 1, 0, 0, min_mn, batchCount, opts.queue);
    
            magma_minproduct_time = magma_minproduct_sync_wtime(0);
    
            info = magma_minproduct_cgeqrf_batched(M, N, dA_array, ldda, dtau_array, dinfo_magma_minproduct, batchCount, opts.queue);

            magma_minproduct_time = magma_minproduct_sync_wtime(0) - magma_minproduct_time;
            magma_minproduct_perf = gflops / magma_minproduct_time;

            magma_minproduct_cgetmatrix( M, column, d_A, ldda, h_Amagma_minproduct, lda);

            if (info != 0)
                printf("magma_minproduct_cgeqrf_batched returned error %d: %s.\n",
                       (int) info, magma_minproduct_strerror( info ));
            
            /* ====================================================================
               Performs operation using CUBLAS
               =================================================================== */

            /* cublasCgeqrfBatched is only available from CUBLAS v6.5 */
            #if CUDA_VERSION >= 6050
            magma_minproduct_csetmatrix( M, column, h_R, lda,  d_A, ldda );
            cset_pointer(dA_array, d_A, 1, 0, 0, ldda*N, batchCount, opts.queue);
            cset_pointer(dtau_array, dtau_cublas, 1, 0, 0, min_mn, batchCount, opts.queue);

            cublas_time = magma_minproduct_sync_wtime(0);
    
            int cublas_info;  // int, not magma_minproduct_int_t
            cublasCgeqrfBatched(opts.handle, M, N, dA_array, ldda, dtau_array, &cublas_info, batchCount);

            cublas_time = magma_minproduct_sync_wtime(0) - cublas_time;
            cublas_perf = gflops / cublas_time;

            if (cublas_info != 0)
                printf("cublasCgeqrfBatched returned error %d: %s.\n",
                       (int) cublas_info, magma_minproduct_strerror( cublas_info ));
            #endif

            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.check ) {
                cpu_time = magma_minproduct_wtime();

                for (int i=0; i < batchCount; i++)
                {
                    lapackf77_cgeqrf(&M, &N, h_A + i*lda*N, &lda, tau + i*min_mn, h_work + i * lwork, &lwork, &info);
                  
                    /* lapackf77_clarft( Magma_minproductForwardStr, Magma_minproductColumnwiseStr,
                                     &M, &N, h_A + i*lda*N, &lda, tau + i*min_mn, h_work + i * lwork, &N); */
                }

                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_cgeqrf returned error %d: %s.\n",
                           (int) info, magma_minproduct_strerror( info ));
                
                /* =====================================================================
                   Check the MAGMA_minproduct CUBLAS result compared to LAPACK
                   =================================================================== */
                magma_minproduct_int_t ldq = M;
                magma_minproduct_int_t ldr = min_mn;
                magma_minproductFloatComplex *Q, *R;
                float *work;

                TESTING_MALLOC_CPU( Q,    magma_minproductFloatComplex, ldq*min_mn );  // M by K
                TESTING_MALLOC_CPU( R,    magma_minproductFloatComplex, ldr*N );       // K by N
                TESTING_MALLOC_CPU( work, float,             min_mn );

                /* check magma_minproduct result */
                magma_minproduct_cgetvector(min_mn*batchCount, dtau_magma_minproduct, 1, tau, 1);
                for (int i=0; i < batchCount; i++)
                {
                    get_QR_error(M, N, min_mn,
                             h_Amagma_minproduct + i*lda*N, h_R + i*lda*N, lda, tau + i*min_mn,
                             Q, ldq, R, ldr, h_work, lwork,
                             work, &error, &error2);

                    if ( isnan(error) || isinf(error) ) {
                        magma_minproduct_error = error;
                        break;
                    }
                    magma_minproduct_error  = max( fabs(error),  magma_minproduct_error  );
                    magma_minproduct_error2 = max( fabs(error2), magma_minproduct_error2 );
                }

                /* check cublas result */
                #if CUDA_VERSION >= 6050
                magma_minproduct_cgetvector(min_mn*batchCount, dtau_magma_minproduct, 1, tau, 1);
                magma_minproduct_cgetmatrix( M, column, d_A, ldda, h_A, lda);
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

                bool okay = (magma_minproduct_error < tol && magma_minproduct_error2 < tol);
                //bool okay_cublas = (cublas_error < tol && cublas_error2 < tol);
                status += ! okay;

                printf("%5d       %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %7.2f (%7.2f)   %15.2e   %15.2e   %15.2e   %15.2e   %s\n",
                       (int)batchCount, (int) M, (int) N,
                       magma_minproduct_perf,  1000.*magma_minproduct_time,
                       cublas_perf, 1000.*cublas_time,
                       cpu_perf,    1000.*cpu_time,
                       magma_minproduct_error, magma_minproduct_error2,
                       cublas_error, cublas_error2,
                       (okay ? "ok" : "failed") );
            }
            else {
                printf("%5d       %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)      ---   (  ---  )   ---  \n",
                       (int)batchCount, (int) M, (int) N,
                       magma_minproduct_perf,  1000.*magma_minproduct_time,
                       cublas_perf, 1000.*cublas_time );
            }
            
            TESTING_FREE_CPU( tau    );
            TESTING_FREE_CPU( h_A    );
            TESTING_FREE_CPU( h_Amagma_minproduct);
            TESTING_FREE_CPU( h_work );
            TESTING_FREE_PIN( h_R    );
        
            TESTING_FREE_DEV( d_A   );
            TESTING_FREE_DEV( dtau_magma_minproduct  );
            TESTING_FREE_DEV( dtau_cublas );

            TESTING_FREE_DEV( dinfo_magma_minproduct );
            TESTING_FREE_DEV( dinfo_cublas );

            TESTING_FREE_DEV( dA_array   );
            TESTING_FREE_DEV( dtau_array  );

            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    
    magma_minproduct_queue_destroy( stream[0] );
    magma_minproduct_queue_destroy( stream[1] );

    TESTING_FINALIZE();
    return status;
}
