/*
    -- MAGMA_tally3 (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @generated from testing_ztrsm_batched.cpp normal z -> s, Fri May  1 21:31:40 2015
       @author Chongxiao Cao
       @author Tingxing Dong
       @author Azzam Haidar
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <cuda.h>  // for CUDA_VERSION

// includes, project
#include "testings.h"  // before magma_tally3.h, to include cublas_v2
#include "flops.h"
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "batched_kernel_param.h"

#define h_A(i,j,s) (h_A + (i) + (j)*lda + (s)*lda*Ak)


//#define PRINTMAT

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing strsm_batched
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, magma_tally3_perf, magma_tally3_time=0, cublas_perf=0, cublas_time=0, cpu_perf=0, cpu_time=0;
    float          magma_tally3_error, cublas_error, lapack_error, work[1];
    magma_tally3_int_t M, N, info;
    magma_tally3_int_t Ak;
    magma_tally3_int_t sizeA, sizeB;
    magma_tally3_int_t lda, ldb, ldda, lddb;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t *ipiv;

    float c_zero = MAGMA_tally3_S_ZERO;
    
    float *h_A, *h_B, *h_Bcublas, *h_Bmagma_tally3, *h_Blapack, *h_X;
    float *d_A, *d_B;
    float **d_A_array = NULL;
    float **d_B_array = NULL;
    
    float **dW1_displ  = NULL;
    float **dW2_displ  = NULL;
    float **dW3_displ  = NULL;
    float **dW4_displ  = NULL;
    float **dinvA_array = NULL;
    float **dwork_array = NULL;

    float c_neg_one = MAGMA_tally3_S_NEG_ONE;
    float c_one = MAGMA_tally3_S_ONE;
    float alpha = MAGMA_tally3_S_MAKE(  0.29, -0.86 );
    magma_tally3_int_t status = 0;
    magma_tally3_int_t batchCount = 1;

    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    batchCount = opts.batchcount;
    
    float tol = opts.tolerance * lapackf77_slamch("E");

    printf("side = %s, uplo = %s, transA = %s, diag = %s \n",
           lapack_side_const_tally3(opts.side), lapack_uplo_const_tally3(opts.uplo),
           lapack_trans_const_tally3(opts.transA), lapack_diag_const_tally3(opts.diag) );
    printf("BatchCount   M     N   MAGMA_tally3 Gflop/s (ms)   CUBLAS Gflop/s (ms)   CPU Gflop/s (ms)      MAGMA_tally3     CUBLAS   LAPACK error\n");
    printf("============================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            gflops = FLOPS_STRSM(opts.side, M, N) / 1e9 * batchCount;

            if ( opts.side == Magma_tally3Left ) {
                lda = M;
                Ak  = M;
            } else {
                lda = N;
                Ak  = N;
            }
            
            ldb = M;
            
            ldda = lda = magma_tally3_roundup( lda, opts.roundup );  // multiple of 32 by default
            lddb = ldb = magma_tally3_roundup( ldb, opts.roundup );  // multiple of 32 by default

            sizeA = lda*Ak*batchCount;
            sizeB = ldb*N*batchCount;
            magma_tally3_int_t NN = ldb*N;

            TESTING_MALLOC_CPU( h_A,       float, sizeA  );
            TESTING_MALLOC_CPU( h_B,       float, sizeB   );
            TESTING_MALLOC_CPU( h_X,       float, sizeB   );
            TESTING_MALLOC_CPU( h_Blapack, float, sizeB   );
            TESTING_MALLOC_CPU( h_Bcublas, float, sizeB   );
            TESTING_MALLOC_CPU( h_Bmagma_tally3,  float, sizeB   );
            TESTING_MALLOC_CPU( ipiv,      magma_tally3_int_t,        Ak      );
            
            TESTING_MALLOC_DEV( d_A,       float, ldda*Ak*batchCount );
            TESTING_MALLOC_DEV( d_B,       float, lddb*N*batchCount  );
            
            magma_tally3_malloc((void**)&d_A_array, batchCount * sizeof(*d_A_array));
            magma_tally3_malloc((void**)&d_B_array, batchCount * sizeof(*d_B_array));

            magma_tally3_malloc((void**)&dW1_displ,  batchCount * sizeof(*dW1_displ));
            magma_tally3_malloc((void**)&dW2_displ,  batchCount * sizeof(*dW2_displ));
            magma_tally3_malloc((void**)&dW3_displ,  batchCount * sizeof(*dW3_displ));
            magma_tally3_malloc((void**)&dW4_displ,  batchCount * sizeof(*dW4_displ));
            magma_tally3_malloc((void**)&dinvA_array, batchCount * sizeof(*dinvA_array));
            magma_tally3_malloc((void**)&dwork_array, batchCount * sizeof(*dwork_array));

            float* dinvA=NULL;
            float* dwork=NULL; // invA and work are workspace in strsm
 
            magma_tally3_int_t dinvA_batchSize = magma_tally3_roundup( Ak, TRI_NB )*TRI_NB;
            magma_tally3_int_t dwork_batchSize = lddb*N;
            magma_tally3_smalloc( &dinvA, dinvA_batchSize * batchCount);
            magma_tally3_smalloc( &dwork, dwork_batchSize * batchCount );
    
            sset_pointer(dwork_array, dwork, lddb, 0, 0, dwork_batchSize, batchCount, opts.queue);
            sset_pointer(dinvA_array, dinvA, magma_tally3_roundup( Ak, TRI_NB ), 0, 0, dinvA_batchSize, batchCount, opts.queue);

            memset(h_Bmagma_tally3, 0, batchCount*ldb*N*sizeof(float));
            magma_tally3blas_slaset( Magma_tally3Full, lddb, N*batchCount, c_zero, c_zero, dwork, lddb);

            /* Initialize the matrices */
            /* Factor A into LU to get well-conditioned triangular matrix.
             * Copy L to U, since L seems okay when used with non-unit diagonal
             * (i.e., from U), while U fails when used with unit diagonal. */
            lapackf77_slarnv( &ione, ISEED, &sizeA, h_A );

            for (int s=0; s < batchCount; s++) {
                lapackf77_sgetrf( &Ak, &Ak, h_A + s * lda * Ak, &lda, ipiv, &info );
                for( int j = 0; j < Ak; ++j ) {
                    for( int i = 0; i < j; ++i ) {
                        *h_A(i,j,s) = *h_A(j,i,s);
                    }
                }
            }

            lapackf77_slarnv( &ione, ISEED, &sizeB, h_B );
            memcpy( h_Blapack, h_B, sizeB*sizeof(float) );
            
            /* =====================================================================
               Performs operation using MAGMA_tally3BLAS
               =================================================================== */
            magma_tally3_ssetmatrix( Ak, Ak*batchCount, h_A, lda, d_A, ldda );
            magma_tally3_ssetmatrix( M,  N*batchCount, h_B, ldb, d_B, lddb );

            sset_pointer(d_A_array, d_A, ldda, 0, 0, ldda*Ak, batchCount, opts.queue);
            sset_pointer(d_B_array, d_B, lddb, 0, 0, lddb*N, batchCount, opts.queue);
            sset_pointer(dwork_array, dwork, lddb, 0, 0, lddb*N, batchCount, opts.queue);

            magma_tally3_time = magma_tally3_sync_wtime( NULL );
            #if 1
                magma_tally3blas_strsm_outofplace_batched(
                    opts.side, opts.uplo, opts.transA, opts.diag, 1,
                    M, N, alpha,
                    d_A_array,    ldda, // dA
                    d_B_array,    lddb, // dB
                    dwork_array,  lddb, // dX output
                    dinvA_array,  dinvA_batchSize,
                    dW1_displ,   dW2_displ,
                    dW3_displ,   dW4_displ,
                    1, batchCount, opts.queue);
                magma_tally3_time = magma_tally3_sync_wtime( NULL ) - magma_tally3_time;
                magma_tally3_perf = gflops / magma_tally3_time;
                magma_tally3_sgetmatrix( M, N*batchCount, dwork, lddb, h_Bmagma_tally3, ldb );
            #else
                magma_tally3blas_strsm_batched(
                    opts.side, opts.uplo, opts.transA, opts.diag,
                    M, N, alpha,
                    d_A_array, ldda,
                    d_B_array, lddb,
                    batchCount, opts.queue );
                magma_tally3_time = magma_tally3_sync_wtime( NULL ) - magma_tally3_time;
                magma_tally3_perf = gflops / magma_tally3_time;
                magma_tally3_sgetmatrix( M, N*batchCount, d_B, lddb, h_Bmagma_tally3, ldb );
            #endif
       
            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            magma_tally3_ssetmatrix( M, N*batchCount, h_B, ldb, d_B, lddb );
            sset_pointer(d_B_array, d_B, lddb, 0, 0, lddb*N, batchCount, opts.queue);

            // CUBLAS version <= 6.0 has float **            dA_array, no cast needed.
            // CUBLAS version    6.5 has float const**       dA_array, requiring cast.
            // Correctly, it should be   float const* const* dA_array, to avoid requiring cast.
            #if CUDA_VERSION >= 6050
                cublas_time = magma_tally3_sync_wtime( NULL );
                cublasStrsmBatched(
                    opts.handle, cublas_side_const_tally3(opts.side), cublas_uplo_const_tally3(opts.uplo),
                    cublas_trans_const_tally3(opts.transA), cublas_diag_const_tally3(opts.diag),
                    M, N, &alpha,
                    (const float**) d_A_array, ldda,
                    d_B_array, lddb, batchCount);
                cublas_time = magma_tally3_sync_wtime( NULL ) - cublas_time;
                cublas_perf = gflops / cublas_time;
            #endif

            magma_tally3_sgetmatrix( M, N*batchCount, d_B, lddb, h_Bcublas, ldb );
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally3_wtime();
                for (int s=0; s < batchCount; s++) {
                    blasf77_strsm(
                        lapack_side_const_tally3(opts.side), lapack_uplo_const_tally3(opts.uplo),
                        lapack_trans_const_tally3(opts.transA), lapack_diag_const_tally3(opts.diag),
                        &M, &N, &alpha,
                        h_A       + s * lda * Ak, &lda,
                        h_Blapack + s * ldb * N,  &ldb );
                }
                cpu_time = magma_tally3_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            // ||b - 1/alpha*A*x|| / (||A||*||x||)
            float alpha2 = MAGMA_tally3_S_DIV( c_one, alpha );
            float normR, normX, normA=0.;
            magma_tally3_error  = 0.0;
            cublas_error = 0.0;

            memcpy( h_X, h_Bmagma_tally3, sizeB*sizeof(float) );

            // check magma_tally3
            for (int s=0; s < batchCount; s++) {
                normA = lapackf77_slange( "M", &Ak, &Ak, h_A + s * lda * Ak, &lda, work );
                blasf77_strmm(
                    lapack_side_const_tally3(opts.side), lapack_uplo_const_tally3(opts.uplo),
                    lapack_trans_const_tally3(opts.transA), lapack_diag_const_tally3(opts.diag),
                    &M, &N, &alpha2,
                    h_A + s * lda * Ak, &lda,
                    h_X + s * ldb * N,  &ldb );

                blasf77_saxpy( &NN, &c_neg_one, h_B + s * ldb * N, &ione, h_X + s * ldb * N, &ione );

                normR = lapackf77_slange( "M", &M, &N, h_X + s * ldb * N,      &ldb, work );
                normX = lapackf77_slange( "M", &M, &N, h_Bmagma_tally3 + s * ldb * N, &ldb, work );
                float magma_tally3_err = normR/(normX*normA);

                if ( isnan(magma_tally3_err) || isinf(magma_tally3_err) ) {
                    printf("error for matrix %d magma_tally3_error = %7.2f where normR=%7.2f normX=%7.2f and normA=%7.2f\n", s, magma_tally3_err, normR, normX, normA);
                    magma_tally3_error = magma_tally3_err;
                    break;
                }
                magma_tally3_error = max(fabs(magma_tally3_err), magma_tally3_error);
            }

            memcpy( h_X, h_Bcublas, sizeB*sizeof(float) );
            // check cublas
            #if CUDA_VERSION >= 6050
            for (int s=0; s < batchCount; s++) {
                normA = lapackf77_slange( "M", &Ak, &Ak, h_A + s * lda * Ak, &lda, work );
                blasf77_strmm(
                    lapack_side_const_tally3(opts.side), lapack_uplo_const_tally3(opts.uplo),
                    lapack_trans_const_tally3(opts.transA), lapack_diag_const_tally3(opts.diag),
                    &M, &N, &alpha2,
                    h_A + s * lda * Ak, &lda,
                    h_X + s * ldb * N, &ldb );

                blasf77_saxpy( &NN, &c_neg_one, h_B + s * ldb * N, &ione, h_X  + s * ldb * N, &ione );
                normR = lapackf77_slange( "M", &M, &N, h_X  + s * ldb * N,       &ldb, work );
                normX = lapackf77_slange( "M", &M, &N, h_Bcublas  + s * ldb * N, &ldb, work );
                float cublas_err = normR/(normX*normA);

                if ( isnan(cublas_err) || isinf(cublas_err) ) {
                    printf("error for matrix %d cublas_error = %7.2f where normR=%7.2f normX=%7.2f and normA=%7.2f\n", s, cublas_err, normR, normX, normA);
                    cublas_error = cublas_err;
                    break;
                }
                cublas_error = max(fabs(cublas_err), cublas_error);
            }
            #endif

            if ( opts.lapack ) {
                // check lapack
                // this verifies that the matrix wasn't so bad that it couldn't be solved accurately.
                lapack_error = 0.0;
                memcpy( h_X, h_Blapack, sizeB*sizeof(float) );
                for (int s=0; s < batchCount; s++) {
                    blasf77_strmm(
                        lapack_side_const_tally3(opts.side), lapack_uplo_const_tally3(opts.uplo),
                        lapack_trans_const_tally3(opts.transA), lapack_diag_const_tally3(opts.diag),
                        &M, &N, &alpha2,
                        h_A + s * lda * Ak, &lda,
                        h_X + s * ldb * N,  &ldb );
    
                    blasf77_saxpy( &NN, &c_neg_one, h_B + s * ldb * N, &ione, h_X + s * ldb * N, &ione );
                    normR = lapackf77_slange( "M", &M, &N, h_X + s * ldb * N,       &ldb, work );
                    normX = lapackf77_slange( "M", &M, &N, h_Blapack + s * ldb * N, &ldb, work );
                    float lapack_err = normR/(normX*normA);

                    if (lapack_error < lapack_err)
                        lapack_error = lapack_err;
                }

                printf("%5d     %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %8.2e   %8.2e   %s\n",
                        (int)batchCount, (int) M, (int) N,
                        magma_tally3_perf,  1000.*magma_tally3_time,
                        cublas_perf, 1000.*cublas_time,
                        cpu_perf,    1000.*cpu_time,
                        magma_tally3_error, cublas_error, lapack_error,
                        (magma_tally3_error < tol && cublas_error < tol? "ok" : "failed"));
                status += ! (magma_tally3_error < tol && cublas_error < tol);
            }
            else {
                printf("%5d     %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)     ---   (  ---  )   %8.2e   %8.2e     ---      %s\n",
                        (int)batchCount, (int) M, (int) N,
                        magma_tally3_perf,  1000.*magma_tally3_time,
                        cublas_perf, 1000.*cublas_time,
                        magma_tally3_error, cublas_error,
                        (magma_tally3_error < tol && cublas_error < tol ? "ok" : "failed"));
                status += ! (magma_tally3_error < tol && cublas_error < tol);
            }
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_B );
            TESTING_FREE_CPU( h_X );
            TESTING_FREE_CPU( h_Blapack );
            TESTING_FREE_CPU( h_Bcublas );
            TESTING_FREE_CPU( h_Bmagma_tally3  );
            TESTING_FREE_CPU( ipiv );
            
            TESTING_FREE_DEV( d_A );
            TESTING_FREE_DEV( d_B );
            magma_tally3_free(d_A_array);
            magma_tally3_free(d_B_array);

            magma_tally3_free(dW1_displ);
            magma_tally3_free(dW2_displ);
            magma_tally3_free(dW3_displ);
            magma_tally3_free(dW4_displ);

            TESTING_FREE_DEV( dinvA );
            TESTING_FREE_DEV( dwork );
            magma_tally3_free(dwork_array);
            magma_tally3_free(dinvA_array);
            
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
