/*
    -- MAGMA_minproduct (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @generated from testing_ztrsm_batched.cpp normal z -> c, Fri May  1 21:31:38 2015
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
#include "testings.h"  // before magma_minproduct.h, to include cublas_v2
#include "flops.h"
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "batched_kernel_param.h"

#define h_A(i,j,s) (h_A + (i) + (j)*lda + (s)*lda*Ak)


//#define PRINTMAT

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing ctrsm_batched
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, magma_minproduct_perf, magma_minproduct_time=0, cublas_perf=0, cublas_time=0, cpu_perf=0, cpu_time=0;
    float          magma_minproduct_error, cublas_error, lapack_error, work[1];
    magma_minproduct_int_t M, N, info;
    magma_minproduct_int_t Ak;
    magma_minproduct_int_t sizeA, sizeB;
    magma_minproduct_int_t lda, ldb, ldda, lddb;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t *ipiv;

    magma_minproductFloatComplex c_zero = MAGMA_minproduct_C_ZERO;
    
    magma_minproductFloatComplex *h_A, *h_B, *h_Bcublas, *h_Bmagma_minproduct, *h_Blapack, *h_X;
    magma_minproductFloatComplex *d_A, *d_B;
    magma_minproductFloatComplex **d_A_array = NULL;
    magma_minproductFloatComplex **d_B_array = NULL;
    
    magma_minproductFloatComplex **dW1_displ  = NULL;
    magma_minproductFloatComplex **dW2_displ  = NULL;
    magma_minproductFloatComplex **dW3_displ  = NULL;
    magma_minproductFloatComplex **dW4_displ  = NULL;
    magma_minproductFloatComplex **dinvA_array = NULL;
    magma_minproductFloatComplex **dwork_array = NULL;

    magma_minproductFloatComplex c_neg_one = MAGMA_minproduct_C_NEG_ONE;
    magma_minproductFloatComplex c_one = MAGMA_minproduct_C_ONE;
    magma_minproductFloatComplex alpha = MAGMA_minproduct_C_MAKE(  0.29, -0.86 );
    magma_minproduct_int_t status = 0;
    magma_minproduct_int_t batchCount = 1;

    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    batchCount = opts.batchcount;
    
    float tol = opts.tolerance * lapackf77_slamch("E");

    printf("side = %s, uplo = %s, transA = %s, diag = %s \n",
           lapack_side_const(opts.side), lapack_uplo_const(opts.uplo),
           lapack_trans_const(opts.transA), lapack_diag_const(opts.diag) );
    printf("BatchCount   M     N   MAGMA_minproduct Gflop/s (ms)   CUBLAS Gflop/s (ms)   CPU Gflop/s (ms)      MAGMA_minproduct     CUBLAS   LAPACK error\n");
    printf("============================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            gflops = FLOPS_CTRSM(opts.side, M, N) / 1e9 * batchCount;

            if ( opts.side == Magma_minproductLeft ) {
                lda = M;
                Ak  = M;
            } else {
                lda = N;
                Ak  = N;
            }
            
            ldb = M;
            
            ldda = lda = magma_minproduct_roundup( lda, opts.roundup );  // multiple of 32 by default
            lddb = ldb = magma_minproduct_roundup( ldb, opts.roundup );  // multiple of 32 by default

            sizeA = lda*Ak*batchCount;
            sizeB = ldb*N*batchCount;
            magma_minproduct_int_t NN = ldb*N;

            TESTING_MALLOC_CPU( h_A,       magma_minproductFloatComplex, sizeA  );
            TESTING_MALLOC_CPU( h_B,       magma_minproductFloatComplex, sizeB   );
            TESTING_MALLOC_CPU( h_X,       magma_minproductFloatComplex, sizeB   );
            TESTING_MALLOC_CPU( h_Blapack, magma_minproductFloatComplex, sizeB   );
            TESTING_MALLOC_CPU( h_Bcublas, magma_minproductFloatComplex, sizeB   );
            TESTING_MALLOC_CPU( h_Bmagma_minproduct,  magma_minproductFloatComplex, sizeB   );
            TESTING_MALLOC_CPU( ipiv,      magma_minproduct_int_t,        Ak      );
            
            TESTING_MALLOC_DEV( d_A,       magma_minproductFloatComplex, ldda*Ak*batchCount );
            TESTING_MALLOC_DEV( d_B,       magma_minproductFloatComplex, lddb*N*batchCount  );
            
            magma_minproduct_malloc((void**)&d_A_array, batchCount * sizeof(*d_A_array));
            magma_minproduct_malloc((void**)&d_B_array, batchCount * sizeof(*d_B_array));

            magma_minproduct_malloc((void**)&dW1_displ,  batchCount * sizeof(*dW1_displ));
            magma_minproduct_malloc((void**)&dW2_displ,  batchCount * sizeof(*dW2_displ));
            magma_minproduct_malloc((void**)&dW3_displ,  batchCount * sizeof(*dW3_displ));
            magma_minproduct_malloc((void**)&dW4_displ,  batchCount * sizeof(*dW4_displ));
            magma_minproduct_malloc((void**)&dinvA_array, batchCount * sizeof(*dinvA_array));
            magma_minproduct_malloc((void**)&dwork_array, batchCount * sizeof(*dwork_array));

            magma_minproductFloatComplex* dinvA=NULL;
            magma_minproductFloatComplex* dwork=NULL; // invA and work are workspace in ctrsm
 
            magma_minproduct_int_t dinvA_batchSize = magma_minproduct_roundup( Ak, TRI_NB )*TRI_NB;
            magma_minproduct_int_t dwork_batchSize = lddb*N;
            magma_minproduct_cmalloc( &dinvA, dinvA_batchSize * batchCount);
            magma_minproduct_cmalloc( &dwork, dwork_batchSize * batchCount );
    
            cset_pointer(dwork_array, dwork, lddb, 0, 0, dwork_batchSize, batchCount, opts.queue);
            cset_pointer(dinvA_array, dinvA, magma_minproduct_roundup( Ak, TRI_NB ), 0, 0, dinvA_batchSize, batchCount, opts.queue);

            memset(h_Bmagma_minproduct, 0, batchCount*ldb*N*sizeof(magma_minproductFloatComplex));
            magma_minproductblas_claset( Magma_minproductFull, lddb, N*batchCount, c_zero, c_zero, dwork, lddb);

            /* Initialize the matrices */
            /* Factor A into LU to get well-conditioned triangular matrix.
             * Copy L to U, since L seems okay when used with non-unit diagonal
             * (i.e., from U), while U fails when used with unit diagonal. */
            lapackf77_clarnv( &ione, ISEED, &sizeA, h_A );

            for (int s=0; s < batchCount; s++) {
                lapackf77_cgetrf( &Ak, &Ak, h_A + s * lda * Ak, &lda, ipiv, &info );
                for( int j = 0; j < Ak; ++j ) {
                    for( int i = 0; i < j; ++i ) {
                        *h_A(i,j,s) = *h_A(j,i,s);
                    }
                }
            }

            lapackf77_clarnv( &ione, ISEED, &sizeB, h_B );
            memcpy( h_Blapack, h_B, sizeB*sizeof(magma_minproductFloatComplex) );
            
            /* =====================================================================
               Performs operation using MAGMA_minproductBLAS
               =================================================================== */
            magma_minproduct_csetmatrix( Ak, Ak*batchCount, h_A, lda, d_A, ldda );
            magma_minproduct_csetmatrix( M,  N*batchCount, h_B, ldb, d_B, lddb );

            cset_pointer(d_A_array, d_A, ldda, 0, 0, ldda*Ak, batchCount, opts.queue);
            cset_pointer(d_B_array, d_B, lddb, 0, 0, lddb*N, batchCount, opts.queue);
            cset_pointer(dwork_array, dwork, lddb, 0, 0, lddb*N, batchCount, opts.queue);

            magma_minproduct_time = magma_minproduct_sync_wtime( NULL );
            #if 1
                magma_minproductblas_ctrsm_outofplace_batched(
                    opts.side, opts.uplo, opts.transA, opts.diag, 1,
                    M, N, alpha,
                    d_A_array,    ldda, // dA
                    d_B_array,    lddb, // dB
                    dwork_array,  lddb, // dX output
                    dinvA_array,  dinvA_batchSize,
                    dW1_displ,   dW2_displ,
                    dW3_displ,   dW4_displ,
                    1, batchCount, opts.queue);
                magma_minproduct_time = magma_minproduct_sync_wtime( NULL ) - magma_minproduct_time;
                magma_minproduct_perf = gflops / magma_minproduct_time;
                magma_minproduct_cgetmatrix( M, N*batchCount, dwork, lddb, h_Bmagma_minproduct, ldb );
            #else
                magma_minproductblas_ctrsm_batched(
                    opts.side, opts.uplo, opts.transA, opts.diag,
                    M, N, alpha,
                    d_A_array, ldda,
                    d_B_array, lddb,
                    batchCount, opts.queue );
                magma_minproduct_time = magma_minproduct_sync_wtime( NULL ) - magma_minproduct_time;
                magma_minproduct_perf = gflops / magma_minproduct_time;
                magma_minproduct_cgetmatrix( M, N*batchCount, d_B, lddb, h_Bmagma_minproduct, ldb );
            #endif
       
            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            magma_minproduct_csetmatrix( M, N*batchCount, h_B, ldb, d_B, lddb );
            cset_pointer(d_B_array, d_B, lddb, 0, 0, lddb*N, batchCount, opts.queue);

            // CUBLAS version <= 6.0 has magma_minproductFloatComplex **            dA_array, no cast needed.
            // CUBLAS version    6.5 has magma_minproductFloatComplex const**       dA_array, requiring cast.
            // Correctly, it should be   magma_minproductFloatComplex const* const* dA_array, to avoid requiring cast.
            #if CUDA_VERSION >= 6050
                cublas_time = magma_minproduct_sync_wtime( NULL );
                cublasCtrsmBatched(
                    opts.handle, cublas_side_const(opts.side), cublas_uplo_const(opts.uplo),
                    cublas_trans_const(opts.transA), cublas_diag_const(opts.diag),
                    M, N, &alpha,
                    (const magma_minproductFloatComplex**) d_A_array, ldda,
                    d_B_array, lddb, batchCount);
                cublas_time = magma_minproduct_sync_wtime( NULL ) - cublas_time;
                cublas_perf = gflops / cublas_time;
            #endif

            magma_minproduct_cgetmatrix( M, N*batchCount, d_B, lddb, h_Bcublas, ldb );
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_minproduct_wtime();
                for (int s=0; s < batchCount; s++) {
                    blasf77_ctrsm(
                        lapack_side_const(opts.side), lapack_uplo_const(opts.uplo),
                        lapack_trans_const(opts.transA), lapack_diag_const(opts.diag),
                        &M, &N, &alpha,
                        h_A       + s * lda * Ak, &lda,
                        h_Blapack + s * ldb * N,  &ldb );
                }
                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            // ||b - 1/alpha*A*x|| / (||A||*||x||)
            magma_minproductFloatComplex alpha2 = MAGMA_minproduct_C_DIV( c_one, alpha );
            float normR, normX, normA=0.;
            magma_minproduct_error  = 0.0;
            cublas_error = 0.0;

            memcpy( h_X, h_Bmagma_minproduct, sizeB*sizeof(magma_minproductFloatComplex) );

            // check magma_minproduct
            for (int s=0; s < batchCount; s++) {
                normA = lapackf77_clange( "M", &Ak, &Ak, h_A + s * lda * Ak, &lda, work );
                blasf77_ctrmm(
                    lapack_side_const(opts.side), lapack_uplo_const(opts.uplo),
                    lapack_trans_const(opts.transA), lapack_diag_const(opts.diag),
                    &M, &N, &alpha2,
                    h_A + s * lda * Ak, &lda,
                    h_X + s * ldb * N,  &ldb );

                blasf77_caxpy( &NN, &c_neg_one, h_B + s * ldb * N, &ione, h_X + s * ldb * N, &ione );

                normR = lapackf77_clange( "M", &M, &N, h_X + s * ldb * N,      &ldb, work );
                normX = lapackf77_clange( "M", &M, &N, h_Bmagma_minproduct + s * ldb * N, &ldb, work );
                float magma_minproduct_err = normR/(normX*normA);

                if ( isnan(magma_minproduct_err) || isinf(magma_minproduct_err) ) {
                    printf("error for matrix %d magma_minproduct_error = %7.2f where normR=%7.2f normX=%7.2f and normA=%7.2f\n", s, magma_minproduct_err, normR, normX, normA);
                    magma_minproduct_error = magma_minproduct_err;
                    break;
                }
                magma_minproduct_error = max(fabs(magma_minproduct_err), magma_minproduct_error);
            }

            memcpy( h_X, h_Bcublas, sizeB*sizeof(magma_minproductFloatComplex) );
            // check cublas
            #if CUDA_VERSION >= 6050
            for (int s=0; s < batchCount; s++) {
                normA = lapackf77_clange( "M", &Ak, &Ak, h_A + s * lda * Ak, &lda, work );
                blasf77_ctrmm(
                    lapack_side_const(opts.side), lapack_uplo_const(opts.uplo),
                    lapack_trans_const(opts.transA), lapack_diag_const(opts.diag),
                    &M, &N, &alpha2,
                    h_A + s * lda * Ak, &lda,
                    h_X + s * ldb * N, &ldb );

                blasf77_caxpy( &NN, &c_neg_one, h_B + s * ldb * N, &ione, h_X  + s * ldb * N, &ione );
                normR = lapackf77_clange( "M", &M, &N, h_X  + s * ldb * N,       &ldb, work );
                normX = lapackf77_clange( "M", &M, &N, h_Bcublas  + s * ldb * N, &ldb, work );
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
                memcpy( h_X, h_Blapack, sizeB*sizeof(magma_minproductFloatComplex) );
                for (int s=0; s < batchCount; s++) {
                    blasf77_ctrmm(
                        lapack_side_const(opts.side), lapack_uplo_const(opts.uplo),
                        lapack_trans_const(opts.transA), lapack_diag_const(opts.diag),
                        &M, &N, &alpha2,
                        h_A + s * lda * Ak, &lda,
                        h_X + s * ldb * N,  &ldb );
    
                    blasf77_caxpy( &NN, &c_neg_one, h_B + s * ldb * N, &ione, h_X + s * ldb * N, &ione );
                    normR = lapackf77_clange( "M", &M, &N, h_X + s * ldb * N,       &ldb, work );
                    normX = lapackf77_clange( "M", &M, &N, h_Blapack + s * ldb * N, &ldb, work );
                    float lapack_err = normR/(normX*normA);

                    if (lapack_error < lapack_err)
                        lapack_error = lapack_err;
                }

                printf("%5d     %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %8.2e   %8.2e   %s\n",
                        (int)batchCount, (int) M, (int) N,
                        magma_minproduct_perf,  1000.*magma_minproduct_time,
                        cublas_perf, 1000.*cublas_time,
                        cpu_perf,    1000.*cpu_time,
                        magma_minproduct_error, cublas_error, lapack_error,
                        (magma_minproduct_error < tol && cublas_error < tol? "ok" : "failed"));
                status += ! (magma_minproduct_error < tol && cublas_error < tol);
            }
            else {
                printf("%5d     %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)     ---   (  ---  )   %8.2e   %8.2e     ---      %s\n",
                        (int)batchCount, (int) M, (int) N,
                        magma_minproduct_perf,  1000.*magma_minproduct_time,
                        cublas_perf, 1000.*cublas_time,
                        magma_minproduct_error, cublas_error,
                        (magma_minproduct_error < tol && cublas_error < tol ? "ok" : "failed"));
                status += ! (magma_minproduct_error < tol && cublas_error < tol);
            }
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_B );
            TESTING_FREE_CPU( h_X );
            TESTING_FREE_CPU( h_Blapack );
            TESTING_FREE_CPU( h_Bcublas );
            TESTING_FREE_CPU( h_Bmagma_minproduct  );
            TESTING_FREE_CPU( ipiv );
            
            TESTING_FREE_DEV( d_A );
            TESTING_FREE_DEV( d_B );
            magma_minproduct_free(d_A_array);
            magma_minproduct_free(d_B_array);

            magma_minproduct_free(dW1_displ);
            magma_minproduct_free(dW2_displ);
            magma_minproduct_free(dW3_displ);
            magma_minproduct_free(dW4_displ);

            TESTING_FREE_DEV( dinvA );
            TESTING_FREE_DEV( dwork );
            magma_minproduct_free(dwork_array);
            magma_minproduct_free(dinvA_array);
            
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
