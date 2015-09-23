/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zpotrf_mgpu.cpp normal z -> s, Fri Jan 30 19:00:24 2015
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing spotrf_mgpu
*/
int main( int argc, char** argv )
{
    TESTING_INIT();

    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf=0, cpu_time=0;
    float           error, work[1];
    float c_neg_one = MAGMA_minproduct_S_NEG_ONE;
    float *h_A, *h_R;
    magma_minproductFloat_ptr d_lA[ Magma_minproductMaxGPUs ];
    magma_minproduct_int_t N, n2, lda, ldda, max_size, ngpu;
    magma_minproduct_int_t info, nb;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t  status = 0;
    
    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    float tol = opts.tolerance * lapackf77_slamch("E");
    
    printf("ngpu = %d, uplo = %s\n", (int) opts.ngpu, lapack_uplo_const(opts.uplo) );
    printf("    N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R||_F / ||A||_F\n");
    printf("=================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            n2     = lda*N;
            nb     = magma_minproduct_get_spotrf_nb( N );
            gflops = FLOPS_SPOTRF( N ) / 1e9;
            
            // ngpu must be at least the number of blocks
            ngpu = min( opts.ngpu, int((N+nb-1)/nb) );
            if ( ngpu < opts.ngpu ) {
                printf( " * too many GPUs for the matrix size, using %d GPUs\n", (int) ngpu );
            }
            
            // Allocate host memory for the matrix
            TESTING_MALLOC_CPU( h_A, float, n2 );
            TESTING_MALLOC_PIN( h_R, float, n2 );
            
            // Allocate device memory
            // matrix is distributed by block-rows or block-columns
            // this is maximum size that any GPU stores;
            // size is rounded up to full blocks in both rows and columns
            max_size = nb*(1+N/(nb*ngpu)) * nb*((N+nb-1)/nb);
            for( int dev=0; dev < ngpu; dev++ ) {
                magma_minproduct_setdevice( dev );
                TESTING_MALLOC_DEV( d_lA[dev], float, max_size );
            }
            
            /* Initialize the matrix */
            lapackf77_slarnv( &ione, ISEED, &n2, h_A );
            magma_minproduct_smake_hpd( N, h_A, lda );
            lapackf77_slacpy( Magma_minproductUpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_minproduct_wtime();
                lapackf77_spotrf( lapack_uplo_const(opts.uplo), &N, h_A, &lda, &info );
                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_spotrf returned error %d: %s.\n",
                           (int) info, magma_minproduct_strerror( info ));
            }
            
            /* ====================================================================
               Performs operation using MAGMA_minproduct
               =================================================================== */
            if ( opts.uplo == Magma_minproductUpper ) {
                ldda = ((N+nb-1)/nb)*nb;
                magma_minproduct_ssetmatrix_1D_col_bcyclic( N, N, h_R, lda, d_lA, ldda, ngpu, nb );
            } else {
                ldda = (1+N/(nb*ngpu))*nb;
                magma_minproduct_ssetmatrix_1D_row_bcyclic( N, N, h_R, lda, d_lA, ldda, ngpu, nb );
            }
            
            gpu_time = magma_minproduct_wtime();
            magma_minproduct_spotrf_mgpu( ngpu, opts.uplo, N, d_lA, ldda, &info );
            gpu_time = magma_minproduct_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_minproduct_spotrf_mgpu returned error %d: %s.\n",
                       (int) info, magma_minproduct_strerror( info ));
            
            if ( opts.uplo == Magma_minproductUpper ) {
                magma_minproduct_sgetmatrix_1D_col_bcyclic( N, N, d_lA, ldda, h_R, lda, ngpu, nb );
            } else {
                magma_minproduct_sgetmatrix_1D_row_bcyclic( N, N, d_lA, ldda, h_R, lda, ngpu, nb );
            }
            
            /* =====================================================================
               Check the result compared to LAPACK
               =================================================================== */
            for( int dev=0; dev < ngpu; dev++ ){
                magma_minproduct_setdevice( dev );
                magma_minproduct_device_sync();
            }
            if ( opts.lapack ) {
                error = lapackf77_slange("f", &N, &N, h_A, &lda, work );
                blasf77_saxpy( &n2, &c_neg_one, h_A, &ione, h_R, &ione );
                error = lapackf77_slange("f", &N, &N, h_R, &lda, work ) / error;
                
                printf("%5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                       (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time,
                       error, (error < tol ? "ok" : "failed") );
                status += ! (error < tol);
            }
            else {
                printf("%5d     ---   (  ---  )   %7.2f (%7.2f)     ---\n",
                       (int) N, gpu_perf, gpu_time );
            }
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_PIN( h_R );
            for( int dev=0; dev < ngpu; dev++ ){
                magma_minproduct_setdevice( dev );
                TESTING_FREE_DEV( d_lA[dev] );
            }
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
