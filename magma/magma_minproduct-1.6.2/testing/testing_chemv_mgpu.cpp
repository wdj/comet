/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zhemv_mgpu.cpp normal z -> c, Fri Jan 30 19:00:23 2015
       
       @author Mark Gates
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "testings.h"  // before magma_minproduct.h, to include cublas_v2
#include "flops.h"
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"

#include "magma_minproduct_operators.h"

#define PRECISION_c


// --------------------
int main(int argc, char **argv)
{
    TESTING_INIT();

    real_Double_t gflops, cpu_time=0, cpu_perf=0, gpu_time, gpu_perf, mgpu_time, mgpu_perf, cuda_time, cuda_perf;
    float      error=0, error2=0, work[1];
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproductFloatComplex c_neg_one = MAGMA_minproduct_C_NEG_ONE;
    magma_minproduct_int_t n_local[Magma_minproductMaxGPUs];

    magma_minproduct_int_t N, Noffset, lda, ldda, blocks, lhwork, ldwork, matsize, vecsize;
    magma_minproduct_int_t incx = 1;

    magma_minproductFloatComplex alpha = MAGMA_minproduct_C_MAKE(  1.5, -2.3 );
    magma_minproductFloatComplex beta  = MAGMA_minproduct_C_MAKE( -0.6,  0.8 );
    magma_minproductFloatComplex *A, *X, *Y, *Ylapack, *Ycublas, *Ymagma_minproduct, *Ymagma_minproduct1, *hwork;
    magma_minproductFloatComplex_ptr dA, dX, dY;
    magma_minproductFloatComplex_ptr d_lA[Magma_minproductMaxGPUs], dwork[Magma_minproductMaxGPUs];

    magma_minproduct_device_t dev;
    magma_minproduct_queue_t queues[Magma_minproductMaxGPUs];
    magma_minproduct_int_t     status = 0;
    
    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");

    magma_minproduct_int_t nb = 64;  // required by magma_minproductblas_chemv_mgpu implementation

    for( dev=0; dev < opts.ngpu; ++dev ) {
        magma_minproduct_setdevice( dev );
        magma_minproduct_queue_create( &queues[dev] );
    }
    
    // currently, tests all offsets in the offsets array;
    // comment out loop below to test a specific offset.
    magma_minproduct_int_t offset = opts.offset;
    magma_minproduct_int_t offsets[] = { 0, 1, 31, 32, 33, 63, 64, 65, 100, 200 };
    magma_minproduct_int_t noffsets = sizeof(offsets) / sizeof(*offsets);
    
    printf("uplo = %s, ngpu %d, block size = %d, offset %d\n",
            lapack_uplo_const(opts.uplo), (int) opts.ngpu, (int) nb, (int) offset );
    printf( "                  BLAS                CUBLAS              MAGMA_minproduct 1 GPU         MAGMA_minproduct MGPU       Error rel  Error rel\n"
            "    N  offset     Gflop/s (msec)      Gflop/s (msec)      Gflop/s (msec)      Gflop/s (msec)   to CUBLAS  to LAPACK\n"
            "===================================================================================================================\n" );
    for( int itest = 0; itest < opts.ntest; ++itest ) {
      
      // comment out these two lines & end of loop to test a specific offset
      for( int ioffset=0; ioffset < noffsets; ioffset += 1 ) {
        offset = offsets[ioffset];
        
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N       = opts.nsize[itest];
            Noffset = N + offset;
            lda     = Noffset;
            ldda    = ((Noffset+31)/32)*32;
            matsize = Noffset*ldda;
            vecsize = (Noffset-1)*incx + 1;
            gflops  = FLOPS_CHEMV( N ) / 1e9;
            
            blocks = (N + (offset % nb) - 1)/nb + 1;
            lhwork = N*opts.ngpu;
            ldwork = ldda*(blocks + 1);

            TESTING_MALLOC_CPU( A,       magma_minproductFloatComplex, matsize );
            TESTING_MALLOC_CPU( Y,       magma_minproductFloatComplex, vecsize );
            TESTING_MALLOC_CPU( Ycublas, magma_minproductFloatComplex, vecsize );
            TESTING_MALLOC_CPU( Ymagma_minproduct,  magma_minproductFloatComplex, vecsize );
            TESTING_MALLOC_CPU( Ymagma_minproduct1, magma_minproductFloatComplex, vecsize );
            TESTING_MALLOC_CPU( Ylapack, magma_minproductFloatComplex, vecsize );

            TESTING_MALLOC_PIN( X,       magma_minproductFloatComplex, vecsize );
            TESTING_MALLOC_PIN( hwork,   magma_minproductFloatComplex, lhwork  );
            
            magma_minproduct_setdevice( opts.device );
            TESTING_MALLOC_DEV( dA, magma_minproductFloatComplex, matsize );
            TESTING_MALLOC_DEV( dX, magma_minproductFloatComplex, vecsize );
            TESTING_MALLOC_DEV( dY, magma_minproductFloatComplex, vecsize );
            
            // TODO make magma_minproduct_cmalloc_bcyclic helper function?
            for( dev=0; dev < opts.ngpu; dev++ ) {
                n_local[dev] = ((Noffset/nb)/opts.ngpu)*nb;
                if (dev < (Noffset/nb) % opts.ngpu)
                    n_local[dev] += nb;
                else if (dev == (Noffset/nb) % opts.ngpu)
                    n_local[dev] += Noffset % nb;
                
                magma_minproduct_setdevice( dev );
                TESTING_MALLOC_DEV( d_lA[dev],  magma_minproductFloatComplex, ldda*n_local[dev] );
                TESTING_MALLOC_DEV( dwork[dev], magma_minproductFloatComplex, ldwork );
            }
            
            //////////////////////////////////////////////////////////////////////////
            
            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &matsize, A );
            magma_minproduct_cmake_hermitian( Noffset, A, lda );
            
            lapackf77_clarnv( &ione, ISEED, &vecsize, X );
            lapackf77_clarnv( &ione, ISEED, &vecsize, Y );
            
            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            magma_minproduct_setdevice( opts.device );
            magma_minproduct_csetmatrix( Noffset, Noffset, A, lda, dA, ldda );
            magma_minproduct_csetvector( Noffset, X, incx, dX, incx );
            magma_minproduct_csetvector( Noffset, Y, incx, dY, incx );
            
            cuda_time = magma_minproduct_sync_wtime(0);
            cublasChemv( opts.handle, cublas_uplo_const(opts.uplo), N,
                         &alpha, dA + offset + offset*ldda, ldda,
                                 dX + offset, incx,
                         &beta,  dY + offset, incx );
            cuda_time = magma_minproduct_sync_wtime(0) - cuda_time;
            cuda_perf = gflops / cuda_time;
            
            magma_minproduct_cgetvector( Noffset, dY, incx, Ycublas, incx );
            
            /* =====================================================================
               Performs operation using MAGMA_minproductBLAS (1 GPU)
               =================================================================== */
            magma_minproduct_setdevice( opts.device );
            magma_minproduct_csetvector( Noffset, Y, incx, dY, incx );
            
            gpu_time = magma_minproduct_sync_wtime( opts.queue );
            
            magma_minproductblas_chemv_work( opts.uplo, N,
                                  alpha, dA + offset + offset*ldda, ldda,
                                         dX + offset, incx,
                                  beta,  dY + offset, incx, dwork[ opts.device ], ldwork,
                                  opts.queue );
            
            gpu_time = magma_minproduct_sync_wtime( opts.queue ) - gpu_time;
            gpu_perf = gflops / gpu_time;
            magma_minproduct_cgetvector( Noffset, dY, incx, Ymagma_minproduct1, incx );
            
            /* =====================================================================
               Performs operation using MAGMA_minproductBLAS (multi-GPU)
               =================================================================== */
            magma_minproduct_csetmatrix_1D_col_bcyclic( Noffset, Noffset, A, lda, d_lA, ldda, opts.ngpu, nb );
            blasf77_ccopy( &Noffset, Y, &incx, Ymagma_minproduct, &incx );
            
            // workspaces do NOT need to be zero -- set to NAN to prove
            for( dev=0; dev < opts.ngpu; ++dev ) {
                magma_minproduct_setdevice( dev );
                magma_minproductblas_claset( Magma_minproductFull, ldwork, 1, MAGMA_minproduct_C_NAN, MAGMA_minproduct_C_NAN, dwork[dev], ldwork );
            }
            lapackf77_claset( "Full", &lhwork, &ione, &MAGMA_minproduct_C_NAN, &MAGMA_minproduct_C_NAN, hwork, &lhwork );
            
            mgpu_time = magma_minproduct_sync_wtime(0);
            
            magma_minproduct_int_t info;
            info = magma_minproductblas_chemv_mgpu(
                opts.uplo, N,
                alpha,
                d_lA, ldda, offset,
                X + offset, incx,
                beta,
                Ymagma_minproduct + offset, incx,
                hwork, lhwork,
                dwork, ldwork,
                opts.ngpu, nb, queues );
            if ( info != 0 )
                printf("magma_minproductblas_chemv_mgpu returned error %d: %s.\n",
                       (int) info, magma_minproduct_strerror( info ));
            
            info = magma_minproductblas_chemv_mgpu_sync(
                opts.uplo, N,
                alpha,
                d_lA, ldda, offset,
                X + offset, incx,
                beta,
                Ymagma_minproduct + offset, incx,
                hwork, lhwork,
                dwork, ldwork,
                opts.ngpu, nb, queues );
            if ( info != 0 )
                printf("magma_minproductblas_chemv_sync returned error %d: %s.\n",
                       (int) info, magma_minproduct_strerror( info ));
            
            mgpu_time = magma_minproduct_sync_wtime(0) - mgpu_time;
            mgpu_perf = gflops / mgpu_time;
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                blasf77_ccopy( &Noffset, Y, &incx, Ylapack, &incx );
                
                cpu_time = magma_minproduct_wtime();
                blasf77_chemv( lapack_uplo_const(opts.uplo), &N,
                               &alpha, A + offset + offset*lda, &lda,
                                       X + offset, &incx,
                               &beta,  Ylapack + offset, &incx );
                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
    
                /* =====================================================================
                   Compute the Difference LAPACK vs. Magma_minproduct
                   =================================================================== */
                error2 = lapackf77_clange( "F", &Noffset, &ione, Ylapack, &Noffset, work );
                blasf77_caxpy( &Noffset, &c_neg_one, Ymagma_minproduct, &incx, Ylapack, &incx );
                error2 = lapackf77_clange( "F", &Noffset, &ione, Ylapack, &Noffset, work ) / error2;
            }
            
            /* =====================================================================
               Compute the Difference Cublas vs. Magma_minproduct
               =================================================================== */            
            error = lapackf77_clange( "F", &Noffset, &ione, Ycublas, &Noffset, work );
            blasf77_caxpy( &Noffset, &c_neg_one, Ymagma_minproduct, &incx, Ycublas, &incx );
            error = lapackf77_clange( "F", &Noffset, &ione, Ycublas, &Noffset, work ) / error;
            
            bool okay = (error < tol && error2 < tol);
            status += ! okay;
            if ( opts.lapack ) {
                printf( "%5d  %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %8.2e   %s\n",
                        (int) N, (int) offset,
                         cpu_perf,  cpu_time*1000.,
                        cuda_perf, cuda_time*1000.,
                         gpu_perf,  gpu_time*1000.,
                        mgpu_perf, mgpu_time*1000.,
                        error, error2, (okay ? "ok" : "failed") );
            }
            else {
                printf( "%5d  %5d     ---   (  ---  )   %7.2f (%7.2f)   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e     ---      %s\n",
                        (int) N, (int) offset,
                        cuda_perf, cuda_time*1000.,
                         gpu_perf,  gpu_time*1000.,
                        mgpu_perf, mgpu_time*1000.,
                        error, (okay ? "ok" : "failed") );
            }
            
            /* Free Memory */
            TESTING_FREE_CPU( A );
            TESTING_FREE_CPU( Y );
            TESTING_FREE_CPU( Ycublas );
            TESTING_FREE_CPU( Ymagma_minproduct  );
            TESTING_FREE_CPU( Ymagma_minproduct1 );
            TESTING_FREE_CPU( Ylapack );

            TESTING_FREE_PIN( X );
            TESTING_FREE_PIN( hwork   );
            
            magma_minproduct_setdevice( opts.device );
            TESTING_FREE_DEV( dA );
            TESTING_FREE_DEV( dX );
            TESTING_FREE_DEV( dY );
            
            for( dev=0; dev < opts.ngpu; dev++ ) {
                magma_minproduct_setdevice( dev );
                TESTING_FREE_DEV( d_lA[dev]  );
                TESTING_FREE_DEV( dwork[dev] );
            }
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
        
      // comment out these two lines line & top of loop test a specific offset
      }  // end for ioffset
      printf( "\n" );
      
    }
    
    for( dev=0; dev < opts.ngpu; ++dev ) {
        magma_minproduct_setdevice( dev );
        magma_minproduct_queue_destroy( queues[dev] );
    }
    
    TESTING_FINALIZE();
    return status;
}
