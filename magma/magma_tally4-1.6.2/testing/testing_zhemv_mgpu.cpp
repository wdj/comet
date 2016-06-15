/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> c d s
       
       @author Mark Gates
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "testings.h"  // before magma_tally4.h, to include cublas_v2
#include "flops.h"
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"

#include "magma_tally4_operators.h"

#define PRECISION_z


// --------------------
int main(int argc, char **argv)
{
    TESTING_INIT();

    real_Double_t gflops, cpu_time=0, cpu_perf=0, gpu_time, gpu_perf, mgpu_time, mgpu_perf, cuda_time, cuda_perf;
    double      error=0, error2=0, work[1];
    magma_tally4_int_t ione     = 1;
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    magma_tally4DoubleComplex c_neg_one = MAGMA_tally4_Z_NEG_ONE;
    magma_tally4_int_t n_local[Magma_tally4MaxGPUs];

    magma_tally4_int_t N, Noffset, lda, ldda, blocks, lhwork, ldwork, matsize, vecsize;
    magma_tally4_int_t incx = 1;

    magma_tally4DoubleComplex alpha = MAGMA_tally4_Z_MAKE(  1.5, -2.3 );
    magma_tally4DoubleComplex beta  = MAGMA_tally4_Z_MAKE( -0.6,  0.8 );
    magma_tally4DoubleComplex *A, *X, *Y, *Ylapack, *Ycublas, *Ymagma_tally4, *Ymagma_tally41, *hwork;
    magma_tally4DoubleComplex_ptr dA, dX, dY;
    magma_tally4DoubleComplex_ptr d_lA[Magma_tally4MaxGPUs], dwork[Magma_tally4MaxGPUs];

    magma_tally4_device_t dev;
    magma_tally4_queue_t queues[Magma_tally4MaxGPUs];
    magma_tally4_int_t     status = 0;
    
    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );
    
    double tol = opts.tolerance * lapackf77_dlamch("E");

    magma_tally4_int_t nb = 64;  // required by magma_tally4blas_zhemv_mgpu implementation

    for( dev=0; dev < opts.ngpu; ++dev ) {
        magma_tally4_setdevice( dev );
        magma_tally4_queue_create( &queues[dev] );
    }
    
    // currently, tests all offsets in the offsets array;
    // comment out loop below to test a specific offset.
    magma_tally4_int_t offset = opts.offset;
    magma_tally4_int_t offsets[] = { 0, 1, 31, 32, 33, 63, 64, 65, 100, 200 };
    magma_tally4_int_t noffsets = sizeof(offsets) / sizeof(*offsets);
    
    printf("uplo = %s, ngpu %d, block size = %d, offset %d\n",
            lapack_uplo_const_tally4(opts.uplo), (int) opts.ngpu, (int) nb, (int) offset );
    printf( "                  BLAS                CUBLAS              MAGMA_tally4 1 GPU         MAGMA_tally4 MGPU       Error rel  Error rel\n"
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
            gflops  = FLOPS_ZHEMV( N ) / 1e9;
            
            blocks = (N + (offset % nb) - 1)/nb + 1;
            lhwork = N*opts.ngpu;
            ldwork = ldda*(blocks + 1);

            TESTING_MALLOC_CPU( A,       magma_tally4DoubleComplex, matsize );
            TESTING_MALLOC_CPU( Y,       magma_tally4DoubleComplex, vecsize );
            TESTING_MALLOC_CPU( Ycublas, magma_tally4DoubleComplex, vecsize );
            TESTING_MALLOC_CPU( Ymagma_tally4,  magma_tally4DoubleComplex, vecsize );
            TESTING_MALLOC_CPU( Ymagma_tally41, magma_tally4DoubleComplex, vecsize );
            TESTING_MALLOC_CPU( Ylapack, magma_tally4DoubleComplex, vecsize );

            TESTING_MALLOC_PIN( X,       magma_tally4DoubleComplex, vecsize );
            TESTING_MALLOC_PIN( hwork,   magma_tally4DoubleComplex, lhwork  );
            
            magma_tally4_setdevice( opts.device );
            TESTING_MALLOC_DEV( dA, magma_tally4DoubleComplex, matsize );
            TESTING_MALLOC_DEV( dX, magma_tally4DoubleComplex, vecsize );
            TESTING_MALLOC_DEV( dY, magma_tally4DoubleComplex, vecsize );
            
            // TODO make magma_tally4_zmalloc_bcyclic helper function?
            for( dev=0; dev < opts.ngpu; dev++ ) {
                n_local[dev] = ((Noffset/nb)/opts.ngpu)*nb;
                if (dev < (Noffset/nb) % opts.ngpu)
                    n_local[dev] += nb;
                else if (dev == (Noffset/nb) % opts.ngpu)
                    n_local[dev] += Noffset % nb;
                
                magma_tally4_setdevice( dev );
                TESTING_MALLOC_DEV( d_lA[dev],  magma_tally4DoubleComplex, ldda*n_local[dev] );
                TESTING_MALLOC_DEV( dwork[dev], magma_tally4DoubleComplex, ldwork );
            }
            
            //////////////////////////////////////////////////////////////////////////
            
            /* Initialize the matrix */
            lapackf77_zlarnv( &ione, ISEED, &matsize, A );
            magma_tally4_zmake_hermitian( Noffset, A, lda );
            
            lapackf77_zlarnv( &ione, ISEED, &vecsize, X );
            lapackf77_zlarnv( &ione, ISEED, &vecsize, Y );
            
            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            magma_tally4_setdevice( opts.device );
            magma_tally4_zsetmatrix( Noffset, Noffset, A, lda, dA, ldda );
            magma_tally4_zsetvector( Noffset, X, incx, dX, incx );
            magma_tally4_zsetvector( Noffset, Y, incx, dY, incx );
            
            cuda_time = magma_tally4_sync_wtime(0);
            cublasZhemv( opts.handle, cublas_uplo_const_tally4(opts.uplo), N,
                         &alpha, dA + offset + offset*ldda, ldda,
                                 dX + offset, incx,
                         &beta,  dY + offset, incx );
            cuda_time = magma_tally4_sync_wtime(0) - cuda_time;
            cuda_perf = gflops / cuda_time;
            
            magma_tally4_zgetvector( Noffset, dY, incx, Ycublas, incx );
            
            /* =====================================================================
               Performs operation using MAGMA_tally4BLAS (1 GPU)
               =================================================================== */
            magma_tally4_setdevice( opts.device );
            magma_tally4_zsetvector( Noffset, Y, incx, dY, incx );
            
            gpu_time = magma_tally4_sync_wtime( opts.queue );
            
            magma_tally4blas_zhemv_work( opts.uplo, N,
                                  alpha, dA + offset + offset*ldda, ldda,
                                         dX + offset, incx,
                                  beta,  dY + offset, incx, dwork[ opts.device ], ldwork,
                                  opts.queue );
            
            gpu_time = magma_tally4_sync_wtime( opts.queue ) - gpu_time;
            gpu_perf = gflops / gpu_time;
            magma_tally4_zgetvector( Noffset, dY, incx, Ymagma_tally41, incx );
            
            /* =====================================================================
               Performs operation using MAGMA_tally4BLAS (multi-GPU)
               =================================================================== */
            magma_tally4_zsetmatrix_1D_col_bcyclic( Noffset, Noffset, A, lda, d_lA, ldda, opts.ngpu, nb );
            blasf77_zcopy( &Noffset, Y, &incx, Ymagma_tally4, &incx );
            
            // workspaces do NOT need to be zero -- set to NAN to prove
            for( dev=0; dev < opts.ngpu; ++dev ) {
                magma_tally4_setdevice( dev );
                magma_tally4blas_zlaset( Magma_tally4Full, ldwork, 1, MAGMA_tally4_Z_NAN, MAGMA_tally4_Z_NAN, dwork[dev], ldwork );
            }
            lapackf77_zlaset( "Full", &lhwork, &ione, &MAGMA_tally4_Z_NAN, &MAGMA_tally4_Z_NAN, hwork, &lhwork );
            
            mgpu_time = magma_tally4_sync_wtime(0);
            
            magma_tally4_int_t info;
            info = magma_tally4blas_zhemv_mgpu(
                opts.uplo, N,
                alpha,
                d_lA, ldda, offset,
                X + offset, incx,
                beta,
                Ymagma_tally4 + offset, incx,
                hwork, lhwork,
                dwork, ldwork,
                opts.ngpu, nb, queues );
            if ( info != 0 )
                printf("magma_tally4blas_zhemv_mgpu returned error %d: %s.\n",
                       (int) info, magma_tally4_strerror( info ));
            
            info = magma_tally4blas_zhemv_mgpu_sync(
                opts.uplo, N,
                alpha,
                d_lA, ldda, offset,
                X + offset, incx,
                beta,
                Ymagma_tally4 + offset, incx,
                hwork, lhwork,
                dwork, ldwork,
                opts.ngpu, nb, queues );
            if ( info != 0 )
                printf("magma_tally4blas_zhemv_sync returned error %d: %s.\n",
                       (int) info, magma_tally4_strerror( info ));
            
            mgpu_time = magma_tally4_sync_wtime(0) - mgpu_time;
            mgpu_perf = gflops / mgpu_time;
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                blasf77_zcopy( &Noffset, Y, &incx, Ylapack, &incx );
                
                cpu_time = magma_tally4_wtime();
                blasf77_zhemv( lapack_uplo_const_tally4(opts.uplo), &N,
                               &alpha, A + offset + offset*lda, &lda,
                                       X + offset, &incx,
                               &beta,  Ylapack + offset, &incx );
                cpu_time = magma_tally4_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
    
                /* =====================================================================
                   Compute the Difference LAPACK vs. Magma_tally4
                   =================================================================== */
                error2 = lapackf77_zlange( "F", &Noffset, &ione, Ylapack, &Noffset, work );
                blasf77_zaxpy( &Noffset, &c_neg_one, Ymagma_tally4, &incx, Ylapack, &incx );
                error2 = lapackf77_zlange( "F", &Noffset, &ione, Ylapack, &Noffset, work ) / error2;
            }
            
            /* =====================================================================
               Compute the Difference Cublas vs. Magma_tally4
               =================================================================== */            
            error = lapackf77_zlange( "F", &Noffset, &ione, Ycublas, &Noffset, work );
            blasf77_zaxpy( &Noffset, &c_neg_one, Ymagma_tally4, &incx, Ycublas, &incx );
            error = lapackf77_zlange( "F", &Noffset, &ione, Ycublas, &Noffset, work ) / error;
            
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
            TESTING_FREE_CPU( Ymagma_tally4  );
            TESTING_FREE_CPU( Ymagma_tally41 );
            TESTING_FREE_CPU( Ylapack );

            TESTING_FREE_PIN( X );
            TESTING_FREE_PIN( hwork   );
            
            magma_tally4_setdevice( opts.device );
            TESTING_FREE_DEV( dA );
            TESTING_FREE_DEV( dX );
            TESTING_FREE_DEV( dY );
            
            for( dev=0; dev < opts.ngpu; dev++ ) {
                magma_tally4_setdevice( dev );
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
        magma_tally4_setdevice( dev );
        magma_tally4_queue_destroy( queues[dev] );
    }
    
    TESTING_FINALIZE();
    return status;
}
