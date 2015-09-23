/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
       
       @author Mark Gates
       @author Azzam Haidar
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// includes, project
#include "flops.h"
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "testings.h"

#include "trace.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing magma_minproduct_zhemm_mgpu
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    magma_minproductDoubleComplex c_neg_one = MAGMA_minproduct_Z_NEG_ONE;
    magma_minproductDoubleComplex calpha    = MAGMA_minproduct_Z_MAKE( 3.456, 5.678 );
    magma_minproductDoubleComplex cbeta     = MAGMA_minproduct_Z_MAKE( 1.234, 2.456 );
    
    real_Double_t    gflops, gpu_perf=0., cpu_perf=0., gpu_time=0., cpu_time=0.;
    real_Double_t    gpu_perf2=0., gpu_time2=0.;
    double           error=0., errorbis=0., work[1];
    magma_minproductDoubleComplex *hA, *hX, *hB, *hR;
    magma_minproductDoubleComplex_ptr dA[Magma_minproductMaxGPUs], dX[Magma_minproductMaxGPUs], dB[Magma_minproductMaxGPUs], dwork[Magma_minproductMaxGPUs], hwork[Magma_minproductMaxGPUs+1];
    magma_minproductDoubleComplex_ptr dA2;
    magma_minproduct_int_t M, N, size, lda, ldda, msize, nb, nstream;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t iseed[4] = {0,0,0,1};
    magma_minproduct_int_t status = 0;
    
    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    // default values
    nb      = (opts.nb      > 0 ? opts.nb      : 64);
    nstream = (opts.nstream > 0 ? opts.nstream :  2);
    
    magma_minproduct_int_t gnode[Magma_minproductMaxGPUs][Magma_minproductMaxGPUs+2];
    magma_minproduct_int_t nbcmplx = 0;
    magma_minproduct_buildconnection_mgpu(gnode, &nbcmplx, opts.ngpu);
    printf("Initializing communication pattern... GPU-ncmplx %d\n\n", (int) nbcmplx);

    for (int i=0; i < nbcmplx; ++i) {
        int myngpu = gnode[i][Magma_minproductMaxGPUs];
        printf("cmplx %d has %d gpu ", i, myngpu);
        for(int j=0; j < myngpu; ++j)
            printf("  %d", (int) gnode[i][j]);
        printf("\n");
    }

    magma_minproduct_int_t nbevents = 2;
    magma_minproduct_queue_t streams[Magma_minproductMaxGPUs][20];
    magma_minproduct_event_t redevents[Magma_minproductMaxGPUs][20];
    magma_minproduct_event_t redevents2[Magma_minproductMaxGPUs][Magma_minproductMaxGPUs*Magma_minproductMaxGPUs+10];
    for( int d = 0; d < opts.ngpu; ++d ) {
        for( magma_minproduct_int_t i = 0; i < nstream; ++i ) {
            magma_minproduct_queue_create( &streams[d][i] );
        }
        for( magma_minproduct_int_t i = 0; i < nbevents; ++i ) {
            cudaEventCreateWithFlags(&redevents[d][i],  cudaEventDisableTiming);
            cudaEventCreateWithFlags(&redevents2[d][i], cudaEventDisableTiming);
        }
    }

    printf( "nb %d, ngpu %d, nstream %d version %d\n", (int) nb, (int) opts.ngpu, (int) nstream, (int) opts.version );
    printf("    M     N    nb offset  CPU GFlop/s (sec)   GPU GFlop/s (sec)   CUBLAS hemm (sec)   ||R|| / ||A||*||X||\n");
    printf("=========================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
      M = opts.msize[itest];
      N = opts.nsize[itest];
      for( int offset = 0; offset < N; offset += min(N,nb) ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            msize = M - offset;
            lda   = M;
            ldda  = ((M + 31)/32)*32;
            size  = lda*M;
            gflops = FLOPS_ZHEMM( Magma_minproductLeft, (double)msize, (double)N ) / 1e9;
            
            magma_minproduct_int_t dworksiz = ldda*N*3;
            magma_minproduct_int_t hworksiz = lda*N;
            
            TESTING_MALLOC_CPU( hA, magma_minproductDoubleComplex, lda*M );
            TESTING_MALLOC_CPU( hX, magma_minproductDoubleComplex, lda*N );
            TESTING_MALLOC_CPU( hB, magma_minproductDoubleComplex, lda*N );
            
            TESTING_MALLOC_PIN( hR, magma_minproductDoubleComplex, lda*N );

            for( int d = 0; d < opts.ngpu; ++d ) {
                magma_minproduct_int_t mlocal = ((M / nb) / opts.ngpu + 1) * nb;
                magma_minproduct_setdevice( d );
                TESTING_MALLOC_DEV( dA[d],    magma_minproductDoubleComplex, ldda*mlocal );
                TESTING_MALLOC_DEV( dX[d],    magma_minproductDoubleComplex, ldda*N      );
                TESTING_MALLOC_DEV( dB[d],    magma_minproductDoubleComplex, ldda*N      );
                TESTING_MALLOC_DEV( dwork[d], magma_minproductDoubleComplex, dworksiz    );
                
                TESTING_MALLOC_PIN( hwork[d], magma_minproductDoubleComplex, hworksiz    );
            }
            TESTING_MALLOC_PIN( hwork[opts.ngpu], magma_minproductDoubleComplex, lda*N );
        
            if ( opts.check ) {
                magma_minproduct_setdevice( 0 );
                TESTING_MALLOC_DEV( dA2, magma_minproductDoubleComplex, ldda*M );
            }

            lapackf77_zlarnv( &ione, iseed, &size, hA );
            magma_minproduct_zmake_hermitian( M, hA, lda );
            
            size = lda*N;
            lapackf77_zlarnv( &ione, iseed, &size, hX );
            lapackf77_zlarnv( &ione, iseed, &size, hB );
            lapackf77_zlacpy( "Full", &M, &N, hB, &lda, hR, &lda );
            
            /* ====================================================================
               Performs operation using MAGMA_minproduct
               =================================================================== */
            magma_minproduct_zsetmatrix_1D_col_bcyclic( M, M, hA, lda, dA, ldda, opts.ngpu, nb );
            for( int d = 0; d < opts.ngpu; ++d ) {
                magma_minproduct_setdevice( d );
                //magma_minproductblasSetKernelStream( streams[ d ][  0 ] );
                magma_minproduct_zsetmatrix( M, N, hX, lda, dX[d], ldda );
                //if (d == 0) magma_minproduct_zsetmatrix( M, N, hB, lda, dB[d], ldda ); // this is wrong coz when offset != 0 the gpu who do the beta*C may be not 0 so this should be related to stdev(starting device who own i=0 first col)
                magma_minproduct_zsetmatrix( M, N, hB, lda, dB[d], ldda );
            }
        
            //memset(hR, 0, lda*N*sizeof(magma_minproductDoubleComplex));
    
            trace_init( 1, opts.ngpu, nstream, (magma_minproduct_queue_t*) streams );
    
            //magma_minproduct_int_t offset = 0; //nb;
    
            gpu_time = magma_minproduct_sync_wtime(0);
        
            magma_minproductblas_zhemm_mgpu_com(
                Magma_minproductLeft, Magma_minproductLower, msize, N,
                calpha,    dA, ldda, offset,
                           dX, ldda,
                cbeta,     dB, ldda, dwork, dworksiz, hR, lda, hwork, hworksiz,
                opts.ngpu, nb, streams, nstream, redevents2, nbevents, gnode, nbcmplx);
           
            gpu_time = magma_minproduct_sync_wtime(0) - gpu_time;
            gpu_perf = gflops / gpu_time;
                
            #ifdef TRACING
            char buf[80];
            snprintf( buf, sizeof(buf), "zhemm-m%d-n%d-nb%d-stream%d-ngpu%d-run%d.svg",
                      (int) M, (int) N, (int) nb, (int) nstream, (int) opts.ngpu, (int) iter );
            trace_finalize( buf, "trace.css" );
            #endif
            
            /* ====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            if ( opts.check && iter == 0 ) {
                magma_minproduct_setdevice( 0 );
                magma_minproductblasSetKernelStream(  0  );
                magma_minproduct_zsetmatrix( M, M, hA, lda, dA2, ldda );
                magma_minproduct_zsetmatrix( M, N, hX, lda, dX[0], ldda );
                magma_minproduct_zsetmatrix( M, N, hB, lda, dwork[0], ldda );
                
                gpu_time2 = magma_minproduct_sync_wtime(0);
                magma_minproduct_zhemm(
                    Magma_minproductLeft, Magma_minproductLower, msize, N,
                    calpha,    dA2+offset*ldda+offset, ldda,
                               dX[0],    ldda,
                    cbeta,     dwork[0], ldda );
                gpu_time2 = magma_minproduct_sync_wtime(0) - gpu_time2;
                gpu_perf2 = gflops / gpu_time2;
            }
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.check ) {
                // store ||A||*||X||
                errorbis  = lapackf77_zlange("fro", &msize, &msize, hA+offset*lda+offset, &lda, work );
                errorbis *= lapackf77_zlange("fro", &msize, &N, hX, &lda, work );
                
                //printf( "A =" ); magma_minproduct_zprint( M, M, hA, lda );
                //printf( "X =" ); magma_minproduct_zprint( M, N, hX, lda );
                //printf( "B =" ); magma_minproduct_zprint( M, N, hB, lda );
                
                cpu_time = magma_minproduct_wtime();
                blasf77_zhemm( "Left", "Lower", &msize, &N,
                                &calpha, hA+offset*lda+offset, &lda,
                                         hX, &lda,
                                &cbeta,  hB, &lda );
                cpu_time = magma_minproduct_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                /*
                trace_file = fopen("AJETE/C", "w");
                for (int j = 0; j < N; j++)
                    for (int i = 0; i < siz; i++)
                        fprintf(trace_file, "%10d%10d%40.30e\n", i+1, j+1, hB[j*lda+i]);
                fclose(trace_file);
                */
                magma_minproduct_int_t firstprint=0;
                for(magma_minproduct_int_t dev=0; dev < opts.ngpu; ++dev) {
                    magma_minproduct_setdevice( dev );
                    magma_minproduct_zgetmatrix( M, N, dB[dev], ldda, hR, lda );
    
                    // compute relative error ||R||/||A||*||X||, where R := B_magma_minproduct - B_lapack = R - B
                    size = lda*N;
                    blasf77_zaxpy( &size, &c_neg_one, hB, &ione, hR, &ione );
                    error = lapackf77_zlange("fro", &msize, &N, hR, &lda, work) / errorbis;
                    
                    //printf( "R ="  ); magma_minproduct_zprint( M, N, hR, lda );
                    if (firstprint == 0) {
                        printf( "%5d %5d %5d %5d   %7.1f (%7.4f)   %7.1f (%7.4f)   %7.1f (%7.4f)   %8.2e   %s\n",
                                (int) M, (int) N, (int) nb, (int) offset,
                                cpu_perf, cpu_time,
                                gpu_perf, gpu_time,
                                gpu_perf2, gpu_time2,
                                error, (error < tol ? "ok" : "failed") );
                    }
                    else {
                        printf( "%89s  %8.2e   %s\n", " ",
                                error, (error < tol ? "ok" : "failed") );
                    }
                    status += ! (error < tol);
                    firstprint =1;
                }
            } else {
                printf( "%5d %5d %5d %5d     ---   (  ---  )   %7.1f (%7.4f)     ---   (  ---  )   ---\n",
                        (int) M, (int) N, (int) nb, (int) offset,
                        gpu_perf, gpu_time );
            }
    
            TESTING_FREE_CPU( hA );
            TESTING_FREE_CPU( hX );
            TESTING_FREE_CPU( hB );
            
            TESTING_FREE_PIN( hR );
        
            for( int d = 0; d < opts.ngpu; ++d ) {
                magma_minproduct_setdevice( d );
                TESTING_FREE_DEV( dA[d]    );
                TESTING_FREE_DEV( dX[d]    );
                TESTING_FREE_DEV( dB[d]    );
                TESTING_FREE_DEV( dwork[d] );
                
                TESTING_FREE_PIN( hwork[d] );
            }
            TESTING_FREE_PIN( hwork[opts.ngpu] );
        
            if ( opts.check ) {
                magma_minproduct_setdevice( 0 );
                TESTING_FREE_DEV( dA2 );
            }
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
      }  // offset
      printf( "\n" );
    }

    for( int d = 0; d < opts.ngpu; ++d ) {
        magma_minproduct_setdevice( d );
        for( magma_minproduct_int_t i = 0; i < nstream; ++i ) {
            magma_minproduct_queue_destroy( streams[d][i] );
        }
        for( magma_minproduct_int_t i = 0; i < nbevents; ++i ) {
            magma_minproduct_event_destroy( redevents[d][i]  );
            magma_minproduct_event_destroy( redevents2[d][i] );
        }
    }
    
    TESTING_FINALIZE();
    return status;
}
