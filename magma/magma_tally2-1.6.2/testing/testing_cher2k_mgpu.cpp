/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zher2k_mgpu.cpp normal z -> c, Fri Jan 30 19:00:24 2015
       
       @author Mark Gates
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// includes, project
#include "flops.h"
#include "magma_tally2.h"
#include "magma_tally2_lapack.h"
#include "testings.h"

// define ICHI to test with Ichi's version, too
#undef ICHI


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing magma_tally2_cher2k_mgpu
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    magma_tally2FloatComplex c_neg_one = MAGMA_tally2_C_NEG_ONE;
    magma_tally2FloatComplex alpha = MAGMA_tally2_C_MAKE( 1.2345, 4.321 );
    float beta = 3.14159;
    
    real_Double_t    gflops, gpu_perf, cpu_perf, gpu_time, cpu_time;
    float           error, work[1];
    magma_tally2FloatComplex *hA, *hR, *hR2, *hV, *hW;
    magma_tally2FloatComplex_ptr dV[Magma_tally2MaxGPUs], dW[Magma_tally2MaxGPUs], dA[Magma_tally2MaxGPUs];
    magma_tally2_int_t n, k, size, lda, ldda, nb, ngpu, nstream;
    magma_tally2_int_t ione     = 1;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};

    magma_tally2_queue_t streams[Magma_tally2MaxGPUs][20];
    magma_tally2_int_t status = 0;
    
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );

    float tol = opts.tolerance * lapackf77_slamch("E");
    
    ngpu    = opts.ngpu;
    nb      = (opts.nb      > 0 ? opts.nb      : 64);
    nstream = (opts.nstream > 0 ? opts.nstream :  2);
    
    printf( "version 1: magma_tally2blas_cher2k_mgpu2     %s\n", (opts.version==1 ? "(enabled)" : ""));
    printf( "version 2: magma_tally2blas_cher2k_mgpu_spec %s\n", (opts.version==2 ? "(enabled)" : ""));
#ifdef ICHI
    printf( "version 3: magma_tally2_cher2k_mgpu (Ichi)   %s\n", (opts.version==3 ? "(enabled)" : ""));
#endif
    printf( "\n" );
    
    printf( "nb %d, ngpu %d, nstream %d\n", (int) nb, (int) ngpu, (int) nstream );
    printf("    n     k    nb offset  CPU GFlop/s (sec)   GPU GFlop/s (sec)   |R|/(|V|*|W|+|A|)\n");
    printf("===================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        n = opts.nsize[itest];
        k = opts.ksize[itest];
        
        for( int offset = 0; offset < n; offset += min(k,nb) ) {
            for( int iter = 0; iter < opts.niter; ++iter ) {
                lda    = n;
                ldda   = ((n + 31)/32)*32;
                gflops = FLOPS_CHER2K( k, n-offset ) / 1e9;
                
                TESTING_MALLOC_CPU( hA,  magma_tally2FloatComplex, lda*n   );
                TESTING_MALLOC_CPU( hR,  magma_tally2FloatComplex, lda*n   );
                TESTING_MALLOC_CPU( hR2, magma_tally2FloatComplex, lda*n   );
                TESTING_MALLOC_CPU( hV,  magma_tally2FloatComplex, lda*k*2 );
                //TESTING_MALLOC_CPU( hW,  magma_tally2FloatComplex, lda*k   );
                for( int d = 0; d < ngpu; ++d ) {
                    magma_tally2_int_t nlocal = ((n / nb) / ngpu + 1) * nb;
                    magma_tally2_setdevice( d );
                    TESTING_MALLOC_DEV( dA[d], magma_tally2FloatComplex, ldda*nlocal );
                    TESTING_MALLOC_DEV( dV[d], magma_tally2FloatComplex, ldda*k*2    );
                    //TESTING_MALLOC_DEV( dW[d], magma_tally2FloatComplex, ldda*k      );
                    for( int i = 0; i < nstream; ++i ) {
                        magma_tally2_queue_create( &streams[d][i] );
                    }
                }
                
                size = lda*n;
                lapackf77_clarnv( &ione, ISEED, &size, hA );
                size = lda*k*2;
                lapackf77_clarnv( &ione, ISEED, &size, hV );
                hW = hV + lda*k;
                //lapackf77_clarnv( &ione, ISEED, &size, hW );
                
                /* ====================================================================
                   Performs operation using MAGMA_tally2
                   =================================================================== */
                magma_tally2_csetmatrix_1D_col_bcyclic( n, n, hA, lda, dA, ldda, ngpu, nb );
                for( int d = 0; d < ngpu; ++d ) {
                    magma_tally2_setdevice( d );
                    dW[d] = dV[d] + ldda*k;
                    magma_tally2_csetmatrix( n, k, hV, lda, dV[d], ldda );
                    magma_tally2_csetmatrix( n, k, hW, lda, dW[d], ldda );
                }
                
                gpu_time = magma_tally2_sync_wtime(0);
                
                if ( opts.version == 1 ) {
                    magma_tally2blas_cher2k_mgpu2(
                        Magma_tally2Lower, Magma_tally2NoTrans, n-offset, k,
                        alpha, dV, ldda, 0,
                               dW, ldda, 0,
                        beta,  dA, ldda, offset,
                        ngpu, nb, streams, nstream );
                }
                else if ( opts.version == 2 ) {
                    magma_tally2blas_cher2k_mgpu_spec(
                        Magma_tally2Lower, Magma_tally2NoTrans, n-offset, k,
                        alpha, dV, ldda, 0,
                               dW, ldda, 0,
                        beta,  dA, ldda, offset,
                        ngpu, nb, streams, nstream );
                }
                else {
#ifdef ICHI
                    magma_tally2_cher2k_mgpu(
                        ngpu, Magma_tally2Lower, Magma_tally2NoTrans, nb, n-offset, k,
                        alpha, dV, ldda,
                               //dW, ldda,
                        beta,  dA, ldda, offset,
                        nstream, streams );
#endif
                }
                
                gpu_time = magma_tally2_sync_wtime(0) - gpu_time;
                gpu_perf = gflops / gpu_time;
                
                // Get dA back to the CPU to compare with the CPU result.
                magma_tally2_cgetmatrix_1D_col_bcyclic( n, n, dA, ldda, hR, lda, ngpu, nb );
                
                /* =====================================================================
                   Performs operation using LAPACK
                   =================================================================== */
                if ( opts.lapack || opts.check ) {
                    // store ||V||*||W|| + ||A||
                    magma_tally2_int_t n_offset = n - offset;
                    error  = lapackf77_clange("f", &n_offset, &k, hV, &lda, work );
                    error *= lapackf77_clange("f", &n_offset, &k, hW, &lda, work );
                    error += lapackf77_clange("f", &n_offset, &n_offset, &hA[offset + offset*lda], &lda, work );
                    
                    cpu_time = magma_tally2_wtime();
                    blasf77_cher2k( "Lower", "NoTrans", &n_offset, &k,
                                    &alpha, hV, &lda,
                                            hW, &lda,
                                    &beta,  &hA[offset + offset*lda], &lda );
                    cpu_time = magma_tally2_wtime() - cpu_time;
                    cpu_perf = gflops / cpu_time;
                    
                    // compute relative error ||R||/||A||, where R := A_magma_tally2 - A_lapack = R - A
                    size = lda*n;
                    blasf77_caxpy( &size, &c_neg_one, hA, &ione, hR, &ione );
                    error = lapackf77_clanhe("fro", "Lower", &n_offset, &hR[offset + offset*lda], &lda, work) / error;
                    
                    printf( "%5d %5d %5d %5d   %7.1f (%7.4f)   %7.1f (%7.4f)   %8.2e   %s\n",
                            (int) n, (int) k, (int) nb, (int) offset,
                            cpu_perf, cpu_time, gpu_perf, gpu_time,
                            error, (error < tol ? "ok" : "failed"));
                            //, gpu_perf2, gpu_time2, error, error2 );
                    status += ! (error < tol);
                }
                else {
                    printf( "%5d %5d %5d %5d     ---   (  ---  )   %7.1f (%7.4f)     ---\n",
                            (int) n, (int) k, (int) nb, (int) offset,
                            gpu_perf, gpu_time );
                }
                
                TESTING_FREE_CPU( hA  );
                TESTING_FREE_CPU( hR  );
                TESTING_FREE_CPU( hR2 );
                TESTING_FREE_CPU( hV  );
                //TESTING_FREE_CPU( hW );
                for( int d = 0; d < ngpu; ++d ) {
                    magma_tally2_setdevice( d );
                    TESTING_FREE_DEV( dA[d] );
                    TESTING_FREE_DEV( dV[d] );
                    //TESTING_FREE_DEV( dW[d] );
                    for( int i = 0; i < nstream; ++i ) {
                        magma_tally2_queue_destroy( streams[d][i] );
                    }
                }
                fflush( stdout );
            }
            if ( opts.niter > 1 ) {
                printf( "\n" );
            }
        } // offset
        printf( "\n" );
    }
    
    TESTING_FINALIZE();
    return status;
}