/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
       
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
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"
#include "testings.h"

// define ICHI to test with Ichi's version, too
#undef ICHI


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing magma_tally4_zher2k_mgpu
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    magma_tally4DoubleComplex c_neg_one = MAGMA_tally4_Z_NEG_ONE;
    magma_tally4DoubleComplex alpha = MAGMA_tally4_Z_MAKE( 1.2345, 4.321 );
    double beta = 3.14159;
    
    real_Double_t    gflops, gpu_perf, cpu_perf, gpu_time, cpu_time;
    double           error, work[1];
    magma_tally4DoubleComplex *hA, *hR, *hR2, *hV, *hW;
    magma_tally4DoubleComplex_ptr dV[Magma_tally4MaxGPUs], dW[Magma_tally4MaxGPUs], dA[Magma_tally4MaxGPUs];
    magma_tally4_int_t n, k, size, lda, ldda, nb, ngpu, nstream;
    magma_tally4_int_t ione     = 1;
    magma_tally4_int_t ISEED[4] = {0,0,0,1};

    magma_tally4_queue_t streams[Magma_tally4MaxGPUs][20];
    magma_tally4_int_t status = 0;
    
    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );

    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    ngpu    = opts.ngpu;
    nb      = (opts.nb      > 0 ? opts.nb      : 64);
    nstream = (opts.nstream > 0 ? opts.nstream :  2);
    
    printf( "version 1: magma_tally4blas_zher2k_mgpu2     %s\n", (opts.version==1 ? "(enabled)" : ""));
    printf( "version 2: magma_tally4blas_zher2k_mgpu_spec %s\n", (opts.version==2 ? "(enabled)" : ""));
#ifdef ICHI
    printf( "version 3: magma_tally4_zher2k_mgpu (Ichi)   %s\n", (opts.version==3 ? "(enabled)" : ""));
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
                gflops = FLOPS_ZHER2K( k, n-offset ) / 1e9;
                
                TESTING_MALLOC_CPU( hA,  magma_tally4DoubleComplex, lda*n   );
                TESTING_MALLOC_CPU( hR,  magma_tally4DoubleComplex, lda*n   );
                TESTING_MALLOC_CPU( hR2, magma_tally4DoubleComplex, lda*n   );
                TESTING_MALLOC_CPU( hV,  magma_tally4DoubleComplex, lda*k*2 );
                //TESTING_MALLOC_CPU( hW,  magma_tally4DoubleComplex, lda*k   );
                for( int d = 0; d < ngpu; ++d ) {
                    magma_tally4_int_t nlocal = ((n / nb) / ngpu + 1) * nb;
                    magma_tally4_setdevice( d );
                    TESTING_MALLOC_DEV( dA[d], magma_tally4DoubleComplex, ldda*nlocal );
                    TESTING_MALLOC_DEV( dV[d], magma_tally4DoubleComplex, ldda*k*2    );
                    //TESTING_MALLOC_DEV( dW[d], magma_tally4DoubleComplex, ldda*k      );
                    for( int i = 0; i < nstream; ++i ) {
                        magma_tally4_queue_create( &streams[d][i] );
                    }
                }
                
                size = lda*n;
                lapackf77_zlarnv( &ione, ISEED, &size, hA );
                size = lda*k*2;
                lapackf77_zlarnv( &ione, ISEED, &size, hV );
                hW = hV + lda*k;
                //lapackf77_zlarnv( &ione, ISEED, &size, hW );
                
                /* ====================================================================
                   Performs operation using MAGMA_tally4
                   =================================================================== */
                magma_tally4_zsetmatrix_1D_col_bcyclic( n, n, hA, lda, dA, ldda, ngpu, nb );
                for( int d = 0; d < ngpu; ++d ) {
                    magma_tally4_setdevice( d );
                    dW[d] = dV[d] + ldda*k;
                    magma_tally4_zsetmatrix( n, k, hV, lda, dV[d], ldda );
                    magma_tally4_zsetmatrix( n, k, hW, lda, dW[d], ldda );
                }
                
                gpu_time = magma_tally4_sync_wtime(0);
                
                if ( opts.version == 1 ) {
                    magma_tally4blas_zher2k_mgpu2(
                        Magma_tally4Lower, Magma_tally4NoTrans, n-offset, k,
                        alpha, dV, ldda, 0,
                               dW, ldda, 0,
                        beta,  dA, ldda, offset,
                        ngpu, nb, streams, nstream );
                }
                else if ( opts.version == 2 ) {
                    magma_tally4blas_zher2k_mgpu_spec(
                        Magma_tally4Lower, Magma_tally4NoTrans, n-offset, k,
                        alpha, dV, ldda, 0,
                               dW, ldda, 0,
                        beta,  dA, ldda, offset,
                        ngpu, nb, streams, nstream );
                }
                else {
#ifdef ICHI
                    magma_tally4_zher2k_mgpu(
                        ngpu, Magma_tally4Lower, Magma_tally4NoTrans, nb, n-offset, k,
                        alpha, dV, ldda,
                               //dW, ldda,
                        beta,  dA, ldda, offset,
                        nstream, streams );
#endif
                }
                
                gpu_time = magma_tally4_sync_wtime(0) - gpu_time;
                gpu_perf = gflops / gpu_time;
                
                // Get dA back to the CPU to compare with the CPU result.
                magma_tally4_zgetmatrix_1D_col_bcyclic( n, n, dA, ldda, hR, lda, ngpu, nb );
                
                /* =====================================================================
                   Performs operation using LAPACK
                   =================================================================== */
                if ( opts.lapack || opts.check ) {
                    // store ||V||*||W|| + ||A||
                    magma_tally4_int_t n_offset = n - offset;
                    error  = lapackf77_zlange("f", &n_offset, &k, hV, &lda, work );
                    error *= lapackf77_zlange("f", &n_offset, &k, hW, &lda, work );
                    error += lapackf77_zlange("f", &n_offset, &n_offset, &hA[offset + offset*lda], &lda, work );
                    
                    cpu_time = magma_tally4_wtime();
                    blasf77_zher2k( "Lower", "NoTrans", &n_offset, &k,
                                    &alpha, hV, &lda,
                                            hW, &lda,
                                    &beta,  &hA[offset + offset*lda], &lda );
                    cpu_time = magma_tally4_wtime() - cpu_time;
                    cpu_perf = gflops / cpu_time;
                    
                    // compute relative error ||R||/||A||, where R := A_magma_tally4 - A_lapack = R - A
                    size = lda*n;
                    blasf77_zaxpy( &size, &c_neg_one, hA, &ione, hR, &ione );
                    error = lapackf77_zlanhe("fro", "Lower", &n_offset, &hR[offset + offset*lda], &lda, work) / error;
                    
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
                    magma_tally4_setdevice( d );
                    TESTING_FREE_DEV( dA[d] );
                    TESTING_FREE_DEV( dV[d] );
                    //TESTING_FREE_DEV( dW[d] );
                    for( int i = 0; i < nstream; ++i ) {
                        magma_tally4_queue_destroy( streams[d][i] );
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
