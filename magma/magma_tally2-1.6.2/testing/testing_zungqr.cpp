/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

       @author Stan Tomov
       @author Mathieu Faverge
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


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zungqr
*/
int main( int argc, char** argv )
{
    TESTING_INIT();

    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double           Anorm, error, work[1];
    magma_tally2DoubleComplex  c_neg_one = MAGMA_tally2_Z_NEG_ONE;
    magma_tally2DoubleComplex *hA, *hR, *tau, *h_work;
    magma_tally2DoubleComplex_ptr dA, dT;
    magma_tally2_int_t m, n, k;
    magma_tally2_int_t n2, lda, ldda, lwork, min_mn, nb, info;
    magma_tally2_int_t ione     = 1;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    magma_tally2_int_t status = 0;
    
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );

    double tol = opts.tolerance * lapackf77_dlamch("E");
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    printf("Running version %d; available are (specified through --version num):\n",
           (int) opts.version);
    printf("1 - uses precomputed zlarft matrices (default)\n");
    printf("2 - recomputes the zlarft matrices on the fly\n\n");

    printf("    m     n     k   CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R|| / ||A||\n");
    printf("=========================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            m = opts.msize[itest];
            n = opts.nsize[itest];
            k = opts.ksize[itest];
            if ( m < n || n < k ) {
                printf( "%5d %5d %5d   skipping because m < n or n < k\n", (int) m, (int) n, (int) k );
                continue;
            }
            
            lda  = m;
            ldda = ((m + 31)/32)*32;
            n2 = lda*n;
            min_mn = min(m, n);
            nb = magma_tally2_get_zgeqrf_nb( m );
            lwork  = (m + 2*n+nb)*nb;
            gflops = FLOPS_ZUNGQR( m, n, k ) / 1e9;
            
            TESTING_MALLOC_PIN( h_work, magma_tally2DoubleComplex, lwork  );
            TESTING_MALLOC_PIN( hR,     magma_tally2DoubleComplex, lda*n  );
            
            TESTING_MALLOC_CPU( hA,     magma_tally2DoubleComplex, lda*n  );
            TESTING_MALLOC_CPU( tau,    magma_tally2DoubleComplex, min_mn );
            
            TESTING_MALLOC_DEV( dA,     magma_tally2DoubleComplex, ldda*n );
            TESTING_MALLOC_DEV( dT,     magma_tally2DoubleComplex, ( 2*min_mn + ((n + 31)/32)*32 )*nb );
            
            lapackf77_zlarnv( &ione, ISEED, &n2, hA );
            lapackf77_zlacpy( Magma_tally2FullStr, &m, &n, hA, &lda, hR, &lda );
            
            Anorm = lapackf77_zlange("f", &m, &n, hA, &lda, work );
            
            /* ====================================================================
               Performs operation using MAGMA_tally2
               =================================================================== */
            // first, get QR factors in both hA and hR
            // okay that magma_tally2_zgeqrf_gpu has special structure for R; R isn't used here.
            magma_tally2_zsetmatrix( m, n, hA, lda, dA, ldda );
            magma_tally2_zgeqrf_gpu( m, n, dA, ldda, tau, dT, &info );
            if (info != 0)
                printf("magma_tally2_zgeqrf_gpu returned error %d: %s.\n",
                       (int) info, magma_tally2_strerror( info ));
            magma_tally2_zgetmatrix( m, n, dA, ldda, hA, lda );
            lapackf77_zlacpy( Magma_tally2FullStr, &m, &n, hA, &lda, hR, &lda );
            
            gpu_time = magma_tally2_wtime();
            if (opts.version == 1)
                magma_tally2_zungqr( m, n, k, hR, lda, tau, dT, nb, &info );
            else
                magma_tally2_zungqr2(m, n, k, hR, lda, tau, &info );
            gpu_time = magma_tally2_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally2_zungqr_gpu returned error %d: %s.\n",
                       (int) info, magma_tally2_strerror( info ));
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally2_wtime();
                lapackf77_zungqr( &m, &n, &k, hA, &lda, tau, h_work, &lwork, &info );
                cpu_time = magma_tally2_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_zungqr returned error %d: %s.\n",
                           (int) info, magma_tally2_strerror( info ));
                
                // compute relative error |R|/|A| := |Q_magma_tally2 - Q_lapack|/|A|
                blasf77_zaxpy( &n2, &c_neg_one, hA, &ione, hR, &ione );
                error = lapackf77_zlange("f", &m, &n, hR, &lda, work) / Anorm;
                
                bool okay = (error < tol);
                status += ! okay;
                printf("%5d %5d %5d   %7.1f (%7.2f)   %7.1f (%7.2f)   %8.2e   %s\n",
                       (int) m, (int) n, (int) k,
                       cpu_perf, cpu_time, gpu_perf, gpu_time,
                       error, (okay ? "ok" : "failed"));
            }
            else {
                printf("%5d %5d %5d     ---   (  ---  )   %7.1f (%7.2f)     ---  \n",
                       (int) m, (int) n, (int) k,
                       gpu_perf, gpu_time );
            }
            
            TESTING_FREE_PIN( h_work );
            TESTING_FREE_PIN( hR     );
            
            TESTING_FREE_CPU( hA  );
            TESTING_FREE_CPU( tau );
            
            TESTING_FREE_DEV( dA );
            TESTING_FREE_DEV( dT );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    
    TESTING_FINALIZE();
    return status;
}