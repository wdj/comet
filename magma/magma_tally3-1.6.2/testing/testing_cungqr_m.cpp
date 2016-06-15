/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zungqr_m.cpp normal z -> c, Fri Jan 30 19:00:25 2015

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
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cungqr
*/
int main( int argc, char** argv )
{
    TESTING_INIT();

    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float           Anorm, error, work[1];
    magma_tally3FloatComplex  c_neg_one = MAGMA_tally3_C_NEG_ONE;
    magma_tally3FloatComplex *hA, *hR, *hT, *tau, *h_work;
    magma_tally3FloatComplex_ptr dA, dT;
    magma_tally3_int_t m, n, k;
    magma_tally3_int_t n2, lda, ldda, lwork, min_mn, nb, info;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t status = 0;
    
    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    float tol = opts.tolerance * lapackf77_slamch("E");
    
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
            nb = magma_tally3_get_cgeqrf_nb( m );
            lwork  = (m + 2*n+nb)*nb;
            gflops = FLOPS_CUNGQR( m, n, k ) / 1e9;
            
            TESTING_MALLOC_PIN( hA,     magma_tally3FloatComplex, lda*n     );
            TESTING_MALLOC_PIN( h_work, magma_tally3FloatComplex, lwork     );
            
            TESTING_MALLOC_CPU( hR,     magma_tally3FloatComplex, lda*n     );
            TESTING_MALLOC_CPU( hT,     magma_tally3FloatComplex, min_mn*nb );
            TESTING_MALLOC_CPU( tau,    magma_tally3FloatComplex, min_mn    );
                                        
            TESTING_MALLOC_DEV( dA,     magma_tally3FloatComplex, ldda*n    );
            TESTING_MALLOC_DEV( dT,     magma_tally3FloatComplex, ( 2*min_mn + ((n + 31)/32)*32 )*nb );
            
            lapackf77_clarnv( &ione, ISEED, &n2, hA );
            lapackf77_clacpy( Magma_tally3UpperLowerStr, &m, &n, hA, &lda, hR, &lda );
            
            Anorm = lapackf77_clange("f", &m, &n, hA, &lda, work );
                
            /* ====================================================================
               Performs operation using MAGMA_tally3
               =================================================================== */
            // first, get QR factors in both hA and hR, and dT in hT
            // okay that magma_tally3_cgeqrf_gpu has special structure for R; R isn't used here.
            magma_tally3_csetmatrix( m, n, hA, lda, dA, ldda );
            magma_tally3_cgeqrf_gpu( m, n, dA, ldda, tau, dT, &info );
            if ( info != 0 )
                printf("magma_tally3_cgeqrf_gpu returned error %d: %s.\n",
                       (int) info, magma_tally3_strerror( info ));
            magma_tally3_cgetmatrix( m, n, dA, ldda, hA, lda );
            lapackf77_clacpy( Magma_tally3FullStr, &m, &n, hA, &lda, hR, &lda );
            magma_tally3_cgetmatrix( nb, min_mn, dT, nb, hT, nb );
            
            gpu_time = magma_tally3_wtime();
            magma_tally3_cungqr_m( m, n, k, hR, lda, tau, hT, nb, &info );
            gpu_time = magma_tally3_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if ( info != 0 )
                printf("magma_tally3_cungqr_gpu returned error %d: %s.\n",
                       (int) info, magma_tally3_strerror( info ));
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_tally3_wtime();
                lapackf77_cungqr( &m, &n, &k, hA, &lda, tau, h_work, &lwork, &info );
                cpu_time = magma_tally3_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if ( info != 0 )
                    printf("lapackf77_cungqr returned error %d: %s.\n",
                           (int) info, magma_tally3_strerror( info ));
                
                // compute relative error |R|/|A| := |Q_magma_tally3 - Q_lapack|/|A|
                blasf77_caxpy( &n2, &c_neg_one, hA, &ione, hR, &ione );
                error = lapackf77_clange("f", &m, &n, hR, &lda, work) / Anorm;
                
                bool okay = (error < tol);
                status += ! okay;
                printf("%5d %5d %5d   %7.1f (%7.2f)   %7.1f (%7.2f)   %8.2e   %s\n",
                       (int) m, (int) n, (int) k,
                       cpu_perf, cpu_time, gpu_perf, gpu_time,
                       error, (okay ? "ok" : "failed") );
            }
            else {
                printf("%5d %5d %5d     ---   (  ---  )   %7.1f (%7.2f)     ---  \n",
                       (int) m, (int) n, (int) k,
                       gpu_perf, gpu_time );
            }
            
            TESTING_FREE_PIN( hA     );
            TESTING_FREE_PIN( h_work );
            
            TESTING_FREE_CPU( hR  );
            TESTING_FREE_CPU( hT  );
            TESTING_FREE_CPU( tau );
            
            TESTING_FREE_DEV( dA  );
            TESTING_FREE_DEV( dT  );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    
    TESTING_FINALIZE();
    return status;
}
