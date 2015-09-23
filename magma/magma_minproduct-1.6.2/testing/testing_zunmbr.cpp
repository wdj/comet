/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @precisions normal z -> c d s
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

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zunmbr
*/
int main( int argc, char** argv )
{
    TESTING_INIT();
    
    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double error, dwork[1];
    magma_minproductDoubleComplex c_neg_one = MAGMA_minproduct_Z_NEG_ONE;
    magma_minproduct_int_t ione = 1;
    magma_minproduct_int_t m, n, k, mi, ni, mm, nn, nq, size, info;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};
    magma_minproduct_int_t nb, ldc, lda, lwork, lwork_max;
    magma_minproductDoubleComplex *C, *R, *A, *work, *tau, *tauq, *taup;
    double *d, *e;
    magma_minproduct_int_t status = 0;
    
    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
    
    // need slightly looser bound (60*eps instead of 30*eps) for some tests
    opts.tolerance = max( 60., opts.tolerance );
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    // test all combinations of input parameters
    magma_minproduct_vect_t  vect [] = { Magma_minproductQ,          Magma_minproductP       };
    magma_minproduct_side_t  side [] = { Magma_minproductLeft,       Magma_minproductRight   };
    magma_minproduct_trans_t trans[] = { Magma_minproduct_ConjTrans, Magma_minproductNoTrans };

    printf("    M     N     K   vect side   trans   CPU GFlop/s (sec)   GPU GFlop/s (sec)   ||R||_F / ||QC||_F\n");
    printf("===============================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
      for( int ivect = 0; ivect < 2; ++ivect ) {
      for( int iside = 0; iside < 2; ++iside ) {
      for( int itran = 0; itran < 2; ++itran ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            m = opts.msize[itest];
            n = opts.nsize[itest];
            k = opts.ksize[itest];
            nb  = magma_minproduct_get_zgebrd_nb( m );
            ldc = m;
            // A is nq x k (vect=Q) or k x nq (vect=P)
            // where nq=m (left) or nq=n (right)
            nq  = (side[iside] == Magma_minproductLeft ? m  : n );
            mm  = (vect[ivect] == Magma_minproductQ    ? nq : k );
            nn  = (vect[ivect] == Magma_minproductQ    ? k  : nq);
            lda = mm;
            
            // MBR calls either MQR or MLQ in various ways
            if ( vect[ivect] == Magma_minproductQ ) {
                if ( nq >= k ) {
                    gflops = FLOPS_ZUNMQR( m, n, k, side[iside] ) / 1e9;
                }
                else {
                    if ( side[iside] == Magma_minproductLeft ) {
                        mi = m - 1;
                        ni = n;
                    }
                    else {
                        mi = m;
                        ni = n - 1;
                    }
                    gflops = FLOPS_ZUNMQR( mi, ni, nq-1, side[iside] ) / 1e9;
                }
            }
            else {
                if ( nq > k ) {
                    gflops = FLOPS_ZUNMLQ( m, n, k, side[iside] ) / 1e9;
                }
                else {
                    if ( side[iside] == Magma_minproductLeft ) {
                        mi = m - 1;
                        ni = n;
                    }
                    else {
                        mi = m;
                        ni = n - 1;
                    }
                    gflops = FLOPS_ZUNMLQ( mi, ni, nq-1, side[iside] ) / 1e9;
                }
            }
            
            // workspace for gebrd is (mm + nn)*nb
            // workspace for unmbr is m*nb or n*nb, depending on side
            lwork_max = max( (mm + nn)*nb, max( m*nb, n*nb ));
            
            TESTING_MALLOC_CPU( C,    magma_minproductDoubleComplex, ldc*n );
            TESTING_MALLOC_CPU( R,    magma_minproductDoubleComplex, ldc*n );
            TESTING_MALLOC_CPU( A,    magma_minproductDoubleComplex, lda*nn );
            TESTING_MALLOC_CPU( work, magma_minproductDoubleComplex, lwork_max );
            TESTING_MALLOC_CPU( d,    double,             min(mm,nn) );
            TESTING_MALLOC_CPU( e,    double,             min(mm,nn) );
            TESTING_MALLOC_CPU( tauq, magma_minproductDoubleComplex, min(mm,nn) );
            TESTING_MALLOC_CPU( taup, magma_minproductDoubleComplex, min(mm,nn) );
            
            // C is full, m x n
            size = ldc*n;
            lapackf77_zlarnv( &ione, ISEED, &size, C );
            lapackf77_zlacpy( "Full", &m, &n, C, &ldc, R, &ldc );
            
            size = lda*nn;
            lapackf77_zlarnv( &ione, ISEED, &size, A );
            
            // compute BRD factorization to get Householder vectors in A, tauq, taup
            //lapackf77_zgebrd( &mm, &nn, A, &lda, d, e, tauq, taup, work, &lwork_max, &info );
            magma_minproduct_zgebrd( mm, nn, A, lda, d, e, tauq, taup, work, lwork_max, &info );
            if (info != 0)
                printf("magma_minproduct_zgebrd returned error %d: %s.\n",
                       (int) info, magma_minproduct_strerror( info ));
            
            if ( vect[ivect] == Magma_minproductQ ) {
                tau = tauq;
            } else {
                tau = taup;
            }
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_minproduct_wtime();
            lapackf77_zunmbr( lapack_vect_const( vect[ivect] ),
                              lapack_side_const( side[iside] ),
                              lapack_trans_const( trans[itran] ),
                              &m, &n, &k,
                              A, &lda, tau, C, &ldc, work, &lwork_max, &info );
            cpu_time = magma_minproduct_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            if (info != 0)
                printf("lapackf77_zunmbr returned error %d: %s.\n",
                       (int) info, magma_minproduct_strerror( info ));
            
            /* ====================================================================
               Performs operation using MAGMA_minproduct
               =================================================================== */
            // query for workspace size
            lwork = -1;
            magma_minproduct_zunmbr( vect[ivect], side[iside], trans[itran],
                          m, n, k,
                          A, lda, tau, R, ldc, work, lwork, &info );
            if (info != 0)
                printf("magma_minproduct_zunmbr (lwork query) returned error %d: %s.\n",
                       (int) info, magma_minproduct_strerror( info ));
            lwork = (magma_minproduct_int_t) MAGMA_minproduct_Z_REAL( work[0] );
            if ( lwork < 0 || lwork > lwork_max ) {
                printf("optimal lwork %d > lwork_max %d\n", (int) lwork, (int) lwork_max );
                lwork = lwork_max;
            }
            
            gpu_time = magma_minproduct_wtime();
            magma_minproduct_zunmbr( vect[ivect], side[iside], trans[itran],
                          m, n, k,
                          A, lda, tau, R, ldc, work, lwork, &info );
            gpu_time = magma_minproduct_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_minproduct_zunmbr returned error %d: %s.\n",
                       (int) info, magma_minproduct_strerror( info ));
            
            /* =====================================================================
               compute relative error |QC_magma_minproduct - QC_lapack| / |QC_lapack|
               =================================================================== */
            error = lapackf77_zlange( "Fro", &m, &n, C, &ldc, dwork );
            size = ldc*n;
            blasf77_zaxpy( &size, &c_neg_one, C, &ione, R, &ione );
            error = lapackf77_zlange( "Fro", &m, &n, R, &ldc, dwork ) / error;
            
            printf( "%5d %5d %5d   %c   %4c   %5c   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                    (int) m, (int) n, (int) k,
                    lapacke_vect_const( vect[ivect] ),
                    lapacke_side_const( side[iside] ),
                    lapacke_trans_const( trans[itran] ),
                    cpu_perf, cpu_time, gpu_perf, gpu_time,
                    error, (error < tol ? "ok" : "failed") );
            status += ! (error < tol);
            
            TESTING_FREE_CPU( C );
            TESTING_FREE_CPU( R );
            TESTING_FREE_CPU( A );
            TESTING_FREE_CPU( work );
            TESTING_FREE_CPU( d );
            TESTING_FREE_CPU( e );
            TESTING_FREE_CPU( taup );
            TESTING_FREE_CPU( tauq );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
      }}}  // end ivect, iside, itran
      printf( "\n" );
    }
    
    TESTING_FINALIZE();
    return status;
}
