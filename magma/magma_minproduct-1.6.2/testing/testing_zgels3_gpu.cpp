/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

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

#define PRECISION_z

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zgels
*/
int main( int argc, char** argv)
{
    TESTING_INIT();
    
    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    double           gpu_error, cpu_error, error, Anorm, work[1];
    magma_minproductDoubleComplex  c_one     = MAGMA_minproduct_Z_ONE;
    magma_minproductDoubleComplex  c_neg_one = MAGMA_minproduct_Z_NEG_ONE;
    magma_minproductDoubleComplex *h_A, *h_A2, *h_B, *h_X, *h_R, *tau, *h_work, tmp[1];
    magma_minproductDoubleComplex_ptr d_A, d_B;
    magma_minproduct_int_t M, N, size, nrhs, lda, ldb, ldda, lddb, min_mn, max_mn, nb, info;
    magma_minproduct_int_t lworkgpu, lhwork, lhwork2;
    magma_minproduct_int_t ione     = 1;
    magma_minproduct_int_t ISEED[4] = {0,0,0,1};

    magma_minproduct_opts opts;
    parse_opts( argc, argv, &opts );
 
    magma_minproduct_int_t status = 0;
    double tol = opts.tolerance * lapackf77_dlamch("E");

    nrhs = opts.nrhs;
    
    printf("                                                            ||b-Ax|| / (N||A||)   ||dx-x||/(N||A||)\n");
    printf("    M     N  NRHS   CPU GFlop/s (sec)   GPU GFlop/s (sec)   CPU        GPU                         \n");
    printf("===================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            if ( M < N ) {
                printf( "%5d %5d %5d   skipping because M < N is not yet supported.\n", (int) M, (int) N, (int) nrhs );
                continue;
            }
            min_mn = min(M, N);
            max_mn = max(M, N);
            lda    = M;
            ldb    = max_mn;
            size   = lda*N;
            ldda   = ((M+31)/32)*32;
            lddb   = ((max_mn+31)/32)*32;
            nb     = magma_minproduct_get_zgeqrf_nb(M);
            gflops = (FLOPS_ZGEQRF( M, N ) + FLOPS_ZGEQRS( M, N, nrhs )) / 1e9;
            
            lworkgpu = (M - N + nb)*(nrhs + nb) + nrhs*nb;
            
            // query for workspace size
            lhwork = -1;
            lapackf77_zgeqrf(&M, &N, NULL, &M, NULL, tmp, &lhwork, &info);
            lhwork2 = (magma_minproduct_int_t) MAGMA_minproduct_Z_REAL( tmp[0] );
            
            lhwork = -1;
            lapackf77_zunmqr( Magma_minproductLeftStr, Magma_minproduct_ConjTransStr,
                              &M, &nrhs, &min_mn, NULL, &lda, NULL,
                              NULL, &ldb, tmp, &lhwork, &info);
            lhwork = (magma_minproduct_int_t) MAGMA_minproduct_Z_REAL( tmp[0] );
            lhwork = max( max( lhwork, lhwork2 ), lworkgpu );
            
            TESTING_MALLOC_CPU( tau,    magma_minproductDoubleComplex, min_mn    );
            TESTING_MALLOC_CPU( h_A,    magma_minproductDoubleComplex, lda*N     );
            TESTING_MALLOC_CPU( h_A2,   magma_minproductDoubleComplex, lda*N     );
            TESTING_MALLOC_CPU( h_B,    magma_minproductDoubleComplex, ldb*nrhs  );
            TESTING_MALLOC_CPU( h_X,    magma_minproductDoubleComplex, ldb*nrhs  );
            TESTING_MALLOC_CPU( h_R,    magma_minproductDoubleComplex, ldb*nrhs  );
            TESTING_MALLOC_CPU( h_work, magma_minproductDoubleComplex, lhwork    );
            
            TESTING_MALLOC_DEV( d_A,    magma_minproductDoubleComplex, ldda*N    );
            TESTING_MALLOC_DEV( d_B,    magma_minproductDoubleComplex, lddb*nrhs );
            
            /* Initialize the matrices */
            lapackf77_zlarnv( &ione, ISEED, &size, h_A );
            lapackf77_zlacpy( Magma_minproductUpperLowerStr, &M, &N, h_A, &lda, h_A2, &lda );
            
            // make random RHS
            size = M*nrhs;
            lapackf77_zlarnv( &ione, ISEED, &size, h_B );
            lapackf77_zlacpy( Magma_minproductUpperLowerStr, &M, &nrhs, h_B, &ldb, h_R, &ldb );
            
            // make consistent RHS
            //size = N*nrhs;
            //lapackf77_zlarnv( &ione, ISEED, &size, h_X );
            //blasf77_zgemm( Magma_minproductNoTransStr, Magma_minproductNoTransStr, &M, &nrhs, &N,
            //               &c_one,  h_A, &lda,
            //                        h_X, &ldb,
            //               &c_zero, h_B, &ldb );
            //lapackf77_zlacpy( Magma_minproductUpperLowerStr, &M, &nrhs, h_B, &ldb, h_R, &ldb );
            
            /* ====================================================================
               Performs operation using MAGMA_minproduct
               =================================================================== */
            magma_minproduct_zsetmatrix( M, N,    h_A, lda, d_A, ldda );
            magma_minproduct_zsetmatrix( M, nrhs, h_B, ldb, d_B, lddb );
            
            gpu_time = magma_minproduct_wtime();
            magma_minproduct_zgels3_gpu( Magma_minproductNoTrans, M, N, nrhs, d_A, ldda,
                              d_B, lddb, h_work, lworkgpu, &info);
            gpu_time = magma_minproduct_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_minproduct_zgels3_gpu returned error %d: %s.\n",
                       (int) info, magma_minproduct_strerror( info ));
            
            // Get the solution in h_X
            magma_minproduct_zgetmatrix( N, nrhs, d_B, lddb, h_X, ldb );
            
            // compute the residual
            blasf77_zgemm( Magma_minproductNoTransStr, Magma_minproductNoTransStr, &M, &nrhs, &N,
                           &c_neg_one, h_A, &lda,
                                       h_X, &ldb,
                           &c_one,     h_R, &ldb);
            Anorm = lapackf77_zlange("f", &M, &N, h_A, &lda, work);
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            lapackf77_zlacpy( Magma_minproductUpperLowerStr, &M, &nrhs, h_B, &ldb, h_X, &ldb );
            
            cpu_time = magma_minproduct_wtime();
            lapackf77_zgels( Magma_minproductNoTransStr, &M, &N, &nrhs,
                             h_A, &lda, h_X, &ldb, h_work, &lhwork, &info);
            cpu_time = magma_minproduct_wtime() - cpu_time;
            cpu_perf = gflops / cpu_time;
            if (info != 0)
                printf("lapackf77_zgels returned error %d: %s.\n",
                       (int) info, magma_minproduct_strerror( info ));
            
            blasf77_zgemm( Magma_minproductNoTransStr, Magma_minproductNoTransStr, &M, &nrhs, &N,
                           &c_neg_one, h_A2, &lda,
                                       h_X,  &ldb,
                           &c_one,     h_B,  &ldb);
            
            cpu_error = lapackf77_zlange("f", &M, &nrhs, h_B, &ldb, work) / (min_mn*Anorm);
            gpu_error = lapackf77_zlange("f", &M, &nrhs, h_R, &ldb, work) / (min_mn*Anorm);
            
            // error relative to LAPACK
            size = M*nrhs;
            blasf77_zaxpy( &size, &c_neg_one, h_B, &ione, h_R, &ione );
            error = lapackf77_zlange("f", &M, &nrhs, h_R, &ldb, work) / (min_mn*Anorm);
            
            printf("%5d %5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %8.2e   %8.2e",
                   (int) M, (int) N, (int) nrhs,
                   cpu_perf, cpu_time, gpu_perf, gpu_time, cpu_error, gpu_error, error );
                        
            if ( M == N ) {
                printf( "   %s\n", (gpu_error < tol && error < tol ? "ok" : "failed"));
                status += ! (gpu_error < tol && error < tol);
            }
            else {
                printf( "   %s\n", (error < tol ? "ok" : "failed"));
                status += ! (error < tol);
            }

            TESTING_FREE_CPU( tau    );
            TESTING_FREE_CPU( h_A    );
            TESTING_FREE_CPU( h_A2   );
            TESTING_FREE_CPU( h_B    );
            TESTING_FREE_CPU( h_X    );
            TESTING_FREE_CPU( h_R    );
            TESTING_FREE_CPU( h_work );
            
            TESTING_FREE_DEV( d_A    );
            TESTING_FREE_DEV( d_B    );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
