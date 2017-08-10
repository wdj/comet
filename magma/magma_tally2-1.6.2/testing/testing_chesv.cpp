/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zhesv.cpp normal z -> c, Fri Jan 30 19:00:24 2015
       @author Ichitaro Yamazaki
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_tally2.h"
#include "magma_tally2_lapack.h"
#include "testings.h"

/* ================================================================================================== */

// Initialize matrix to random & symmetrize.
// Having this in separate function ensures the same ISEED is always used,
// so we can re-generate the identical matrix.
void init_matrix( int m, int n, magma_tally2FloatComplex *h_A, magma_tally2_int_t lda )
{
    assert( m == n );
    magma_tally2_int_t ione = 1;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    magma_tally2_int_t n2 = lda*n;
    lapackf77_clarnv( &ione, ISEED, &n2, h_A );
    magma_tally2_cmake_hermitian( n, h_A, lda );
}


// On input, A and ipiv is LU factorization of A. On output, A is overwritten.
// Requires m == n.
// Uses init_matrix() to re-generate original A as needed.
// Generates random RHS b and solves Ax=b.
// Returns residual, |Ax - b| / (n |A| |x|).
float get_residual(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2FloatComplex *x, magma_tally2_int_t ldx,
    magma_tally2FloatComplex *b, magma_tally2_int_t ldb)
{
    const magma_tally2FloatComplex c_one     = MAGMA_tally2_C_ONE;
    const magma_tally2FloatComplex c_neg_one = MAGMA_tally2_C_NEG_ONE;
    const magma_tally2_int_t ione = 1;
    
    // reset to original A
    init_matrix( n, n, A, lda );
    
    // compute r = Ax - b, saved in b
    blasf77_cgemv( "Notrans", &n, &n, &c_one, A, &lda, x, &ione, &c_neg_one, b, &ione );
    
    // compute residual |Ax - b| / (n*|A|*|x|)
    float norm_x, norm_A, norm_r, work[1];
    norm_A = lapackf77_clange( Magma_tally2FullStr, &n, &n, A, &lda, work );
    norm_r = lapackf77_clange( Magma_tally2FullStr, &n, &ione, b, &n, work );
    norm_x = lapackf77_clange( Magma_tally2FullStr, &n, &ione, x, &n, work );
    
    //printf( "r=\n" ); magma_tally2_cprint( 1, n, b, 1 );
    //printf( "r=%.2e, A=%.2e, x=%.2e, n=%d\n", norm_r, norm_A, norm_x, n );
    return norm_r / (n * norm_A * norm_x);
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing chesv
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    magma_tally2FloatComplex *h_A, *h_B, *h_X, *work, temp;
    real_Double_t   gflops, gpu_perf, gpu_time = 0.0, cpu_perf=0, cpu_time=0;
    float          error, error_lapack = 0.0;
    magma_tally2_int_t     *ipiv;
    magma_tally2_int_t     N, n2, lda, ldb, sizeB, lwork, info;
    magma_tally2_int_t     status = 0, ione = 1;
    magma_tally2_int_t     ISEED[4] = {0,0,0,1};

    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");

    printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |Ax-b|/(N*|A|*|x|)\n");
    printf("=========================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            ldb    = N;
            lda    = N;
            n2     = lda*N;
            sizeB  = ldb*opts.nrhs;
            gflops = ( FLOPS_CPOTRF( N ) + FLOPS_CPOTRS( N, opts.nrhs ) ) / 1e9;
            
            TESTING_MALLOC_CPU( ipiv, magma_tally2_int_t, N );
            TESTING_MALLOC_PIN( h_A,  magma_tally2FloatComplex, n2 );
            TESTING_MALLOC_PIN( h_B,  magma_tally2FloatComplex, sizeB );
            TESTING_MALLOC_PIN( h_X,  magma_tally2FloatComplex, sizeB );
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                lwork = -1;
                lapackf77_chesv(lapack_uplo_const_tally2(opts.uplo), &N, &opts.nrhs, 
                                h_A, &lda, ipiv, h_X, &ldb, &temp, &lwork, &info);
                lwork = (int)MAGMA_tally2_C_REAL(temp);
                TESTING_MALLOC_CPU( work, magma_tally2FloatComplex, lwork );

                init_matrix( N, N, h_A, lda );
                lapackf77_clarnv( &ione, ISEED, &sizeB, h_B );
                lapackf77_clacpy( Magma_tally2UpperLowerStr, &N, &opts.nrhs, h_B, &ldb, h_X, &ldb );

                cpu_time = magma_tally2_wtime();
                lapackf77_chesv(lapack_uplo_const_tally2(opts.uplo), &N, &opts.nrhs,
                                h_A, &lda, ipiv, h_X, &ldb, work, &lwork, &info);
                cpu_time = magma_tally2_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_chesv returned error %d: %s.\n",
                           (int) info, magma_tally2_strerror( info ));
                error_lapack = get_residual( opts.uplo, N, opts.nrhs, h_A, lda, ipiv, h_X, ldb, h_B, ldb );

                TESTING_FREE_CPU( work );
            }
           
            /* ====================================================================
               Performs operation using MAGMA_tally2
               =================================================================== */
            init_matrix( N, N, h_A, lda );
            lapackf77_clarnv( &ione, ISEED, &sizeB, h_B );
            lapackf77_clacpy( Magma_tally2UpperLowerStr, &N, &opts.nrhs, h_B, &ldb, h_X, &ldb );

            magma_tally2_setdevice(0);
            gpu_time = magma_tally2_wtime();
            magma_tally2_chesv( opts.uplo, N, opts.nrhs, h_A, lda, ipiv, h_X, ldb, &info);
            gpu_time = magma_tally2_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally2_chesv returned error %d: %s.\n",
                       (int) info, magma_tally2_strerror( info ));
            
            /* =====================================================================
               Check the factorization
               =================================================================== */
            if ( opts.lapack ) {
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)",
                       (int) N, (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time );
            }
            else {
                printf("%5d %5d     ---   (  ---  )   %7.2f (%7.2f)",
                       (int) N, (int) N, gpu_perf, gpu_time );
            }
            if ( opts.check == 0 ) {
                printf("     ---   \n");
            } else {
                error = get_residual( opts.uplo, N, opts.nrhs, h_A, lda, ipiv, h_X, ldb, h_B, ldb );
                printf("   %8.2e   %s", error, (error < tol ? "ok" : "failed"));
                if (opts.lapack)
                    printf(" (lapack rel.res. = %8.2e)", error_lapack);
                printf("\n");
                status += ! (error < tol);
            }
            
            TESTING_FREE_CPU( ipiv );
            TESTING_FREE_PIN( h_X  );
            TESTING_FREE_PIN( h_B  );
            TESTING_FREE_PIN( h_A  );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}