/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zgetrf_gpu.cpp normal z -> c, Fri Jan 30 19:00:25 2015
       @author Mark Gates
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"
#include "testings.h"


// Initialize matrix to random.
// Having this in separate function ensures the same ISEED is always used,
// so we can re-generate the identical matrix.
void init_matrix( int m, int n, magma_tally4FloatComplex *h_A, magma_tally4_int_t lda )
{
    magma_tally4_int_t ione = 1;
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    magma_tally4_int_t n2 = lda*n;
    lapackf77_clarnv( &ione, ISEED, &n2, h_A );
}


// On input, A and ipiv is LU factorization of A. On output, A is overwritten.
// Requires m == n.
// Uses init_matrix() to re-generate original A as needed.
// Generates random RHS b and solves Ax=b.
// Returns residual, |Ax - b| / (n |A| |x|).
float get_residual(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *ipiv )
{
    if ( m != n ) {
        printf( "\nERROR: residual check defined only for square matrices\n" );
        return -1;
    }
    
    const magma_tally4FloatComplex c_one     = MAGMA_tally4_C_ONE;
    const magma_tally4FloatComplex c_neg_one = MAGMA_tally4_C_NEG_ONE;
    const magma_tally4_int_t ione = 1;
    
    // this seed should be DIFFERENT than used in init_matrix
    // (else x is column of A, so residual can be exactly zero)
    magma_tally4_int_t ISEED[4] = {0,0,0,2};
    magma_tally4_int_t info = 0;
    magma_tally4FloatComplex *x, *b;
    
    // initialize RHS
    TESTING_MALLOC_CPU( x, magma_tally4FloatComplex, n );
    TESTING_MALLOC_CPU( b, magma_tally4FloatComplex, n );
    lapackf77_clarnv( &ione, ISEED, &n, b );
    blasf77_ccopy( &n, b, &ione, x, &ione );
    
    // solve Ax = b
    lapackf77_cgetrs( "Notrans", &n, &ione, A, &lda, ipiv, x, &n, &info );
    if (info != 0)
        printf("lapackf77_cgetrs returned error %d: %s.\n",
               (int) info, magma_tally4_strerror( info ));
    
    // reset to original A
    init_matrix( m, n, A, lda );
    
    // compute r = Ax - b, saved in b
    blasf77_cgemv( "Notrans", &m, &n, &c_one, A, &lda, x, &ione, &c_neg_one, b, &ione );
    
    // compute residual |Ax - b| / (n*|A|*|x|)
    float norm_x, norm_A, norm_r, work[1];
    norm_A = lapackf77_clange( "F", &m, &n, A, &lda, work );
    norm_r = lapackf77_clange( "F", &n, &ione, b, &n, work );
    norm_x = lapackf77_clange( "F", &n, &ione, x, &n, work );
    
    //printf( "r=\n" ); magma_tally4_cprint( 1, n, b, 1 );
    
    TESTING_FREE_CPU( x );
    TESTING_FREE_CPU( b );
    
    //printf( "r=%.2e, A=%.2e, x=%.2e, n=%d\n", norm_r, norm_A, norm_x, n );
    return norm_r / (n * norm_A * norm_x);
}


// On input, LU and ipiv is LU factorization of A. On output, LU is overwritten.
// Works for any m, n.
// Uses init_matrix() to re-generate original A as needed.
// Returns error in factorization, |PA - LU| / (n |A|)
// This allocates 3 more matrices to store A, L, and U.
float get_LU_error(magma_tally4_int_t M, magma_tally4_int_t N,
                    magma_tally4FloatComplex *LU, magma_tally4_int_t lda,
                    magma_tally4_int_t *ipiv)
{
    magma_tally4_int_t min_mn = min(M,N);
    magma_tally4_int_t ione   = 1;
    magma_tally4_int_t i, j;
    magma_tally4FloatComplex alpha = MAGMA_tally4_C_ONE;
    magma_tally4FloatComplex beta  = MAGMA_tally4_C_ZERO;
    magma_tally4FloatComplex *A, *L, *U;
    float work[1], matnorm, residual;
    
    TESTING_MALLOC_CPU( A, magma_tally4FloatComplex, lda*N    );
    TESTING_MALLOC_CPU( L, magma_tally4FloatComplex, M*min_mn );
    TESTING_MALLOC_CPU( U, magma_tally4FloatComplex, min_mn*N );
    memset( L, 0, M*min_mn*sizeof(magma_tally4FloatComplex) );
    memset( U, 0, min_mn*N*sizeof(magma_tally4FloatComplex) );

    // set to original A
    init_matrix( M, N, A, lda );
    lapackf77_claswp( &N, A, &lda, &ione, &min_mn, ipiv, &ione);
    
    // copy LU to L and U, and set diagonal to 1
    lapackf77_clacpy( Magma_tally4LowerStr, &M, &min_mn, LU, &lda, L, &M      );
    lapackf77_clacpy( Magma_tally4UpperStr, &min_mn, &N, LU, &lda, U, &min_mn );
    for(j=0; j<min_mn; j++)
        L[j+j*M] = MAGMA_tally4_C_MAKE( 1., 0. );
    
    matnorm = lapackf77_clange("f", &M, &N, A, &lda, work);

    blasf77_cgemm("N", "N", &M, &N, &min_mn,
                  &alpha, L, &M, U, &min_mn, &beta, LU, &lda);

    for( j = 0; j < N; j++ ) {
        for( i = 0; i < M; i++ ) {
            LU[i+j*lda] = MAGMA_tally4_C_SUB( LU[i+j*lda], A[i+j*lda] );
        }
    }
    residual = lapackf77_clange("f", &M, &N, LU, &lda, work);

    TESTING_FREE_CPU( A );
    TESTING_FREE_CPU( L );
    TESTING_FREE_CPU( U );

    return residual / (matnorm * N);
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cgetrf
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gflops, gpu_perf, gpu_time, cpu_perf=0, cpu_time=0;
    float          error;
    magma_tally4FloatComplex *h_A;
    magma_tally4FloatComplex_ptr d_A;
    magma_tally4_int_t     *ipiv;
    magma_tally4_int_t M, N, n2, lda, ldda, info, min_mn;
    magma_tally4_int_t status   = 0;

    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );

    float tol = opts.tolerance * lapackf77_slamch("E");
    
    if ( opts.check == 2 ) {
        printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |Ax-b|/(N*|A|*|x|)\n");
    }
    else {
        printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |PA-LU|/(N*|A|)\n");
    }
    printf("=========================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N;
            ldda   = ((M+31)/32)*32;
            gflops = FLOPS_CGETRF( M, N ) / 1e9;
            
            TESTING_MALLOC_CPU( ipiv, magma_tally4_int_t,        min_mn );
            TESTING_MALLOC_CPU( h_A,  magma_tally4FloatComplex, n2     );
            TESTING_MALLOC_DEV( d_A,  magma_tally4FloatComplex, ldda*N );
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                init_matrix( M, N, h_A, lda );
                
                cpu_time = magma_tally4_wtime();
                lapackf77_cgetrf(&M, &N, h_A, &lda, ipiv, &info);
                cpu_time = magma_tally4_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_cgetrf returned error %d: %s.\n",
                           (int) info, magma_tally4_strerror( info ));
            }
            
            /* ====================================================================
               Performs operation using MAGMA_tally4
               =================================================================== */
            init_matrix( M, N, h_A, lda );
            magma_tally4_csetmatrix( M, N, h_A, lda, d_A, ldda );
            
            gpu_time = magma_tally4_wtime();
            magma_tally4_cgetrf_gpu( M, N, d_A, ldda, ipiv, &info);
            gpu_time = magma_tally4_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally4_cgetrf_gpu returned error %d: %s.\n",
                       (int) info, magma_tally4_strerror( info ));
            
            /* =====================================================================
               Check the factorization
               =================================================================== */
            if ( opts.lapack ) {
                printf("%5d %5d   %7.2f (%7.2f)   %7.2f (%7.2f)",
                       (int) M, (int) N, cpu_perf, cpu_time, gpu_perf, gpu_time );
            }
            else {
                printf("%5d %5d     ---   (  ---  )   %7.2f (%7.2f)",
                       (int) M, (int) N, gpu_perf, gpu_time );
            }
            if ( opts.check == 2 ) {
                magma_tally4_cgetmatrix( M, N, d_A, ldda, h_A, lda );
                error = get_residual( M, N, h_A, lda, ipiv );
                printf("   %8.2e   %s\n", error, (error < tol ? "ok" : "failed"));
                status += ! (error < tol);
            }
            else if ( opts.check ) {
                magma_tally4_cgetmatrix( M, N, d_A, ldda, h_A, lda );
                error = get_LU_error( M, N, h_A, lda, ipiv );
                printf("   %8.2e   %s\n", error, (error < tol ? "ok" : "failed"));
                status += ! (error < tol);
            }
            else {
                printf("     ---  \n");
            }
            
            TESTING_FREE_CPU( ipiv );
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_DEV( d_A );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
