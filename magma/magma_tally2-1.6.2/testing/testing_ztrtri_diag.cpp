/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> c d s
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

#define     h_A(i_, j_) (h_A     + (i_) + (j_)*lda)
#define h_dinvA(i_, j_) (h_dinvA + (i_) + (j_)*nb)


/* ////////////////////////////////////////////////////////////////////////////
   -- like axpy for matrices: B += alpha*A.
*/
void zgeadd(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    const magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex       *B, magma_tally2_int_t ldb )
{
    #define A(i_, j_) (A + (i_) + (j_)*lda)
    #define B(i_, j_) (B + (i_) + (j_)*ldb)
    
    const magma_tally2_int_t ione = 1;
    
    for( int j=0; j < n; ++j ) {
        blasf77_zaxpy( &m, &alpha, A(0,j), &ione, B(0,j), &ione );
    }
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing ztrtri
*/
int main( int argc, char** argv )
{
    TESTING_INIT();

    real_Double_t   gflops, magma_tally2_perf, magma_tally2_time=0;  //, cpu_perf=0, cpu_time=0;
    double          magma_tally2_error, norm_invA, work[1];
    magma_tally2_int_t N, lda, ldda, info;
    magma_tally2_int_t jb, nb, nblock, sizeA, size_inv;
    magma_tally2_int_t ione     = 1;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};
    magma_tally2_int_t *ipiv;

    magma_tally2DoubleComplex *h_A, *h_dinvA;
    magma_tally2DoubleComplex_ptr d_A, d_dinvA;
    magma_tally2DoubleComplex c_neg_one = MAGMA_tally2_Z_NEG_ONE;
    magma_tally2_int_t status = 0;
    
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    double tol = opts.tolerance * lapackf77_dlamch("E");
    const char *uplo_ = lapack_uplo_const_tally2(opts.uplo);

    // this is the NB hard coded into ztrtri_diag.
    nb = 128;
    
    printf("uplo = %s, diag = %s\n",
           lapack_uplo_const_tally2(opts.uplo), lapack_diag_const_tally2(opts.diag) );
    printf("    N  MAGMA_tally2 Gflop/s (ms)   MAGMA_tally2 error\n");
    printf("=======================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda = N;
            ldda = ((lda+31)/32)*32;
            nblock = (N+nb-1)/nb;
            gflops = nblock * FLOPS_ZTRTRI( nb ) / 1e9;
            
            TESTING_MALLOC_CPU( h_A,    magma_tally2DoubleComplex, lda*N );
            TESTING_MALLOC_CPU( ipiv,   magma_tally2_int_t,        N     );
            
            size_inv = nblock*nb*nb;
            TESTING_MALLOC_DEV( d_A,    magma_tally2DoubleComplex, ldda*N );
            TESTING_MALLOC_DEV( d_dinvA, magma_tally2DoubleComplex, size_inv );
            TESTING_MALLOC_CPU( h_dinvA, magma_tally2DoubleComplex, size_inv );
            
            /* Initialize the matrices */
            /* Factor A into LU to get well-conditioned triangular matrix.
             * Copy L to U, since L seems okay when used with non-unit diagonal
             * (i.e., from U), while U fails when used with unit diagonal. */
            sizeA = lda*N;            
            lapackf77_zlarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_zgetrf( &N, &N, h_A, &lda, ipiv, &info );
            for( int j = 0; j < N; ++j ) {
                for( int i = 0; i < j; ++i ) {
                    *h_A(i,j) = *h_A(j,i);
                }
            }
            
            /* =====================================================================
               Performs operation using MAGMA_tally2BLAS
               =================================================================== */
            magma_tally2_zsetmatrix( N, N, h_A, lda, d_A, ldda );
            
            magma_tally2_time = magma_tally2_sync_wtime( NULL );
            magma_tally2blas_ztrtri_diag( opts.uplo, opts.diag, N, d_A, ldda, d_dinvA );
            magma_tally2_time = magma_tally2_sync_wtime( NULL ) - magma_tally2_time;
            magma_tally2_perf = gflops / magma_tally2_time;
            
            magma_tally2_zgetvector( size_inv, d_dinvA, 1, h_dinvA, 1 );
            
            if ( opts.verbose ) {
                printf( "A%d=", (int) N );
                magma_tally2_zprint( N, N, h_A, lda );
                printf( "d_dinvA%d=", (int) N );
                magma_tally2_zprint( min(N+4, nb), min(N+4, nblock*nb), h_dinvA, nb );
            }
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                //cpu_time = magma_tally2_wtime();
                lapackf77_ztrtri(
                    lapack_uplo_const_tally2(opts.uplo), lapack_diag_const_tally2(opts.diag),
                    &N, h_A, &lda, &info );
                //cpu_time = magma_tally2_wtime() - cpu_time;
                //cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.check ) {
                // |invA - invA_magma_tally2| / |invA|, accumulated over all diagonal blocks
                magma_tally2_error = 0;
                norm_invA   = 0;
                for( int i=0; i < N; i += nb ) {
                    jb = min( nb, N-i );
                    zgeadd( jb, jb, c_neg_one, h_A(i, i), lda, h_dinvA(0, i), nb );
                    magma_tally2_error = max( magma_tally2_error, lapackf77_zlantr( "M", uplo_, Magma_tally2NonUnitStr, &jb, &jb, h_dinvA(0, i), &nb,  work ));
                    norm_invA   = max( norm_invA,   lapackf77_zlantr( "M", uplo_, Magma_tally2NonUnitStr, &jb, &jb, h_A(i, i),     &lda, work ));
                }
                magma_tally2_error /= norm_invA;
                
                // CPU is doing N-by-N inverse, while GPU is doing (N/NB) NB-by-NB inverses.
                // So don't compare performance.
                printf("%5d   %7.2f (%7.2f)   %8.2e   %s\n",
                        (int) N,
                        magma_tally2_perf,  1000.*magma_tally2_time,
                        //cpu_perf,    1000.*cpu_time,
                        magma_tally2_error,
                        (magma_tally2_error < tol ? "ok" : "failed"));
                status += ! (magma_tally2_error < tol);
            }
            else {
                printf("%5d   %7.2f (%7.2f)      ---\n",
                        (int) N,
                        magma_tally2_perf,  1000.*magma_tally2_time );
            }
            
            TESTING_FREE_CPU( h_A     );
            TESTING_FREE_CPU( ipiv    );
            
            TESTING_FREE_DEV( d_A     );
            TESTING_FREE_DEV( d_dinvA );
            TESTING_FREE_CPU( h_dinvA );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
