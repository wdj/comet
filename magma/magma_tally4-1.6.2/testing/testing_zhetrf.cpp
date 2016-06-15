/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> c d s
       @author Ichitaro Yamazaki
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

/* ================================================================================================== */

// Initialize matrix to random & symmetrize. If nopiv, make positive definite.
// Having this in separate function ensures the same ISEED is always used,
// so we can re-generate the identical matrix.
void init_matrix( int nopiv, int m, int n, magma_tally4DoubleComplex *h_A, magma_tally4_int_t lda )
{
    assert( m == n );
    magma_tally4_int_t ione = 1;
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    magma_tally4_int_t n2 = lda*n;
    lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
    if (nopiv) {
        magma_tally4_zmake_hpd( n, h_A, lda );
    }
    else {
        magma_tally4_zmake_hermitian( n, h_A, lda );
    }
}


// On input, A and ipiv is LU factorization of A. On output, A is overwritten.
// Requires m == n.
// Uses init_matrix() to re-generate original A as needed.
// Generates random RHS b and solves Ax=b.
// Returns residual, |Ax - b| / (n |A| |x|).
double get_residual(
    int nopiv, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *ipiv )
{
    const magma_tally4DoubleComplex c_one     = MAGMA_tally4_Z_ONE;
    const magma_tally4DoubleComplex c_neg_one = MAGMA_tally4_Z_NEG_ONE;
    const magma_tally4_int_t ione = 1;
    magma_tally4_int_t upper = (uplo == Magma_tally4Upper);
    
    // this seed should be DIFFERENT than used in init_matrix
    // (else x is column of A, so residual can be exactly zero)
    magma_tally4_int_t ISEED[4] = {0,0,0,2};
    magma_tally4_int_t info = 0;
    magma_tally4DoubleComplex *x, *b;
    
    // initialize RHS
    TESTING_MALLOC_CPU( x, magma_tally4DoubleComplex, n );
    TESTING_MALLOC_CPU( b, magma_tally4DoubleComplex, n );
    lapackf77_zlarnv( &ione, ISEED, &n, b );
    blasf77_zcopy( &n, b, &ione, x, &ione );
    
    // solve Ax = b
    if (nopiv) {
        if (upper) {
            blasf77_ztrsm( Magma_tally4LeftStr, Magma_tally4UpperStr, 
                           Magma_tally4ConjTransStr, Magma_tally4UnitStr, 
                           &n, &ione, &c_one,
                           A, &lda, x, &n );
            for (int i=0; i<n; i++) x[i] = MAGMA_tally4_Z_DIV( x[i], A[i+i*lda] );
            blasf77_ztrsm( Magma_tally4LeftStr, Magma_tally4UpperStr, 
                           Magma_tally4NoTransStr, Magma_tally4UnitStr, 
                           &n, &ione, &c_one,
                           A, &lda, x, &n );
        } else {
            blasf77_ztrsm( Magma_tally4LeftStr, Magma_tally4LowerStr, 
                           Magma_tally4NoTransStr, Magma_tally4UnitStr, 
                           &n, &ione, &c_one,
                           A, &lda, x, &n );
            for (int i=0; i<n; i++) x[i] = MAGMA_tally4_Z_DIV( x[i], A[i+i*lda] );
            blasf77_ztrsm( Magma_tally4LeftStr, Magma_tally4LowerStr, 
                           Magma_tally4ConjTransStr, Magma_tally4UnitStr, 
                           &n, &ione, &c_one,
                           A, &lda, x, &n );
        }
    }else {
        lapackf77_zhetrs( lapack_uplo_const_tally4(uplo), &n, &ione, A, &lda, ipiv, x, &n, &info );
    }
    if (info != 0)
        printf("lapackf77_zhetrs returned error %d: %s.\n",
               (int) info, magma_tally4_strerror( info ));
    // reset to original A
    init_matrix( nopiv, n, n, A, lda );
    
    // compute r = Ax - b, saved in b
    blasf77_zgemv( "Notrans", &n, &n, &c_one, A, &lda, x, &ione, &c_neg_one, b, &ione );
    
    // compute residual |Ax - b| / (n*|A|*|x|)
    double norm_x, norm_A, norm_r, work[1];
    norm_A = lapackf77_zlange( Magma_tally4FullStr, &n, &n, A, &lda, work );
    norm_r = lapackf77_zlange( Magma_tally4FullStr, &n, &ione, b, &n, work );
    norm_x = lapackf77_zlange( Magma_tally4FullStr, &n, &ione, x, &n, work );
    
    //printf( "r=\n" ); magma_tally4_zprint( 1, n, b, 1 );
    
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
double get_LDLt_error(int nopiv, magma_tally4_uplo_t uplo, magma_tally4_int_t N,
                      magma_tally4DoubleComplex *LD, magma_tally4_int_t lda,
                      magma_tally4_int_t *ipiv)
{
    magma_tally4_int_t i, j;
    magma_tally4DoubleComplex c_one  = MAGMA_tally4_Z_ONE;
    magma_tally4DoubleComplex c_zero = MAGMA_tally4_Z_ZERO;
    magma_tally4DoubleComplex *A, *L, *D;
    double work[1], matnorm, residual;
    #define LD(i,j) (LD[(i) + (j)*lda])
    #define  A(i,j) ( A[(i) + (j)*N])
    #define  L(i,j) ( L[(i) + (j)*N])
    #define  D(i,j) ( D[(i) + (j)*N])

    TESTING_MALLOC_CPU( A, magma_tally4DoubleComplex, N*N );
    TESTING_MALLOC_CPU( L, magma_tally4DoubleComplex, N*N );
    TESTING_MALLOC_CPU( D, magma_tally4DoubleComplex, N*N );
    memset( L, 0, N*N*sizeof(magma_tally4DoubleComplex) );
    memset( D, 0, N*N*sizeof(magma_tally4DoubleComplex) );

    // set to original A, and apply pivoting
    init_matrix( nopiv, N, N, A, N );
    if (uplo == Magma_tally4Upper) {
        for (j=N-1; j>=0; j--) {
            int piv = (nopiv ? j+1 : ipiv[j]);
            if (piv < 0) {
                piv = -(piv+1);
                // extract 2-by-2 pivot
                D(j,j)     = LD(j,j);
                D(j,j-1)   = MAGMA_tally4_Z_CNJG(LD(j-1,j));
                D(j-1,j)   = LD(j-1,j);
                D(j-1,j-1) = LD(j-1,j-1);
                // exract L
                L(j,j) = c_one;
                for (i=0; i<j-1; i++) {
                    L(i,j) = LD(i,j);
                }
                j--;
                L(j,j) = c_one;
                for (i=0; i<j; i++) {
                    L(i,j) = LD(i,j);
                }
                if (piv != j) {
                    // apply row-pivoting to previous L
                    for (i=j+2; i<N; i++) {
                        magma_tally4DoubleComplex val = L(j,i);
                        L(j,i) = L(piv,i);
                        L(piv,i) = val;
                    }
                    // apply row-pivoting to A
                    for (i=0; i<N; i++) {
                        magma_tally4DoubleComplex val = A(j,i);
                        A(j,i) = A(piv,i);
                        A(piv,i) = val;
                    }
                    // apply col-pivoting to A
                    for (i=0; i<N; i++) {
                        magma_tally4DoubleComplex val = A(i,j);
                        A(i,j) = A(i,piv);
                        A(i,piv) = val;
                    }
                }
            } else {
                piv = piv-1;
                // extract 1-by-1 pivot
                D(j,j) = LD(j,j);
                // exract L
                L(j,j) = c_one;
                for (i=0; i<j; i++) {
                    L(i,j) = LD(i,j);
                }
                if (piv != j) {
                    // apply row-pivoting to previous L
                    for (i=j+1; i<N; i++) {
                        magma_tally4DoubleComplex val = L(j,i);
                        L(j,i) = L(piv,i);
                        L(piv,i) = val;
                    }
                    // apply row-pivoting to A
                    for (i=0; i<N; i++) {
                        magma_tally4DoubleComplex val = A(j,i);
                        A(j,i) = A(piv,i);
                        A(piv,i) = val;
                    }
                    // apply col-pivoting to A
                    for (i=0; i<N; i++) {
                        magma_tally4DoubleComplex val = A(i,j);
                        A(i,j) = A(i,piv);
                        A(i,piv) = val;
                    }
                }
            }
        }
        if (nopiv) {
            // compute W = D*U
            blasf77_zgemm(Magma_tally4NoTransStr, Magma_tally4NoTransStr, &N, &N, &N,
                          &c_one, D, &N, L, &N, &c_zero, LD, &lda);
            // compute D = U'*W
            blasf77_zgemm(Magma_tally4ConjTransStr, Magma_tally4NoTransStr, &N, &N, &N,
                          &c_one, L, &N, LD, &lda, &c_zero, D, &N);
        } else {
            // compute W = U*D
            blasf77_zgemm(Magma_tally4NoTransStr, Magma_tally4NoTransStr, &N, &N, &N,
                          &c_one, L, &N, D, &N, &c_zero, LD, &lda);
            // compute D = W*U'
            blasf77_zgemm(Magma_tally4NoTransStr, Magma_tally4ConjTransStr, &N, &N, &N,
                          &c_one, LD, &lda, L, &N, &c_zero, D, &N);
        }
    } else {
        for (j=0; j<N; j++) {
            int piv = (nopiv ? j+1 : ipiv[j]);
            if (piv < 0) {
                piv = -(piv+1);
                // extract 2-by-2 pivot
                D(j,j)     = LD(j,j);
                D(j,j+1)   = MAGMA_tally4_Z_CNJG(LD(j+1,j));
                D(j+1,j)   = LD(j+1,j);
                D(j+1,j+1) = LD(j+1,j+1);
                // exract L
                L(j,j) = c_one;
                for (i=j+2; i<N; i++) {
                    L(i,j) = LD(i,j);
                }
                j++;
                L(j,j) = c_one;
                for (i=j+1; i<N; i++) {
                    L(i,j) = LD(i,j);
                }
                if (piv != j) {
                    // apply row-pivoting to previous L
                    for (i=0; i<j-1; i++) {
                        magma_tally4DoubleComplex val = L(j,i);
                        L(j,i) = L(piv,i);
                        L(piv,i) = val;
                    }
                    // apply row-pivoting to A
                    for (i=0; i<N; i++) {
                        magma_tally4DoubleComplex val = A(j,i);
                        A(j,i) = A(piv,i);
                        A(piv,i) = val;
                    }
                    // apply col-pivoting to A
                    for (i=0; i<N; i++) {
                        magma_tally4DoubleComplex val = A(i,j);
                        A(i,j) = A(i,piv);
                        A(i,piv) = val;
                    }
                }
            } else {
                piv = piv-1;
                // extract 1-by-1 pivot
                D(j,j) = LD(j,j);
                // exract L
                L(j,j) = c_one;
                for (i=j+1; i<N; i++) {
                    L(i,j) = LD(i,j);
                }
                if (piv != j) {
                    // apply row-pivoting to previous L
                    for (i=0; i<j; i++) {
                        magma_tally4DoubleComplex val = L(j,i);
                        L(j,i) = L(piv,i);
                        L(piv,i) = val;
                    }
                    // apply row-pivoting to A
                    for (i=0; i<N; i++) {
                        magma_tally4DoubleComplex val = A(j,i);
                        A(j,i) = A(piv,i);
                        A(piv,i) = val;
                    }
                    // apply col-pivoting to A
                    for (i=0; i<N; i++) {
                        magma_tally4DoubleComplex val = A(i,j);
                        A(i,j) = A(i,piv);
                        A(i,piv) = val;
                    }
                }
            }
        }
        // compute W = L*D
        blasf77_zgemm(Magma_tally4NoTransStr, Magma_tally4NoTransStr, &N, &N, &N,
                      &c_one, L, &N, D, &N, &c_zero, LD, &lda);
        // compute D = W*L'
        blasf77_zgemm(Magma_tally4NoTransStr, Magma_tally4ConjTransStr, &N, &N, &N,
                      &c_one, LD, &lda, L, &N, &c_zero, D, &N);
    }
    // compute norm of A
    matnorm = lapackf77_zlange(Magma_tally4FullStr, &N, &N, A, &lda, work);

    for( j = 0; j < N; j++ ) {
        for( i = 0; i < N; i++ ) {
            D(i,j) = MAGMA_tally4_Z_SUB( D(i,j), A(i,j) );
        }
    }
    residual = lapackf77_zlange(Magma_tally4FullStr, &N, &N, D, &N, work);

    TESTING_FREE_CPU( A );
    TESTING_FREE_CPU( L );
    TESTING_FREE_CPU( D );

    return residual / (matnorm * N);
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zhetrf
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    magma_tally4DoubleComplex *h_A, *work, temp;
    real_Double_t   gflops, gpu_perf, gpu_time = 0.0, cpu_perf=0, cpu_time=0;
    double          error, error_lapack = 0.0;
    magma_tally4_int_t     *ipiv;
    magma_tally4_int_t     N, n2, lda, lwork, info;
    magma_tally4_int_t     status = 0;
    magma_tally4_int_t     cpu = 0, nopiv = 0, nopiv_gpu = 0, row = 0;
    
    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );
    switch (opts.version) {
        case 1:
            cpu = 1;
            printf( "\nCPU-Interface to Bunch-Kauffman on GPU" );
            break;
        case 2:
            //gpu = 1;
            printf( "\nGPU-Interface to Bunch-Kauffman on GPU" );
            printf( "\n not yet..\n\n" );
            return 0;
            break;
        case 3:
            nopiv = 1;
            printf( "\nCPU-Interface to hybrid Non-pivoted LDLt (A is SPD)" );
            break;
        case 4:
            nopiv_gpu = 1;
            printf( "\nGPU-Interface to hybrid Non-pivoted LDLt (A is SPD)" );
            break;
            break;
        //case 5:
        //    row = 1;
        //    printf( "\n Bunch-Kauffman: GPU-only version (row-major)" );
        //    break;
        default:
        //  printf( " hybrid CPU-GPU version" );
            printf( " version = %d not supported\n\n", (int) opts.version );
            return 0;
    }

    printf( " (%s)\n", lapack_uplo_const_tally4(opts.uplo) );
    printf( " (--version: 1 = Bunch-Kauffman (CPU), 2 = Bunch-Kauffman (GPU), 3 = No-piv (CPU), 4 = No-piv (GPU))\n\n" );
    
    double tol = opts.tolerance * lapackf77_dlamch("E");

    if ( opts.check == 2 ) {
        printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |Ax-b|/(N*|A|*|x|)\n");
    }
    else {
        printf("    M     N   CPU GFlop/s (sec)   GPU GFlop/s (sec)   |PAP'-LDL'|/(N*|A|)\n");
    }
    printf("=========================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            n2     = lda*N;
            gflops = FLOPS_ZPOTRF( N ) / 1e9;
            
            TESTING_MALLOC_CPU( ipiv, magma_tally4_int_t, N );
            TESTING_MALLOC_PIN( h_A,  magma_tally4DoubleComplex, n2 );
            

            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                lwork = -1;
                lapackf77_zhetrf( lapack_uplo_const_tally4(opts.uplo), &N, h_A, &lda, ipiv, &temp, &lwork, &info);
                lwork = (int)MAGMA_tally4_Z_REAL(temp);
                TESTING_MALLOC_CPU( work, magma_tally4DoubleComplex, lwork );

                init_matrix( nopiv, N, N, h_A, lda );
                cpu_time = magma_tally4_wtime();
                lapackf77_zhetrf( lapack_uplo_const_tally4(opts.uplo), &N, h_A, &lda, ipiv, work, &lwork, &info);
                cpu_time = magma_tally4_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0)
                    printf("lapackf77_zhetrf returned error %d: %s.\n",
                           (int) info, magma_tally4_strerror( info ));
                error_lapack = get_residual( nopiv, opts.uplo, N, h_A, lda, ipiv );

                TESTING_FREE_CPU( work );
            }
           
            /* ====================================================================
               Performs operation using MAGMA_tally4
               =================================================================== */
            init_matrix( (nopiv | nopiv_gpu), N, N, h_A, lda );

            if (nopiv) {
                // CPU-interface to non-piv LDLt
                magma_tally4_setdevice(0);
                gpu_time = magma_tally4_wtime();
                magma_tally4_zhetrf_nopiv_tally4( opts.uplo, N, h_A, lda, &info);
                gpu_time = magma_tally4_wtime() - gpu_time;
            } else if (cpu) {
                // CPU-interface to Bunch-Kauffman LDLt
                magma_tally4_setdevice(0);
                gpu_time = magma_tally4_wtime();
                magma_tally4_zhetrf( opts.uplo, N, h_A, lda, ipiv, &info);
                gpu_time = magma_tally4_wtime() - gpu_time;
            } else if (nopiv_gpu) {
                // GPU-interface to non-piv LDLt
                magma_tally4_setdevice(0);
                magma_tally4_int_t ldda = 32*((N+31)/32);
                magma_tally4DoubleComplex_ptr d_A;
                if (MAGMA_tally4_SUCCESS != magma_tally4_zmalloc( &d_A, N*ldda  )) {
                    printf( " failed to allocate d_A(%dx%d)\n", (int) N, (int) ldda );
                    return 0;
                }
                magma_tally4_zsetmatrix(N, N, h_A, lda, d_A, ldda);
                gpu_time = magma_tally4_wtime();
                magma_tally4_zhetrf_nopiv_tally4_gpu( opts.uplo, N, d_A, ldda, &info);
                gpu_time = magma_tally4_wtime() - gpu_time;
                magma_tally4_zgetmatrix(N, N, d_A, ldda, h_A, lda);
                magma_tally4_free( d_A );
            } else if (row) {
                magma_tally4_setdevice(0);
                //magma_tally4_zhetrf_gpu_row( opts.uplo, N, h_A, lda, ipiv, work, lwork, &info);
            } else {
                magma_tally4_setdevice(0);
                //magma_tally4_zhetrf_hybrid( opts.uplo, N, h_A, lda, ipiv, work, lwork, &info);
            }
            gpu_perf = gflops / gpu_time;
            if (info != 0)
                printf("magma_tally4_zhetrf returned error %d: %s.\n",
                       (int) info, magma_tally4_strerror( info ));
            
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
            if ( opts.check == 2 ) {
                error = get_residual( (nopiv | nopiv_gpu), opts.uplo, N, h_A, lda, ipiv );
                printf("   %8.2e   %s", error, (error < tol ? "ok" : "failed"));
                if (opts.lapack)
                    printf(" (lapack rel.res. = %8.2e)", error_lapack);
                printf("\n");
                status += ! (error < tol);
            }
            else if ( opts.check ) {
                error = get_LDLt_error( (nopiv | nopiv_gpu), opts.uplo, N, h_A, lda, ipiv );
                printf("   %8.2e   %s\n", error, (error < tol ? "ok" : "failed"));
                status += ! (error < tol);
            }
            else {
                printf("     ---   \n");
            }
            
            TESTING_FREE_CPU( ipiv );
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