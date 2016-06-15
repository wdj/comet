/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from testing_zswap.cpp normal z -> c, Fri Jan 30 19:00:25 2015
       @author Mark Gates
*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// includes, project
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "testings.h"

#define PRECISION_c

// if ( A==B ) return 0, else return 1
static int diff_matrix( magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3FloatComplex *A, magma_tally3_int_t lda, magma_tally3FloatComplex *B, magma_tally3_int_t ldb )
{
    for( magma_tally3_int_t j = 0; j < n; j++ ) {
        for( magma_tally3_int_t i = 0; i < m; i++ ) {
            if ( ! MAGMA_tally3_C_EQUAL( A[lda*j+i], B[ldb*j+i] ) )
                return 1;
        }
    }
    return 0;
}

// fill matrix with entries Aij = offset + (i+1) + (j+1)/10000,
// which makes it easy to identify which rows & cols have been swapped.
static void init_matrix( magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3FloatComplex *A, magma_tally3_int_t lda, magma_tally3_int_t offset )
{
    assert( lda >= m );
    for( magma_tally3_int_t j = 0; j < n; ++j ) {
        for( magma_tally3_int_t i=0; i < m; ++i ) {
            A[i + j*lda] = MAGMA_tally3_C_MAKE( offset + (i+1) + (j+1)/10000., 0 );
        }
    }
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cswap, cswapblk, claswp, claswpx
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    magma_tally3FloatComplex *h_A1, *h_A2;
    magma_tally3FloatComplex *h_R1, *h_R2;
    magma_tally3FloatComplex_ptr d_A1, d_A2;
    
    // row-major and column-major performance
    real_Double_t row_perf0 = MAGMA_tally3_D_NAN, col_perf0 = MAGMA_tally3_D_NAN;
    real_Double_t row_perf1 = MAGMA_tally3_D_NAN, col_perf1 = MAGMA_tally3_D_NAN;
    real_Double_t row_perf2 = MAGMA_tally3_D_NAN, col_perf2 = MAGMA_tally3_D_NAN;
    real_Double_t row_perf4 = MAGMA_tally3_D_NAN;
    real_Double_t row_perf5 = MAGMA_tally3_D_NAN, col_perf5 = MAGMA_tally3_D_NAN;
    real_Double_t row_perf6 = MAGMA_tally3_D_NAN, col_perf6 = MAGMA_tally3_D_NAN;
    real_Double_t row_perf7 = MAGMA_tally3_D_NAN;
    real_Double_t cpu_perf  = MAGMA_tally3_D_NAN;

    real_Double_t time, gbytes;

    magma_tally3_int_t N, lda, ldda, nb, j;
    magma_tally3_int_t ione = 1;
    magma_tally3_int_t *ipiv, *ipiv2;
    magma_tally3Int_ptr d_ipiv;
    magma_tally3_int_t status = 0;
    
    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );

    magma_tally3_queue_t queue = 0;
    
    printf("            %8s cswap    cswap             cswapblk          claswp   claswp2  claswpx           ccopymatrix      CPU      (all in )\n", g_platform_str );
    printf("    N   nb  row-maj/col-maj   row-maj/col-maj   row-maj/col-maj   row-maj  row-maj  row-maj/col-maj   row-blk/col-blk  claswp   (GByte/s)\n");
    printf("=========================================================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            // For an N x N matrix, swap nb rows or nb columns using various methods.
            // Each test is assigned one bit in the 'check' bitmask; bit=1 indicates failure.
            // The variable 'shift' keeps track of which bit is for current test
            int shift = 1;
            int check = 0;
            N = opts.nsize[itest];
            lda    = N;
            ldda   = ((N+31)/32)*32;
            nb     = (opts.nb > 0 ? opts.nb : magma_tally3_get_cgetrf_nb( N ));
            nb     = min( N, nb );
            // each swap does 2N loads and 2N stores, for nb swaps
            gbytes = sizeof(magma_tally3FloatComplex) * 4.*N*nb / 1e9;
            
            TESTING_MALLOC_PIN( h_A1, magma_tally3FloatComplex, lda*N );
            TESTING_MALLOC_PIN( h_A2, magma_tally3FloatComplex, lda*N );
            TESTING_MALLOC_PIN( h_R1, magma_tally3FloatComplex, lda*N );
            TESTING_MALLOC_PIN( h_R2, magma_tally3FloatComplex, lda*N );
            
            TESTING_MALLOC_CPU( ipiv,  magma_tally3_int_t, nb );
            TESTING_MALLOC_CPU( ipiv2, magma_tally3_int_t, nb );
            
            TESTING_MALLOC_DEV( d_ipiv, magma_tally3_int_t, nb );
            TESTING_MALLOC_DEV( d_A1, magma_tally3FloatComplex, ldda*N );
            TESTING_MALLOC_DEV( d_A2, magma_tally3FloatComplex, ldda*N );
            
            // getrf always makes ipiv[j] >= j+1, where ipiv is one based and j is zero based
            // some implementations (e.g., MacOS dlaswp) assume this
            for( j=0; j < nb; j++ ) {
                ipiv[j] = (rand() % (N-j)) + j + 1;
                assert( ipiv[j] >= j+1 );
                assert( ipiv[j] <= N   );
            }
            
            /* =====================================================================
             * cublas / clBLAS / Xeon Phi cswap, row-by-row (2 matrices)
             */
            
            /* Row Major */
            init_matrix( N, N, h_A1, lda, 0 );
            init_matrix( N, N, h_A2, lda, 100 );
            magma_tally3_csetmatrix( N, N, h_A1, lda, d_A1, ldda );
            magma_tally3_csetmatrix( N, N, h_A2, lda, d_A2, ldda );
            
            time = magma_tally3_sync_wtime( queue );
            for( j=0; j < nb; j++) {
                if ( j != (ipiv[j]-1)) {
                    #ifdef HAVE_CUBLAS
                        cublasCswap( opts.handle, N, d_A1+ldda*j, 1, d_A2+ldda*(ipiv[j]-1), 1 );
                    #else
                        magma_tally3_cswap( N, d_A1, ldda*j, 1, d_A2, ldda*(ipiv[j]-1), 1, opts.queue );
                    #endif
                }
            }
            time = magma_tally3_sync_wtime( queue ) - time;
            row_perf0 = gbytes / time;
            
            for( j=0; j < nb; j++) {
                if ( j != (ipiv[j]-1)) {
                    blasf77_cswap( &N, h_A1+lda*j, &ione, h_A2+lda*(ipiv[j]-1), &ione);
                }
            }
            magma_tally3_cgetmatrix( N, N, d_A1, ldda, h_R1, lda );
            magma_tally3_cgetmatrix( N, N, d_A2, ldda, h_R2, lda );
            check += (diff_matrix( N, N, h_A1, lda, h_R1, lda ) ||
                      diff_matrix( N, N, h_A2, lda, h_R2, lda ))*shift;
            shift *= 2;
            
            /* Column Major */
            init_matrix( N, N, h_A1, lda, 0 );
            init_matrix( N, N, h_A2, lda, 100 );
            magma_tally3_csetmatrix( N, N, h_A1, lda, d_A1, ldda );
            magma_tally3_csetmatrix( N, N, h_A2, lda, d_A2, ldda );
            
            time = magma_tally3_sync_wtime( queue );
            for( j=0; j < nb; j++) {
                if ( j != (ipiv[j]-1)) {
                    #ifdef HAVE_CUBLAS
                        cublasCswap( opts.handle, N, d_A1+j, ldda, d_A2+ipiv[j]-1, ldda );
                    #else
                        magma_tally3_cswap( N, d_A1, j, ldda, d_A2, ipiv[j]-1, ldda, opts.queue );
                    #endif
                }
            }
            time = magma_tally3_sync_wtime( queue ) - time;
            col_perf0 = gbytes / time;
            
            for( j=0; j < nb; j++) {
                if ( j != (ipiv[j]-1)) {
                    blasf77_cswap( &N, h_A1+j, &lda, h_A2+(ipiv[j]-1), &lda);
                }
            }
            magma_tally3_cgetmatrix( N, N, d_A1, ldda, h_R1, lda );
            magma_tally3_cgetmatrix( N, N, d_A2, ldda, h_R2, lda );
            check += (diff_matrix( N, N, h_A1, lda, h_R1, lda ) ||
                      diff_matrix( N, N, h_A2, lda, h_R2, lda ))*shift;
            shift *= 2;

            /* =====================================================================
             * cswap, row-by-row (2 matrices)
             */
            
            /* Row Major */
            init_matrix( N, N, h_A1, lda, 0 );
            init_matrix( N, N, h_A2, lda, 100 );
            magma_tally3_csetmatrix( N, N, h_A1, lda, d_A1, ldda );
            magma_tally3_csetmatrix( N, N, h_A2, lda, d_A2, ldda );
            
            time = magma_tally3_sync_wtime( queue );
            for( j=0; j < nb; j++) {
                if ( j != (ipiv[j]-1)) {
                    magma_tally3blas_cswap( N, d_A1+ldda*j, 1, d_A2+ldda*(ipiv[j]-1), 1);
                }
            }
            time = magma_tally3_sync_wtime( queue ) - time;
            row_perf1 = gbytes / time;
            
            for( j=0; j < nb; j++) {
                if ( j != (ipiv[j]-1)) {
                    blasf77_cswap( &N, h_A1+lda*j, &ione, h_A2+lda*(ipiv[j]-1), &ione);
                }
            }
            magma_tally3_cgetmatrix( N, N, d_A1, ldda, h_R1, lda );
            magma_tally3_cgetmatrix( N, N, d_A2, ldda, h_R2, lda );
            check += (diff_matrix( N, N, h_A1, lda, h_R1, lda ) ||
                      diff_matrix( N, N, h_A2, lda, h_R2, lda ))*shift;
            shift *= 2;
            
            /* Column Major */
            init_matrix( N, N, h_A1, lda, 0 );
            init_matrix( N, N, h_A2, lda, 100 );
            magma_tally3_csetmatrix( N, N, h_A1, lda, d_A1, ldda );
            magma_tally3_csetmatrix( N, N, h_A2, lda, d_A2, ldda );
            
            time = magma_tally3_sync_wtime( queue );
            for( j=0; j < nb; j++) {
                if ( j != (ipiv[j]-1)) {
                    magma_tally3blas_cswap( N, d_A1+j, ldda, d_A2+ipiv[j]-1, ldda );
                }
            }
            time = magma_tally3_sync_wtime( queue ) - time;
            col_perf1 = gbytes / time;
            
            for( j=0; j < nb; j++) {
                if ( j != (ipiv[j]-1)) {
                    blasf77_cswap( &N, h_A1+j, &lda, h_A2+(ipiv[j]-1), &lda);
                }
            }
            magma_tally3_cgetmatrix( N, N, d_A1, ldda, h_R1, lda );
            magma_tally3_cgetmatrix( N, N, d_A2, ldda, h_R2, lda );
            check += (diff_matrix( N, N, h_A1, lda, h_R1, lda ) ||
                      diff_matrix( N, N, h_A2, lda, h_R2, lda ))*shift;
            shift *= 2;

            /* =====================================================================
             * cswapblk, blocked version (2 matrices)
             */
            
            #ifdef HAVE_CUBLAS
            /* Row Major */
            init_matrix( N, N, h_A1, lda, 0 );
            init_matrix( N, N, h_A2, lda, 100 );
            magma_tally3_csetmatrix( N, N, h_A1, lda, d_A1, ldda );
            magma_tally3_csetmatrix( N, N, h_A2, lda, d_A2, ldda );
            
            time = magma_tally3_sync_wtime( queue );
            magma_tally3blas_cswapblk( Magma_tally3RowMajor, N, d_A1, ldda, d_A2, ldda, 1, nb, ipiv, 1, 0);
            time = magma_tally3_sync_wtime( queue ) - time;
            row_perf2 = gbytes / time;
            
            for( j=0; j < nb; j++) {
                if ( j != (ipiv[j]-1)) {
                    blasf77_cswap( &N, h_A1+lda*j, &ione, h_A2+lda*(ipiv[j]-1), &ione);
                }
            }
            magma_tally3_cgetmatrix( N, N, d_A1, ldda, h_R1, lda );
            magma_tally3_cgetmatrix( N, N, d_A2, ldda, h_R2, lda );
            check += (diff_matrix( N, N, h_A1, lda, h_R1, lda ) ||
                      diff_matrix( N, N, h_A2, lda, h_R2, lda ))*shift;
            shift *= 2;
            
            /* Column Major */
            init_matrix( N, N, h_A1, lda, 0 );
            init_matrix( N, N, h_A2, lda, 100 );
            magma_tally3_csetmatrix( N, N, h_A1, lda, d_A1, ldda );
            magma_tally3_csetmatrix( N, N, h_A2, lda, d_A2, ldda );
            
            time = magma_tally3_sync_wtime( queue );
            magma_tally3blas_cswapblk( Magma_tally3ColMajor, N, d_A1, ldda, d_A2, ldda, 1, nb, ipiv, 1, 0);
            time = magma_tally3_sync_wtime( queue ) - time;
            col_perf2 = gbytes / time;
            
            for( j=0; j < nb; j++) {
                if ( j != (ipiv[j]-1)) {
                    blasf77_cswap( &N, h_A1+j, &lda, h_A2+(ipiv[j]-1), &lda);
                }
            }
            magma_tally3_cgetmatrix( N, N, d_A1, ldda, h_R1, lda );
            magma_tally3_cgetmatrix( N, N, d_A2, ldda, h_R2, lda );
            check += (diff_matrix( N, N, h_A1, lda, h_R1, lda ) ||
                      diff_matrix( N, N, h_A2, lda, h_R2, lda ))*shift;
            shift *= 2;
            #endif

            /* =====================================================================
             * LAPACK-style claswp (1 matrix)
             */
            
            #ifdef HAVE_CUBLAS
            /* Row Major */
            init_matrix( N, N, h_A1, lda, 0 );
            magma_tally3_csetmatrix( N, N, h_A1, lda, d_A1, ldda );
            
            time = magma_tally3_sync_wtime( queue );
            magma_tally3blas_claswp( N, d_A1, ldda, 1, nb, ipiv, 1);
            time = magma_tally3_sync_wtime( queue ) - time;
            row_perf4 = gbytes / time;
            
            for( j=0; j < nb; j++) {
                if ( j != (ipiv[j]-1)) {
                    blasf77_cswap( &N, h_A1+lda*j, &ione, h_A1+lda*(ipiv[j]-1), &ione);
                }
            }
            magma_tally3_cgetmatrix( N, N, d_A1, ldda, h_R1, lda );
            check += diff_matrix( N, N, h_A1, lda, h_R1, lda )*shift;
            shift *= 2;
            #endif

            /* =====================================================================
             * LAPACK-style claswp (1 matrix) - d_ipiv on GPU
             */
            
            #ifdef HAVE_CUBLAS
            /* Row Major */
            init_matrix( N, N, h_A1, lda, 0 );
            magma_tally3_csetmatrix( N, N, h_A1, lda, d_A1, ldda );
            
            time = magma_tally3_sync_wtime( queue );
            magma_tally3_setvector( nb, sizeof(magma_tally3_int_t), ipiv, 1, d_ipiv, 1 );
            magma_tally3blas_claswp2( N, d_A1, ldda, 1, nb, d_ipiv, 1 );
            time = magma_tally3_sync_wtime( queue ) - time;
            row_perf7 = gbytes / time;
            
            for( j=0; j < nb; j++) {
                if ( j != (ipiv[j]-1)) {
                    blasf77_cswap( &N, h_A1+lda*j, &ione, h_A1+lda*(ipiv[j]-1), &ione);
                }
            }
            magma_tally3_cgetmatrix( N, N, d_A1, ldda, h_R1, lda );
            check += diff_matrix( N, N, h_A1, lda, h_R1, lda )*shift;
            shift *= 2;
            #endif

            /* =====================================================================
             * LAPACK-style claswpx (extended for row- and col-major) (1 matrix)
             */
            
            #ifdef HAVE_CUBLAS
            /* Row Major */
            init_matrix( N, N, h_A1, lda, 0 );
            magma_tally3_csetmatrix( N, N, h_A1, lda, d_A1, ldda );
            
            time = magma_tally3_sync_wtime( queue );
            magma_tally3blas_claswpx( N, d_A1, ldda, 1, 1, nb, ipiv, 1);
            time = magma_tally3_sync_wtime( queue ) - time;
            row_perf5 = gbytes / time;
            
            for( j=0; j < nb; j++) {
                if ( j != (ipiv[j]-1)) {
                    blasf77_cswap( &N, h_A1+lda*j, &ione, h_A1+lda*(ipiv[j]-1), &ione);
                }
            }
            magma_tally3_cgetmatrix( N, N, d_A1, ldda, h_R1, lda );
            check += diff_matrix( N, N, h_A1, lda, h_R1, lda )*shift;
            shift *= 2;
            
            /* Col Major */
            init_matrix( N, N, h_A1, lda, 0 );
            magma_tally3_csetmatrix( N, N, h_A1, lda, d_A1, ldda );
            
            time = magma_tally3_sync_wtime( queue );
            magma_tally3blas_claswpx( N, d_A1, 1, ldda, 1, nb, ipiv, 1);
            time = magma_tally3_sync_wtime( queue ) - time;
            col_perf5 = gbytes / time;
            #endif
            
            /* LAPACK swap on CPU for comparison */
            time = magma_tally3_wtime();
            lapackf77_claswp( &N, h_A1, &lda, &ione, &nb, ipiv, &ione);
            time = magma_tally3_wtime() - time;
            cpu_perf = gbytes / time;
            
            #ifdef HAVE_CUBLAS
            magma_tally3_cgetmatrix( N, N, d_A1, ldda, h_R1, lda );
            check += diff_matrix( N, N, h_A1, lda, h_R1, lda )*shift;
            shift *= 2;
            #endif

            /* =====================================================================
             * Copy matrix.
             */
            
            time = magma_tally3_sync_wtime( queue );
            magma_tally3_ccopymatrix( N, nb, d_A1, ldda, d_A2, ldda );
            time = magma_tally3_sync_wtime( queue ) - time;
            // copy reads 1 matrix and writes 1 matrix, so has half gbytes of swap
            col_perf6 = 0.5 * gbytes / time;
            
            time = magma_tally3_sync_wtime( queue );
            magma_tally3_ccopymatrix( nb, N, d_A1, ldda, d_A2, ldda );
            time = magma_tally3_sync_wtime( queue ) - time;
            // copy reads 1 matrix and writes 1 matrix, so has half gbytes of swap
            row_perf6 = 0.5 * gbytes / time;

            printf("%5d  %3d  %6.2f%c/ %6.2f%c  %6.2f%c/ %6.2f%c  %6.2f%c/ %6.2f%c  %6.2f%c  %6.2f%c  %6.2f%c/ %6.2f%c  %6.2f / %6.2f  %6.2f  %10s\n",
                   (int) N, (int) nb,
                   row_perf0, ((check & 0x001) != 0 ? '*' : ' '),
                   col_perf0, ((check & 0x002) != 0 ? '*' : ' '),
                   row_perf1, ((check & 0x004) != 0 ? '*' : ' '),
                   col_perf1, ((check & 0x008) != 0 ? '*' : ' '),
                   row_perf2, ((check & 0x010) != 0 ? '*' : ' '),
                   col_perf2, ((check & 0x020) != 0 ? '*' : ' '),
                   row_perf4, ((check & 0x040) != 0 ? '*' : ' '),
                   row_perf7, ((check & 0x080) != 0 ? '*' : ' '),
                   row_perf5, ((check & 0x100) != 0 ? '*' : ' '),
                   col_perf5, ((check & 0x200) != 0 ? '*' : ' '),
                   row_perf6,
                   col_perf6,
                   cpu_perf,
                   (check == 0 ? "ok" : "* failed") );
            status += ! (check == 0);
            
            TESTING_FREE_PIN( h_A1 );
            TESTING_FREE_PIN( h_A2 );
            TESTING_FREE_PIN( h_R1 );
            TESTING_FREE_PIN( h_R2 );
            
            TESTING_FREE_CPU( ipiv  );
            TESTING_FREE_CPU( ipiv2 );
            
            TESTING_FREE_DEV( d_ipiv );
            TESTING_FREE_DEV( d_A1 );
            TESTING_FREE_DEV( d_A2 );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }
    
    TESTING_FINALIZE();
    return status;
}
