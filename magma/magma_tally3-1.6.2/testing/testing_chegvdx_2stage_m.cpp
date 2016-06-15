/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Azzam Haidar
       @author Mark Gates

       @generated from testing_zhegvdx_2stage_m.cpp normal z -> c, Fri Jan 30 19:00:26 2015

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "testings.h"
#include "magma_tally3_cbulge.h"
#include "magma_tally3_threadsetting.h"

#define PRECISION_c

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing chegvdx
*/
int main( int argc, char** argv)
{

    TESTING_INIT();

    real_Double_t   mgpu_time;
    magma_tally3FloatComplex *h_A, *h_Ainit, *h_B, *h_Binit, *h_work;

#if defined(PRECISION_z) || defined(PRECISION_c)
    float *rwork;
    magma_tally3_int_t lrwork;
#endif

    float *w1, result=0;
    magma_tally3_int_t *iwork;
    magma_tally3_int_t N, n2, info, lwork, liwork;
    magma_tally3FloatComplex c_zero    = MAGMA_tally3_C_ZERO;
    magma_tally3FloatComplex c_one     = MAGMA_tally3_C_ONE;
    magma_tally3FloatComplex c_neg_one = MAGMA_tally3_C_NEG_ONE;
    magma_tally3_int_t ione     = 1;
    magma_tally3_int_t ISEED[4] = {0,0,0,1};
    magma_tally3_int_t status = 0;

    magma_tally3_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol = opts.tolerance * lapackf77_slamch("E");

    magma_tally3_range_t range = Magma_tally3RangeAll;
    if (opts.fraction != 1)
        range = Magma_tally3RangeI;

    if ( opts.check && opts.jobz == Magma_tally3NoVec ) {
        fprintf( stderr, "checking results requires vectors; setting jobz=V (option -JV)\n" );
        opts.jobz = Magma_tally3Vec;
    }

    printf("using: ngpu = %d, itype = %d, jobz = %s, range = %s, uplo = %s, opts.check = %d, fraction = %6.4f\n",
           (int) opts.ngpu, (int) opts.itype,
           lapack_vec_const_tally3(opts.jobz), lapack_range_const_tally3(range), lapack_uplo_const_tally3(opts.uplo),
           (int) opts.check, opts.fraction);
    
    printf("    N     M   ngpu   MGPU Time (sec)\n");
    printf("====================================\n");
    magma_tally3_int_t threads = magma_tally3_get_parallel_numthreads();
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            n2     = N*N;
            #if defined(PRECISION_z) || defined(PRECISION_c)
            lwork  = magma_tally3_cbulge_get_lq2(N, threads) + 2*N + N*N;
            lrwork = 1 + 5*N +2*N*N;
            #else
            lwork  = magma_tally3_cbulge_get_lq2(N, threads) + 1 + 6*N + 2*N*N;
            #endif
            liwork = 3 + 5*N;


            //magma_tally3_int_t NB = 96;//magma_tally3_bulge_get_nb(N);
            //magma_tally3_int_t sizvblg = magma_tally3_cbulge_get_lq2(N, threads);        
            //magma_tally3_int_t siz = max(sizvblg,n2)+2*(N*NB+N)+24*N; 
            /* Allocate host memory for the matrix */
            TESTING_MALLOC_PIN( h_A,    magma_tally3FloatComplex, n2 );
            TESTING_MALLOC_PIN( h_B,    magma_tally3FloatComplex, n2 );
            TESTING_MALLOC_PIN( h_work, magma_tally3FloatComplex, lwork );
            #if defined(PRECISION_z) || defined(PRECISION_c)
            TESTING_MALLOC_PIN( rwork,  float, lrwork);
            #endif

            TESTING_MALLOC_CPU( w1,     float, N );
            TESTING_MALLOC_CPU( iwork,  magma_tally3_int_t, liwork);
            
            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &n2, h_A );
            lapackf77_clarnv( &ione, ISEED, &n2, h_B );
            magma_tally3_cmake_hpd( N, h_B, N );
            magma_tally3_cmake_hermitian( N, h_A, N );

            if ( opts.warmup || opts.check ) {
                TESTING_MALLOC_CPU( h_Ainit, magma_tally3FloatComplex, n2 );
                TESTING_MALLOC_CPU( h_Binit, magma_tally3FloatComplex, n2 );
                lapackf77_clacpy( Magma_tally3UpperLowerStr, &N, &N, h_A, &N, h_Ainit, &N );
                lapackf77_clacpy( Magma_tally3UpperLowerStr, &N, &N, h_B, &N, h_Binit, &N );
            }



            magma_tally3_int_t m1 = 0;
            float vl = 0;
            float vu = 0;
            magma_tally3_int_t il = 0;
            magma_tally3_int_t iu = 0;

            if (range == Magma_tally3RangeI) {
                il = 1;
                iu = (int) (opts.fraction*N);
            }

            if ( opts.warmup ) {

                // ==================================================================
                // Warmup using MAGMA_tally3. I prefer to use smalltest to warmup A-
                // ==================================================================
                magma_tally3_chegvdx_2stage_m(opts.ngpu, opts.itype, opts.jobz, range, opts.uplo,
                                       N, h_A, N, h_B, N, vl, vu, il, iu, &m1, w1,
                                       h_work, lwork,
                                       #if defined(PRECISION_z) || defined(PRECISION_c)
                                       rwork, lrwork,
                                       #endif
                                       iwork, liwork,
                                       &info);
                lapackf77_clacpy( Magma_tally3UpperLowerStr, &N, &N, h_Ainit, &N, h_A, &N );
                lapackf77_clacpy( Magma_tally3UpperLowerStr, &N, &N, h_Binit, &N, h_B, &N );
            }

            // ===================================================================
            // Performs operation using MAGMA_tally3
            // ===================================================================

            mgpu_time = magma_tally3_wtime();
            magma_tally3_chegvdx_2stage_m(opts.ngpu, opts.itype, opts.jobz, range, opts.uplo,
                                   N, h_A, N, h_B, N, vl, vu, il, iu, &m1, w1,
                                   h_work, lwork,
                                       #if defined(PRECISION_z) || defined(PRECISION_c)
                                   rwork, lrwork,
                                       #endif
                                   iwork, liwork,
                                   &info);
            mgpu_time = magma_tally3_wtime() - mgpu_time;

            if ( opts.check && opts.jobz != Magma_tally3NoVec ) {
                // ===================================================================
                // Check the results following the LAPACK's [zc]hegvdx routine.
                // A x = lambda B x is solved
                // and the following 3 tests computed:
                // (1)    | A Z - B Z D | / ( |A||Z| N )  (itype = 1)
                //        | A B Z - Z D | / ( |A||Z| N )  (itype = 2)
                //        | B A Z - Z D | / ( |A||Z| N )  (itype = 3)
                // ===================================================================
                #if defined(PRECISION_d) || defined(PRECISION_s)
                float *rwork = h_work + N*N;
                #endif
                result = 1.;
                result /= lapackf77_clanhe("1", lapack_uplo_const_tally3(opts.uplo), &N, h_Ainit, &N, rwork);
                result /= lapackf77_clange("1", &N , &m1, h_A, &N, rwork);

                if (opts.itype == 1) {
                    blasf77_chemm("L", lapack_uplo_const_tally3(opts.uplo), &N, &m1, &c_one, h_Ainit, &N, h_A, &N, &c_zero, h_work, &N);
                    for(int i=0; i<m1; ++i)
                        blasf77_csscal(&N, &w1[i], &h_A[i*N], &ione);
                    blasf77_chemm("L", lapack_uplo_const_tally3(opts.uplo), &N, &m1, &c_neg_one, h_Binit, &N, h_A, &N, &c_one, h_work, &N);
                    result *= lapackf77_clange("1", &N, &m1, h_work, &N, rwork)/N;
                }
                else if (opts.itype == 2) {
                    blasf77_chemm("L", lapack_uplo_const_tally3(opts.uplo), &N, &m1, &c_one, h_Binit, &N, h_A, &N, &c_zero, h_work, &N);
                    for(int i=0; i<m1; ++i)
                        blasf77_csscal(&N, &w1[i], &h_A[i*N], &ione);
                    blasf77_chemm("L", lapack_uplo_const_tally3(opts.uplo), &N, &m1, &c_one, h_Ainit, &N, h_work, &N, &c_neg_one, h_A, &N);
                    result *= lapackf77_clange("1", &N, &m1, h_A, &N, rwork)/N;
                }
                else if (opts.itype == 3) {
                    blasf77_chemm("L", lapack_uplo_const_tally3(opts.uplo), &N, &m1, &c_one, h_Ainit, &N, h_A, &N, &c_zero, h_work, &N);
                    for(int i=0; i<m1; ++i)
                        blasf77_csscal(&N, &w1[i], &h_A[i*N], &ione);
                    blasf77_chemm("L", lapack_uplo_const_tally3(opts.uplo), &N, &m1, &c_one, h_Binit, &N, h_work, &N, &c_neg_one, h_A, &N);
                    result *= lapackf77_clange("1", &N, &m1, h_A, &N, rwork)/N;
                }
            }

            // ===================================================================
            // Print execution time
            // ===================================================================
            printf("%5d %5d   %4d   %7.2f\n",
                   (int) N, (int) m1, (int) opts.ngpu, mgpu_time);
            if ( opts.check && opts.jobz != Magma_tally3NoVec ) {
                printf("Testing the eigenvalues and eigenvectors for correctness:\n");
                if (opts.itype==1) {
                    printf("    | A Z - B Z D | / (|A| |Z| N) = %8.2e   %s\n", result, (result < tol ? "ok" : "failed") );
                }
                else if (opts.itype==2) {
                    printf("    | A B Z - Z D | / (|A| |Z| N) = %8.2e   %s\n", result, (result < tol ? "ok" : "failed") );
                }
                else if (opts.itype==3) {
                    printf("    | B A Z - Z D | / (|A| |Z| N) = %8.2e   %s\n", result, (result < tol ? "ok" : "failed") );
                }
                printf("\n");
                status += ! (result < tol);
            }

            TESTING_FREE_PIN( h_A    );
            TESTING_FREE_PIN( h_B    );
            TESTING_FREE_PIN( h_work );
            #if defined(PRECISION_z) || defined(PRECISION_c)
            TESTING_FREE_PIN( rwork  );
            #endif

            TESTING_FREE_CPU( w1    );
            TESTING_FREE_CPU( iwork );
            if ( opts.warmup || opts.check ) {
                TESTING_FREE_CPU( h_Ainit );
                TESTING_FREE_CPU( h_Binit );
            }
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    /* Shutdown */
    TESTING_FINALIZE();
    return status;
}
