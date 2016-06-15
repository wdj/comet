/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Azzam Haidar
       @author Mark Gates

       @generated from testing_zhegvdx_2stage.cpp normal z -> s, Fri Jan 30 19:00:26 2015

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"
#include "testings.h"
#include "magma_tally4_sbulge.h"
#include "magma_tally4_threadsetting.h"

#define PRECISION_s


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing ssygvdx
*/
int main( int argc, char** argv)
{

    TESTING_INIT();

    real_Double_t gpu_time;

    float *h_A, *h_R, *h_B, *h_S, *h_work;

    #if defined(PRECISION_z) || defined(PRECISION_c)
    float *rwork;
    magma_tally4_int_t lrwork;
    #endif

    /* Matrix size */
    float *w1, *w2, result[2]={0,0};
    magma_tally4_int_t *iwork;
    magma_tally4_int_t N, n2, info, lwork, liwork;
    float c_zero    = MAGMA_tally4_S_ZERO;
    float c_one     = MAGMA_tally4_S_ONE;
    float c_neg_one = MAGMA_tally4_S_NEG_ONE;
    magma_tally4_int_t ione     = 1;
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    magma_tally4_int_t status = 0;

    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol    = opts.tolerance * lapackf77_slamch("E");
    float tolulp = opts.tolerance * lapackf77_slamch("P");

    magma_tally4_range_t range = Magma_tally4RangeAll;
    if (opts.fraction != 1)
        range = Magma_tally4RangeI;

    if ( opts.check && opts.jobz == Magma_tally4NoVec ) {
        fprintf( stderr, "checking results requires vectors; setting jobz=V (option -JV)\n" );
        opts.jobz = Magma_tally4Vec;
    }

    printf("using: itype = %d, jobz = %s, range = %s, uplo = %s, opts.check = %d, fraction = %6.4f\n",
           (int) opts.itype, lapack_vec_const_tally4(opts.jobz), lapack_range_const_tally4(range), lapack_uplo_const_tally4(opts.uplo),
           (int) opts.check, opts.fraction);

    printf("    N     M   GPU Time (sec)\n");
    printf("============================\n");
    magma_tally4_int_t threads = magma_tally4_get_parallel_numthreads();
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            n2     = N*N;
            #if defined(PRECISION_z) || defined(PRECISION_c)
            lwork  = magma_tally4_sbulge_get_lq2(N, threads) + 2*N + N*N;
            lrwork = 1 + 5*N +2*N*N;
            #else
            lwork  = magma_tally4_sbulge_get_lq2(N, threads) + 1 + 6*N + 2*N*N;
            #endif
            liwork = 3 + 5*N;

            /* Allocate host memory for the matrix */
            TESTING_MALLOC_CPU( h_A,    float, n2 );
            TESTING_MALLOC_CPU( h_B,    float, n2 );
            TESTING_MALLOC_CPU( w1,     float, N );
            TESTING_MALLOC_CPU( w2,     float, N );
            TESTING_MALLOC_CPU( iwork,  magma_tally4_int_t, liwork );
            
            TESTING_MALLOC_PIN( h_R,    float, n2 );
            TESTING_MALLOC_PIN( h_S,    float, n2 );
            TESTING_MALLOC_PIN( h_work, float, lwork );
            #if defined(PRECISION_z) || defined(PRECISION_c)
            TESTING_MALLOC_PIN( rwork,  float, lrwork);
            #endif

            /* Initialize the matrix */
            lapackf77_slarnv( &ione, ISEED, &n2, h_A );
            lapackf77_slarnv( &ione, ISEED, &n2, h_B );
            magma_tally4_smake_hpd( N, h_B, N );
            magma_tally4_smake_symmetric( N, h_A, N );

            magma_tally4_int_t m1 = 0;
            float vl = 0;
            float vu = 0;
            magma_tally4_int_t il = 0;
            magma_tally4_int_t iu = 0;

            if (range == Magma_tally4RangeI) {
                il = 1;
                iu = (int) (opts.fraction*N);
            }

            // ==================================================================
            // Warmup using MAGMA_tally4
            // ==================================================================
            if (opts.warmup) {
                lapackf77_slacpy( Magma_tally4UpperLowerStr, &N, &N, h_A, &N, h_R, &N );
                lapackf77_slacpy( Magma_tally4UpperLowerStr, &N, &N, h_B, &N, h_S, &N );

                magma_tally4_ssygvdx_2stage(opts.itype, opts.jobz, range, opts.uplo,
                                     N, h_R, N, h_S, N, vl, vu, il, iu, &m1, w1,
                                     h_work, lwork,
                                     #if defined(PRECISION_z) || defined(PRECISION_c)
                                     rwork, lrwork,
                                     #endif
                                     iwork, liwork,
                                     &info);
            }
            // ===================================================================
            // Performs operation using MAGMA_tally4
            // ===================================================================
            lapackf77_slacpy( Magma_tally4UpperLowerStr, &N, &N, h_A, &N, h_R, &N );
            lapackf77_slacpy( Magma_tally4UpperLowerStr, &N, &N, h_B, &N, h_S, &N );

            gpu_time = magma_tally4_wtime();
            magma_tally4_ssygvdx_2stage(opts.itype, opts.jobz, range, opts.uplo,
                                 N, h_R, N, h_S, N, vl, vu, il, iu, &m1, w1,
                                 h_work, lwork,
                                 #if defined(PRECISION_z) || defined(PRECISION_c)
                                 rwork, lrwork,
                                 #endif
                                 iwork, liwork,
                                 &info);
            gpu_time = magma_tally4_wtime() - gpu_time;


            if ( opts.check && opts.jobz != Magma_tally4NoVec ) {
                /* =====================================================================
                   Check the results following the LAPACK's [zc]hegvdx routine.
                   A x = lambda B x is solved
                   and the following 3 tests computed:
                   (1)    | A Z - B Z D | / ( |A||Z| N )  (itype = 1)
                          | A B Z - Z D | / ( |A||Z| N )  (itype = 2)
                          | B A Z - Z D | / ( |A||Z| N )  (itype = 3)
                   (2)    | S(with V) - S(w/o V) | / | S |
                   =================================================================== */
                #if defined(PRECISION_d) || defined(PRECISION_s)
                float *rwork = h_work + N*N;
                #endif

                result[0] = 1.;
                result[0] /= lapackf77_slansy("1", lapack_uplo_const_tally4(opts.uplo), &N, h_A, &N, rwork);
                result[0] /= lapackf77_slange("1", &N, &m1, h_R, &N, rwork);

                if (opts.itype == 1) {
                    blasf77_ssymm("L", lapack_uplo_const_tally4(opts.uplo), &N, &m1, &c_one, h_A, &N, h_R, &N, &c_zero, h_work, &N);
                    for(int i=0; i<m1; ++i)
                        blasf77_sscal(&N, &w1[i], &h_R[i*N], &ione);
                    blasf77_ssymm("L", lapack_uplo_const_tally4(opts.uplo), &N, &m1, &c_neg_one, h_B, &N, h_R, &N, &c_one, h_work, &N);
                    result[0] *= lapackf77_slange("1", &N, &m1, h_work, &N, rwork)/N;
                }
                else if (opts.itype == 2) {
                    blasf77_ssymm("L", lapack_uplo_const_tally4(opts.uplo), &N, &m1, &c_one, h_B, &N, h_R, &N, &c_zero, h_work, &N);
                    for(int i=0; i<m1; ++i)
                        blasf77_sscal(&N, &w1[i], &h_R[i*N], &ione);
                    blasf77_ssymm("L", lapack_uplo_const_tally4(opts.uplo), &N, &m1, &c_one, h_A, &N, h_work, &N, &c_neg_one, h_R, &N);
                    result[0] *= lapackf77_slange("1", &N, &m1, h_R, &N, rwork)/N;
                }
                else if (opts.itype == 3) {
                    blasf77_ssymm("L", lapack_uplo_const_tally4(opts.uplo), &N, &m1, &c_one, h_A, &N, h_R, &N, &c_zero, h_work, &N);
                    for(int i=0; i<m1; ++i)
                        blasf77_sscal(&N, &w1[i], &h_R[i*N], &ione);
                    blasf77_ssymm("L", lapack_uplo_const_tally4(opts.uplo), &N, &m1, &c_one, h_B, &N, h_work, &N, &c_neg_one, h_R, &N);
                    result[0] *= lapackf77_slange("1", &N, &m1, h_R, &N, rwork)/N;
                }

                lapackf77_slacpy( Magma_tally4UpperLowerStr, &N, &N, h_A, &N, h_R, &N );
                lapackf77_slacpy( Magma_tally4UpperLowerStr, &N, &N, h_B, &N, h_S, &N );

                magma_tally4_int_t m2 = m1;
                lapackf77_ssygvd(&opts.itype, "N", lapack_uplo_const_tally4(opts.uplo), &N,
                              h_R, &N, h_S, &N, w2,
                              h_work, &lwork,
                              #if defined(PRECISION_z) || defined(PRECISION_c)
                              rwork, &lrwork,
                              #endif
                              iwork, &liwork,
                              &info);

                float maxw=0, diff=0;
                for(int j=0; j<m2; j++) {
                    maxw = max(maxw, fabs(w1[j]));
                    maxw = max(maxw, fabs(w2[j]));
                    diff = max(diff, fabs(w1[j] - w2[j]));
                }
                result[1] = diff / (m2*maxw);
            }


            /* =====================================================================
               Print execution time
               =================================================================== */
            printf("%5d %5d   %7.2f\n",
                   (int) N, (int) m1, gpu_time);
            if ( opts.check && opts.jobz != Magma_tally4NoVec ) {
                printf("Testing the eigenvalues and eigenvectors for correctness:\n");
                if (opts.itype==1) {
                    printf("    | A Z - B Z D | / (|A| |Z| N) = %8.2e   %s\n",   result[0], (result[0] < tol    ? "ok" : "failed"));
                }
                else if (opts.itype==2) {
                    printf("    | A B Z - Z D | / (|A| |Z| N) = %8.2e   %s\n",   result[0], (result[0] < tol    ? "ok" : "failed"));
                }
                else if (opts.itype==3) {
                    printf("    | B A Z - Z D | / (|A| |Z| N) = %8.2e   %s\n",   result[0], (result[0] < tol    ? "ok" : "failed"));
                }
                printf(    "    | D(w/ Z) - D(w/o Z) | / |D|  = %8.2e   %s\n\n", result[1], (result[1] < tolulp ? "ok" : "failed"));
                status += ! (result[0] < tol && result[1] < tolulp);
            }

            TESTING_FREE_CPU( h_A   );
            TESTING_FREE_CPU( h_B   );
            TESTING_FREE_CPU( w1    );
            TESTING_FREE_CPU( w2    );
            TESTING_FREE_CPU( iwork );
            
            TESTING_FREE_PIN( h_R );
            TESTING_FREE_PIN( h_S );
            TESTING_FREE_PIN( h_work );
            #if defined(PRECISION_z) || defined(PRECISION_c)
            TESTING_FREE_PIN( rwork );
            #endif
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
