/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Azzam Haidar
       @author Mark Gates

       @generated from testing_zhegvd_m.cpp normal z -> c, Fri Jan 30 19:00:26 2015

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

#define COMPLEX

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing chegvd
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    magma_tally4FloatComplex *h_A, *h_Ainit, *h_B, *h_Binit, *h_work;
    #ifdef COMPLEX
    float *rwork;
    #endif
    float *w1, *w2, result[2]={0, 0};
    magma_tally4_int_t *iwork;
    real_Double_t mgpu_time, gpu_time, cpu_time;

    /* Matrix size */
    magma_tally4_int_t N, n2, nb;

    magma_tally4_int_t info;
    magma_tally4_int_t ione = 1;

    magma_tally4FloatComplex c_zero    = MAGMA_tally4_C_ZERO;
    magma_tally4FloatComplex c_one     = MAGMA_tally4_C_ONE;
    magma_tally4FloatComplex c_neg_one = MAGMA_tally4_C_NEG_ONE;

    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    magma_tally4_int_t status = 0;

    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );
    
    float tol    = opts.tolerance * lapackf77_slamch("E");
    float tolulp = opts.tolerance * lapackf77_slamch("P");

    // checking NoVec requires LAPACK
    opts.lapack |= (opts.check && opts.jobz == Magma_tally4NoVec);
    
    printf("using: ngpu = %d, itype = %d, jobz = %s, uplo = %s, check = %d\n",
           (int) opts.ngpu, (int) opts.itype,
           lapack_vec_const_tally4(opts.jobz), lapack_uplo_const_tally4(opts.uplo), (int) opts.check);

    printf("    N   CPU Time (sec)   GPU Time (sec)   MGPU Time (sec)\n");
    printf("=========================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            // TODO define lda
            N = opts.nsize[itest];
            n2     = N*N;
            nb     = magma_tally4_get_chetrd_nb(N);
            #ifdef COMPLEX
                magma_tally4_int_t lwork  = max( N + N*nb, 2*N + N*N );
                magma_tally4_int_t lrwork = 1 + 5*N +2*N*N;
            #else
                magma_tally4_int_t lwork  = max( 2*N + N*nb, 1 + 6*N + 2*N*N );
            #endif
            magma_tally4_int_t liwork = 3 + 5*N;

            TESTING_MALLOC_PIN( h_A,    magma_tally4FloatComplex, n2    );
            TESTING_MALLOC_PIN( h_B,    magma_tally4FloatComplex, n2    );
            TESTING_MALLOC_PIN( h_work, magma_tally4FloatComplex, lwork );
            #ifdef COMPLEX
            TESTING_MALLOC_PIN( rwork, float, lrwork );
            #endif

            TESTING_MALLOC_CPU( w1,    float, N );
            TESTING_MALLOC_CPU( w2,    float, N );
            TESTING_MALLOC_CPU( iwork, magma_tally4_int_t, liwork );

            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &n2, h_A );
            lapackf77_clarnv( &ione, ISEED, &n2, h_B );
            magma_tally4_cmake_hpd( N, h_B, N );
            magma_tally4_cmake_hermitian( N, h_A, N );

            if ( opts.warmup || opts.check ) {
                TESTING_MALLOC_CPU( h_Ainit, magma_tally4FloatComplex, n2 );
                TESTING_MALLOC_CPU( h_Binit, magma_tally4FloatComplex, n2 );
                lapackf77_clacpy( Magma_tally4FullStr, &N, &N, h_A, &N, h_Ainit, &N );
                lapackf77_clacpy( Magma_tally4FullStr, &N, &N, h_B, &N, h_Binit, &N );
            }

            if (opts.warmup) {
                // ==================================================================
                // Warmup using MAGMA_tally4.
                // ==================================================================
                magma_tally4_chegvd_m( opts.ngpu, opts.itype, opts.jobz, opts.uplo,
                                N, h_A, N, h_B, N, w1,
                                h_work, lwork,
                                #ifdef COMPLEX
                                rwork, lrwork,
                                #endif
                                iwork, liwork,
                                &info);
                lapackf77_clacpy( Magma_tally4FullStr, &N, &N, h_Ainit, &N, h_A, &N );
                lapackf77_clacpy( Magma_tally4FullStr, &N, &N, h_Binit, &N, h_B, &N );
            }

            // ===================================================================
            // Performs operation using MAGMA_tally4
            // ===================================================================
            mgpu_time = magma_tally4_wtime();
            magma_tally4_chegvd_m( opts.ngpu, opts.itype, opts.jobz, opts.uplo,
                            N, h_A, N, h_B, N, w1,
                            h_work, lwork,
                            #ifdef COMPLEX
                            rwork, lrwork,
                            #endif
                            iwork, liwork,
                            &info);
            mgpu_time = magma_tally4_wtime() - mgpu_time;

            if (info != 0)
                printf("magma_tally4_chegvd_m returned error %d: %s.\n",
                       (int) info, magma_tally4_strerror( info ));

            if ( opts.check && opts.jobz != Magma_tally4NoVec ) {
                /* =====================================================================
                   Check the results following the LAPACK's [zc]hegvd routine.
                   A x = lambda B x is solved
                   and the following 3 tests computed:
                   (1)    | A Z - B Z D | / ( |A||Z| N )  (itype = 1)
                          | A B Z - Z D | / ( |A||Z| N )  (itype = 2)
                          | B A Z - Z D | / ( |A||Z| N )  (itype = 3)
                   =================================================================== */

                #ifdef REAL
                float *rwork = h_work + N*N;
                #endif

                result[0] = 1.;
                result[0] /= lapackf77_clanhe("1", lapack_uplo_const_tally4(opts.uplo), &N, h_Ainit, &N, rwork);
                result[0] /= lapackf77_clange("1", &N, &N, h_A, &N, rwork);

                if (opts.itype == 1) {
                    blasf77_chemm("L", lapack_uplo_const_tally4(opts.uplo), &N, &N, &c_one, h_Ainit, &N, h_A, &N, &c_zero, h_work, &N);
                    for(int i=0; i < N; ++i)
                        blasf77_csscal(&N, &w1[i], &h_A[i*N], &ione);
                    blasf77_chemm("L", lapack_uplo_const_tally4(opts.uplo), &N, &N, &c_neg_one, h_Binit, &N, h_A, &N, &c_one, h_work, &N);
                    result[0] *= lapackf77_clange("1", &N, &N, h_work, &N, rwork)/N;
                }
                else if (opts.itype == 2) {
                    blasf77_chemm("L", lapack_uplo_const_tally4(opts.uplo), &N, &N, &c_one, h_Binit, &N, h_A, &N, &c_zero, h_work, &N);
                    for(int i=0; i < N; ++i)
                        blasf77_csscal(&N, &w1[i], &h_A[i*N], &ione);
                    blasf77_chemm("L", lapack_uplo_const_tally4(opts.uplo), &N, &N, &c_one, h_Ainit, &N, h_work, &N, &c_neg_one, h_A, &N);
                    result[0] *= lapackf77_clange("1", &N, &N, h_A, &N, rwork)/N;
                }
                else if (opts.itype == 3) {
                    blasf77_chemm("L", lapack_uplo_const_tally4(opts.uplo), &N, &N, &c_one, h_Ainit, &N, h_A, &N, &c_zero, h_work, &N);
                    for(int i=0; i < N; ++i)
                        blasf77_csscal(&N, &w1[i], &h_A[i*N], &ione);
                    blasf77_chemm("L", lapack_uplo_const_tally4(opts.uplo), &N, &N, &c_one, h_Binit, &N, h_work, &N, &c_neg_one, h_A, &N);
                    result[0] *= lapackf77_clange("1", &N, &N, h_A, &N, rwork)/N;
                }
            }
            if ( opts.lapack ) {
                lapackf77_clacpy( Magma_tally4FullStr, &N, &N, h_Ainit, &N, h_A, &N );
                lapackf77_clacpy( Magma_tally4FullStr, &N, &N, h_Binit, &N, h_B, &N );
                
                /* ====================================================================
                   Performs operation using MAGMA_tally4
                   =================================================================== */
                gpu_time = magma_tally4_wtime();
                magma_tally4_chegvd(opts.itype, opts.jobz, opts.uplo,
                             N, h_A, N, h_B, N, w2,
                             h_work, lwork,
                             #ifdef COMPLEX
                             rwork, lrwork,
                             #endif
                             iwork, liwork,
                             &info);
                gpu_time = magma_tally4_wtime() - gpu_time;

                if (info != 0)
                    printf("magma_tally4_chegvd returned error %d: %s.\n",
                           (int) info, magma_tally4_strerror( info ));

                /* =====================================================================
                   Performs operation using LAPACK
                   =================================================================== */
                cpu_time = magma_tally4_wtime();
                lapackf77_chegvd(&opts.itype, lapack_vec_const_tally4(opts.jobz), lapack_uplo_const_tally4(opts.uplo),
                                 &N, h_Ainit, &N, h_Binit, &N, w2,
                                 h_work, &lwork,
                                 #ifdef COMPLEX
                                 rwork, &lrwork,
                                 #endif
                                 iwork, &liwork,
                                 &info);
                cpu_time = magma_tally4_wtime() - cpu_time;
                if (info != 0)
                    printf("lapackf77_chegvd returned error %d: %s.\n",
                           (int) info, magma_tally4_strerror( info ));

                float maxw=0, diff=0;
                for(int j=0; j < N; j++) {
                    maxw = max(maxw, fabs(w1[j]));
                    maxw = max(maxw, fabs(w2[j]));
                    diff = max(diff, fabs(w1[j] - w2[j]));
                }
                result[1] = diff / (N*maxw);

                /* =====================================================================
                   Print execution time
                   =================================================================== */
                printf("%5d   %7.2f          %7.2f          %7.2f\n",
                       (int) N, cpu_time, gpu_time, mgpu_time);
            }
            else {
                printf("%5d     ---              ---            %7.2f\n",
                       (int) N, mgpu_time);
            }
            if ( opts.check && opts.jobz != Magma_tally4NoVec ) {
                printf("Testing the eigenvalues and eigenvectors for correctness:\n");
                if (opts.itype == 1) {
                    printf("    | A Z - B Z D | / (|A| |Z| N) = %8.2e   %s\n",   result[0],  (result[0] < tol ? "ok" : "failed") );
                }
                else if (opts.itype == 2) {
                    printf("    | A B Z - Z D | / (|A| |Z| N) = %8.2e   %s\n",   result[0],  (result[0] < tol ? "ok" : "failed") );
                }
                else if (opts.itype == 3) {
                    printf("    | B A Z - Z D | / (|A| |Z| N) = %8.2e   %s\n",   result[0],  (result[0] < tol ? "ok" : "failed") );
                }
                status += ! (result[0] < tol);
            }
            if ( opts.lapack ) {
                printf(    "    | D_mgpu - D_lapack | / |D|   = %8.2e   %s\n\n", result[1], (result[1] < tolulp ? "ok" : "failed") );
                status += ! (result[1] < tolulp);
            }

            /* Memory clean up */
            TESTING_FREE_PIN( h_A    );
            TESTING_FREE_PIN( h_B    );
            TESTING_FREE_PIN( h_work );
            #ifdef COMPLEX
            TESTING_FREE_PIN( rwork  );
            #endif
            
            TESTING_FREE_CPU( w1    );
            TESTING_FREE_CPU( w2    );
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

    TESTING_FINALIZE();
    return status;
}
