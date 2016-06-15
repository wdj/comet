/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Azzam Haidar
       @author Mark Gates

       @generated from testing_zhegvdx.cpp normal z -> c, Fri Jan 30 19:00:26 2015

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
   -- Testing chegvdx
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t   gpu_time /*cpu_time*/;
    magma_tally4FloatComplex *h_A, *h_R, *h_B, *h_S, *h_work;
    float *w1, *w2, vl=0, vu=0;
    float result[2] = {0};
    magma_tally4_int_t *iwork;
    magma_tally4_int_t N, n2, info, il, iu, m1, nb, lwork, liwork;
    magma_tally4FloatComplex c_zero    = MAGMA_tally4_C_ZERO;
    magma_tally4FloatComplex c_one     = MAGMA_tally4_C_ONE;
    magma_tally4FloatComplex c_neg_one = MAGMA_tally4_C_NEG_ONE;
    #ifdef COMPLEX
    float *rwork;
    magma_tally4_int_t lrwork;
    #endif
    //float d_one         =  1.;
    //float d_ten         = 10.;
    magma_tally4_int_t ione     = 1;
    magma_tally4_int_t ISEED[4] = {0,0,0,1};
    magma_tally4_int_t status = 0;

    magma_tally4_opts opts;
    parse_opts( argc, argv, &opts );

    float tol    = opts.tolerance * lapackf77_slamch("E");
    //float tolulp = opts.tolerance * lapackf77_slamch("P");
    
    // TODO: how to check NoVec? This doesn't have an equivalent call to LAPACK.
    if ( opts.check && opts.jobz == Magma_tally4NoVec ) {
        fprintf( stderr, "WARNING: cannot currently check results when not computing vectors (option -JN).\n" );
    }
    
    printf("using: itype = %d, jobz = %s, uplo = %s, check = %d, fraction = %6.4f\n",
           (int) opts.itype, lapack_vec_const_tally4(opts.jobz), lapack_uplo_const_tally4(opts.uplo),
           (int) opts.check, opts.fraction);

    printf("    N     M   GPU Time (sec)\n");
    printf("============================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            n2     = N*N;
            nb     = magma_tally4_get_chetrd_nb(N);
            #ifdef COMPLEX
                lwork  = max( N + N*nb, 2*N + N*N );
                lrwork = 1 + 5*N +2*N*N;
            #else
                lwork  = max( 2*N + N*nb, 1 + 6*N + 2*N*N );
            #endif
            liwork = 3 + 5*N;

            if ( opts.fraction == 0 ) {
                il = N / 10;
                iu = N / 5+il;
            }
            else {
                il = 1;
                iu = (int) (opts.fraction*N);
                if (iu < 1) iu = 1;
            }

            TESTING_MALLOC_CPU( h_A,    magma_tally4FloatComplex, n2     );
            TESTING_MALLOC_CPU( h_B,    magma_tally4FloatComplex, n2     );
            TESTING_MALLOC_CPU( w1,     float,             N      );
            TESTING_MALLOC_CPU( w2,     float,             N      );
            TESTING_MALLOC_CPU( iwork,  magma_tally4_int_t,        liwork );
            
            TESTING_MALLOC_PIN( h_R,    magma_tally4FloatComplex, n2     );
            TESTING_MALLOC_PIN( h_S,    magma_tally4FloatComplex, n2     );
            TESTING_MALLOC_PIN( h_work, magma_tally4FloatComplex, lwork  );
            #ifdef COMPLEX
            TESTING_MALLOC_PIN( rwork, float, lrwork);
            #endif
            
            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &n2, h_A );
            lapackf77_clarnv( &ione, ISEED, &n2, h_B );
            magma_tally4_cmake_hpd( N, h_B, N );
            magma_tally4_cmake_hermitian( N, h_A, N );

            // ==================================================================
            // Warmup using MAGMA_tally4
            // ==================================================================
            if ( opts.warmup ) {
                lapackf77_clacpy( Magma_tally4FullStr, &N, &N, h_A, &N, h_R, &N );
                lapackf77_clacpy( Magma_tally4FullStr, &N, &N, h_B, &N, h_S, &N );
                
                magma_tally4_chegvdx( opts.itype, opts.jobz, Magma_tally4RangeI, opts.uplo,
                               N, h_R, N, h_S, N, vl, vu, il, iu, &m1, w1,
                               h_work, lwork,
                               #ifdef COMPLEX
                               rwork, lrwork,
                               #endif      
                               iwork, liwork,
                               &info );
                if (info != 0)
                    printf("magma_tally4_chegvdx returned error %d: %s.\n",
                           (int) info, magma_tally4_strerror( info ));
            }
            /* ====================================================================
               Performs operation using MAGMA_tally4
               =================================================================== */
            lapackf77_clacpy( Magma_tally4FullStr, &N, &N, h_A, &N, h_R, &N );
            lapackf77_clacpy( Magma_tally4FullStr, &N, &N, h_B, &N, h_S, &N );

            gpu_time = magma_tally4_wtime();
            magma_tally4_chegvdx( opts.itype, opts.jobz, Magma_tally4RangeI, opts.uplo,
                           N, h_R, N, h_S, N, vl, vu, il, iu, &m1, w1,
                           h_work, lwork,
                           #ifdef COMPLEX
                           rwork, lrwork,
                           #endif
                           iwork, liwork,
                           &info );
            gpu_time = magma_tally4_wtime() - gpu_time;
            if (info != 0)
                printf("magma_tally4_chegvdx returned error %d: %s.\n",
                       (int) info, magma_tally4_strerror( info ));
            
            if ( opts.check && opts.jobz != Magma_tally4NoVec ) {
                /* =====================================================================
                   Check the results following the LAPACK's [zc]hegvdx routine.
                   A x = lambda B x is solved
                   and the following 3 tests computed:
                   (1)    | A Z - B Z D | / ( |A||Z| N )  (itype = 1)
                          | A B Z - Z D | / ( |A||Z| N )  (itype = 2)
                          | B A Z - Z D | / ( |A||Z| N )  (itype = 3)
                   (2)    | D(with V) - D(w/o V) | / | D |
                   =================================================================== */
                #ifdef REAL
                float *rwork = h_work + N*N;
                #endif
                
                result[0] = 1.;
                result[0] /= lapackf77_clanhe("1", lapack_uplo_const_tally4(opts.uplo), &N, h_A, &N, rwork);
                result[0] /= lapackf77_clange("1", &N, &m1, h_R, &N, rwork);
                
                if (opts.itype == 1) {
                    blasf77_chemm("L", lapack_uplo_const_tally4(opts.uplo), &N, &m1, &c_one, h_A, &N, h_R, &N, &c_zero, h_work, &N);
                    for(int i=0; i < m1; ++i)
                        blasf77_csscal(&N, &w1[i], &h_R[i*N], &ione);
                    blasf77_chemm("L", lapack_uplo_const_tally4(opts.uplo), &N, &m1, &c_neg_one, h_B, &N, h_R, &N, &c_one, h_work, &N);
                    result[0] *= lapackf77_clange("1", &N, &m1, h_work, &N, rwork)/N;
                }
                else if (opts.itype == 2) {
                    blasf77_chemm("L", lapack_uplo_const_tally4(opts.uplo), &N, &m1, &c_one, h_B, &N, h_R, &N, &c_zero, h_work, &N);
                    for(int i=0; i < m1; ++i)
                        blasf77_csscal(&N, &w1[i], &h_R[i*N], &ione);
                    blasf77_chemm("L", lapack_uplo_const_tally4(opts.uplo), &N, &m1, &c_one, h_A, &N, h_work, &N, &c_neg_one, h_R, &N);
                    result[0] *= lapackf77_clange("1", &N, &m1, h_R, &N, rwork)/N;
                }
                else if (opts.itype == 3) {
                    blasf77_chemm("L", lapack_uplo_const_tally4(opts.uplo), &N, &m1, &c_one, h_A, &N, h_R, &N, &c_zero, h_work, &N);
                    for(int i=0; i < m1; ++i)
                        blasf77_csscal(&N, &w1[i], &h_R[i*N], &ione);
                    blasf77_chemm("L", lapack_uplo_const_tally4(opts.uplo), &N, &m1, &c_one, h_B, &N, h_work, &N, &c_neg_one, h_R, &N);
                    result[0] *= lapackf77_clange("1", &N, &m1, h_R, &N, rwork)/N;
                }
                
                // Disable eigenvalue check which calls routine again --
                // it obscures whether error occurs in first call above or in this call.
                // TODO: add comparison to LAPACK, as in testing_chegvd.cpp?
                //
                //lapackf77_clacpy( Magma_tally4FullStr, &N, &N, h_A, &N, h_R, &N );
                //lapackf77_clacpy( Magma_tally4FullStr, &N, &N, h_B, &N, h_S, &N );
                //
                //magma_tally4_int_t m2;
                //magma_tally4_chegvdx( opts.itype, Magma_tally4NoVec, Magma_tally4RangeI, opts.uplo,
                //               N, h_R, N, h_S, N, vl, vu, il, iu, &m2, w2,
                //               h_work, lwork,
                //               #ifdef COMPLEX
                //               rwork, lrwork,
                //               #endif
                //               iwork, liwork,
                //               &info );
                //if (info != 0)
                //    printf("magma_tally4_chegvdx returned error %d: %s.\n",
                //           (int) info, magma_tally4_strerror( info ));
                //
                //float maxw=0, diff=0;
                //for(int j=0; j < m2; j++) {
                //    maxw = max(maxw, fabs(w1[j]));
                //    maxw = max(maxw, fabs(w2[j]));
                //    diff = max(diff, fabs(w1[j]-w2[j]));
                //}
                //result[1] = diff / (m2*maxw);
            }
            
            /* =====================================================================
               Print execution time
               =================================================================== */
            printf("%5d %5d   %7.2f\n",
                   (int) N, (int) m1, gpu_time);
            if ( opts.check && opts.jobz != Magma_tally4NoVec ) {
                printf("Testing the eigenvalues and eigenvectors for correctness:\n");
                if (opts.itype == 1) {
                    printf("    | A Z - B Z D | / (|A| |Z| N) = %8.2e   %s\n",   result[0], (result[0] < tol    ? "ok" : "failed"));
                }
                else if (opts.itype == 2) {
                    printf("    | A B Z - Z D | / (|A| |Z| N) = %8.2e   %s\n",   result[0], (result[0] < tol    ? "ok" : "failed"));
                }
                else if (opts.itype == 3) {
                    printf("    | B A Z - Z D | / (|A| |Z| N) = %8.2e   %s\n",   result[0], (result[0] < tol    ? "ok" : "failed"));
                }
                //printf(    "    | D(w/ Z) - D(w/o Z) | / |D|  = %8.2e   %s\n\n", result[1], (result[1] < tolulp ? "ok" : "failed"));
                status += ! (result[0] < tol);  // && result[1] < tolulp);
            }
            
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( h_B );
            TESTING_FREE_CPU( w1  );
            TESTING_FREE_CPU( w2  );
            TESTING_FREE_CPU( iwork );
            
            TESTING_FREE_PIN( h_R    );
            TESTING_FREE_PIN( h_S    );
            TESTING_FREE_PIN( h_work );
            #ifdef COMPLEX
            TESTING_FREE_PIN( rwork );
            #endif
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return status;
}
