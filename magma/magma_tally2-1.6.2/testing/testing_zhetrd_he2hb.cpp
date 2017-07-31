/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Azzam Haidar

       @precisions normal z -> s d c

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
#include "magma_tally2_zbulge.h"
#include "testings.h"
#include "magma_tally2_threadsetting.h"

#if defined(USEMKL)
#include <mkl_service.h>
#endif

#if defined(USEACML)
#include <omp.h>
#endif

#define PRECISION_z


#if defined(PRECISION_z) || defined(PRECISION_d)
extern "C" void cmp_vals(int n, double *wr1, double *wr2, double *nrmI, double *nrm1, double *nrm2);
extern "C" void zcheck_eig_(char *JOBZ, int  *MATYPE, int  *N, int  *NB,
                       magma_tally2DoubleComplex* A, int  *LDA, double *AD, double *AE, double *D1, double *EIG,
                    magma_tally2DoubleComplex *Z, int  *LDZ, magma_tally2DoubleComplex *WORK, double *RWORK, double *RESU);
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zhetrd_he2hb
*/
int main( int argc, char** argv)
{
    TESTING_INIT();

    real_Double_t    gpu_time, gpu_perf, gflops;
    magma_tally2DoubleComplex *h_A, *h_R, *h_work, *dT1;
    magma_tally2DoubleComplex *tau;
    double *D, *E;

    /* Matrix size */
    magma_tally2_int_t N, n2, lda, lwork, lwork0;  //ldt

    magma_tally2_int_t info;
    magma_tally2_int_t ione     = 1;
    magma_tally2_int_t ISEED[4] = {0,0,0,1};

#if defined(CHECKEIG)
#if defined(PRECISION_z)  || defined(PRECISION_d)
    magma_tally2_int_t WANTZ=0;
    magma_tally2_int_t THREADS=1;
#endif
#endif

    magma_tally2_int_t NE = 0;
    magma_tally2_int_t NB = 0;
    magma_tally2_int_t ngpu = 1;
    
    magma_tally2_opts opts;
    parse_opts( argc, argv, &opts );
    
    NB = opts.nb;
    if (NB < 1)
        NB  = 64; //64; //magma_tally2_get_zhetrd_he2hb_nb(N);

    // what is NE ?
    if (NE < 1)
        NE  = 64; //N;  //magma_tally2_get_zhetrd_he2hb_nb(N);  // N not yet initialized

    printf("  N    GPU GFlop/s   \n");
    printf("=====================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda = N;
            //ldt = N;
            n2  = N*lda;
            gflops = FLOPS_ZHETRD( N ) / 1e9;
            
            /* We suppose the magma_tally2 NB is bigger than lapack NB */
            lwork0 = N*NB;
        
            /* Allocate host memory for the matrix */
            TESTING_MALLOC_CPU( h_A,    magma_tally2DoubleComplex, lda*N  );
            TESTING_MALLOC_CPU( tau,    magma_tally2DoubleComplex, N-1    );
            
            TESTING_MALLOC_PIN( h_R,    magma_tally2DoubleComplex, lda*N  );
            TESTING_MALLOC_PIN( h_work, magma_tally2DoubleComplex, lwork0 );
            TESTING_MALLOC_PIN( D, double, N );
            TESTING_MALLOC_PIN( E, double, N );
            
            //TESTING_MALLOC_DEV( dT1, magma_tally2DoubleComplex, (2*min(N,N)+(N+31)/32*32)*NB );
            TESTING_MALLOC_DEV( dT1, magma_tally2DoubleComplex, (N*NB) );
        
            // if (WANTZ) gflops = 2.0*gflops;
    
            /* ====================================================================
               Initialize the matrix
               =================================================================== */
            lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
            magma_tally2_zmake_hermitian( N, h_A, lda );
            
            lapackf77_zlacpy( Magma_tally2UpperLowerStr, &N, &N, h_A, &lda, h_R, &lda );
    
    
            /* ====================================================================
               Performs operation using MAGMA_tally2
               =================================================================== */
            magma_tally2_device_t cdev;
            magma_tally2_getdevice( &cdev );
    
            gpu_time = magma_tally2_wtime();
            /*
            magma_tally2_zhetrd_he2hb( opts.uplo, N, NB, h_R, lda, tau, h_work, lwork0, dT1, THREADS, &info);
            tband = magma_tally2_wtime - gpu_time();
            printf("  Finish BAND  N %d  NB %d  ngpu %d timing= %f\n", N, NB, ngpu, tband);
            magma_tally2_zhetrd_bhe2trc_v5(THREADS, WANTZ, opts.uplo, NE, N, NB, h_R, lda, D, E, dT1, ldt);
            */
    
            /*
            magma_tally2_zhetrd_he2hb( opts.uplo, N, NB, h_R, lda, tau, h_work, lwork, dT1, THREADS, &info);
            tband = magma_tally2_wtime - gpu_time();
            printf("  Finish BAND  N %d  NB %d  ngpu %d timing= %f\n", N, NB, ngpu, tband);
            magma_tally2_zhetrd_bhe2trc(THREADS, WANTZ, opts.uplo, NE, N, NB, h_R, lda, D, E, dT1, ldt);
            */

            magma_tally2_range_t range = Magma_tally2RangeAll;
            magma_tally2_int_t fraction_ev = 100;
            magma_tally2_int_t il, iu, m1;
            double vl=0., vu=0.;
    
            if (fraction_ev == 0) {
                il = N / 10;
                iu = N / 5+il;
            }
            else {
                il = 1;
                iu = (int)(fraction_ev*N);
                if (iu < 1) iu = 1;
            }
            magma_tally2DoubleComplex *hh_work;
            magma_tally2_int_t *iwork;
            magma_tally2_int_t /*nb,*/ /*lwork,*/ liwork;
            magma_tally2_int_t threads = magma_tally2_get_parallel_numthreads();
            #if defined(PRECISION_z) || defined(PRECISION_c)
                double *rwork;
                magma_tally2_int_t lrwork;
                lwork  = magma_tally2_zbulge_get_lq2(N, threads) + 2*N + N*N;
                lrwork = 1 + 5*N +2*N*N;
                TESTING_MALLOC_PIN( rwork, double, lrwork );
            #else
                lwork  = magma_tally2_zbulge_get_lq2(N, threads) + 1 + 6*N + 2*N*N;
            #endif
            liwork = 3 + 5*N;
            //nb = magma_tally2_get_zhetrd_nb(N);
            TESTING_MALLOC_PIN( hh_work, magma_tally2DoubleComplex, lwork  );
            TESTING_MALLOC_CPU( iwork,   magma_tally2_int_t,        liwork );
    
            if (ngpu == 1) {
                printf("calling zheevdx_2stage 1 GPU\n");
                magma_tally2_zheevdx_2stage( opts.jobz, range, opts.uplo, N,
                                h_R, lda,
                                vl, vu, il, iu,
                                &m1, D,
                                hh_work, lwork,
                                #if defined(PRECISION_z) || defined(PRECISION_c)
                                rwork, lrwork,
                                #endif
                                iwork, liwork,
                                &info);
    
            } else {
                printf("calling zheevdx_2stage_m %d GPU\n", (int) ngpu);
                magma_tally2_zheevdx_2stage_m(ngpu, opts.jobz, range, opts.uplo, N,
                                h_R, lda,
                                vl, vu, il, iu,
                                &m1, D,
                                hh_work, lwork,
                                #if defined(PRECISION_z) || defined(PRECISION_c)
                                rwork, lrwork,
                                #endif
                                iwork, liwork,
                                &info);
            }
    
    
            magma_tally2_setdevice( cdev );
            gpu_time = magma_tally2_wtime() - gpu_time;
            gpu_perf = gflops / gpu_time;
    
            /* =====================================================================
               Check the factorization
               =================================================================== */
            /*
            if ( opts.check ) {
                FILE        *fp ;
    
                printf("Writing input matrix in matlab_i_mat.txt ...\n");
                fp = fopen ("matlab_i_mat.txt", "w") ;
                if ( fp == NULL ) { printf("Couldn't open output file\n"); exit(1); }
    
                for (j=0; j < N; j++) {
                    for (k=0; k < N; k++) {
                        #if defined(PRECISION_z) || defined(PRECISION_c)
                        fprintf(fp, "%5d %5d %11.8f %11.8f\n", k+1, j+1,
                                h_A[k+j*lda].x, h_A[k+j*lda].y);
                        #else
                        fprintf(fp, "%5d %5d %11.8f\n", k+1, j+1, h_A[k+j*lda]);
                        #endif
                    }
                }
                fclose( fp ) ;
    
                printf("Writing output matrix in matlab_o_mat.txt ...\n");
                fp = fopen ("matlab_o_mat.txt", "w") ;
                if ( fp == NULL ) { printf("Couldn't open output file\n"); exit(1); }
    
                for (j=0; j < N; j++) {
                    for (k=0; k < N; k++) {
                        #if defined(PRECISION_z) || defined(PRECISION_c)
                        fprintf(fp, "%5d %5d %11.8f %11.8f\n", k+1, j+1,
                                h_R[k+j*lda].x, h_R[k+j*lda].y);
                        #else
                        fprintf(fp, "%5d %5d %11.8f\n", k+1, j+1, h_R[k+j*lda]);
                        #endif
                    }
                }
                fclose( fp ) ;
            }
            */
    
    
    
            /* =====================================================================
               Print performance and error.
               =================================================================== */
#if defined(CHECKEIG)
#if defined(PRECISION_z)  || defined(PRECISION_d)
            if ( opts.check ) {
                printf("  Total N %5d  gflops %6.2f  timing %6.2f seconds\n", (int) N, gpu_perf, gpu_time );
                char JOBZ;
                if (WANTZ == 0)
                    JOBZ = 'N';
                else
                    JOBZ = 'V';
                double nrmI=0.0, nrm1=0.0, nrm2=0.0;
                int    lwork2 = 256*N;
                magma_tally2DoubleComplex *work2, *AINIT;
                double *rwork2, *D2;
                // TODO free this memory !
                magma_tally2_zmalloc_cpu( &work2, lwork2 );
                magma_tally2_dmalloc_cpu( &rwork2, N );
                magma_tally2_dmalloc_cpu( &D2, N );
                magma_tally2_zmalloc_cpu( &AINIT, N*lda );
                memcpy(AINIT, h_A, N*lda*sizeof(magma_tally2DoubleComplex));
                /* compute the eigenvalues using lapack routine to be able to compare to it and used as ref */
                cpu_time = magma_tally2_wtime();
                i= min(12, THREADS);
    
                #if defined(USEMKL)
                mkl_set_num_threads( i );
                #endif
                #if defined(USEACML)
                omp_set_num_threads(i);
                #endif
    
                lapackf77_zheev( "N", "L", &N, h_A, &lda, D2, work2, &lwork2,
                    #if defined(PRECISION_z) || defined (PRECISION_c)
                    rwork2,
                    #endif
                    &info );
                
                ///* call eigensolver for our resulting tridiag [D E] and for Q */
                //dstedc_withZ('V', N, D, E, h_R, lda);
                ////dsterf_( &N, D, E, &info);
                ////
                cpu_time = magma_tally2_wtime() - cpu_time;
                printf("  Finish CHECK - EIGEN   timing= %f  threads %d\n", cpu_time, i);
    
                /*
                for (i=0; i < 10; i++)
                    printf(" voici lpk D[%d] %8.2e\n", i, D2[i]);
                */
    
                //magma_tally2DoubleComplex mydz=0.0, mydo=1.0;
                //magma_tally2DoubleComplex *Z;
                // magma_tally2_zmalloc_cpu( &Z, N*lda );
                // dgemm_("N", "N", &N, &N, &N, &mydo, h_R, &lda, h_A, &lda, &mydz, Z, &lda);
    
    
                /* compare result */
                cmp_vals(N, D2, D, &nrmI, &nrm1, &nrm2);
    
    
                magma_tally2DoubleComplex *WORKAJETER;
                double *RWORKAJETER, *RESU;
                // TODO free this memory !
                magma_tally2_zmalloc_cpu( &WORKAJETER, (2* N * N + N)  );
                magma_tally2_dmalloc_cpu( &RWORKAJETER, N  );
                magma_tally2_dmalloc_cpu( &RESU, 10 );
                int MATYPE;
                memset(RESU, 0, 10*sizeof(double));
    
     
                MATYPE=3;
                double NOTHING=0.0;
                cpu_time = magma_tally2_wtime();
                // check results
                zcheck_eig_(&JOBZ, &MATYPE, &N, &NB, AINIT, &lda, &NOTHING, &NOTHING, D2, D, h_R, &lda, WORKAJETER, RWORKAJETER, RESU );
                cpu_time = magma_tally2_wtime() - cpu_time;
                printf("  Finish CHECK - results timing= %f\n", cpu_time);
                #if defined(USEMKL)
                mkl_set_num_threads( 1 );
                #endif
                #if defined(USEACML)
                omp_set_num_threads(1);
                #endif
    
                printf("\n");
                printf(" ================================================================================================================\n");
                printf("   ==> INFO voici  threads=%d    N=%d    NB=%d   WANTZ=%d\n", (int) THREADS, (int) N, (int) NB, (int) WANTZ);
                printf(" ================================================================================================================\n");
                printf("            DSBTRD                : %15s \n", "STATblgv9withQ    ");
                printf(" ================================================================================================================\n");
                if (WANTZ > 0)
                    printf(" | A - U S U' | / ( |A| n ulp )   : %15.3E   \n", RESU[0]);
                if (WANTZ > 0)
                    printf(" | I - U U' | / ( n ulp )         : %15.3E   \n", RESU[1]);
                printf(" | D1 - EVEIGS | / (|D| ulp)      : %15.3E   \n",  RESU[2]);
                printf(" max | D1 - EVEIGS |              : %15.3E   \n",  RESU[6]);
                printf(" ================================================================================================================\n\n\n");
                
                printf(" ****************************************************************************************************************\n");
                printf(" * Hello here are the norm  Infinite (max)=%8.2e  norm one (sum)=%8.2e   norm2(sqrt)=%8.2e *\n", nrmI, nrm1, nrm2);
                printf(" ****************************************************************************************************************\n\n");
            }
#endif
#endif
            
            printf("  Total N %5d  gflops %6.2f        timing %6.2f seconds\n", (int) N, gpu_perf, gpu_time );
            printf("============================================================================\n\n\n");
            
            /* Memory clean up */
            TESTING_FREE_CPU( h_A );
            TESTING_FREE_CPU( tau );
            
            TESTING_FREE_PIN( h_R    );
            TESTING_FREE_PIN( h_work );
            TESTING_FREE_PIN( D      );
            TESTING_FREE_PIN( E      );
            
            TESTING_FREE_DEV( dT1 );
            
            /* TODO - not all memory has been freed inside loop */
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    TESTING_FINALIZE();
    return EXIT_SUCCESS;
}
