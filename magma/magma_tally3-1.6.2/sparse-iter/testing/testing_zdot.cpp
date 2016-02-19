/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s
       @author Hartwig Anzt
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_tally3.h"
#include "magma_tally3_lapack.h"
#include "testings.h"
#include "common_magma_tally3sparse.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- testing zdot
*/
int main(  int argc, char** argv )
{
    magma_tally3_int_t info = 0;
    // set queue for old dense routines
    magma_tally3_queue_t queue=NULL;
    magma_tally3_queue_create( /*devices[ opts->device ],*/ &queue );
    magma_tally3blasGetKernelStream( &queue );

    TESTING_INIT();


    magma_tally3_z_matrix a={Magma_tally3_CSR}, b={Magma_tally3_CSR}, x={Magma_tally3_CSR}, y={Magma_tally3_CSR}, skp={Magma_tally3_CSR};

        printf("#================================================================================================================================================\n");
        printf("\n");
        printf("            |                            runtime                             |                              GFLOPS\n");
        printf("#n num_vecs |  CUDOT       CUGEMV       MAGMA_tally3GEMV       MDOT       MDGM      |      CUDOT       CUGEMV      MAGMA_tally3GEMV       MDOT       MDGM      \n");
        printf("#------------------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("\n");

    for( magma_tally3_int_t num_vecs=5; num_vecs<6; num_vecs+=1 ) {
        for( magma_tally3_int_t n=10000; n<100000001; n=n+10000 ) {
            int iters = 10;
            double computations = (2.* n * iters * num_vecs);

            magma_tally3DoubleComplex one = MAGMA_tally3_Z_MAKE(1.0, 0.0);
            magma_tally3DoubleComplex zero = MAGMA_tally3_Z_MAKE(0.0, 0.0);
            magma_tally3DoubleComplex alpha;

            #define ENABLE_TIMER
            #ifdef ENABLE_TIMER
            real_Double_t mdot1, mdot2, mdgm1, mdgm2, magma_tally3gemv1, magma_tally3gemv2, cugemv1, cugemv2, cudot1, cudot2;
            real_Double_t mdot_time, mdgm_time, magma_tally3gemv_time, cugemv_time, cudot_time;
            #endif

            CHECK( magma_tally3_zvinit( &a, Magma_tally3_DEV, n, num_vecs, one, queue ));
            CHECK( magma_tally3_zvinit( &b, Magma_tally3_DEV, num_vecs, 1, one, queue ));
            int min_ten = min(num_vecs, 15);
            CHECK( magma_tally3_zvinit( &x, Magma_tally3_DEV, min_ten, n, one, queue ));
            CHECK( magma_tally3_zvinit( &y, Magma_tally3_DEV, min_ten, n, one, queue ));
            CHECK( magma_tally3_zvinit( &skp, Magma_tally3_DEV, num_vecs, 1, zero, queue ));

            // warm up
            CHECK( magma_tally3_zgemvmdot( n, num_vecs, a.dval, b.dval, x.dval, y.dval, skp.dval, queue ));

            // CUDOT
            #ifdef ENABLE_TIMER
            cudot1 = magma_tally3_sync_wtime( queue );
            #endif
            for( int h=0; h<iters; h++) {
                for( int l=0; l<num_vecs; l++)
                    alpha = magma_tally3_zdotc(n, a.dval, 1, b.dval, 1);
            }
            #ifdef ENABLE_TIMER
            cudot2 = magma_tally3_sync_wtime( queue );
            cudot_time=cudot2-cudot1;
            #endif
            // CUGeMV
            #ifdef ENABLE_TIMER
            cugemv1 = magma_tally3_sync_wtime( queue );
            #endif
            for( int h=0; h<iters; h++) {
                magma_tally3_zgemv(Magma_tally3Trans, n, num_vecs, one, a.dval, n, b.dval, 1, zero, skp.dval, 1);
                //h++;
            }
            #ifdef ENABLE_TIMER
            cugemv2 = magma_tally3_sync_wtime( queue );
            cugemv_time=cugemv2-cugemv1;
            #endif
            // MAGMA_tally3GeMV
            #ifdef ENABLE_TIMER
            magma_tally3gemv1 = magma_tally3_sync_wtime( queue );
            #endif
            for( int h=0; h<iters; h++) {
                magma_tally3blas_zgemv(Magma_tally3Trans, n, num_vecs, one, a.dval, n, b.dval, 1, zero, skp.dval, 1);
                //h++;
            }
            #ifdef ENABLE_TIMER
            magma_tally3gemv2 = magma_tally3_sync_wtime( queue );
            magma_tally3gemv_time=magma_tally3gemv2-magma_tally3gemv1;
            #endif
            // MDOT
            #ifdef ENABLE_TIMER
            mdot1 = magma_tally3_sync_wtime( queue );
            #endif
            for( int h=0; h<iters; h++) {
                //magma_tally3_zmdotc( n, num_vecs, a.dval, b.dval, x.dval, y.dval, skp.dval, queue );
                CHECK( magma_tally3_zmdotc( n, 2, a.dval, b.dval, x.dval, y.dval, skp.dval, queue ));
                CHECK( magma_tally3_zmdotc( n, 2, a.dval, b.dval, x.dval, y.dval, skp.dval, queue ));
                CHECK( magma_tally3_zmdotc( n, 1, a.dval, b.dval, x.dval, y.dval, skp.dval, queue ));
                //h++;
            }
            #ifdef ENABLE_TIMER
            mdot2 = magma_tally3_sync_wtime( queue );
            mdot_time=mdot2-mdot1;
            #endif
            // MDGM
            #ifdef ENABLE_TIMER
            mdgm1 = magma_tally3_sync_wtime( queue );
            #endif
            for( int h=0; h<iters; h++) {
                CHECK( magma_tally3_zgemvmdot( n, num_vecs, a.dval, b.dval, x.dval, y.dval, skp.dval, queue ));
                //h++;
            }
            #ifdef ENABLE_TIMER
            mdgm2 = magma_tally3_sync_wtime( queue );
            mdgm_time=mdgm2-mdgm1;
            #endif

            //magma_tally3_zprint_gpu(num_vecs,1,skp.dval,num_vecs);

            //Chronometry
            #ifdef ENABLE_TIMER
            printf("%d  %d  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n",
                    n, num_vecs,
                    cudot_time/iters,
                    (cugemv_time)/iters,
                    (magma_tally3gemv_time)/iters,
                    (mdot_time)/iters,
                    (mdgm_time)/iters,
                    (double)(computations)/(cudot_time*(1.e+09)),
                    (double)(computations)/(cugemv_time*(1.e+09)),
                    (double)(computations)/(magma_tally3gemv_time*(1.e+09)),
                    (double)(computations)/(mdot_time*(1.e+09)),
                    (double)(computations)/(mdgm_time*(1.e+09)) );
            #endif

            magma_tally3_zmfree(&a, queue );
            magma_tally3_zmfree(&b, queue );
            magma_tally3_zmfree(&x, queue );
            magma_tally3_zmfree(&y, queue );
            magma_tally3_zmfree(&skp, queue );
        }

        printf("#================================================================================================================================================\n");
        printf("\n");
        printf("\n");
    }

cleanup:
    magma_tally3_zmfree(&a, queue );
    magma_tally3_zmfree(&b, queue );
    magma_tally3_zmfree(&x, queue );
    magma_tally3_zmfree(&y, queue );
    magma_tally3_zmfree(&skp, queue );
    magma_tally3_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
