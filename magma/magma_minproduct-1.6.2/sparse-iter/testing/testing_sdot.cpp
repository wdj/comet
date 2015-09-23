/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from testing_zdot.cpp normal z -> s, Sun May  3 11:23:02 2015
       @author Hartwig Anzt
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "testings.h"
#include "common_magma_minproductsparse.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- testing zdot
*/
int main(  int argc, char** argv )
{
    magma_minproduct_int_t info = 0;
    // set queue for old dense routines
    magma_minproduct_queue_t queue=NULL;
    magma_minproduct_queue_create( /*devices[ opts->device ],*/ &queue );
    magma_minproductblasGetKernelStream( &queue );

    TESTING_INIT();


    magma_minproduct_s_matrix a={Magma_minproduct_CSR}, b={Magma_minproduct_CSR}, x={Magma_minproduct_CSR}, y={Magma_minproduct_CSR}, skp={Magma_minproduct_CSR};

        printf("#================================================================================================================================================\n");
        printf("\n");
        printf("            |                            runtime                             |                              GFLOPS\n");
        printf("#n num_vecs |  CUDOT       CUGEMV       MAGMA_minproductGEMV       MDOT       MDGM      |      CUDOT       CUGEMV      MAGMA_minproductGEMV       MDOT       MDGM      \n");
        printf("#------------------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("\n");

    for( magma_minproduct_int_t num_vecs=5; num_vecs<6; num_vecs+=1 ) {
        for( magma_minproduct_int_t n=10000; n<100000001; n=n+10000 ) {
            int iters = 10;
            float computations = (2.* n * iters * num_vecs);

            float one = MAGMA_minproduct_S_MAKE(1.0, 0.0);
            float zero = MAGMA_minproduct_S_MAKE(0.0, 0.0);
            float alpha;

            #define ENABLE_TIMER
            #ifdef ENABLE_TIMER
            real_Double_t mdot1, mdot2, mdgm1, mdgm2, magma_minproductgemv1, magma_minproductgemv2, cugemv1, cugemv2, cudot1, cudot2;
            real_Double_t mdot_time, mdgm_time, magma_minproductgemv_time, cugemv_time, cudot_time;
            #endif

            CHECK( magma_minproduct_svinit( &a, Magma_minproduct_DEV, n, num_vecs, one, queue ));
            CHECK( magma_minproduct_svinit( &b, Magma_minproduct_DEV, num_vecs, 1, one, queue ));
            int min_ten = min(num_vecs, 15);
            CHECK( magma_minproduct_svinit( &x, Magma_minproduct_DEV, min_ten, n, one, queue ));
            CHECK( magma_minproduct_svinit( &y, Magma_minproduct_DEV, min_ten, n, one, queue ));
            CHECK( magma_minproduct_svinit( &skp, Magma_minproduct_DEV, num_vecs, 1, zero, queue ));

            // warm up
            CHECK( magma_minproduct_sgemvmdot( n, num_vecs, a.dval, b.dval, x.dval, y.dval, skp.dval, queue ));

            // CUDOT
            #ifdef ENABLE_TIMER
            cudot1 = magma_minproduct_sync_wtime( queue );
            #endif
            for( int h=0; h<iters; h++) {
                for( int l=0; l<num_vecs; l++)
                    alpha = magma_minproduct_sdot(n, a.dval, 1, b.dval, 1);
            }
            #ifdef ENABLE_TIMER
            cudot2 = magma_minproduct_sync_wtime( queue );
            cudot_time=cudot2-cudot1;
            #endif
            // CUGeMV
            #ifdef ENABLE_TIMER
            cugemv1 = magma_minproduct_sync_wtime( queue );
            #endif
            for( int h=0; h<iters; h++) {
                magma_minproduct_sgemv(Magma_minproductTrans, n, num_vecs, one, a.dval, n, b.dval, 1, zero, skp.dval, 1);
                //h++;
            }
            #ifdef ENABLE_TIMER
            cugemv2 = magma_minproduct_sync_wtime( queue );
            cugemv_time=cugemv2-cugemv1;
            #endif
            // MAGMA_minproductGeMV
            #ifdef ENABLE_TIMER
            magma_minproductgemv1 = magma_minproduct_sync_wtime( queue );
            #endif
            for( int h=0; h<iters; h++) {
                magma_minproductblas_sgemv(Magma_minproductTrans, n, num_vecs, one, a.dval, n, b.dval, 1, zero, skp.dval, 1);
                //h++;
            }
            #ifdef ENABLE_TIMER
            magma_minproductgemv2 = magma_minproduct_sync_wtime( queue );
            magma_minproductgemv_time=magma_minproductgemv2-magma_minproductgemv1;
            #endif
            // MDOT
            #ifdef ENABLE_TIMER
            mdot1 = magma_minproduct_sync_wtime( queue );
            #endif
            for( int h=0; h<iters; h++) {
                //magma_minproduct_smdotc( n, num_vecs, a.dval, b.dval, x.dval, y.dval, skp.dval, queue );
                CHECK( magma_minproduct_smdotc( n, 2, a.dval, b.dval, x.dval, y.dval, skp.dval, queue ));
                CHECK( magma_minproduct_smdotc( n, 2, a.dval, b.dval, x.dval, y.dval, skp.dval, queue ));
                CHECK( magma_minproduct_smdotc( n, 1, a.dval, b.dval, x.dval, y.dval, skp.dval, queue ));
                //h++;
            }
            #ifdef ENABLE_TIMER
            mdot2 = magma_minproduct_sync_wtime( queue );
            mdot_time=mdot2-mdot1;
            #endif
            // MDGM
            #ifdef ENABLE_TIMER
            mdgm1 = magma_minproduct_sync_wtime( queue );
            #endif
            for( int h=0; h<iters; h++) {
                CHECK( magma_minproduct_sgemvmdot( n, num_vecs, a.dval, b.dval, x.dval, y.dval, skp.dval, queue ));
                //h++;
            }
            #ifdef ENABLE_TIMER
            mdgm2 = magma_minproduct_sync_wtime( queue );
            mdgm_time=mdgm2-mdgm1;
            #endif

            //magma_minproduct_sprint_gpu(num_vecs,1,skp.dval,num_vecs);

            //Chronometry
            #ifdef ENABLE_TIMER
            printf("%d  %d  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n",
                    n, num_vecs,
                    cudot_time/iters,
                    (cugemv_time)/iters,
                    (magma_minproductgemv_time)/iters,
                    (mdot_time)/iters,
                    (mdgm_time)/iters,
                    (float)(computations)/(cudot_time*(1.e+09)),
                    (float)(computations)/(cugemv_time*(1.e+09)),
                    (float)(computations)/(magma_minproductgemv_time*(1.e+09)),
                    (float)(computations)/(mdot_time*(1.e+09)),
                    (float)(computations)/(mdgm_time*(1.e+09)) );
            #endif

            magma_minproduct_smfree(&a, queue );
            magma_minproduct_smfree(&b, queue );
            magma_minproduct_smfree(&x, queue );
            magma_minproduct_smfree(&y, queue );
            magma_minproduct_smfree(&skp, queue );
        }

        printf("#================================================================================================================================================\n");
        printf("\n");
        printf("\n");
    }

cleanup:
    magma_minproduct_smfree(&a, queue );
    magma_minproduct_smfree(&b, queue );
    magma_minproduct_smfree(&x, queue );
    magma_minproduct_smfree(&y, queue );
    magma_minproduct_smfree(&skp, queue );
    magma_minproduct_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
