/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from testing_zblas.cpp normal z -> d, Sun May  3 11:23:02 2015
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
   -- testing any solver
*/
int main(  int argc, char** argv )
{
    magma_minproduct_int_t info = 0;
    /* Initialize */
    TESTING_INIT();
    magma_minproduct_queue_t queue=NULL;
    magma_minproduct_queue_create( &queue );
    magma_minproductblasSetKernelStream( queue );

    magma_minproduct_int_t j, n=1000000, FLOPS;
    
    double one = MAGMA_minproduct_D_MAKE( 1.0, 0.0 );
    double two = MAGMA_minproduct_D_MAKE( 2.0, 0.0 );

    magma_minproduct_d_matrix a={Magma_minproduct_CSR}, ad={Magma_minproduct_CSR}, bd={Magma_minproduct_CSR}, cd={Magma_minproduct_CSR};
    CHECK( magma_minproduct_dvinit( &a, Magma_minproduct_CPU, n, 1, one, queue ));
    CHECK( magma_minproduct_dvinit( &bd, Magma_minproduct_DEV, n, 1, two, queue ));
    CHECK( magma_minproduct_dvinit( &cd, Magma_minproduct_DEV, n, 1, one, queue ));
    
    CHECK( magma_minproduct_dmtransfer( a, &ad, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));

    real_Double_t start, end, res;
    
    FLOPS = 2*n;
    start = magma_minproduct_sync_wtime( queue );
    for (j=0; j<100; j++)
        res = magma_minproduct_dnrm2(n, ad.dval, 1);
    end = magma_minproduct_sync_wtime( queue );
    printf( " > MAGMA_minproduct nrm2: %.2e seconds %.2e GFLOP/s\n",
                                    (end-start)/100, FLOPS*100/1e9/(end-start) );
    FLOPS = n;
    start = magma_minproduct_sync_wtime( queue );
    for (j=0; j<100; j++)
        magma_minproduct_dscal( n, two, ad.dval, 1 );
    end = magma_minproduct_sync_wtime( queue );
    printf( " > MAGMA_minproduct scal: %.2e seconds %.2e GFLOP/s\n",
                                    (end-start)/100, FLOPS*100/1e9/(end-start) );
    FLOPS = 2*n;
    start = magma_minproduct_sync_wtime( queue );
    for (j=0; j<100; j++)
        magma_minproduct_daxpy( n, one, ad.dval, 1, bd.dval, 1 );
    end = magma_minproduct_sync_wtime( queue );
    printf( " > MAGMA_minproduct axpy: %.2e seconds %.2e GFLOP/s\n",
                                    (end-start)/100, FLOPS*100/1e9/(end-start) );
    FLOPS = n;
    start = magma_minproduct_sync_wtime( queue );
    for (j=0; j<100; j++)
        magma_minproduct_dcopy( n, bd.dval, 1, ad.dval, 1 );
    end = magma_minproduct_sync_wtime( queue );
    printf( " > MAGMA_minproduct copy: %.2e seconds %.2e GFLOP/s\n",
                                    (end-start)/100, FLOPS*100/1e9/(end-start) );
    FLOPS = 2*n;
    start = magma_minproduct_sync_wtime( queue );
    for (j=0; j<100; j++)
        res = MAGMA_minproduct_D_REAL( magma_minproduct_ddot(n, ad.dval, 1, bd.dval, 1) );
    end = magma_minproduct_sync_wtime( queue );
    printf( " > MAGMA_minproduct dotc: %.2e seconds %.2e GFLOP/s\n",
                                    (end-start)/100, FLOPS*100/1e9/(end-start) );

    printf("# tester BLAS:  ok\n");


    magma_minproduct_dmfree( &a, queue);
    magma_minproduct_dmfree(&ad, queue);
    magma_minproduct_dmfree(&bd, queue);
    magma_minproduct_dmfree(&cd, queue);

    
cleanup:
    magma_minproduct_dmfree( &a, queue);
    magma_minproduct_dmfree(&ad, queue);
    magma_minproduct_dmfree(&bd, queue);
    magma_minproduct_dmfree(&cd, queue);
    magma_minproductblasSetKernelStream( NULL );
    magma_minproduct_queue_destroy( queue );
    magma_minproduct_finalize();
    return info;
}
