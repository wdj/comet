/*
    -- MAGMA_tally2 (version 1.6.2) --
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
#include "magma_tally2.h"
#include "magma_tally2_lapack.h"
#include "testings.h"
#include "common_magma_tally2sparse.h"



/* ////////////////////////////////////////////////////////////////////////////
   -- testing any solver
*/
int main(  int argc, char** argv )
{
    magma_tally2_int_t info = 0;
    /* Initialize */
    TESTING_INIT();
    magma_tally2_queue_t queue=NULL;
    magma_tally2_queue_create( &queue );
    magma_tally2blasSetKernelStream( queue );

    magma_tally2_int_t j, n=1000000, FLOPS;
    
    double one = MAGMA_tally2_D_MAKE( 1.0, 0.0 );
    double two = MAGMA_tally2_D_MAKE( 2.0, 0.0 );

    magma_tally2_d_matrix a={Magma_tally2_CSR}, ad={Magma_tally2_CSR}, bd={Magma_tally2_CSR}, cd={Magma_tally2_CSR};
    CHECK( magma_tally2_dvinit( &a, Magma_tally2_CPU, n, 1, one, queue ));
    CHECK( magma_tally2_dvinit( &bd, Magma_tally2_DEV, n, 1, two, queue ));
    CHECK( magma_tally2_dvinit( &cd, Magma_tally2_DEV, n, 1, one, queue ));
    
    CHECK( magma_tally2_dmtransfer( a, &ad, Magma_tally2_CPU, Magma_tally2_DEV, queue ));

    real_Double_t start, end, res;
    
    FLOPS = 2*n;
    start = magma_tally2_sync_wtime( queue );
    for (j=0; j<100; j++)
        res = magma_tally2_dnrm2(n, ad.dval, 1);
    end = magma_tally2_sync_wtime( queue );
    printf( " > MAGMA_tally2 nrm2: %.2e seconds %.2e GFLOP/s\n",
                                    (end-start)/100, FLOPS*100/1e9/(end-start) );
    FLOPS = n;
    start = magma_tally2_sync_wtime( queue );
    for (j=0; j<100; j++)
        magma_tally2_dscal( n, two, ad.dval, 1 );
    end = magma_tally2_sync_wtime( queue );
    printf( " > MAGMA_tally2 scal: %.2e seconds %.2e GFLOP/s\n",
                                    (end-start)/100, FLOPS*100/1e9/(end-start) );
    FLOPS = 2*n;
    start = magma_tally2_sync_wtime( queue );
    for (j=0; j<100; j++)
        magma_tally2_daxpy( n, one, ad.dval, 1, bd.dval, 1 );
    end = magma_tally2_sync_wtime( queue );
    printf( " > MAGMA_tally2 axpy: %.2e seconds %.2e GFLOP/s\n",
                                    (end-start)/100, FLOPS*100/1e9/(end-start) );
    FLOPS = n;
    start = magma_tally2_sync_wtime( queue );
    for (j=0; j<100; j++)
        magma_tally2_dcopy( n, bd.dval, 1, ad.dval, 1 );
    end = magma_tally2_sync_wtime( queue );
    printf( " > MAGMA_tally2 copy: %.2e seconds %.2e GFLOP/s\n",
                                    (end-start)/100, FLOPS*100/1e9/(end-start) );
    FLOPS = 2*n;
    start = magma_tally2_sync_wtime( queue );
    for (j=0; j<100; j++)
        res = MAGMA_tally2_D_REAL( magma_tally2_ddot(n, ad.dval, 1, bd.dval, 1) );
    end = magma_tally2_sync_wtime( queue );
    printf( " > MAGMA_tally2 dotc: %.2e seconds %.2e GFLOP/s\n",
                                    (end-start)/100, FLOPS*100/1e9/(end-start) );

    printf("# tester BLAS:  ok\n");


    magma_tally2_dmfree( &a, queue);
    magma_tally2_dmfree(&ad, queue);
    magma_tally2_dmfree(&bd, queue);
    magma_tally2_dmfree(&cd, queue);

    
cleanup:
    magma_tally2_dmfree( &a, queue);
    magma_tally2_dmfree(&ad, queue);
    magma_tally2_dmfree(&bd, queue);
    magma_tally2_dmfree(&cd, queue);
    magma_tally2blasSetKernelStream( NULL );
    magma_tally2_queue_destroy( queue );
    magma_tally2_finalize();
    return info;
}
