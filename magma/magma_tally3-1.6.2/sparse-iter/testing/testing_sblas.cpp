/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from testing_zblas.cpp normal z -> s, Sun May  3 11:23:02 2015
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
   -- testing any solver
*/
int main(  int argc, char** argv )
{
    magma_tally3_int_t info = 0;
    /* Initialize */
    TESTING_INIT();
    magma_tally3_queue_t queue=NULL;
    magma_tally3_queue_create( &queue );
    magma_tally3blasSetKernelStream( queue );

    magma_tally3_int_t j, n=1000000, FLOPS;
    
    float one = MAGMA_tally3_S_MAKE( 1.0, 0.0 );
    float two = MAGMA_tally3_S_MAKE( 2.0, 0.0 );

    magma_tally3_s_matrix a={Magma_tally3_CSR}, ad={Magma_tally3_CSR}, bd={Magma_tally3_CSR}, cd={Magma_tally3_CSR};
    CHECK( magma_tally3_svinit( &a, Magma_tally3_CPU, n, 1, one, queue ));
    CHECK( magma_tally3_svinit( &bd, Magma_tally3_DEV, n, 1, two, queue ));
    CHECK( magma_tally3_svinit( &cd, Magma_tally3_DEV, n, 1, one, queue ));
    
    CHECK( magma_tally3_smtransfer( a, &ad, Magma_tally3_CPU, Magma_tally3_DEV, queue ));

    real_Double_t start, end, res;
    
    FLOPS = 2*n;
    start = magma_tally3_sync_wtime( queue );
    for (j=0; j<100; j++)
        res = magma_tally3_snrm2(n, ad.dval, 1);
    end = magma_tally3_sync_wtime( queue );
    printf( " > MAGMA_tally3 nrm2: %.2e seconds %.2e GFLOP/s\n",
                                    (end-start)/100, FLOPS*100/1e9/(end-start) );
    FLOPS = n;
    start = magma_tally3_sync_wtime( queue );
    for (j=0; j<100; j++)
        magma_tally3_sscal( n, two, ad.dval, 1 );
    end = magma_tally3_sync_wtime( queue );
    printf( " > MAGMA_tally3 scal: %.2e seconds %.2e GFLOP/s\n",
                                    (end-start)/100, FLOPS*100/1e9/(end-start) );
    FLOPS = 2*n;
    start = magma_tally3_sync_wtime( queue );
    for (j=0; j<100; j++)
        magma_tally3_saxpy( n, one, ad.dval, 1, bd.dval, 1 );
    end = magma_tally3_sync_wtime( queue );
    printf( " > MAGMA_tally3 axpy: %.2e seconds %.2e GFLOP/s\n",
                                    (end-start)/100, FLOPS*100/1e9/(end-start) );
    FLOPS = n;
    start = magma_tally3_sync_wtime( queue );
    for (j=0; j<100; j++)
        magma_tally3_scopy( n, bd.dval, 1, ad.dval, 1 );
    end = magma_tally3_sync_wtime( queue );
    printf( " > MAGMA_tally3 copy: %.2e seconds %.2e GFLOP/s\n",
                                    (end-start)/100, FLOPS*100/1e9/(end-start) );
    FLOPS = 2*n;
    start = magma_tally3_sync_wtime( queue );
    for (j=0; j<100; j++)
        res = MAGMA_tally3_S_REAL( magma_tally3_sdot(n, ad.dval, 1, bd.dval, 1) );
    end = magma_tally3_sync_wtime( queue );
    printf( " > MAGMA_tally3 dotc: %.2e seconds %.2e GFLOP/s\n",
                                    (end-start)/100, FLOPS*100/1e9/(end-start) );

    printf("# tester BLAS:  ok\n");


    magma_tally3_smfree( &a, queue);
    magma_tally3_smfree(&ad, queue);
    magma_tally3_smfree(&bd, queue);
    magma_tally3_smfree(&cd, queue);

    
cleanup:
    magma_tally3_smfree( &a, queue);
    magma_tally3_smfree(&ad, queue);
    magma_tally3_smfree(&bd, queue);
    magma_tally3_smfree(&cd, queue);
    magma_tally3blasSetKernelStream( NULL );
    magma_tally3_queue_destroy( queue );
    magma_tally3_finalize();
    return info;
}
