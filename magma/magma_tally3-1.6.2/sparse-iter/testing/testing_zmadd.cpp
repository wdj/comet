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
   -- testing csr matrix add
*/
int main(  int argc, char** argv )
{
    magma_tally3_int_t info = 0;
    TESTING_INIT();
    
    magma_tally3_queue_t queue=NULL;
    magma_tally3_queue_create( /*devices[ opts->device ],*/ &queue );

    real_Double_t res;
    magma_tally3_z_matrix A={Magma_tally3_CSR}, B={Magma_tally3_CSR}, B2={Magma_tally3_CSR}, 
    A_d={Magma_tally3_CSR}, B_d={Magma_tally3_CSR}, C_d={Magma_tally3_CSR};

    magma_tally3DoubleComplex one = MAGMA_tally3_Z_MAKE(1.0, 0.0);
    magma_tally3DoubleComplex mone = MAGMA_tally3_Z_MAKE(-1.0, 0.0);

    magma_tally3_int_t i=1;

    if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
        i++;
        magma_tally3_int_t laplace_size = atoi( argv[i] );
        CHECK( magma_tally3_zm_5stencil(  laplace_size, &A, queue ));
    } else {                        // file-matrix test
        CHECK( magma_tally3_z_csr_mtx( &A,  argv[i], queue ));
    }
    printf( "# matrix info: %d-by-%d with %d nonzeros\n",
                        (int) A.num_rows,(int) A.num_cols,(int) A.nnz );
    i++;

    if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
        i++;
        magma_tally3_int_t laplace_size = atoi( argv[i] );
        CHECK( magma_tally3_zm_5stencil(  laplace_size, &B, queue ));
    } else {                        // file-matrix test
        CHECK( magma_tally3_z_csr_mtx( &B,  argv[i], queue ));
    }
    printf( "# matrix info: %d-by-%d with %d nonzeros\n",
                        (int) B.num_rows,(int) B.num_cols,(int) B.nnz );


    CHECK( magma_tally3_zmtransfer( A, &A_d, Magma_tally3_CPU, Magma_tally3_DEV, queue ));
    CHECK( magma_tally3_zmtransfer( B, &B_d, Magma_tally3_CPU, Magma_tally3_DEV, queue ));

    CHECK( magma_tally3_zcuspaxpy( &one, A_d, &one, B_d, &C_d, queue ));

    magma_tally3_zmfree(&B_d, queue );

    CHECK( magma_tally3_zcuspaxpy( &mone, A_d, &one, C_d, &B_d, queue ));
    
    CHECK( magma_tally3_zmtransfer( B_d, &B2, Magma_tally3_DEV, Magma_tally3_CPU, queue ));

    magma_tally3_zmfree(&A_d, queue );
    magma_tally3_zmfree(&B_d, queue );
    magma_tally3_zmfree(&C_d, queue );

    // check difference
    CHECK( magma_tally3_zmdiff( B, B2, &res, queue ));
    printf("# ||A-B||_F = %8.2e\n", res);
    if ( res < .000001 )
        printf("# tester matrix add:  ok\n");
    else
        printf("# tester matrix add:  failed\n");

    magma_tally3_zmfree(&A, queue );
    magma_tally3_zmfree(&B, queue );
    magma_tally3_zmfree(&B2, queue );

cleanup:
    magma_tally3_zmfree(&A_d, queue );
    magma_tally3_zmfree(&B_d, queue );
    magma_tally3_zmfree(&C_d, queue );
    magma_tally3_zmfree(&A, queue );
    magma_tally3_zmfree(&B, queue );
    magma_tally3_zmfree(&B2, queue );
    magma_tally3_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
