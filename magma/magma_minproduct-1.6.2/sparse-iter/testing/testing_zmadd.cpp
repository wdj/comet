/*
    -- MAGMA_minproduct (version 1.6.2) --
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
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "testings.h"
#include "common_magma_minproductsparse.h"



/* ////////////////////////////////////////////////////////////////////////////
   -- testing csr matrix add
*/
int main(  int argc, char** argv )
{
    magma_minproduct_int_t info = 0;
    TESTING_INIT();
    
    magma_minproduct_queue_t queue=NULL;
    magma_minproduct_queue_create( /*devices[ opts->device ],*/ &queue );

    real_Double_t res;
    magma_minproduct_z_matrix A={Magma_minproduct_CSR}, B={Magma_minproduct_CSR}, B2={Magma_minproduct_CSR}, 
    A_d={Magma_minproduct_CSR}, B_d={Magma_minproduct_CSR}, C_d={Magma_minproduct_CSR};

    magma_minproductDoubleComplex one = MAGMA_minproduct_Z_MAKE(1.0, 0.0);
    magma_minproductDoubleComplex mone = MAGMA_minproduct_Z_MAKE(-1.0, 0.0);

    magma_minproduct_int_t i=1;

    if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
        i++;
        magma_minproduct_int_t laplace_size = atoi( argv[i] );
        CHECK( magma_minproduct_zm_5stencil(  laplace_size, &A, queue ));
    } else {                        // file-matrix test
        CHECK( magma_minproduct_z_csr_mtx( &A,  argv[i], queue ));
    }
    printf( "# matrix info: %d-by-%d with %d nonzeros\n",
                        (int) A.num_rows,(int) A.num_cols,(int) A.nnz );
    i++;

    if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
        i++;
        magma_minproduct_int_t laplace_size = atoi( argv[i] );
        CHECK( magma_minproduct_zm_5stencil(  laplace_size, &B, queue ));
    } else {                        // file-matrix test
        CHECK( magma_minproduct_z_csr_mtx( &B,  argv[i], queue ));
    }
    printf( "# matrix info: %d-by-%d with %d nonzeros\n",
                        (int) B.num_rows,(int) B.num_cols,(int) B.nnz );


    CHECK( magma_minproduct_zmtransfer( A, &A_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
    CHECK( magma_minproduct_zmtransfer( B, &B_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));

    CHECK( magma_minproduct_zcuspaxpy( &one, A_d, &one, B_d, &C_d, queue ));

    magma_minproduct_zmfree(&B_d, queue );

    CHECK( magma_minproduct_zcuspaxpy( &mone, A_d, &one, C_d, &B_d, queue ));
    
    CHECK( magma_minproduct_zmtransfer( B_d, &B2, Magma_minproduct_DEV, Magma_minproduct_CPU, queue ));

    magma_minproduct_zmfree(&A_d, queue );
    magma_minproduct_zmfree(&B_d, queue );
    magma_minproduct_zmfree(&C_d, queue );

    // check difference
    CHECK( magma_minproduct_zmdiff( B, B2, &res, queue ));
    printf("# ||A-B||_F = %8.2e\n", res);
    if ( res < .000001 )
        printf("# tester matrix add:  ok\n");
    else
        printf("# tester matrix add:  failed\n");

    magma_minproduct_zmfree(&A, queue );
    magma_minproduct_zmfree(&B, queue );
    magma_minproduct_zmfree(&B2, queue );

cleanup:
    magma_minproduct_zmfree(&A_d, queue );
    magma_minproduct_zmfree(&B_d, queue );
    magma_minproduct_zmfree(&C_d, queue );
    magma_minproduct_zmfree(&A, queue );
    magma_minproduct_zmfree(&B, queue );
    magma_minproduct_zmfree(&B2, queue );
    magma_minproduct_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
