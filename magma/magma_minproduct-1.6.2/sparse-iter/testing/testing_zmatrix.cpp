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
   -- testing any solver
*/
int main(  int argc, char** argv )
{
    magma_minproduct_int_t info = 0;
    TESTING_INIT();

    magma_minproduct_zopts zopts;
    magma_minproduct_queue_t queue=NULL;
    magma_minproduct_queue_create( /*devices[ opts->device ],*/ &queue );
    
    real_Double_t res;
    magma_minproduct_z_matrix Z={Magma_minproduct_CSR}, A={Magma_minproduct_CSR}, AT={Magma_minproduct_CSR}, 
    A2={Magma_minproduct_CSR}, B={Magma_minproduct_CSR}, B_d={Magma_minproduct_CSR};
    
    int i=1;
    CHECK( magma_minproduct_zparse_opts( argc, argv, &zopts, &i, queue ));

    B.blocksize = zopts.blocksize;
    B.alignment = zopts.alignment;

    while(  i < argc ) {

        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_minproduct_int_t laplace_size = atoi( argv[i] );
            CHECK( magma_minproduct_zm_5stencil(  laplace_size, &Z, queue ));
        } else {                        // file-matrix test
            CHECK( magma_minproduct_z_csr_mtx( &Z,  argv[i], queue ));
        }

        printf( "# matrix info: %d-by-%d with %d nonzeros\n",
                            (int) Z.num_rows,(int) Z.num_cols,(int) Z.nnz );

        // scale matrix
        CHECK( magma_minproduct_zmscale( &Z, zopts.scaling, queue ));

        // remove nonzeros in matrix
        CHECK( magma_minproduct_zmcsrcompressor( &Z, queue ));
        
        // convert to be non-symmetric
        CHECK( magma_minproduct_zmconvert( Z, &A, Magma_minproduct_CSR, Magma_minproduct_CSRL, queue ));
        
        // transpose
        CHECK( magma_minproduct_zmtranspose( A, &AT, queue ));

        // convert, copy back and forth to check everything works

        CHECK( magma_minproduct_zmconvert( AT, &B, Magma_minproduct_CSR, zopts.output_format, queue ));
        magma_minproduct_zmfree(&AT, queue );
        CHECK( magma_minproduct_zmtransfer( B, &B_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
        magma_minproduct_zmfree(&B, queue );
        CHECK( magma_minproduct_zmcsrcompressor_gpu( &B_d, queue ));
        CHECK( magma_minproduct_zmtransfer( B_d, &B, Magma_minproduct_DEV, Magma_minproduct_CPU, queue ));
        magma_minproduct_zmfree(&B_d, queue );
        CHECK( magma_minproduct_zmconvert( B, &AT, zopts.output_format,Magma_minproduct_CSR, queue ));
        magma_minproduct_zmfree(&B, queue );

        // transpose back
        CHECK( magma_minproduct_zmtranspose( AT, &A2, queue ));
        magma_minproduct_zmfree(&AT, queue );
        CHECK( magma_minproduct_zmdiff( A, A2, &res, queue));
        printf("# ||A-B||_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# tester:  ok\n");
        else
            printf("# tester:  failed\n");

        magma_minproduct_zmfree(&A, queue );
        magma_minproduct_zmfree(&A2, queue );
        magma_minproduct_zmfree(&Z, queue );

        i++;
    }

cleanup:
    magma_minproduct_zmfree(&AT, queue );
    magma_minproduct_zmfree(&A, queue );
    magma_minproduct_zmfree(&B, queue );
    magma_minproduct_zmfree(&B_d, queue );
    magma_minproduct_zmfree(&A2, queue );
    magma_minproduct_zmfree(&Z, queue );
    magma_minproduct_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
