/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from testing_zmcompressor.cpp normal z -> s, Sun May  3 11:23:02 2015
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

    magma_minproduct_sopts zopts;
    magma_minproduct_queue_t queue=NULL;
    magma_minproduct_queue_create( /*devices[ opts->device ],*/ &queue );

    real_Double_t res;
    magma_minproduct_s_matrix A={Magma_minproduct_CSR}, AT={Magma_minproduct_CSR}, A2={Magma_minproduct_CSR}, 
    B={Magma_minproduct_CSR}, B_d={Magma_minproduct_CSR};
    
    int i=1;
    real_Double_t start, end;
    CHECK( magma_minproduct_sparse_opts( argc, argv, &zopts, &i, queue ));

    B.blocksize = zopts.blocksize;
    B.alignment = zopts.alignment;

    while(  i < argc ) {

        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_minproduct_int_t laplace_size = atoi( argv[i] );
            CHECK( magma_minproduct_sm_5stencil(  laplace_size, &A, queue ));
        } else {                        // file-matrix test
            CHECK( magma_minproduct_s_csr_mtx( &A,  argv[i], queue ));
        }

        printf( "\n# matrix info: %d-by-%d with %d nonzeros\n\n",
                            (int) A.num_rows,(int) A.num_cols,(int) A.nnz );

        // scale matrix
        CHECK( magma_minproduct_smscale( &A, zopts.scaling, queue ));

        // remove nonzeros in matrix
        start = magma_minproduct_sync_wtime( queue );
        for (int j=0; j<10; j++)
            CHECK( magma_minproduct_smcsrcompressor( &A, queue ));
        end = magma_minproduct_sync_wtime( queue );
        printf( " > MAGMA_minproduct CPU: %.2e seconds.\n", (end-start)/10 );
        // transpose
        CHECK( magma_minproduct_smtranspose( A, &AT, queue ));

        // convert, copy back and forth to check everything works
        CHECK( magma_minproduct_smconvert( AT, &B, Magma_minproduct_CSR, Magma_minproduct_CSR, queue ));
        magma_minproduct_smfree(&AT, queue );
        CHECK( magma_minproduct_smtransfer( B, &B_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
        magma_minproduct_smfree(&B, queue );

        start = magma_minproduct_sync_wtime( queue );
        for (int j=0; j<10; j++)
            CHECK( magma_minproduct_smcsrcompressor_gpu( &B_d, queue ));
        end = magma_minproduct_sync_wtime( queue );
        printf( " > MAGMA_minproduct GPU: %.2e seconds.\n", (end-start)/10 );


        CHECK( magma_minproduct_smtransfer( B_d, &B, Magma_minproduct_DEV, Magma_minproduct_CPU, queue ));
        magma_minproduct_smfree(&B_d, queue );
        CHECK( magma_minproduct_smconvert( B, &AT, Magma_minproduct_CSR, Magma_minproduct_CSR, queue ));
        magma_minproduct_smfree(&B, queue );

        // transpose back
        CHECK( magma_minproduct_smtranspose( AT, &A2, queue ));
        magma_minproduct_smfree(&AT, queue );
        CHECK( magma_minproduct_smdiff( A, A2, &res, queue ));
        printf("# ||A-B||_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# tester matrix compressor:  ok\n");
        else
            printf("# tester matrix compressor:  failed\n");

        magma_minproduct_smfree(&A, queue );
        magma_minproduct_smfree(&A2, queue );

        i++;
    }
    
cleanup:
    magma_minproduct_smfree(&AT, queue );
    magma_minproduct_smfree(&B, queue );
    magma_minproduct_smfree(&A, queue );
    magma_minproduct_smfree(&A2, queue );
    magma_minproduct_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
