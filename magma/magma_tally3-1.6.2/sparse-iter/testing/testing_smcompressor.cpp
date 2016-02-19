/*
    -- MAGMA_tally3 (version 1.6.2) --
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
    TESTING_INIT();

    magma_tally3_sopts zopts;
    magma_tally3_queue_t queue=NULL;
    magma_tally3_queue_create( /*devices[ opts->device ],*/ &queue );

    real_Double_t res;
    magma_tally3_s_matrix A={Magma_tally3_CSR}, AT={Magma_tally3_CSR}, A2={Magma_tally3_CSR}, 
    B={Magma_tally3_CSR}, B_d={Magma_tally3_CSR};
    
    int i=1;
    real_Double_t start, end;
    CHECK( magma_tally3_sparse_opts( argc, argv, &zopts, &i, queue ));

    B.blocksize = zopts.blocksize;
    B.alignment = zopts.alignment;

    while(  i < argc ) {

        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_tally3_int_t laplace_size = atoi( argv[i] );
            CHECK( magma_tally3_sm_5stencil(  laplace_size, &A, queue ));
        } else {                        // file-matrix test
            CHECK( magma_tally3_s_csr_mtx( &A,  argv[i], queue ));
        }

        printf( "\n# matrix info: %d-by-%d with %d nonzeros\n\n",
                            (int) A.num_rows,(int) A.num_cols,(int) A.nnz );

        // scale matrix
        CHECK( magma_tally3_smscale( &A, zopts.scaling, queue ));

        // remove nonzeros in matrix
        start = magma_tally3_sync_wtime( queue );
        for (int j=0; j<10; j++)
            CHECK( magma_tally3_smcsrcompressor( &A, queue ));
        end = magma_tally3_sync_wtime( queue );
        printf( " > MAGMA_tally3 CPU: %.2e seconds.\n", (end-start)/10 );
        // transpose
        CHECK( magma_tally3_smtranspose( A, &AT, queue ));

        // convert, copy back and forth to check everything works
        CHECK( magma_tally3_smconvert( AT, &B, Magma_tally3_CSR, Magma_tally3_CSR, queue ));
        magma_tally3_smfree(&AT, queue );
        CHECK( magma_tally3_smtransfer( B, &B_d, Magma_tally3_CPU, Magma_tally3_DEV, queue ));
        magma_tally3_smfree(&B, queue );

        start = magma_tally3_sync_wtime( queue );
        for (int j=0; j<10; j++)
            CHECK( magma_tally3_smcsrcompressor_gpu( &B_d, queue ));
        end = magma_tally3_sync_wtime( queue );
        printf( " > MAGMA_tally3 GPU: %.2e seconds.\n", (end-start)/10 );


        CHECK( magma_tally3_smtransfer( B_d, &B, Magma_tally3_DEV, Magma_tally3_CPU, queue ));
        magma_tally3_smfree(&B_d, queue );
        CHECK( magma_tally3_smconvert( B, &AT, Magma_tally3_CSR, Magma_tally3_CSR, queue ));
        magma_tally3_smfree(&B, queue );

        // transpose back
        CHECK( magma_tally3_smtranspose( AT, &A2, queue ));
        magma_tally3_smfree(&AT, queue );
        CHECK( magma_tally3_smdiff( A, A2, &res, queue ));
        printf("# ||A-B||_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# tester matrix compressor:  ok\n");
        else
            printf("# tester matrix compressor:  failed\n");

        magma_tally3_smfree(&A, queue );
        magma_tally3_smfree(&A2, queue );

        i++;
    }
    
cleanup:
    magma_tally3_smfree(&AT, queue );
    magma_tally3_smfree(&B, queue );
    magma_tally3_smfree(&A, queue );
    magma_tally3_smfree(&A2, queue );
    magma_tally3_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
