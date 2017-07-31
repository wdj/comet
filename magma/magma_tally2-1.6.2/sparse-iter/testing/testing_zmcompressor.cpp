/*
    -- MAGMA_tally2 (version 1.6.2) --
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
    TESTING_INIT();

    magma_tally2_zopts zopts;
    magma_tally2_queue_t queue=NULL;
    magma_tally2_queue_create( /*devices[ opts->device ],*/ &queue );

    real_Double_t res;
    magma_tally2_z_matrix A={Magma_tally2_CSR}, AT={Magma_tally2_CSR}, A2={Magma_tally2_CSR}, 
    B={Magma_tally2_CSR}, B_d={Magma_tally2_CSR};
    
    int i=1;
    real_Double_t start, end;
    CHECK( magma_tally2_zparse_opts( argc, argv, &zopts, &i, queue ));

    B.blocksize = zopts.blocksize;
    B.alignment = zopts.alignment;

    while(  i < argc ) {

        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_tally2_int_t laplace_size = atoi( argv[i] );
            CHECK( magma_tally2_zm_5stencil(  laplace_size, &A, queue ));
        } else {                        // file-matrix test
            CHECK( magma_tally2_z_csr_mtx( &A,  argv[i], queue ));
        }

        printf( "\n# matrix info: %d-by-%d with %d nonzeros\n\n",
                            (int) A.num_rows,(int) A.num_cols,(int) A.nnz );

        // scale matrix
        CHECK( magma_tally2_zmscale( &A, zopts.scaling, queue ));

        // remove nonzeros in matrix
        start = magma_tally2_sync_wtime( queue );
        for (int j=0; j<10; j++)
            CHECK( magma_tally2_zmcsrcompressor( &A, queue ));
        end = magma_tally2_sync_wtime( queue );
        printf( " > MAGMA_tally2 CPU: %.2e seconds.\n", (end-start)/10 );
        // transpose
        CHECK( magma_tally2_zmtranspose( A, &AT, queue ));

        // convert, copy back and forth to check everything works
        CHECK( magma_tally2_zmconvert( AT, &B, Magma_tally2_CSR, Magma_tally2_CSR, queue ));
        magma_tally2_zmfree(&AT, queue );
        CHECK( magma_tally2_zmtransfer( B, &B_d, Magma_tally2_CPU, Magma_tally2_DEV, queue ));
        magma_tally2_zmfree(&B, queue );

        start = magma_tally2_sync_wtime( queue );
        for (int j=0; j<10; j++)
            CHECK( magma_tally2_zmcsrcompressor_gpu( &B_d, queue ));
        end = magma_tally2_sync_wtime( queue );
        printf( " > MAGMA_tally2 GPU: %.2e seconds.\n", (end-start)/10 );


        CHECK( magma_tally2_zmtransfer( B_d, &B, Magma_tally2_DEV, Magma_tally2_CPU, queue ));
        magma_tally2_zmfree(&B_d, queue );
        CHECK( magma_tally2_zmconvert( B, &AT, Magma_tally2_CSR, Magma_tally2_CSR, queue ));
        magma_tally2_zmfree(&B, queue );

        // transpose back
        CHECK( magma_tally2_zmtranspose( AT, &A2, queue ));
        magma_tally2_zmfree(&AT, queue );
        CHECK( magma_tally2_zmdiff( A, A2, &res, queue ));
        printf("# ||A-B||_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# tester matrix compressor:  ok\n");
        else
            printf("# tester matrix compressor:  failed\n");

        magma_tally2_zmfree(&A, queue );
        magma_tally2_zmfree(&A2, queue );

        i++;
    }
    
cleanup:
    magma_tally2_zmfree(&AT, queue );
    magma_tally2_zmfree(&B, queue );
    magma_tally2_zmfree(&A, queue );
    magma_tally2_zmfree(&A2, queue );
    magma_tally2_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
