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
   -- testing any solver
*/
int main(  int argc, char** argv )
{
    magma_tally3_int_t info = 0;
    TESTING_INIT();

    magma_tally3_zopts zopts;
    magma_tally3_queue_t queue=NULL;
    magma_tally3_queue_create( /*devices[ opts->device ],*/ &queue );
    
    real_Double_t res;
    magma_tally3_z_matrix Z={Magma_tally3_CSR}, A={Magma_tally3_CSR}, AT={Magma_tally3_CSR}, 
    A2={Magma_tally3_CSR}, B={Magma_tally3_CSR}, B_d={Magma_tally3_CSR};
    
    int i=1;
    CHECK( magma_tally3_zparse_opts( argc, argv, &zopts, &i, queue ));

    B.blocksize = zopts.blocksize;
    B.alignment = zopts.alignment;

    while(  i < argc ) {

        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_tally3_int_t laplace_size = atoi( argv[i] );
            CHECK( magma_tally3_zm_5stencil(  laplace_size, &Z, queue ));
        } else {                        // file-matrix test
            CHECK( magma_tally3_z_csr_mtx( &Z,  argv[i], queue ));
        }

        printf( "# matrix info: %d-by-%d with %d nonzeros\n",
                            (int) Z.num_rows,(int) Z.num_cols,(int) Z.nnz );

        // scale matrix
        CHECK( magma_tally3_zmscale( &Z, zopts.scaling, queue ));

        // remove nonzeros in matrix
        CHECK( magma_tally3_zmcsrcompressor( &Z, queue ));
        
        // convert to be non-symmetric
        CHECK( magma_tally3_zmconvert( Z, &A, Magma_tally3_CSR, Magma_tally3_CSRL, queue ));
        
        // transpose
        CHECK( magma_tally3_zmtranspose( A, &AT, queue ));

        // convert, copy back and forth to check everything works

        CHECK( magma_tally3_zmconvert( AT, &B, Magma_tally3_CSR, zopts.output_format, queue ));
        magma_tally3_zmfree(&AT, queue );
        CHECK( magma_tally3_zmtransfer( B, &B_d, Magma_tally3_CPU, Magma_tally3_DEV, queue ));
        magma_tally3_zmfree(&B, queue );
        CHECK( magma_tally3_zmcsrcompressor_gpu( &B_d, queue ));
        CHECK( magma_tally3_zmtransfer( B_d, &B, Magma_tally3_DEV, Magma_tally3_CPU, queue ));
        magma_tally3_zmfree(&B_d, queue );
        CHECK( magma_tally3_zmconvert( B, &AT, zopts.output_format,Magma_tally3_CSR, queue ));
        magma_tally3_zmfree(&B, queue );

        // transpose back
        CHECK( magma_tally3_zmtranspose( AT, &A2, queue ));
        magma_tally3_zmfree(&AT, queue );
        CHECK( magma_tally3_zmdiff( A, A2, &res, queue));
        printf("# ||A-B||_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# tester:  ok\n");
        else
            printf("# tester:  failed\n");

        magma_tally3_zmfree(&A, queue );
        magma_tally3_zmfree(&A2, queue );
        magma_tally3_zmfree(&Z, queue );

        i++;
    }

cleanup:
    magma_tally3_zmfree(&AT, queue );
    magma_tally3_zmfree(&A, queue );
    magma_tally3_zmfree(&B, queue );
    magma_tally3_zmfree(&B_d, queue );
    magma_tally3_zmfree(&A2, queue );
    magma_tally3_zmfree(&Z, queue );
    magma_tally3_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
