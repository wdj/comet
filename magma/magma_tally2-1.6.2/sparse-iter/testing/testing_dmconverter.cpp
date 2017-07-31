/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from testing_zmconverter.cpp normal z -> d, Sun May  3 11:23:02 2015
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

    magma_tally2_dopts zopts;
    magma_tally2_queue_t queue=NULL;
    magma_tally2_queue_create( /*devices[ opts->device ],*/ &queue );

    real_Double_t res;
    magma_tally2_d_matrix Z={Magma_tally2_CSR}, Z2={Magma_tally2_CSR}, A={Magma_tally2_CSR}, A2={Magma_tally2_CSR}, 
    AT={Magma_tally2_CSR}, AT2={Magma_tally2_CSR}, B={Magma_tally2_CSR};
    printf("check1\n");
    int i=1;
    CHECK( magma_tally2_dparse_opts( argc, argv, &zopts, &i, queue ));

    B.blocksize = zopts.blocksize;
    B.alignment = zopts.alignment;

    while(  i < argc ) {

        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_tally2_int_t laplace_size = atoi( argv[i] );
            CHECK( magma_tally2_dm_5stencil(  laplace_size, &Z, queue ));
        } else {                        // file-matrix test
            CHECK( magma_tally2_d_csr_mtx( &Z,  argv[i], queue ));
        }

        printf( "# matrix info: %d-by-%d with %d nonzeros\n",
                            (int) Z.num_rows,(int) Z.num_cols,(int) Z.nnz );
        
        // convert to be non-symmetric
        CHECK( magma_tally2_dmconvert( Z, &A, Magma_tally2_CSR, Magma_tally2_CSRL, queue ));
        CHECK( magma_tally2_dmconvert( Z, &B, Magma_tally2_CSR, Magma_tally2_CSRU, queue ));

        // transpose
        CHECK( magma_tally2_dmtranspose( A, &AT, queue ));

        // quite some conversions
                    
        //ELL
        CHECK( magma_tally2_dmconvert( AT, &AT2, Magma_tally2_CSR, Magma_tally2_ELL, queue ));
        magma_tally2_dmfree(&AT, queue );
        CHECK( magma_tally2_dmconvert( AT2, &AT, Magma_tally2_ELL, Magma_tally2_CSR, queue ));
        magma_tally2_dmfree(&AT2, queue );
        //ELLPACKT
        CHECK( magma_tally2_dmconvert( AT, &AT2, Magma_tally2_CSR, Magma_tally2_ELLPACKT, queue ));
        magma_tally2_dmfree(&AT, queue );
        CHECK( magma_tally2_dmconvert( AT2, &AT, Magma_tally2_ELLPACKT, Magma_tally2_CSR, queue ));
        magma_tally2_dmfree(&AT2, queue );
        //ELLRT
        AT2.blocksize = 8;
        AT2.alignment = 8;
        CHECK( magma_tally2_dmconvert( AT, &AT2, Magma_tally2_CSR, Magma_tally2_ELLRT, queue ));
        magma_tally2_dmfree(&AT, queue );
        CHECK( magma_tally2_dmconvert( AT2, &AT, Magma_tally2_ELLRT, Magma_tally2_CSR, queue ));
        magma_tally2_dmfree(&AT2, queue );
        //SELLP
        AT2.blocksize = 8;
        AT2.alignment = 8;
        CHECK( magma_tally2_dmconvert( AT, &AT2, Magma_tally2_CSR, Magma_tally2_SELLP, queue ));
        magma_tally2_dmfree(&AT, queue );
        CHECK( magma_tally2_dmconvert( AT2, &AT, Magma_tally2_SELLP, Magma_tally2_CSR, queue ));
        magma_tally2_dmfree(&AT2, queue );
        //ELLD
        CHECK( magma_tally2_dmconvert( AT, &AT2, Magma_tally2_CSR, Magma_tally2_ELLD, queue ));
        magma_tally2_dmfree(&AT, queue );
        CHECK( magma_tally2_dmconvert( AT2, &AT, Magma_tally2_ELLD, Magma_tally2_CSR, queue ));
        magma_tally2_dmfree(&AT2, queue );
        //CSRCOO
        CHECK( magma_tally2_dmconvert( AT, &AT2, Magma_tally2_CSR, Magma_tally2_CSRCOO, queue ));
        magma_tally2_dmfree(&AT, queue );
        CHECK( magma_tally2_dmconvert( AT2, &AT, Magma_tally2_CSRCOO, Magma_tally2_CSR, queue ));
        magma_tally2_dmfree(&AT2, queue );
        //CSRD
        CHECK( magma_tally2_dmconvert( AT, &AT2, Magma_tally2_CSR, Magma_tally2_CSRD, queue ));
        magma_tally2_dmfree(&AT, queue );
        CHECK( magma_tally2_dmconvert( AT2, &AT, Magma_tally2_CSRD, Magma_tally2_CSR, queue ));
        magma_tally2_dmfree(&AT2, queue );
        //BCSR
        CHECK( magma_tally2_dmconvert( AT, &AT2, Magma_tally2_CSR, Magma_tally2_BCSR, queue ));
        magma_tally2_dmfree(&AT, queue );
        CHECK( magma_tally2_dmconvert( AT2, &AT, Magma_tally2_BCSR, Magma_tally2_CSR, queue ));
        magma_tally2_dmfree(&AT2, queue );
        
        // transpose
        CHECK( magma_tally2_dmtranspose( AT, &A2, queue ));
        
        CHECK( magma_tally2_dmdiff( A, A2, &res, queue));
        printf("# ||A-A2||_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# conversion tester:  ok\n");
        else
            printf("# conversion tester:  failed\n");
        
        CHECK( magma_tally2_dmlumerge( A2, B, &Z2, queue ));

        
        CHECK( magma_tally2_dmdiff( Z, Z2, &res, queue));
        printf("# ||Z-Z2||_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# LUmerge tester:  ok\n");
        else
            printf("# LUmerge tester:  failed\n");


        magma_tally2_dmfree(&A, queue );
        magma_tally2_dmfree(&A2, queue );
        magma_tally2_dmfree(&AT, queue );
        magma_tally2_dmfree(&AT2, queue );
        magma_tally2_dmfree(&B, queue );
        magma_tally2_dmfree(&Z2, queue );
        magma_tally2_dmfree(&Z, queue );

        i++;
    }

cleanup:
    magma_tally2_dmfree(&A, queue );
    magma_tally2_dmfree(&A2, queue );
    magma_tally2_dmfree(&AT, queue );
    magma_tally2_dmfree(&AT2, queue );
    magma_tally2_dmfree(&B, queue );
    magma_tally2_dmfree(&Z2, queue );
    magma_tally2_dmfree(&Z, queue );
    
    magma_tally2_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
