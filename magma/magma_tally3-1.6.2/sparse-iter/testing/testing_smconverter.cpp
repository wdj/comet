/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from testing_zmconverter.cpp normal z -> s, Sun May  3 11:23:02 2015
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
    magma_tally3_s_matrix Z={Magma_tally3_CSR}, Z2={Magma_tally3_CSR}, A={Magma_tally3_CSR}, A2={Magma_tally3_CSR}, 
    AT={Magma_tally3_CSR}, AT2={Magma_tally3_CSR}, B={Magma_tally3_CSR};
    printf("check1\n");
    int i=1;
    CHECK( magma_tally3_sparse_opts( argc, argv, &zopts, &i, queue ));

    B.blocksize = zopts.blocksize;
    B.alignment = zopts.alignment;

    while(  i < argc ) {

        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_tally3_int_t laplace_size = atoi( argv[i] );
            CHECK( magma_tally3_sm_5stencil(  laplace_size, &Z, queue ));
        } else {                        // file-matrix test
            CHECK( magma_tally3_s_csr_mtx( &Z,  argv[i], queue ));
        }

        printf( "# matrix info: %d-by-%d with %d nonzeros\n",
                            (int) Z.num_rows,(int) Z.num_cols,(int) Z.nnz );
        
        // convert to be non-symmetric
        CHECK( magma_tally3_smconvert( Z, &A, Magma_tally3_CSR, Magma_tally3_CSRL, queue ));
        CHECK( magma_tally3_smconvert( Z, &B, Magma_tally3_CSR, Magma_tally3_CSRU, queue ));

        // transpose
        CHECK( magma_tally3_smtranspose( A, &AT, queue ));

        // quite some conversions
                    
        //ELL
        CHECK( magma_tally3_smconvert( AT, &AT2, Magma_tally3_CSR, Magma_tally3_ELL, queue ));
        magma_tally3_smfree(&AT, queue );
        CHECK( magma_tally3_smconvert( AT2, &AT, Magma_tally3_ELL, Magma_tally3_CSR, queue ));
        magma_tally3_smfree(&AT2, queue );
        //ELLPACKT
        CHECK( magma_tally3_smconvert( AT, &AT2, Magma_tally3_CSR, Magma_tally3_ELLPACKT, queue ));
        magma_tally3_smfree(&AT, queue );
        CHECK( magma_tally3_smconvert( AT2, &AT, Magma_tally3_ELLPACKT, Magma_tally3_CSR, queue ));
        magma_tally3_smfree(&AT2, queue );
        //ELLRT
        AT2.blocksize = 8;
        AT2.alignment = 8;
        CHECK( magma_tally3_smconvert( AT, &AT2, Magma_tally3_CSR, Magma_tally3_ELLRT, queue ));
        magma_tally3_smfree(&AT, queue );
        CHECK( magma_tally3_smconvert( AT2, &AT, Magma_tally3_ELLRT, Magma_tally3_CSR, queue ));
        magma_tally3_smfree(&AT2, queue );
        //SELLP
        AT2.blocksize = 8;
        AT2.alignment = 8;
        CHECK( magma_tally3_smconvert( AT, &AT2, Magma_tally3_CSR, Magma_tally3_SELLP, queue ));
        magma_tally3_smfree(&AT, queue );
        CHECK( magma_tally3_smconvert( AT2, &AT, Magma_tally3_SELLP, Magma_tally3_CSR, queue ));
        magma_tally3_smfree(&AT2, queue );
        //ELLD
        CHECK( magma_tally3_smconvert( AT, &AT2, Magma_tally3_CSR, Magma_tally3_ELLD, queue ));
        magma_tally3_smfree(&AT, queue );
        CHECK( magma_tally3_smconvert( AT2, &AT, Magma_tally3_ELLD, Magma_tally3_CSR, queue ));
        magma_tally3_smfree(&AT2, queue );
        //CSRCOO
        CHECK( magma_tally3_smconvert( AT, &AT2, Magma_tally3_CSR, Magma_tally3_CSRCOO, queue ));
        magma_tally3_smfree(&AT, queue );
        CHECK( magma_tally3_smconvert( AT2, &AT, Magma_tally3_CSRCOO, Magma_tally3_CSR, queue ));
        magma_tally3_smfree(&AT2, queue );
        //CSRD
        CHECK( magma_tally3_smconvert( AT, &AT2, Magma_tally3_CSR, Magma_tally3_CSRD, queue ));
        magma_tally3_smfree(&AT, queue );
        CHECK( magma_tally3_smconvert( AT2, &AT, Magma_tally3_CSRD, Magma_tally3_CSR, queue ));
        magma_tally3_smfree(&AT2, queue );
        //BCSR
        CHECK( magma_tally3_smconvert( AT, &AT2, Magma_tally3_CSR, Magma_tally3_BCSR, queue ));
        magma_tally3_smfree(&AT, queue );
        CHECK( magma_tally3_smconvert( AT2, &AT, Magma_tally3_BCSR, Magma_tally3_CSR, queue ));
        magma_tally3_smfree(&AT2, queue );
        
        // transpose
        CHECK( magma_tally3_smtranspose( AT, &A2, queue ));
        
        CHECK( magma_tally3_smdiff( A, A2, &res, queue));
        printf("# ||A-A2||_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# conversion tester:  ok\n");
        else
            printf("# conversion tester:  failed\n");
        
        CHECK( magma_tally3_smlumerge( A2, B, &Z2, queue ));

        
        CHECK( magma_tally3_smdiff( Z, Z2, &res, queue));
        printf("# ||Z-Z2||_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# LUmerge tester:  ok\n");
        else
            printf("# LUmerge tester:  failed\n");


        magma_tally3_smfree(&A, queue );
        magma_tally3_smfree(&A2, queue );
        magma_tally3_smfree(&AT, queue );
        magma_tally3_smfree(&AT2, queue );
        magma_tally3_smfree(&B, queue );
        magma_tally3_smfree(&Z2, queue );
        magma_tally3_smfree(&Z, queue );

        i++;
    }

cleanup:
    magma_tally3_smfree(&A, queue );
    magma_tally3_smfree(&A2, queue );
    magma_tally3_smfree(&AT, queue );
    magma_tally3_smfree(&AT2, queue );
    magma_tally3_smfree(&B, queue );
    magma_tally3_smfree(&Z2, queue );
    magma_tally3_smfree(&Z, queue );
    
    magma_tally3_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
