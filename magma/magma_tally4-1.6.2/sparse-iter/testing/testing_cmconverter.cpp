/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from testing_zmconverter.cpp normal z -> c, Sun May  3 11:23:02 2015
       @author Hartwig Anzt
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"
#include "testings.h"
#include "common_magma_tally4sparse.h"



/* ////////////////////////////////////////////////////////////////////////////
   -- testing any solver
*/
int main(  int argc, char** argv )
{
    magma_tally4_int_t info = 0;
    TESTING_INIT();

    magma_tally4_copts zopts;
    magma_tally4_queue_t queue=NULL;
    magma_tally4_queue_create( /*devices[ opts->device ],*/ &queue );

    real_Double_t res;
    magma_tally4_c_matrix Z={Magma_tally4_CSR}, Z2={Magma_tally4_CSR}, A={Magma_tally4_CSR}, A2={Magma_tally4_CSR}, 
    AT={Magma_tally4_CSR}, AT2={Magma_tally4_CSR}, B={Magma_tally4_CSR};
    printf("check1\n");
    int i=1;
    CHECK( magma_tally4_cparse_opts( argc, argv, &zopts, &i, queue ));

    B.blocksize = zopts.blocksize;
    B.alignment = zopts.alignment;

    while(  i < argc ) {

        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_tally4_int_t laplace_size = atoi( argv[i] );
            CHECK( magma_tally4_cm_5stencil(  laplace_size, &Z, queue ));
        } else {                        // file-matrix test
            CHECK( magma_tally4_c_csr_mtx( &Z,  argv[i], queue ));
        }

        printf( "# matrix info: %d-by-%d with %d nonzeros\n",
                            (int) Z.num_rows,(int) Z.num_cols,(int) Z.nnz );
        
        // convert to be non-symmetric
        CHECK( magma_tally4_cmconvert( Z, &A, Magma_tally4_CSR, Magma_tally4_CSRL, queue ));
        CHECK( magma_tally4_cmconvert( Z, &B, Magma_tally4_CSR, Magma_tally4_CSRU, queue ));

        // transpose
        CHECK( magma_tally4_cmtranspose( A, &AT, queue ));

        // quite some conversions
                    
        //ELL
        CHECK( magma_tally4_cmconvert( AT, &AT2, Magma_tally4_CSR, Magma_tally4_ELL, queue ));
        magma_tally4_cmfree(&AT, queue );
        CHECK( magma_tally4_cmconvert( AT2, &AT, Magma_tally4_ELL, Magma_tally4_CSR, queue ));
        magma_tally4_cmfree(&AT2, queue );
        //ELLPACKT
        CHECK( magma_tally4_cmconvert( AT, &AT2, Magma_tally4_CSR, Magma_tally4_ELLPACKT, queue ));
        magma_tally4_cmfree(&AT, queue );
        CHECK( magma_tally4_cmconvert( AT2, &AT, Magma_tally4_ELLPACKT, Magma_tally4_CSR, queue ));
        magma_tally4_cmfree(&AT2, queue );
        //ELLRT
        AT2.blocksize = 8;
        AT2.alignment = 8;
        CHECK( magma_tally4_cmconvert( AT, &AT2, Magma_tally4_CSR, Magma_tally4_ELLRT, queue ));
        magma_tally4_cmfree(&AT, queue );
        CHECK( magma_tally4_cmconvert( AT2, &AT, Magma_tally4_ELLRT, Magma_tally4_CSR, queue ));
        magma_tally4_cmfree(&AT2, queue );
        //SELLP
        AT2.blocksize = 8;
        AT2.alignment = 8;
        CHECK( magma_tally4_cmconvert( AT, &AT2, Magma_tally4_CSR, Magma_tally4_SELLP, queue ));
        magma_tally4_cmfree(&AT, queue );
        CHECK( magma_tally4_cmconvert( AT2, &AT, Magma_tally4_SELLP, Magma_tally4_CSR, queue ));
        magma_tally4_cmfree(&AT2, queue );
        //ELLD
        CHECK( magma_tally4_cmconvert( AT, &AT2, Magma_tally4_CSR, Magma_tally4_ELLD, queue ));
        magma_tally4_cmfree(&AT, queue );
        CHECK( magma_tally4_cmconvert( AT2, &AT, Magma_tally4_ELLD, Magma_tally4_CSR, queue ));
        magma_tally4_cmfree(&AT2, queue );
        //CSRCOO
        CHECK( magma_tally4_cmconvert( AT, &AT2, Magma_tally4_CSR, Magma_tally4_CSRCOO, queue ));
        magma_tally4_cmfree(&AT, queue );
        CHECK( magma_tally4_cmconvert( AT2, &AT, Magma_tally4_CSRCOO, Magma_tally4_CSR, queue ));
        magma_tally4_cmfree(&AT2, queue );
        //CSRD
        CHECK( magma_tally4_cmconvert( AT, &AT2, Magma_tally4_CSR, Magma_tally4_CSRD, queue ));
        magma_tally4_cmfree(&AT, queue );
        CHECK( magma_tally4_cmconvert( AT2, &AT, Magma_tally4_CSRD, Magma_tally4_CSR, queue ));
        magma_tally4_cmfree(&AT2, queue );
        //BCSR
        CHECK( magma_tally4_cmconvert( AT, &AT2, Magma_tally4_CSR, Magma_tally4_BCSR, queue ));
        magma_tally4_cmfree(&AT, queue );
        CHECK( magma_tally4_cmconvert( AT2, &AT, Magma_tally4_BCSR, Magma_tally4_CSR, queue ));
        magma_tally4_cmfree(&AT2, queue );
        
        // transpose
        CHECK( magma_tally4_cmtranspose( AT, &A2, queue ));
        
        CHECK( magma_tally4_cmdiff( A, A2, &res, queue));
        printf("# ||A-A2||_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# conversion tester:  ok\n");
        else
            printf("# conversion tester:  failed\n");
        
        CHECK( magma_tally4_cmlumerge( A2, B, &Z2, queue ));

        
        CHECK( magma_tally4_cmdiff( Z, Z2, &res, queue));
        printf("# ||Z-Z2||_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# LUmerge tester:  ok\n");
        else
            printf("# LUmerge tester:  failed\n");


        magma_tally4_cmfree(&A, queue );
        magma_tally4_cmfree(&A2, queue );
        magma_tally4_cmfree(&AT, queue );
        magma_tally4_cmfree(&AT2, queue );
        magma_tally4_cmfree(&B, queue );
        magma_tally4_cmfree(&Z2, queue );
        magma_tally4_cmfree(&Z, queue );

        i++;
    }

cleanup:
    magma_tally4_cmfree(&A, queue );
    magma_tally4_cmfree(&A2, queue );
    magma_tally4_cmfree(&AT, queue );
    magma_tally4_cmfree(&AT2, queue );
    magma_tally4_cmfree(&B, queue );
    magma_tally4_cmfree(&Z2, queue );
    magma_tally4_cmfree(&Z, queue );
    
    magma_tally4_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
