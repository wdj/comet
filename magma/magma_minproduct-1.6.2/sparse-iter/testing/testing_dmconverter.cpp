/*
    -- MAGMA_minproduct (version 1.6.2) --
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

    magma_minproduct_dopts zopts;
    magma_minproduct_queue_t queue=NULL;
    magma_minproduct_queue_create( /*devices[ opts->device ],*/ &queue );

    real_Double_t res;
    magma_minproduct_d_matrix Z={Magma_minproduct_CSR}, Z2={Magma_minproduct_CSR}, A={Magma_minproduct_CSR}, A2={Magma_minproduct_CSR}, 
    AT={Magma_minproduct_CSR}, AT2={Magma_minproduct_CSR}, B={Magma_minproduct_CSR};
    printf("check1\n");
    int i=1;
    CHECK( magma_minproduct_dparse_opts( argc, argv, &zopts, &i, queue ));

    B.blocksize = zopts.blocksize;
    B.alignment = zopts.alignment;

    while(  i < argc ) {

        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_minproduct_int_t laplace_size = atoi( argv[i] );
            CHECK( magma_minproduct_dm_5stencil(  laplace_size, &Z, queue ));
        } else {                        // file-matrix test
            CHECK( magma_minproduct_d_csr_mtx( &Z,  argv[i], queue ));
        }

        printf( "# matrix info: %d-by-%d with %d nonzeros\n",
                            (int) Z.num_rows,(int) Z.num_cols,(int) Z.nnz );
        
        // convert to be non-symmetric
        CHECK( magma_minproduct_dmconvert( Z, &A, Magma_minproduct_CSR, Magma_minproduct_CSRL, queue ));
        CHECK( magma_minproduct_dmconvert( Z, &B, Magma_minproduct_CSR, Magma_minproduct_CSRU, queue ));

        // transpose
        CHECK( magma_minproduct_dmtranspose( A, &AT, queue ));

        // quite some conversions
                    
        //ELL
        CHECK( magma_minproduct_dmconvert( AT, &AT2, Magma_minproduct_CSR, Magma_minproduct_ELL, queue ));
        magma_minproduct_dmfree(&AT, queue );
        CHECK( magma_minproduct_dmconvert( AT2, &AT, Magma_minproduct_ELL, Magma_minproduct_CSR, queue ));
        magma_minproduct_dmfree(&AT2, queue );
        //ELLPACKT
        CHECK( magma_minproduct_dmconvert( AT, &AT2, Magma_minproduct_CSR, Magma_minproduct_ELLPACKT, queue ));
        magma_minproduct_dmfree(&AT, queue );
        CHECK( magma_minproduct_dmconvert( AT2, &AT, Magma_minproduct_ELLPACKT, Magma_minproduct_CSR, queue ));
        magma_minproduct_dmfree(&AT2, queue );
        //ELLRT
        AT2.blocksize = 8;
        AT2.alignment = 8;
        CHECK( magma_minproduct_dmconvert( AT, &AT2, Magma_minproduct_CSR, Magma_minproduct_ELLRT, queue ));
        magma_minproduct_dmfree(&AT, queue );
        CHECK( magma_minproduct_dmconvert( AT2, &AT, Magma_minproduct_ELLRT, Magma_minproduct_CSR, queue ));
        magma_minproduct_dmfree(&AT2, queue );
        //SELLP
        AT2.blocksize = 8;
        AT2.alignment = 8;
        CHECK( magma_minproduct_dmconvert( AT, &AT2, Magma_minproduct_CSR, Magma_minproduct_SELLP, queue ));
        magma_minproduct_dmfree(&AT, queue );
        CHECK( magma_minproduct_dmconvert( AT2, &AT, Magma_minproduct_SELLP, Magma_minproduct_CSR, queue ));
        magma_minproduct_dmfree(&AT2, queue );
        //ELLD
        CHECK( magma_minproduct_dmconvert( AT, &AT2, Magma_minproduct_CSR, Magma_minproduct_ELLD, queue ));
        magma_minproduct_dmfree(&AT, queue );
        CHECK( magma_minproduct_dmconvert( AT2, &AT, Magma_minproduct_ELLD, Magma_minproduct_CSR, queue ));
        magma_minproduct_dmfree(&AT2, queue );
        //CSRCOO
        CHECK( magma_minproduct_dmconvert( AT, &AT2, Magma_minproduct_CSR, Magma_minproduct_CSRCOO, queue ));
        magma_minproduct_dmfree(&AT, queue );
        CHECK( magma_minproduct_dmconvert( AT2, &AT, Magma_minproduct_CSRCOO, Magma_minproduct_CSR, queue ));
        magma_minproduct_dmfree(&AT2, queue );
        //CSRD
        CHECK( magma_minproduct_dmconvert( AT, &AT2, Magma_minproduct_CSR, Magma_minproduct_CSRD, queue ));
        magma_minproduct_dmfree(&AT, queue );
        CHECK( magma_minproduct_dmconvert( AT2, &AT, Magma_minproduct_CSRD, Magma_minproduct_CSR, queue ));
        magma_minproduct_dmfree(&AT2, queue );
        //BCSR
        CHECK( magma_minproduct_dmconvert( AT, &AT2, Magma_minproduct_CSR, Magma_minproduct_BCSR, queue ));
        magma_minproduct_dmfree(&AT, queue );
        CHECK( magma_minproduct_dmconvert( AT2, &AT, Magma_minproduct_BCSR, Magma_minproduct_CSR, queue ));
        magma_minproduct_dmfree(&AT2, queue );
        
        // transpose
        CHECK( magma_minproduct_dmtranspose( AT, &A2, queue ));
        
        CHECK( magma_minproduct_dmdiff( A, A2, &res, queue));
        printf("# ||A-A2||_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# conversion tester:  ok\n");
        else
            printf("# conversion tester:  failed\n");
        
        CHECK( magma_minproduct_dmlumerge( A2, B, &Z2, queue ));

        
        CHECK( magma_minproduct_dmdiff( Z, Z2, &res, queue));
        printf("# ||Z-Z2||_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# LUmerge tester:  ok\n");
        else
            printf("# LUmerge tester:  failed\n");


        magma_minproduct_dmfree(&A, queue );
        magma_minproduct_dmfree(&A2, queue );
        magma_minproduct_dmfree(&AT, queue );
        magma_minproduct_dmfree(&AT2, queue );
        magma_minproduct_dmfree(&B, queue );
        magma_minproduct_dmfree(&Z2, queue );
        magma_minproduct_dmfree(&Z, queue );

        i++;
    }

cleanup:
    magma_minproduct_dmfree(&A, queue );
    magma_minproduct_dmfree(&A2, queue );
    magma_minproduct_dmfree(&AT, queue );
    magma_minproduct_dmfree(&AT2, queue );
    magma_minproduct_dmfree(&B, queue );
    magma_minproduct_dmfree(&Z2, queue );
    magma_minproduct_dmfree(&Z, queue );
    
    magma_minproduct_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
