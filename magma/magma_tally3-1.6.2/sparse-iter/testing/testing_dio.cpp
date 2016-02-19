/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from testing_zio.cpp normal z -> d, Sun May  3 11:23:02 2015
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

    magma_tally3_dopts zopts;
    magma_tally3_queue_t queue=NULL;
    magma_tally3_queue_create( /*devices[ opts->device ],*/ &queue );
    
    real_Double_t res;
    magma_tally3_d_matrix A={Magma_tally3_CSR}, A2={Magma_tally3_CSR}, 
    A3={Magma_tally3_CSR}, A4={Magma_tally3_CSR}, A5={Magma_tally3_CSR};
    
    int i=1;
    CHECK( magma_tally3_dparse_opts( argc, argv, &zopts, &i, queue ));

    while(  i < argc ) {

        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_tally3_int_t laplace_size = atoi( argv[i] );
            CHECK( magma_tally3_dm_5stencil(  laplace_size, &A, queue ));
        } else {                        // file-matrix test
            CHECK( magma_tally3_d_csr_mtx( &A,  argv[i], queue ));
        }

        printf( "# matrix info: %d-by-%d with %d nonzeros\n",
                            (int) A.num_rows,(int) A.num_cols,(int) A.nnz );

        // filename for temporary matrix storage
        const char *filename = "testmatrix.mtx";

        // write to file
        CHECK( magma_tally3_dwrite_csrtomtx( A, filename, queue ));
        // read from file
        CHECK( magma_tally3_d_csr_mtx( &A2, filename, queue ));

        // delete temporary matrix
        unlink( filename );
                
        //visualize
        printf("A2:\n");
        CHECK( magma_tally3_dprint_matrix( A2, queue ));
        
        //visualize
        CHECK( magma_tally3_dmconvert(A2, &A4, Magma_tally3_CSR, Magma_tally3_CSRL, queue ));
        printf("A4:\n");
        CHECK( magma_tally3_dprint_matrix( A4, queue ));
        CHECK( magma_tally3_dmconvert(A4, &A5, Magma_tally3_CSR, Magma_tally3_ELL, queue ));
        printf("A5:\n");
        CHECK( magma_tally3_dprint_matrix( A5, queue ));

        // pass it to another application and back
        magma_tally3_int_t m, n;
        magma_tally3_index_t *row, *col;
        double *val=NULL;
        CHECK( magma_tally3_dcsrget( A2, &m, &n, &row, &col, &val, queue ));
        CHECK( magma_tally3_dcsrset( m, n, row, col, val, &A3, queue ));

        CHECK( magma_tally3_dmdiff( A, A2, &res, queue ));
        printf("# ||A-B||_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# tester IO:  ok\n");
        else
            printf("# tester IO:  failed\n");

        CHECK( magma_tally3_dmdiff( A, A3, &res, queue ));
        printf("# ||A-B||_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# tester matrix interface:  ok\n");
        else
            printf("# tester matrix interface:  failed\n");

        magma_tally3_dmfree(&A, queue );
        magma_tally3_dmfree(&A2, queue );
        magma_tally3_dmfree(&A4, queue );
        magma_tally3_dmfree(&A5, queue );


        i++;
    }
    
cleanup:
    magma_tally3_dmfree(&A, queue );
    magma_tally3_dmfree(&A2, queue );
    magma_tally3_dmfree(&A4, queue );
    magma_tally3_dmfree(&A5, queue );
    magma_tally3_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
