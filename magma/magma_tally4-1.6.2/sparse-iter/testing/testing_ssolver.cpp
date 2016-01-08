/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from testing_zsolver.cpp normal z -> s, Sun May  3 11:23:02 2015
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

    magma_tally4_sopts zopts;
    magma_tally4_queue_t queue=NULL;
    magma_tally4_queue_create( /*devices[ opts->device ],*/ &queue );
    
    float one = MAGMA_tally4_S_MAKE(1.0, 0.0);
    float zero = MAGMA_tally4_S_MAKE(0.0, 0.0);
    magma_tally4_s_matrix A={Magma_tally4_CSR}, B={Magma_tally4_CSR}, B_d={Magma_tally4_CSR};
    magma_tally4_s_matrix x={Magma_tally4_CSR}, b={Magma_tally4_CSR};
    
    int i=1;
    CHECK( magma_tally4_sparse_opts( argc, argv, &zopts, &i, queue ));

    B.blocksize = zopts.blocksize;
    B.alignment = zopts.alignment;

    if ( zopts.solver_par.solver != Magma_tally4_PCG &&
         zopts.solver_par.solver != Magma_tally4_PGMRES &&
         zopts.solver_par.solver != Magma_tally4_PBICGSTAB &&
         zopts.solver_par.solver != Magma_tally4_ITERREF  &&
         zopts.solver_par.solver != Magma_tally4_LOBPCG )
        zopts.precond_par.solver = Magma_tally4_NONE;

    CHECK( magma_tally4_ssolverinfo_init( &zopts.solver_par, &zopts.precond_par, queue ));

    while(  i < argc ) {

        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_tally4_int_t laplace_size = atoi( argv[i] );
            CHECK( magma_tally4_sm_5stencil(  laplace_size, &A, queue ));
        } else {                        // file-matrix test
            CHECK( magma_tally4_s_csr_mtx( &A,  argv[i], queue ));
        }

        printf( "\n# matrix info: %d-by-%d with %d nonzeros\n\n",
                            (int) A.num_rows,(int) A.num_cols,(int) A.nnz );


        // for the eigensolver case
        zopts.solver_par.ev_length = A.num_rows;
        CHECK( magma_tally4_seigensolverinfo_init( &zopts.solver_par, queue ));

        // scale matrix
        CHECK( magma_tally4_smscale( &A, zopts.scaling, queue ));

        CHECK( magma_tally4_smconvert( A, &B, Magma_tally4_CSR, zopts.output_format, queue ));
        CHECK( magma_tally4_smtransfer( B, &B_d, Magma_tally4_CPU, Magma_tally4_DEV, queue ));

        // vectors and initial guess
        CHECK( magma_tally4_svinit( &b, Magma_tally4_DEV, A.num_cols, 1, one, queue ));
        //magma_tally4_svinit( &x, Magma_tally4_DEV, A.num_cols, 1, one, queue );
        //magma_tally4_s_spmv( one, B_d, x, zero, b, queue );                 //  b = A x
        //magma_tally4_smfree(&x, queue );
        CHECK( magma_tally4_svinit( &x, Magma_tally4_DEV, A.num_cols, 1, zero, queue ));
        
        info = magma_tally4_s_solver( B_d, b, &x, &zopts, queue );
        if( info != 0 ){
            printf("error: solver returned: %s (%d).\n",
                magma_tally4_strerror( info ), info );
        }
        magma_tally4_ssolverinfo( &zopts.solver_par, &zopts.precond_par, queue );

        magma_tally4_smfree(&B_d, queue );
        magma_tally4_smfree(&B, queue );
        magma_tally4_smfree(&A, queue );
        magma_tally4_smfree(&x, queue );
        magma_tally4_smfree(&b, queue );

        i++;
    }


    


cleanup:
    magma_tally4_smfree(&B_d, queue );
    magma_tally4_smfree(&B, queue );
    magma_tally4_smfree(&A, queue );
    magma_tally4_smfree(&x, queue );
    magma_tally4_smfree(&b, queue );
    magma_tally4_ssolverinfo_free( &zopts.solver_par, &zopts.precond_par, queue );
    magma_tally4_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
