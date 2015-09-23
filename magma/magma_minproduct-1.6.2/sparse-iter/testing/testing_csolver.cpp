/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from testing_zsolver.cpp normal z -> c, Sun May  3 11:23:02 2015
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

    magma_minproduct_copts zopts;
    magma_minproduct_queue_t queue=NULL;
    magma_minproduct_queue_create( /*devices[ opts->device ],*/ &queue );
    
    magma_minproductFloatComplex one = MAGMA_minproduct_C_MAKE(1.0, 0.0);
    magma_minproductFloatComplex zero = MAGMA_minproduct_C_MAKE(0.0, 0.0);
    magma_minproduct_c_matrix A={Magma_minproduct_CSR}, B={Magma_minproduct_CSR}, B_d={Magma_minproduct_CSR};
    magma_minproduct_c_matrix x={Magma_minproduct_CSR}, b={Magma_minproduct_CSR};
    
    int i=1;
    CHECK( magma_minproduct_cparse_opts( argc, argv, &zopts, &i, queue ));

    B.blocksize = zopts.blocksize;
    B.alignment = zopts.alignment;

    if ( zopts.solver_par.solver != Magma_minproduct_PCG &&
         zopts.solver_par.solver != Magma_minproduct_PGMRES &&
         zopts.solver_par.solver != Magma_minproduct_PBICGSTAB &&
         zopts.solver_par.solver != Magma_minproduct_ITERREF  &&
         zopts.solver_par.solver != Magma_minproduct_LOBPCG )
        zopts.precond_par.solver = Magma_minproduct_NONE;

    CHECK( magma_minproduct_csolverinfo_init( &zopts.solver_par, &zopts.precond_par, queue ));

    while(  i < argc ) {

        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_minproduct_int_t laplace_size = atoi( argv[i] );
            CHECK( magma_minproduct_cm_5stencil(  laplace_size, &A, queue ));
        } else {                        // file-matrix test
            CHECK( magma_minproduct_c_csr_mtx( &A,  argv[i], queue ));
        }

        printf( "\n# matrix info: %d-by-%d with %d nonzeros\n\n",
                            (int) A.num_rows,(int) A.num_cols,(int) A.nnz );


        // for the eigensolver case
        zopts.solver_par.ev_length = A.num_rows;
        CHECK( magma_minproduct_ceigensolverinfo_init( &zopts.solver_par, queue ));

        // scale matrix
        CHECK( magma_minproduct_cmscale( &A, zopts.scaling, queue ));

        CHECK( magma_minproduct_cmconvert( A, &B, Magma_minproduct_CSR, zopts.output_format, queue ));
        CHECK( magma_minproduct_cmtransfer( B, &B_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));

        // vectors and initial guess
        CHECK( magma_minproduct_cvinit( &b, Magma_minproduct_DEV, A.num_cols, 1, one, queue ));
        //magma_minproduct_cvinit( &x, Magma_minproduct_DEV, A.num_cols, 1, one, queue );
        //magma_minproduct_c_spmv( one, B_d, x, zero, b, queue );                 //  b = A x
        //magma_minproduct_cmfree(&x, queue );
        CHECK( magma_minproduct_cvinit( &x, Magma_minproduct_DEV, A.num_cols, 1, zero, queue ));
        
        info = magma_minproduct_c_solver( B_d, b, &x, &zopts, queue );
        if( info != 0 ){
            printf("error: solver returned: %s (%d).\n",
                magma_minproduct_strerror( info ), info );
        }
        magma_minproduct_csolverinfo( &zopts.solver_par, &zopts.precond_par, queue );

        magma_minproduct_cmfree(&B_d, queue );
        magma_minproduct_cmfree(&B, queue );
        magma_minproduct_cmfree(&A, queue );
        magma_minproduct_cmfree(&x, queue );
        magma_minproduct_cmfree(&b, queue );

        i++;
    }


    


cleanup:
    magma_minproduct_cmfree(&B_d, queue );
    magma_minproduct_cmfree(&B, queue );
    magma_minproduct_cmfree(&A, queue );
    magma_minproduct_cmfree(&x, queue );
    magma_minproduct_cmfree(&b, queue );
    magma_minproduct_csolverinfo_free( &zopts.solver_par, &zopts.precond_par, queue );
    magma_minproduct_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
