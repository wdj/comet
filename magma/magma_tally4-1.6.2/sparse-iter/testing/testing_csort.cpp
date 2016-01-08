/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from testing_zsort.cpp normal z -> c, Sun May  3 11:23:02 2015
       @author Hartwig Anzt
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>


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
    /* Initialize */
    TESTING_INIT();
    magma_tally4_queue_t queue=NULL;
    magma_tally4_queue_create( &queue );
    magma_tally4blasSetKernelStream( queue );

    magma_tally4_int_t i=1, n=100;
    magma_tally4_index_t *x=NULL;
    
    magma_tally4_c_matrix A={Magma_tally4_CSR};

    CHECK( magma_tally4_index_malloc_cpu( &x, n ));
    printf("unsorted:\n");
    srand(time(NULL));
    for(magma_tally4_int_t i = 0; i<n; i++ ){
        int r = rand()%100;
        x[i] = r;
        printf("%d  ", x[i]);
    }
    printf("\n\n");
    
    printf("sorting...");
    CHECK( magma_tally4_cindexsort(x, 0, n-1, queue ));
    printf("done.\n\n");
    
    printf("sorted:\n");
    for(magma_tally4_int_t i = 0; i<n; i++ ){
        printf("%d  ", x[i]);
    }
    printf("\n\n");

    magma_tally4_free_cpu( x );
    
    while(  i < argc ) {

        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_tally4_int_t laplace_size = atoi( argv[i] );
            CHECK( magma_tally4_cm_5stencil(  laplace_size, &A, queue ));
        } else {                        // file-matrix test
            CHECK( magma_tally4_c_csr_mtx( &A,  argv[i], queue ));
        }

        printf( "\n# matrix info: %d-by-%d with %d nonzeros\n\n",
                            (int) A.num_rows,(int) A.num_cols,(int) A.nnz );
    
        CHECK( magma_tally4_index_malloc_cpu( &x, A.num_rows*10 ));
        magma_tally4_int_t num_ind = 0;

        CHECK( magma_tally4_cdomainoverlap( A.num_rows, &num_ind, A.row, A.col, x, queue ));
                printf("domain overlap indices:\n");
        for(magma_tally4_int_t j = 0; j<num_ind; j++ ){
            printf("%d  ", x[j]);
        }
        printf("\n\n");
        magma_tally4_free_cpu( x );
        magma_tally4_cmfree(&A, queue);
        
        i++;
        
    }

cleanup:
    magma_tally4_cmfree(&A, queue );
    magma_tally4blasSetKernelStream( NULL );
    magma_tally4_queue_destroy( queue );
    magma_tally4_finalize();
    return info;
}
