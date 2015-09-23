/*
    -- MAGMA_minproduct (version 1.6.2) --
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
#include <time.h>


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
    /* Initialize */
    TESTING_INIT();
    magma_minproduct_queue_t queue=NULL;
    magma_minproduct_queue_create( &queue );
    magma_minproductblasSetKernelStream( queue );

    magma_minproduct_int_t i=1, n=100;
    magma_minproduct_index_t *x=NULL;
    
    magma_minproduct_z_matrix A={Magma_minproduct_CSR};

    CHECK( magma_minproduct_index_malloc_cpu( &x, n ));
    printf("unsorted:\n");
    srand(time(NULL));
    for(magma_minproduct_int_t i = 0; i<n; i++ ){
        int r = rand()%100;
        x[i] = r;
        printf("%d  ", x[i]);
    }
    printf("\n\n");
    
    printf("sorting...");
    CHECK( magma_minproduct_zindexsort(x, 0, n-1, queue ));
    printf("done.\n\n");
    
    printf("sorted:\n");
    for(magma_minproduct_int_t i = 0; i<n; i++ ){
        printf("%d  ", x[i]);
    }
    printf("\n\n");

    magma_minproduct_free_cpu( x );
    
    while(  i < argc ) {

        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_minproduct_int_t laplace_size = atoi( argv[i] );
            CHECK( magma_minproduct_zm_5stencil(  laplace_size, &A, queue ));
        } else {                        // file-matrix test
            CHECK( magma_minproduct_z_csr_mtx( &A,  argv[i], queue ));
        }

        printf( "\n# matrix info: %d-by-%d with %d nonzeros\n\n",
                            (int) A.num_rows,(int) A.num_cols,(int) A.nnz );
    
        CHECK( magma_minproduct_index_malloc_cpu( &x, A.num_rows*10 ));
        magma_minproduct_int_t num_ind = 0;

        CHECK( magma_minproduct_zdomainoverlap( A.num_rows, &num_ind, A.row, A.col, x, queue ));
                printf("domain overlap indices:\n");
        for(magma_minproduct_int_t j = 0; j<num_ind; j++ ){
            printf("%d  ", x[j]);
        }
        printf("\n\n");
        magma_minproduct_free_cpu( x );
        magma_minproduct_zmfree(&A, queue);
        
        i++;
        
    }

cleanup:
    magma_minproduct_zmfree(&A, queue );
    magma_minproductblasSetKernelStream( NULL );
    magma_minproduct_queue_destroy( queue );
    magma_minproduct_finalize();
    return info;
}
