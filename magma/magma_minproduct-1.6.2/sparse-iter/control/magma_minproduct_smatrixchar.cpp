/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_minproduct_zmatrixchar.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/
#include "common_magma_minproductsparse.h"

#define THRESHOLD 10e-99



/**
    Purpose
    -------

    Checks the maximal number of nonzeros in a row of matrix A.
    Inserts the data into max_nnz_row.


    Arguments
    ---------

    @param[in,out]
    A           magma_minproduct_s_matrix*
                sparse matrix
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_saux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_srowentries(
    magma_minproduct_s_matrix *A,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_index_t *length=NULL;
    magma_minproduct_index_t i,j, maxrowlength=0;
    
    // check whether matrix on CPU
    if ( A->memory_location == Magma_minproduct_CPU ) {
        // CSR
        if ( A->storage_type == Magma_minproduct_CSR ) {
            CHECK( magma_minproduct_index_malloc_cpu( &length, A->num_rows));
            for( i=0; i<A->num_rows; i++ ) {
                length[i] = A->row[i+1]-A->row[i];
                if (length[i] > maxrowlength)
                     maxrowlength = length[i];
            }
            A->max_nnz_row = maxrowlength;
        }
        // Dense
        else if ( A->storage_type == Magma_minproduct_DENSE ) {
            CHECK( magma_minproduct_index_malloc_cpu( &length, A->num_rows));

            for( i=0; i<A->num_rows; i++ ) {
                length[i] = 0;
                for( j=0; j<A->num_cols; j++ ) {
                    if ( MAGMA_minproduct_S_REAL( A->val[i*A->num_cols + j] ) != 0. )
                        length[i]++;
                    }
                if (length[i] > maxrowlength)
                     maxrowlength = length[i];
            }
            A->max_nnz_row = maxrowlength;
        }
    } // end CPU case

    else {
        printf("error: matrix not on CPU.\n");
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
    }
cleanup:
    magma_minproduct_free( length );
    return info;
}


/**
    Purpose
    -------

    Computes the diameter of a sparse matrix and stores the value in diameter.


    Arguments
    ---------

    @param[in,out]
    A           magma_minproduct_s_matrix*
                sparse matrix
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_saux
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_sdiameter(
    magma_minproduct_s_matrix *A,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_index_t i, j, tmp,  *dim=NULL, maxdim=0;
    
    // check whether matrix on CPU
    if ( A->memory_location == Magma_minproduct_CPU ) {
        // CSR
        if ( A->storage_type == Magma_minproduct_CSR ) {
            CHECK( magma_minproduct_index_malloc_cpu( &dim, A->num_rows));
            for( i=0; i<A->num_rows; i++ ) {
                dim[i] = 0;
                for( j=A->row[i]; j<A->row[i+1]; j++ ) {
                   // if ( MAGMA_minproduct_S_REAL(A->val[j]) > THRESHOLD ) {
                        tmp = abs( i - A->col[j] );
                        if ( tmp > dim[i] )
                            dim[i] = tmp;
                   // }
                }
                if ( dim[i] > maxdim )
                     maxdim = dim[i];
            }
            A->diameter = maxdim;
        }
        // Dense
        else if ( A->storage_type == Magma_minproduct_DENSE ) {
            magma_minproduct_index_t i, j, tmp,  *dim, maxdim=0;
            CHECK( magma_minproduct_index_malloc_cpu( &dim, A->num_rows));
            for( i=0; i<A->num_rows; i++ ) {
                dim[i] = 0;
                for( j=0; j<A->num_cols; j++ ) {
                    if ( MAGMA_minproduct_S_REAL( A->val[i*A->num_cols + j] ) !=  0.0 ) {
                        tmp = abs( i -j );
                        if ( tmp > dim[i] )
                            dim[i] = tmp;
                    }
                }
                if ( dim[i] > maxdim )
                     maxdim = dim[i];
            }
            A->diameter = maxdim;
        }
        // ELLPACK
        else if ( A->storage_type == Magma_minproduct_ELL ) {
            CHECK( magma_minproduct_index_malloc_cpu( &dim, A->num_rows));
            for( i=0; i<A->num_rows; i++ ) {
                dim[i] = 0;
                for( j=i*A->max_nnz_row; j<(i+1)*A->max_nnz_row; j++ ) {
                    if ( MAGMA_minproduct_S_REAL( A->val[j] ) > THRESHOLD ) {
                        tmp = abs( i - A->col[j] );
                        if ( tmp > dim[i] )
                            dim[i] = tmp;
                    }
                }
                if ( dim[i] > maxdim )
                     maxdim = dim[i];
            }
            A->diameter = maxdim;
        }
        // ELL
        else if ( A->storage_type == Magma_minproduct_ELL ) {
            printf("error:format not supported.\n");
            info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
        }
    } // end CPU case

    else {
        printf("error: matrix not on CPU.\n");
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
    }
cleanup:
    magma_minproduct_free( &dim );
    return info;
}
