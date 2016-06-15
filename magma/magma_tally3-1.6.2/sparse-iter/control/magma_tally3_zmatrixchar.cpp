/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt
*/
#include "common_magma_tally3sparse.h"

#define THRESHOLD 10e-99



/**
    Purpose
    -------

    Checks the maximal number of nonzeros in a row of matrix A.
    Inserts the data into max_nnz_row.


    Arguments
    ---------

    @param[in,out]
    A           magma_tally3_z_matrix*
                sparse matrix
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zaux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_zrowentries(
    magma_tally3_z_matrix *A,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    magma_tally3_index_t *length=NULL;
    magma_tally3_index_t i,j, maxrowlength=0;
    
    // check whether matrix on CPU
    if ( A->memory_location == Magma_tally3_CPU ) {
        // CSR
        if ( A->storage_type == Magma_tally3_CSR ) {
            CHECK( magma_tally3_index_malloc_cpu( &length, A->num_rows));
            for( i=0; i<A->num_rows; i++ ) {
                length[i] = A->row[i+1]-A->row[i];
                if (length[i] > maxrowlength)
                     maxrowlength = length[i];
            }
            A->max_nnz_row = maxrowlength;
        }
        // Dense
        else if ( A->storage_type == Magma_tally3_DENSE ) {
            CHECK( magma_tally3_index_malloc_cpu( &length, A->num_rows));

            for( i=0; i<A->num_rows; i++ ) {
                length[i] = 0;
                for( j=0; j<A->num_cols; j++ ) {
                    if ( MAGMA_tally3_Z_REAL( A->val[i*A->num_cols + j] ) != 0. )
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
        info = MAGMA_tally3_ERR_NOT_SUPPORTED;
    }
cleanup:
    magma_tally3_free( length );
    return info;
}


/**
    Purpose
    -------

    Computes the diameter of a sparse matrix and stores the value in diameter.


    Arguments
    ---------

    @param[in,out]
    A           magma_tally3_z_matrix*
                sparse matrix
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zaux
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_zdiameter(
    magma_tally3_z_matrix *A,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    magma_tally3_index_t i, j, tmp,  *dim=NULL, maxdim=0;
    
    // check whether matrix on CPU
    if ( A->memory_location == Magma_tally3_CPU ) {
        // CSR
        if ( A->storage_type == Magma_tally3_CSR ) {
            CHECK( magma_tally3_index_malloc_cpu( &dim, A->num_rows));
            for( i=0; i<A->num_rows; i++ ) {
                dim[i] = 0;
                for( j=A->row[i]; j<A->row[i+1]; j++ ) {
                   // if ( MAGMA_tally3_Z_REAL(A->val[j]) > THRESHOLD ) {
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
        else if ( A->storage_type == Magma_tally3_DENSE ) {
            magma_tally3_index_t i, j, tmp,  *dim, maxdim=0;
            CHECK( magma_tally3_index_malloc_cpu( &dim, A->num_rows));
            for( i=0; i<A->num_rows; i++ ) {
                dim[i] = 0;
                for( j=0; j<A->num_cols; j++ ) {
                    if ( MAGMA_tally3_Z_REAL( A->val[i*A->num_cols + j] ) !=  0.0 ) {
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
        else if ( A->storage_type == Magma_tally3_ELL ) {
            CHECK( magma_tally3_index_malloc_cpu( &dim, A->num_rows));
            for( i=0; i<A->num_rows; i++ ) {
                dim[i] = 0;
                for( j=i*A->max_nnz_row; j<(i+1)*A->max_nnz_row; j++ ) {
                    if ( MAGMA_tally3_Z_REAL( A->val[j] ) > THRESHOLD ) {
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
        else if ( A->storage_type == Magma_tally3_ELL ) {
            printf("error:format not supported.\n");
            info = MAGMA_tally3_ERR_NOT_SUPPORTED;
        }
    } // end CPU case

    else {
        printf("error: matrix not on CPU.\n");
        info = MAGMA_tally3_ERR_NOT_SUPPORTED;
    }
cleanup:
    magma_tally3_free( &dim );
    return info;
}