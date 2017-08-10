/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally2_zdomainoverlap.cpp normal z -> d, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/
#include "common_magma_tally2sparse.h"


extern "C"
magma_tally2_int_t
magma_tally2_dindexcopy(
    magma_tally2_int_t num_copy,
    magma_tally2_int_t offset,
    magma_tally2_index_t *tmp_x,
    magma_tally2_index_t *x,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    CHECK( magma_tally2_dindexsort( tmp_x, 0, num_copy-1, queue ));
    for( magma_tally2_int_t j=0; j<num_copy; j++ ){
        x[ j+offset ] = tmp_x[ j ];
        tmp_x[ j ] = -1;
    }
        
cleanup:    
    return info;
}


/**
    Purpose
    -------

    Generates the update list.

    Arguments
    ---------

    @param[in]
    x           magma_tally2_index_t*
                array to sort

    @param[in]
    num_rows    magma_tally2_int_t
                number of rows in matrix
                
    @param[out]
    num_indices magma_tally2_int_t*
                number of indices in array

    @param[in]
    rowptr      magma_tally2_index_t*
                rowpointer of matrix
                
    @param[in]
    colidx      magma_tally2_index_t*
                colindices of matrix
                
    @param[in]
    x           magma_tally2_index_t*
                array containing indices for domain overlap

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_daux
    ********************************************************************/

extern "C"
magma_tally2_int_t
magma_tally2_ddomainoverlap(
    magma_tally2_index_t num_rows,
    magma_tally2_index_t *num_indices,
    magma_tally2_index_t *rowptr,
    magma_tally2_index_t *colidx,
    magma_tally2_index_t *x,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    magma_tally2_int_t blocksize=128;

    magma_tally2_int_t row=0, col=0, num_ind=0, offset=0;

    magma_tally2_index_t *tmp_x;
    CHECK( magma_tally2_index_malloc_cpu( &tmp_x, blocksize ));

    for(magma_tally2_int_t i=0; i<blocksize; i++ ){
        tmp_x[i] = -1;
    }
    
    for(magma_tally2_int_t i=0; i<num_rows; i++){
        row = i;
        for(magma_tally2_int_t j=rowptr[row]; j<rowptr[row+1]; j++){
            col = colidx[j];
            int doubleitem = 0;
            for(magma_tally2_int_t k=0; k<blocksize; k++){
              if( tmp_x[k] == col )
                  doubleitem = 1;
            }
            if( doubleitem == 0 ){
                tmp_x[num_ind] = col;
                num_ind++;
                (*num_indices)++;
            }
            if( num_ind == blocksize || j == rowptr[num_rows]-1 ){
                magma_tally2_dindexcopy( num_ind, offset, tmp_x, x, queue );
                offset=offset+num_ind;
                num_ind = 0;
                break;
            }
        }
    }

cleanup:
    magma_tally2_free_cpu( tmp_x );
    return info;

}
