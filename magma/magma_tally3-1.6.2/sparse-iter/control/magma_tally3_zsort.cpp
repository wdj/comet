/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from
//  the IO functions provided by MatrixMarket

#include "common_magma_tally3sparse.h"


/**
    Purpose
    -------

    Sorts an array of integers.

    Arguments
    ---------

    @param[in,out]
    x           magma_tally3_index_t*
                array to sort

    @param[in]
    first       magma_tally3_int_t
                pointer to first element

    @param[in]
    last        magma_tally3_int_t
                pointer to last element

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zaux
    ********************************************************************/

extern "C"
magma_tally3_int_t
magma_tally3_zindexsort(
    magma_tally3_index_t *x,
    magma_tally3_int_t first,
    magma_tally3_int_t last,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    magma_tally3_index_t pivot,j,temp,i;

    if(first<last){
         pivot=first;
         i=first;
         j=last;

        while(i<j){
            while( x[i]<=x[pivot] && i<last )
                i++;
            while( x[j]>x[pivot] )
                j--;
            if( i<j ){
                temp = x[i];
                x[i] = x[j];
                x[j] = temp;
            }
        }

        temp=x[pivot];
        x[pivot]=x[j];
        x[j]=temp;
        CHECK( magma_tally3_zindexsort( x, first, j-1, queue ));
        CHECK( magma_tally3_zindexsort( x, j+1, last, queue ));

    }
cleanup:
    return info;

}

