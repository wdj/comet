/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally4_zsort.cpp normal z -> c, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from
//  the IO functions provided by MatrixMarket

#include "common_magma_tally4sparse.h"


/**
    Purpose
    -------

    Sorts an array of integers.

    Arguments
    ---------

    @param[in,out]
    x           magma_tally4_index_t*
                array to sort

    @param[in]
    first       magma_tally4_int_t
                pointer to first element

    @param[in]
    last        magma_tally4_int_t
                pointer to last element

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_caux
    ********************************************************************/

extern "C"
magma_tally4_int_t
magma_tally4_cindexsort(
    magma_tally4_index_t *x,
    magma_tally4_int_t first,
    magma_tally4_int_t last,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    magma_tally4_index_t pivot,j,temp,i;

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
        CHECK( magma_tally4_cindexsort( x, first, j-1, queue ));
        CHECK( magma_tally4_cindexsort( x, j+1, last, queue ));

    }
cleanup:
    return info;

}

