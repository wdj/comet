/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
       @author Mark Gates
*/

#ifndef MAGMA_tally4_ZGEHRD_H
#define MAGMA_tally4_ZGEHRD_H

#include "magma_tally4_types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct zgehrd_data
{
    magma_tally4_int_t ngpu;
    
    magma_tally4_int_t ldda;
    magma_tally4_int_t ldv;
    magma_tally4_int_t ldvd;
    
    magma_tally4DoubleComplex *A    [ Magma_tally4MaxGPUs ];  // ldda*nlocal
    magma_tally4DoubleComplex *V    [ Magma_tally4MaxGPUs ];  // ldv *nb, whole panel
    magma_tally4DoubleComplex *Vd   [ Magma_tally4MaxGPUs ];  // ldvd*nb, block-cyclic
    magma_tally4DoubleComplex *Y    [ Magma_tally4MaxGPUs ];  // ldda*nb
    magma_tally4DoubleComplex *W    [ Magma_tally4MaxGPUs ];  // ldda*nb
    magma_tally4DoubleComplex *Ti   [ Magma_tally4MaxGPUs ];  // nb*nb
    
    magma_tally4_queue_t streams[ Magma_tally4MaxGPUs ];
};

#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_tally4_ZGEHRD_H
