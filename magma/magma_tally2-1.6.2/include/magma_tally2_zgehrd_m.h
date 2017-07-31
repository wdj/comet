/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
       @author Mark Gates
*/

#ifndef MAGMA_tally2_ZGEHRD_H
#define MAGMA_tally2_ZGEHRD_H

#include "magma_tally2_types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct zgehrd_data_tally2
{
    magma_tally2_int_t ngpu;
    
    magma_tally2_int_t ldda;
    magma_tally2_int_t ldv;
    magma_tally2_int_t ldvd;
    
    magma_tally2DoubleComplex *A    [ Magma_tally2MaxGPUs ];  // ldda*nlocal
    magma_tally2DoubleComplex *V    [ Magma_tally2MaxGPUs ];  // ldv *nb, whole panel
    magma_tally2DoubleComplex *Vd   [ Magma_tally2MaxGPUs ];  // ldvd*nb, block-cyclic
    magma_tally2DoubleComplex *Y    [ Magma_tally2MaxGPUs ];  // ldda*nb
    magma_tally2DoubleComplex *W    [ Magma_tally2MaxGPUs ];  // ldda*nb
    magma_tally2DoubleComplex *Ti   [ Magma_tally2MaxGPUs ];  // nb*nb
    
    magma_tally2_queue_t streams[ Magma_tally2MaxGPUs ];
};

#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_tally2_ZGEHRD_H
