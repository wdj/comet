/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
       @author Mark Gates
*/

#ifndef MAGMA_tally3_ZGEHRD_H
#define MAGMA_tally3_ZGEHRD_H

#include "magma_tally3_types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct zgehrd_data_tally3
{
    magma_tally3_int_t ngpu;
    
    magma_tally3_int_t ldda;
    magma_tally3_int_t ldv;
    magma_tally3_int_t ldvd;
    
    magma_tally3DoubleComplex *A    [ Magma_tally3MaxGPUs ];  // ldda*nlocal
    magma_tally3DoubleComplex *V    [ Magma_tally3MaxGPUs ];  // ldv *nb, whole panel
    magma_tally3DoubleComplex *Vd   [ Magma_tally3MaxGPUs ];  // ldvd*nb, block-cyclic
    magma_tally3DoubleComplex *Y    [ Magma_tally3MaxGPUs ];  // ldda*nb
    magma_tally3DoubleComplex *W    [ Magma_tally3MaxGPUs ];  // ldda*nb
    magma_tally3DoubleComplex *Ti   [ Magma_tally3MaxGPUs ];  // nb*nb
    
    magma_tally3_queue_t streams[ Magma_tally3MaxGPUs ];
};

#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_tally3_ZGEHRD_H
