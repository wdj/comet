/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
       @author Mark Gates
*/

#ifndef MAGMA_minproduct_ZGEHRD_H
#define MAGMA_minproduct_ZGEHRD_H

#include "magma_minproduct_types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct zgehrd_data
{
    magma_minproduct_int_t ngpu;
    
    magma_minproduct_int_t ldda;
    magma_minproduct_int_t ldv;
    magma_minproduct_int_t ldvd;
    
    magma_minproductDoubleComplex *A    [ Magma_minproductMaxGPUs ];  // ldda*nlocal
    magma_minproductDoubleComplex *V    [ Magma_minproductMaxGPUs ];  // ldv *nb, whole panel
    magma_minproductDoubleComplex *Vd   [ Magma_minproductMaxGPUs ];  // ldvd*nb, block-cyclic
    magma_minproductDoubleComplex *Y    [ Magma_minproductMaxGPUs ];  // ldda*nb
    magma_minproductDoubleComplex *W    [ Magma_minproductMaxGPUs ];  // ldda*nb
    magma_minproductDoubleComplex *Ti   [ Magma_minproductMaxGPUs ];  // nb*nb
    
    magma_minproduct_queue_t streams[ Magma_minproductMaxGPUs ];
};

#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_minproduct_ZGEHRD_H
