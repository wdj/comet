/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproduct_zgehrd_m.h normal z -> c, Fri Jan 30 19:00:05 2015
       @author Mark Gates
*/

#ifndef MAGMA_minproduct_CGEHRD_H
#define MAGMA_minproduct_CGEHRD_H

#include "magma_minproduct_types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct cgehrd_data
{
    magma_minproduct_int_t ngpu;
    
    magma_minproduct_int_t ldda;
    magma_minproduct_int_t ldv;
    magma_minproduct_int_t ldvd;
    
    magma_minproductFloatComplex *A    [ Magma_minproductMaxGPUs ];  // ldda*nlocal
    magma_minproductFloatComplex *V    [ Magma_minproductMaxGPUs ];  // ldv *nb, whole panel
    magma_minproductFloatComplex *Vd   [ Magma_minproductMaxGPUs ];  // ldvd*nb, block-cyclic
    magma_minproductFloatComplex *Y    [ Magma_minproductMaxGPUs ];  // ldda*nb
    magma_minproductFloatComplex *W    [ Magma_minproductMaxGPUs ];  // ldda*nb
    magma_minproductFloatComplex *Ti   [ Magma_minproductMaxGPUs ];  // nb*nb
    
    magma_minproduct_queue_t streams[ Magma_minproductMaxGPUs ];
};

#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_minproduct_CGEHRD_H
