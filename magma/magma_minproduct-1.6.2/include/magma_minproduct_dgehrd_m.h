/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproduct_zgehrd_m.h normal z -> d, Fri Jan 30 19:00:05 2015
       @author Mark Gates
*/

#ifndef MAGMA_minproduct_DGEHRD_H
#define MAGMA_minproduct_DGEHRD_H

#include "magma_minproduct_types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct dgehrd_data
{
    magma_minproduct_int_t ngpu;
    
    magma_minproduct_int_t ldda;
    magma_minproduct_int_t ldv;
    magma_minproduct_int_t ldvd;
    
    double *A    [ Magma_minproductMaxGPUs ];  // ldda*nlocal
    double *V    [ Magma_minproductMaxGPUs ];  // ldv *nb, whole panel
    double *Vd   [ Magma_minproductMaxGPUs ];  // ldvd*nb, block-cyclic
    double *Y    [ Magma_minproductMaxGPUs ];  // ldda*nb
    double *W    [ Magma_minproductMaxGPUs ];  // ldda*nb
    double *Ti   [ Magma_minproductMaxGPUs ];  // nb*nb
    
    magma_minproduct_queue_t streams[ Magma_minproductMaxGPUs ];
};

#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_minproduct_DGEHRD_H
