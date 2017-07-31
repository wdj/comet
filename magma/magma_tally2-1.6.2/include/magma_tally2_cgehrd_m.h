/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally2_zgehrd_m.h normal z -> c, Fri Jan 30 19:00:05 2015
       @author Mark Gates
*/

#ifndef MAGMA_tally2_CGEHRD_H
#define MAGMA_tally2_CGEHRD_H

#include "magma_tally2_types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct cgehrd_data_tally2
{
    magma_tally2_int_t ngpu;
    
    magma_tally2_int_t ldda;
    magma_tally2_int_t ldv;
    magma_tally2_int_t ldvd;
    
    magma_tally2FloatComplex *A    [ Magma_tally2MaxGPUs ];  // ldda*nlocal
    magma_tally2FloatComplex *V    [ Magma_tally2MaxGPUs ];  // ldv *nb, whole panel
    magma_tally2FloatComplex *Vd   [ Magma_tally2MaxGPUs ];  // ldvd*nb, block-cyclic
    magma_tally2FloatComplex *Y    [ Magma_tally2MaxGPUs ];  // ldda*nb
    magma_tally2FloatComplex *W    [ Magma_tally2MaxGPUs ];  // ldda*nb
    magma_tally2FloatComplex *Ti   [ Magma_tally2MaxGPUs ];  // nb*nb
    
    magma_tally2_queue_t streams[ Magma_tally2MaxGPUs ];
};

#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_tally2_CGEHRD_H
