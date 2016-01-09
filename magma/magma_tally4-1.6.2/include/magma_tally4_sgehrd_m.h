/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally4_zgehrd_m.h normal z -> s, Fri Jan 30 19:00:05 2015
       @author Mark Gates
*/

#ifndef MAGMA_tally4_SGEHRD_H
#define MAGMA_tally4_SGEHRD_H

#include "magma_tally4_types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct sgehrd_data_tally4
{
    magma_tally4_int_t ngpu;
    
    magma_tally4_int_t ldda;
    magma_tally4_int_t ldv;
    magma_tally4_int_t ldvd;
    
    float *A    [ Magma_tally4MaxGPUs ];  // ldda*nlocal
    float *V    [ Magma_tally4MaxGPUs ];  // ldv *nb, whole panel
    float *Vd   [ Magma_tally4MaxGPUs ];  // ldvd*nb, block-cyclic
    float *Y    [ Magma_tally4MaxGPUs ];  // ldda*nb
    float *W    [ Magma_tally4MaxGPUs ];  // ldda*nb
    float *Ti   [ Magma_tally4MaxGPUs ];  // nb*nb
    
    magma_tally4_queue_t streams[ Magma_tally4MaxGPUs ];
};

#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_tally4_SGEHRD_H
