/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally3_zgehrd_m.h normal z -> s, Fri Jan 30 19:00:05 2015
       @author Mark Gates
*/

#ifndef MAGMA_tally3_SGEHRD_H
#define MAGMA_tally3_SGEHRD_H

#include "magma_tally3_types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct sgehrd_data_tally3
{
    magma_tally3_int_t ngpu;
    
    magma_tally3_int_t ldda;
    magma_tally3_int_t ldv;
    magma_tally3_int_t ldvd;
    
    float *A    [ Magma_tally3MaxGPUs ];  // ldda*nlocal
    float *V    [ Magma_tally3MaxGPUs ];  // ldv *nb, whole panel
    float *Vd   [ Magma_tally3MaxGPUs ];  // ldvd*nb, block-cyclic
    float *Y    [ Magma_tally3MaxGPUs ];  // ldda*nb
    float *W    [ Magma_tally3MaxGPUs ];  // ldda*nb
    float *Ti   [ Magma_tally3MaxGPUs ];  // nb*nb
    
    magma_tally3_queue_t streams[ Magma_tally3MaxGPUs ];
};

#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_tally3_SGEHRD_H
