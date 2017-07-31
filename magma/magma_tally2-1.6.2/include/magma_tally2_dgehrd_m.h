/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally2_zgehrd_m.h normal z -> d, Fri Jan 30 19:00:05 2015
       @author Mark Gates
*/

#ifndef MAGMA_tally2_DGEHRD_H
#define MAGMA_tally2_DGEHRD_H

#include "magma_tally2_types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct dgehrd_data_tally2
{
    magma_tally2_int_t ngpu;
    
    magma_tally2_int_t ldda;
    magma_tally2_int_t ldv;
    magma_tally2_int_t ldvd;
    
    double *A    [ Magma_tally2MaxGPUs ];  // ldda*nlocal
    double *V    [ Magma_tally2MaxGPUs ];  // ldv *nb, whole panel
    double *Vd   [ Magma_tally2MaxGPUs ];  // ldvd*nb, block-cyclic
    double *Y    [ Magma_tally2MaxGPUs ];  // ldda*nb
    double *W    [ Magma_tally2MaxGPUs ];  // ldda*nb
    double *Ti   [ Magma_tally2MaxGPUs ];  // nb*nb
    
    magma_tally2_queue_t streams[ Magma_tally2MaxGPUs ];
};

#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_tally2_DGEHRD_H
