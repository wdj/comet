/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally4_zgehrd_m.h normal z -> d, Fri Jan 30 19:00:05 2015
       @author Mark Gates
*/

#ifndef MAGMA_tally4_DGEHRD_H
#define MAGMA_tally4_DGEHRD_H

#include "magma_tally4_types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct dgehrd_data_tally4
{
    magma_tally4_int_t ngpu;
    
    magma_tally4_int_t ldda;
    magma_tally4_int_t ldv;
    magma_tally4_int_t ldvd;
    
    double *A    [ Magma_tally4MaxGPUs ];  // ldda*nlocal
    double *V    [ Magma_tally4MaxGPUs ];  // ldv *nb, whole panel
    double *Vd   [ Magma_tally4MaxGPUs ];  // ldvd*nb, block-cyclic
    double *Y    [ Magma_tally4MaxGPUs ];  // ldda*nb
    double *W    [ Magma_tally4MaxGPUs ];  // ldda*nb
    double *Ti   [ Magma_tally4MaxGPUs ];  // nb*nb
    
    magma_tally4_queue_t streams[ Magma_tally4MaxGPUs ];
};

#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_tally4_DGEHRD_H
