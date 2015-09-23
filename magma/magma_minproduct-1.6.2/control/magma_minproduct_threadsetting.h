/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Azzam Haidar
*/

#ifndef MAGMA_minproduct_THREADSETTING_H
#define MAGMA_minproduct_THREADSETTING_H

#ifdef __cplusplus
extern "C" {
#endif
/***************************************************************************//**
 *  Internal routines
 **/
void magma_minproduct_set_lapack_numthreads(magma_minproduct_int_t numthreads);
magma_minproduct_int_t magma_minproduct_get_lapack_numthreads();
magma_minproduct_int_t magma_minproduct_get_parallel_numthreads();
/***************************************************************************/
#ifdef __cplusplus
}
#endif

#endif  // MAGMA_minproduct_THREADSETTING_H
