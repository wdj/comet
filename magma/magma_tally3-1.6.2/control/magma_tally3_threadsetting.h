/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Azzam Haidar
*/

#ifndef MAGMA_tally3_THREADSETTING_H
#define MAGMA_tally3_THREADSETTING_H

#ifdef __cplusplus
extern "C" {
#endif
/***************************************************************************//**
 *  Internal routines
 **/
void magma_tally3_set_lapack_numthreads(magma_tally3_int_t numthreads);
magma_tally3_int_t magma_tally3_get_lapack_numthreads();
magma_tally3_int_t magma_tally3_get_parallel_numthreads();
/***************************************************************************/
#ifdef __cplusplus
}
#endif

#endif  // MAGMA_tally3_THREADSETTING_H
