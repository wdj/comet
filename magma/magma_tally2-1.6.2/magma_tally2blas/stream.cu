/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Mark Gates
*/

#include "common_magma_tally2.h"

magma_tally2_queue_t magma_tally2_stream = 0;


/**
    Purpose
    -------
    magma_tally2blasSetKernelStream sets the CUDA stream that MAGMA_tally2 BLAS and
    CUBLAS (v1) routines use (unless explicitly given a stream).
    
    In a multi-threaded application, be careful to avoid race conditions
    when using this. For instance, if calls are executed in this order:
    
    @verbatim
        thread 1                            thread 2
        ------------------------------      ------------------------------
    1.  magma_tally2blasSetKernelStream( s1 )         
    2.                                      magma_tally2blasSetKernelStream( s2 )
    3.  magma_tally2_dgemm( ... )
    4.                                      magma_tally2_dgemm( ... )
    @endverbatim
    
    both magma_tally2_dgemm would occur on stream s2. A lock should be used to prevent
    this, so the dgemm in thread 1 uses stream s1, and the dgemm in thread 2
    uses s2:
    
    @verbatim
        thread 1                            thread 2
        ------------------------------      ------------------------------
    1.  lock()                                  
    2.  magma_tally2blasSetKernelStream( s1 )          
    3.  magma_tally2_dgemm( ... )                      
    4.  unlock()                                
    5.                                      lock()
    6.                                      magma_tally2blasSetKernelStream( s2 )
    7.                                      magma_tally2_dgemm( ... )
    8.                                      unlock()
    @endverbatim
    
    Most BLAS calls in MAGMA_tally2, such as magma_tally2_dgemm, are asynchronous, so the lock
    will only have to wait until dgemm is queued, not until it is finished.
    
    Arguments
    ---------
    @param[in]
    stream  magma_tally2_queue_t
            The CUDA stream.

    @ingroup magma_tally2_util
    ********************************************************************/
extern "C"
cublasStatus_t magma_tally2blasSetKernelStream( magma_tally2_queue_t stream )
{
    magma_tally2_stream = stream;
    return cublasSetKernelStream( stream );
}


/**
    Purpose
    -------
    magma_tally2blasGetKernelStream gets the CUDA stream that MAGMA_tally2 BLAS
    routines use.

    Arguments
    ---------
    @param[out]
    stream  magma_tally2_queue_t
            The CUDA stream.

    @ingroup magma_tally2_util
    ********************************************************************/
extern "C"
cublasStatus_t magma_tally2blasGetKernelStream( magma_tally2_queue_t *stream )
{
    *stream = magma_tally2_stream;
    return CUBLAS_STATUS_SUCCESS;
}
