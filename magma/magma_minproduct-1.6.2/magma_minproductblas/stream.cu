/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Mark Gates
*/

#include "common_magma_minproduct.h"

magma_minproduct_queue_t magma_minproduct_stream = 0;


/**
    Purpose
    -------
    magma_minproductblasSetKernelStream sets the CUDA stream that MAGMA_minproduct BLAS and
    CUBLAS (v1) routines use (unless explicitly given a stream).
    
    In a multi-threaded application, be careful to avoid race conditions
    when using this. For instance, if calls are executed in this order:
    
    @verbatim
        thread 1                            thread 2
        ------------------------------      ------------------------------
    1.  magma_minproductblasSetKernelStream( s1 )         
    2.                                      magma_minproductblasSetKernelStream( s2 )
    3.  magma_minproduct_dgemm( ... )
    4.                                      magma_minproduct_dgemm( ... )
    @endverbatim
    
    both magma_minproduct_dgemm would occur on stream s2. A lock should be used to prevent
    this, so the dgemm in thread 1 uses stream s1, and the dgemm in thread 2
    uses s2:
    
    @verbatim
        thread 1                            thread 2
        ------------------------------      ------------------------------
    1.  lock()                                  
    2.  magma_minproductblasSetKernelStream( s1 )          
    3.  magma_minproduct_dgemm( ... )                      
    4.  unlock()                                
    5.                                      lock()
    6.                                      magma_minproductblasSetKernelStream( s2 )
    7.                                      magma_minproduct_dgemm( ... )
    8.                                      unlock()
    @endverbatim
    
    Most BLAS calls in MAGMA_minproduct, such as magma_minproduct_dgemm, are asynchronous, so the lock
    will only have to wait until dgemm is queued, not until it is finished.
    
    Arguments
    ---------
    @param[in]
    stream  magma_minproduct_queue_t
            The CUDA stream.

    @ingroup magma_minproduct_util
    ********************************************************************/
extern "C"
cublasStatus_t magma_minproductblasSetKernelStream( magma_minproduct_queue_t stream )
{
    magma_minproduct_stream = stream;
    return cublasSetKernelStream( stream );
}


/**
    Purpose
    -------
    magma_minproductblasGetKernelStream gets the CUDA stream that MAGMA_minproduct BLAS
    routines use.

    Arguments
    ---------
    @param[out]
    stream  magma_minproduct_queue_t
            The CUDA stream.

    @ingroup magma_minproduct_util
    ********************************************************************/
extern "C"
cublasStatus_t magma_minproductblasGetKernelStream( magma_minproduct_queue_t *stream )
{
    *stream = magma_minproduct_stream;
    return CUBLAS_STATUS_SUCCESS;
}
