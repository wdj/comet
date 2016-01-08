/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Mark Gates
*/

#include <cuda.h>  // for CUDA_VERSION

#include "common_magma_tally4sparse.h"

/**
    Maps a cuSPARSE error to a MAGMA_tally4 error.
    
    @param[in]
    status      cuSPARSE error
    
    @return MAGMA_tally4 error
    
    @ingroup magma_tally4_util
    ********************************************************************/
extern "C" magma_tally4_int_t
cusparse2magma_tally4_error( cusparseStatus_t status )
{
    switch( status ) {
        case CUSPARSE_STATUS_SUCCESS:                   return MAGMA_tally4_SUCCESS;                                break;
        case CUSPARSE_STATUS_NOT_INITIALIZED:           return MAGMA_tally4_ERR_CUSPARSE_NOT_INITIALIZED;           break;
        case CUSPARSE_STATUS_ALLOC_FAILED:              return MAGMA_tally4_ERR_CUSPARSE_ALLOC_FAILED;              break;
        case CUSPARSE_STATUS_INVALID_VALUE:             return MAGMA_tally4_ERR_CUSPARSE_INVALID_VALUE;             break;
        case CUSPARSE_STATUS_ARCH_MISMATCH:             return MAGMA_tally4_ERR_CUSPARSE_ARCH_MISMATCH;             break;
        case CUSPARSE_STATUS_MAPPING_ERROR:             return MAGMA_tally4_ERR_CUSPARSE_MAPPING_ERROR;             break;
        case CUSPARSE_STATUS_EXECUTION_FAILED:          return MAGMA_tally4_ERR_CUSPARSE_EXECUTION_FAILED;          break;
        case CUSPARSE_STATUS_INTERNAL_ERROR:            return MAGMA_tally4_ERR_CUSPARSE_INTERNAL_ERROR;            break;
        case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED: return MAGMA_tally4_ERR_CUSPARSE_MATRIX_TYPE_NOT_SUPPORTED; break;
        
        // added in CUDA 6.0
        #if CUDA_VERSION >= 6000
        case CUSPARSE_STATUS_ZERO_PIVOT:                return MAGMA_tally4_ERR_CUSPARSE_ZERO_PIVOT;                break;
        #endif
        
        default:
            return MAGMA_tally4_ERR_UNKNOWN;
            break;
    }
}
