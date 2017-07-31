#ifndef ERROR_H
#define ERROR_H

#include "common_magma_tally2.h"

// cuGetErrorString requires cuda.h, which we don't include elsewhere
// since we don't use the CUDA driver routines, only the CUDA runtime routines.
#include <cuda.h>

// overloaded C++ functions to deal with errors
void magma_tally2_xerror( cudaError_t    err, const char* func, const char* file, int line );
void magma_tally2_xerror( CUresult       err, const char* func, const char* file, int line );
void magma_tally2_xerror( cublasStatus_t err, const char* func, const char* file, int line );
void magma_tally2_xerror( magma_tally2_int_t    err, const char* func, const char* file, int line );

#ifdef __cplusplus
extern "C" {
#endif

// cuda provides cudaGetErrorString,
// but not cuGetErrorString or cublasGetErrorString, so provide our own.
// In magma_tally2.h, we also provide magma_tally2_strerror.
const char* magma_tally2_cuGetErrorString( CUresult error );
const char* magma_tally2_cublasGetErrorString( cublasStatus_t error );

#ifdef __cplusplus
}
#endif

#ifdef NDEBUG
#define check_error( err )                     ((void)0)
#define check_xerror( err, func, file, line )  ((void)0)
#else
#define check_error( err )                     magma_tally2_xerror( err, __func__, __FILE__, __LINE__ )
#define check_xerror( err, func, file, line )  magma_tally2_xerror( err, func, file, line )
#endif

#endif        //  #ifndef ERROR_H
