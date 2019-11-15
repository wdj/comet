//-----------------------------------------------------------------------------
/*!
 * \file   linalg_tc.hh
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  CUDA code, primarily for using tensor cores.
 * \note   Copyright (C) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_linalg_tc_hh_
#define _gm_linalg_tc_hh_

#if defined USE_CUDA
#include "cublas_v2.h"
#elif defined USE_HIP
#include "rocblas.h"
#endif

#include "env.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

#if defined USE_CUDA
//typedef cublasStatus_t accelblasStatus_t;
typedef cublasHandle_t accelblasHandle_t;
#elif defined USE_HIP
//typedef rocblas_status accelblasStatus_t;
typedef rocblas_handle accelblasHandle_t;
#else
//typedef int accelblasStatus_t;
typedef int accelblasHandle_t;
#endif

struct TCBufs {
  size_t tc_buf_size;
  void* tc_buf_left;
  void* tc_buf_right;
  accelblasHandle_t accelblas_handle;
};

//-----------------------------------------------------------------------------

void gm_tc_gemm_start(int m, int n, int k,
                      void* dA, int ldda,
                      void* dB, int lddb,
                      void* dC, int lddc,
                      TCBufs& tc_bufs,
                      GMEnv* env);

void gm_tc_bufs_malloc(int num_vector_local,
                       int num_field_local,
                       int num_packedval_field_local,
                       TCBufs& tc_bufs,
                       GMEnv* env);

void gm_tc_bufs_free(TCBufs& tc_bufs,
                     GMEnv* env);

size_t gm_gemm_divisibility_required(GMEnv* const env);
size_t gm_gemm_size_required(size_t size_requested, GMEnv* const env);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _gm_linalg_tc_hh_

//-----------------------------------------------------------------------------