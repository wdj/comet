//-----------------------------------------------------------------------------
/*!
 * \file   tc.hh
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  CUDA code, primarily for using tensor cores.
 * \note   Copyright (C) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_tc_hh_
#define _comet_tc_hh_

#if defined COMET_USE_CUDA
#include "cublas_v2.h"
#elif defined COMET_USE_HIP
#include "rocblas.h"
#endif

#include "env.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

#if defined COMET_USE_CUDA
//typedef cublasStatus_t accelblasStatus_t;
typedef cublasHandle_t accelblasHandle_t;
#elif defined COMET_USE_HIP
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

void tc_gemm_start(
  int m, int n, int k,
  const void* matA1, int ldda1, const void* matA2, int ldda2,
  const void* matB, int lddb, void* matC, int lddc,
  int step_2way, TCBufs& tc_bufs, CEnv& env);

void tc_bufs_malloc(
  int num_vector_local,
  int num_field_local,
  int num_packedval_field_local,
  TCBufs& tc_bufs,
  CEnv& env);

void tc_bufs_free(TCBufs& tc_bufs, CEnv& env);

size_t tc_gemm_divisibility_required(const CEnv& env);
size_t tc_gemm_size_required(size_t size_requested, const CEnv& env);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_tc_hh_

//-----------------------------------------------------------------------------
