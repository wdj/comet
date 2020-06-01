//-----------------------------------------------------------------------------
/*!
 * \file   tc.hh
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  CUDA code, primarily for using tensor cores.
 */
//-----------------------------------------------------------------------------
/*-----------------------------------------------------------------------------

Copyright 2020, UT-Battelle, LLC

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-----------------------------------------------------------------------------*/

#ifndef _COMET_TC_HH_
#define _COMET_TC_HH_

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
  GMFloat* sums_I, GMFloat* sums_J, GMFloat* sums_K,
  GMFloat* counts_I, GMFloat* counts_J, GMFloat* counts_K, int J,
  int nfal, int step_2way, TCBufs& tc_bufs, CEnv& env);

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

#endif // _COMET_TC_HH_

//-----------------------------------------------------------------------------
