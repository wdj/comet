//-----------------------------------------------------------------------------
/*!
 * \file   tc.cc
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

#include "cstdint"
#include "cstdio"

#if defined COMET_USE_CUDA
#  include "cublas_v2.h"
#  include "cuda_fp16.h"
#elif defined COMET_USE_HIP
//#  include "hip/hip_runtime_api.h"
//#  pragma GCC diagnostic ignored "-Wc99-designator"
//#  include "hip/hip_runtime.h"
#  include "rocblas.h"
#endif

#if defined COMET_USE_CPUBLAS
#if defined COMET_USE_HIP
#include "blis.h"
#else
#include BLAS_H
#endif
#endif

#include "env.hh"

#include "tc.hh"
#include "tc_helpers.i.hh"
#include "tc_in.i.hh"
#include "tc_solve.i.hh"
#include "tc_out.i.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Set matrix to zero, possibly asynchronously.
//
// This is needed because GEMM functions with beta=0 can fail to initialize
// C to zero if C is filled with trash (e.g. NaNs).

template<int TC_METHOD>
static void set_matrix_zero_start(void* matC, int lddc, int m, CEnv& env) {
  COMET_INSIST(matC);

  typedef typename TCTraits<TC_METHOD>::GemmOut_t GemmOut_t;

  const size_t size = lddc * (size_t)m * sizeof(GemmOut_t) * 4;

  if (env.is_compute_method_gpu()) {
#   ifdef COMET_USE_CUDA
      cudaMemsetAsync(matC, 0, size, env.stream_compute());
#   endif
#   ifdef COMET_USE_HIP
      hipMemsetAsync(matC, 0, size, env.stream_compute());
#    endif
    COMET_INSIST(System::accel_last_call_succeeded());
  } else {
    memset(matC, 0, size);
  }
}

//-----------------------------------------------------------------------------
/// \brief Use a standard GEMM to compute bitwise result: implementation.
///
///        This is the main function to perform the relevant
///        bitwise modified GEMM operation by use of standard GEMM
///        computations, typically using reduced precision arithmetic
///        and associated hardware features.
///
///        This is composed of three steps:
///        1. copy the input matrices into the required matrix format
///        2. apply the GEMM
///        3. adjust the results in-place to the required format.
///        To save on memory, this 3-step process is broken into
///        a sequence of steps as an outer loop.
///        All of these operations are pipelined in a (CUDA) execution
///        stream.
///

template<int TC_METHOD>
static void tc_gemm_start_impl_(
  int m, int n, int k,
  const void* matA1, const void* matA2, const void* matB, void* matC, int lddc,
  GMFloat* sums_I, GMFloat* sums_J, GMFloat* sums_K,
  GMFloat* counts_I, GMFloat* counts_J, GMFloat* counts_K, int J,
  TCBufs& tc_bufs, int nfal, int step_2way, CEnv& env) {

  const int nvl = n;
  const int npvfl = k;
  const int I_max = m;
  const int I_max_dim = lddc;
  COMET_INSIST(I_max <= I_max_dim && I_max_dim <= nvl);
  // nvll is the effective nvl (column dim) for the left matrix
  // We only really only need up to I_max, but must compute to I_max_dim
  // to satisfy cublas divisibility requirements.
  // Note nvl is always the column dim for the right matrix (CHECK).
  const int nvll = I_max_dim;
  COMET_INSIST((size_t)nvll == tc_gemm_size_required(nvll, env));

  // Get matX counts if needed.

  tc_compute_matX_counts(I_max, I_max_dim, nvl, npvfl, nfal,
    (uint32_t*)matA1, (uint32_t*)matA2, tc_bufs, step_2way, env);

  set_matrix_zero_start<TC_METHOD>(matC, lddc, m, env);

  const int num_tc_steps = env.num_tc_steps();

  // Loop over steps of algorithm.
  for (int tc_step_num = 0; tc_step_num < num_tc_steps; ++tc_step_num) {

    // Select the block row of the left and right matrices for this step.
    const int pvfl_min = ((tc_step_num+0) * npvfl) / num_tc_steps;
    const int pvfl_max = ((tc_step_num+1) * npvfl) / num_tc_steps;
    const int npvfl_thisstep = pvfl_max - pvfl_min;

    const bool is_empty_block_row = 0 == npvfl_thisstep;
    if (is_empty_block_row)
      continue;

    // Convert the input matrices of packed bit values into matrices
    // of a type suitable for the GEMM.
    const bool left_matrix = false; // A
    const bool right_matrix = true; // B
    tc_buf_write_<TC_METHOD>(left_matrix, I_max, I_max_dim, nvl, npvfl,
      npvfl_thisstep, pvfl_min, nfal, (uint32_t*)matA1, (uint32_t*)matA2,
      tc_bufs, step_2way, env);
    tc_buf_write_<TC_METHOD>(right_matrix, I_max, I_max_dim, nvl, npvfl,
      npvfl_thisstep, pvfl_min, nfal, (uint32_t*)matB, (uint32_t*)matB,
      tc_bufs, step_2way, env);

    // Perform the GEMM for this pair of block rows; accumulate.
    const bool is_first = 0 == pvfl_min;
    tc_solve_<TC_METHOD>(is_first, nvll, nvl, npvfl_thisstep,
      matC, tc_bufs, env);

  } // for

  // Postprocess GEMM results.

  if (env.is_threshold_tc()) {

    tc_out_<TC_METHOD, MetricFormat::SINGLE>(nvll, nvl, matC,
      sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, J, step_2way, env);

  } else {

    tc_out_<TC_METHOD, MetricFormat::PACKED_DOUBLE>(nvll, nvl, matC,
      sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, J, step_2way, env);

  }
}

//=============================================================================
// EXTERNALLY VISIBLE FUNCTIONS: GENERAL
//=============================================================================

//-----------------------------------------------------------------------------
/// \brief Use a standard GEMM to compute CoMet metrics bitwise result.

void tc_gemm_start(
  int m, int n, int k,
  const void* matA1, int ldda1, const void* matA2, int ldda2,
  const void* matB, int lddb, void* matC, int lddc,
  GMFloat* sums_I, GMFloat* sums_J, GMFloat* sums_K,
  GMFloat* counts_I, GMFloat* counts_J, GMFloat* counts_K, int J,
  int nfal, int step_2way, TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST(matA1 && matA2 && matB && matC);
  COMET_INSIST(m >= 0 && n >= 0 && k >= 0);
  COMET_INSIST(ldda1 >= 0 && ldda2 >= 0 && lddb >= 0 && lddc >= 0);
  COMET_INSIST(k <= ldda1 && k <= ldda2 && k <= lddb && m <= lddc);
  COMET_INSIST(env.tc_eff() != TC::NO);
  COMET_INSIST(env.is_metric_type_bitwise());
  COMET_INSIST(nfal <= 64 * k);

  COMET_INSIST(tc_bufs.tc_buf_left);

  COMET_INSIST(ldda1 == k && ldda2 == k && lddb == k); // always true here

  // Select required template function instance.

  switch (env.tc_eff()) {
    // --------------
    case TC::FP32: {
      tc_gemm_start_impl_<TC::FP32>(
        m, n, k, matA1, matA2, matB, matC, lddc,
        sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, J,
        tc_bufs, nfal, step_2way, env);
    } break;
    // --------------
    case TC::FP16: {
      tc_gemm_start_impl_<TC::FP16>(
        m, n, k, matA1, matA2, matB, matC, lddc,
        sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, J,
        tc_bufs, nfal, step_2way, env);
    } break;
    // --------------
    case TC::INT8: {
      tc_gemm_start_impl_<TC::INT8>(
        m, n, k, matA1, matA2, matB, matC, lddc,
        sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, J,
        tc_bufs, nfal, step_2way, env);
    } break;
    // --------------
    case TC::B1: {
      tc_gemm_start_impl_<TC::B1>(
        m, n, k, matA1, matA2, matB, matC, lddc,
        sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, J,
        tc_bufs, nfal, step_2way, env);
    } break;
    // --------------
    default:
      COMET_INSIST(false && "Invalid tc type.");
  } // switch
}

//-----------------------------------------------------------------------------
/// \brief Divisibility requirement for GEMM.

size_t tc_gemm_divisibility_required(const CEnv& env) {

  const bool need_divisible_by_4 = env.tc_eff() != TC::NO;

  return need_divisible_by_4 ? 4 : 1;
}

//-----------------------------------------------------------------------------
/// \brief Size requirement for GEMM.

size_t tc_gemm_size_required(size_t size_requested, const CEnv& env) {

  const size_t factor = tc_gemm_divisibility_required(env);

  return utils::ceil(size_requested, factor)*factor;
}

//=============================================================================
// EXTERNALLY VISIBLE FUNCTIONS: BUFFER MANAGEMENT
//=============================================================================

//-----------------------------------------------------------------------------
/// \brief Initialize TCBufs object by allocating memory etc.

void tc_bufs_malloc(int num_vector_local,
                    int num_field_local,
                    int num_packedval_field_local,
                    TCBufs& tc_bufs,
                    CEnv& env) {
  COMET_INSIST(num_vector_local >= 0);
  COMET_INSIST(num_packedval_field_local >= 0);
  COMET_INSIST(!tc_bufs.tc_buf_left);
  COMET_INSIST(!tc_bufs.tc_buf_right);

  if (!env.is_metric_type_bitwise() || env.tc_eff() == TC::NO)
    return;

  // Calculate sizes.

  const size_t nvl = num_vector_local;
  const size_t npvfl = num_packedval_field_local;
  const size_t npvfl_thisstep_max = utils::ceil(npvfl, (size_t)env.num_tc_steps());

  const int sizeof_gemm_in_t =
     env.tc_eff() == TC::FP32 ?
       sizeof(typename TCTraits<TC::FP32>::GemmIn_t) :
     env.tc_eff() == TC::FP16 ?
       sizeof(typename TCTraits<TC::FP16>::GemmIn_t) :
     env.tc_eff() == TC::INT8 ?
       sizeof(typename TCTraits<TC::INT8>::GemmIn_t) :
     env.tc_eff() == TC::B1 ?
       sizeof(typename TCTraits<TC::B1>::GemmIn_t) :
     0;
  COMET_INSIST(TC::is_valid(env.tc_eff())); // this code must be updated if new method

  const size_t nvlX2 = nvl * 2;

  tc_bufs.tc_buf_size = nvlX2 * (npvfl_thisstep_max * 64) * sizeof_gemm_in_t;
  tc_bufs.tc_buf_size = tc_bufs.tc_buf_size ? tc_bufs.tc_buf_size : 1;

  tc_bufs.matX_counts_size = env.is_using_xor() && env.num_way() == NumWay::_3 ?
    nvlX2 * sizeof(uint32_t) : 1;

  if (env.is_compute_method_gpu()) {

    // Allocate buffers.

#if defined COMET_USE_CUDA
    cudaMalloc(&tc_bufs.tc_buf_left, tc_bufs.tc_buf_size);
#elif defined COMET_USE_HIP
    hipMalloc(&tc_bufs.tc_buf_left, tc_bufs.tc_buf_size);
#endif
    COMET_INSIST(System::accel_last_call_succeeded());
    env.gpu_mem_local_inc(tc_bufs.tc_buf_size);

#if defined COMET_USE_CUDA
    cudaMalloc(&tc_bufs.tc_buf_right, tc_bufs.tc_buf_size);
#elif defined COMET_USE_HIP
    hipMalloc(&tc_bufs.tc_buf_right, tc_bufs.tc_buf_size);
#endif
    COMET_INSIST(System::accel_last_call_succeeded());
    env.gpu_mem_local_inc(tc_bufs.tc_buf_size);

#if defined COMET_USE_CUDA
    cudaMalloc(&tc_bufs.matX_counts, tc_bufs.matX_counts_size);
#elif defined COMET_USE_HIP
    hipMalloc(&tc_bufs.matX_counts, tc_bufs.matX_counts_size);
#endif
    COMET_INSIST(System::accel_last_call_succeeded());
    env.gpu_mem_local_inc(tc_bufs.matX_counts_size);

    // Set up accel blas handle.

#if defined COMET_USE_CUDA
    cublasStatus_t status = cublasCreate(&tc_bufs.accelblas_handle);
    COMET_INSIST(status == CUBLAS_STATUS_SUCCESS && "Error in cublasCreate.");
    COMET_INSIST(System::accel_last_call_succeeded());

    status = cublasSetStream(tc_bufs.accelblas_handle, env.stream_compute());
    COMET_INSIST(status == CUBLAS_STATUS_SUCCESS && "Error in cublasSetStream.");

    status = cublasSetMathMode(tc_bufs.accelblas_handle, CUBLAS_TENSOR_OP_MATH);
    COMET_INSIST(status == CUBLAS_STATUS_SUCCESS &&
                 "Error in cublasSetMathMode.");
    COMET_INSIST(System::accel_last_call_succeeded());
#elif defined COMET_USE_HIP
    int status = rocblas_create_handle(&tc_bufs.accelblas_handle);
    COMET_INSIST(status == rocblas_status_success &&
             "Error in rocblas_create_handle.");
    COMET_INSIST(System::accel_last_call_succeeded());

    status = rocblas_set_stream(tc_bufs.accelblas_handle, env.stream_compute());
    COMET_INSIST(status == rocblas_status_success &&
             "Error in rocblas_set_stream.");
    COMET_INSIST(System::accel_last_call_succeeded());

    //FIX - will this be needed for AMD gpu?
    //  status = cublasSetMathMode(tc_bufs.accelblas_handle,
    //                             CUBLAS_TENSOR_OP_MATH);
    //  COMET_INSIST(status == CUBLAS_STATUS_SUCCESS &&
    //           "Error in cublasSetMathMode.");
#endif

  } else { // compute_method

    // Allocate buffers.

    tc_bufs.tc_buf_left = malloc(tc_bufs.tc_buf_size);
    COMET_INSIST(tc_bufs.tc_buf_left);
    env.cpu_mem_local_inc(tc_bufs.tc_buf_size);

    tc_bufs.tc_buf_right = malloc(tc_bufs.tc_buf_size);
    COMET_INSIST(tc_bufs.tc_buf_right);
    env.cpu_mem_local_inc(tc_bufs.tc_buf_size);
//memset((void*)tc_bufs.tc_buf_left, 0, tc_bufs.tc_buf_size);
//memset((void*)tc_bufs.tc_buf_right, 0, tc_bufs.tc_buf_size);

  } // compute_method
}

//-----------------------------------------------------------------------------
/// \brief Terminate TCBufs object by deallocating memory etc.

void tc_bufs_free(TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST((tc_bufs.tc_buf_left != 0) == (tc_bufs.tc_buf_right != 0));

  if (!tc_bufs.tc_buf_left) {
    return;
  }

  if (env.is_compute_method_gpu()) {

    // Free buffers.

#if defined COMET_USE_CUDA
    cudaFree(tc_bufs.tc_buf_left);
#elif defined COMET_USE_HIP
    hipFree(tc_bufs.tc_buf_left);
#endif
    COMET_INSIST(System::accel_last_call_succeeded());
    tc_bufs.tc_buf_left = NULL;
    env.gpu_mem_local_dec(tc_bufs.tc_buf_size);

#if defined COMET_USE_CUDA
    cudaFree(tc_bufs.tc_buf_right);
#elif defined COMET_USE_HIP
    hipFree(tc_bufs.tc_buf_right);
#endif
    COMET_INSIST(System::accel_last_call_succeeded());
    tc_bufs.tc_buf_right = NULL;
    env.gpu_mem_local_dec(tc_bufs.tc_buf_size);

#if defined COMET_USE_CUDA
    cudaFree(tc_bufs.matX_counts);
#elif defined COMET_USE_HIP
    hipFree(tc_bufs.matX_counts);
#endif
    COMET_INSIST(System::accel_last_call_succeeded());
    tc_bufs.matX_counts = NULL;
    env.gpu_mem_local_dec(tc_bufs.matX_counts_size);

    // Free accel blas handle.

#if defined COMET_USE_CUDA
    cublasStatus_t status = cublasDestroy(tc_bufs.accelblas_handle);
    COMET_INSIST(status == CUBLAS_STATUS_SUCCESS && "Error in cublasDestroy.");
#elif defined COMET_USE_HIP
    int status = rocblas_destroy_handle(tc_bufs.accelblas_handle);
    COMET_INSIST(status == rocblas_status_success &&
             "Error in rocblas_destroy_handle.");
#endif
    COMET_INSIST(System::accel_last_call_succeeded());

  } else { // compute_method

    // Free buffers.

    free(tc_bufs.tc_buf_left);
    tc_bufs.tc_buf_left = NULL;
    env.cpu_mem_local_dec(tc_bufs.tc_buf_size);

    free(tc_bufs.tc_buf_right);
    tc_bufs.tc_buf_right = NULL;
    env.cpu_mem_local_dec(tc_bufs.tc_buf_size);

  } // compute_method
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
