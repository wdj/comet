//-----------------------------------------------------------------------------
/*!
 * \file   tc_solve.i.hh
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  CUDA code, tc package: gemm operation.
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

#ifndef _COMET_TC_SOLVE_I_HH_
#define _COMET_TC_SOLVE_I_HH_

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Call cublas to perform required GEMM.

template<int TC_METHOD>
static void tc_solve_impl(bool is_first, int m, int n, int k,
  void* matC, TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST(matC);
  COMET_INSIST(m >= 0 && n >= 0 && k >= 0);

  // NOTE: from https://devblogs.nvidia.com/programming-tensor-cores-cuda-9/
  // "Invoke the GEMM, ensuring k, lda, ldb, and ldc are all multiples of 8, 
  //  and m is a multiple of 4"
  // "GEMMs that do not satisfy the above rules will fall back
  //  to a non-Tensor Core implementation"
  // See also https://docs.nvidia.com/cuda/cublas/index.html#cublas-gemmEx

  // k (=nfl) is derived from padded-up npvfl (multiple of 64), so always ok.
  COMET_INSIST(k % 8 == 0 && "Failed divisibility condition for tc gemm.");
  // since I_max_dim % 4 == 0; see tc_gemm_divisibility_required()
  COMET_INSIST(m % 8 == 0 && "Failed divisibility condition for tc gemm.");
  // since nvl % 4 == 0; see tc_gemm_divisibility_required()
  COMET_INSIST(n % 8 == 0 && "Failed divisibility condition for tc gemm.");

  // Make BLAS call.

  if (env.is_compute_method_gpu()) {

    // Make accelerator BLAS call.

#   ifdef COMET_USE_ACCEL

      const typename TCSelector<TC_METHOD>::GemmOut_t alpha = 1;
      const typename TCSelector<TC_METHOD>::GemmOut_t beta = is_first ? 0 : 1;

      // GPU BLAS call.

#     ifdef COMET_USE_CUDA
        const cublasStatus_t status = cublasGemmEx(
#     else
        //int status = rocblas_gemm_ex(
        const rocblas_status status = rocblas_gemm_ex(
#     endif
        tc_bufs.accelblas_handle
#     ifdef COMET_USE_CUDA
        , CUBLAS_OP_N, CUBLAS_OP_T
#     else
        , rocblas_operation_none, rocblas_operation_transpose
#     endif
        , m, n, k
        , (void*)&alpha
        , tc_bufs.tc_buf_left, TCSelector<TC_METHOD>::gemm_type_in(), m
        , tc_bufs.tc_buf_right, TCSelector<TC_METHOD>::gemm_type_in(), n
        , (void*)&beta
        , matC, TCSelector<TC_METHOD>::gemm_type_out(), m
#     ifdef COMET_USE_HIP
        , matC, TCSelector<TC_METHOD>::gemm_type_out(), m
#     endif
        , TCSelector<TC_METHOD>::gemm_type_out()
#     ifdef COMET_USE_CUDA
        //, CUBLAS_GEMM_ALGO3_TENSOR_OP // best timing for cuda 9.1.85 transpose
        //, CUBLAS_GEMM_DFALT_TENSOR_OP // good timing for cuda 9.2.88 transpose
        , CUBLAS_GEMM_ALGO4_TENSOR_OP // best timing for cuda 9.2.88 transpose
#     else
        , rocblas_gemm_algo_standard
        , 0, 0  // solution_index, flags, workspace_size, workspace
#     endif
      );
      // TODO: use CUDA 10 autotuning capability here (later).

#     ifdef COMET_USE_CUDA
        if (CUBLAS_STATUS_SUCCESS != status) {
          // Decode error message.
          fprintf(stderr, "Error: %s\n",
                   CUBLAS_STATUS_NOT_INITIALIZED == status ?
                  "CUBLAS_STATUS_NOT_INITIALIZED" :
                   CUBLAS_STATUS_ARCH_MISMATCH == status ?
                  "CUBLAS_STATUS_ARCH_MISMATCH" :
                   CUBLAS_STATUS_NOT_SUPPORTED == status ?
                  "CUBLAS_STATUS_NOT_SUPPORTED" :
                   CUBLAS_STATUS_INVALID_VALUE == status ?
                  "CUBLAS_STATUS_INVALID_VALUE" :
                   CUBLAS_STATUS_EXECUTION_FAILED == status ?
                  "CUBLAS_STATUS_EXECUTION_FAILED" : "");
        }
        COMET_INSIST(CUBLAS_STATUS_SUCCESS == status &&
                     "Failure in call to cublasGemmEx.");
#     else
        COMET_INSIST(status == rocblas_status_success &&
                     "Failure in call to rocblas_gemm_ex.");
#     endif

#   else // COMET_USE_ACCEL

      COMET_INSIST(false && "Failure to call GEMM function.");

#   endif // COMET_USE_ACCEL

    System::accel_last_call_succeeded();

  } else { // (!env.is_compute_method_gpu()) {

#   ifdef COMET_USE_CPUBLAS

      COMET_INSIST(env.tc_eff() == TC::FP32);

      const float alpha = 1;
      const float beta = is_first ? 0 : 1;

      // Make CPU BLAS call.

      cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans,
        m, n, k, alpha, (float*)tc_bufs.tc_buf_left, m,
        (float*)tc_bufs.tc_buf_right, n, beta, (float*)matC, m);

#   else // COMET_USE_CPUBLAS

      COMET_INSIST(false && "Failure to call GEMM function.");

#   endif // COMET_USE_CPUBLAS

  } // if compute_method

  env.ops_local_inc(2 * m * (double)n * (double)k);
}

//-----------------------------------------------------------------------------
/// \brief Call to perform required GEMM.

template<int TC_METHOD>
void tc_solve_(bool is_first, int nvll, int nvl, int npvfl_thisstep,
               void* matC, TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST(matC);
  COMET_INSIST(nvll >= 0 && nvl >= 0 && nvll <= nvl);
  COMET_INSIST(npvfl_thisstep >= 0);
  COMET_INSIST(env.tc_eff() != TC::NO);

  const int nfl_thisstep = npvfl_thisstep * 64;

  const int m = 2 * nvll; // metrics array dim
  const int n = 2 * nvl; // metrics array dim
  const int k = nfl_thisstep; // vectors array (as GemmIn_t) dim

  tc_solve_impl<TC_METHOD>(is_first, m, n, k, matC, tc_bufs, env);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TC_SOLVE_I_HH_

//-----------------------------------------------------------------------------
