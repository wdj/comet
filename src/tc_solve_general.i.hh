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

#ifndef _COMET_TC_SOLVE_GENERAL_I_HH_
#define _COMET_TC_SOLVE_GENERAL_I_HH_

#include "cstdlib"
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include <mma.h>

#include "tc_solve_cutlass_general.i.hh"

using namespace nvcuda;

// Tensor core GEMM defines
// 1-bit int/int tensor core blocks
//#define WMMA1B_M   8
//#define WMMA1B_N   8
//#define WMMA1B_K 128

// Number of bits in a uint8_t
//#define NBITS      8

//=============================================================================

namespace comet {

#if 0
//-----------------------------------------------------------------------------
/// \brief 1-bit xor gemm, cpu version (not high performance).

void b1_xor_gemm_cpu(size_t m, size_t n, size_t k,
  uint8_t* a, uint8_t* b, bool beta, int32_t* c) {
  COMET_INSIST(a && b && c);
  //COMET_INSIST(m == n);

  for (size_t ind_i = 0; ind_i < m; ++ind_i) {
    for (size_t ind_j = 0; ind_j < n; ++ind_j) {
      for (size_t ind_k = 0; ind_k < k; ++ind_k) {

        const uint8_t ai = a[ind_i+m*ind_k];
        const uint8_t bj = b[ind_j+n*ind_k];
        int32_t& ci = c[ind_i+m*ind_j];

        const int32_t v = utils::popc8(ai ^ bj);
        ci = beta ? ci + v : v;
      }
    }
  }
}
#endif

//-----------------------------------------------------------------------------
/// \brief 1-bit xor gemm, gpu version (not high performance).

__global__
void b1_xor_gemm_gpu(size_t m, size_t n, size_t k,
  uint8_t* a, uint8_t* b, bool beta, int32_t* c) {
  //COMET_INSIST(a && b && c);
  //COMET_INSIST(m == n);

  const int ind_m = threadIdx_x_() + blockIdx_x_() * blockDim_x_();
  const int ind_n = blockIdx_y_();

  if (ind_m >= m || ind_n >= n)
    return;

  for (size_t ind_k = 0; ind_k < k; ++ind_k) {

    const uint8_t aik = a[ind_k + k*ind_m];
    const uint8_t bjk = b[ind_k + k*ind_n];

    int32_t& cij = c[ind_m + m*ind_n];

    const int32_t v = utils::popc8(aik ^ bjk);
    cij = beta || ind_k ? cij + v : v;
  }
}

//-----------------------------------------------------------------------------
/// \brief Simple WMMA tensor core 1-bit xor gemm

__global__
void b1_xor_gemm_gpu_tc_simple(size_t m, size_t n, size_t k, uint8_t* a,
                               uint8_t* b, bool beta, int32_t* c) {
  // Block and thread indices
  //int tx = threadIdx.x, ty = threadIdx.y;
  int bx = blockIdx.x, by = blockIdx.y;

  // Index of the first sub-matrix of A processed by the block
  // Index of the last sub-matrix of A processed by the block
  // Step size used to iterate through the sub-matrices of A
  // Loop over each block in a row
  int aBegin = k * WMMA1B_M * bx;
  int aStep  = WMMA1B_K/NBITS;

  // Index of the first sub-matrix of B processed by the block
  // Step size used to iterate through the sub-matrices of B
  // Loop over each block in a column
  int bBegin = k * WMMA1B_N * by;
  int bStep  = WMMA1B_K/NBITS;

  //printf("b=%d,%d t=%d,%d mnk=%u,%u,%u\n",bx,by,tx,ty,(unsigned int)m,
  //       (unsigned int)n,(unsigned int)k);

  // Declare fragments
  wmma::fragment<wmma::matrix_a, WMMA1B_M, WMMA1B_N, WMMA1B_K, wmma::experimental::precision::b1,
                 wmma::row_major> a_frag;
  wmma::fragment<wmma::matrix_b, WMMA1B_M, WMMA1B_N, WMMA1B_K, wmma::experimental::precision::b1,
                 wmma::col_major> b_frag;

  wmma::fragment<wmma::accumulator, WMMA1B_M, WMMA1B_N, WMMA1B_K, int> c_frag;
  wmma::fill_fragment(c_frag, 0);

  // Loop over all sub-matrices of A and B to compute block sub-matrix
  int nblocks=((k*NBITS)+WMMA1B_K-1)/WMMA1B_K;
  for(int block=0; block<nblocks; block++) {

    // Load the inputs
    wmma::load_matrix_sync(a_frag, a+aBegin+block*aStep, k*NBITS);
    wmma::load_matrix_sync(b_frag, b+bBegin+block*bStep, k*NBITS);

    // Perform the matrix multiplication
    wmma::bmma_sync(c_frag, a_frag, b_frag, c_frag);
  }

  // Store the output
  int cBegin = n*WMMA1B_M*bx + WMMA1B_N*by;
  wmma::store_matrix_sync(c+cBegin, c_frag, n, wmma::mem_row_major);

  // Print individual c values
  //__syncthreads();
  //int cShift = tx*n+ty;
  //int cInd = cBegin + cShift;
  //printf("b=%d,%d t=%d,%d mnk=%u,%u,%u cBegin=%d cShift=%d cInd=%d val=%d\n",
  //       bx,by,tx,ty,(unsigned int)m,(unsigned int)n,(unsigned int)k,cBegin,cShift,cInd,c[cInd]);
}

//-----------------------------------------------------------------------------
/// \brief Simple WMMA tensor core 1-bit xor gemm that first loads data into
///        shared memory

__global__
void b1_xor_gemm_gpu_tc_sm(size_t m, size_t n, size_t k, uint8_t* a,
                           uint8_t *b, bool beta, int32_t* c) {
  // Block and thread indices
  int tx = threadIdx.x, ty = threadIdx.y;
  int bx = blockIdx.x, by = blockIdx.y;

  // Index of the first sub-matrix of A processed by the block
  // Index of the last sub-matrix of A processed by the block
  // Step size used to iterate through the sub-matrices of A
  // Loop over each block in a row
  int aBegin = k * WMMA1B_M * bx;
  int aStep  = WMMA1B_K/NBITS;

  // Index of the first sub-matrix of B processed by the block
  // Step size used to iterate through the sub-matrices of B
  // Loop over each block in a column
  int bBegin = k * WMMA1B_N * by;
  int bStep  = WMMA1B_K/NBITS;

  // Declare fragments
  wmma::fragment<wmma::matrix_a, WMMA1B_M, WMMA1B_N, WMMA1B_K, wmma::experimental::precision::b1,
                 wmma::row_major> a_frag;
  wmma::fragment<wmma::matrix_b, WMMA1B_M, WMMA1B_N, WMMA1B_K, wmma::experimental::precision::b1,
                 wmma::col_major> b_frag;
  wmma::fragment<wmma::accumulator, WMMA1B_M, WMMA1B_N, WMMA1B_K, int> c_frag;
  wmma::fill_fragment(c_frag, 0);

  // Loop over all sub-matrices of A and B to compute block sub-matrix
  int nblocks=((k*NBITS)+WMMA1B_K-1)/WMMA1B_K;
  for(int block=0; block<nblocks; block++) {

    // Load into shared memory
    __shared__ uint8_t As[WMMA1B_M][WMMA1B_K/NBITS];
    __shared__ uint8_t Bs[WMMA1B_N][WMMA1B_K/NBITS];
    for(int i=0; i<2; i++) As[tx][ty*2+i] = a[aBegin+block*aStep+tx*k+ty*2+i];
    for(int i=0; i<2; i++) Bs[ty][tx*2+i] = b[bBegin+block*bStep+ty*k+tx*2+i];
    __syncthreads();

    // Load the inputs
    wmma::load_matrix_sync(a_frag, *As, WMMA1B_K);
    wmma::load_matrix_sync(b_frag, *Bs, WMMA1B_K);

    // Perform the matrix multiplication
    wmma::bmma_sync(c_frag, a_frag, b_frag, c_frag);
    __syncthreads();
  }

  // Store the output
  int cBegin = n*WMMA1B_M*bx + WMMA1B_N*by;
  wmma::store_matrix_sync(c+cBegin, c_frag, n, wmma::mem_row_major);
}

//-----------------------------------------------------------------------------
/// \brief Perform required GEMM.

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

  if(env.print_details()) printf("In tc_solve_impl mnk=%d,%d,%d\n",m,n,k);
  double tbegin = env.synced_time();

  // Make BLAS call.

  const bool is_timing_gemm = false; // true;

  if (is_timing_gemm)
    env.stream_synchronize(env.stream_compute());
  double t1 = !is_timing_gemm ? 0 : System::time();

  if (env.is_compute_method_gpu() && TC_METHOD == TC::B1) {

#   ifdef COMET_USE_ACCEL

      COMET_INSIST(TCTraits<TC_METHOD>::IS_B_FIELD_MAJOR);

      // Original test 1-bit GEMM kernel
      if(env.num_kernel()==0) {
        // Original b1 xor GEMM
        if(env.print_details()) printf("Launching b1_xor_gemm_gpu\n");
        enum {NUM_FL_PER_PVFL = 64};
        COMET_INSIST(k % NUM_FL_PER_PVFL == 0 && "Failed divisibility condition for tc gemm.");

        const bool beta = is_first ? 0 : 1;

        // 8 == number of uint8_t values used to store NUM_FL_PER_PVFL fields
        // in the tc buf.
        enum {BYTES_PER_PVFL_FIELDS = 8};

        const int bytes_per_gi = sizeof(typename TCTraits<TC_METHOD>::GemmIn_t);
        const size_t k_eff = (k / NUM_FL_PER_PVFL) *
                             (BYTES_PER_PVFL_FIELDS / bytes_per_gi);

        const int threadblocksize = 256;
        COMET_INSIST((threadblocksize <= 256 || ! BuildHas::HIP) &&
                     "Current HIP limitation.");
        const int num_threadblocks_0 = utils::ceil(m, threadblocksize);
        const int num_threadblocks_1 = n;

        COMET_LAUNCH_KERNEL(b1_xor_gemm_gpu,
          dim3(num_threadblocks_0, num_threadblocks_1, 1),
          dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
          m, n, k_eff, (uint8_t*)tc_bufs.tc_buf_left,
          (uint8_t*)tc_bufs.tc_buf_right, beta, (int32_t*)matC);

        System::accel_last_call_succeeded();
        env.ops_local_inc(2 * m * (double)n * (double)k_eff);
      }
      // Call WMMA 1-bit GEMM kernels
      else if(env.num_kernel()>=1 && env.num_kernel()<20) {
        enum {NUM_FL_PER_PVFL = 64};
        COMET_INSIST(k % NUM_FL_PER_PVFL == 0 && "Failed divisibility condition for tc gemm.");

        const bool beta = is_first ? 0 : 1;

        // 8 == number of uint8_t values used to store NUM_FL_PER_PVFL fields
        // in the tc buf.
        enum {BYTES_PER_PVFL_FIELDS = 8};

        const int bytes_per_gi = sizeof(typename TCTraits<TC_METHOD>::GemmIn_t);
        const size_t k_eff = (k / NUM_FL_PER_PVFL) *
                             (BYTES_PER_PVFL_FIELDS / bytes_per_gi);

        const int threadblockx = 8, threadblocky = 8;

        int gridblockx = m/threadblockx;
        int gridblocky = n/threadblocky;

        if(env.print_details()) 
          printf("Launching 1-bit general GEMM kernel m=%d n=%d k_eff=%zu k=%d beta=%d "
                 "bytes_per_gi=%d NUM_FL_PER_PVFL=%d gridDim=%d,%d threadDim=%d,%d\n",
                 m,n,k_eff,k,(int)beta,bytes_per_gi,NUM_FL_PER_PVFL,
                 gridblockx,gridblocky,threadblockx,threadblocky);

        switch(env.num_kernel()) {
          // Basic GEMM
          case 1: {
            if(env.print_details()) printf("Using simple tensor core kernel\n");
            COMET_LAUNCH_KERNEL(b1_xor_gemm_gpu_tc_simple,
              dim3(gridblockx, gridblocky, 1),
              dim3(threadblockx, threadblocky, 1), 0, env.stream_compute(),
              n, m, k_eff, (uint8_t*)tc_bufs.tc_buf_right,
              (uint8_t*)tc_bufs.tc_buf_left, beta, (int32_t*)matC);
          } break;
          // Simple shared memory GEMM
          case 2: {
            if(env.print_details()) printf("Using shared memory tensor core kernel\n");
            COMET_LAUNCH_KERNEL(b1_xor_gemm_gpu_tc_sm,
              dim3(gridblockx, gridblocky, 1),
              dim3(threadblockx, threadblocky, 1), 0, env.stream_compute(),
              n, m, k_eff, (uint8_t*)tc_bufs.tc_buf_right,
              (uint8_t*)tc_bufs.tc_buf_left, beta, (int32_t*)matC);
          } break;
          // Cutlass kernels
          case 10: {
            if(env.print_details()) printf("Using Cutlass kernel 256x128\n");
            CutlassTCGemm1B_256x128(n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k,
              (uint8_t*)tc_bufs.tc_buf_left, k, (int32_t*)matC, n);
          } break;
          case 11: {
            if(env.print_details()) printf("Using Cutlass kernel 128x256\n");
            CutlassTCGemm1B_128x256(n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k,
              (uint8_t*)tc_bufs.tc_buf_left, k, (int32_t*)matC, n);
          } break;
          case 12: {
            if(env.print_details()) printf("Using Cutlass kernel 128x128\n");
            CutlassTCGemm1B_128x128(n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k,
              (uint8_t*)tc_bufs.tc_buf_left, k, (int32_t*)matC, n);
          } break;
          case 13: {
            if(env.print_details()) printf("Using Cutlass kernel 128x64\n");
            CutlassTCGemm1B_128x64(n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k,
              (uint8_t*)tc_bufs.tc_buf_left, k, (int32_t*)matC, n);
          } break;
          case 14: {
            if(env.print_details()) printf("Using Cutlass kernel 64x128\n");
            CutlassTCGemm1B_64x128(n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k,
              (uint8_t*)tc_bufs.tc_buf_left, k, (int32_t*)matC, n);
          } break;
          case 15: {
            if(env.print_details()) printf("Using Cutlass kernel 64x64\n");
            CutlassTCGemm1B_64x64(n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k,
              (uint8_t*)tc_bufs.tc_buf_left, k, (int32_t*)matC, n);
          } break;
          case 16: {
            if(env.print_details()) printf("Using Cutlass WMMA kernel 64x64\n");
            CutlassTCGemm1BWmma_64x64(n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k,
              (uint8_t*)tc_bufs.tc_buf_left, k, (int32_t*)matC, n);
          } break;
          default: {
            printf("Failed to call appropriate 1-bit GEMM kernel for num_kernel=%d\n",
               env.num_kernel());
            COMET_INSIST(false && "Failure to call GEMM function.");
          }
        }
        System::accel_last_call_succeeded();
        env.ops_local_inc(2 * m * (double)n * (double)k_eff);
      }
      // Failed to call appropriate kernel number
      else {
        printf("Failed to call appropriate 1-bit GEMM kernel for num_kernel=%d\n",
               env.num_kernel());
        COMET_INSIST(false && "Failure to call GEMM function.");
      }
#   else // COMET_USE_ACCEL

      COMET_INSIST(false && "Failure to call GEMM function.");

#   endif // COMET_USE_ACCEL

  } else if (env.is_compute_method_gpu()) { // && TC_METHOD != TC::B1

    // Make accelerator BLAS call.

#   ifdef COMET_USE_ACCEL

      const typename TCTraits<TC_METHOD>::GemmOut_t alpha = 1;
      const typename TCTraits<TC_METHOD>::GemmOut_t beta = is_first ? 0 : 1;

      enum {IS_B_FIELD_MAJOR = TCTraits<TC_METHOD>::IS_B_FIELD_MAJOR};

      // GPU BLAS call.

#     ifdef COMET_USE_CUDA
        const cublasStatus_t status = cublasGemmEx(
#     else
        //int status = rocblas_gemm_ex(
        const rocblas_status status = rocblas_gemm_ex(
#     endif
        tc_bufs.accelblas_handle
#     ifdef COMET_USE_CUDA
        , IS_B_FIELD_MAJOR ? CUBLAS_OP_T : CUBLAS_OP_N
        , IS_B_FIELD_MAJOR ? CUBLAS_OP_N : CUBLAS_OP_T
#     else
        , IS_B_FIELD_MAJOR ? rocblas_operation_transpose : rocblas_operation_none
        , IS_B_FIELD_MAJOR ? rocblas_operation_none : rocblas_operation_transpose
#     endif
        , m, n, k
        , (void*)&alpha
        , tc_bufs.tc_buf_left, TCTraits<TC_METHOD>::gemm_type_in(), m
        , tc_bufs.tc_buf_right, TCTraits<TC_METHOD>::gemm_type_in(), n
        , (void*)&beta
        , matC, TCTraits<TC_METHOD>::gemm_type_out(), m
#     ifdef COMET_USE_HIP
        , matC, TCTraits<TC_METHOD>::gemm_type_out(), m
#     endif
        , TCTraits<TC_METHOD>::gemm_type_out()
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
        if (CUBLAS_STATUS_SUCCESS != status)
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
        COMET_INSIST(CUBLAS_STATUS_SUCCESS == status &&
                     "Failure in call to cublasGemmEx.");
#     else
        if (status != rocblas_status_success)
          // Decode error message.
          fprintf(stderr, "Error: %s\n",
                  rocblas_status_invalid_handle      == status ?
                  "handle not initialized, invalid or null" :
                  rocblas_status_not_implemented     == status ?
                  "function is not implemented" :
                  rocblas_status_invalid_pointer     == status ?
                  "invalid pointer argument" :
                  rocblas_status_invalid_size        == status ?
                  "invalid size argument" :
                  rocblas_status_memory_error        == status ?
                  "failed internal memory allocation, copy or dealloc" :
                  rocblas_status_internal_error      == status ?
                  "other internal library failure" :
                  rocblas_status_perf_degraded       == status ?
                  "performance degraded due to low device memory" :
                  rocblas_status_size_query_mismatch == status ?
                  "unmatched start/stop size query" :
                  rocblas_status_size_increased      == status ?
                  "queried device memory size increased" :
                  rocblas_status_size_unchanged      == status ?
                  "queried device memory size unchanged" :
                  rocblas_status_invalid_value       == status ?
                  "passed argument not valid" :
                  rocblas_status_continue            == status ?
                  "nothing preventing function to proceed" : "");
        COMET_INSIST(status == rocblas_status_success &&
                     "Failure in call to rocblas_gemm_ex.");
#     endif

#   else // COMET_USE_ACCEL

      COMET_INSIST(false && "Failure to call GEMM function.");

#   endif // COMET_USE_ACCEL

    if (! BuildHas::HIP) // FIX - this is a bug workaround
    COMET_INSIST(System::accel_last_call_succeeded());
    //System::accel_last_call_succeeded();

  } else { // (!env.is_compute_method_gpu()) {

    if (env.tc_eff() == TC::FP32) {

#     ifdef COMET_USE_CPUBLAS

        const float alpha = 1;
        const float beta = is_first ? 0 : 1;

        // Make CPU BLAS call.

        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans,
          m, n, k, alpha, (float*)tc_bufs.tc_buf_left, m,
          (float*)tc_bufs.tc_buf_right, n, beta, (float*)matC, m);

#     else // COMET_USE_CPUBLAS

        COMET_INSIST(false && "Failure to call GEMM function.");

#     endif // COMET_USE_CPUBLAS

#if 0
    } else if (env.tc_eff() == TC::B1) {

      enum {NUM_FL_PER_PVFL = 64};
      COMET_INSIST(k % NUM_FL_PER_PVFL == 0 && "Failed divisibility condition for tc gemm.");

      const bool beta = is_first ? 0 : 1;

      const int size_gi = sizeof(typename TCTraits<TC_METHOD>::GemmIn_t);
      const size_t k_eff = (k / NUM_FL_PER_PVFL) * (8 / size_gi);

      b1_xor_gemm_cpu(m, n, k_eff, (uint8_t*)tc_bufs.tc_buf_left,
        (uint8_t*)tc_bufs.tc_buf_right, beta, (int32_t*)matC);
#endif

    } else { // if env.tc_eff()

      COMET_INSIST(false && "Failure to call GEMM function.");

    } // if env.tc_eff()

  } // if compute_method

  env.gemmtime_inc(env.synced_time() - tbegin);

  if (is_timing_gemm) {
    env.stream_synchronize(env.stream_compute());
    double t2 = System::time();
    const double t = t2 - t1;
    printf("%i %i %i   time %f TOP/s %f\n", (int)m, (int)n, (int)k, t,
           (2 * m * (double)n * (double)k) / (t * 1e12));
  }
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

  if(env.print_details())
    printf("Calling tc_solve_impl with m=2*nvll=2*%d=%d n=2*nvl=2*%d=%d "
           "k=npvfl_thisstep*64=%d*64=%d\n",nvll,m,nvl,n,npvfl_thisstep,k);
  tc_solve_impl<TC_METHOD>(is_first, m, n, k, matC, tc_bufs, env);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TC_SOLVE_I_HH_

//-----------------------------------------------------------------------------
