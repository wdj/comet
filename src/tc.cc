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
#  include "cuda_fp16.h"
#endif

#if defined COMET_USE_CPUBLAS
#if defined COMET_USE_HIP
//#include "blis.h"
#include "cblas.h"
#else
#include BLAS_H
#endif
#endif

#include "env.hh"
#include "histograms.hh"

#include "tc.hh"
#include "tc_helpers.i.hh"
#include "tc_in.i.hh"
#include "tc_solve.i.hh"
#include "tc_out.i.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Divisibility requirement for GEMM: vector axis.

// The units are "vector"s (as counted by num_vector_local).

size_t tc_gemm_vaxis_divisibility_required(const CEnv& env) {

  const size_t result = env.tc_eff() == TC::NO ? 1 :
    // NOTE: for tc, must be >= 2 -- see tc_in.
    //!env.is_compute_method_gpu() ? 1 :
    //!env.is_metric_type_bitwise() ? 1 :
    env.is_using_cutlass() ? 256 :
    // Curent requirement >= 4 - see tc_in.
    env.is_using_cutlass_mockup() ? 4 :
    8; // cuBLAS/rocBLAS

    //// Curent requirement >= 4 - see tc_in.
    //tc_solve_use_mockup(env) &&
    //  (env.tc_eff() == TC::B1 || env.tc_eff() == TC::INT4) ? 4 :
    //env.tc_eff() == TC::B1 ? 256 :
    //env.tc_eff() == TC::INT4 ? 256 : 8;

  return result;
}

//-----------------------------------------------------------------------------
/// \brief Divisibility requirement for GEMM: field axis.

size_t tc_gemm_faxis_divisibility_required(const CEnv& env) {

  // The units here are "packed field"s -- sizeof = sizeof(double[2]) = 16.

  const size_t result = env.tc_eff() == TC::NO ? 1 :
    //!env.is_compute_method_gpu() ? 1 :
    //!env.is_metric_type_bitwise() ? 1 :
    env.is_using_cutlass() ? 64 :
    1;

    //!(env.tc_eff() == TC::B1 || env.tc_eff() == TC::INT4) ? 1 : 4;
    //!(env.tc_eff() == TC::B1 || env.tc_eff() == TC::INT4) ? 1 : 256;
    //!(env.tc_eff() == TC::B1 || env.tc_eff() == TC::INT4) ? 1 : 64;

  return result;
}

//-----------------------------------------------------------------------------
/// \brief Divisibility requirement for num_vector_local, due to GEMMs.

size_t tc_nvl_divisibility_required_for_gemm(const CEnv& env) {

  const size_t factor_gemm = tc_gemm_vaxis_divisibility_required(env);

  // tc methods require 2 columns per vector, exploit this if possible..
  const size_t factor = env.tc_eff() == 0 ? factor_gemm :
                        factor_gemm % 2 == 0 ? factor_gemm / 2 :
                        factor_gemm;

  return factor;
}

//-----------------------------------------------------------------------------
/// \brief Size requirement for GEMM: vector axis.

size_t tc_gemm_vaxis_size_required(size_t size_requested, const CEnv& env) {

  const size_t factor = tc_gemm_vaxis_divisibility_required(env);

  // Pad up.
  return utils::ceil(size_requested, factor)*factor;
}

//-----------------------------------------------------------------------------
/// \brief Size requirement for GEMM: field axis.

size_t tc_gemm_faxis_size_required(size_t size_requested, const CEnv& env) {

  const size_t factor = tc_gemm_faxis_divisibility_required(env);

  // Pad up.
  return utils::ceil(size_requested, factor)*factor;
}

//-----------------------------------------------------------------------------
/// \brief Size requirement for num_vectror_loca, due to GEMMs.

size_t tc_nvl_size_required_for_gemm(size_t size_requested, const CEnv& env) {

  const size_t factor = tc_nvl_divisibility_required_for_gemm(env);

  // Pad up.
  return utils::ceil(size_requested, factor)*factor;
}

//-----------------------------------------------------------------------------
/// \brief Set matrix to zero, possibly asynchronously.
//
// This is needed because GEMM functions with beta=0 can fail to initialize
// C to zero if C is filled with trash (e.g. NaNs).

template<int TC_METHOD>
static void tc_set_matrix_zero_start(void* matC, int lddc, int m, CEnv& env) {
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
  int nfal, int step_2way, TCBufs& tc_bufs, Histograms& histograms,
  GemmShapes& gemm_shapes, CEnv& env, TCDebug tc_debug) {

  const int nvl = n;
  const int npfl = k;
  const int d = tc_gemm_faxis_divisibility_required(env);
  COMET_INSIST(npfl % d == 0);
  const int I_max = m;
  const int I_max_dim = lddc;
  COMET_INSIST(I_max <= I_max_dim && I_max_dim <= nvl);
  // nvll is the effective nvl (column dim) for the left matrix
  // We only really only need up to I_max, but must compute to I_max_dim
  // to satisfy cublas divisibility requirements.
  // Note nvl is always the column dim for the right matrix (CHECK).
  const int nvll = I_max_dim;
  COMET_INSIST((size_t)nvll == tc_nvl_size_required_for_gemm(nvll, env));

  // Get matX counts if needed.

  tc_compute_matX_counts(I_max, I_max_dim, nvl, npfl, nfal,
    (TCWord_t*)matA1, (TCWord_t*)matA2, tc_bufs, step_2way, env);

  tc_set_matrix_zero_start<TC_METHOD>(matC, lddc, m, env);

  const int num_tc_steps = env.num_tc_steps();

//double t1, t2;
  // Loop over steps of algorithm.
  for (int tc_step_num = 0; tc_step_num < num_tc_steps; ++tc_step_num) {

//t1 = env.synced_time();

    // Select the block row of the left and right matrices for this step.
    const int pfl_min = (((tc_step_num+0) * (npfl/d)) / num_tc_steps) * d;
    const int pfl_max = (((tc_step_num+1) * (npfl/d)) / num_tc_steps) * d;
    const int npfl_thisstep = pfl_max - pfl_min;

    const bool is_empty_block_row = 0 == npfl_thisstep;
    if (is_empty_block_row)
      continue;

    // Convert the input matrices of packed bit values into matrices
    // of a type suitable for the GEMM.
    enum {IS_LEFT = true};
    tc_buf_write_<TC_METHOD, IS_LEFT>(I_max, I_max_dim, nvl, npfl,
      npfl_thisstep, pfl_min, nfal, (TCWord_t*)matA1, (TCWord_t*)matA2,
      tc_bufs, step_2way, env);
    tc_buf_write_<TC_METHOD, !IS_LEFT>(I_max, I_max_dim, nvl, npfl,
      npfl_thisstep, pfl_min, nfal, (TCWord_t*)matB, (TCWord_t*)matB,
      tc_bufs, step_2way, env);

//t2 = env.synced_time();
//if (System::is_proc_num_0()) printf("1 %.20f %.20f\n", (t2-t1), (double)t2);
//t1 = env.synced_time();
    // Perform the GEMM for this pair of block rows; accumulate.
    const bool is_first = 0 == pfl_min;
//for (int i=0; i<10; ++i)
    tc_solve_<TC_METHOD>(is_first, nvll, nvl, npfl_thisstep,
      matC, tc_bufs, env);

//t2 = env.synced_time();
//if (System::is_proc_num_0()) printf("2 %.20f %.20f\n", (t2-t1), (double)t2);
  } // for

  // Postprocess GEMM results.

//t1 = env.synced_time();
  if (env.is_threshold_tc()) {

    tc_out_<TC_METHOD, MetricFormat::SINGLE>(nvll, nvl, matC,
      sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, tc_bufs.matX_counts,
      J, step_2way, histograms, gemm_shapes, env, tc_debug);

  } else {

    tc_out_<TC_METHOD, MetricFormat::PACKED_DOUBLE>(nvll, nvl, matC,
      sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, tc_bufs.matX_counts,
      J, step_2way, histograms, gemm_shapes, env, tc_debug);

  }
//t2 = env.synced_time();
//if (System::is_proc_num_0()) printf("3 %.20f %.20f\n", (t2-t1), (double)t2);
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
  int nfal, int step_2way, TCBufs& tc_bufs, Histograms& histograms,
  GemmShapes& gemm_shapes,
  CEnv& env, TCDebug tc_debug) {
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
        nfal, step_2way, tc_bufs, histograms, gemm_shapes, env, tc_debug);
    } break;
    // --------------
    case TC::FP16: {
      tc_gemm_start_impl_<TC::FP16>(
        m, n, k, matA1, matA2, matB, matC, lddc,
        sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, J,
        nfal, step_2way, tc_bufs, histograms, gemm_shapes, env, tc_debug);
    } break;
    // --------------
    case TC::INT8: {
      tc_gemm_start_impl_<TC::INT8>(
        m, n, k, matA1, matA2, matB, matC, lddc,
        sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, J,
        nfal, step_2way, tc_bufs, histograms, gemm_shapes, env, tc_debug);
    } break;
    // --------------
    case TC::B1: {
      tc_gemm_start_impl_<TC::B1>(
        m, n, k, matA1, matA2, matB, matC, lddc,
        sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, J,
        nfal, step_2way, tc_bufs, histograms, gemm_shapes, env, tc_debug);
    } break;
    // --------------
    case TC::INT4: {
      tc_gemm_start_impl_<TC::INT4>(
        m, n, k, matA1, matA2, matB, matC, lddc,
        sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, J,
        nfal, step_2way, tc_bufs, histograms, gemm_shapes, env, tc_debug);
    } break;
    // --------------
    default:
      COMET_INSIST(false && "Invalid tc type.");
  } // switch
}

//=============================================================================
// EXTERNALLY VISIBLE FUNCTIONS: BUFFER MANAGEMENT
//=============================================================================

//-----------------------------------------------------------------------------
/// \brief Initialize TCBufs object by allocating memory etc.

void TCBufs::malloc(int num_vector_local,
                    //int num_field_local,
                    int num_packedfield_local,
                    TCBufs& tc_bufs,
                    CEnv& env) {
  COMET_INSIST(num_vector_local >= 0);
  COMET_INSIST(num_packedfield_local >= 0);
  COMET_INSIST(!tc_bufs.tc_buf_left);
  COMET_INSIST(!tc_bufs.tc_buf_right);

  if (!env.is_metric_type_bitwise() || env.tc_eff() == TC::NO)
    return;

  // Calculate sizes.

  const size_t nvl = num_vector_local;
  const size_t npfl = num_packedfield_local;
  const int d = tc_gemm_faxis_divisibility_required(env);
  const size_t npfl_thisstep_max =
    utils::ceil(npfl/d, (size_t)env.num_tc_steps()) * d;

  const int sizeof_gemm_in_t =
     env.tc_eff() == TC::FP32 ? sizeof(typename TCTraits<TC::FP32>::GemmIn_t) :
     env.tc_eff() == TC::FP16 ? sizeof(typename TCTraits<TC::FP16>::GemmIn_t) :
     env.tc_eff() == TC::INT8 ? sizeof(typename TCTraits<TC::INT8>::GemmIn_t) :
     env.tc_eff() == TC::B1 ? sizeof(typename TCTraits<TC::B1>::GemmIn_t) :
     env.tc_eff() == TC::INT4 ? sizeof(typename TCTraits<TC::INT4>::GemmIn_t) :
     0;

  const int nfpgi =
     env.tc_eff() == TC::FP32 ? TCTraits<TC::FP32>::NFPGI :
     env.tc_eff() == TC::FP16 ? TCTraits<TC::FP16>::NFPGI :
     env.tc_eff() == TC::INT8 ? TCTraits<TC::INT8>::NFPGI :
     env.tc_eff() == TC::B1 ? TCTraits<TC::B1>::NFPGI :
     env.tc_eff() == TC::INT4 ? TCTraits<TC::INT4>::NFPGI :
     0;
  COMET_INSIST(TC::is_valid(env.tc_eff())); // must be updated if new method

  const size_t nvlX2 = nvl * 2;

  tc_bufs.tc_buf_size = nvlX2 * (npfl_thisstep_max * 64 / nfpgi) *
                        sizeof_gemm_in_t;
  tc_bufs.tc_buf_size = tc_bufs.tc_buf_size ? tc_bufs.tc_buf_size : 1;

  tc_bufs.matX_counts_size = env.is_using_xor() && env.num_way() == NumWay::_3 ?
    nvlX2 * sizeof(uint32_t) : 1;

  if (env.is_compute_method_gpu()) {

    // Allocate buffers.

#   if defined COMET_USE_CUDA
//printf("3a %zu %i\n", tc_bufs.tc_buf_size, System::proc_num());
      cudaMalloc(&tc_bufs.tc_buf_left, tc_bufs.tc_buf_size);
#   elif defined COMET_USE_HIP
      hipMalloc(&tc_bufs.tc_buf_left, tc_bufs.tc_buf_size);
#   endif
    COMET_INSIST(System::accel_last_call_succeeded());
    env.gpu_mem_local_inc(tc_bufs.tc_buf_size);

#   if defined COMET_USE_CUDA
//printf("3b %zu %i\n", tc_bufs.tc_buf_size, System::proc_num());
      cudaMalloc(&tc_bufs.tc_buf_right, tc_bufs.tc_buf_size);
#   elif defined COMET_USE_HIP
      hipMalloc(&tc_bufs.tc_buf_right, tc_bufs.tc_buf_size);
#   endif
    COMET_INSIST(System::accel_last_call_succeeded());
    env.gpu_mem_local_inc(tc_bufs.tc_buf_size);

#   if defined COMET_USE_CUDA
//printf("3c %zu %i\n", tc_bufs.matX_counts_size, System::proc_num());
      cudaMalloc(&tc_bufs.matX_counts, tc_bufs.matX_counts_size);
#   elif defined COMET_USE_HIP
      hipMalloc(&tc_bufs.matX_counts, tc_bufs.matX_counts_size);
#   endif
    COMET_INSIST(System::accel_last_call_succeeded());
    env.gpu_mem_local_inc(tc_bufs.matX_counts_size);

    // Set up accel blas handle.

#   if defined COMET_USE_CUDA
      cublasStatus_t status = cublasCreate(&tc_bufs.accelblas_handle);
      COMET_INSIST(status == CUBLAS_STATUS_SUCCESS && "Error in cublasCreate.");
      COMET_INSIST(System::accel_last_call_succeeded());

      status = cublasSetStream(tc_bufs.accelblas_handle, env.stream_compute());
      COMET_INSIST(status == CUBLAS_STATUS_SUCCESS &&
                   "Error in cublasSetStream.");

      status = cublasSetMathMode(tc_bufs.accelblas_handle,
                                 CUBLAS_TENSOR_OP_MATH);
      COMET_INSIST(status == CUBLAS_STATUS_SUCCESS &&
                   "Error in cublasSetMathMode.");
      COMET_INSIST(System::accel_last_call_succeeded());
#   elif defined COMET_USE_HIP
      int status = rocblas_create_handle(&tc_bufs.accelblas_handle);
      COMET_INSIST(status == rocblas_status_success &&
               "Error in rocblas_create_handle.");
      COMET_INSIST(System::accel_last_call_succeeded());

      status = rocblas_set_stream(tc_bufs.accelblas_handle,
                                  env.stream_compute());
      COMET_INSIST(status == rocblas_status_success &&
               "Error in rocblas_set_stream.");
      COMET_INSIST(System::accel_last_call_succeeded());

      //FIX - will this be needed for AMD gpu?
      //  status = cublasSetMathMode(tc_bufs.accelblas_handle,
      //                             CUBLAS_TENSOR_OP_MATH);
      //  COMET_INSIST(status == CUBLAS_STATUS_SUCCESS &&
      //           "Error in cublasSetMathMode.");
#   endif

  } else { // compute_method

    // Allocate buffers.

    tc_bufs.tc_buf_left = ::malloc(tc_bufs.tc_buf_size);
    COMET_INSIST(tc_bufs.tc_buf_left);
    env.cpu_mem_local_inc(tc_bufs.tc_buf_size);
    //memset((void*)tc_bufs.tc_buf_left, 0, tc_bufs.tc_buf_size);

    tc_bufs.tc_buf_right = ::malloc(tc_bufs.tc_buf_size);
    COMET_INSIST(tc_bufs.tc_buf_right);
    env.cpu_mem_local_inc(tc_bufs.tc_buf_size);
    //memset((void*)tc_bufs.tc_buf_right, 0, tc_bufs.tc_buf_size);

  } // compute_method
}

//-----------------------------------------------------------------------------
/// \brief Terminate TCBufs object by deallocating memory etc.

void TCBufs::free(TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST((tc_bufs.tc_buf_left != 0) == (tc_bufs.tc_buf_right != 0));

  if (!tc_bufs.tc_buf_left)
    return;

  if (env.is_compute_method_gpu()) {

    // Free buffers.

#   if defined COMET_USE_CUDA
      cudaFree(tc_bufs.tc_buf_left);
#   elif defined COMET_USE_HIP
      hipFree(tc_bufs.tc_buf_left);
#   endif
    COMET_INSIST(System::accel_last_call_succeeded());
    tc_bufs.tc_buf_left = NULL;
    env.gpu_mem_local_dec(tc_bufs.tc_buf_size);

#   if defined COMET_USE_CUDA
      cudaFree(tc_bufs.tc_buf_right);
#   elif defined COMET_USE_HIP
      hipFree(tc_bufs.tc_buf_right);
#   endif
    COMET_INSIST(System::accel_last_call_succeeded());
    tc_bufs.tc_buf_right = NULL;
    env.gpu_mem_local_dec(tc_bufs.tc_buf_size);

#   if defined COMET_USE_CUDA
      cudaFree(tc_bufs.matX_counts);
#   elif defined COMET_USE_HIP
      hipFree(tc_bufs.matX_counts);
#   endif
    COMET_INSIST(System::accel_last_call_succeeded());
    tc_bufs.matX_counts = NULL;
    env.gpu_mem_local_dec(tc_bufs.matX_counts_size);

    // Free accel blas handle.

#   if defined COMET_USE_CUDA
      cublasStatus_t status = cublasDestroy(tc_bufs.accelblas_handle);
      COMET_INSIST(status == CUBLAS_STATUS_SUCCESS &&
                   "Error in cublasDestroy.");
#   elif defined COMET_USE_HIP
      int status = rocblas_destroy_handle(tc_bufs.accelblas_handle);
      COMET_INSIST(status == rocblas_status_success &&
             "Error in rocblas_destroy_handle.");
#   endif
    COMET_INSIST(System::accel_last_call_succeeded());

  } else { // compute_method CPU

    // Free buffers.

    ::free(tc_bufs.tc_buf_left);
    tc_bufs.tc_buf_left = NULL;
    env.cpu_mem_local_dec(tc_bufs.tc_buf_size);

    ::free(tc_bufs.tc_buf_right);
    tc_bufs.tc_buf_right = NULL;
    env.cpu_mem_local_dec(tc_bufs.tc_buf_size);

  } // compute_method
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
