//-----------------------------------------------------------------------------
/*!
 * \file   linalg_tc.cc
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  CUDA code, primarily for using tensor cores.
 * \note   Copyright (C) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "cstdint"

#if defined USE_CUDA
#  include "cublas_v2.h"
#  include "cuda_fp16.h"
#elif defined USE_HIP
#  include "hip/hip_runtime_api.h"
#  include "hip/hip_runtime.h"
#  include "rocblas.h"
#else
#  define __device__
#  define __global__
#endif

#if defined USE_CPUBLAS
#include BLAS_H
#endif // USE_CPUBLAS

#include "env.hh"
#include "linalg_tc.hh"

//=============================================================================

namespace comet {

//=============================================================================
// HELPERS
//=============================================================================

//-----------------------------------------------------------------------------
/// \brief Abstracted thread indexing/dimensions functions.

#if defined USE_CUDA
__device__ static int threadIdx_x_() { return threadIdx.x; }

__device__ static int blockIdx_x_() { return blockIdx.x; }
__device__ static int blockIdx_y_() { return blockIdx.y; }
__device__ static int blockIdx_z_() { return blockIdx.z; }

__device__ static int blockDim_x_() { return blockDim.x; }

__device__ static int gridDim_y_() { return gridDim.y; }
#elif defined USE_HIP
__device__ static int threadIdx_x_() { return hipThreadIdx_x; }

__device__ static int blockIdx_x_() { return hipBlockIdx_x; }
__device__ static int blockIdx_y_() { return hipBlockIdx_y; }
__device__ static int blockIdx_z_() { return hipBlockIdx_z; }

__device__ static int blockDim_x_() { return hipBlockDim_x; }

__device__ static int gridDim_y_() { return hipGridDim_y; }
#else
__device__ static int threadIdx_x_() { return 0; }

__device__ static int blockIdx_x_() { return 0; }
__device__ static int blockIdx_y_() { return 0; }
__device__ static int blockIdx_z_() { return 0; }

__device__ static int blockDim_x_() { return 0; }

__device__ static int gridDim_y_() { return 0; }
#endif

//-----------------------------------------------------------------------------
/// \brief Provide needed constants of GemmIn_t type.
///
///        Provide constants "0", "1" and "2" in the datatype
///        required to input to the reduced precision GEMM.

// Note: we could use __half here instead of uint16_t.  The intent here
// was to use a type based on standard C/C++.  No actual computations
// are done in this code based on the specifics of the type, so it doesn't
// matter.  Important thing is that sizeof(uint16_t) == sizeof(__half) == 2.

template<typename GemmIn_t> struct TCBufTypes;

template<> struct TCBufTypes<uint16_t> {
  static __host__ __device__ uint16_t zero() {return (uint16_t)0x0000;}
  static __host__ __device__ uint16_t one() {return (uint16_t)0x3c00;}
  static __host__ __device__ uint16_t two() {return (uint16_t)0x4000;}
                                        // = *(uint16_t*)&__float2half(0.);
                                        // = *(uint16_t*)&__float2half(1.);
                                        // = *(uint16_t*)&__float2half(2.);
};

//----------

template<> struct TCBufTypes<int8_t> {
  static __host__ __device__ int8_t zero() {return (int8_t)0;}
  static __host__ __device__ int8_t one() {return (int8_t)1;}
  static __host__ __device__ int8_t two() {return (int8_t)2;}
};

//----------

template<> struct TCBufTypes<GMFp32> {
  static __host__ __device__ GMFp32 zero() {return (GMFp32)0;}
  static __host__ __device__ GMFp32 one() {return (GMFp32)1;}
  static __host__ __device__ GMFp32 two() {return (GMFp32)2;}
};

//-----------------------------------------------------------------------------
/// \brief Select types etc. based on the setting of the tc param.

template<int TC_METHOD> struct TCSelector;

template<> struct TCSelector<TC::INT8> {
  // types.
  typedef int8_t GemmIn_t;
  typedef int32_t GemmOut_t;
#if defined USE_CUDA
  // type selector parameters.
  static cudaDataType __host__ __device__ gemm_type_in() {return CUDA_R_8I;}
  static cudaDataType __host__ __device__ gemm_type_out() {return CUDA_R_32I;}
#elif defined USE_HIP
  // type selector parameters.
  static rocblas_datatype __host__ __device__ gemm_type_in() {
   return rocblas_datatype_u8_r;
  }
  static rocblas_datatype __host__ __device__ gemm_type_out() {
   return rocblas_datatype_u32_r;
  }
#endif
  enum { COUNT = 1 };
};

//----------

template<> struct TCSelector<TC::FP16> {
  // types.
  typedef uint16_t GemmIn_t;
  typedef GMFp32 GemmOut_t;
#if defined USE_CUDA
  // type selector parameters.
  static cudaDataType __host__ __device__ gemm_type_in() {return CUDA_R_16F;}
  static cudaDataType __host__ __device__ gemm_type_out() {return CUDA_R_32F;}
#elif defined USE_HIP
  // type selector parameters.
  static rocblas_datatype __host__ __device__ gemm_type_in() {
    return rocblas_datatype_f16_r;
  }
  static rocblas_datatype __host__ __device__ gemm_type_out() {
    return rocblas_datatype_f32_r;
  }
#endif
  enum { COUNT = 1 };
};

//----------

template<> struct TCSelector<TC::FP32> {
  // types.
  typedef GMFp32 GemmIn_t;
  typedef GMFp32 GemmOut_t;
#if defined USE_CUDA
  // type selector parameters.
  static cudaDataType __host__ __device__ gemm_type_in() {return CUDA_R_32F;}
  static cudaDataType __host__ __device__ gemm_type_out() {return CUDA_R_32F;}
#elif defined USE_HIP
  // type selector parameters.
  static rocblas_datatype __host__ __device__ gemm_type_in() {
    return rocblas_datatype_f32_r;
  }
  static rocblas_datatype __host__ __device__ gemm_type_out() {
    return rocblas_datatype_f32_r;
  }
#endif
  enum { COUNT = 1 };
};

//=============================================================================
// FILE_LOCAL (STATIC) FUNCTIONS
//=============================================================================

//-----------------------------------------------------------------------------
/// \brief Write individual elements to buf.

//CHANGE: for non-kernel functions, struct to manage these variables

template<typename GemmIn_t>
__host__ __device__ static void gm_tc_buf_write_kernel_elt_(
  GemmIn_t* vo,
  const uint32_t* vi32,
  int vi32_dim0,
  int num_way,
  bool is_sparse,
  bool is_right,
  bool is_duo,
  int nvlea,
  int nvle,
  int nvleD2,
  int nvleX2,
  int nfl,
  int nflD2,
  int nflD2_thisstep,
  int flD2_min,
  int vlX2,
  int flD2_thisstep) {

  // Two fields (seminibbles) map to two halves of (2*sizeof(GemmIn_t))-bit word

  const int i01 = vlX2 % 2; // count either 0 bits or 1 bits.
  const int vl = vlX2 / 2;

  const int flD2 = flD2_min + flD2_thisstep;

  // Output array interpreted as having GemmIn_t scalars has nfl rows.

  const uint32_t* const vi32_col = vi32 + vl * (size_t)vi32_dim0;

  // Pick up two consecutive field values:
  // first field seminibble0, second field seminibble1
  // Set to zero if outside of active range.

  const int nibble = vl<nvlea ? (vi32_col[flD2/8] >> (4*(flD2%8))) & 15 : 0;

  const int seminibble0 = nibble & 3;
  const int seminibble1 = (nibble>>2) & 3;

  // Count number of 0 (or 1) bits in respective seminibble.
  // Determine whether to skip (1,0) null indicator value.

  const bool skip_10 = is_sparse || (num_way == 3 && ! is_right);

  // Possible counts, represented in target type.
  const GemmIn_t zero = TCBufTypes<GemmIn_t>::zero();
  const GemmIn_t one  = TCBufTypes<GemmIn_t>::one();
  const GemmIn_t two  = TCBufTypes<GemmIn_t>::two();

  const GemmIn_t out0 = is_duo ? (
                          seminibble0 == 2         ? zero :
                          (seminibble0 & 1) == i01 ? one :
                                                     zero
                        ) : (
                          seminibble0 == 3*i01     ? two :
                          seminibble0 == 3*(1-i01) ? zero :
                                         !skip_10  ? one :
                          seminibble0 == 1         ? one :
                                                     zero
                        );

  const GemmIn_t out1 = is_duo ? (
                          seminibble1 == 2         ? zero :
                          (seminibble1 & 1) == i01 ? one :
                                                     zero
                        ) : (
                          seminibble1 == 3*i01     ? two :
                          seminibble1 == 3*(1-i01) ? zero :
                                         !skip_10  ? one :
                          seminibble1 == 1         ? one :
                                                     zero
                        );
//if (flD2_thisstep == 0) printf("%f %f\n", (float)out0, (float)out1);

  // Always keep pair of cols together, corresponding to the two i01 values.
  // Right case: straight copy of cols to cols in sequence.
  // Left case: interleave to make later swizzling of metrics array work:
  // [ A A B B C C D D E E F F ] -> [ A A D D B B E E C C F F]

  const int vl_index = is_right ? vl : vl < nvleD2 ? 2*vl : 2*vl - nvle + 1;
  const int vlX2_index = i01 + 2*vl_index;

  const int flD2_index = flD2_thisstep;

  const int fl_index_0 = 0 + 2 * flD2_index;
  const int fl_index_1 = 1 + 2 * flD2_index;

  const int vlX2_dim = nvleX2;

  vo[vlX2_index + vlX2_dim * (size_t)fl_index_0] = out0;
  vo[vlX2_index + vlX2_dim * (size_t)fl_index_1] = out1;
}

//-----------------------------------------------------------------------------
/// \brief GPU kernel to support gm_tc_buf_write_.

template<typename GemmIn_t>
__global__ static void gm_tc_buf_write_kernel_(
  GemmIn_t* vo,
  const uint32_t* vi32,
  int vi32_dim0,
  int num_way,
  bool is_sparse,
  bool is_right,
  bool is_duo,
  int nvlea,
  int nvle,
  int nvleD2,
  int nvleX2,
  int nfl,
  int nflD2,
  int nflD2_thisstep,
  int flD2_min) {

  // Two fields (seminibbles) map to two halves of (2*sizeof(GemmIn_t))-bit word

  const int vlX2 = threadIdx_x_() + blockIdx_x_() * blockDim_x_();
  const int flD2_thisstep = blockIdx_y_() + gridDim_y_() * blockIdx_z_();

  if (vlX2 >= nvleX2 || flD2_thisstep >= nflD2_thisstep) {
    return;
  }

  gm_tc_buf_write_kernel_elt_<GemmIn_t>(vo, vi32, vi32_dim0,
    num_way, is_sparse, is_right, is_duo,
    nvlea, nvle, nvleD2, nvleX2, nfl, nflD2, nflD2_thisstep, flD2_min,
    vlX2, flD2_thisstep);
}

//-----------------------------------------------------------------------------
/// \brief Convert bitwise matrix to required format for GEMM.

template<int TC_METHOD>
static void gm_tc_buf_write_(
  bool is_right,
  int I_max,
  int I_max_dim,
  int nvl,
  int npvfl,
  int npvfl_thisstep,
  int pvfl_min,
  const void* vi,
  TCBufs& tc_bufs,
  bool is_duo,
  GMEnv* env) {

  COMET_INSIST(env && vi);
  COMET_INSIST(I_max_dim >= 0 && I_max_dim <= nvl);
  COMET_INSIST(I_max >= 0 && I_max <= I_max_dim);
  COMET_INSIST(nvl >= 0);
  COMET_INSIST(npvfl >= 0);
  COMET_INSIST(tc_bufs.tc_buf_left);
  COMET_INSIST(tc_bufs.tc_buf_right);
  COMET_INSIST(npvfl >= 0);
  COMET_INSIST(npvfl_thisstep >= 0 && npvfl_thisstep <= npvfl);
  COMET_INSIST(pvfl_min >= 0 && pvfl_min + npvfl_thisstep <= npvfl);

  // num_vector-related dimensions.

  const int nvle = is_right ? nvl : I_max_dim; // effective nvl dimension
  const int nvleD2 = nvle / 2;
  const int nvleX2 = nvle * 2;
  const int nvlea = is_right ? nvl : I_max; // num active nvle; others zeroed
  // NOTE: we are ignoring the issue from decomp_mgr that
  // num_vector_active_local may be strictly less than num_vector_local;
  // doesn't matter: just compute on fake values that will later be ignored.

  COMET_INSIST(nvle % 2 == 0 && nvl % 2 == 0 &&
           "tc method here requires num_vector_local multiple of 2.");

  // num_field-related dimensions.

  const int nfl = npvfl * 64;
  const int nflD2 = nfl / 2;
  const int nfl_thisstep = npvfl_thisstep * 64;
  const int nflD2_thisstep = nfl_thisstep / 2;
  const int fl_min = pvfl_min * 64;
  const int flD2_min = fl_min / 2;
  // Remember: end padding is set to zero; will correct zero counts later.

  // accelerator thread dims.

  const int threadblocksize = 256;
  const int blockdim_y = 32768;
  const int num_threadblocks_0 = utils::ceil(nvleX2, threadblocksize);
  const int num_threadblocks_1 = utils::min(nflD2_thisstep, blockdim_y);
  const int num_threadblocks_2 = utils::ceil(nflD2_thisstep, blockdim_y);

  // Arrays.

  typedef typename TCSelector<TC_METHOD>::GemmIn_t GemmIn_t;
  uint32_t* vi32 = (uint32_t*)vi;
  const int vi32_dim0 = npvfl * 4; // 4 = sizeof(doublecomplex) / sizeof(int32)
  GemmIn_t* const tc_buf = is_right ? (GemmIn_t*)tc_bufs.tc_buf_right :
                                      (GemmIn_t*)tc_bufs.tc_buf_left;
  COMET_INSIST(nvleX2 * (size_t)(2*nflD2_thisstep) *
           sizeof(typename TCSelector<TC_METHOD>::GemmIn_t)
           <= tc_bufs.tc_buf_size &&
           "Subscriptrange error on tc buf.");

  if (env->compute_method() == ComputeMethod::GPU) {

    // Kernel call.

#ifdef USE_ACCEL
#  ifdef USE_HIP
    hipLaunchKernelGGL(
#  endif
    gm_tc_buf_write_kernel_<GemmIn_t>
#  ifdef USE_CUDA
        <<<
#  else
        ,
#  endif
        dim3(num_threadblocks_0, num_threadblocks_1, num_threadblocks_2),
        dim3(threadblocksize, 1, 1),
        0,
        env->stream_compute()
#  ifdef USE_CUDA
        >>> (
#  else
        ,
#  endif
      tc_buf, vi32, vi32_dim0,
      env->num_way(), env->sparse(), is_right, is_duo,
      nvlea, nvle, nvleD2, nvleX2, nfl, nflD2, nflD2_thisstep, flD2_min);

    System::accel_last_call_succeeded();

#else // USE_ACCEL

      int dummy = 0;
      dummy += num_threadblocks_0 + num_threadblocks_1 + num_threadblocks_2;

#endif // USE_ACCEL

  } else { // if (env->compute_method() != ComputeMethod::GPU)

    for (int flD2_thisstep=0; flD2_thisstep<nflD2_thisstep; ++flD2_thisstep) {
      for (int vlX2=0; vlX2<nvleX2; ++vlX2) {

        gm_tc_buf_write_kernel_elt_<GemmIn_t>(
          tc_buf, vi32, vi32_dim0,
          env->num_way(), env->sparse(), is_right, is_duo,
          nvlea, nvle, nvleD2, nvleX2, nfl, nflD2, nflD2_thisstep, flD2_min,
          vlX2, flD2_thisstep);

      }
    }

  } // if compute_method
}

//-----------------------------------------------------------------------------
/// \brief Call cublas to perform required GEMM.

template<int TC_METHOD>
static void gm_tc_solve_accelblasgemmex_(
  bool is_first,
  int m,
  int n,
  int k,
  const void* matA,
  const void* matB,
  void* matC,
  TCBufs& tc_bufs,
  GMEnv* env) {

  COMET_INSIST(env && matA && matB && matC);

  // See https://devblogs.nvidia.com/programming-tensor-cores-cuda-9/
  // "Invoke the GEMM, ensuring k, lda, ldb, and ldc are all multiples of 8, 
  // and m is a multiple of 4"
  // "GEMMs that do not satisfy the above rules will fall back
  // to a non-Tensor Core implementation"
  // See also https://docs.nvidia.com/cuda/cublas/index.html#cublas-gemmEx

  // nfl is derived from padded-up npvfl, so always ok.
  COMET_INSIST(k % 8 == 0 && "Failed divisibility condition for tc gemm.");
  // since I_max_dim % 4 == 0; see gm_gemm_divisibility_required()
  COMET_INSIST(m % 8 == 0 && "Failed divisibility condition for tc gemm.");
  // since nvl % 4 == 0; see gm_gemm_divisibility_required()
  COMET_INSIST(n % 8 == 0 && "Failed divisibility condition for tc gemm.");

  // Make BLAS call.

  if (env->compute_method() == ComputeMethod::GPU) {

#ifdef USE_ACCEL

    const typename TCSelector<TC_METHOD>::GemmOut_t alpha = 1;
    const typename TCSelector<TC_METHOD>::GemmOut_t beta = is_first ? 0 : 1;

    // GPU BLAS call.

#  ifdef USE_CUDA
    const cublasStatus_t status = cublasGemmEx(
#  else
    //int status = rocblas_gemm_ex(
    const rocblas_status status = rocblas_gemm_ex(
#  endif
      tc_bufs.accelblas_handle
#  ifdef USE_CUDA
      , CUBLAS_OP_N, CUBLAS_OP_T
#  else
      , rocblas_operation_none, rocblas_operation_transpose
#  endif
      , m, n, k
      , (void*)&alpha
      , tc_bufs.tc_buf_left, TCSelector<TC_METHOD>::gemm_type_in(), m
      , tc_bufs.tc_buf_right, TCSelector<TC_METHOD>::gemm_type_in(), n
      , (void*)&beta
      , matC, TCSelector<TC_METHOD>::gemm_type_out(), m
#  ifdef USE_HIP
      , matC, TCSelector<TC_METHOD>::gemm_type_out(), m
#  endif
      , TCSelector<TC_METHOD>::gemm_type_out()
#  ifdef USE_CUDA
      //, CUBLAS_GEMM_ALGO3_TENSOR_OP // best timing, for cuda 9.1.85, transpose
      //, CUBLAS_GEMM_DFALT_TENSOR_OP // good timing, for cuda 9.2.88, transpose
      , CUBLAS_GEMM_ALGO4_TENSOR_OP // best timing, for cuda 9.2.88, transpose
#  else
      , rocblas_gemm_algo_standard
      , 0, 0  // solution_index, flags, workspace_size, workspace
#  endif
    );
    // TODO: use CUDA 10 autotuning capability here (later).

#  ifdef USE_CUDA
    if (status == CUBLAS_STATUS_NOT_INITIALIZED) {
      printf("Error: CUBLAS_STATUS_NOT_INITIALIZED\n");
    } else if (status == CUBLAS_STATUS_ARCH_MISMATCH) {
      printf("Error: CUBLAS_STATUS_ARCH_MISMATCH\n");
    } else if (status == CUBLAS_STATUS_NOT_SUPPORTED) {
      printf("Error: CUBLAS_STATUS_NOT_SUPPORTED\n");
    } else if (status == CUBLAS_STATUS_INVALID_VALUE) {
      printf("Error: CUBLAS_STATUS_INVALID_VALUE\n");
    } else if (status == CUBLAS_STATUS_EXECUTION_FAILED) {
      printf("Error: CUBLAS_STATUS_EXECUTION_FAILED\n");
    }
    COMET_INSIST(status == CUBLAS_STATUS_SUCCESS &&
             "Failure in call to cublasGemmEx.");
#  else
    COMET_INSIST(status == rocblas_status_success &&
             "Failure in call to rocblas_gemm_ex.");
#  endif

#else // USE_ACCEL

    COMET_INSIST(false && "Failure to call GEMM function.");

#endif // USE_ACCEL

    System::accel_last_call_succeeded();

  } else { // if (env->compute_method() != ComputeMethod::GPU) {

#ifdef USE_CPUBLAS

    COMET_INSIST(env->tc_eff() == TC::FP32);

    const float alpha = 1;
    const float beta = is_first ? 0 : 1;

    // CPU BLAS call.

    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans,
      m, n, k, alpha, (float*)tc_bufs.tc_buf_left, m,
      (float*)tc_bufs.tc_buf_right, n, beta, (float*)matC, m);

//for (int i=0; i< m*n; ++i)
//  printf("%i %f\n", i, ((float*)matC)[i]);

#else // USE_CPUBLAS

    COMET_INSIST(false && "Failure to call GEMM function.");

#endif // USE_CPUBLAS

  } // if compute_method

  env->ops_local_inc(2 * m * (double)n * (double)k);
}

//-----------------------------------------------------------------------------
/// \brief Call to perform required GEMM.

template<int TC_METHOD>
static void gm_tc_solve_(
  bool is_first,
  int nvll,
  int nvl,
  int npvfl_thisstep,
  const void* matA,
  const void* matB,
  void* matC,
  TCBufs& tc_bufs,
  GMEnv* env) {

  COMET_INSIST(env && matA && matB && matC);
  COMET_INSIST(nvll >= 0);
  COMET_INSIST(nvl >= 0);
  COMET_INSIST(nvll <= nvl);
  COMET_INSIST(npvfl_thisstep >= 0);
  COMET_INSIST(env->tc_eff() != TC::NO);

  const int nfl_thisstep = npvfl_thisstep * 64;

  const int m = 2 * nvll; // metrics array dim
  const int n = 2 * nvl; // metrics array dim
  const int k = nfl_thisstep; // vectors array (as GemmIn_t) dim

  gm_tc_solve_accelblasgemmex_<TC_METHOD>(is_first, m, n, k,
                                          matA, matB, matC, tc_bufs, env);
}

//-----------------------------------------------------------------------------
/// \brief Swizzle individual elements in buf.

template<typename GemmOut_t>
__host__ __device__ static void gm_tc_repair_metrics_kernel_elt_(
  int nvl, int nvll, int nvllD2, void* vo,
  int thread_r, int thread_c) { 

  // Considered as an array of floats, array is 2*nvl rows X 2*nvl cols.
  // Each thread manipulates a block of 4 rows and 2 cols.
  // Thus the dimensions of the metrics array in blocks is nvllD2 X nvl.
  // Each block viewed as an array of doubles is 2 X 2.

  // Two col numbers being processed of this (float) array.

  // ISSUE: does the compiler need to / understand that the pointers are aliased

  const size_t fcr_offset0 = 4*thread_r + thread_c * (size_t)(4*nvll);
  const size_t fcr_offset1 = 4*thread_r + thread_c * (size_t)(4*nvll) + 2*nvll;

  // Read the 8 values.

  GemmOut_t* const fvo = (GemmOut_t*)vo;

  const GemmOut_t f00 = fvo[fcr_offset0+0];
  const GemmOut_t f01 = fvo[fcr_offset0+1];
  const GemmOut_t f02 = fvo[fcr_offset0+2];
  const GemmOut_t f03 = fvo[fcr_offset0+3];

  const GemmOut_t f10 = fvo[fcr_offset1+0];
  const GemmOut_t f11 = fvo[fcr_offset1+1];
  const GemmOut_t f12 = fvo[fcr_offset1+2];
  const GemmOut_t f13 = fvo[fcr_offset1+3];
//printf("%f %f %f %f %f %f %f %f\n", (float)f00, (float)f01, (float)f02, (float)f03, (float)f10, (float)f11, (float)f12, (float)f13);

  // Apply the permutation:

  // [ f00  f10 ]  ->  [ f00  f02 ]
  // [ f01  f11 ]  ->  [ f01  f03 ]
  // [ f02  f12 ]  ->  [ f10  f12 ]
  // [ f03  f13 ]  ->  [ f11  f13 ]

  const GemmOut_t f00p = f00;
  const GemmOut_t f01p = f01;

  const GemmOut_t f02p = f10;
  const GemmOut_t f03p = f11;

  const GemmOut_t f10p = f02;
  const GemmOut_t f11p = f03;

  const GemmOut_t f12p = f12;
  const GemmOut_t f13p = f13;

  // Use "shifter" to move a value to the upper half of the mantissa.

  const double shifter = (((uint32_t)1) << GM_TALLY1_MAX_VALUE_BITS);

  // Pack two 26-bit integers into mantissa of double.

  const double d00 = (double)f00p + (double)f02p * shifter;
  const double d01 = (double)f01p + (double)f03p * shifter;

  const double d10 = (double)f10p + (double)f12p * shifter;
  const double d11 = (double)f11p + (double)f13p * shifter;

  // Overwrite block with the new values.
  // All is isolated to a single thread, should be thread safe.

  const size_t dc_offset0 = thread_c * (size_t)(2*nvll);
  const size_t dc_offset1 = thread_c * (size_t)(2*nvll) + nvll;

  const size_t dcr_offset0 = dc_offset0 + 2*thread_r;
  const size_t dcr_offset1 = dc_offset1 + 2*thread_r;

  double* const dvo = (double*)vo;

  dvo[dcr_offset0+0] = d00;
  dvo[dcr_offset0+1] = d01;

  dvo[dcr_offset1+0] = d10;
  dvo[dcr_offset1+1] = d11;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
/// \brief GPU kernel to support gm_tc_repair_metrics_.
///
///        This function has two purposes:
///        1. Convert the 2X2 table from each pair of compared vectors
///        from 4 32-bit (int32 or float32) values to the required
///        16-byte double complex packed format.
///        2. Permute the table elements to the required places.
///
///        The reason for the permutation is as follows.
///        For the output matrix of this function, each single 2X2 matrix
///        is arranged contiguously in memory as a double complex value.
///        However, the input matrices to the GEMM do not give a result
///        matrix that is consistent with this ordering.
///        Thus there needs to be a copy to rearrange.  Furthermore,
///        we want to make this an in-place rearrangement to save
///        space, and additionally we want to assign work to threads
///        with no race conditions and with coalesced memory accesses.
///
///        The method can be explained as follows.
///        1. The input "left" and "right" matrices to the modified GEMM
///        can be thought of each as a group of column vectors.
///        2. Each column (of 2-bit entries) is converted into two columns,
///        with entries being the counts of 0 bits and 1 bits of the
///        original vectors.  Each pair of vectors is kept together
///        side-by-side in these new left and right matrices L and R.
///        3. The columns of L are permuted, to give L' = L P
///        Example:
///          R  = [ G, G, H, H, I, I, J, J, K, K, L, L ]
///          L  = [ A, A, B, B, C, C, D, D, E, E, F, F ]
///          L' = [ A, A, D, D, B, B, E, E, C, C, F, F ]
///        (note L is used in 2 different senses here)
///        4. The GEMM is computed, M = (L')^T R = P^T L^T R.  Because of
///        the permutation of L, the rows of M are permuted.
///        Here, for brevity we drop the transpose, writing A^T G as AG, etc.
///          M = [ AG, AG, AH, AH, . . . ]
///              [ AG, AG, AH, AH, . . . ]
///              [ DG, DG, DH, DH, . . . ]
///              [ DG, DG, DH, DH, . . . ]
///              [ BG, BG, BH, BH, . . . ]
///              [ BG, BG, BH, BH, . . . ]
///              [ EG, EG, EH, EH, . . . ]
///              [ EG, EG, EH, EH, . . . ]
///              [ CG, CG, CH, CH, . . . ]
///              [ CG, CG, CH, CH, . . . ]
///              [ FG, FG, FH, FH, . . . ]
///              [ FG, FG, FH, FH, . . . ]
///        Here we are considering M to be stored in column-major order.
///        5. Next we consider this as composed of size 4X2 blocks,
///        assign a CUDA thread to each block and do an in-block
///        permutation. Note each thread loads 2 16-byte (double) words,
///        with stride between threads of 16 bytes.
///        (need to check on efficiency of this w.r.t. coalescing etc.)
///          [ AG, AG ] -> [ AG, DG ]
///          [ AG, AG ] -> [ AG, DG ]
///          [ DG, DG ] -> [ AG, DG ]
///          [ DG, DG ] -> [ AG, DG ]
///        As can be seen, all four entries AG of the table are now
///        contiguous in memory.

template<typename GemmOut_t>
__global__ static void gm_tc_repair_metrics_kernel_(
  int nvl, int nvll, int nvllD2, void* vo) { 

  // Row and column threads of metrics array.
  const int thread_r = threadIdx_x_() + blockIdx_x_() * blockDim_x_();
  const int thread_c = blockIdx_y_();

  if (thread_r >= nvllD2 || thread_c >= nvl) {
    return;
  }

  gm_tc_repair_metrics_kernel_elt_<GemmOut_t>(
    nvl, nvll, nvllD2, vo,
    thread_r, thread_c);
}

//-----------------------------------------------------------------------------
/// \brief Swizzle/cast values from cublas call into double complex format.
///
///        The cublas gemm poduces a matrix of scalars of 32 bit size
///        (int32 or float).  However the required format of the metrics
///        is a matrix of double complex values, with each double
///        containing two packed 26-bit integers.
///        This code does an in-place transformation from one to the other.

template<int TC_METHOD>
static void gm_tc_repair_metrics_(
  int nvll,
  int nvl,
  void* vo,
  TCBufs& tc_bufs,
  GMEnv* env) {

  COMET_INSIST(env && vo);
  COMET_INSIST(nvll >= 0);
  COMET_INSIST(nvl >= 0);
  COMET_INSIST(nvll <= nvl);

  // always true, because of gm_gemm_divisibility_required()
  COMET_INSIST(nvll % 2 == 0 && "Failed divisibility condition for tc gemm.");
  const int nvllD2 = nvll / 2;

  const int threadblocksize = 256;
  const int vll2_threadblocks = utils::ceil(nvllD2, threadblocksize);

  typedef typename TCSelector<TC_METHOD>::GemmOut_t GemmOut_t;

  if (env->compute_method() == ComputeMethod::GPU) {

    // Kernel call.

#ifdef USE_ACCEL

#  ifdef USE_HIP
    hipLaunchKernelGGL(
#  endif
    gm_tc_repair_metrics_kernel_<GemmOut_t>
#  ifdef USE_CUDA
        <<<
#  else
        ,
#  endif
        dim3(vll2_threadblocks, nvl, 1),
        dim3(threadblocksize, 1, 1),
        0,
        env->stream_compute()
#  ifdef USE_CUDA
        >>> (
#  else
        ,
#  endif
        nvl, nvll, nvllD2, vo);

    System::accel_last_call_succeeded();

#else // USE_ACCEL

  int dummy = 0;
  dummy += vll2_threadblocks + threadblocksize;

#endif // USE_ACCEL

  } else { // if (env->compute_method() != ComputeMethod::GPU)

//    for (int thread_c=0; thread_c<nvllD2; ++thread_c) {
//      for (int thread_r=0; thread_r<nvl; ++thread_r) {
    for (int thread_c=0; thread_c<nvl; ++thread_c) {
      for (int thread_r=0; thread_r<nvllD2; ++thread_r) {

        gm_tc_repair_metrics_kernel_elt_<GemmOut_t>(
          nvl, nvll, nvllD2, vo, thread_r, thread_c);

      }
    }

  } // if compute_method
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
static void gm_tc_gemm_start_impl_(
  int m, int n, int k,
  const void* matA, int ldda,
  const void* matB, int lddb,
  void* matC, int lddc,
  TCBufs& tc_bufs,
  GMEnv* env) {

  COMET_INSIST(ldda == k && lddb == k); // For our purposes, always true

  const int nvl = n;
  const int npvfl = k;
  const int I_max = m;
  const int I_max_dim = lddc;
  COMET_INSIST(I_max <= I_max_dim);
  COMET_INSIST(I_max_dim <= nvl);
  // nvll is the effective nvl (column dim) for left matrix
  // only really need to compute up to I_max, but need to compute to I_max_dim
  // to satisfy divisibility requirements.
  // Note nvl is always the column dim for the right matrix.
  const int nvll = I_max_dim;
  COMET_INSIST((size_t)nvll == gm_gemm_size_required(nvll, env));

  const int num_steps = env->num_tc_steps();

  // Loop over steps of algorithm.
  for (int step_num = 0; step_num < num_steps; ++step_num) {

    // Select the block row of the left and right matrices for this step.
    const int pvfl_min = ((step_num+0) * npvfl) / num_steps;
    const int pvfl_max = ((step_num+1) * npvfl) / num_steps;
    const int npvfl_thisstep = pvfl_max - pvfl_min;

    if (npvfl_thisstep == 0) {  // empty block row
      continue;
    }

    // Convert the input matrices of packed bit values into matrices
    // of values of a type suitable for the GEMM.
    const bool left_matrix = false; // A
    const bool right_matrix = true; // B
    const bool is_duo = env->metric_type() == MetricType::DUO;
    gm_tc_buf_write_<TC_METHOD>(left_matrix, I_max, I_max_dim, nvl, npvfl,
                     npvfl_thisstep, pvfl_min, matA, tc_bufs, is_duo, env);
    gm_tc_buf_write_<TC_METHOD>(right_matrix, I_max, I_max_dim, nvl, npvfl,
                     npvfl_thisstep, pvfl_min, matB, tc_bufs, is_duo, env);

    // Perform the GEMM for this pair of block rows; accumulate.
    gm_tc_solve_<TC_METHOD>(
      pvfl_min==0, nvll, nvl, npvfl_thisstep, matA, matB, matC, tc_bufs, env);
  }

  // Revise the results of the GEMMs to be in the needed double complex format.
  gm_tc_repair_metrics_<TC_METHOD>(nvll, nvl, matC, tc_bufs, env);
}

//=============================================================================
// "PUBLIC" FUNCTIONS
//=============================================================================

//-----------------------------------------------------------------------------
/// \brief Divisibility requirement for GEMM.

size_t gm_gemm_divisibility_required(GMEnv* const env) {
  COMET_INSIST(env);

  const bool need_divisible_by_4 = env->tc_eff() != TC::NO;

  return need_divisible_by_4 ? 4 : 1;
}

//-----------------------------------------------------------------------------
/// \brief Size requirement for GEMM.

size_t gm_gemm_size_required(size_t size_requested, GMEnv* const env) {
  COMET_INSIST(env);

  const size_t factor = gm_gemm_divisibility_required(env);

  return utils::ceil(size_requested, factor)*factor;
}

//-----------------------------------------------------------------------------
/// \brief Use a standard GEMM to compute CoMet metrics bitwise result.

//CHANGE: pass in vectors_I, vectors_J_col pointers (allow both this and older option?)

//CHANGE: ? eliminate ldda, lddb, since redundant with other vars.

void gm_tc_gemm_start(int m, int n, int k,
                      const void* matA, int ldda,
                      const void* matB, int lddb,
                      void* matC, int lddc,
                      TCBufs& tc_bufs,
                      GMEnv* env) {
  COMET_INSIST(matA && matB && matC && env);
  COMET_INSIST(m >= 0 && n >= 0 && k >= 0);
  COMET_INSIST(ldda >= 0 && lddb >= 0 && lddc >= 0);
  COMET_INSIST(k <= ldda);
  COMET_INSIST(k <= lddb);
  COMET_INSIST(m <= lddc);
  COMET_INSIST(env->tc_eff() != TC::NO);
  COMET_INSIST(env->is_metric_type_bitwise());

  COMET_INSIST(tc_bufs.tc_buf_left);

  // Select required template function instance.

  switch (env->tc_eff()) {
    // --------------
    case TC::INT8: {
      gm_tc_gemm_start_impl_<TC::INT8>(
        m, n, k, matA, ldda, matB, lddb, matC, lddc, tc_bufs,  env);
    } break;
    // --------------
    case TC::FP16: {
      gm_tc_gemm_start_impl_<TC::FP16>(
        m, n, k, matA, ldda, matB, lddb, matC, lddc, tc_bufs,  env);
    } break;
    // --------------
    case TC::FP32: {
      gm_tc_gemm_start_impl_<TC::FP32>(
        m, n, k, matA, ldda, matB, lddb, matC, lddc, tc_bufs,  env);
    } break;
    // --------------
    default:
      COMET_INSIST(false && "Invalid tc type.");
  } // switch
}

//-----------------------------------------------------------------------------
/// \brief Initialize TCBufs object by allocating memory etc.

void gm_tc_bufs_malloc(int num_vector_local,
                       int num_field_local,
                       int num_packedval_field_local,
                       TCBufs& tc_bufs,
                       GMEnv* env) {
  COMET_INSIST(env);
  COMET_INSIST(num_vector_local >= 0);
  COMET_INSIST(num_packedval_field_local >= 0);
  COMET_INSIST(!tc_bufs.tc_buf_left);
  COMET_INSIST(!tc_bufs.tc_buf_right);

  if (!env->is_metric_type_bitwise() || env->tc_eff() == TC::NO)
    return;

  // Calculate sizes.

  const size_t nvl = num_vector_local;
  const size_t npvfl = num_packedval_field_local;
  const size_t npvfl_thisstep_max = utils::ceil(npvfl, (size_t)env->num_tc_steps());

  const int sizeof_gemm_in_t =
     env->tc_eff() == TC::INT8 ?
       sizeof(typename TCSelector<TC::INT8>::GemmIn_t) :
     env->tc_eff() == TC::FP16 ?
       sizeof(typename TCSelector<TC::FP16>::GemmIn_t) :
     env->tc_eff() == TC::FP32 ?
       sizeof(typename TCSelector<TC::FP32>::GemmIn_t) :
     0;
  COMET_INSIST(TC::is_valid(env->tc_eff())); // this code must be updated if new method

  const size_t nvlX2 = nvl * 2;

  tc_bufs.tc_buf_size = nvlX2 * (npvfl_thisstep_max * 64) * sizeof_gemm_in_t;
  tc_bufs.tc_buf_size = tc_bufs.tc_buf_size ? tc_bufs.tc_buf_size : 1;

  if (env->compute_method() == ComputeMethod::GPU) {

    // Allocate buffers.

#if defined USE_CUDA
    cudaMalloc(&tc_bufs.tc_buf_left, tc_bufs.tc_buf_size);
#elif defined USE_HIP
    hipMalloc(&tc_bufs.tc_buf_left, tc_bufs.tc_buf_size);
#endif
    System::accel_last_call_succeeded();
    env->gpu_mem_local_inc(tc_bufs.tc_buf_size);

#if defined USE_CUDA
    cudaMalloc(&tc_bufs.tc_buf_right, tc_bufs.tc_buf_size);
#elif defined USE_HIP
    hipMalloc(&tc_bufs.tc_buf_right, tc_bufs.tc_buf_size);
#endif
    System::accel_last_call_succeeded();
    env->gpu_mem_local_inc(tc_bufs.tc_buf_size);

    // Set up accel blas handle.

#if defined USE_CUDA
    cublasStatus_t status = cublasCreate(&tc_bufs.accelblas_handle);
    COMET_INSIST(status == CUBLAS_STATUS_SUCCESS && "Error in cublasCreate.");

    status = cublasSetStream(tc_bufs.accelblas_handle, env->stream_compute());
    COMET_INSIST(status == CUBLAS_STATUS_SUCCESS && "Error in cublasSetStream.");

    status = cublasSetMathMode(tc_bufs.accelblas_handle, CUBLAS_TENSOR_OP_MATH);
    COMET_INSIST(status == CUBLAS_STATUS_SUCCESS && "Error in cublasSetMathMode.");
#elif defined USE_HIP
    //rocbas_status_ status = rocblas_create_handle(&tc_bufs.accelblas_handle);
    int status = rocblas_create_handle(&tc_bufs.accelblas_handle);
    COMET_INSIST(status == rocblas_status_success &&
             "Error in rocblas_create_handle.");

    status = rocblas_set_stream(tc_bufs.accelblas_handle, env->stream_compute());
    COMET_INSIST(status == rocblas_status_success &&
             "Error in rocblas_set_stream.");

    //FIX - will this be needed?
    //  status = cublasSetMathMode(tc_bufs.accelblas_handle,
    //                             CUBLAS_TENSOR_OP_MATH);
    //  COMET_INSIST(status == CUBLAS_STATUS_SUCCESS &&
    //           "Error in cublasSetMathMode.");
#endif

  } else { // compute_method

    // Allocate buffers.

    tc_bufs.tc_buf_left = malloc(tc_bufs.tc_buf_size);
    COMET_INSIST(tc_bufs.tc_buf_left);
    env->cpu_mem_local_inc(tc_bufs.tc_buf_size);

    tc_bufs.tc_buf_right = malloc(tc_bufs.tc_buf_size);
    COMET_INSIST(tc_bufs.tc_buf_right);
    env->cpu_mem_local_inc(tc_bufs.tc_buf_size);
//memset((void*)tc_bufs.tc_buf_left, 0, tc_bufs.tc_buf_size);
//memset((void*)tc_bufs.tc_buf_right, 0, tc_bufs.tc_buf_size);

  } // compute_method
}

//-----------------------------------------------------------------------------
/// \brief Terminate TCBufs object by deallocating memory etc.

void gm_tc_bufs_free(TCBufs& tc_bufs, GMEnv* env) {
  COMET_INSIST(env);
  COMET_INSIST((tc_bufs.tc_buf_left != 0) == (tc_bufs.tc_buf_right != 0));

  if (!tc_bufs.tc_buf_left) {
    return;
  }

  if (env->compute_method() == ComputeMethod::GPU) {

    // Free buffers.

#if defined USE_CUDA
    cudaFree(tc_bufs.tc_buf_left);
#elif defined USE_HIP
    hipFree(tc_bufs.tc_buf_left);
#endif
    System::accel_last_call_succeeded();
    tc_bufs.tc_buf_left = NULL;
    env->gpu_mem_local_dec(tc_bufs.tc_buf_size);

#if defined USE_CUDA
    cudaFree(tc_bufs.tc_buf_right);
#elif defined USE_HIP
    hipFree(tc_bufs.tc_buf_right);
#endif
    System::accel_last_call_succeeded();
    tc_bufs.tc_buf_right = NULL;
    env->gpu_mem_local_dec(tc_bufs.tc_buf_size);

    // Free accel blas handle.

#if defined USE_CUDA
    cublasStatus_t status = cublasDestroy(tc_bufs.accelblas_handle);
    COMET_INSIST(status == CUBLAS_STATUS_SUCCESS && "Error in cublasDestroy.");
#elif defined USE_HIP
    //rocblas_status status = rocblas_destroy_handle(tc_bufs.accelblas_handle);
    int status = rocblas_destroy_handle(tc_bufs.accelblas_handle);
    COMET_INSIST(status == rocblas_status_success &&
             "Error in rocblas_destroy_handle.");
#endif

  } else { // compute_method

    // Free buffers.

    free(tc_bufs.tc_buf_left);
    tc_bufs.tc_buf_left = NULL;
    env->cpu_mem_local_dec(tc_bufs.tc_buf_size);

    free(tc_bufs.tc_buf_right);
    tc_bufs.tc_buf_right = NULL;
    env->cpu_mem_local_dec(tc_bufs.tc_buf_size);

  } // compute_method
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
