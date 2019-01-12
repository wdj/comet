//-----------------------------------------------------------------------------
/*!
 * \file   linalg_cuda.cu
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  CUDA code to support linear algebra operations.
 * \note   Copyright (C) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "stdint.h"

#include "linalg_cuda.cuh"

//-----------------------------------------------------------------------------
#ifdef USE_TC
//-----------------------------------------------------------------------------

#include "cuda_fp16.h"

//-----------------------------------------------------------------------------
/// \brief Specialized class to support gm_tc_buf_write_kernel_ type seletion.

// Note: we could use __half here instead of GMUInt16.  The intent here
// was to use a type based on standard C/C++.  No actual computations
// are done in this code based on the specifics of the type, so it doesn't
// matter.  Important thing is that sizeof(GMUInt16) == sizeof(__half) == 2.

template<typename Buf_t> struct TCBuf_types;

template<> struct TCBuf_types<GMUInt16> {
  static __device__ GMUInt16 zero() {return (GMUInt16)0x0000;}
                                        // = *(GMUInt16*)&__float2half(0.);
  static __device__ GMUInt16 one() {return (GMUInt16)0x3c00;}
                                        // = *(GMUInt16*)&__float2half(1.);
  static __device__ GMUInt16 two() {return (GMUInt16)0x4000;}
                                        // = *(GMUInt16*)&__float2half(2.);
};

template<> struct TCBuf_types<GMUInt8> {
  static __device__ GMUInt8 zero() {return (GMUInt8)0;}
  static __device__ GMUInt8 one() {return (GMUInt8)1;}
  static __device__ GMUInt8 two() {return (GMUInt8)2;}
};

//-----------------------------------------------------------------------------
/// \brief GPU kernel to support gm_tc_buf_write_.

template<typename Buf_t>
__global__ void gm_tc_buf_write_kernel_(
  int num_way,
  bool is_sparse,
  bool is_right,
  GMUInt32* vi32,
  int vi32_dim0,
  int nvla,
  int nvl,
  int nvl2,
  int nvlX2,
  int nfl,
  int nfl2,
  int nfl2_step,
  int fl2_min,
  void* vo) {

  // Two fields (seminibbles) map to two halves of (2*sizeof(Buf_t))-bit word

  const int vlX2 = threadIdx.x + blockIdx.x * blockDim.x;
  const int fl2_step = blockIdx.y + gridDim.y * blockIdx.z;

  if (vlX2 >= nvlX2 || fl2_step >= nfl2_step) {
    return;
  }

  const int i01 = vlX2 % 2; // count either 0 bits or 1 bits.
  const int vl = vlX2 / 2;

  const int fl2 = fl2_min + fl2_step;

  // Output array interpreted as having Buf_t scalars has nfl rows.

  const GMUInt32* const vi32_col = vi32 + vl * (size_t)vi32_dim0;

  // Pick up two consecutive field values:
  // first field seminibble0, second field seminibble1
  // Set to zero if outside of active range.

  const int nibble = vl<nvla ? (vi32_col[fl2/8] >> (4*(fl2%8))) & 15 : 0;

  const int seminibble0 = nibble & 3;
  const int seminibble1 = (nibble>>2) & 3;

  // Count number of 0 (or 1) bits in respective seminibble.
  // Determine whether to skip (1,0) null indicator value.

  const bool skip_10 = is_sparse || (num_way == 3 && ! is_right);

  // Possible counts, represented in target type.
  const Buf_t zero = TCBuf_types<Buf_t>::zero();
  const Buf_t one = TCBuf_types<Buf_t>::one();
  const Buf_t two = TCBuf_types<Buf_t>::two();

  const Buf_t out0 = seminibble0 == 3*i01     ? two :
                     seminibble0 == 3*(1-i01) ? zero :
                                    !skip_10  ? one :
                     seminibble0 == 1         ? one :
                                                zero;

  const Buf_t out1 = seminibble1 == 3*i01     ? two :
                     seminibble1 == 3*(1-i01) ? zero :
                                    !skip_10  ? one :
                     seminibble1 == 1         ? one :
                                                zero;
  // Always keep pair of cols together, corresponding to the two i01 values.
  // Right case: straight copy of cols to cols in sequence.
  // Left case: interleave to make later swizzling of metrics array work:
  // [ A A B B C C D D E E F F ] -> [ A A D D B B E E C C F F]

  const int vl_index = is_right ? vl : vl < nvl2 ? 2*vl : 2*vl - nvl + 1;
  const int vlX2_index = i01 + 2*vl_index;

  const int fl2_index = fl2_step;

  const int fl_index_0 = 0 + 2 * fl2_index;
  const int fl_index_1 = 1 + 2 * fl2_index;

  const int vlX2_dim = nvlX2;

  Buf_t* vo_typed = (Buf_t*)vo;

  vo_typed[vlX2_index + vlX2_dim * (size_t)fl_index_0] = out0;
  vo_typed[vlX2_index + vlX2_dim * (size_t)fl_index_1] = out1;
}

//-----------------------------------------------------------------------------
/// \brief Convert bitwise matrix to required format for GEMM.

void gm_tc_buf_write_(
  bool is_right,
  int I_max,
  int I_max_dim,
  int nvl,
  int npvfl,
  int npvfl_step,
  int pvfl_min,
  void* vi,
  TCBufs& tc_bufs,
  GMEnv* env) {

  GMInsist(env && vi);
  GMInsist(I_max_dim >= 0 && I_max_dim <= nvl);
  GMInsist(I_max >= 0 && I_max <= I_max_dim);
  GMInsist(nvl >= 0);
  GMInsist(npvfl >= 0);
  GMInsist(tc_bufs.tc_buf_left);
  GMInsist(tc_bufs.tc_buf_right);
//TODO: more assertions

  const bool is_int8 = env->tc == 2;

  const int nvle = is_right ? nvl : I_max_dim; // effective nvl dimension
  const int nvle2 = nvle / 2;
  const int nvleX2 = 2 * nvle;
  const int nvlea = is_right ? nvl : I_max; // num active nvle; others zeroed

  const int nfl = npvfl * 64;
  const int nfl2 = nfl / 2;
  const int nfl_step = npvfl_step * 64;
  const int nfl2_step = nfl_step / 2;
  const int fl_min = pvfl_min * 64;
  const int fl2_min = fl_min / 2;

  const int vi_dim0 = npvfl * 4; // 4 = sizeof(doublecomplex) / sizeof(int32)

  GMInsistInterface(env, nvle % 2 == 0 && nvl % 2 == 0 &&
                    "tc method here requires num_vector_local multiple of 2.");

  const int threadblocksize = 256;
  const int blockdim_y = 32768;
  const int num_threadblocks_0 = gm_ceil_i8(nvleX2, threadblocksize);
  const int num_threadblocks_1 = gm_min_i8(nfl2_step, blockdim_y);
  const int num_threadblocks_2 = gm_ceil_i8(nfl2_step, blockdim_y);

  // TODO: check dims against tc_bufs.tc_buf_size

  void* const tc_buf = is_right ? tc_bufs.tc_buf_right : tc_bufs.tc_buf_left;
  GMUInt32* vi32 = (GMUInt32*)vi;

  if (! is_int8) {

    gm_tc_buf_write_kernel_<GMUInt16><<<
        dim3(num_threadblocks_0, num_threadblocks_1, num_threadblocks_2),
        dim3(threadblocksize, 1, 1),
        0,
        env->stream_compute_>>>(
      GMEnv_num_way(env), env->sparse, is_right,
      vi32, vi_dim0, nvlea, nvle, nvle2, nvleX2,
      nfl, nfl2, nfl2_step, fl2_min, tc_buf);

  } else {

    gm_tc_buf_write_kernel_<GMUInt8><<<
        dim3(num_threadblocks_0, num_threadblocks_1, num_threadblocks_2),
        dim3(threadblocksize, 1, 1),
        0,
        env->stream_compute_>>>(
      GMEnv_num_way(env), env->sparse, is_right,
      vi32, vi_dim0, nvlea, nvle, nvle2, nvleX2,
      nfl, nfl2, nfl2_step, fl2_min, tc_buf);

  }

  GMEnv_cuda_last_call_succeeded(env);
}

//-----------------------------------------------------------------------------
/// \brief Call cublas to perform required GEMM.

void gm_tc_solve_(
  bool is_first,
  int nvll,
  int nvl,
  int npvfl_step,
  void* dA,
  void* dB,
  void* dC,
  TCBufs& tc_bufs,
  GMEnv* env) {

  GMInsist(env && dA && dB && dC);
  GMInsist(nvll >= 0);
  GMInsist(nvl >= 0);
  GMInsist(nvll <= nvl);
  GMInsist(npvfl_step >= 0);
  GMInsist(env->tc == 1 || env->tc == 2);

  const int nfl_step = npvfl_step * 64;

  const int m = 2 * nvll; // metrics array dim
  const int n = 2 * nvl; // metrics array dim
  const int k = nfl_step; // vectors array (as halfs/bytes) dim

  const float alpha_f32 = 1;
  const float beta_f32 = is_first ? 0 : 1;

  const int alpha_i32 = 1;
  const int beta_i32 = is_first ? 0 : 1;

  const bool is_int8 = env->tc == 2;

  // See https://devblogs.nvidia.com/programming-tensor-cores-cuda-9/
  // "Invoke the GEMM, ensuring k, lda, ldb, and ldc are all multiples of 8, 
  // and m is a multiple of 4"
  // "GEMMs that do not satisfy the above rules will fall back
  // to a non-Tensor Core implementation"
  // See also https://docs.nvidia.com/cuda/cublas/index.html#cublas-gemmEx

  GMInsist(k % 8 == 0); // nfl is derived from padded-up npvfl, so always ok.
  GMInsist(m % 8 == 0); // need I_max_dim % 4 == 0
  GMInsist(n % 8 == 0); // need nvl % 4 == 0

  // Make BLAS call.

  cublasStatus_t status = cublasGemmEx(
    tc_bufs.cublas_handle,
    CUBLAS_OP_N, CUBLAS_OP_T,
    m, n, k,
    is_int8 ? (void*)&alpha_i32 : (void*)&alpha_f32,
    tc_bufs.tc_buf_left,
    is_int8 ? CUDA_R_8I : CUDA_R_16F,
    m,
    tc_bufs.tc_buf_right,
    is_int8 ? CUDA_R_8I : CUDA_R_16F,
    n,
    is_int8 ? (void*)&beta_i32 : (void*)&beta_f32,
    //(void*)&beta_f32,
    dC,
    is_int8 ? CUDA_R_32I : CUDA_R_32F,
    m,
    is_int8 ? CUDA_R_32I : CUDA_R_32F,
    //CUBLAS_GEMM_ALGO3_TENSOR_OP // best timing, for cuda 9.1.85, transpose
    //CUBLAS_GEMM_DFALT_TENSOR_OP // good timing, for cuda 9.2.88, transpose
    CUBLAS_GEMM_ALGO4_TENSOR_OP // best timing, for cuda 9.2.88, transpose
    //CUBLAS_GEMM_DFALT_TENSOR_OP
  );

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

  GMInsist(status == CUBLAS_STATUS_SUCCESS);

  env->ops_local += 2 * m * (double)n * (double)k;
}

//-----------------------------------------------------------------------------
/// \brief GPU kernel to support gm_tc_repair_metrics_.

template<typename GemmOut_t>
__global__ void gm_tc_repair_metrics_kernel_(
  int nvl, int nvll, int nvll2, void* vo) { 
  // Row and column of metrics array.

  const int thread_r = threadIdx.x + blockIdx.x * blockDim.x;
  const int thread_c = blockIdx.y;

  if (thread_r >= nvll2 || thread_c >= nvl) {
    return;
  }

  // Considered as an array of floats, array is 2*nvl rows X 2*nvl cols.
  // Each thread manipulates a block of 4 rows and 2 cols.
  // Thus the dimensions of the metrics array in blocks is nvll2 X nvl.
  // Each block viewed as an array of doubles is 2 X 2.

  // Two col numbers being processed of this (float) array.

  // ISSUE: does the compiler need to / understand that the pointers are aliased

  const size_t fc_offset0 = thread_c * (size_t)(4*nvll);
  const size_t fc_offset1 = thread_c * (size_t)(4*nvll) + 2*nvll;

  const size_t fcr_offset0 = fc_offset0 + 4*thread_r;
  const size_t fcr_offset1 = fc_offset1 + 4*thread_r;

  // Read the 8 floats.

  GemmOut_t* fvo = (GemmOut_t*)vo;

  const GemmOut_t f00 = fvo[fcr_offset0+0];
  const GemmOut_t f01 = fvo[fcr_offset0+1];
  const GemmOut_t f02 = fvo[fcr_offset0+2];
  const GemmOut_t f03 = fvo[fcr_offset0+3];

  const GemmOut_t f10 = fvo[fcr_offset1+0];
  const GemmOut_t f11 = fvo[fcr_offset1+1];
  const GemmOut_t f12 = fvo[fcr_offset1+2];
  const GemmOut_t f13 = fvo[fcr_offset1+3];

  // Apply the permutation:

  // [ A  A ]  ->  [ A  B ]
  // [ A  A ]  ->  [ A  B ]
  // [ B  B ]  ->  [ A  B ]
  // [ B  B ]  ->  [ A  B ]

  const GemmOut_t f00p = f00;
  const GemmOut_t f01p = f01;

  const GemmOut_t f02p = f10;
  const GemmOut_t f03p = f11;

  const GemmOut_t f10p = f02;
  const GemmOut_t f11p = f03;

  const GemmOut_t f12p = f12;
  const GemmOut_t f13p = f13;

  // Use helper value to move value to upper half of mantissa.

  const double shifter = (((GMUInt32)1) << GM_TALLY1_MAX_VALUE_BITS);

  // Pack two 25-bit integers into mantissa of double.

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

  double* dvo = (double*)vo;

  dvo[dcr_offset0+0] = d00;
  dvo[dcr_offset0+1] = d01;

  dvo[dcr_offset1+0] = d10;
  dvo[dcr_offset1+1] = d11;
}

//-----------------------------------------------------------------------------
/// \brief Swizzle/cast values from cublas call into double complex format.
///
///        The cublas gemm poduces a matrix of scalars of 32 bit size
///        (int32 or float).  However the required format of the metrics
///        is a matrix of double complex values, with each double
///        containing two packed 25-bit integers.
///        This code does an in-place transformation from one to the other.

void gm_tc_repair_metrics_(
  int nvll,
  int nvl,
  void* vo_ptr,
  TCBufs& tc_bufs,
  GMEnv* env) {

  GMInsist(env && vo_ptr);
  GMInsist(nvll >= 0);
  GMInsist(nvl >= 0);
  GMInsist(nvll <= nvl);

  const int nvll2 = nvll / 2;

  const int threadblocksize = 256;
  const int vll2_threadblocks = gm_ceil_i8(nvll2, threadblocksize);

  const bool is_int8 = env->tc == 2;

  if (is_int8) {
    gm_tc_repair_metrics_kernel_<int32_t><<<
        dim3(vll2_threadblocks, nvl, 1),
        dim3(threadblocksize, 1, 1),
      0,
        env->stream_compute_>>>(nvl, nvll, nvll2, vo_ptr);
  } else {
    gm_tc_repair_metrics_kernel_<float><<<
        dim3(vll2_threadblocks, nvl, 1),
        dim3(threadblocksize, 1, 1),
      0,
        env->stream_compute_>>>(nvl, nvll, nvll2, vo_ptr);
  }

  GMEnv_cuda_last_call_succeeded(env);
}

//-----------------------------------------------------------------------------
/// \brief Use a standard GEMM to compute bitwise result: implementation.

void gm_tc_gemm_start_impl_(int m, int n, int k,
                            void* dA, int ldda,
                            void* dB, int lddb,
                            void* dC, int lddc,
                            TCBufs& tc_bufs,
                            GMEnv* env) {

  // Ensure tensor core hardware is available.
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, 0);
  GMInsist(deviceProp.major >= 7);

  const int I_max = m;
  const int I_max_dim = lddc;
  const int nvll = I_max_dim; // effective nvl for left matrix
  const int nvl = n;
  const int npvfl = k;
  const int num_steps = env->num_tc_steps;

  // Loop over steps of algorithm.
  for (int step_num = 0; step_num < num_steps; ++step_num) {

    // Select the block row of the left and right matrices for this step.
    const int pvfl_min = ((step_num+0) * npvfl) / num_steps;
    const int pvfl_max = ((step_num+1) * npvfl) / num_steps;
    const int npvfl_step = pvfl_max - pvfl_min;

    if (npvfl_step == 0) {  // empty block row
      continue;
    }

    // Convert the input matrices of packed bit values into matrices
    // of values of a type suitable for the GEMM.
    const bool left_matrix = false; // A
    const bool right_matrix = true; // B
    gm_tc_buf_write_(left_matrix, I_max, I_max_dim, nvl, npvfl,
                     npvfl_step, pvfl_min, dA, tc_bufs, env);
    gm_tc_buf_write_(right_matrix, I_max, I_max_dim, nvl, npvfl,
                     npvfl_step, pvfl_min, dB, tc_bufs, env);

    // Perform the GEMM for this pair of block rows; accumulate.
    gm_tc_solve_(pvfl_min==0, nvll, nvl, npvfl_step, dA, dB, dC, tc_bufs, env);
  }

  // Revise the results of the GEMMs to be in the needed double complex format.
  gm_tc_repair_metrics_(nvll, nvl, dC, tc_bufs, env);
}

//-----------------------------------------------------------------------------
#endif // USE_TC
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/// \brief Use a standard GEMM to compute bitwise result.

void gm_tc_gemm_start(int m, int n, int k,
                      void* dA, int ldda,
                      void* dB, int lddb,
                      void* dC, int lddc,
                      TCBufs& tc_bufs,
                      GMEnv* env) {
  GMInsist(dA && dB && dC && env);
  GMInsist(m >= 0 && n >= 0 && k >= 0);
  GMInsist(ldda >= 0 && lddb >= 0 && lddc >= 0);
  GMInsist(k <= ldda);
  GMInsist(k <= lddb);
  GMInsist(m <= lddc);
  GMInsist(env->tc);
  GMInsist(GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC);
  GMInsist(GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU);
#ifndef USE_TC
  GMInsistInterface(env,
                    false && "TC option unavailable for this platform/build.");
#endif

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, 0);
  GMInsist(deviceProp.major >= 7);

  gm_tc_gemm_start_impl_(m, n, k, dA, ldda, dB, lddb, dC, lddc, tc_bufs,  env);
}

//-----------------------------------------------------------------------------
/// \brief Initialize TCBufs object by allocating memory etc.

void gm_tc_bufs_malloc(int num_vector_local,
                       int num_field_local,
                       int num_packedval_field_local,
                       TCBufs& tc_bufs,
                       GMEnv* env) {
  GMInsist(env);
  GMInsist(num_vector_local >= 0);
  GMInsist(num_packedval_field_local >= 0);
  GMInsist(!tc_bufs.tc_buf_left);
  GMInsist(!tc_bufs.tc_buf_right);

  if (!env->tc) {
    return;
  }

  if (GMEnv_metric_type(env) != GM_METRIC_TYPE_CCC) {
    return;
  }

  // Calculate sizes.

  const size_t nvl = num_vector_local;
  const size_t npvfl = num_packedval_field_local;
  const size_t npvfl_step_max = gm_ceil_i8(npvfl, env->num_tc_steps);

  const bool is_int8 = env->tc == 2;
  const int sizeof_scalar = is_int8 ? 1 : 2;

  const size_t nvlX2 = nvl * 2;

  tc_bufs.tc_buf_size = nvlX2 * (npvfl_step_max * 64) * sizeof_scalar;
  tc_bufs.tc_buf_size = tc_bufs.tc_buf_size ? tc_bufs.tc_buf_size : 1;

  // Allocate buffers.

  cudaMalloc(&tc_bufs.tc_buf_left, tc_bufs.tc_buf_size);
  GMEnv_cuda_last_call_succeeded(env);
  env->gpu_mem += tc_bufs.tc_buf_size;
  env->gpu_mem_max = gm_max_i8(env->gpu_mem_max, env->gpu_mem);

  cudaMalloc(&tc_bufs.tc_buf_right, tc_bufs.tc_buf_size);
  GMEnv_cuda_last_call_succeeded(env);
  env->gpu_mem += tc_bufs.tc_buf_size;
  env->gpu_mem_max = gm_max_i8(env->gpu_mem_max, env->gpu_mem);

  // Set up cublas handle.

  cublasStatus_t status_cb = cublasCreate(&tc_bufs.cublas_handle);
  GMInsist(status_cb == CUBLAS_STATUS_SUCCESS);

  status_cb = cublasSetStream(tc_bufs.cublas_handle, env->stream_compute_);
  GMInsist(status_cb == CUBLAS_STATUS_SUCCESS);

#ifdef USE_TC
  status_cb = cublasSetMathMode(tc_bufs.cublas_handle, CUBLAS_TENSOR_OP_MATH);
  GMInsist(status_cb == CUBLAS_STATUS_SUCCESS);
#else
  GMInsistInterface(env,
                    false && "TC option unavailable for this platform/build.");
#endif
}

//-----------------------------------------------------------------------------
/// \brief Terminate TCBufs object by deallocating memory etc.

void gm_tc_bufs_free(TCBufs& tc_bufs,
                     GMEnv* env) {
  GMInsist(env);
  GMInsist((tc_bufs.tc_buf_left != 0) == (tc_bufs.tc_buf_right != 0));

  if (!tc_bufs.tc_buf_left) {
    return;
  }

  // Free buffers.

  cudaFree(tc_bufs.tc_buf_left);
  GMEnv_cuda_last_call_succeeded(env);
  tc_bufs.tc_buf_left = NULL;
  env->gpu_mem -= tc_bufs.tc_buf_size;

  cudaFree(tc_bufs.tc_buf_right);
  GMEnv_cuda_last_call_succeeded(env);
  tc_bufs.tc_buf_right = NULL;
  env->gpu_mem -= tc_bufs.tc_buf_size;

  // Free cublas handle.

  cublasStatus_t status_cb = cublasDestroy(tc_bufs.cublas_handle);
  GMInsist(status_cb == CUBLAS_STATUS_SUCCESS);
}

//-----------------------------------------------------------------------------
