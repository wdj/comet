//-----------------------------------------------------------------------------
/*!
 * \file   linalg_cuda.cu
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  Supporting CUDA functions.
 * \note   Copyright (C) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "linalg_cuda.hh"

//-----------------------------------------------------------------------------

#define TRANSPOSE

// If TRANSPOSE, then copied matrices are stored as A, B^T, for speed on
// the tensor cores.  Otherwise, store as the rest of the code: A^T, B.

#ifdef USE_TC

#include "cuda_fp16.h"

//TODO: get rid of multiplier
#define MULTIPLIER 1

//-----------------------------------------------------------------------------

__global__ void gm_tc_buf_write_fp16_kernel_(
  int num_way,
  bool is_sparse,
  bool is_right,
  GMUInt32* vi,
  int vi_dim0,
  int nvlea,
  int nvle,
  int nvle2,
  int nvleX2,
  int nfl,
  int nfl2,
  int nfl2_step,
  int fl2_min,
  GMUInt32* vo) {
//TODO: rename nvle->nvl, nvlea->nvla . . .
//TODO: vi32

  // Two fields (seminibbles) map to two halves of 32-bit word

#ifdef TRANSPOSE
  const int vlX2 = threadIdx.x + blockIdx.x * blockDim.x;
  const int fl2_step = blockIdx.y;
#else
  const int fl2_step = threadIdx.x + blockIdx.x * blockDim.x;
  const int vlX2 = blockIdx.y;
#endif

  if (vlX2 >= nvleX2 || fl2_step >= nfl2_step) {
    return;
  }

  const int i01 = vlX2 % 2; // count either 0 bits or 1 bits.
  const int vl = vlX2 / 2;

  const int fl2 = fl2_min + fl2_step;

  // Output array as floats has nfl/2 rows, as halfs has nfl rows.

//TODO: cast vi_dim0 to size_t
  const GMUInt32* const vi_col = vi + vl * vi_dim0;

  // Pick up two consecutive field values:
  // first field seminibble0, second field seminibble1
  // Set to zero if outside of active range.

  const int nibble = vl<nvlea ? (vi_col[fl2/8] >> (4*(fl2%8))) & 15 : 0;

  const int seminibble0 = nibble & 3;
  const int seminibble1 = (nibble>>2) & 3;

  // Count number of 0 (or 1) bits in respective seminibble.
  // Determine whether to skip (1,0) null indicator value.

  //CHECK
  const bool skip_10 = is_sparse || (num_way == 3 && ! is_right);

  // Possible counts, represented as FP16.
  const GMUInt16 zero = 0x0000;
  const GMUInt16 one = 0x3c00;
  const GMUInt16 two = 0x4000;
  //const GMUInt16 zero = *(GMUInt16*)&__float2half(0.);
  //const GMUInt16 one = *(GMUInt16*)&__float2half(1.);
  //const GMUInt16 two = *(GMUInt16*)&__float2half(2.);

  const GMUInt16 out0 = seminibble0 == 3*i01     ? two :
                        seminibble0 == 3*(1-i01) ? zero :
                                       !skip_10  ? one :
                        seminibble0 == 1         ? one :
                                                   zero;

  const GMUInt16 out1 = seminibble1 == 3*i01     ? two :
                        seminibble1 == 3*(1-i01) ? zero :
                                       !skip_10  ? one :
                        seminibble1 == 1         ? one :
                                                   zero;
  // Always keep pair of cols together, corresponding to the two i01 values.
  // Right case: straight copy of cols to cols in sequence.
  // Left case: interleave to make later swizzling of metrics array work:
  // [ A A B B C C D D E E F F ] -> [ A A D D B B E E C C F F]

//TODO: vl_index = is_right ? vl : vl < nvle2 ? 2*vl : 2*vl - nvle + 1
  const int vlX2_index = is_right ? i01 + 2*vl :
                  i01 + 2*( vl < nvle2 ? 2*vl : 2*vl - nvle + 1 );

  const int fl2_index = fl2_step;

#ifdef TRANSPOSE
  const int fl_index_0 = 0 + 2 * fl2_index;
  const int fl_index_1 = 1 + 2 * fl2_index;

  const int vlX2_dim = nvleX2;

//TODO: void* vo, GMUInt16* vo16;
  ((GMUInt16*)vo)[vlX2_index + vlX2_dim * (size_t)fl_index_0] = out0;
  ((GMUInt16*)vo)[vlX2_index + vlX2_dim * (size_t)fl_index_1] = out1;

#else
  const int fl2_dim = nfl2_step;

  // Combine two halfs into one 32-bit value.

  const GMUInt32 out01 = ((GMUInt32)out0) + ( ((GMUInt32)out1) << 16 );

//TODO: GMUInt32* vo32
//TODO: cast vlX2_index to size_t
  vo[fl2_index + fl2_dim * vlX2_index] = out01;
#endif
}

//-----------------------------------------------------------------------------

__global__ void gm_tc_buf_write_int8_kernel_(
  int num_way,
  bool is_sparse,
  bool is_right,
  GMUInt32* vi,
  int vi_dim0,
  int nvlea,
  int nvle,
  int nvle2,
  int nvleX2,
  int nfl,
  int nfl2,
  int nfl2_step,
  int fl2_min,
  GMUInt16* vo) {

  // Two fields (seminibbles) map to two halves of 16-bit word

#ifdef TRANSPOSE
  const int vlX2 = threadIdx.x + blockIdx.x * blockDim.x;
  const int fl2_step = blockIdx.y;
#else
  const int fl2_step = threadIdx.x + blockIdx.x * blockDim.x;
  const int vlX2 = blockIdx.y;
#endif

  if (vlX2 >= nvleX2 || fl2_step >= nfl2_step) {
    return;
  }

  const int i01 = vlX2 % 2; // count either 0 bits or 1 bits.
  const int vl = vlX2 / 2;

  const int fl2 = fl2_min + fl2_step;

  // Output array as shorts has nfl/2 rows, as chars has nfl rows.

// ISSUE; int32 vs int16
  const GMUInt32* const vi_col = vi + vl * vi_dim0;

  // Pick up two consecutive field values:
  // first field seminibble0, second field seminibble1
  // Set to zero if outside of active range.

  const int nibble = vl<nvlea ? (vi_col[fl2/8] >> (4*(fl2%8))) & 15 : 0; 

  const int seminibble0 = nibble & 3;
  const int seminibble1 = (nibble>>2) & 3;

  // Count number of 0 (or 1) bits in respective seminibble.
  // Determine whether to skip (1,0) null indicator value.

  //CHECK
  const bool skip_10 = is_sparse || (num_way == 3 && ! is_right);

  // Possibe counts, represented as Int8.
  const GMUInt8 zero = 0;
  const GMUInt8 one = 1;
  const GMUInt8 two = 2;

  const GMUInt8 out0 = seminibble0 == 3*i01     ? two :
                       seminibble0 == 3*(1-i01) ? zero :
                                      !skip_10  ? one :
                       seminibble0 == 1         ? one :
                                                  zero;

  const GMUInt8 out1 = seminibble1 == 3*i01     ? two :
                       seminibble1 == 3*(1-i01) ? zero :
                                      !skip_10  ? one :
                       seminibble1 == 1         ? one :
                                                  zero;
  // Always keep pair of cols together, corresponding to the two i01 values.
  // Right case: straight copy of cols to cols in sequence.
  // Left case: interleave to make later swizzling of metrics array work:
  // [ A A B B C C D D E E F F ] -> [ A A D D B B E E C C F F]

  const int vlX2_index = is_right ? i01 + 2*vl :
                  i01 + 2*( vl < nvle2 ? 2*vl : 2*vl - nvle + 1 );

  const int fl2_index = fl2_step;

#ifdef TRANSPOSE
  //CHECK
  const int fl_index_0 = 0 + 2 * fl2_index;
  const int fl_index_1 = 1 + 2 * fl2_index;

  const int vlX2_dim = nvleX2;

  ((GMUInt8*)vo)[ vlX2_index + vlX2_dim * fl_index_0 ] = out0;
  ((GMUInt8*)vo)[ vlX2_index + vlX2_dim * fl_index_1 ] = out1;
#else
  const int fl2_dim = nfl2_step;

  // Combine two chars into one short int.

  const GMUInt16 out01 = ((GMUInt16)out0) + ( ((GMUInt16)out1) << 8 );

  vo[ fl2_index + fl2_dim * vlX2_index ] = out01;
#endif
}

//-----------------------------------------------------------------------------
// Convert matrix stored as packed 2-bit values into matrix of FP16 (or Int8).

void gm_tc_buf_write_(
  bool is_right,
  int I_max,
  int I_max_dim,
  int nvl,
  int npvfl,
  int npvfl_step,
  int pvfl_min,
  void* vi_ptr,
  GMEnv* env) {
  GMInsist(env && vi_ptr);
  GMInsist(I_max_dim >= 0 && I_max_dim <= nvl);
  GMInsist(I_max >= 0 && I_max <= I_max_dim);
  GMInsist(nvl >= 0);
  GMInsist(npvfl >= 0);
//TODO: more assertions
//TODO: void: vi, int32* vi32

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

// ISSUE: is this right for int16 case
  const int vi_dim0 = npvfl * 4; // 4 = sizeof(doublecomplex) / sizeof(int32)

  GMInsistInterface(env, nvle % 2 == 0 && nvl % 2 == 0 &&
                    "tc method here requires num_vector_local multiple of 2.");

  const int threadblocksize = 256;
#ifdef TRANSPOSE
  const int num_threadblocks_0 = gm_ceil_i8(nvleX2, threadblocksize);
  const int num_threadblocks_1 = nfl2_step;
#else
  const int num_threadblocks_0 = gm_ceil_i8(nfl2_step, threadblocksize);
  const int num_threadblocks_1 = nvleX2;
#endif

  const void* tc_buf = is_right ? env->tc_buf_right : env->tc_buf_left;

  if (! is_int8) {

    gm_tc_buf_write_fp16_kernel_<<<
        dim3(num_threadblocks_0, num_threadblocks_1, 1),
        dim3(threadblocksize, 1, 1),
        0,
        env->stream_compute_>>>(
      GMEnv_num_way(env),
      env->sparse,
      is_right,
      (GMUInt32*)vi_ptr,
      vi_dim0,
      nvlea,
      nvle,
      nvle2,
      nvleX2,
      nfl,
      nfl2,
      nfl2_step,
      fl2_min,
      (GMUInt32*)tc_buf);

  } else {

    gm_tc_buf_write_int8_kernel_<<<
        dim3(num_threadblocks_0, num_threadblocks_1, 1),
        dim3(threadblocksize, 1, 1),
        0,
        env->stream_compute_>>>(
      GMEnv_num_way(env),
      env->sparse,
      is_right,
      (GMUInt32*)vi_ptr,
      vi_dim0,
      nvlea,
      nvle,
      nvle2,
      nvleX2,
      nfl,
      nfl2,
      nfl2_step,
      fl2_min,
      (GMUInt16*)tc_buf);

  }

  GMEnv_cuda_last_call_succeeded(env);
}

//-----------------------------------------------------------------------------
// Call tensor core enabled cuBLAS function to tally bits for CCC.

void gm_tc_solve_(
  bool is_first,
  int nvll,
  int nvl,
  int npvfl_step,
  void* dA,
  void* dB,
  void* dC,
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

  const float alpha = 1;
  const float beta = is_first ? 0 : 1;

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
    env->cublas_handle,
#ifdef TRANSPOSE
    CUBLAS_OP_N, CUBLAS_OP_T,
#else
    CUBLAS_OP_T, CUBLAS_OP_N,
#endif
    m, n, k,
    &alpha,
    env->tc_buf_left, is_int8 ? CUDA_R_8I : CUDA_R_16F,
#ifdef TRANSPOSE
    m,
#else
    k,
#endif
    env->tc_buf_right, is_int8 ? CUDA_R_8I : CUDA_R_16F,
#ifdef TRANSPOSE
    n,
#else
    k,
#endif
    &beta,
    dC, CUDA_R_32F, m,
    CUDA_R_32F,
#ifdef TRANSPOSE
    //CUBLAS_GEMM_ALGO3_TENSOR_OP // best timing, for cuda 9.1.85, transpose
    //CUBLAS_GEMM_DFALT_TENSOR_OP // good timing, for cuda 9.2.88, transpose
    CUBLAS_GEMM_ALGO4_TENSOR_OP // best timing, for cuda 9.2.88, transpose
#else
    CUBLAS_GEMM_ALGO4_TENSOR_OP // best timing, for cuda 9.1.85, non-transpose
#endif
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

__global__ void gm_tc_fix_metrics_kernel_(
  int nvl,
  int nvll,
  int nvll2,
  float* vo,
  float multiplier) {

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

//TODO: use only single pointer instead . . . ?

//TODO: make size_t; cast 4*nvll to size_t
  const int fc_offset0 = thread_c * (4*nvll);
  const int fc_offset1 = thread_c * (4*nvll) + 2*nvll;
//TODO: fcr_offset0/1

  // Read the 8 floats.

//TODO: variable fvo
  const float f00 = vo[fc_offset0+0+4*thread_r] * multiplier;
  const float f01 = vo[fc_offset0+1+4*thread_r] * multiplier;
  const float f02 = vo[fc_offset0+2+4*thread_r] * multiplier;
  const float f03 = vo[fc_offset0+3+4*thread_r] * multiplier;

  const float f10 = vo[fc_offset1+0+4*thread_r] * multiplier;
  const float f11 = vo[fc_offset1+1+4*thread_r] * multiplier;
  const float f12 = vo[fc_offset1+2+4*thread_r] * multiplier;
  const float f13 = vo[fc_offset1+3+4*thread_r] * multiplier;

  // Apply the permutation:

  // [ A  A ]  ->  [ A  B ]
  // [ A  A ]  ->  [ A  B ]
  // [ B  B ]  ->  [ A  B ]
  // [ B  B ]  ->  [ A  B ]

  const float f00p = f00;
  const float f01p = f01;

  const float f02p = f10;
  const float f03p = f11;

  const float f10p = f02;
  const float f11p = f03;

  const float f12p = f12;
  const float f13p = f13;

  // Use helper value to move value to upper half of mantissa.

  const double shifter = (((GMUInt32)1) << GM_TALLY1_MAX_VALUE_BITS);

  // Pack two 25-bit integers into mantissa of double.

// TODO: explicitly cast float values to double
  const double d00 = f00p + f02p * shifter;
  const double d01 = f01p + f03p * shifter;

  const double d10 = f10p + f12p * shifter;
  const double d11 = f11p + f13p * shifter;

  // Overwrite block with the new values.
  // All is isolated to a single thread, should be thread safe.

//TODO: make size_t; cast 2*nvll to size_t
  const int dc_offset0 = thread_c * (2*nvll);
  const int dc_offset1 = thread_c * (2*nvll) + nvll;
//TODO: dcr_offset0/1

//TODO: variable dvo
  ((double*)vo)[dc_offset0+0+2*thread_r] = d00;
  ((double*)vo)[dc_offset0+1+2*thread_r] = d01;

  ((double*)vo)[dc_offset1+0+2*thread_r] = d10;
  ((double*)vo)[dc_offset1+1+2*thread_r] = d11;
}

//-----------------------------------------------------------------------------
// Swizzle/cast values from the CUBLAS call into required double complex format.

void gm_tc_fix_metrics_(
  int nvll,
  int nvl,
  void* vo_ptr,
  GMEnv* env) {
  GMInsist(env && vo_ptr);
  GMInsist(nvll >= 0);
  GMInsist(nvl >= 0);
  GMInsist(nvll <= nvl);

  const int nvll2 = nvll / 2;

  const int threadblocksize = 256;
  const int vll2_threadblocks = gm_ceil_i8(nvll2, threadblocksize);

  const bool is_int8 = env->tc == 2;

  gm_tc_fix_metrics_kernel_<<<
      dim3(vll2_threadblocks, nvl, 1),
      dim3(threadblocksize, 1, 1),
      0,
      env->stream_compute_>>>(
    nvl,
    nvll,
    nvll2,
    (float*)vo_ptr,
    is_int8 ? 1 : MULTIPLIER * MULTIPLIER
  );

  GMEnv_cuda_last_call_succeeded(env);
}
#endif

//-----------------------------------------------------------------------------

void gm_tc_gemm_start(int m, int n, int k,
                      void* dA, int ldda,
                      void* dB, int lddb,
                      void* dC, int lddc,
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

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, 0);
  GMInsist(deviceProp.major >= 7);

#ifdef USE_TC
  const int I_max = m;
  const int I_max_dim = lddc;
  const int nvll = I_max_dim; // effective nvl for left matrix
  const int nvl = n;
  const int npvfl = k;

  const int num_steps = env->num_tc_steps;

  for (int step_num = 0; step_num < num_steps; ++step_num) {
    const int pvfl_min = ((step_num+0) * npvfl) / num_steps;
    const int pvfl_max = ((step_num+1) * npvfl) / num_steps;
    const int npvfl_step = pvfl_max - pvfl_min;
    GMAssert(npvfl_step <= env->npvfl_step_max);

    if (npvfl_step == 0) {
      continue;
    }

    const bool left_matrix = false; // A
    const bool right_matrix = true; // B
    gm_tc_buf_write_(left_matrix, I_max, I_max_dim, nvl, npvfl,
                     npvfl_step, pvfl_min, dA, env);
    gm_tc_buf_write_(right_matrix, I_max, I_max_dim, nvl, npvfl,
                     npvfl_step, pvfl_min, dB, env);

    //for (int i=0; i<20; ++i)
    gm_tc_solve_(pvfl_min==0, nvll, nvl, npvfl_step, dA, dB, dC, env);
  }

  gm_tc_fix_metrics_(nvll, nvl, dC, env);

#else
  GMInsistInterface(env,
                    false && "TC option not implemented for this platform.");
#endif
}

//-----------------------------------------------------------------------------
