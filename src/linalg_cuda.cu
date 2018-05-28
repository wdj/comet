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

#if 0
__device__ static const GMUInt32 MGEMM2_table1[] = {
  0x00000000,
  0x00003c00,
  0x00003c00,
  0x00004000,

  0x3c000000,
  0x3c003c00,
  0x3c003c00,
  0x3c004000,

  0x3c000000,
  0x3c003c00,
  0x3c003c00,
  0x3c004000,

  0x40000000,
  0x40003c00,
  0x40003c00,
  0x40004000
  };

__device__ static const GMUInt32 MGEMM2_table0[] = {
  0x40004000,
  0x40003c00,
  0x40003c00,
  0x40000000,

  0x3c004000,
  0x3c003c00,
  0x3c003c00,
  0x3c000000,

  0x3c004000,
  0x3c003c00,
  0x3c003c00,
  0x3c000000,

  0x00004000,
  0x00003c00,
  0x00003c00,
  0x00000000
  };

struct MGEMM2 {
  static const GMUInt32* table[2];
};

const GMUInt32* MGEMM2::table[] = {MGEMM2_table0, MGEMM2_table1};

#endif

#define TRANSPOSE

//-----------------------------------------------------------------------------

// https://docs.nvidia.com/cuda/cublas/index.html#cublas-gemmEx

//-----------------------------------------------------------------------------

__global__ void gm_tc_buf_write_fp16_(
  int num_way,
  bool is_sparse,
  int left_right,
  GMUInt32* vi,
  int vi_dim0,
  int nvl_in,
  int nvl,
  int nvl2,
  int nfl,
  int nfl2,
  GMUInt32* vo) {

  // Two fields (seminibbles) map to two half words of 32-bit word

  const int fl2 = threadIdx.x + blockIdx.x * blockDim.x;
  const int i01 = blockIdx.y; // count either 0 bits or 1 bits.
  const int vl = blockIdx.z;
  const int vl_in = vl <= nvl_in-1 ? vl : nvl_in-1;

  if (fl2 >= nfl2) {
    return;
  }

  if (vl < 0 || vl >= nvl) {
    return;
  }

  // Output array as floats has nfl/2 rows, as halfs has nfl rows.

  // Column offset in input array.

  const GMUInt32 * const vi_col = vi + vl_in * vi_dim0;

  // Pick up two consecutive field values.
  // NOTE: first field seminibble0, second field seminibble1

  const int nibble = (vi_col[fl2/8] >> (4 * (fl2%8))) & 15;

  const int seminibble0 = nibble & 3;
  const int seminibble1 = (nibble>>2) & 3;

  // Count number of 0 (or 1) bits in respective seminibble.

  const bool is_right = left_right != 0;

  const bool skip_10 = is_sparse || (num_way == 3 && ! is_right);

  const GMUInt16 zero = 0x0000;
  const GMUInt16 one = 0x3c00;
  const GMUInt16 two = 0x4000;

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
  // Combine two halfs into one float.

  const GMUInt32 out01 = ((GMUInt32)out0 + ( ((GMUInt32out1) << 16);

  // Always keep pair of cols together, corresponding to the two i01 values.
  // Right case: straight copy of cols to cols in sequence.
  // Left case: interleave to make swizzling of result array work:
  // [ A A B B C C D D E E F F ] -> [ A A D D B B E E C C F F]

  const int vlX2_index = is_right ? i01 + 2*vl :
                  i01 + 2*( vl < nvl2 ? 2*vl : 2*vl - nvl + 1 );

  const int fl2_index = fl2;
  const int fl2_dim = nfl2;

  vo[ fl2_index + fl2_dim * vlX2_index ] = out01;

#ifdef TRANSPOSE
// switch?
  const int fl_index_0 = 0 + 2 * fl2_index;
  const int fl_index_1 = 1 + 2 * fl2_index;

  const int vlX2_dim = 2 * nvl;

  ((GMUInt16*)vo)[ vlX2_index + vlX2_dim * fl_index_0 ] = out0;
  ((GMUInt16*)vo)[ vlX2_index + vlX2_dim * fl_index_1 ] = out1;
#endif

//if (fl2==0) printf("%s vec %i field %i  %i\n", left_right ? "r" : "l", vl, 2*fl2+0, seminibble0);
//if (fl2==0) printf("lr %i vec %i field %i  %i\n", left_right, vl, 2*fl2+1, seminibble1);
}

//-----------------------------------------------------------------------------

__global__ void gm_tc_buf_write_int8_(
  int num_way,
  bool is_sparse,
  int left_right,
  GMUInt32* vi,
  int vi_dim0,
  int nvl_in,
  int nvl,
  int nvl2,
  int nfl,
  int nfl2,
  GMUInt16* vo) {

  // Two fields (seminibbles) map to two half words of 32-bit word

  const int fl2 = threadIdx.x + blockIdx.x * blockDim.x;
  const int i01 = blockIdx.y; // count either 0 bits or 1 bits.
  const int vl = blockIdx.z;
  const int vl_in = vl <= nvl_in-1 ? vl : nvl_in-1;

  if (fl2 >= nfl2) {
    return;
  }

  if (vl < 0 || vl >= nvl) {
    return;
  }

  // Output array as shorts has nfl/2 rows, as chars has nfl rows.

  // Column offset in input array.

  const GMUInt32 * const vi_col = vi + vl_in * vi_dim0;

  // Pick up two consecutive field values.
  // NOTE: first field seminibble0, second field seminibble1

  const int nibble = (vi_col[fl2/8] >> (4 * (fl2%8))) & 15;

  const int seminibble0 = nibble & 3;
  const int seminibble1 = (nibble>>2) & 3;

  // Count number of 0 (or 1) bits in respective seminibble.

  const bool is_right = left_right != 0;

  // this needs work
  const bool skip_10 = is_sparse || (num_way == 3 && ! is_right);

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
  // Combine two chars into one short int.

  const GMUInt16 out01 = ((GMUInt16)out0) + ( ((GMUInt16)out1) << 8);

  // Always keep pair of cols together, corresponding to the two i01 values.
  // Right case: straight copy of cols to cols in sequence.
  // Left case: interleave to make swizzling of result array work:
  // [ A A B B C C D D E E F F ] -> [ A A D D B B E E C C F F]

  const int vlX2_index = is_right ? i01 + 2*vl :
                  i01 + 2*( vl < nvl2 ? 2*vl : 2*vl - nvl + 1 );

  const int fl2_index = fl2;
  const int fl2_dim = nfl2;

  vo[ fl2_index + fl2_dim * vlX2_index ] = out01;

// byte address {0,1} + 2 * fl2 

#ifdef TRANSPOSE
// switch?
  const int fl_index_0 = 0 + 2 * fl2_index;
  const int fl_index_1 = 1 + 2 * fl2_index;

  const int vlX2_dim = 2 * nvl;

  ((GMUInt8*)vo)[ vlX2_index + vlX2_dim * fl_index_0 ] = out0;
  ((GMUInt8*)vo)[ vlX2_index + vlX2_dim * fl_index_1 ] = out1;
#endif

//  if (fl2==0 && vl==0 && i01==0) printf("//////  %i\n", (int) *(GMUInt16*)vo);

//  if (fl2==0)
//    printf("////// %s %i %i  %i %i %i %i  %i %i %i  %i %i\n", left_right ? "r" : "l", (int)nvl, (int)nfl,  fl0, fl1, i01, vl,  fl0, fl1, col,  (int)(col+2*nvl*fl0), (int)(col+2*nvl*fl1));

//printf("%i %i\n", (int)out0, (int)out1);

//if (fl2==0) printf("%s i01 %i vec %i field %i  %i  %i\n", left_right ? "r" : "l", i01, vl, 2*fl2+0, seminibble0, (int)out0);

//if (seminibble0) printf("vec %i field %i  %i\n", vl, 2*fl2+0, seminibble0);
//if (seminibble1) printf("vec %i field %i  %i\n", vl, 2*fl2+1, seminibble1);
}

//-----------------------------------------------------------------------------
// Convert packed-bits fields of input vectors into halfs (or chars).

void gm_tc_buf_write(
  int left_right,
  int I_max,
  int num_vector_local,
  int num_packedval_field_local,
  void* bufd,
  GMEnv* env) {
  GMInsist(left_right == 0 || left_right == 1);
  GMInsist(env && bufd);
  GMInsist(I_max >= 0);
  GMInsist(num_vector_local >= 0);
  GMInsist(num_packedval_field_local >= 0);
  GMInsist(GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC);
  GMInsist(GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU);

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, 0);
  GMInsist(deviceProp.major >= 7);

  //GMInsistInterface(env, GMEnv_num_way(env) == GM_NUM_WAY_2 &&
  //                  "Not yet implemented.");
  //GMInsistInterface(env, !env->sparse && "Not yet implemented.");

  const int nvl = num_vector_local;
  const int nvl2 = nvl / 2;
  const int npvfl = num_packedval_field_local;
  const int nfl = npvfl * 64;
  const int nfl2 = nfl / 2;

  GMInsistInterface(env, num_vector_local % 2 == 0 &&
                    "tc method requires num_vector_local multiple of 2.");

  const int threadblocksize = 256;
  const int fl2_threadblocks = (nfl2+threadblocksize-1) / threadblocksize;

  const bool is_right = left_right != 0;
  const bool is_int8 = env->tc == 2;

  if (is_int8) {

    gm_tc_buf_write_int8_<<<
      dim3(fl2_threadblocks, 2, nvl),
      dim3(threadblocksize, 1, 1),
      0,
      env->stream_compute_>>>(
      GMEnv_num_way(env),
      env->sparse,
      left_right,
      (GMUInt32*)bufd,
      npvfl * 4,
      left_right ? nvl : I_max,
      nvl,
      nvl2,
      nfl,
      nfl2,
      is_right ? (GMUInt16*)env->tc_buf_right : (GMUInt16*)env->tc_buf_left);

  } else {

    gm_tc_buf_write_fp16_<<<
      dim3(fl2_threadblocks, 2, nvl),
      dim3(threadblocksize, 1, 1),
      0,
      env->stream_compute_>>>(
      GMEnv_num_way(env),
      env->sparse,
      left_right,
      (GMUInt32*)bufd,
      npvfl * 4,
      left_right ? nvl : I_max,
      nvl,
      nvl2,
      nfl,
      nfl2,
      is_right ? (GMUInt32*)env->tc_buf_right : (GMUInt32*)env->tc_buf_left);

  }

  GMEnv_cuda_last_call_succeeded(env);
}

//-----------------------------------------------------------------------------
// Call tensor core enabled cuBLAS function to tally bits for CCC.

void gm_tc_solve(
  int I_max,
  int num_vector_local,
  int num_packedval_field_local,
  void* dA,
  int ldda,
  void* dB,
  int lddb,
  void* dC,
  int lddc,
  GMEnv* env) {
  GMInsist(env && dA && dB && dC);
  GMInsist(I_max >= 0);
  GMInsist(num_vector_local >= 0);
  GMInsist(num_packedval_field_local >= 0);
  GMInsist(ldda >= 0 && lddb >= 0 && lddc >= 0);
  GMInsist(env->tc == 1 || env->tc == 2);

  const int nvl = num_vector_local;
  const int npvfl = num_packedval_field_local;
  const int nfl = npvfl * 64;

  const int m = 2 * nvl; // metrics array dim
  const int n = 2 * nvl; // metrics array dim
  const int k = nfl; // vectors array (as halfs) dim

  const float alpha = 1;
  const float beta = 0;

  const bool is_int8 = env->tc == 2;

  GMInsist(k % 8 == 0); // nfl is derived from padded-up npvfl, so ok.
  GMInsist(m % 8 == 0); // need nvl % 4 == 0

#if 0
  cublasStatus_t status = cublasSgemmEx(
    env->cublas_handle,
    CUBLAS_OP_T, CUBLAS_OP_N,
    m, n, k,
    &alpha,
    env->tc_buf_left, CUDA_R_16F, k,
    env->tc_buf_right, CUDA_R_16F, k,
    &beta,
    dC, CUDA_R_32F, m);
#endif

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
    CUBLAS_GEMM_DFALT_TENSOR_OP);
    //CUBLAS_GEMM_ALGO4_TENSOR_OP); // best, for cuda 9.1.85 non-transpose

  if (status == CUBLAS_STATUS_NOT_INITIALIZED) {
    printf("Error: CUBLAS_STATUS_NOT_INITIALIZED\n");
  }
  if (status == CUBLAS_STATUS_ARCH_MISMATCH) {
    printf("Error: CUBLAS_STATUS_ARCH_MISMATCH\n");
  }
  if (status == CUBLAS_STATUS_NOT_SUPPORTED) {
    printf("Error: CUBLAS_STATUS_NOT_SUPPORTED\n");
  }
  if (status == CUBLAS_STATUS_INVALID_VALUE) {
    printf("Error: CUBLAS_STATUS_INVALID_VALUE\n");
  }
  if (status == CUBLAS_STATUS_EXECUTION_FAILED) {
    printf("Error: CUBLAS_STATUS_EXECUTION_FAILED\n");
  }

  GMInsist(status == CUBLAS_STATUS_SUCCESS);

  env->ops_local += 2 * m * (double)n * (double)k;
}

//-----------------------------------------------------------------------------

__global__ void gm_tc_fix_metrics_(
  int nvl,
  int nvl2,
  float* bufd) {

  // Row and column of metrics array.

  const int thread_r = threadIdx.x + blockIdx.x * blockDim.x;
  const int thread_c = blockIdx.y;

  if (thread_r >= nvl2 || thread_c >= nvl) {
    return;
  }

  // Considered as an array of floats, array is 2*nvl rows X 2*nvl cols.
  // Each thread manipulates a block of 4 rows and 2 cols.
  // Thus the dimensions of the metrics array in blocks is nvl2 X nvl.
  // Each block viewed as an array of doubles is 2 X 2.

  // Two col numbers being processed of this (float) array.

  const int fc0 = thread_c * (4*nvl);
  const int fc1 = thread_c * (4*nvl) + 2*nvl;

  // Read the 8 floats.

  const float f00 = bufd[fc0+0+4*thread_r];
  const float f01 = bufd[fc0+1+4*thread_r];
  const float f02 = bufd[fc0+2+4*thread_r];
  const float f03 = bufd[fc0+3+4*thread_r];

  const float f10 = bufd[fc1+0+4*thread_r];
  const float f11 = bufd[fc1+1+4*thread_r];
  const float f12 = bufd[fc1+2+4*thread_r];
  const float f13 = bufd[fc1+3+4*thread_r];

  // Apply the permutation.
  //
  // [ A  A ]      [ A  B ]
  // [ A  A ]  ->  [ A  B ]
  // [ B  B ]  ->  [ A  B ]
  // [ B  B ]      [ A  B ]

  const float f00p = f00;
  const float f01p = f01;

  const float f02p = f10;
  const float f03p = f11;

  const float f10p = f02;
  const float f11p = f03;

  const float f12p = f12;
  const float f13p = f13;

  // Pack two 25-bit integers into mantissa of double.

  const double shifter = (((GMUInt32)1)<<GM_TALLY1_MAX_VALUE_BITS);

  const double d00 = f00p + f02p * shifter;
  const double d01 = f01p + f03p * shifter;

  const double d10 = f10p + f12p * shifter;
  const double d11 = f11p + f13p * shifter;

  // Overwrite block with the new values.

  const int dc0 = thread_c * (2*nvl);
  const int dc1 = thread_c * (2*nvl) + nvl;

  ((double*)bufd)[dc0+0+2*thread_r] = d00;
  ((double*)bufd)[dc0+1+2*thread_r] = d01;

  ((double*)bufd)[dc1+0+2*thread_r] = d10;
  ((double*)bufd)[dc1+1+2*thread_r] = d11;

//printf("%i %i %f %f %f %f\n", 2*thread_r+0, thread_c, f00p, f01p, f02p, f03p);
//printf("%i %i %f %f %f %f\n", 2*thread_r+1, thread_c, f10p, f11p, f12p, f13p);
}

//-----------------------------------------------------------------------------
// Swizzle/cast values from the CUBLAS call into required double complex format.

void gm_tc_fix_metrics(
  int I_max,
  int num_vector_local,
  void* bufd,
  GMEnv* env) {
  GMInsist(env && bufd);
  GMInsist(I_max >= 0);
  GMInsist(num_vector_local >= 0);

  const int nvl = num_vector_local;
  const int nvl2 = nvl / 2;

  const int threadblocksize = 256;
  const int vl2_threadblocks = (nvl2+threadblocksize-1) / threadblocksize;

  gm_tc_fix_metrics_<<<
    dim3(vl2_threadblocks, nvl, 1),
    dim3(threadblocksize, 1, 1),
    0,
    env->stream_compute_>>>(
    nvl,
    nvl2,
    (float*)bufd
  );

  GMEnv_cuda_last_call_succeeded(env);
}

//-----------------------------------------------------------------------------
