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
#include "cstdio"

#if defined COMET_USE_CUDA
#  include "cublas_v2.h"
#  include "cuda_fp16.h"
#elif defined COMET_USE_HIP
#  include "hip/hip_runtime_api.h"
//#  pragma GCC diagnostic ignored "-Wc99-designator"
#  include "hip/hip_runtime.h"
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
#include "linalg_tc.hh"

//=============================================================================

namespace comet {

//=============================================================================
// HELPERS
//=============================================================================

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
  static __host__ __device__ uint16_t four() {return (uint16_t)0x4400;}
                                        // = *(uint16_t*)&__float2half(0.);
                                        // = *(uint16_t*)&__float2half(1.);
                                        // = *(uint16_t*)&__float2half(2.);
                                        // = *(uint16_t*)&__float2half(4.);
};

//----------

template<> struct TCBufTypes<int8_t> {
  static __host__ __device__ int8_t zero() {return (int8_t)0;}
  static __host__ __device__ int8_t one() {return (int8_t)1;}
  static __host__ __device__ int8_t two() {return (int8_t)2;}
  static __host__ __device__ int8_t four() {return (int8_t)4;}
};

//----------

template<> struct TCBufTypes<GMFp32> {
  static __host__ __device__ GMFp32 zero() {return (GMFp32)0;}
  static __host__ __device__ GMFp32 one() {return (GMFp32)1;}
  static __host__ __device__ GMFp32 two() {return (GMFp32)2;}
  static __host__ __device__ GMFp32 four() {return (GMFp32)4;}
};

//-----------------------------------------------------------------------------
/// \brief Select types etc. based on the setting of the tc param.

template<int TC_METHOD> struct TCSelector;

template<> struct TCSelector<TC::INT8> {
  // types.
  typedef int8_t GemmIn_t;
  typedef int32_t GemmOut_t;
#if defined COMET_USE_CUDA
  // type selector parameters.
  static cudaDataType __host__ __device__ gemm_type_in() {return CUDA_R_8I;}
  static cudaDataType __host__ __device__ gemm_type_out() {return CUDA_R_32I;}
#elif defined COMET_USE_HIP
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
#if defined COMET_USE_CUDA
  // type selector parameters.
  static cudaDataType __host__ __device__ gemm_type_in() {return CUDA_R_16F;}
  static cudaDataType __host__ __device__ gemm_type_out() {return CUDA_R_32F;}
#elif defined COMET_USE_HIP
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
#if defined COMET_USE_CUDA
  // type selector parameters.
  static cudaDataType __host__ __device__ gemm_type_in() {return CUDA_R_32F;}
  static cudaDataType __host__ __device__ gemm_type_out() {return CUDA_R_32F;}
#elif defined COMET_USE_HIP
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
// TC_BUF_WRITE FILE-LOCAL (STATIC) FUNCTIONS
//=============================================================================

//-----------------------------------------------------------------------------
/// \brief Helper function: is a nonnegative integer a power of 2.

__host__ __device__ static bool is_po2(int x) {
  return x && (!(x&(x-1))); 
}

//-----------------------------------------------------------------------------
/// \brief Write individual elements to buf.
///
/// Description of is_bitwise_3way_2step option:
///
/// For this method two passes are made (instead of three), each of which
/// calculates exactly 4 of the required 8 metrics values.
/// The main idea is this: the left matrix entries are calculated as
/// the combined number of paths through corresponding elements of the I
/// and J matrices.  For the first pass (step_2way == 0), paths 0-0 and
/// 0-1 (corresponding to I-J element values) are counted; for the second
/// pass, paths 1-0 and 1-1.
/// Below is a tabular representation of the values for the non-sparse CCC case.
/// The "10" cases are omtted here because they are the same as the "01" cases.
/// The 3-way CCC sparse case is easliy adapted from this.
///
/// step_2way ------>   0    0    1    1
/// i01       ------>   0    1    0    1
/// output    ------>  0-0  0-1  1-0  1-1
/// --------------------------------------
///  m (=I)   c (=J)
///  |        |
///  v        v
///  00       00        4    0    0    0
///  00       01        2    2    0    0
///  00       11        2    2    0    0
/// --------------------------------------
///  01       00        2    0    2    0
///  01       01        1    1    1    1
///  01       11        0    2    0    2
/// --------------------------------------
///  11       00        0    0    4    0
///  11       01        0    0    2    2
///  11       11        0    0    0    4
///
/// The 3-way non-sparse DUO method is similar, but we only
/// look at 1 bit of each seminibble, not 2:
///
/// step_2way ------>   0    0    1    1
/// i01       ------>   0    1    0    1
/// output    ------>  0-0  0-1  1-0  1-1
/// --------------------------------------
///  m (=I)   c (=J)
///  |        |
///  v        v
///  *0       *0        1    0    0    0
///  *0       *1        0    1    0    0
/// --------------------------------------
///  *1       *0        0    0    1    0
///  *1       *1        0    0    0    1

template<typename GemmIn_t>
__host__ __device__ static void gm_tc_buf_write_kernel_elt_(
  GemmIn_t* vo,
  const uint32_t* vim,
  const uint32_t* vic,
  int vi_dim0,
  int num_way,
  bool is_sparse,
  bool is_right,
  bool is_duo,
  bool form_matX_on_accel,
  int step_2way,
  bool is_bitwise_3way_2step,
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
  //COMET_ASSERT(i01 == 0 || i01 == 1);
  const int vl = vlX2 / 2;

  const int flD2 = flD2_min + flD2_thisstep;

  // Output array interpreted as having GemmIn_t scalars has nfl rows.

  const uint32_t* const vim_col = vim + vl * (size_t)vi_dim0;

  // Pick up two consecutive field values:
  // first field seminibble0, second field seminibble1
  // Set to zero if outside of active range.

  const int nibblem = vl<nvlea ? (vim_col[flD2/8] >> (4*(flD2%8))) & 15 : 0;
  const int snm0 = nibblem & 3;
  const int snm1 = (nibblem>>2) & 3;

  const int nibblec = vl<nvlea ? (vic[flD2/8] >> (4*(flD2%8))) & 15 : 0;
  const int snc0 = nibblec & 3;
  const int snc1 = (nibblec>>2) & 3;

  // Count number of 0 (or 1) bits in respective seminibble.
  // Determine whether to skip (1,0) null indicator value.
  // NOTE: does not work for all cases.

  const bool is_left = ! is_right;
  const bool skip_10 = is_sparse || (num_way == 3 && is_left);

  // Possible counts, represented in target type.

  const GemmIn_t zero = TCBufTypes<GemmIn_t>::zero();
  const GemmIn_t one  = TCBufTypes<GemmIn_t>::one();
  const GemmIn_t two  = TCBufTypes<GemmIn_t>::two();
  const GemmIn_t four  = TCBufTypes<GemmIn_t>::four();

  // Possible seminibble bit patterns.

  const int _00 = 0;
  const int _01 = 1;
  const int _10 = 2;
  const int _11 = 3;

  // Unimplemented cases:
  //COMET_ASSERT( ! (is_duo && !is_sparse) );
  //COMET_ASSERT( ! (is_duo && !is_bitwise_3way_2step) );
  //COMET_ASSERT( ! (is_bitwise_3way_2step && !form_matX_on_accel) );

  const GemmIn_t out0 = 3 == num_way && is_left && is_duo ? (
                          snm0 == _10                              ? zero :
                          snc0 == _10                              ? zero :
                          (snm0&1) == step_2way && (snc0&1) == i01 ? one :
                                                                     zero
                        ) : //====================
                        is_duo ? (
                           snm0 == _10         ? zero :
                          (snm0 & 1) == i01    ? one :
                       /* (snm0 & 1) == 1-i01 */ zero
                        ) : //====================
                        3 == num_way && is_left &&
                        is_bitwise_3way_2step /* && is_ccc */ ? (
                          snm0 == _10 && is_sparse               ? zero :
                          snc0 == _10 && is_sparse               ? zero :
                          snm0 == _11*(1-step_2way)              ? zero :
                          snc0 == _11*(1-i01)                    ? zero :
                          is_po2(snm0) && is_po2(snc0)           ? one :
                          is_po2(snm0) && snc0 == _11*i01        ? two :
                          is_po2(snc0) && snm0 == _11*step_2way  ? two :
                        /* snm0*(3-snm0) + (snc0)*(3-snc0) == 0 */ four
                        ) : //====================
                        3 == num_way && is_left &&
                        form_matX_on_accel /* && is_ccc */ ? (
                          snm0 == _10 && is_sparse      ? zero :
                          snm0 == _00 && step_2way != 0 ? zero :
                          snm0 == _01 && step_2way != 1 ? zero :
                          snm0 == _10 && step_2way != 1 ? zero :
                          snm0 == _11 && step_2way != 2 ? zero :
                          snc0 == _11*i01               ? two :
                          snc0 == _11*(1-i01)           ? zero :
                          snc0 == _01                   ? one :
                                  is_sparse             ? zero :
                       /* snc0 == _10 */                  one
                        ) : //====================
                        /* is_ccc ... */ (
                          snm0 == _11*i01      ? two :
                          snm0 == _11*(1-i01)  ? zero :
                                  !skip_10     ? one :
                          snm0 == _01          ? one :
                       /* snm0 == _10 */         zero
                        );

  const GemmIn_t out1 = 3 == num_way && is_left && is_duo ? (
                          snm1 == _10                              ? zero :
                          snc1 == _10                              ? zero :
                          (snm1&1) == step_2way && (snc1&1) == i01 ? one :
                                                                     zero
                        ) : //====================
                        is_duo ? (
                           snm1 == _10         ? zero :
                          (snm1 & 1) == i01    ? one :
                       /* (snm1 & 1) == 1-i01 */ zero
                        ) : //====================
                        3 == num_way && is_left &&
                        is_bitwise_3way_2step /* && is_ccc */ ? (
                          snm1 == _10 && is_sparse               ? zero :
                          snc1 == _10 && is_sparse               ? zero :
                          snm1 == _11*(1-step_2way)              ? zero :
                          snc1 == _11*(1-i01)                    ? zero :
                          is_po2(snm1) && is_po2(snc1)           ? one :
                          is_po2(snm1) && snc1 == _11*i01        ? two :
                          is_po2(snc1) && snm1 == _11*step_2way  ? two :
                        /* snm1*(3-snm1) + (snc1)*(3-snc1) == 0 */ four
                        ) : //====================
                        3 == num_way && is_left &&
                        form_matX_on_accel /* && is_ccc */ ? (
                          snm1 == _10 && is_sparse      ? zero :
                          snm1 == _00 && step_2way != 0 ? zero :
                          snm1 == _01 && step_2way != 1 ? zero :
                          snm1 == _10 && step_2way != 1 ? zero :
                          snm1 == _11 && step_2way != 2 ? zero :
                          snc1 == _11*i01               ? two :
                          snc1 == _11*(1-i01)           ? zero :
                          snc1 == _01                   ? one :
                                  is_sparse             ? zero :
                       /* snc1 == _10 */                  one
                        ) : //====================
                        /* is_ccc ... */ (
                          snm1 == _11*i01      ? two :
                          snm1 == _11*(1-i01)  ? zero :
                                  !skip_10     ? one :
                          snm1 == _01          ? one :
                       /* snm1 == _10 */         zero
                        );

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
  const uint32_t* vim,
  const uint32_t* vic,
  int vi_dim0,
  int num_way,
  bool is_sparse,
  bool is_right,
  bool is_duo,
  bool form_matX_on_accel,
  int step_2way,
  bool is_bitwise_3way_2step,
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

  gm_tc_buf_write_kernel_elt_<GemmIn_t>(vo, vim, vic, vi_dim0,
    num_way, is_sparse, is_right, is_duo, form_matX_on_accel, step_2way,
    is_bitwise_3way_2step,
    nvlea, nvle, nvleD2, nvleX2, nfl, nflD2, nflD2_thisstep, flD2_min,
    vlX2, flD2_thisstep);
}

//-----------------------------------------------------------------------------
/// \brief Convert bitwise matrix to required format for GEMM.

template<int TC_METHOD>
static void gm_tc_buf_write_(
  bool is_right, int I_max, int I_max_dim, int nvl,
  int npvfl, int npvfl_thisstep, int pvfl_min,
  const uint32_t* vi1, const uint32_t* vi2, TCBufs& tc_bufs, int step_2way,
  GMEnv& env) {

  COMET_INSIST(vi1 && vi2);
  COMET_INSIST(I_max_dim >= 0 && I_max_dim <= nvl);
  COMET_INSIST(I_max >= 0 && I_max <= I_max_dim);
  COMET_INSIST(nvl >= 0 && npvfl >= 0);
  COMET_INSIST(tc_bufs.tc_buf_left && tc_bufs.tc_buf_right);
  COMET_INSIST(npvfl_thisstep >= 0 && npvfl_thisstep <= npvfl);
  COMET_INSIST(pvfl_min >= 0 && pvfl_min + npvfl_thisstep <= npvfl);

  // num_vector-related dimensions.

  const int nvle = is_right ? nvl : I_max_dim; // effective nvl dimension
  const int nvleD2 = nvle / 2;
  const int nvleX2 = nvle * 2;
  const int nvlea = is_right ? nvl : I_max; // num active nvle; others zeroed
  // NOTE: ignoring here the issue from decomp_mgr that
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

  // Arrays.

  typedef typename TCSelector<TC_METHOD>::GemmIn_t GemmIn_t;
  const int vi_dim0 = npvfl * 4; // 4 = sizeof(doublecomplex) / sizeof(int32)
  GemmIn_t* const tc_buf = is_right ? (GemmIn_t*)tc_bufs.tc_buf_right :
                                      (GemmIn_t*)tc_bufs.tc_buf_left;
  COMET_INSIST(nvleX2 * (size_t)(2*nflD2_thisstep) *
           sizeof(typename TCSelector<TC_METHOD>::GemmIn_t)
           <= tc_bufs.tc_buf_size &&
           "Subscriptrange error on tc buf.");

  const bool is_duo = env.metric_type() == MetricType::DUO;
  const bool form_matX_on_accel = env.form_matX_on_accel();
  const bool is_bitwise_3way_2step = env.is_bitwise_3way_2step();

  const uint32_t* unused_col = form_matX_on_accel ? NULL : vi1; // dummy
  const uint32_t* vim = form_matX_on_accel ? vi2 : vi1; // matrix
  const uint32_t* vic = form_matX_on_accel ? vi1 : unused_col; // column

  if (env.is_compute_method_gpu()) {

    // Kernel call.

#   ifdef COMET_USE_ACCEL

      const int threadblocksize = 256;
      const int blockdim_y = 32768;
      const int num_threadblocks_0 = utils::ceil(nvleX2, threadblocksize);
      const int num_threadblocks_1 = utils::min(nflD2_thisstep, blockdim_y);
      const int num_threadblocks_2 = utils::ceil(nflD2_thisstep, blockdim_y);

#     ifdef COMET_USE_HIP
        hipLaunchKernelGGL(
#     endif
        gm_tc_buf_write_kernel_<GemmIn_t>
#     ifdef COMET_USE_CUDA
        <<<
#     else
        ,
#     endif
        dim3(num_threadblocks_0, num_threadblocks_1, num_threadblocks_2),
        dim3(threadblocksize, 1, 1), 0, env.stream_compute()
#     ifdef COMET_USE_CUDA
        >>> (
#     else
        ,
#     endif
        tc_buf, vim, vic, vi_dim0, env.num_way(), env.sparse(), is_right,
        is_duo, form_matX_on_accel, step_2way, is_bitwise_3way_2step,
        nvlea, nvle, nvleD2, nvleX2, nfl, nflD2, nflD2_thisstep, flD2_min);

      System::accel_last_call_succeeded();

#   endif // COMET_USE_ACCEL

  } else { // (!env.is_compute_method_gpu())

    for (int flD2_thisstep=0; flD2_thisstep<nflD2_thisstep; ++flD2_thisstep) {
      for (int vlX2=0; vlX2<nvleX2; ++vlX2) {

        gm_tc_buf_write_kernel_elt_<GemmIn_t>(
          tc_buf, vim, vic, vi_dim0, env.num_way(), env.sparse(), is_right,
          is_duo, form_matX_on_accel, step_2way, is_bitwise_3way_2step,
          nvlea, nvle, nvleD2, nvleX2, nfl, nflD2, nflD2_thisstep, flD2_min,
          vlX2, flD2_thisstep);

      }
//printf("========================== %i\n", flD2_thisstep);
    }

  } // if (env.is_compute_method_gpu())
}

//=============================================================================
// BLAS GEMM FILE-LOCAL (STATIC) FUNCTIONS
//=============================================================================

//-----------------------------------------------------------------------------
/// \brief Call cublas to perform required GEMM.

template<int TC_METHOD>
static void gm_tc_solve_impl(bool is_first, int m, int n, int k,
  void* matC, TCBufs& tc_bufs, GMEnv& env) {
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
  // since I_max_dim % 4 == 0; see gm_gemm_divisibility_required()
  COMET_INSIST(m % 8 == 0 && "Failed divisibility condition for tc gemm.");
  // since nvl % 4 == 0; see gm_gemm_divisibility_required()
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
          printf("Error: %s\n", CUBLAS_STATUS_NOT_INITIALIZED == status ?
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
static void gm_tc_solve_(bool is_first, int nvll, int nvl, int npvfl_thisstep,
                         void* matC, TCBufs& tc_bufs, GMEnv& env) {
  COMET_INSIST(matC);
  COMET_INSIST(nvll >= 0 && nvl >= 0 && nvll <= nvl);
  COMET_INSIST(npvfl_thisstep >= 0);
  COMET_INSIST(env.tc_eff() != TC::NO);

  const int nfl_thisstep = npvfl_thisstep * 64;

  const int m = 2 * nvll; // metrics array dim
  const int n = 2 * nvl; // metrics array dim
  const int k = nfl_thisstep; // vectors array (as GemmIn_t) dim

  gm_tc_solve_impl<TC_METHOD>(is_first, m, n, k, matC, tc_bufs, env);
}

//=============================================================================
// REPAIR METRICS FILE-LOCAL (STATIC) FUNCTIONS
//=============================================================================

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
  GMEnv& env) {

  COMET_INSIST(vo);
  COMET_INSIST(nvll >= 0);
  COMET_INSIST(nvl >= 0);
  COMET_INSIST(nvll <= nvl);

  // always true, because of gm_gemm_divisibility_required()
  COMET_INSIST(nvll % 2 == 0 && "Failed divisibility condition for tc gemm.");
  const int nvllD2 = nvll / 2;

  typedef typename TCSelector<TC_METHOD>::GemmOut_t GemmOut_t;

  if (env.is_compute_method_gpu()) {

    // Kernel call.

#ifdef COMET_USE_ACCEL

    const int threadblocksize = 256;
    const int vll2_threadblocks = utils::ceil(nvllD2, threadblocksize);

#  ifdef COMET_USE_HIP
    hipLaunchKernelGGL(
#  endif
    gm_tc_repair_metrics_kernel_<GemmOut_t>
#  ifdef COMET_USE_CUDA
        <<<
#  else
        ,
#  endif
        dim3(vll2_threadblocks, nvl, 1),
        dim3(threadblocksize, 1, 1),
        0,
        env.stream_compute()
#  ifdef COMET_USE_CUDA
        >>> (
#  else
        ,
#  endif
        nvl, nvll, nvllD2, vo);

    System::accel_last_call_succeeded();

#endif // COMET_USE_ACCEL

  } else { // (!env.is_compute_method_gpu())

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

//=============================================================================
// TOP-LEVEL FILE-LOCAL (STATIC) FUNCTIONS
//=============================================================================

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
  const void* matA1, const void* matA2, const void* matB, void* matC, int lddc,
  TCBufs& tc_bufs, int step_2way, GMEnv& env) {

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
  COMET_INSIST((size_t)nvll == gm_gemm_size_required(nvll, env));

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
    gm_tc_buf_write_<TC_METHOD>(left_matrix, I_max, I_max_dim, nvl, npvfl,
      npvfl_thisstep, pvfl_min, (uint32_t*)matA1, (uint32_t*)matA2, tc_bufs,
      step_2way, env);
    gm_tc_buf_write_<TC_METHOD>(right_matrix, I_max, I_max_dim, nvl, npvfl,
      npvfl_thisstep, pvfl_min, (uint32_t*)matB, (uint32_t*)matB, tc_bufs,
      step_2way, env);

    // Perform the GEMM for this pair of block rows; accumulate.
    const bool is_first = 0 == pvfl_min;
    gm_tc_solve_<TC_METHOD>(is_first, nvll, nvl, npvfl_thisstep,
      matC, tc_bufs, env);
  }

  // Revise the results of the GEMMs to be in the needed double complex format.
  gm_tc_repair_metrics_<TC_METHOD>(nvll, nvl, matC, tc_bufs, env);
}

//=============================================================================
// EXTERNALLY VISIBLE FUNCTIONS: GENERAL
//=============================================================================

//-----------------------------------------------------------------------------
/// \brief Use a standard GEMM to compute CoMet metrics bitwise result.

void gm_tc_gemm_start(
  int m, int n, int k,
  const void* matA1, int ldda1, const void* matA2, int ldda2,
  const void* matB, int lddb, void* matC, int lddc,
  int step_2way, TCBufs& tc_bufs, GMEnv& env) {
  COMET_INSIST(matA1 && matA2 && matB && matC);
  COMET_INSIST(m >= 0 && n >= 0 && k >= 0);
  COMET_INSIST(ldda1 >= 0 && ldda2 >= 0 && lddb >= 0 && lddc >= 0);
  COMET_INSIST(k <= ldda1 && k <= ldda2 && k <= lddb && m <= lddc);
  COMET_INSIST(env.tc_eff() != TC::NO);
  COMET_INSIST(env.is_metric_type_bitwise());

  COMET_INSIST(tc_bufs.tc_buf_left);

  COMET_INSIST(ldda1 == k && ldda2 == k && lddb == k); // always true here

  // Select required template function instance.

  switch (env.tc_eff()) {
    // --------------
    case TC::INT8: {
      gm_tc_gemm_start_impl_<TC::INT8>(
        m, n, k, matA1, matA2, matB, matC, lddc, tc_bufs, step_2way, env);
    } break;
    // --------------
    case TC::FP16: {
      gm_tc_gemm_start_impl_<TC::FP16>(
        m, n, k, matA1, matA2, matB, matC, lddc, tc_bufs, step_2way, env);
    } break;
    // --------------
    case TC::FP32: {
      gm_tc_gemm_start_impl_<TC::FP32>(
        m, n, k, matA1, matA2, matB, matC, lddc, tc_bufs, step_2way, env);
    } break;
    // --------------
    default:
      COMET_INSIST(false && "Invalid tc type.");
  } // switch
}

//-----------------------------------------------------------------------------
/// \brief Divisibility requirement for GEMM.

size_t gm_gemm_divisibility_required(const GMEnv& env) {

  const bool need_divisible_by_4 = env.tc_eff() != TC::NO;

  return need_divisible_by_4 ? 4 : 1;
}

//-----------------------------------------------------------------------------
/// \brief Size requirement for GEMM.

size_t gm_gemm_size_required(size_t size_requested, const GMEnv& env) {

  const size_t factor = gm_gemm_divisibility_required(env);

  return utils::ceil(size_requested, factor)*factor;
}

//=============================================================================
// EXTERNALLY VISIBLE FUNCTIONS: BUFFER MANAGEMENT
//=============================================================================

//-----------------------------------------------------------------------------
/// \brief Initialize TCBufs object by allocating memory etc.

void gm_tc_bufs_malloc(int num_vector_local,
                       int num_field_local,
                       int num_packedval_field_local,
                       TCBufs& tc_bufs,
                       GMEnv& env) {
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
     env.tc_eff() == TC::INT8 ?
       sizeof(typename TCSelector<TC::INT8>::GemmIn_t) :
     env.tc_eff() == TC::FP16 ?
       sizeof(typename TCSelector<TC::FP16>::GemmIn_t) :
     env.tc_eff() == TC::FP32 ?
       sizeof(typename TCSelector<TC::FP32>::GemmIn_t) :
     0;
  COMET_INSIST(TC::is_valid(env.tc_eff())); // this code must be updated if new method

  const size_t nvlX2 = nvl * 2;

  tc_bufs.tc_buf_size = nvlX2 * (npvfl_thisstep_max * 64) * sizeof_gemm_in_t;
  tc_bufs.tc_buf_size = tc_bufs.tc_buf_size ? tc_bufs.tc_buf_size : 1;

  if (env.is_compute_method_gpu()) {

    // Allocate buffers.

#if defined COMET_USE_CUDA
    cudaMalloc(&tc_bufs.tc_buf_left, tc_bufs.tc_buf_size);
#elif defined COMET_USE_HIP
    hipMalloc(&tc_bufs.tc_buf_left, tc_bufs.tc_buf_size);
#endif
    System::accel_last_call_succeeded();
    env.gpu_mem_local_inc(tc_bufs.tc_buf_size);

#if defined COMET_USE_CUDA
    cudaMalloc(&tc_bufs.tc_buf_right, tc_bufs.tc_buf_size);
#elif defined COMET_USE_HIP
    hipMalloc(&tc_bufs.tc_buf_right, tc_bufs.tc_buf_size);
#endif
    System::accel_last_call_succeeded();
    env.gpu_mem_local_inc(tc_bufs.tc_buf_size);

    // Set up accel blas handle.

#if defined COMET_USE_CUDA
    cublasStatus_t status = cublasCreate(&tc_bufs.accelblas_handle);
    COMET_INSIST(status == CUBLAS_STATUS_SUCCESS && "Error in cublasCreate.");

    status = cublasSetStream(tc_bufs.accelblas_handle, env.stream_compute());
    COMET_INSIST(status == CUBLAS_STATUS_SUCCESS && "Error in cublasSetStream.");

    status = cublasSetMathMode(tc_bufs.accelblas_handle, CUBLAS_TENSOR_OP_MATH);
    COMET_INSIST(status == CUBLAS_STATUS_SUCCESS && "Error in cublasSetMathMode.");
#elif defined COMET_USE_HIP
    //rocbas_status_ status = rocblas_create_handle(&tc_bufs.accelblas_handle);
    int status = rocblas_create_handle(&tc_bufs.accelblas_handle);
    COMET_INSIST(status == rocblas_status_success &&
             "Error in rocblas_create_handle.");

    status = rocblas_set_stream(tc_bufs.accelblas_handle, env.stream_compute());
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

void gm_tc_bufs_free(TCBufs& tc_bufs, GMEnv& env) {
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
    System::accel_last_call_succeeded();
    tc_bufs.tc_buf_left = NULL;
    env.gpu_mem_local_dec(tc_bufs.tc_buf_size);

#if defined COMET_USE_CUDA
    cudaFree(tc_bufs.tc_buf_right);
#elif defined COMET_USE_HIP
    hipFree(tc_bufs.tc_buf_right);
#endif
    System::accel_last_call_succeeded();
    tc_bufs.tc_buf_right = NULL;
    env.gpu_mem_local_dec(tc_bufs.tc_buf_size);

    // Free accel blas handle.

#if defined COMET_USE_CUDA
    cublasStatus_t status = cublasDestroy(tc_bufs.accelblas_handle);
    COMET_INSIST(status == CUBLAS_STATUS_SUCCESS && "Error in cublasDestroy.");
#elif defined COMET_USE_HIP
    //rocblas_status status = rocblas_destroy_handle(tc_bufs.accelblas_handle);
    int status = rocblas_destroy_handle(tc_bufs.accelblas_handle);
    COMET_INSIST(status == rocblas_status_success &&
             "Error in rocblas_destroy_handle.");
#endif

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
