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

#include "cstdlib"
#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include <mma.h>

#include "tc_solve_cutlass_general.i.hh"
#include "tc_solve_comet_xor.i.hh"
#include "tc_solve_comet_xor_int.i.hh"
#include "tc_solve_comet_mult.i.hh"
#include "tc_solve_comet_mult_int.i.hh"

#ifdef COMET_USE_CUTLASS
#include "cutlass/gemm/device/gemm.h"
#endif

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

struct TCSubmethod {
  enum {SIMPLE = 0,
        _256_128 = 1,
        _128_256 = 2,
        _128_128 = 3,
        _128_64 = 4,
        _64_128 = 5,
        _64_64 = 6,
        _64_64_WMMA = 7,
        _128_256_1024 = 8,
        _INT4 = 9,
        _INT8 = 10
  };
};

struct TCGemmOpXorPopc {
  template<typename GemmIn_t>
  __host__ __device__
  static GemmIn_t op(GemmIn_t a, GemmIn_t b) {return a ^ b;}
# ifdef COMET_USE_CUTLASS
    typedef cutlass::arch::OpXorPopc Value; 
#endif
};

struct TCGemmOpMultiplyAdd {
  template<typename GemmIn_t>
  __host__ __device__
  static GemmIn_t op(GemmIn_t a, GemmIn_t b) {return a & b;}
# ifdef COMET_USE_CUTLASS
#   if COMET_COMPUTE_CAPABILITY != 750
      typedef cutlass::arch::OpMultiplyAdd Value; 
#   else
      // UNUSED - fix for template compilation problem on Turing.
      typedef cutlass::arch::OpXorPopc Value; 
#   endif
# endif
};

struct TCGemmOpMultiplyAddSaturate {
  template<typename GemmIn_t>
  __host__ __device__
  static GemmIn_t op(GemmIn_t a, GemmIn_t b) {return a & b;}
# ifdef COMET_USE_CUTLASS
#   if COMET_COMPUTE_CAPABILITY != 750
      typedef cutlass::arch::OpMultiplyAddSaturate Value; 
#   else
      // Fix for template compilation problem on Turing.
      typedef cutlass::arch::OpXorPopc Value; 
#   endif
# endif
};

// Mixin class.
struct CutlassOpClassTensorOp {
# ifdef COMET_USE_CUTLASS
  typedef typename cutlass::arch::OpClassTensorOp OpClass_t;
# endif
};

// Mixin class.
struct CutlassOpClassWmmaTensorOp {
# ifdef COMET_USE_CUTLASS
  typedef typename cutlass::arch::OpClassWmmaTensorOp OpClass_t;
# endif
};

// Mixin class.
struct CutlassArch {
# ifdef COMET_USE_CUTLASS
#   if COMET_COMPUTE_CAPABILITY != 750
      typedef typename cutlass::arch::Sm80 CutlassArch_t;
      enum {STAGES = 3};
#   else
      typedef typename cutlass::arch::Sm75 CutlassArch_t;
      enum {STAGES = 2};
#   endif
# endif
};

// Mixin class.
template<int _0, int _1, int _2>
struct CutlassThreadblockShape {
  enum {ThreadBlockShape0 = _0,
        ThreadBlockShape1 = _1,
        ThreadBlockShape2 = _2
  };
};

// Mixin class.
template<int _0, int _1, int _2>
struct CutlassWarpShape {
  enum {WarpShape0 = _0,
        WarpShape1 = _1,
        WarpShape2 = _2
  };
};

// Mixin class.
template<int _0, int _1, int _2>
struct CutlassInstructionShape {
  enum {InstructionShape0 = _0,
        InstructionShape1 = _1,
        InstructionShape2 = _2
  };
};

//-----------------------------------------------------------------------------

template<int TC_SUBMETHOD> struct CutlassSettings;

template<> struct CutlassSettings<TCSubmethod::_256_128>
  : public CutlassArch,
           CutlassOpClassTensorOp,
           CutlassThreadblockShape<256,128,512>,
           CutlassWarpShape<64,64,512>,
           CutlassInstructionShape<16,8,256> {};

template<> struct CutlassSettings<TCSubmethod::_128_256>
  : public CutlassArch,
           CutlassOpClassTensorOp,
           CutlassThreadblockShape<128,256,512>,
           CutlassWarpShape<64,64,512>,
           CutlassInstructionShape<16,8,256> {};

template<> struct CutlassSettings<TCSubmethod::_128_256_1024>
  : public CutlassArch,
           CutlassOpClassTensorOp,
           CutlassThreadblockShape<128,256,1024>,
           CutlassWarpShape<64,64,1024>,
           CutlassInstructionShape<16,8,256> {};

template<> struct CutlassSettings<TCSubmethod::_128_128>
  : public CutlassArch,
           CutlassOpClassTensorOp,
           CutlassThreadblockShape<128,128,512>,
           CutlassWarpShape<64,64,512>,
           CutlassInstructionShape<16,8,256> {};

template<> struct CutlassSettings<TCSubmethod::_128_64>
  : public CutlassArch,
           CutlassOpClassTensorOp,
           CutlassThreadblockShape<128,64,512>,
           CutlassWarpShape<64,64,512>,
           CutlassInstructionShape<16,8,256> {};

template<> struct CutlassSettings<TCSubmethod::_64_128>
  : public CutlassArch,
           CutlassOpClassTensorOp,
           CutlassThreadblockShape<64,128,512>,
           CutlassWarpShape<64,64,512>,
           CutlassInstructionShape<16,8,256> {};

template<> struct CutlassSettings<TCSubmethod::_64_64>
  : public CutlassArch,
           CutlassOpClassTensorOp,
           CutlassThreadblockShape<64,64,512>,
           CutlassWarpShape<64,64,512>,
           CutlassInstructionShape<16,8,256> {};

template<> struct CutlassSettings<TCSubmethod::_64_64_WMMA>
  : public CutlassArch,
           CutlassOpClassWmmaTensorOp,
           CutlassThreadblockShape<64,64,512>,
           CutlassWarpShape<64,64,512>,
           CutlassInstructionShape<16,8,256> {};

template<> struct CutlassSettings<TCSubmethod::_INT4>
  : public CutlassArch,
           CutlassOpClassTensorOp,
           CutlassThreadblockShape<128,256,256>,
           CutlassWarpShape<64,64,256>,
           CutlassInstructionShape<16,8,64> {};

template<> struct CutlassSettings<TCSubmethod::_INT8>
  : public CutlassArch,
           CutlassOpClassTensorOp,
           CutlassThreadblockShape<128,256,64>,
           CutlassWarpShape<64,64,64>,
           CutlassInstructionShape<16,8,32> {};

           //CutlassThreadblockShape<128,256,256>, CutlassWarpShape<64,64,256>, CutlassInstructionShape<16,8,64> {}; // 959 TOps
           //CutlassThreadblockShape<256,128,256>, CutlassWarpShape<64,64,256>, CutlassInstructionShape<16,8,64> {}; // 575 TOps
           //CutlassThreadblockShape<128,128,256>, CutlassWarpShape<64,64,256>, CutlassInstructionShape<16,8,64> {}; // 600 TOps
           //CutlassThreadblockShape<256,64,256>, CutlassWarpShape<64,64,256>, CutlassInstructionShape<16,8,64> {}; // 316 TOps
           //CutlassThreadblockShape<64,256,256>, CutlassWarpShape<64,64,256>, CutlassInstructionShape<16,8,64> {}; // 849 TOps
           //CutlassThreadblockShape<64,128,256>, CutlassWarpShape<32,64,256>, CutlassInstructionShape<16,8,64> {}; // 549 TOps
           //CutlassThreadblockShape<128,64,256>, CutlassWarpShape<64,32,256>, CutlassInstructionShape<16,8,64> {}; // 309 TOps
           //CutlassThreadblockShape<64,64,256>, CutlassWarpShape<32,32,256>, CutlassInstructionShape<16,8,64> {}; // 307 TOps

           //CutlassThreadblockShape<128,256,128>, CutlassWarpShape<64,64,128>, CutlassInstructionShape<16,8,64> {}; // 560 TOps
           //CutlassThreadblockShape<256,128,128>, CutlassWarpShape<64,64,128>, CutlassInstructionShape<16,8,64> {}; // 362 TOps
           //CutlassThreadblockShape<128,128,128>, CutlassWarpShape<64,64,128>, CutlassInstructionShape<16,8,64> {}; // 344 TOps
           //CutlassThreadblockShape<256,64,128>, CutlassWarpShape<64,64,128>, CutlassInstructionShape<16,8,64> {}; // 198 TOps
           //CutlassThreadblockShape<64,256,128>, CutlassWarpShape<64,64,128>, CutlassInstructionShape<16,8,64> {}; // 460 TOps
           //CutlassThreadblockShape<64,128,128>, CutlassWarpShape<32,64,128>, CutlassInstructionShape<16,8,64> {}; // 319 TOps
           //CutlassThreadblockShape<128,64,128>, CutlassWarpShape<64,32,128>, CutlassInstructionShape<16,8,64> {}; // 191 TOps
           //CutlassThreadblockShape<64,64,128>, CutlassWarpShape<32,32,128>, CutlassInstructionShape<16,8,64> {}; // 189 TOps

//-----------------------------------------------------------------------------

template<int TC_METHOD>
struct CutlassElementInputType {};

template<>
struct CutlassElementInputType<TC::B1> {
# ifdef COMET_USE_CUTLASS
  typedef cutlass::uint1b_t Value;
# endif
};

template<>
struct CutlassElementInputType<TC::INT4> {
# ifdef COMET_USE_CUTLASS
  typedef cutlass::int4b_t Value;
# endif
};

template<>
struct CutlassElementInputType<TC::INT8> {
# ifdef COMET_USE_CUTLASS
  typedef uint8_t Value;
# endif
};

//-----------------------------------------------------------------------------

template<int TC_METHOD, int TC_SUBMETHOD, typename TCGemmOp>
void tc_solve_impl_cutlass(
  bool is_first, int m, int n, int k,
  uint8_t const *A, int lda,
  uint8_t const *B, int ldb,
  int32_t *C, int ldc,
  AccelStream_t accel_stream) {

# ifdef COMET_USE_CUTLASS

  // Extra checks, may be too strict - dimensions, alignment.

  COMET_INSIST(m % 256 == 0);
  COMET_INSIST(n % 256 == 0);
  COMET_INSIST(k % 256 == 0);

  COMET_INSIST(((size_t)A) % 256 == 0);
  COMET_INSIST(((size_t)B) % 256 == 0);
  COMET_INSIST(((size_t)C) % 256 == 0);

  using ElementInput = typename CutlassElementInputType<TC_METHOD>::Value;
  using ElementCompute = int32_t;
  using ElementAccumulator = int32_t;
  using ElementOutput = int32_t;

  // see https://github.com/NVIDIA/cutlass/blob/master/include/cutlass/gemm/device/gemm.h
  // https://github.com/NVIDIA/cutlass/blob/master/include/cutlass/gemm/device/default_gemm_configuration.h
  // https://github.com/NVIDIA/cutlass/blob/master/test/unit/gemm/device/gemm_b1t_b1n_s32t_tensor_op_s32_sm80.cu
  // https://github.com/NVIDIA/cutlass/blob/master/test/unit/gemm/device/gemm_s4t_s4n_s32t_tensor_op_s32_sm80.cu

  using Gemm = cutlass::gemm::device::Gemm<
      ElementInput, cutlass::layout::RowMajor,
      ElementInput, cutlass::layout::ColumnMajor,
      ElementOutput, cutlass::layout::RowMajor,
      ElementAccumulator,
      typename CutlassSettings<TC_SUBMETHOD>::OpClass_t,
      typename CutlassSettings<TC_SUBMETHOD>::CutlassArch_t,
      cutlass::gemm::GemmShape< // ThreadblockShape_
        CutlassSettings<TC_SUBMETHOD>::ThreadBlockShape0,
        CutlassSettings<TC_SUBMETHOD>::ThreadBlockShape1,
        CutlassSettings<TC_SUBMETHOD>::ThreadBlockShape2>,
      cutlass::gemm::GemmShape< // WarpShape_
        CutlassSettings<TC_SUBMETHOD>::WarpShape0,
        CutlassSettings<TC_SUBMETHOD>::WarpShape1,
        CutlassSettings<TC_SUBMETHOD>::WarpShape2>,
      cutlass::gemm::GemmShape< // InstructionShape_
        CutlassSettings<TC_SUBMETHOD>::InstructionShape0,
        CutlassSettings<TC_SUBMETHOD>::InstructionShape1,
        CutlassSettings<TC_SUBMETHOD>::InstructionShape2>,
      cutlass::epilogue::thread::LinearCombinationClamp<
          ElementOutput,
          128 / cutlass::sizeof_bits<ElementOutput>::value,
          ElementAccumulator,
          ElementCompute>,
      cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
      CutlassSettings<TC_SUBMETHOD>::STAGES, // Stages
      128 / cutlass::sizeof_bits<ElementInput>::value, // AlignmentA
      128 / cutlass::sizeof_bits<ElementInput>::value, // AlignmentB
      false, // SplitKSerial
      typename TCGemmOp::Value>;

  Gemm gemm_operator;

  const ElementCompute alpha = 1;
  const ElementCompute beta = is_first ? 0 : 1;

  typename Gemm::Arguments args(
    {m, n, k}, // Gemm Problem dimensions
    {(ElementInput const*)A, lda}, // Tensor-ref for source matrix A
    {(ElementInput const*)B, ldb}, // Tensor-ref for source matrix B
    {(ElementOutput *)C, ldc}, // Tensor-ref for source matrix C
    {(ElementOutput *)C, ldc}, // Tensor-ref for destination matrix C
    {alpha,beta}); // Scalars used in the Epilogue


  // Perform GEMM.
  const cutlass::Status status = gemm_operator(args, nullptr, accel_stream);
  COMET_INSIST(status == cutlass::Status::kSuccess);
  System::accel_last_call_succeeded();

# else

  COMET_INSIST(false && "Failure to call GEMM function.");

# endif
}

//=============================================================================
/// \brief 1-bit gemm kernel (mockup version, not high performance).

template<typename GemmIn_t, typename GemmOut_t, typename TCGemmOp>
__global__ void tc_solve_b1_gemm_mockup_kernel(
  size_t m, size_t n, size_t k,
  GemmIn_t* a, GemmIn_t* b, bool beta, GemmOut_t* c) {
  //COMET_INSIST(a && b && c);

  // k axis serialized; m, n axes threaded.

  const size_t ind_m = threadIdx_x_() + blockIdx_x_() * blockDim_x_();
  const size_t ind_n = blockIdx_y_();

  if (ind_m >= m || ind_n >= n)
    return;

  for (size_t ind_k = 0; ind_k < k; ++ind_k) {

    const GemmIn_t aik = a[ind_k + k*ind_m];
    const GemmIn_t bjk = b[ind_k + k*ind_n];

    GemmOut_t& cij = c[ind_m + m*ind_n];

    // Use "xor" or "bitwise and"; count 1-bits with popcount.
    const GemmOut_t v = utils::popc<GemmIn_t>(TCGemmOp::op(aik, bjk));
    //const GemmOut_t v = utils::popc<GemmIn_t>(aik & bjk);
    if (ind_m < 1024 * 1024 * 1024) // WORKAROUND for undiagnosed error.
      cij = beta || ind_k ? cij + v : v;

  } // for ind_k
}

//=============================================================================
/// \brief INT4 gemm kernel (mockup version, not high performance).

template<typename GemmIn_t, typename GemmOut_t>
__global__ void tc_solve_int4_gemm_mockup_kernel(
  size_t m, size_t n, size_t k,
  GemmIn_t* a, GemmIn_t* b, bool beta, GemmOut_t* c) {
  //COMET_INSIST(a && b && c);

  // k axis serialized; m, n axes threaded.
  
  const size_t ind_m = threadIdx_x_() + blockIdx_x_() * blockDim_x_();
  const size_t ind_n = blockIdx_y_();

  if (ind_m >= m || ind_n >= n)
    return;

  for (size_t ind_k = 0; ind_k < k; ++ind_k) {

    const GemmIn_t aik = a[ind_k + k*ind_m];
    const GemmIn_t bjk = b[ind_k + k*ind_n];

    GemmOut_t& cij = c[ind_m + m*ind_n];

    for (int nibble = 0; nibble < sizeof(GemmIn_t) * 2 ; ++nibble) {
      const GemmOut_t v = ((aik >> (4*nibble)) & (GemmIn_t)(16-1)) *
                          ((bjk >> (4*nibble)) & (GemmIn_t)(16-1));
      //if (ind_m < 1024 * 1024 * 1024) // WORKAROUND for undiagnosed error.
      cij = beta || ind_k || nibble ? cij + v : v;
    } // for i

  } // for ind_k
}

////-----------------------------------------------------------------------------
///// \brief 1-bit xor gemm.
//
//static bool tc_solve_use_mockup(const CEnv& env) {
//  //return true;
//  return ((env.tc_eff() == TC::B1 || env.tc_eff() == TC::INT4) &&
//    ! (BuildHas::CUDA && BuildHas::CUTLASS &&
//       System::compute_capability() > 700));
//}

//-----------------------------------------------------------------------------
/// \brief Simple WMMA tensor core 1-bit xor gemm

__global__
void b1_xor_gemm_gpu_tc_simple(size_t m, size_t n, size_t k, uint8_t* a,
                               uint8_t* b, bool beta, int32_t* c) {

  using namespace nvcuda;

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

  wmma::fragment<wmma::accumulator, WMMA1B_M, WMMA1B_N, WMMA1B_K, int> acc_frag;
  wmma::fragment<wmma::accumulator, WMMA1B_M, WMMA1B_N, WMMA1B_K, int> c_frag;
  wmma::fill_fragment(acc_frag, 0);

  // Loop over all sub-matrices of A and B to compute block sub-matrix
  int nblocks=((k*NBITS)+WMMA1B_K-1)/WMMA1B_K;
  for(int block=0; block<nblocks; block++) {

    // Load the inputs
    wmma::load_matrix_sync(a_frag, a+aBegin+block*aStep, k*NBITS);
    wmma::load_matrix_sync(b_frag, b+bBegin+block*bStep, k*NBITS);

    // Perform the matrix multiplication
    wmma::bmma_sync(acc_frag, a_frag, b_frag, acc_frag);
  }

  int cBegin = n*WMMA1B_M*bx + WMMA1B_N*by;
  if(beta) {
    // Load C fragment
    wmma::load_matrix_sync(c_frag, c+cBegin, n, wmma::mem_row_major);

    // Add acc_frag to c_frag
    for(int i=0; i<c_frag.num_elements; i++) {
      c_frag.x[i] += acc_frag.x[i];
    }

    // Store the output
    wmma::store_matrix_sync(c+cBegin, c_frag, n, wmma::mem_row_major);
  }
  else {
    // Store the output
    wmma::store_matrix_sync(c+cBegin, acc_frag, n, wmma::mem_row_major);
  }

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
  using namespace nvcuda;

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

__global__
void b1_print_matrix(int m, int n, int32_t* c) {

  // Block and thread indices
  int tx = threadIdx.x, ty = threadIdx.y;
  int bx = blockIdx.x, by = blockIdx.y;

  int gridx = bx*BLOCK_SIZE + tx;
  int gridy = by*BLOCK_SIZE + ty;

  if(gridy>=m || gridx>=n) return;

  int cBegin = by*n*2*BLOCK_SIZE*2;
  int rind1 = ty * n*4;
  int rind2 = ty * n*4 + n*2;
  int cind = ((tx + bx*BLOCK_SIZE) % (n/2))*4 + ((tx + bx*BLOCK_SIZE) / (n/2))*2;
  int cInd1 = cBegin + rind1 + cind;
  int cInd2 = cBegin + rind2 + cind;

  int32_t c0, c1, c2, c3;
  c0 = c[cInd1]; c1 = c[cInd1+1]; c2 = c[cInd2]; c3 = c[cInd2+1];

  //int px0=0, px1=1, py0=0, py1=1;
  //int px0=48, px1=49, py0=48, py1=49;
  //int px0=47, px1=49, py0=47, py1=49;
  int px0=95, px1=96, py0=95, py1=96;

  if(bx>=px0 && bx<px1 && by>=py0 && by<py1)
    printf("b=%d,%d t=%d,%d cb=%d ci1=%d=%d,%d ci2=%d=%d,%d c0123=%d,%d,%d,%d\n",
           bx,by,tx,ty,cBegin,cInd1,rind1,cind,cInd2,rind2,cind,c0,c1,c2,c3);
}

//-----------------------------------------------------------------------------

template<int TC_METHOD>
static void tc_solve_impl_subbyte(bool is_first, int m, int n, int k,
  void* matC, TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST(matC);
  COMET_INSIST(m >= 0 && n >= 0 && k >= 0);

  if(env.print_details()) printf("In tc_solve_impl_b1 num_kernel=%d use_mockup=%d\n",
    env.num_kernel(),env.is_using_cutlass_mockup());
  if (env.is_using_cutlass_mockup()) {
  //if (tc_solve_use_mockup(env)) {

    //COMET_INSIST(TCTraits<TC_METHOD>::IS_B_FIELD_MAJOR);

    enum {NUM_FL_PER_PFL = 64};
    COMET_INSIST(k % NUM_FL_PER_PFL == 0 &&
                 "Failed divisibility condition for tc gemm.");

    const bool beta = is_first ? 0 : 1;
    no_unused_variable_warning(beta);

    // 8 == number of uint8_t values used to store each chunk of
    // NUM_FL_PER_PFL fields in the tc buf.
    enum {BITS_PER_BYTE = 8};
    enum {BYTES_PER_PFL_FIELDS =
      (NUM_FL_PER_PFL*TCTraits<TC_METHOD>::NUM_BITS_PER_FIELD) / BITS_PER_BYTE};

    typedef typename TCTraits<TC_METHOD>::GemmIn_t GemmIn_t;
    typedef typename TCTraits<TC_METHOD>::GemmOut_t GemmOut_t;
    no_unused_type_warning<GemmOut_t>();

    const int bytes_per_gi = sizeof(GemmIn_t);
    const size_t k_eff = (k / NUM_FL_PER_PFL) *
                         (BYTES_PER_PFL_FIELDS / bytes_per_gi);
    no_unused_variable_warning(k_eff);

    const int threadblocksize = 256;
    COMET_INSIST((threadblocksize <= 256 || ! BuildHas::HIP) &&
                 "Current HIP limitation.");
    const int num_threadblocks_0 = utils::ceil(m, threadblocksize);
    const int num_threadblocks_1 = n;

    if (TC_METHOD == TC::INT4) {

      COMET_LAUNCH_KERNEL((tc_solve_int4_gemm_mockup_kernel
                           <GemmIn_t, GemmOut_t>),
        dim3(num_threadblocks_0, num_threadblocks_1, 1),
        dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
        m, n, k_eff, (GemmIn_t*)tc_bufs.tc_buf_left,
        (GemmIn_t*)tc_bufs.tc_buf_right, beta, (GemmOut_t*)matC);

    } else if (env.is_using_xor()) { // && TC_METHOD == TC::B1

      COMET_LAUNCH_KERNEL((tc_solve_b1_gemm_mockup_kernel
                           <GemmIn_t, GemmOut_t, TCGemmOpXorPopc>),
        dim3(num_threadblocks_0, num_threadblocks_1, 1),
        dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
        m, n, k_eff, (GemmIn_t*)tc_bufs.tc_buf_left,
        (GemmIn_t*)tc_bufs.tc_buf_right, beta, (GemmOut_t*)matC);

    } else { // !env.is_using_xor() // && TC_METHOD == TC::B1

      if(env.print_details()) printf("Launching b1_xor_gemm_gpu m=%d n=%d k_eff=%zu k=%d beta=%d "
          "bytes_per_gi=%d NUM_FL_PER_PVFL=%d gridDim=%d,%d threadDim=%d,1\n",
          m,n,k_eff,k,(int)beta,bytes_per_gi,NUM_FL_PER_PFL,
          num_threadblocks_0,num_threadblocks_1,threadblocksize);
      COMET_LAUNCH_KERNEL((tc_solve_b1_gemm_mockup_kernel
                           <GemmIn_t, GemmOut_t, TCGemmOpMultiplyAdd>),
        dim3(num_threadblocks_0, num_threadblocks_1, 1),
        dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
        m, n, k_eff, (GemmIn_t*)tc_bufs.tc_buf_left,
        (GemmIn_t*)tc_bufs.tc_buf_right, beta, (GemmOut_t*)matC);

    } // if (TC_METHOD == TC::INT4)

    System::accel_last_call_succeeded();

  } else if(env.num_kernel()==10) {

    COMET_INSIST(env.is_using_cutlass());

#   if COMET_COMPUTE_CAPABILITY != 750
      enum {TC_SUBMETHOD_B1 = TCSubmethod::_128_256_1024};
#   else
      enum {TC_SUBMETHOD_B1 = TCSubmethod::_128_128};
#   endif

    if (TC_METHOD == TC::INT8) {

      tc_solve_impl_cutlass<TC::INT8, TCSubmethod::_INT8, TCGemmOpMultiplyAddSaturate>(
        is_first, n, m, k, // NOTE: switching order of A, B.
        (uint8_t*)tc_bufs.tc_buf_right, k,
        (uint8_t*)tc_bufs.tc_buf_left, k,
        (int32_t*)matC, m,
        env.stream_compute());

    } else if (TC_METHOD == TC::INT4) {
      if(env.print_details()) printf("Using INT4 Mult-Add Cutlass GEMM\n");
      tc_solve_impl_cutlass<TC::INT4, TCSubmethod::_INT4, TCGemmOpMultiplyAddSaturate>(
        is_first, n, m, k, // NOTE: switching order of A, B.
        (uint8_t*)tc_bufs.tc_buf_right, k,
        (uint8_t*)tc_bufs.tc_buf_left, k,
        (int32_t*)matC, m,
        env.stream_compute());

    } else if (env.is_using_xor()) { // && TC_METHOD == TC::B1

      if(env.print_details()) printf("Using Xor Cutlass GEMM\n");
      tc_solve_impl_cutlass<TC::B1, TC_SUBMETHOD_B1, TCGemmOpXorPopc>(
        is_first, n, m, k, // NOTE: switching order of A, B.
        (uint8_t*)tc_bufs.tc_buf_right, k,
        (uint8_t*)tc_bufs.tc_buf_left, k,
        (int32_t*)matC, m,
        env.stream_compute());

    } else { // !env.is_using_xor() // && TC_METHOD == TC::B1

      if(env.print_details()) printf("Using Mult-Add Cutlass GEMM\n");
      tc_solve_impl_cutlass<TC::B1, TC_SUBMETHOD_B1, TCGemmOpMultiplyAdd>(
        is_first, n, m, k, // NOTE: switching order of A, B.
        (uint8_t*)tc_bufs.tc_buf_right, k,
        (uint8_t*)tc_bufs.tc_buf_left, k,
        (int32_t*)matC, m,
        env.stream_compute());

    } // if (TC_METHOD == TC::INT4)

    System::accel_last_call_succeeded(); // extra precaution.

  }
  // Call WMMA 1-bit GEMM kernels
  else if(env.num_kernel()>=1 && env.num_kernel()<100) {
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
        /*case 11: {
          if(env.print_details()) printf("Using Cutlass kernel 128x256x512\n");
          CutlassTCGemm1B<TCTBlockType::_128_256_512,
                          TCWarpType::_64_64_512,
			  TCInstType::_8_8_128,
			  TCOpType::Xor, 2>(
	    n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k,
            (uint8_t*)tc_bufs.tc_buf_left, k, beta, (int32_t*)matC, n, env.stream_compute());
        } break;
	case 12: {
          if(env.print_details()) printf("Using Cutlass kernel 256x128x512\n");
          CutlassTCGemm1B<TCTBlockType::_256_128_512,
		          TCWarpType::_64_64_512,
			  TCInstType::_8_8_128,
			  TCOpType::Xor, 2>(
	    n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k,
            (uint8_t*)tc_bufs.tc_buf_left, k, beta, (int32_t*)matC, n, env.stream_compute());
        } break;
        case 13: {
          if(env.print_details()) printf("Using Cutlass kernel 128x128\n");
          CutlassTCGemm1B<TCTBlockType::_128_128_512,
		          TCWarpType::_64_64_512,
			  TCInstType::_8_8_128,
			  TCOpType::Xor, 2>(
	    n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k,
            (uint8_t*)tc_bufs.tc_buf_left, k, beta, (int32_t*)matC, n, env.stream_compute());
        } break;
        case 14: {
          if(env.print_details()) printf("Using Cutlass WMMA kernel 64x64\n");
          CutlassTCGemm1BWmma<TCTBlockType::_64_64_512,
		              TCWarpType::_32_32_512,
			      TCInstType::_8_8_128,
			      TCOpType::Xor, 2>(
	    n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k,
            (uint8_t*)tc_bufs.tc_buf_left, k, beta, (int32_t*)matC, n, env.stream_compute());
        } break;*/
#if defined COMET_USE_AMPERE
        // Stages=2
	/*case 15: {
          if(env.print_details()) printf("Using Cutlass kernel TB=128x256x1024 W=64x64x1024 I=16x8x256 NStages=2\n");
          CutlassTCGemm1B<TCTBlockType::_128_256_1024,
                          TCWarpType::_64_64_1024,
                          TCInstType::_16_8_256,
			  TCOpType::Xor, 2>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;
        case 16: {
          if(env.print_details()) printf("Using Cutlass kernel TB=256x128x1024 W=64x64x1024 I=16x8x256 NStages=2\n");
          CutlassTCGemm1B<TCTBlockType::_256_128_1024,
                          TCWarpType::_64_64_1024,
                          TCInstType::_16_8_256,
			  TCOpType::Xor, 2>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;*/
	// Stages=3
	case 17: {
          if(env.print_details()) printf("Using Cutlass kernel TB=128x256x1024 W=64x64x1024 I=16x8x256 NStages=3\n");
          CutlassTCGemm1B<TCTBlockType::_128_256_1024,
                          TCWarpType::_64_64_1024,
                          TCInstType::_16_8_256,
			  TCOpType::Xor, 3>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;
	/*case 18: {
          if(env.print_details()) printf("Using Cutlass kernel TB=256x128x1024 W=64x64x1024 I=16x8x256 NStages=3\n");
          CutlassTCGemm1B<TCTBlockType::_256_128_1024,
                          TCWarpType::_64_64_1024,
                          TCInstType::_16_8_256,
			  TCOpType::Xor, 3>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;
	case 19: {
          if(env.print_details()) printf("Using Cutlass kernel TB=128x256x1024 W=32x64x1024 I=16x8x256 NStages=3\n");
          CutlassTCGemm1B<TCTBlockType::_128_256_1024,
                          TCWarpType::_64_32_1024,
                          TCInstType::_16_8_256,
			  TCOpType::Xor, 3>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;
	case 20: {
          if(env.print_details()) printf("Using Cutlass kernel TB=128x256x1024 W=32x32x1024 I=16x8x256 NStages=3\n");
          CutlassTCGemm1B<TCTBlockType::_128_128_1024,
                          TCWarpType::_64_64_1024,
                          TCInstType::_16_8_256,
			  TCOpType::Xor, 3>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;
	// Stages=4
	case 21: {
          if(env.print_details()) printf("Using Cutlass kernel TB=64x64x1024 W=32x32x1024 I=16x8x256 NStages=4\n");
          CutlassTCGemm1B<TCTBlockType::_64_64_1024,
                          TCWarpType::_32_32_1024,
                          TCInstType::_16_8_256,
			  TCOpType::Xor, 4>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;
	// Stages=6
	case 22: {
	  if(env.print_details()) printf("Using Cutlass xor kernel TB=64x64x512 W=32x32x512 I=16x8x256 NStages=6\n");
          CutlassTCGemm1B<TCTBlockType::_64_64_512,
                          TCWarpType::_32_32_512,
                          TCInstType::_16_8_256,
			  TCOpType::Xor, 6>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;
	// Larger unlisted sizes
	case 23: {
	  if(env.print_details()) printf("Using Cutlass kernel TB=128x256x1024 W=64x128x1024\n");
          CutlassTCGemm1BTest<TCTBlockType::_128_256_1024,
                          TCWarpType::_64_128_1024,
                          TCInstType::_16_8_256,
                          TCOpType::Xor, 3>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;
	case 24: {
          if(env.print_details()) printf("Using Cutlass kernel TB=128x256x1024 W=128x64x1024\n");
          CutlassTCGemm1BTest<TCTBlockType::_128_256_1024,
                          TCWarpType::_128_64_1024,
                          TCInstType::_16_8_256,
                          TCOpType::Xor, 3>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;
	case 25: {
          if(env.print_details()) printf("Using Cutlass kernel TB=128x256x1024 W=128x128x1024\n");
          CutlassTCGemm1BTest<TCTBlockType::_128_256_1024,
                          TCWarpType::_128_128_1024,
                          TCInstType::_16_8_256,
                          TCOpType::Xor, 3>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;
	case 26: {
          if(env.print_details()) printf("Using Cutlass kernel TB=128x288x1024 W=64x96x1024\n");
          CutlassTCGemm1BTest<TCTBlockType::_128_288_1024,
                          TCWarpType::_64_96_1024,
                          TCInstType::_16_8_256,
                          TCOpType::Xor, 3>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;*/
        // Reasonable performance but not as good for large problems
	/*case 27: {
          if(env.print_details()) printf("Using Cutlass kernel TB=160x256x1024 W=80x64x1024\n");
          CutlassTCGemm1BTest<TCTBlockType::_160_256_1024,
                          TCWarpType::_80_64_1024,
                          TCInstType::_16_8_256,
                          TCOpType::Xor, 3>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;
	case 28: {
          if(env.print_details()) printf("Using Cutlass kernel TB=160x288x1024 W=80x72x1024\n");
          CutlassTCGemm1BTest<TCTBlockType::_192_224_1024,
                          TCWarpType::_96_56_1024,
                          TCInstType::_16_8_256,
                          TCOpType::Xor, 3>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;*/
	/*case 29: {
          if(env.print_details()) printf("Using Cutlass kernel TB=160x288x1024 W=80x72x1024\n");
          CutlassTCGemm1BTest<TCTBlockType::_128_256_1024,
                          TCWarpType::_32_128_1024,
                          TCInstType::_16_8_256,
                          TCOpType::Xor, 3>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;*/

	// Other settings
	/*case 48: {
          if(env.print_details()) printf("Using Cutlass kernel TB=128x256x1024 W=64x64x1024 I=16x8x256 NStages=3 Alt\n");
          CutlassTCGemm1BTest<TCTBlockType::_128_256_1024,
                          TCWarpType::_64_64_1024,
                          TCInstType::_16_8_256,
                          TCOpType::Xor, 3>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;
        case 49: {
          if(env.print_details()) printf("Using Cutlass kernel TB=128x256x1024 W=64x64x1024 I=16x8x256 NStages=3 Alt\n");
          CutlassTCGemm1BSplitK<TCTBlockType::_128_256_1024,
                          TCWarpType::_64_64_1024,
                          TCInstType::_16_8_256,
                          TCOpType::Xor, 3>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;*/
#endif
        // Multiply GEMMs
#if defined COMET_USE_AMPERE
	// Stages=3
        case 50: {
          if(env.print_details()) printf("Using Cutlass mult kernel TB=128x256x1024 W=64x64x1024 I=16x8x256 NStages=3\n");
          CutlassTCGemm1B<TCTBlockType::_128_256_1024,
                          TCWarpType::_64_64_1024,
                          TCInstType::_16_8_256,
                          TCOpType::Mult, 3>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;

	case 51: {
          if(env.print_details()) printf("Using Cutlass kernel TB=256x128x1024 W=64x64x1024 I=16x8x256 NStages=3\n");
          CutlassTCGemm1B<TCTBlockType::_256_128_1024,
                          TCWarpType::_64_64_1024,
                          TCInstType::_16_8_256,
                          TCOpType::Mult, 3>(
            n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;
#endif


	default: {
          printf("Failed to call appropriate 1-bit GEMM kernel for num_kernel=%d\n",
             env.num_kernel());
          COMET_INSIST(false && "Failure to call GEMM function.");
        }
      }
      System::accel_last_call_succeeded();
      //if(env.print_details()) printf("Number of ops = 2*%d*%d*%d\n",m,n,k);
      //env.ops_local_inc(2 * m * (double)n * (double)k);

      // Print matrix contents
      cudaStreamSynchronize(env.stream_compute());
      printf("Printing matrix info\n");
      int m2=m/2, n2=n/2;
      const int threadblockx2 = BLOCK_SIZE, threadblocky2 = BLOCK_SIZE;
      int gridblockx2 = (int)ceil((double)m2/threadblockx2);
      int gridblocky2 = (int)ceil((double)n2/threadblocky2);
      printf("Printing matrix info mn=%d,%d m2n2=%d,%d\n",m,n,m2,n2);
      COMET_LAUNCH_KERNEL(b1_print_matrix,
        dim3(gridblockx2, gridblocky2, 1),
        dim3(threadblockx2, threadblocky2, 1), 0, env.stream_compute(),
        n2, m2, (int32_t*)matC);
      System::accel_last_call_succeeded();
  }
  else {
    printf("Failed to call appropriate 1-bit GEMM kernel for num_kernel=%d\n",
               env.num_kernel());
    COMET_INSIST(false && "Failure to call GEMM function.");
  } // if
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

  int rank = 0;
  COMET_MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

  if(env.print_details()) printf("\nrank=%d In tc_solve_impl mnk=%d,%d,%d num_kernel=%d\n",
    rank,m,n,k,env.num_kernel());
  if(env.print_details()) printf("rank=%d Starting timer\n",rank);

  // Make the appropriate BLAS call.

  const bool is_timing_gemm = false; //System::is_proc_num_0(); //false; // true;

  if (is_timing_gemm)
    env.stream_synchronize(env.stream_compute());
  double t1 = !is_timing_gemm ? 0 : System::time();

  //if (env.is_compute_method_gpu() &&
  //  (TC_METHOD == TC::B1 || TC_METHOD == TC::INT4 ||
  //   (TC_METHOD == TC::INT8 && BuildHas::CUDA && BuildHas::CUTLASS &&
  //    System::compute_capability() > 700))) {

  if (env.is_using_cutlass() || env.is_using_cutlass_mockup()) {

    //------------------------------
    // CASE: GPU, CUTLASS OR MOCKUP.
    //------------------------------

    enum {IS_B_FIELD_MAJOR = true};
    COMET_INSIST(IS_B_FIELD_MAJOR == tc_is_b_field_major(env));

#   ifdef COMET_USE_ACCEL

      if(env.print_details()) printf("rank=%d Calling tc_solve_impl_subbyte B1\n",rank);
      env.gemm_timer.record();

      env.gemm_timer.start();

      tc_solve_impl_subbyte<TC_METHOD>(is_first, m, n, k, matC, tc_bufs, env);

      env.gemm_timer.end();
      if(env.print_details()) printf("rank=%d Done calling tc_solve_impl_subbyte B1\n",rank);

#   else // COMET_USE_ACCEL

      COMET_INSIST(false && "Failure to call GEMM function.");

#   endif // COMET_USE_ACCEL

  } else if (env.is_compute_method_gpu()) { // && not cutlass or mockup

    //---------------------------
    // CASE: GPU, CUBLAS/ROCBLAS.
    //---------------------------

    // Make accelerator BLAS call.

#   ifdef COMET_USE_ACCEL

      // NOTE: from https://devblogs.nvidia.com/programming-tensor-cores-cuda-9/
      // "Invoke the GEMM, ensuring k, lda, ldb, and ldc are all multiples of 8, 
      //  and m is a multiple of 4"
      // "GEMMs that do not satisfy the above rules will fall back
      //  to a non-Tensor Core implementation"
      // See also https://docs.nvidia.com/cuda/cublas/index.html#cublas-gemmEx
      // NOTE: this may be relaxed for later versions of CUDA.

      // k (=nfl) is derived from padded-up npfl (multiple of 64), so always ok.
      COMET_INSIST(k % 8 == 0 && "Failed divisibility condition for tc gemm.");
      // since I_max_dim % 4 == 0; see tc_gemm_vaxis_divisibility_required()
      COMET_INSIST(m % 8 == 0 && "Failed divisibility condition for tc gemm.");
      // since nvl % 4 == 0; see tc_gemm_vaxis_divisibility_required()
      COMET_INSIST(n % 8 == 0 && "Failed divisibility condition for tc gemm.");

      const typename TCTraits<TC_METHOD>::GemmOut_t alpha = 1;
      const typename TCTraits<TC_METHOD>::GemmOut_t beta = is_first ? 0 : 1;

      //enum {IS_B_FIELD_MAJOR = TCTraits<TC_METHOD>::IS_B_FIELD_MAJOR};
      enum {IS_B_FIELD_MAJOR = false};
      COMET_INSIST(IS_B_FIELD_MAJOR == tc_is_b_field_major(env));

      if(env.print_details())
          printf("Launching Cublas/Rocblas GEMM kernel m=%d n=%d k=%d beta=%d\n",
                 m,n,k,(int)beta);

      env.gemm_timer.record();
      env.gemm_timer.start();

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

      env.gemm_timer.end();

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

  } else { // (!env.is_compute_method_gpu() ...) {

    //-----------------------
    // CASE: CPU.
    //-----------------------

    if (env.tc_eff() == TC::FP32) {

#     ifdef COMET_USE_CPUBLAS

        const float alpha = 1;
        const float beta = is_first ? 0 : 1;

        // Make CPU BLAS call.

        env.gemm_timer.record();
	env.gemm_timer.start();

        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans,
          m, n, k, alpha, (float*)tc_bufs.tc_buf_left, m,
          (float*)tc_bufs.tc_buf_right, n, beta, (float*)matC, m);

	env.gemm_timer.end();

#     else // COMET_USE_CPUBLAS

        COMET_INSIST(false && "Failure to call GEMM function.");

#     endif // COMET_USE_CPUBLAS

    } else { // if env.tc_eff()

      COMET_INSIST(false && "Failure to call GEMM function.");

    } // if env.tc_eff()

  } // if compute_method

  const double ops = 2.0 * (double)m * (double)n * (double)k;

  env.ops_gemm_local_inc(ops);
  env.ops_local_inc(ops);
  
  if (is_timing_gemm) {
    env.stream_synchronize(env.stream_compute());
    double t2 = System::time();
    const double t = t2 - t1;
    printf("%i %i %i   time %f TOP/s %f   tbeg %f tend %f\n",
           (int)m, (int)n, (int)k, t,
           (2 * m * (double)n * (double)k) / (t * 1e12), t1, t2);
  }

  if(env.print_details()) printf("rank=%d Done in tc_solve_impl\n",rank);  
}

//-----------------------------------------------------------------------------
/// \brief Call to perform required GEMM.

template<int TC_METHOD>
void tc_solve_(bool is_first, int nvll, int nvl, int npfl_thisstep,
               void* matC, TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST(matC);
  COMET_INSIST(nvll >= 0 && nvl >= 0 && nvll <= nvl);
  COMET_INSIST(npfl_thisstep >= 0);
  COMET_INSIST(env.tc_eff() != TC::NO);

  const int nfl_thisstep = npfl_thisstep * 64;

  const int m = 2 * nvll; // metrics array dim
  const int n = 2 * nvl; // metrics array dim
  const int k = nfl_thisstep; // vectors array (as GemmIn_t) dim

  if(env.print_details())
    printf("Calling tc_solve_impl with m=2*nvll=2*%d=%d n=2*nvl=2*%d=%d "
           "k=npfl_thisstep*64=%d*64=%d\n",nvll,m,nvl,n,npfl_thisstep,k);
  tc_solve_impl<TC_METHOD>(is_first, m, n, k, matC, tc_bufs, env);
}

//-----------------------------------------------------------------------------
/// \brief Call to perform required GEMM.

template<int TC_METHOD>
void tc_solve_comet_(bool is_first, int nvll, int nvl, int npvfl_thisstep,
                     const void *matA, const void *matB,void* matC, TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST(matC);
  COMET_INSIST(nvll >= 0 && nvl >= 0 && nvll <= nvl);
  COMET_INSIST(npvfl_thisstep >= 0);
  COMET_INSIST(env.tc_eff() != TC::NO);

  const int nfl_thisstep = npvfl_thisstep;
  //const int nfl_thisstep = npvfl_thisstep * 64;

  const int m = nvll; // metrics array dim
  const int n = nvl; // metrics array dim
  const int k = nfl_thisstep; // vectors array (as GemmIn_t) dim

  if(env.print_details()) printf("In tc_solve_comet_ calling tc_solve_comet_impl with mnk=%d,%d,%d nvll=%d nvl=%d\n",m,n,k,nvll,nvl);
  if(env.num_kernel() >= 100 && env.num_kernel() < 125) {
    tc_solve_comet_xor_impl<TC_METHOD>(is_first, m, n, k, matA, matB, matC, tc_bufs, env);
  } else if(env.num_kernel() >= 150 && env.num_kernel() < 175) {
    tc_solve_comet_mult_impl<TC_METHOD>(is_first, m, n, k, matA, matB, matC, tc_bufs, env);
  } else {
    printf("Failed to call appropriate 1-bit CoMet double GEMM kernel for num_kernel=%d\n",
         env.num_kernel());
      COMET_INSIST(false && "Failure to call GEMM function.");
  }
}

//-----------------------------------------------------------------------------
/// \brief Call to perform required GEMM.

template<int TC_METHOD>
void tc_solve_comet_int_(bool is_first, int nvll, int nvl, int npfl_thisstep, int pfl_min, int step_2way, int nfal,
                         const void *matA1, const void *matA2, const void *matB,void* matC, TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST(matC);
  COMET_INSIST(nvll >= 0 && nvl >= 0 && nvll <= nvl);
  COMET_INSIST(npfl_thisstep >= 0);
  COMET_INSIST(env.tc_eff() != TC::NO);

  const int nfl_thisstep = npfl_thisstep;

  const int m = nvll; // metrics array dim
  const int n = nvl; // metrics array dim
  const int k = nfl_thisstep; // vectors array (as GemmIn_t) dim

  if(env.num_kernel() >= 125 && env.num_kernel() < 150) {
    if(env.print_details()) printf("Calling tc_solve_comet_int_impl with mnk=%d,%d,%d\n",m,n,k);
    tc_solve_comet_int_impl<TC_METHOD>(is_first, m, n, k, matA1, matB, matC, tc_bufs, env);
  } else if(env.num_kernel() >= 175 && env.num_kernel() < 200) {
    tc_solve_comet_mult_int_impl<TC_METHOD>(is_first, m, n, k, pfl_min, step_2way, nfal, matA1, matA2, matB, matC, tc_bufs, env);
  } else {
    printf("Failed to call appropriate 1-bit CoMet int GEMM kernel for num_kernel=%d\n",
         env.num_kernel());
      COMET_INSIST(false && "Failure to call GEMM function.");
  }
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TC_SOLVE_I_HH_

//-----------------------------------------------------------------------------
