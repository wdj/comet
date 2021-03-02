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
        _64_64_WMMA = 7
  };
};

struct TCGemmOpXorPopc {
  template<typename GemmIn_t>
  __host__ __device__
  static GemmIn_t op(GemmIn_t a, GemmIn_t b) {return a ^ b;}
# ifdef COMET_USE_CUTLASS
  typedef cutlass::arch::OpXorPopc Value; 
#else
  typedef int Value; 
#endif
};

struct TCGemmOpMultiplyAdd {
  template<typename GemmIn_t>
  __host__ __device__
  static GemmIn_t op(GemmIn_t a, GemmIn_t b) {return a & b;}
# ifdef COMET_USE_CUTLASS
  typedef cutlass::arch::OpMultiplyAdd Value; 
#else
  typedef int Value; 
#endif
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

template<int TC_SUBMETHOD> struct CutlassSettings;

template<> struct CutlassSettings<TCSubmethod::_256_128>
  : public CutlassOpClassTensorOp {
  enum {ThreadBlockShape0 = 256,
        ThreadBlockShape1 = 128,
        WarpShape0 = 64,
        WarpShape1 = 64
  };
};

template<> struct CutlassSettings<TCSubmethod::_128_256>
  : public CutlassOpClassTensorOp {
  enum {ThreadBlockShape0 = 128,
        ThreadBlockShape1 = 256,
        WarpShape0 = 64,
        WarpShape1 = 64
  };
};

template<> struct CutlassSettings<TCSubmethod::_128_128>
  : public CutlassOpClassTensorOp {
  enum {ThreadBlockShape0 = 128,
        ThreadBlockShape1 = 128,
        WarpShape0 = 64,
        WarpShape1 = 64
  };
};

template<> struct CutlassSettings<TCSubmethod::_128_64>
  : public CutlassOpClassTensorOp {
  enum {ThreadBlockShape0 = 128,
        ThreadBlockShape1 = 64,
        WarpShape0 = 64,
        WarpShape1 = 32
  };
};

template<> struct CutlassSettings<TCSubmethod::_64_128>
  : public CutlassOpClassTensorOp {
  enum {ThreadBlockShape0 = 64,
        ThreadBlockShape1 = 128,
        WarpShape0 = 32,
        WarpShape1 = 64
  };
};

template<> struct CutlassSettings<TCSubmethod::_64_64>
  : public CutlassOpClassTensorOp {
  enum {ThreadBlockShape0 = 64,
        ThreadBlockShape1 = 64,
        WarpShape0 = 32,
        WarpShape1 = 32
  };
};

template<> struct CutlassSettings<TCSubmethod::_64_64_WMMA>
  : public CutlassOpClassWmmaTensorOp {
  enum {ThreadBlockShape0 = 64,
        ThreadBlockShape1 = 64,
        WarpShape0 = 32,
        WarpShape1 = 32
  };
};

//-----------------------------------------------------------------------------

template<int TC_SUBMETHOD, typename TCGemmOp>
void tc_solve_impl_b1_cutlass(
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

  using ElementInput = cutlass::uint1b_t;
  using ElementOutput = int32_t;
  using ElementAccumulator = int32_t;
  using ElementCompute = int32_t;

  // NOTE: COMET_CUTLASS_ARCH is a #define
  typedef typename cutlass::arch::COMET_CUTLASS_ARCH CutlassArch_t;

  // see https://github.com/NVIDIA/cutlass/blob/master/include/cutlass/gemm/device/gemm.h
  // https://github.com/NVIDIA/cutlass/blob/master/include/cutlass/gemm/device/default_gemm_configuration.h

  using Gemm = cutlass::gemm::device::Gemm<
      ElementInput, cutlass::layout::RowMajor,
      ElementInput, cutlass::layout::ColumnMajor,
      ElementOutput, cutlass::layout::RowMajor,
      ElementAccumulator,
      typename CutlassSettings<TC_SUBMETHOD>::OpClass_t,
      CutlassArch_t,
      cutlass::gemm::GemmShape< // ThreadblockShape_
        CutlassSettings<TC_SUBMETHOD>::ThreadBlockShape0,
        CutlassSettings<TC_SUBMETHOD>::ThreadBlockShape1,
        512>,
      cutlass::gemm::GemmShape< // WarpShape_
        CutlassSettings<TC_SUBMETHOD>::WarpShape0,
        CutlassSettings<TC_SUBMETHOD>::WarpShape1,
        512>,
      //cutlass::gemm::GemmShape<8, 8, 128>, // InstructionShape_
      cutlass::gemm::GemmShape<16, 8, 256>, // InstructionShape_
      cutlass::epilogue::thread::LinearCombination<
          ElementOutput,
          128 / cutlass::sizeof_bits<ElementOutput>::value,
          ElementAccumulator,
          ElementCompute>,
      cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
      //2, // Stages
      3, // Stages
      128, // AlignmentA
      128, // AlignmentB
      false, // SplitKSerial
      typename TCGemmOp::Value>;
      //cutlass::arch::OpXorPopc>;
      //cutlass::arch::OpMultiplyAdd>;

  Gemm gemm_operator;

  const int32_t alpha = 1;
  const int32_t beta = is_first ? 0 : 1;

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

#if 0
//=============================================================================
/// \brief 1-bit xor gemm kernel (mockup version, not high performance).

template<typename GemmIn_t, typename GemmOut_t>
__global__ void tc_solve_b1_xor_gemm_mockup_kernel(
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

    // Use xor; count 1-bits with popcount.
    const GemmOut_t v = utils::popc<GemmIn_t>(aik ^ bjk);
    if (ind_m < 1024 * 1024 * 1024) // WORKAROUND for undiagnosed error.
      cij = beta || ind_k ? cij + v : v;

  } // for ind_k
}
#endif

//-----------------------------------------------------------------------------
/// \brief 1-bit xor gemm.

static bool tc_solve_b1_use_mockup(const CEnv& env) {
  //return true;
  return env.tc_eff() == TC::B1 &&
    ! (BuildHas::CUDA && System::compute_capability() > 700);
}

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

template<int TC_METHOD>
static void tc_solve_impl_b1(bool is_first, int m, int n, int k,
  void* matC, TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST(matC);
  COMET_INSIST(m >= 0 && n >= 0 && k >= 0);

  if(env.print_details()) printf("In tc_solve_impl_b1 num_kernel=%d use_mockup=%d\n",
    env.num_kernel(),tc_solve_b1_use_mockup(env));
  if (env.num_kernel()==0) { // && tc_solve_b1_use_mockup(env)) {

    COMET_INSIST(TCTraits<TC_METHOD>::IS_B_FIELD_MAJOR);

    enum {NUM_FL_PER_PFL = 64};
    COMET_INSIST(k % NUM_FL_PER_PFL == 0 &&
                 "Failed divisibility condition for tc gemm.");

    const bool beta = is_first ? 0 : 1;

    // 8 == number of uint8_t values used to store each chunk of
    // NUM_FL_PER_PFL fields in the tc buf.
    enum {BYTES_PER_PFL_FIELDS = 8};

    typedef typename TCTraits<TC_METHOD>::GemmIn_t GemmIn_t;
    typedef typename TCTraits<TC_METHOD>::GemmOut_t GemmOut_t;

    const int bytes_per_gi = sizeof(GemmIn_t);
    const size_t k_eff = (k / NUM_FL_PER_PFL) *
                         (BYTES_PER_PFL_FIELDS / bytes_per_gi);

    const int threadblocksize = 256;
    COMET_INSIST((threadblocksize <= 256 || ! BuildHas::HIP) &&
                 "Current HIP limitation.");
    const int num_threadblocks_0 = utils::ceil(m, threadblocksize);
    const int num_threadblocks_1 = n;

    if (env.is_using_xor()) {

      COMET_LAUNCH_KERNEL((tc_solve_b1_gemm_mockup_kernel
                           <GemmIn_t, GemmOut_t, TCGemmOpXorPopc>),
        dim3(num_threadblocks_0, num_threadblocks_1, 1),
        dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
        m, n, k_eff, (GemmIn_t*)tc_bufs.tc_buf_left,
        (GemmIn_t*)tc_bufs.tc_buf_right, beta, (GemmOut_t*)matC);

    } else { // !env.is_using_xor()

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

    } // if (env.is_using_xor())

    System::accel_last_call_succeeded();

  } else if(env.num_kernel()==10) {

    if (env.is_using_xor()) {

      if(env.print_details()) printf("Using Xor Cutlass GEMM\n");
      tc_solve_impl_b1_cutlass<TCSubmethod::_128_256, TCGemmOpXorPopc>(
        is_first, n, m, k, // NOTE: switching order of A, B.
        (uint8_t*)tc_bufs.tc_buf_right, k,
        (uint8_t*)tc_bufs.tc_buf_left, k,
        (int32_t*)matC, m,
        env.stream_compute());

    } else { // !env.is_using_xor()

      if(env.print_details()) printf("Using MultAdd Cutlass GEMM\n");
      tc_solve_impl_b1_cutlass<TCSubmethod::_128_256, TCGemmOpMultiplyAdd>(
        is_first, n, m, k, // NOTE: switching order of A, B.
        (uint8_t*)tc_bufs.tc_buf_right, k,
        (uint8_t*)tc_bufs.tc_buf_left, k,
        (int32_t*)matC, m,
        env.stream_compute());


    } // if (env.is_using_xor())

    System::accel_last_call_succeeded(); // extra precaution.

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
        case 11: {
          if(env.print_details()) printf("Using Cutlass kernel 128x256x512\n");
          CutlassTCGemm1B<TCTBlockType::_128_256_512,
                          TCWarpType::_64_64_512,
			  TCInstType::_8_8_128,2>(
	    n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k,
            (uint8_t*)tc_bufs.tc_buf_left, k, beta, (int32_t*)matC, n, env.stream_compute());
        } break;
	case 12: {
          if(env.print_details()) printf("Using Cutlass kernel 256x128x512\n");
          CutlassTCGemm1B<TCTBlockType::_256_128_512,
		          TCWarpType::_64_64_512,
			  TCInstType::_8_8_128,2>(
	    n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k,
            (uint8_t*)tc_bufs.tc_buf_left, k, beta, (int32_t*)matC, n, env.stream_compute());
        } break;
        case 13: {
          if(env.print_details()) printf("Using Cutlass kernel 128x128\n");
          CutlassTCGemm1B<TCTBlockType::_128_128_512,
		          TCWarpType::_64_64_512,
			  TCInstType::_8_8_128,2>(
	    n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k,
            (uint8_t*)tc_bufs.tc_buf_left, k, beta, (int32_t*)matC, n, env.stream_compute());
        } break;
        case 14: {
          if(env.print_details()) printf("Using Cutlass WMMA kernel 64x64\n");
          CutlassTCGemm1BWmma<TCTBlockType::_64_64_512,
		              TCWarpType::_32_32_512,
			      TCInstType::_8_8_128,2>(
	    n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k,
            (uint8_t*)tc_bufs.tc_buf_left, k, beta, (int32_t*)matC, n, env.stream_compute());
        } break;
        case 17: {
          if(env.print_details()) printf("Using Cutlass kernel TB=128x256x1024 W=64x64x1024 I=16x8x256 NStages=2\n");
          CutlassTCGemm1B<TCTBlockType::_128_256_1024,
		          TCWarpType::_64_64_1024,
			  TCInstType::_16_8_256,2>(
	    n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;
        case 18: {
          if(env.print_details()) printf("Using Cutlass kernel TB=128x256x1024 W=64x64x1024 I=16x8x256 NStages=3\n");
          CutlassTCGemm1B<TCTBlockType::_128_256_1024,
		          TCWarpType::_64_64_1024,
			  TCInstType::_16_8_256,3>(
	    n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;
	// Stages = 4 is invalid argument
	/*case 19: {
          if(env.print_details()) printf("Using Cutlass kernel TB=128x256x1024 W=64x64x1024 I=16x8x256 NStages=4\n");
          CutlassTCGemm1B<TCTBlockType::_128_256_1024,
		          TCWarpType::_64_64_1024,
			  TCInstType::_16_8_256,4>(
	    n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k, (uint8_t*)tc_bufs.tc_buf_left, k, beta,
            (int32_t*)matC, n, env.stream_compute());
        } break;*/
	/*case 30: {
          if(env.print_details()) printf("Using Cutlass kernel 128x256\n");
          CutlassTCGemm1B_128x256(n, m, k, (uint8_t*)tc_bufs.tc_buf_right, k,
            (uint8_t*)tc_bufs.tc_buf_left, k, (int32_t*)matC, n, env.stream_compute());
        } break;*/
        default: {
          printf("Failed to call appropriate 1-bit GEMM kernel for num_kernel=%d\n",
             env.num_kernel());
          COMET_INSIST(false && "Failure to call GEMM function.");
        }
      }
      System::accel_last_call_succeeded();
      //if(env.print_details()) printf("Number of ops = 2*%d*%d*%d\n",m,n,k);
      //env.ops_local_inc(2 * m * (double)n * (double)k);
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
  double tbegin = env.get_time();
  if(env.print_details()) printf("rank=%d Starting timer\n",rank);

  // Make the appropriate BLAS call.

  const bool is_timing_gemm = false; //true;

  if (is_timing_gemm)
    env.stream_synchronize(env.stream_compute());
  double t1 = !is_timing_gemm ? 0 : System::time();

  if (env.is_compute_method_gpu() && TC_METHOD == TC::B1) {

    //-------------------
    // CASE: GPU, TC::B1.
    //-------------------

#   ifdef COMET_USE_ACCEL

      //env.gemmtime_record();

      //env.gemmtime_start();
      if(env.print_details()) printf("rank=%d Calling tc_solve_impl_b1\n",rank);
      tc_solve_impl_b1<TC_METHOD>(is_first, m, n, k, matC, tc_bufs, env);
      if(env.print_details()) printf("rank=%d Done calling tc_solve_impl_b1\n",rank);

      //env.gemmtime_end();

#   else // COMET_USE_ACCEL

      COMET_INSIST(false && "Failure to call GEMM function.");

#   endif // COMET_USE_ACCEL

  } else if (env.is_compute_method_gpu()) { // && TC_METHOD != TC::B1

    //-----------------------
    // CASE: GPU, non-TC::B1.
    //-----------------------

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

      enum {IS_B_FIELD_MAJOR = TCTraits<TC_METHOD>::IS_B_FIELD_MAJOR};

      if(env.print_details())
          printf("Launching Cublas/Rocblas GEMM kernel m=%d n=%d k=%d beta=%d\n",
                 m,n,k,(int)beta);

      env.gemmtime_record();
      env.gemmtime_start();

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

      env.gemmtime_end();

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

    //-----------------------
    // CASE: CPU.
    //-----------------------

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

    } else { // if env.tc_eff()

      COMET_INSIST(false && "Failure to call GEMM function.");

    } // if env.tc_eff()

  } // if compute_method

  const double ops = 2.0 * (double)m * (double)n * (double)k;

  env.ops_gemm_local_inc(ops);
  env.ops_local_inc(ops);
  env.gemmtime_inc(env.get_time() - tbegin);

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

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TC_SOLVE_I_HH_

//-----------------------------------------------------------------------------
