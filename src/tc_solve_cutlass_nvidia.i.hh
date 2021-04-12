//-----------------------------------------------------------------------------
/*!
 * \file   tc_solve_cutlass_nvidia.i.hh
 * \author Paul Eller
 * \date   Tue Nov  3 08:26:29 EST 2020
 * \brief  CUDA code, gemm operation, 
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

#ifndef _COMET_TC_SOLVE_CUTLASS_NVIDIA_I_HH_
#define _COMET_TC_SOLVE_CUTLASS_NVIDIA_I_HH_

#ifdef COMET_USE_CUTLASS

#define WMMA1B_M 8
#define WMMA1B_N 8
#define WMMA1B_K 128

#include <cstdint>

#include <cutlass/arch/wmma.h>

#include <cutlass/gemm/warp/default_mma_wmma_tensor_op.h>
#include <cutlass/gemm/warp/default_mma_tensor_op.h>

#include <cutlass/matrix_shape.h>

#include <cutlass/gemm/gemm.h>
#include <cutlass/gemm/warp/mma.h>

#include <cutlass/gemm/warp/mma_tensor_op_policy.h>

#include <cutlass/gemm/warp/mma_tensor_op_tile_iterator.h>
#include <cutlass/gemm/warp/mma_tensor_op_tile_iterator_sm80.h>

#include <cutlass/cutlass.h>
#include <cutlass/aligned_buffer.h>
#include <cutlass/subbyte_reference.h>
#include <cutlass/platform/platform.h>

// GPU bit count routine
#define gm_popcount64(x) __popcll(x)

// Cuda block size
#define BLOCK_SIZE 8

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// 

__device__ inline void process_bits(GMBits2x64 vi, uint64_t &nvi, uint64_t &pvi)
{ 
  //--------------------
  // Nomenclature:
  // 
  // ( )v(i)(0)_(0)
  // (n)v(j)(1)_(1)
  //  ^   ^  ^   ^
  //  |   |  |   |--- lower or upper bit of each seminibble
  //  |   |  |--- lower or upper word
  //  |   |--- left or right vector 
  //  |---test for value or for its negative/complement
  //--------------------
  const uint64_t vi0 = vi.data[0];
  const uint64_t vi1 = vi.data[1];
    
  // Compute masks to sample the single needed bit from each seminibble,
  // and to ignore undefined vector entries.
  const uint64_t oddbits = 0x5555555555555555; 
  const uint64_t vi0mask = (vi0 | ~(vi0 >> 1)) & oddbits;
  const uint64_t vi1mask = (vi1 | ~(vi1 >> 1)) & oddbits;
    
  // Extract elts that are a "1" bit (=01).
  const uint64_t pvi0 = vi0 & vi0mask;
  const uint64_t pvi1 = vi1 & vi1mask;
    
  // Extract elts that are an "0" bit (=00).
  const uint64_t nvi0 = ~vi0 & vi0mask;
  const uint64_t nvi1 = ~vi1 & vi1mask;
    
  // Combine lower, upper words - each only uses odd bits - make packed.
  pvi = pvi0 | (pvi1 << 1);
  nvi = nvi0 | (nvi1 << 1);

  //printf("pvi0=%lu pvi1=%lu pvi=%lu\n",pvi0,pvi1,pvi);
}

//-----------------------------------------------------------------------------
/// 

__device__ inline void combine_bits(int nn, int np, int pn, int pp, double &c0, double &c1)
{
  enum { GM_TALLY1_MAX_VALUE_BITS = 26 };

  const uint64_t r00 = nn;
  const uint64_t r01 = np;
  const uint64_t r10 = pn;
  const uint64_t r11 = pp;

  c0 = r00 | (r10 << GM_TALLY1_MAX_VALUE_BITS);
  c1 = r01 | (r11 << GM_TALLY1_MAX_VALUE_BITS);
  //c0 = r00 | (r01 << GM_TALLY1_MAX_VALUE_BITS);
  //c1 = r10 | (r11 << GM_TALLY1_MAX_VALUE_BITS);
}

//-----------------------------------------------------------------------------
/// 

__device__ inline void combine_bits_sum(int nn, int np, int pn, int pp, double &c0, double &c1,
                                        double &cc0, double &cc1)
{
  enum { GM_TALLY1_MAX_VALUE_BITS = 26 };

  const uint64_t r00 = nn;
  const uint64_t r01 = np;
  const uint64_t r10 = pn;
  const uint64_t r11 = pp;

  c0 = r00 | (r10 << GM_TALLY1_MAX_VALUE_BITS);
  c1 = r01 | (r11 << GM_TALLY1_MAX_VALUE_BITS);
  //c0 = r00 | (r01 << GM_TALLY1_MAX_VALUE_BITS);
  //c1 = r10 | (r11 << GM_TALLY1_MAX_VALUE_BITS);

  cc0 += c0;
  cc1 += c1;
}

//-----------------------------------------------------------------------------
///

template <int pM, int pK, int block_x, int block_y>
__device__ inline void g2r(GMBits2x64 *gmem, int begin, int k, uint64_t frag[])
{
  constexpr int iter_y = pM / block_y;
  constexpr int iter_x = pK / block_x;

  constexpr int iter_xy = iter_x * iter_y;

#pragma unroll
  for (int y = 0; y < iter_y; y++) {
#pragma unroll
    for (int x = 0; x < iter_x; x++) {

      int nrow = threadIdx.y + y * block_y;
      int kval = threadIdx.x + x * block_x;

      // A row major and B col major
      int ind = begin + k * nrow + kval;

      //---Extract input values to process---
      // GMBits2x64 vi = a[aInd];
      auto gmem_u2 = reinterpret_cast<ulonglong2 *>(gmem);
      auto u2 = gmem_u2[ind];
      GMBits2x64 vi;
      vi.data[0] = u2.x;
      vi.data[1] = u2.y;

      uint64_t nvi;
      uint64_t pvi;
      process_bits(vi, nvi, pvi);

      frag[y * iter_x + x] = nvi;
      frag[y * iter_x + x + iter_xy] = pvi;
    }
  }
}

//-----------------------------------------------------------------------------
/// 

template <int pM, int pK, int block_x, int block_y, bool is_a, class Layout>
__device__ inline void r2s(uint64_t frag[], uint64_t *smem, Layout layout)
{
  constexpr int iter_y = pM / block_y;
  constexpr int iter_x = pK / block_x;

  constexpr int iter_xy = iter_x * iter_y;

#pragma unroll
  for (int y = 0; y < iter_y; y++) {
#pragma unroll
    for (int x = 0; x < iter_x; x++) {

      int nrow = threadIdx.y + y * block_y;
      int kval = threadIdx.x + x * block_x;

      int nbrow = nrow;
      int nbcol = kval;

      const int patch = 64; // sizeof(uint64_t) / sizeof(cutlass::uint1b_t)
      if (is_a) {
        // operand A
        smem[layout({nbrow, nbcol * patch}) / patch] = frag[y * iter_x + x];
        smem[layout({nbrow + pM, nbcol * patch}) / patch] = frag[y * iter_x + x + iter_xy];
      } else {
        // operand B
        smem[layout({nbcol * patch, nbrow}) / patch] = frag[y * iter_x + x];
        smem[layout({nbcol * patch, nbrow + pM}) / patch] = frag[y * iter_x + x + iter_xy];
      }
    }
  }
}

//-----------------------------------------------------------------------------
/// 

template <int pM, int pK, int block_x, int block_y, bool is_a, class Layout>
__device__ inline void g2s(uint64_t *smem, GMBits2x64 *gmem, int begin, int k, Layout layout)
{
#pragma unroll
  for (int y = 0; y < pM; y += block_y) {
#pragma unroll
    for (int x = 0; x < pK; x += block_x) {

      int nrow = threadIdx.y + y;
      int kval = threadIdx.x + x;

      // A row major and B col major
      int ind = begin + k * nrow + kval;

      //---Extract input values to process---
      // GMBits2x64 vi = a[aInd];
      auto gmem_u2 = reinterpret_cast<ulonglong2 *>(gmem);
      auto u2 = gmem_u2[ind];
      GMBits2x64 vi;
      vi.data[0] = u2.x;
      vi.data[1] = u2.y;

      int nbrow = nrow;
      int nbcol = kval;

      const int patch = 64; // sizeof(uint64_t) / sizeof(cutlass::uint1b_t)

      uint64_t nvi;
      uint64_t pvi;
      process_bits(vi, nvi, pvi);

      if (is_a) {
        // operand A
        smem[layout({nbrow, nbcol * patch}) / patch] = nvi;
        smem[layout({nbrow + pM, nbcol * patch}) / patch] = pvi;
      } else {
        // operand B
        smem[layout({nbcol * patch, nbrow}) / patch] = nvi;
        smem[layout({nbcol * patch, nbrow + pM}) / patch] = pvi;
      }
    }
  }
}

//-----------------------------------------------------------------------------
/// 

template <int bM, int bN, int bK, int wM, int wN, int wK, int block_x, int block_y>
__global__ void __launch_bounds__(block_x *block_y, 1)
b1_comet_xor_gemm_gpu_cutlass(int m, int n, int k, GMBits2x64 *a, GMBits2x64 *b, bool beta, GMTally2x2 *c)
{
  /**
    bM, bN, bK - threadblock MMA shape: bK is expressed in bits
    pM, pN, pK - threadblock MMA shape: bK is expressed in uint64_t's, and pM and pN are the shapes for loading GMBits2x64

    When loading data from gmem to smem, the "n" and "p" parts are stored to separate parts of smem, and then
    each warp works on and accumulates on the "nn", "np", "pn" and "pp" results separately. At the end each
    warp combines results stored in the 4 accumualtes and then write to gmem, which does not need to go through
    smem anymore.

    In other words, in each threadblock, this specific MMA problem is decomposed into 4 natural sub-problems.
    */

  // Block indices
  int bx = blockIdx.x, by = blockIdx.y;

  //if(bx==0 && by==0 && threadIdx.x==0 && threadIdx.y==0)
  //  printf("In b1_comet_xor_gemm_gpu_cutlass\n");

  int thread_idx = threadIdx.y * blockDim.x + threadIdx.x;

  int warp_id = thread_idx / 32;
  int lane_id = thread_idx % 32;

  constexpr int pM = bM / 2;
  constexpr int pN = bN / 2;
  constexpr int pK = bK / 64;

  // Matrix block location
  int aBegin = k * pM * bx;
  int bBegin = k * pN * by;

#if (__CUDA_ARCH__ < 800)
  // mma for sm75
  constexpr int alignment = 512;
  using InstructionShape = cutlass::gemm::GemmShape<8, 8, 128>;
#else
  // mma for sm80
  constexpr int alignment = 1024;
  using InstructionShape = cutlass::gemm::GemmShape<16, 8, 256>;
#endif

  static_assert(wK % alignment == 0, "alignment");

  using WarpShape = cutlass::gemm::GemmShape<wM, wN, wK>;
  using ElementA = cutlass::uint1b_t;
  using ElementB = cutlass::uint1b_t;
  using ElementC = int;
  using LayoutA
    = cutlass::layout::RowMajorTensorOpMultiplicandCrosswise<cutlass::sizeof_bits<ElementA>::value, alignment>;
  using LayoutB
    = cutlass::layout::ColumnMajorTensorOpMultiplicandCrosswise<cutlass::sizeof_bits<ElementB>::value, alignment>;
  using LayoutC = cutlass::layout::RowMajor;

  using Mma =
    typename cutlass::gemm::warp::DefaultMmaTensorOp<WarpShape, InstructionShape, ElementA, LayoutA, ElementB,
                                                     LayoutB, ElementC, LayoutC, cutlass::arch::OpXorPopc>::Type;

  using ThreadblockShape = cutlass::gemm::GemmShape<bM, bN, bK>;

  //
  // Distribute each sub-problem to the warps.
  //
  constexpr int warp_count_m = pM / WarpShape::kM;
  constexpr int warp_count_n = pN / WarpShape::kN;

  int warp_idx_mn = warp_id % (warp_count_m * warp_count_n);
  int warp_idx_m = warp_idx_mn % warp_count_m;
  int warp_idx_n = warp_idx_mn / warp_count_m;

  extern __shared__ char smem_ptr[];

  char *smem_ptr_A = smem_ptr;
  char *smem_ptr_B = &smem_ptr_A[ThreadblockShape::kM * ThreadblockShape::kK * 1 / 8];

  uint64_t *smem_buffer_A = reinterpret_cast<uint64_t *>(smem_ptr_A);
  uint64_t *smem_buffer_B = reinterpret_cast<uint64_t *>(smem_ptr_B);

  using FragmentA = typename Mma::FragmentA;
  using FragmentB = typename Mma::FragmentB;
  using FragmentC = typename Mma::FragmentC;

  typename Mma::LayoutA layout_A = Mma::LayoutA::packed({ThreadblockShape::kM, ThreadblockShape::kK});
  typename Mma::LayoutB layout_B = Mma::LayoutB::packed({ThreadblockShape::kK, ThreadblockShape::kN});
  typename Mma::LayoutC layout_C = Mma::LayoutC::packed({m, n});

  FragmentC accum_nn, accum_np;
  FragmentC accum_pn, accum_pp;

  accum_nn.clear();
  accum_np.clear();
  accum_pn.clear();
  accum_pp.clear();

  static_assert(pM % block_y == 0, "block-global-memory-loader needs to be in shape.");
  static_assert(pN % block_y == 0, "block-global-memory-loader needs to be in shape.");
  static_assert(pK % block_x == 0, "block-global-memory-loader needs to be in shape.");

  constexpr int iter_x = pK / block_x;
  constexpr int iter_y_a = pM / block_y;
  constexpr int iter_xy_a = iter_x * iter_y_a;
  constexpr int iter_y_b = pN / block_y;
  constexpr int iter_xy_b = iter_x * iter_y_b;

  uint64_t frag_a[2 * iter_xy_a];
  uint64_t frag_b[2 * iter_xy_b];

  for (int l = 0; l < k; l += pK) {
    // Here gmem -> register == "g2r" and then register -> smem == "r2s".
    // Operand A
    g2r<pM, pK, block_x, block_y>(a, aBegin + l, k, frag_a);
    // Operand B
    g2r<pN, pK, block_x, block_y>(b, bBegin + l, k, frag_b);
    // Operand A
    r2s<pM, pK, block_x, block_y, 1>(frag_a, smem_buffer_A, layout_A);
    // Operand B
    r2s<pN, pK, block_x, block_y, 0>(frag_b, smem_buffer_B, layout_B);

    //g2s<pM, pK, block_x, block_y, 1>(smem_buffer_A, a, aBegin+l, k, layout_A);
    //g2s<pN, pK, block_x, block_y, 0>(smem_buffer_B, b, bBegin+l, k, layout_B);

    __syncthreads();

    //
    // Construct warp-level matrix product
    //

    typename Mma::IteratorA iter_A({reinterpret_cast<typename cutlass::uint1b_t *>(smem_buffer_A), layout_A}, lane_id);
    typename Mma::IteratorB iter_B({reinterpret_cast<typename cutlass::uint1b_t *>(smem_buffer_B), layout_B}, lane_id);

    iter_A.add_tile_offset({warp_idx_m, 0});
    iter_B.add_tile_offset({0, warp_idx_n});

    FragmentA frag_A_n, frag_A_p;
    FragmentB frag_B_n, frag_B_p;

    if (warp_id < warp_count_n * warp_count_m) {

      Mma mma;

      CUTLASS_PRAGMA_UNROLL
      for (int warp_k = 0; warp_k < ThreadblockShape::kK; warp_k += Mma::Policy::MmaShape::kK) {
        iter_A.load(frag_A_n);
        iter_A.add_tile_offset({+warp_count_m, 0});
        iter_A.load(frag_A_p);
        iter_A.add_tile_offset({-warp_count_m, 0});

        iter_B.load(frag_B_n);
        iter_B.add_tile_offset({0, +warp_count_n});
        iter_B.load(frag_B_p);
        iter_B.add_tile_offset({0, -warp_count_n});

        ++iter_A;
        ++iter_B;

        mma(accum_nn, frag_A_n, frag_B_n, accum_nn);
        mma(accum_np, frag_A_n, frag_B_p, accum_np);
        mma(accum_pn, frag_A_p, frag_B_n, accum_pn);
        mma(accum_pp, frag_A_p, frag_B_p, accum_pp);
      }
    }

    if (l + pK < k) { __syncthreads(); }
  }

  if (warp_id < warp_count_n * warp_count_m) {
    // use "double2" instead of "GMTally2x2" to make sure compiler issue 128-bit store instructions.
    using output_type = double2;

    using IteratorTallyC = typename cutlass::gemm::warp::MmaTensorOpAccumulatorTileIterator<
      typename cutlass::MatrixShape<WarpShape::kM, WarpShape::kN>, output_type, LayoutC, InstructionShape,
      typename Mma::Policy::OpDelta>;

    typename IteratorTallyC::Fragment accum_tally;
    IteratorTallyC iter_tally_C({reinterpret_cast<output_type *>(c), layout_C}, lane_id);

    if(true) { //!beta) {
#pragma unroll
      for (int idx = 0; idx < FragmentC::kElements; idx++) {
        // Combine the results in the sub-problems.
        combine_bits(accum_nn[idx], accum_np[idx], accum_pn[idx], accum_pp[idx], accum_tally[idx].x, accum_tally[idx].y);
      }

      iter_tally_C.add_tile_offset({(pM / wM) * bx + warp_idx_m, (pN / wN) * by + warp_idx_n});
      //iter_tally_C.add_tile_offset({(pN / wN) * by + warp_idx_n, (pM / wM) * bx + warp_idx_m});
      // The following code does not translates into the most efficient instructions, even with 128-bit stores.
      // With current version of CUTLASS this is as far as we can go.
      iter_tally_C.store(accum_tally);
    }
    else {
      // Add Accum to C then store
      using IteratorC = typename cutlass::gemm::warp::MmaTensorOpAccumulatorTileIterator<
        typename cutlass::MatrixShape<WarpShape::kM, WarpShape::kN>, output_type, LayoutC, InstructionShape,
        typename Mma::Policy::OpDelta>;

      typename IteratorC::Fragment frag_C;
      IteratorC iter_C({reinterpret_cast<output_type *>(c), layout_C}, lane_id);
      iter_C.add_tile_offset({(pM / wM) * bx + warp_idx_m, (pN / wN) * by + warp_idx_n});
      iter_C.load(frag_C);

#pragma unroll
      for (int idx = 0; idx < FragmentC::kElements; idx++) {
        // Combine the results in the sub-problems.
        combine_bits_sum(accum_nn[idx], accum_np[idx], accum_pn[idx], accum_pp[idx],
                         accum_tally[idx].x, accum_tally[idx].y, frag_C[idx].x, frag_C[idx].y);
      }

      iter_C.store(frag_C);
    }
  }
}

//-----------------------------------------------------------------------------
/// 

void set_max_shared_bytes(const void *func)
{
  cudaFuncSetAttribute(func, cudaFuncAttributePreferredSharedMemoryCarveout, (int)cudaSharedmemCarveoutMaxShared);
  int max_shared_bytes;
  cudaDeviceGetAttribute(&max_shared_bytes, cudaDevAttrMaxSharedMemoryPerBlockOptin, 0);
  cudaFuncSetAttribute(func, cudaFuncAttributeMaxDynamicSharedMemorySize, max_shared_bytes);
  //printf("max_shared_bytes=%d\n",max_shared_bytes);
}

//-----------------------------------------------------------------------------
/// \brief Perform required GEMM.

void tc_solve_comet_impl_cutlass(int m, int n, int k, const void *matA, const void *matB, bool beta, void *matC)
{
#if defined COMET_USE_TURING
  // Use following for Turing
  constexpr int block_x = 8;
  constexpr int block_y = 16;

  constexpr int threadblock_m = 128;
  constexpr int threadblock_n = 128;
  constexpr int threadblock_k = 512; // 64 * 8

  constexpr int warp_m = 32;
  constexpr int warp_n = 32;
  constexpr int warp_k = 512; // 64 * 8

#elif defined COMET_USE_AMPERE
  // Use following for Ampere
  constexpr int block_x = 16;
  constexpr int block_y = 8;

  constexpr int threadblock_m = 128;
  constexpr int threadblock_n = 128;
  constexpr int threadblock_k = 2048; // 64 * 32

  constexpr int warp_m = 32;
  constexpr int warp_n = 32;
  constexpr int warp_k = 1024; // 64 * 16
#endif
  int grid_x = m / (threadblock_m / 2);
  int grid_y = n / (threadblock_n / 2);

  auto gemm_kernel
    = b1_comet_xor_gemm_gpu_cutlass<threadblock_m, threadblock_n, threadblock_k, warp_m, warp_n, warp_k, block_x, block_y>;

  int shared_bytes = (threadblock_m + threadblock_n) * threadblock_k / 8;
  set_max_shared_bytes((const void *)gemm_kernel);

  /*printf("Calling b1_comet_xor_gemm_gpu_cutlass kernel mnk = (%d,%d,%d) gridDim = (%d,%d,1) threadDim = (%d,%d,1) threadblock = (%d,%d,%d) warp = (%d,%d,%d) shared_bytes = %d beta=%d\n",
    m, n, k, grid_x, grid_y, block_x, block_y, threadblock_m, threadblock_n, threadblock_k,
    warp_m, warp_n, warp_k, shared_bytes, beta);*/

  gemm_kernel<<<dim3(grid_x, grid_y, 1), dim3(block_x, block_y, 1), shared_bytes>>>(
    m, n, k, (GMBits2x64 *)matA, (GMBits2x64 *)matB, beta, (GMTally2x2 *)matC);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif

#endif // _COMET_TC_SOLVE_CUTLASS_NVIDIA_I_HH_

//-----------------------------------------------------------------------------
