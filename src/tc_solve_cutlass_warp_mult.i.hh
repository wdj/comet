//-----------------------------------------------------------------------------
/*!
 * \file   tc_solve_cutlass_nvidia.i.hh
 * \author Paul Eller
 * \date   Wed Nov  11 16:44:20 EST 2020
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

#ifndef _COMET_TC_SOLVE_CUTLASS_WARP_I_HH_
#define _COMET_TC_SOLVE_CUTLASS_WARP_I_HH_

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

__device__ inline void process_bits_int(GMBits2x64 vi, uint64_t &nvi, uint64_t &pvi)
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

__device__ inline void combine_bits_int(int nn, int np, int pn, int pp, double &c0, double &c1)
{
  enum { GM_TALLY1_MAX_VALUE_BITS = 26 };

  const uint64_t r00 = nn;
  const uint64_t r01 = np;
  const uint64_t r10 = pn;
  const uint64_t r11 = pp;

  c0 = r00 | (r01 << GM_TALLY1_MAX_VALUE_BITS);
  c1 = r10 | (r11 << GM_TALLY1_MAX_VALUE_BITS);
}

//-----------------------------------------------------------------------------
///

template <int pM, int pK, int block_x, int block_y>
__device__ inline void g2r_int(GMBits2x64 *gmem, int begin, int k, uint64_t frag[])
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
      process_bits_int(vi, nvi, pvi);

      frag[y * iter_x + x] = nvi;
      frag[y * iter_x + x + iter_xy] = pvi;
    }
  }
}

//-----------------------------------------------------------------------------
/// 

template <int pM, int pK, int block_x, int block_y, bool is_a, class Layout>
__device__ inline void r2s_int(uint64_t frag[], uint64_t *smem, Layout layout)
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
__device__ inline void g2s_int(uint64_t *smem, GMBits2x64 *gmem, int begin, int k, Layout layout)
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

      int nbrow = nrow*2;
      int nbcol = kval;

      printf("nrow=%d kval=%d ind=%d nbrow=%d nbcol=%d x=%d xr=0:blx=%d:pK=%d y=%d yr=0:bly=%d:pM=%d\n",
             nrow,kval,ind,nbrow,nbcol,x,block_x,pK,y,block_y,pM);

      const int patch = 64; // sizeof(uint64_t) / sizeof(cutlass::uint1b_t)

      uint64_t nvi;
      uint64_t pvi;
      process_bits_int(vi, nvi, pvi);

      if (is_a) {
        // operand A
        smem[layout({nbrow, nbcol * patch}) / patch] = nvi;
        smem[layout({nbrow + 1, nbcol * patch}) / patch] = pvi;
        printf("A smem[%d]=%lu [%d]=%lu\n",(int)layout({nbrow, nbcol * patch}) / patch,nvi,
               (int)layout({nbrow + 1, nbcol * patch}) / patch,pvi);
      } else {
        // operand B
        smem[layout({nbcol * patch, nbrow}) / patch] = nvi;
        smem[layout({nbcol * patch, nbrow + 1}) / patch] = pvi;
        printf("B smem[%d]=%lu [%d]=%lu\n",(int)layout({nbcol * patch, nbrow}) / patch,nvi,
               (int)layout({nbcol * patch, nbrow + 1}) / patch,pvi);
      }
    }
  }
}

//-----------------------------------------------------------------------------
/// 

template <int bM, int bN, int bK, int wM, int wN, int wK, int block_x, int block_y>
__global__ void __launch_bounds__(block_x *block_y, 1)
b1_comet_xor_gemm_gpu_cutlass_int(int m, int n, int k, GMBits2x64 *a, GMBits2x64 *b, int *c)
{
  /**
    bM, bN, bK - threadblock MMA shape: bK is expressed in bits
    pM, pN, pK - threadblock MMA shape: bK is expressed in uint64_t's, and pM and pN are the shapes for loading GMBits2x64
    */

  // Block indices
  int bx = blockIdx.x, by = blockIdx.y;

  int thread_idx = threadIdx.y * blockDim.x + threadIdx.x;

  int warp_id = thread_idx / 32;
  int lane_id = thread_idx % 32;

  constexpr int pM = bM / 2;
  constexpr int pN = bN / 2;
  constexpr int pK = bK / 64;

  printf("b=%d,%d tid=%d warp=%d lane=%d block_x=%d block_y=%d bMNK=%d,%d,%d pMNK=%d,%d,%d wMNK=%d,%d,%d\n",
         bx,by,thread_idx,warp_id,lane_id,block_x,block_y,bM,bN,bK,pM,pN,pK,wM,wN,wK);

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

  //static_assert(wK % alignment == 0, "alignment");

  using WarpShape = cutlass::gemm::GemmShape<wM*2, wN*2, wK>;
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
  constexpr int warp_count_m = pM * 2 / WarpShape::kM;
  constexpr int warp_count_n = pN * 2 / WarpShape::kN;

  int warp_idx_mn = warp_id % (warp_count_m * warp_count_n);
  int warp_idx_m = warp_idx_mn % warp_count_m;
  int warp_idx_n = warp_idx_mn / warp_count_m;

  extern __shared__ char smem_ptr[];

  //char *smem_ptr_A = smem_ptr;
  //char *smem_ptr_B = &smem_ptr_A[ThreadblockShape::kM*2 * ThreadblockShape::kK * 1 / 8];

  //uint64_t *smem_buffer_A = reinterpret_cast<uint64_t *>(smem_ptr_A);
  //uint64_t *smem_buffer_B = reinterpret_cast<uint64_t *>(smem_ptr_B);

  using FragmentA = typename Mma::FragmentA;
  using FragmentB = typename Mma::FragmentB;
  using FragmentC = typename Mma::FragmentC;

  typename Mma::LayoutA layout_A = Mma::LayoutA::packed({ThreadblockShape::kM*2, ThreadblockShape::kK});
  typename Mma::LayoutB layout_B = Mma::LayoutB::packed({ThreadblockShape::kK,   ThreadblockShape::kN*2});
  typename Mma::LayoutC layout_C = Mma::LayoutC::packed({m*2, n*2});

  FragmentC accum;
  accum.clear();

  static_assert(pM % block_y == 0, "block-global-memory-loader needs to be in shape.");
  static_assert(pN % block_y == 0, "block-global-memory-loader needs to be in shape.");
  static_assert(pK % block_x == 0, "block-global-memory-loader needs to be in shape.");

  printf("b=%d,%d tid=%d Asize=%dx%d Astart=0 aBegin=%d Bsize=%dx%d Bstart=%d bBegin=%d ABTotal=%d "
         "warpshape=%d,%d warp_id=%d warp_count=%d,%d warp_idx=%d,%d=%d\n",
         bx,by,thread_idx,ThreadblockShape::kM*2,ThreadblockShape::kK,aBegin,
         ThreadblockShape::kK,ThreadblockShape::kN*2,ThreadblockShape::kM*2 * ThreadblockShape::kK * 1 / 8,bBegin,
         (ThreadblockShape::kM*2 + ThreadblockShape::kN) * ThreadblockShape::kK / 8,
         wM*2,wN*2,warp_id,warp_count_m,warp_count_n,warp_idx_m,warp_idx_n,warp_idx_mn);

  /*for (int l = 0; l < k; l += pK) {
    // Operand A
    g2s_int<pM,pK,block_x,block_y,1>(smem_buffer_A, a, aBegin+l, k, layout_A);
    // Operand B
    g2s_int<pN,pK,block_x,block_y,0>(smem_buffer_B, b, bBegin+l, k, layout_B);

    __syncthreads();

    //
    // Construct warp-level matrix product
    //

    typename Mma::IteratorA iter_A({reinterpret_cast<typename cutlass::uint1b_t *>(smem_buffer_A), layout_A}, lane_id);
    typename Mma::IteratorB iter_B({reinterpret_cast<typename cutlass::uint1b_t *>(smem_buffer_B), layout_B}, lane_id);

    iter_A.add_tile_offset({warp_idx_m, 0});
    iter_B.add_tile_offset({0, warp_idx_n});

    FragmentA frag_A;
    FragmentB frag_B;

    if (warp_id < warp_count_n * warp_count_m) {

      Mma mma;

      CUTLASS_PRAGMA_UNROLL
      for (int warp_k = 0; warp_k < ThreadblockShape::kK; warp_k += Mma::Policy::MmaShape::kK) {
        iter_A.load(frag_A);
        iter_B.load(frag_B);

        ++iter_A;
        ++iter_B;

        mma(accum, frag_A, frag_B, accum);
      }
    }

    if (l + pK < k) { __syncthreads(); }
  }

  if (warp_id < warp_count_n * warp_count_m) {
    using output_type = int;

    using IteratorC = typename cutlass::gemm::warp::MmaTensorOpAccumulatorTileIterator<
      typename cutlass::MatrixShape<WarpShape::kM, WarpShape::kN>, output_type, LayoutC, InstructionShape,
      typename Mma::Policy::OpDelta>;

    IteratorC iter_C({reinterpret_cast<output_type *>(c), layout_C}, lane_id);

    iter_C.add_tile_offset({(pM / wM) * bx + warp_idx_m, (pN / wN) * by + warp_idx_n});
    iter_C.store(accum);
    printf("b=%d,%d tid=%d mnk=%d,%d,%d warp_id=%d warp_count_m=%d warp_count_n=%d offset=%d,%d warpshape=%d,%d c=%d\n",
           bx,by,thread_idx,m,n,k,warp_id,warp_count_m,warp_count_n,
           (pM / wM) * bx + warp_idx_m, (pN / wN) * by + warp_idx_n,
           wM*2,wN*2,c[0]);
  } else {
    printf("b=%d,%d tid=%d warp_id=%d warp_count_m=%d warp_count_n=%d\n",
           bx,by,thread_idx,warp_id,warp_count_m,warp_count_n);
  }
  __syncthreads();
  int cStart = bx*blockDim.x*n + blockDim.y*by;
  printf("b=%d,%d t=%d,%d bD=%d,%d c[%d]=%d\n",bx,by,
         threadIdx.x,threadIdx.y,blockDim.x,blockDim.y,
         cStart+threadIdx.x*n+threadIdx.y,c[cStart+threadIdx.x*n+threadIdx.y]);*/
}

//-----------------------------------------------------------------------------
/// 

void set_max_shared_bytes_int(const void *func)
{
  cudaFuncSetAttribute(func, cudaFuncAttributePreferredSharedMemoryCarveout, (int)cudaSharedmemCarveoutMaxShared);
  int max_shared_bytes;
  cudaDeviceGetAttribute(&max_shared_bytes, cudaDevAttrMaxSharedMemoryPerBlockOptin, 0);
  cudaFuncSetAttribute(func, cudaFuncAttributeMaxDynamicSharedMemorySize, max_shared_bytes);
  //printf("max_shared_bytes=%d\n",max_shared_bytes);
}

//-----------------------------------------------------------------------------
/// \brief Perform required GEMM.

void tc_solve_comet_impl_cutlass_int(int m, int n, int k, const void *matA, const void *matB, void *matC)
{
#if 1
  // Use following for Turing
  /*constexpr int block_x = 8;
  constexpr int block_y = 16;

  constexpr int threadblock_m = 128;
  constexpr int threadblock_n = 128;
  constexpr int threadblock_k = 512; // 64 * 8

  constexpr int warp_m = 32;
  constexpr int warp_n = 32;
  constexpr int warp_k = 512;*/ // 64 * 8
#else
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

  // Smaller test settings
  constexpr int block_x = 8;
  constexpr int block_y = 8;

  constexpr int threadblock_m = 16; // 4 blocks per thread block
  constexpr int threadblock_n = 16;
  constexpr int threadblock_k = 512; // 64 * 8

  constexpr int warp_m = 8; // 64 threads per block
  constexpr int warp_n = 8;
  constexpr int warp_k = 512; // 64 * 8

  // Smallest test settings
  constexpr int block_x = 8;
  constexpr int block_y = 8;

  constexpr int threadblock_m = 8; // 1 block per thread block
  constexpr int threadblock_n = 8;
  constexpr int threadblock_k = 512; // 64 * 8

  constexpr int warp_m = 8; // 64 threads per block
  constexpr int warp_n = 8;
  constexpr int warp_k = 512; // 64 * 8

  int grid_x = m / (threadblock_m / 2);
  int grid_y = n / (threadblock_n / 2);

  auto gemm_kernel
    = b1_comet_xor_gemm_gpu_cutlass_int<threadblock_m, threadblock_n, threadblock_k, warp_m, warp_n, warp_k, block_x, block_y>;

  int shared_bytes = (threadblock_m + threadblock_n) * threadblock_k / 8;
  set_max_shared_bytes_int((const void *)gemm_kernel);

  printf("Calling b1_comet_xor_gemm_gpu_cutlass_int kernel gridDim = (%d,%d,1) threadDim = (%d,%d,1) threadblock = (%d,%d,%d) warp = (%d,%d,%d) shared_bytes = %d\n",
    grid_x, grid_y, block_x, block_y, threadblock_m, threadblock_n, threadblock_k,
    warp_m, warp_n, warp_k, shared_bytes);

  gemm_kernel<<<dim3(grid_x, grid_y, 1), dim3(block_x, block_y, 1), shared_bytes>>>(
    m, n, k, (GMBits2x64 *)matA, (GMBits2x64 *)matB, (int*)matC);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif

//-----------------------------------------------------------------------------
