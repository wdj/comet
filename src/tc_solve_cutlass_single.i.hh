//-----------------------------------------------------------------------------
/*!
 * \file   tc_solve_cutlass_single.i.hh
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

#ifndef _COMET_TC_SOLVE_CUTLASS_SINGLE_I_HH_
#define _COMET_TC_SOLVE_CUTLASS_SINGLE_I_HH_

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

__device__ inline void process_bits_single(GMBits2x64 vi, uint64_t &nvi, uint64_t &pvi)
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

__device__ inline void combine_bits_single(int nn, int np, int pn, int pp, double &c0, double &c1)
{
  enum { GM_TALLY1_MAX_VALUE_BITS = 26 };

  const uint64_t r00 = nn;
  const uint64_t r01 = np;
  const uint64_t r10 = pn;
  const uint64_t r11 = pp;

  //c0 = r00 | (r10 << GM_TALLY1_MAX_VALUE_BITS);
  //c1 = r01 | (r11 << GM_TALLY1_MAX_VALUE_BITS);
  c0 = r00 | (r01 << GM_TALLY1_MAX_VALUE_BITS);
  c1 = r10 | (r11 << GM_TALLY1_MAX_VALUE_BITS);
}

//-----------------------------------------------------------------------------
/// 

__device__ inline void combine_bits_sum_single(int nn, int np, int pn, int pp, double &c0, double &c1,
                                              double &cc0, double &cc1)
{
  enum { GM_TALLY1_MAX_VALUE_BITS = 26 };

  const uint64_t r00 = nn;
  const uint64_t r01 = np;
  const uint64_t r10 = pn;
  const uint64_t r11 = pp;

  //c0 = r00 | (r10 << GM_TALLY1_MAX_VALUE_BITS);
  //c1 = r01 | (r11 << GM_TALLY1_MAX_VALUE_BITS);
  c0 = r00 | (r01 << GM_TALLY1_MAX_VALUE_BITS);
  c1 = r10 | (r11 << GM_TALLY1_MAX_VALUE_BITS);

  cc0 += c0;
  cc1 += c1;
}

//-----------------------------------------------------------------------------
///

template <int pM, int pK, int block_x, int block_y>
__device__ inline void g2r_single(GMBits2x64 *gmem, int begin, int k, uint64_t frag[])
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
      process_bits_single(vi, nvi, pvi);

      frag[y * iter_x + x] = nvi;
      frag[y * iter_x + x + iter_xy] = pvi;
    }
  }
}

//-----------------------------------------------------------------------------
/// 

template <int pM, int pK, int block_x, int block_y, bool is_a, class Layout>
__device__ inline void r2s_single(uint64_t frag[], uint64_t *smem, Layout layout)
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
__device__ inline void g2s_single(uint64_t *smem, GMBits2x64 *gmem, int begin, int m, int n, int k, Layout layout)
{
  // Block indices
  int bx = blockIdx.x, by = blockIdx.y;
  int thread_idx = threadIdx.y * blockDim.x + threadIdx.x;

  uint64_t nvi, pvi;

#pragma unroll
  for (int y = 0; y < 1; y += 100) { //block_y) { //pM
#pragma unroll
    for (int x = 0; x < 1; x += 100) { //block_x) { //pK

      int nrow = threadIdx.y + y;
      int kval = threadIdx.x + x;

      int nbrown = nrow*2;
      int nbrowp = nrow*2 + 1;
      int nbcol  = kval;
      const int patch = 64; // sizeof(uint64_t) / sizeof(cutlass::uint1b_t)

      // A row major and B col major
      int ind = begin + k * nrow + kval;

      //---Extract input values to process---
      // GMBits2x64 vi = a[aInd];
      //auto gmem_u2 = reinterpret_cast<ulonglong2 *>(gmem);
      //auto u2 = gmem_u2[ind];
      GMBits2x64 vi;

        if(kval<k*2) {
	//if(kval<k*2 && nrow<1) { //m) {
          //vi.data[0] = u2.x;
          //vi.data[1] = u2.y;

	  //process_bits_single(vi, nvi, pvi);

	  // Try setting as uint32_t values
          /*uint32_t v0[2], v1[2];
	  v0[0]=nrow+1;   v0[1]=nrow+101;
	  v1[0]=kval+201; v1[1]=kval+301;

          vi.data[0] = reinterpret_cast<uint64_t>(v0);
	  vi.data[1] = reinterpret_cast<uint64_t>(v1);*/

          //vi.data[0]=nrow+1;
	  //vi.data[1]=kval+101;
        /*if(nrow<1) { 
          vi.data[0] = UINT32_MAX;
          vi.data[1] = UINT32_MAX;
          
	  nvi = vi.data[0];
	  pvi = vi.data[1];

        } else {*/
	  vi.data[0] = 0;
	  vi.data[1] = 0;

          nvi = 0;
          pvi = 0;
        //}

        if (is_a) {
          // operand A
	  uint64_t ind1 = layout({nbrown, nbcol * patch}) / patch;
	  uint64_t ind2 = layout({nbrowp, nbcol * patch}) / patch;
          smem[ind1] = nvi;
          smem[ind2] = pvi;
          //if(kval<k && nrow<m) {
	  uint32_t* nvi2 = reinterpret_cast<uint32_t*>(&nvi);
	  uint32_t* pvi2 = reinterpret_cast<uint32_t*>(&pvi);
	  printf("b=%d,%d tidx=%d k=%d isA=%d xy=%d,%d xyblock=%d,%d pMK=%d,%d nrow=%d kval=%d patch=%d ind=%d "
		 "vi=%lu,%lu nvi=%lu=%u,%u pvi=%lu=%u,%u nbrown/rowp/col=%d,%d,%d "
		 "smemn=%d,%d=%lu=%lu smemp=%d,%d=%lu=%lu\n",
                 bx,by,thread_idx,k,is_a,x,y,block_x,block_y,pM,pK,nrow,kval,patch,ind,
		 vi.data[0],vi.data[1],nvi,nvi2[0],nvi2[1],pvi,pvi2[0],pvi2[1],nbrown,nbrowp,nbcol,
                 nbrown,nbcol*patch,ind1*patch,ind1,
                 nbrowp,nbcol*patch,ind2*patch,ind2);
	  printf("b=%d,%d tidx=%d k=%d isA=%d smemn[%lu]=%lu smemp[%lu]=%lu\n",
	         bx,by,thread_idx,k,is_a,
		 ind1,smem[ind1],ind2,smem[ind2]);
          //}
        } else {
          // operand B
	  uint64_t ind1 = layout({nbcol * patch, nbrown}) / patch;
	  uint64_t ind2 = layout({nbcol * patch, nbrowp}) / patch;
          smem[ind1] = nvi;
          smem[ind2] = pvi;
          //if(kval<k && nrow<m) {
          uint32_t* nvi2 = reinterpret_cast<uint32_t*>(&nvi);
          uint32_t* pvi2 = reinterpret_cast<uint32_t*>(&pvi);
	  printf("b=%d,%d tidx=%d k=%d isA=%d xy=%d,%d xyblock=%d,%d nrow=%d kval=%d patch=%d ind=%d vi=%lu,%lu nvi=%lu=%u,%u pvi=%lu=%u,%u "
               "nbrown/rowp/col=%d,%d,%d smemn=%d,%d=%lu=%lu smemp=%d,%d=%lu=%lu\n",
               bx,by,thread_idx,k,is_a,x,y,block_x,block_y,nrow,kval,patch,ind,vi.data[0],vi.data[1],nvi,nvi2[0],nvi2[1],pvi,pvi2[0],pvi2[1],nbrown,nbrowp,nbcol,
               nbcol*patch,nbrown,ind1*patch,ind1,
               nbcol*patch,nbrowp,ind2*patch,ind2);
          printf("b=%d,%d tidx=%d k=%d isA=%d smemn[%lu]=%lu smemp[%lu]=%lu\n",
                 bx,by,thread_idx,k,is_a,
                 ind1,smem[ind1],
	         ind2,smem[ind2]);
        //}
        }
	}
    }
  }
}

//-----------------------------------------------------------------------------
/// 

template <int bM, int bN, int bK, int wM, int wN, int wK, int block_x, int block_y>
__global__ void __launch_bounds__(block_x *block_y, 1)
b1_comet_mult_gemm_gpu_cutlass_single(int m, int n, int k, GMBits2x64 *a, GMBits2x64 *b, bool beta, int *c)
{
  /**
    bM, bN, bK - threadblock MMA shape: bK is expressed in bits
    pM, pN, pK - threadblock MMA shape: bK is expressed in uint64_t's, and pM and pN are the shapes for loading GMBits2x64

    When loading data from gmem to smem, the "n" and "p" parts are stored to separate parts of smem, and then
    each warp works on and accumulates on the "nn", "np", "pn" and "pp" results separately. At the end each
    t
    warp combines results stored in the 4 accumualtes and then write to gmem, which does not need to go through
    smem anymore.

    In other words, in each threadblock, this specific MMA problem is decomposed into 4 natural sub-problems.
    */

  // Block indices
  int bx = blockIdx.x, by = blockIdx.y;

  //if(bx==0 && by==0 && threadIdx.x==0 && threadIdx.y==0)
  //  printf("In b1_comet_mult_gemm_gpu_cutlass\n");

  int thread_idx = threadIdx.y * blockDim.x + threadIdx.x;

  int warp_id = thread_idx / 32;
  int lane_id = thread_idx % 32;

  constexpr int pM = bM;
  constexpr int pN = bN;
  constexpr int pK = bK / 64;

  // Matrix block location
  int aBegin = k * pM * bx;
  int bBegin = k * pN * by;

/*#if (__CUDA_ARCH__ < 800)
  // mma for sm75
  constexpr int alignment = 512;
  using InstructionShape = cutlass::gemm::GemmShape<8, 8, 128>;
#else
  // mma for sm80
  constexpr int alignment = 1024;
  using InstructionShape = cutlass::gemm::GemmShape<16, 8, 256>;
#endif*/
  //constexpr int alignment = 512;
  //using InstructionShape = cutlass::gemm::GemmShape<16, 8, 256>;
  using InstructionShape = cutlass::gemm::GemmShape<8, 8, 128>;

  //static_assert(wK % alignment == 0, "alignment");

  using WarpShape = cutlass::gemm::GemmShape<wM, wN, wK>;
  using ElementA = cutlass::uint1b_t;
  using ElementB = cutlass::uint1b_t;
  using ElementC = int;
  using LayoutA = cutlass::layout::RowMajor;
  //using LayoutA
  //  = cutlass::layout::RowMajorTensorOpMultiplicandCrosswise<cutlass::sizeof_bits<ElementA>::value, alignment>;
  using LayoutB = cutlass::layout::ColumnMajor;
  //using LayoutB
  //  = cutlass::layout::ColumnMajorTensorOpMultiplicandCrosswise<cutlass::sizeof_bits<ElementB>::value, alignment>;
  using LayoutC = cutlass::layout::RowMajor;

  using Mma =
    typename cutlass::gemm::warp::DefaultMmaTensorOp<WarpShape, InstructionShape, ElementA, LayoutA, ElementB,
                                                     LayoutB, ElementC, LayoutC, cutlass::arch::OpMultiplyAdd>::Type;

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
  int smem_B_start = ThreadblockShape::kM * ThreadblockShape::kK * 1 / 8;
  char *smem_ptr_B = &smem_ptr_A[smem_B_start];

  printf("b=%d,%d t=%d warp=%d lane=%d TB=%d,%d,%d TB-P=%d,%d,%d W=%d,%d,%d block=%d,%d "
	 "aBegin=%d bBegin=%d warp_count=%d,%d warp_idx=%d,%d=%d smemA=0 smemB=%d TBmnk=%d,%d,%d\n",
         bx,by,thread_idx,warp_id,lane_id,bM,bN,bK,pM,pN,pK,wM,wN,wK,block_x,block_y,
	 aBegin,bBegin,warp_count_m,warp_count_n,warp_idx_m,warp_idx_n,warp_idx_mn,
	 smem_B_start,ThreadblockShape::kM,ThreadblockShape::kN,ThreadblockShape::kK);

  uint64_t *smem_buffer_A = reinterpret_cast<uint64_t *>(smem_ptr_A);
  uint64_t *smem_buffer_B = reinterpret_cast<uint64_t *>(smem_ptr_B);

  using FragmentA = typename Mma::FragmentA;
  using FragmentB = typename Mma::FragmentB;
  using FragmentC = typename Mma::FragmentC;

  int lM = ThreadblockShape::kM, lN = ThreadblockShape::kN, lK = ThreadblockShape::kK;
  typename Mma::LayoutA layout_A = Mma::LayoutA::packed({lM, lK});
  typename Mma::LayoutB layout_B = Mma::LayoutB::packed({lK, lN});
  typename Mma::LayoutC layout_C = Mma::LayoutC::packed({m*2, n*2});
  __syncthreads();

  FragmentC accum;
  accum.clear();

  printf("b=%d,%d t=%d layoutA=%d,%d layoutB=%d,%d layoutC=%d,%d\n",
         bx,by,thread_idx,lM,lK,lK,lN,m*2,n*2);
  __syncthreads();

  //static_assert(pM % block_y == 0, "block-global-memory-loader needs to be in shape.");
  //static_assert(pN % block_y == 0, "block-global-memory-loader needs to be in shape.");
  //static_assert(pK % block_x == 0, "block-global-memory-loader needs to be in shape.");

  //constexpr int iter_x = pK / block_x;
  //constexpr int iter_y_a = pM / block_y;
  //constexpr int iter_xy_a = iter_x * iter_y_a;
  //constexpr int iter_y_b = pN / block_y;
  //constexpr int iter_xy_b = iter_x * iter_y_b;

  //uint64_t frag_a[2 * iter_xy_a];
  //uint64_t frag_b[2 * iter_xy_b];

  for (int l = 0; l < 1; l += 1) { //k; l += pK) {
    FragmentA frag_A;
    FragmentB frag_B;
    frag_A.clear();
    frag_B.clear();
    __syncthreads();

    // Here gmem -> register == "g2r" and then register -> smem == "r2s".
    // Operand A
    //g2r_single<pM, pK, block_x, block_y>(a, aBegin + l, k, frag_a);
    // Operand B
    //g2r_single<pN, pK, block_x, block_y>(b, bBegin + l, k, frag_b);
    // Operand A
    //r2s_single<pM, pK, block_x, block_y, 1>(frag_a, smem_buffer_A, layout_A);
    // Operand B
    //r2s_single<pN, pK, block_x, block_y, 0>(frag_b, smem_buffer_B, layout_B);

    //printf("b=%d,%d t=%d l=%d/%d a=%d=%d b=%d=%d\n",
    //       bx,by,thread_idx,l,k,aBegin,aBegin+l,bBegin,bBegin+l);

    g2s_single<pM, pK, block_x, block_y, 1>(smem_buffer_A, a, aBegin+l, m, n, k, layout_A);
    __syncthreads();
    g2s_single<pN, pK, block_x, block_y, 0>(smem_buffer_B, b, bBegin+l, m, n, k, layout_B);
    __syncthreads();

    // Construct warp-level matrix product
    typename Mma::IteratorA iter_A({reinterpret_cast<typename cutlass::uint1b_t *>(smem_buffer_A), layout_A}, lane_id);
    typename Mma::IteratorB iter_B({reinterpret_cast<typename cutlass::uint1b_t *>(smem_buffer_B), layout_B}, lane_id);

    iter_A.add_tile_offset({0, 0});
    iter_B.add_tile_offset({0, 0});
    //iter_A.add_tile_offset({warp_idx_m, 0});
    //iter_B.add_tile_offset({0, warp_idx_n});
    __syncthreads();

    Mma mma;
    __syncthreads();

    int nwarps = 1; //max(k*64,ThreadblockShape::kK);

//CUTLASS_PRAGMA_UNROLL
    for (int warp_k = 0; warp_k < nwarps; warp_k += Mma::Policy::MmaShape::kK) {

      printf("b=%d,%d t=%d Loading A\n",bx,by,thread_idx);
      iter_A.load(frag_A);
      iter_A.add_tile_offset({+warp_count_m, 0});
      __syncthreads();

      printf("b=%d,%d t=%d Loading B\n",bx,by,thread_idx);
      iter_B.load(frag_B);
      iter_B.add_tile_offset({0, +warp_count_n});
      __syncthreads();

      //cutlass::uint1b_t av0 = frag_A.at(0);
      //cutlass::uint1b_t* av = frag_A.data();
      uint32_t* av32 = reinterpret_cast<uint32_t*>(frag_A.data());
      uint32_t* bv32 = reinterpret_cast<uint32_t*>(frag_B.data());
      printf("b=%d,%d t=%d size=%d A=%u B=%u\n",
             bx,by,thread_idx,(int)frag_A.size(),av32[0],bv32[0]);
      __syncthreads();

      //printf("b=%d,%d t=%d l=%d/%d warp=%d:%d:%d k=%d kK=%d Aoffset=%d,0=%ld Boffset=0,%d=%ld A+shift=%d,0=%ld A-shift=%d,0=%ld B+shift=0,%d=%ld B-shift=0,%d=%ld\n",
      //       bx,by,thread_idx,l,k,warp_k,nwarps,Mma::Policy::MmaShape::kK,k,ThreadblockShape::kK,warp_idx_m,layout_A({warp_idx_m, 0}),warp_idx_n,layout_B({0,warp_idx_n}),
      //       warp_count_m,layout_A({+warp_count_m, 0}),-warp_count_m,layout_A({-warp_count_m, 0}),
      //       warp_count_n,layout_B({0, +warp_count_n}),-warp_count_n,layout_B({0, -warp_count_n}));

      /*++iter_A;
      ++iter_B;

      mma(accum, frag_A, frag_B, accum);
      __syncthreads();
      int* cv = accum.data();
      printf("b=%d,%d t=%d warp_k=%d size=%d accum=%d,%d\n",bx,by,thread_idx,warp_k,(int)accum.size(),cv[0],cv[1]);
      __syncthreads();*/
    }

    //if (l + pK < k) 
    //__syncthreads();
  }

  /*int* cv2 = accum.data();
  printf("b=%d,%d t=%d mid accum=%u,%u\n",
         bx,by,thread_idx,cv2[0],cv2[1]);
  __syncthreads();

  using output_type = int;

  using IteratorC = typename cutlass::gemm::warp::MmaTensorOpAccumulatorTileIterator<
    typename cutlass::MatrixShape<WarpShape::kM, WarpShape::kN>, output_type, LayoutC, InstructionShape,
    typename Mma::Policy::OpDelta>;

  IteratorC iter_C({reinterpret_cast<output_type *>(c), layout_C}, lane_id);

  iter_C.add_tile_offset({(pM / wM) * bx + warp_idx_m, (pN / wN) * by + warp_idx_n});
  iter_C.store(accum);
  printf("b=%d,%d tid=%d mnk=%d,%d,%d kMN=%d,%d warp_id=%d warp_count_m=%d warp_count_n=%d offset=%d,%d warpshape=%d,%d\n",
         bx,by,thread_idx,m,n,k,WarpShape::kM,WarpShape::kN,warp_id,warp_count_m,warp_count_n,
         (pM / wM) * bx + warp_idx_m, (pN / wN) * by + warp_idx_n,
         wM*2,wN*2);
  __syncthreads();

  //int cBegin = n*bx*block_x+by*block_y;
  //int ind = cBegin + 4*thread_idx;
  int* cv3 = accum.data();
  printf("b=%d,%d t=%d C[%d]=%d accum=%u,%u\n",
         bx,by,thread_idx,thread_idx,c[thread_idx],cv3[0],cv3[1]);
  __syncthreads();*/

  /*if (warp_id < warp_count_n * warp_count_m) {
    // use "double2" instead of "GMTally2x2" to make sure compiler issue 128-bit store instructions.
    using output_type = double2;

    using IteratorTallyC = typename cutlass::gemm::warp::MmaTensorOpAccumulatorTileIterator<
      typename cutlass::MatrixShape<WarpShape::kM, WarpShape::kN>, output_type, LayoutC, InstructionShape,
      typename Mma::Policy::OpDelta>;

    typename IteratorTallyC::Fragment accum_tally;
    IteratorTallyC iter_tally_C({reinterpret_cast<output_type *>(c), layout_C}, lane_id);*/

    //if(true) { //!beta) {
//#pragma unroll
      //for (int idx = 0; idx < FragmentC::kElements; idx++) {
        // Combine the results in the sub-problems.
        //combine_bits_single(accum_nn[idx], accum_np[idx], accum_pn[idx], accum_pp[idx], accum_tally[idx].x, accum_tally[idx].y);
      //}

      //iter_tally_C.add_tile_offset({(pM / wM) * bx + warp_idx_m, (pN / wN) * by + warp_idx_n});
      ////iter_tally_C.add_tile_offset({(pN / wN) * by + warp_idx_n, (pM / wM) * bx + warp_idx_m});
      // The following code does not translates into the most efficient instructions, even with 128-bit stores.
      // With current version of CUTLASS this is as far as we can go.
      //iter_tally_C.store(accum_tally);
    /*}
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
        combine_bits_sum_single(accum_nn[idx], accum_np[idx], accum_pn[idx], accum_pp[idx],
                               accum_tally[idx].x, accum_tally[idx].y, frag_C[idx].x, frag_C[idx].y);
      }

      iter_C.store(frag_C);
    }*/
  //}
}

template <int pM, int pK, bool is_a, class Layout>
__device__ inline void g2s_single2(uint64_t *smem, GMBits2x64 *gmem, int begin, int m, int n, int k, Layout layout)
{
  // Block indices
  int bx = blockIdx.x, by = blockIdx.y;
  int thread_idx = threadIdx.y * blockDim.x + threadIdx.x;

  uint64_t nvi, pvi;

  int nrow = thread_idx/(k*2);
  int kval = thread_idx%(k*2);

  // A row major and B col major
  int ind = begin + k*2 * nrow + kval;

  int nbrown = nrow*2;
  int nbrowp = nrow*2 + 1;
  int nbcol  = kval;
  const int patch = 64; // sizeof(uint64_t) / sizeof(cutlass::uint1b_t)

  //---Extract input values to process---
  auto gmem_u2 = reinterpret_cast<ulonglong2 *>(gmem);
  auto u2 = gmem_u2[ind];
  GMBits2x64 vi;

  if(nrow<m) {
    if(nrow==0) {
    vi.data[0] = u2.x;
    vi.data[1] = u2.y;
    } else {
    vi.data[0] = 0;
    vi.data[1] = 0;
    }

    //process_bits_single(vi, nvi, pvi);
    nvi = vi.data[0];
    pvi = vi.data[1];

    if (is_a) {
      // operand A
      int ind1 = thread_idx; //layout({nbrown, nbcol * patch}) / patch;
      int ind2 = thread_idx+32; //layout({nbrowp, nbcol * patch}) / patch;
      smem[ind1] = nvi;
      smem[ind2] = pvi;
      uint32_t* nvi2 = reinterpret_cast<uint32_t*>(&nvi);
      uint32_t* pvi2 = reinterpret_cast<uint32_t*>(&pvi);
      printf("b=%d,%d tidx=%d k=%d isA=%d pMK=%d,%d nrow=%d kval=%d ind=%d patch=%d "
             "nvi=%lu=%u,%u pvi=%lu=%u,%u n r/c=%d,%d p r/c=%d,%d "
             "smemn=%d,%d=%d=%d smemp=%d,%d=%d=%d "
	     "smemn[%d]=%lu smemp[%d]=%lu\n",
             bx,by,thread_idx,k,is_a,pM,pK,nrow,kval,ind,patch,
             nvi,nvi2[0],nvi2[1],pvi,pvi2[0],pvi2[1],nbrown,nbcol,nbrowp,nbcol,
             nbrown,nbcol*patch,ind1*patch,ind1,
             nbrowp,nbcol*patch,ind2*patch,ind2,
             ind1,smem[ind1],ind2,smem[ind2]);
    }
    else {
      // operand B
      int ind1 = thread_idx; //layout({nbcol * patch, nbrown}) / patch;
      int ind2 = thread_idx+32; //layout({nbcol * patch, nbrowp}) / patch;
      smem[ind1] = nvi;
      smem[ind2] = pvi;
      uint32_t* nvi2 = reinterpret_cast<uint32_t*>(&nvi);
      uint32_t* pvi2 = reinterpret_cast<uint32_t*>(&pvi);
      printf("b=%d,%d tidx=%d k=%d isA=%d pMK=%d,%d nrow=%d kval=%d ind=%d patch=%d "
             "nvi=%lu=%u,%u pvi=%lu=%u,%u n r/c==%d,%d p r/c=%d,%d "
	     "smemn=%d,%d=%d=%d smemp=%d,%d=%d=%d "
	     "smemn[%d]=%lu smemp[%d]=%lu\n",
             bx,by,thread_idx,k,is_a,pM,pK,nrow,kval,ind,patch,
	     nvi,nvi2[0],nvi2[1],pvi,pvi2[0],pvi2[1],nbrown,nbcol,nbrowp,nbcol,
             nbcol*patch,nbrown,ind1*patch,ind1,
             nbcol*patch,nbrowp,ind2*patch,ind2,
             ind1,smem[ind1],ind2,smem[ind2]);
    }
  }
}


//-----------------------------------------------------------------------------
/// 

template <int bM, int bN, int bK, int wM, int wN, int wK>
__global__ void
b1_comet_mult_gemm_gpu_cutlass_single2(int m, int n, int k, GMBits2x64 *a, GMBits2x64 *b, bool beta, int *c)
{
  /**
    bM, bN, bK - threadblock MMA shape: bK is expressed in bits
    pM, pN, pK - threadblock MMA shape: bK is expressed in uint64_t's, and pM and pN are the shapes for loading GMBits2x64

    When loading data from gmem to smem, the "n" and "p" parts are stored to separate parts of smem, and then
    each warp works on and accumulates on the "nn", "np", "pn" and "pp" results separately. At the end each
    t
    warp combines results stored in the 4 accumualtes and then write to gmem, which does not need to go through
    smem anymore.

    In other words, in each threadblock, this specific MMA problem is decomposed into 4 natural sub-problems.
    */

  // Block indices
  int bx = blockIdx.x, by = blockIdx.y;
  int thread_idx = threadIdx.y * blockDim.x + threadIdx.x;
  int warp_id = thread_idx / 32;
  int lane_id = thread_idx % 32;

  constexpr int pM = bM;
  constexpr int pN = bN;
  constexpr int pK = bK / 64;

  // Matrix block location
  int aBegin = k * pM * bx;
  int bBegin = k * pN * by;
  __syncthreads();

  using InstructionShape = cutlass::gemm::GemmShape<8, 8, 128>;

  using WarpShape = cutlass::gemm::GemmShape<wM, wN, wK>;
  using ElementA = cutlass::uint1b_t;
  using ElementB = cutlass::uint1b_t;
  using ElementC = int;
  using LayoutA = cutlass::layout::RowMajor;
  using LayoutB = cutlass::layout::ColumnMajor;
  using LayoutC = cutlass::layout::RowMajor;
  __syncthreads();

  using Mma =
    typename cutlass::gemm::warp::DefaultMmaTensorOp<WarpShape, InstructionShape, ElementA, LayoutA, ElementB,
                                                     LayoutB, ElementC, LayoutC, cutlass::arch::OpMultiplyAdd>::Type;
  using ThreadblockShape = cutlass::gemm::GemmShape<bM, bN, bK>;
  __syncthreads();

  //extern __shared__ char smem_ptr[];
  //char *smem_ptr_A = smem_ptr;
  int smem_B_start = ThreadblockShape::kM * ThreadblockShape::kK * 1 / 8;
  //char *smem_ptr_B = &smem_ptr_A[smem_B_start];


  __shared__ uint64_t smem_buffer_A[128]; 
  __shared__ uint64_t smem_buffer_B[128];
  __syncthreads();

  smem_buffer_A[lane_id] = lane_id;
  smem_buffer_B[lane_id] = lane_id;
  smem_buffer_A[lane_id+32] = lane_id+32;
  smem_buffer_B[lane_id+32] = lane_id+32;

  printf("b=%d,%d t=%d warp=%d lane=%d TB=%d,%d,%d TB-P=%d,%d,%d W=%d,%d,%d "
         "aBegin=%d bBegin=%d smemA=0 smemB=%d TBmnk=%d,%d,%d ind=%d,%d smemA=%lu,%lu smemB=%lu,%lu\n",
         bx,by,thread_idx,warp_id,lane_id,bM,bN,bK,pM,pN,pK,wM,wN,wK,
         aBegin,bBegin,smem_B_start,ThreadblockShape::kM,ThreadblockShape::kN,ThreadblockShape::kK,
	 lane_id,lane_id+32,smem_buffer_A[lane_id],smem_buffer_A[lane_id+32],smem_buffer_B[lane_id],smem_buffer_B[lane_id+32]);
  __syncthreads();

  //uint64_t *smem_buffer_A = reinterpret_cast<uint64_t *>(smem_ptr_A);
  //uint64_t *smem_buffer_B = reinterpret_cast<uint64_t *>(smem_ptr_B);

  using FragmentA = typename Mma::FragmentA;
  using FragmentB = typename Mma::FragmentB;
  using FragmentC = typename Mma::FragmentC;

  int lM = ThreadblockShape::kM, lN = ThreadblockShape::kN, lK = ThreadblockShape::kK;
  typename Mma::LayoutA layout_A = Mma::LayoutA::packed({lM, lK});
  typename Mma::LayoutB layout_B = Mma::LayoutB::packed({lK, lN});
  typename Mma::LayoutC layout_C = Mma::LayoutC::packed({m*2, n*2});
  __syncthreads();

  FragmentC accum;
  accum.clear();

  printf("b=%d,%d t=%d layoutA=%d,%d layoutB=%d,%d layoutC=%d,%d\n",
         bx,by,thread_idx,lM,lK,lK,lN,m*2,n*2);
  __syncthreads();

  // Setup arrays
    int l=0;
    FragmentA frag_A;
    FragmentB frag_B;
    __syncthreads();

    g2s_single2<pM, pK, 1>(smem_buffer_A, a, aBegin+l, m, n, k, layout_A);
    __syncthreads();
    g2s_single2<pN, pK, 0>(smem_buffer_B, b, bBegin+l, m, n, k, layout_B);
    __syncthreads();

    printf("b=%d,%d t=%d Done computing g2s\n",bx,by,thread_idx);
    __syncthreads();

    printf("b=%d,%d t=%d ind=%d,%d smemA=%lu,%lu smemB=%lu,%lu\n",
           bx,by,thread_idx,lane_id,lane_id+32,smem_buffer_A[lane_id],smem_buffer_A[lane_id+32],smem_buffer_B[lane_id],smem_buffer_B[lane_id+32]);
    __syncthreads();

    // Construct warp-level matrix product
    typename Mma::IteratorA iter_A({reinterpret_cast<typename cutlass::uint1b_t *>(smem_buffer_A), layout_A}, lane_id);
    typename Mma::IteratorB iter_B({reinterpret_cast<typename cutlass::uint1b_t *>(smem_buffer_B), layout_B}, lane_id);
    __syncthreads();

    iter_A.add_tile_offset({0, 0});
    iter_B.add_tile_offset({0, 0});

    Mma mma;
    __syncthreads();

    // Multiply two blocks
    printf("b=%d,%d t=%d Loading A\n",bx,by,thread_idx);
    iter_A.load(frag_A);
    iter_A.add_tile_offset({0, 0});
    __syncthreads();

    printf("b=%d,%d t=%d Loading B\n",bx,by,thread_idx);
    iter_B.load(frag_B);
    iter_B.add_tile_offset({0, 0});
    __syncthreads();

    uint32_t* av32 = reinterpret_cast<uint32_t*>(frag_A.data());
    uint32_t* bv32 = reinterpret_cast<uint32_t*>(frag_B.data());
    printf("b=%d,%d t=%d size=%d A=%u,%u B=%u,%u\n",
           bx,by,thread_idx,(int)frag_A.size(),av32[0],av32[1],bv32[0],bv32[1]);
    __syncthreads();

  //printf("b=%d,%d t=%d l=%d/%d warp=%d:%d:%d k=%d kK=%d Aoffset=%d,0=%ld Boffset=0,%d=%ld A+shift=%d,0=%ld A-shift=%d,0=%ld B+shift=0,%d=%ld B-shift=0,%d=%ld\n",
      //       bx,by,thread_idx,l,k,warp_k,nwarps,Mma::Policy::MmaShape::kK,k,ThreadblockShape::kK,warp_idx_m,layout_A({warp_idx_m, 0}),warp_idx_n,layout_B({0,warp_idx_n}),
      //       warp_count_m,layout_A({+warp_count_m, 0}),-warp_count_m,layout_A({-warp_count_m, 0}),
      //       warp_count_n,layout_B({0, +warp_count_n}),-warp_count_n,layout_B({0, -warp_count_n}));

      ++iter_A;
      ++iter_B;

      mma(accum, frag_A, frag_B, accum);
      __syncthreads();

      // Check
      uint64_t check = gm_popcount64(smem_buffer_A[0] & smem_buffer_B[0]);

      int* cv = accum.data();
      printf("b=%d,%d t=%d size=%d accum=%d,%d check=%lu\n",bx,by,thread_idx,(int)accum.size(),cv[0],cv[1],check);
      __syncthreads();
    
    //if (l + pK < k) 
    __syncthreads();
  

  int* cv2 = accum.data();
  printf("b=%d,%d t=%d mid accum=%u,%u\n",
         bx,by,thread_idx,cv2[0],cv2[1]);
  __syncthreads();

  using output_type = int;

  using IteratorC = typename cutlass::gemm::warp::MmaTensorOpAccumulatorTileIterator<
    typename cutlass::MatrixShape<WarpShape::kM, WarpShape::kN>, output_type, LayoutC, InstructionShape,
    typename Mma::Policy::OpDelta>;

  IteratorC iter_C({reinterpret_cast<output_type *>(c), layout_C}, lane_id);

  iter_C.add_tile_offset({0,0});
  //iter_C.add_tile_offset({(pM / wM) * bx + warp_idx_m, (pN / wN) * by + warp_idx_n});
  iter_C.store(accum);
  printf("b=%d,%d tid=%d mnk=%d,%d,%d kMN=%d,%d warp_id=%d warpshape=%d,%d\n",
         bx,by,thread_idx,m,n,k,WarpShape::kM,WarpShape::kN,warp_id,wM*2,wN*2);
  __syncthreads();

  int* cv3 = accum.data();
  printf("b=%d,%d t=%d C[%d]=%d accum=%u,%u\n",
         bx,by,thread_idx,thread_idx,c[thread_idx],cv3[0],cv3[1]);
  __syncthreads();
}

template <int pM, int pK, bool is_a, class Layout>
__device__ inline void g2s_single3(uint64_t *smem, GMBits2x64 *gmem, int begin, int m, int n, int k, Layout layout)
{
  // Block indices
  int bx = blockIdx.x, by = blockIdx.y;
  int thread_idx = threadIdx.y * blockDim.x + threadIdx.x;

  uint64_t nvi, pvi;

  int nrow = thread_idx/(k*2);
  int kval = thread_idx%(k*2);

  // A row major and B col major
  int ind = begin + k*2 * nrow + kval;

  int nbrown = nrow*2;
  int nbrowp = nrow*2 + 1;
  int nbcol  = kval;
  const int patch = 64; // sizeof(uint64_t) / sizeof(cutlass::uint1b_t)

  //---Extract input values to process---
  auto gmem_u2 = reinterpret_cast<ulonglong2 *>(gmem);
  auto u2 = gmem_u2[ind];
  GMBits2x64 vi;

  if(nrow<m) {
    if(nrow==0 && kval==0) {
      vi.data[0] = u2.x;
      vi.data[1] = u2.y;
    } else {
      vi.data[0] = 0;
      vi.data[1] = 0;
    }

    if(nrow==0 && kval==0) {
      process_bits_single(vi, nvi, pvi);
    } else {
      nvi = vi.data[0];
      pvi = vi.data[1];
    }

    //uint32_t* nvi2 = reinterpret_cast<uint32_t*>(&nvi);
    //uint32_t* pvi2 = reinterpret_cast<uint32_t*>(&pvi);

    if (is_a) {
      // operand A
      int ind1 = layout({nbrown, nbcol * patch}) / patch;
      int ind2 = layout({nbrowp, nbcol * patch}) / patch;
      //int ind1 = thread_idx;
      //int ind2 = thread_idx+32;
      smem[ind1] = nvi;
      smem[ind2] = pvi;
      printf("b=%d,%d tidx=%d k=%d isA=%d pMK=%d,%d nrow=%d kval=%d ind=%d patch=%d "
             "nvi=%lu pvi=%lu n r/c=%d,%d p r/c=%d,%d "
             "smemn=%d,%d=%d=%d smemp=%d,%d=%d=%d "
             "smemn[%d]=%lu smemp[%d]=%lu\n",
             bx,by,thread_idx,k,is_a,pM,pK,nrow,kval,ind,patch,
             nvi,pvi,nbrown,nbcol,nbrowp,nbcol,
             nbrown,nbcol*patch,ind1*patch,ind1,
             nbrowp,nbcol*patch,ind2*patch,ind2,
             ind1,smem[ind1],ind2,smem[ind2]);
    }
    else {
      // operand B
      int ind1 = layout({nbcol * patch, nbrown}) / patch;
      int ind2 = layout({nbcol * patch, nbrowp}) / patch;
      //int ind1 = thread_idx;
      //int ind2 = thread_idx+32;
      smem[ind1] = nvi;
      smem[ind2] = pvi;
      printf("b=%d,%d tidx=%d k=%d isA=%d pMK=%d,%d nrow=%d kval=%d ind=%d patch=%d "
             "nvi=%lu pvi=%lu n r/c==%d,%d p r/c=%d,%d "
             "smemn=%d,%d=%d=%d smemp=%d,%d=%d=%d "
             "smemn[%d]=%lu smemp[%d]=%lu\n",
             bx,by,thread_idx,k,is_a,pM,pK,nrow,kval,ind,patch,
             nvi,pvi,nbrown,nbcol,nbrowp,nbcol,
             nbcol*patch,nbrown,ind1*patch,ind1,
             nbcol*patch,nbrowp,ind2*patch,ind2,
             ind1,smem[ind1],ind2,smem[ind2]);
    }
  }
}


template <int bM, int bN, int bK, int wM, int wN, int wK>
__global__ void //__launch_bounds__(block_x *block_y, 1)
b1_comet_mult_gemm_gpu_cutlass_single3(int m, int n, int k, GMBits2x64 *a, GMBits2x64 *b, bool beta, int *c)
{
  /**
    bM, bN, bK - threadblock MMA shape: bK is expressed in bits
    pM, pN, pK - threadblock MMA shape: bK is expressed in uint64_t's, and pM and pN are the shapes for loading GMBits2x64

    When loading data from gmem to smem, the "n" and "p" parts are stored to separate parts of smem, and then
    each warp works on and accumulates on the "nn", "np", "pn" and "pp" results separately. At the end each
    t
    warp combines results stored in the 4 accumualtes and then write to gmem, which does not need to go through
    smem anymore.

    In other words, in each threadblock, this specific MMA problem is decomposed into 4 natural sub-problems.
    */

  // Block indices
  int bx = blockIdx.x, by = blockIdx.y;
  int thread_idx = threadIdx.y * blockDim.x + threadIdx.x;
  int warp_id = thread_idx / 32;
  int lane_id = thread_idx % 32;

  constexpr int pM = bM;
  constexpr int pN = bN;
  constexpr int pK = bK / 64;

  // Matrix block location
  int aBegin = k * pM * bx;
  int bBegin = k * pN * by;
  __syncthreads();

  using InstructionShape = cutlass::gemm::GemmShape<8, 8, 128>;

  using WarpShape = cutlass::gemm::GemmShape<wM, wN, wK>;
  using ElementA = cutlass::uint1b_t;
  using ElementB = cutlass::uint1b_t;
  using ElementC = int;
  using LayoutA = cutlass::layout::RowMajor;
  using LayoutB = cutlass::layout::ColumnMajor;
  using LayoutC = cutlass::layout::RowMajor;
  __syncthreads();

  using Mma =
    typename cutlass::gemm::warp::DefaultMmaTensorOp<WarpShape, InstructionShape, ElementA, LayoutA, ElementB,
                                                     LayoutB, ElementC, LayoutC, cutlass::arch::OpMultiplyAdd>::Type;
  using ThreadblockShape = cutlass::gemm::GemmShape<bM, bN, bK>;
  __syncthreads();

  int lM = ThreadblockShape::kM, lN = ThreadblockShape::kN, lK = ThreadblockShape::kK;
  typename Mma::LayoutA layout_A = Mma::LayoutA::packed({lM, lK});
  typename Mma::LayoutB layout_B = Mma::LayoutB::packed({lK, lN});
  typename Mma::LayoutC layout_C = Mma::LayoutC::packed({m*2, n*2});

  using FragmentA = typename Mma::FragmentA;
  using FragmentB = typename Mma::FragmentB;
  using FragmentC = typename Mma::FragmentC;
  __syncthreads();

  extern __shared__ char smem_ptr[];
  char *smem_ptr_A = smem_ptr;
  int smem_B_start = 1024; //bM*bK * 1 / 8; //ThreadblockShape::kM * ThreadblockShape::kK * 1 / 8;
  char *smem_ptr_B = &smem_ptr_A[smem_B_start];

  uint64_t *smem_buffer_A = reinterpret_cast<uint64_t *>(smem_ptr_A);
  uint64_t *smem_buffer_B = reinterpret_cast<uint64_t *>(smem_ptr_B);

  smem_buffer_A[lane_id] = 0; //lane_id;
  smem_buffer_B[lane_id] = 0; //lane_id;
  smem_buffer_A[lane_id+32] = 0; //lane_id+32;
  smem_buffer_B[lane_id+32] = 0; //lane_id+32;

  printf("b=%d,%d t=%d warp=%d lane=%d TB=%d,%d,%d TB-P=%d,%d,%d W=%d,%d,%d "
         "aBegin=%d bBegin=%d smemA=0 smemB=%d TBmnk=%d,%d,%d ind=%d,%d smemA=%lu,%lu smemB=%lu,%lu\n",
         bx,by,thread_idx,warp_id,lane_id,bM,bN,bK,pM,pN,pK,wM,wN,wK,
         aBegin,bBegin,smem_B_start,ThreadblockShape::kM,ThreadblockShape::kN,ThreadblockShape::kK,
         lane_id,lane_id+32,smem_buffer_A[lane_id],smem_buffer_A[lane_id+32],smem_buffer_B[lane_id],smem_buffer_B[lane_id+32]);
  __syncthreads();

  FragmentC accum;
  accum.clear();

  printf("b=%d,%d t=%d layoutA=%d,%d layoutB=%d,%d layoutC=%d,%d\n",
         bx,by,thread_idx,lM,lK,lK,lN,m*2,n*2);
  __syncthreads();

  // Setup arrays
    int l=0;
    FragmentA frag_A;
    FragmentB frag_B;
    __syncthreads();

  g2s_single3<pM, pK, 1>(smem_buffer_A, a, aBegin+l, m, n, k, layout_A);
  __syncthreads();
  g2s_single3<pN, pK, 0>(smem_buffer_B, b, bBegin+l, m, n, k, layout_B);
  __syncthreads();

  printf("b=%d,%d t=%d ind=%d,%d smemA=%lu,%lu smemB=%lu,%lu\n",
         bx,by,thread_idx,lane_id,lane_id+32,smem_buffer_A[lane_id],smem_buffer_A[lane_id+32],smem_buffer_B[lane_id],smem_buffer_B[lane_id+32]);
  __syncthreads();

    // Construct warp-level matrix product
    typename Mma::IteratorA iter_A({reinterpret_cast<typename cutlass::uint1b_t *>(smem_buffer_A), layout_A}, lane_id);
    typename Mma::IteratorB iter_B({reinterpret_cast<typename cutlass::uint1b_t *>(smem_buffer_B), layout_B}, lane_id);

    iter_A.add_tile_offset({1, 0});
    iter_B.add_tile_offset({0, 1});

    Mma mma;
    __syncthreads();

    // Multiply two blocks
    printf("b=%d,%d t=%d Loading A\n",bx,by,thread_idx);
    iter_A.load(frag_A);
    iter_A.add_tile_offset({1, 0});
    __syncthreads();

    printf("b=%d,%d t=%d Loading B\n",bx,by,thread_idx);
    iter_B.load(frag_B);
    iter_B.add_tile_offset({0, 1});
    __syncthreads();

    uint32_t* av32 = (uint32_t*)(frag_A.data());
    uint32_t* bv32 = (uint32_t*)(frag_B.data());
    printf("b=%d,%d t=%d size=%d A=%u,%u,%u,%u,%u,%u,%u,%u B=%u,%u,%u,%u,%u,%u,%u,%u\n",
           bx,by,thread_idx,(int)frag_A.size(),av32[0],av32[1],av32[2],av32[3],av32[4],av32[5],av32[6],av32[7],
	   bv32[0],bv32[1],bv32[2],bv32[3],bv32[4],bv32[5],bv32[6],bv32[7]);
    __syncthreads();

  //printf("b=%d,%d t=%d l=%d/%d warp=%d:%d:%d k=%d kK=%d Aoffset=%d,0=%ld Boffset=0,%d=%ld A+shift=%d,0=%ld A-shift=%d,0=%ld B+shift=0,%d=%ld B-shift=0,%d=%ld\n",
      //       bx,by,thread_idx,l,k,warp_k,nwarps,Mma::Policy::MmaShape::kK,k,ThreadblockShape::kK,warp_idx_m,layout_A({warp_idx_m, 0}),warp_idx_n,layout_B({0,warp_idx_n}),
      //       warp_count_m,layout_A({+warp_count_m, 0}),-warp_count_m,layout_A({-warp_count_m, 0}),
      //       warp_count_n,layout_B({0, +warp_count_n}),-warp_count_n,layout_B({0, -warp_count_n}));

      ++iter_A;
      ++iter_B;

      mma(accum, frag_A, frag_B, accum);
      __syncthreads();

      // Check
      int check00 = gm_popcount64(smem_buffer_A[0] & smem_buffer_B[0]);
      int check01 = gm_popcount64(smem_buffer_A[0] & smem_buffer_B[1]);
      int check10 = gm_popcount64(smem_buffer_A[1] & smem_buffer_B[0]);
      int check11 = gm_popcount64(smem_buffer_A[1] & smem_buffer_B[1]);
      int check0 = gm_popcount32(av32[0] & bv32[0]);
      int check1 = gm_popcount32(av32[1] & bv32[1]);
      int check2 = gm_popcount32(av32[2] & bv32[2]);
      int check3 = gm_popcount32(av32[3] & bv32[3]);

      int* cv = accum.data();
      printf("b=%d,%d t=%d size=%d accum=%d,%d,%d,%d,%d,%d,%d,%d check=%d,%d,%d,%d check2=%d,%d,%d,%d\n",
	     bx,by,thread_idx,(int)accum.size(),
	     cv[0],cv[1],cv[2],cv[3],cv[4],cv[5],cv[6],cv[7],check00,check01,check10,check11,check0,check1,check2,check3);
      __syncthreads();
}

//-----------------------------------------------------------------------------
/// 

void set_max_shared_bytes_single(const void *func)
{
  cudaFuncSetAttribute(func, cudaFuncAttributePreferredSharedMemoryCarveout, (int)cudaSharedmemCarveoutMaxShared);
  int max_shared_bytes;
  cudaDeviceGetAttribute(&max_shared_bytes, cudaDevAttrMaxSharedMemoryPerBlockOptin, 0);
  cudaFuncSetAttribute(func, cudaFuncAttributeMaxDynamicSharedMemorySize, max_shared_bytes);
  //printf("max_shared_bytes=%d\n",max_shared_bytes);
}

//-----------------------------------------------------------------------------
/// \brief Perform required GEMM.

void tc_solve_comet_impl_cutlass_single(int m, int n, int k, const void *matA, const void *matB, bool beta, void *matC, CEnv& env)
{
/*#if defined COMET_USE_TURING
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
#endif*/

  // Smaller test settings
  // Cutlass block parameters
  constexpr int threadblock_m = 8; // 4 blocks per thread block
  constexpr int threadblock_n = 8;
  constexpr int threadblock_k = 128; //256; // 64 * 

  constexpr int warp_m = 8; // 64 threads per block
  constexpr int warp_n = 8;
  constexpr int warp_k = 128; //256; // 64 * 8

  // Kernel launch parameters
  constexpr int threadblockx = 8;//16;
  constexpr int threadblocky = 4;//8;
  
  int gridblockx = (int)ceil((double)m/threadblockx);
  int gridblocky = (int)ceil((double)n/threadblocky);

  auto gemm_kernel
    = b1_comet_mult_gemm_gpu_cutlass_single3<threadblock_m, threadblock_n, threadblock_k, warp_m, warp_n, warp_k>; //, threadblockx, threadblocky>;

  int shared_bytes = 8192; //(threadblock_m + threadblock_n) * threadblock_k / 8;
  set_max_shared_bytes_single((const void *)gemm_kernel);

  printf("Calling b1_comet_mult_gemm_gpu_cutlass_single kernel mnk = (%d,%d,%d) gridDim = (%d,%d,1) threadDim = (%d,%d,1) threadblock = (%d,%d,%d) warp = (%d,%d,%d) shared_bytes = %d beta=%d\n",
    m, n, k, gridblockx, gridblocky, threadblockx, threadblocky, threadblock_m, threadblock_n, threadblock_k,
    warp_m, warp_n, warp_k, shared_bytes, beta);

  cudaStreamSynchronize(env.stream_compute());

  COMET_LAUNCH_KERNEL(gemm_kernel,dim3(gridblockx, gridblocky, 1), dim3(threadblockx, threadblocky, 1), shared_bytes, env.stream_compute(), m, n, k, (GMBits2x64 *)matA, (GMBits2x64 *)matB, beta, (int*) matC);

  //gemm_kernel<<<dim3(gridblockx, gridblocky, 1), dim3(threadblockx, threadblocky, 1), shared_bytes>>>(
    //m, n, k, (GMBits2x64 *)matA, (GMBits2x64 *)matB, beta, (int*) matC); // (GMTally2x2 *)matC);
  cudaStreamSynchronize(env.stream_compute());
  printf("Done with b1_comet_mult_gemm_gpu_cutlass_single kernel\n");
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TC_SOLVE_CUTLASS_SINGLE_I_HH_

//-----------------------------------------------------------------------------
