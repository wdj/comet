//-----------------------------------------------------------------------------
/*!
 * \file   tc_solve_comet_xor.i.hh
 * \author Paul Eller
 * \date   Tue Nov  3 08:26:29 EST 2020
 * \brief  CUDA code, gemm operation, native CoMet fmt, tc.
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

#ifndef _COMET_TC_SOLVE_XOR_I_HH_
#define _COMET_TC_SOLVE_XOR_I_HH_

// Includes
#include "cstdlib"
//#include <stdlib.h>
#include <cuda_runtime.h>
#include "cuda_fp16.h"
#include "types.hh"

#include "tc_solve_cutlass_nvidia.i.hh"
#include "tc_solve_cutlass_warp.i.hh"

// GPU bit count routine
#define gm_popcount64(x) __popcll(x)

// Cuda block size
#define BLOCK_SIZE 8

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief GPU kernel for custom 1-bit WMMA GEMM with bitcount and 
/// int output array

/*__global__ void b1_comet_gemm_gpu_int_simple(int m, int n, int k, GMBits2x64* a,
                                               GMBits2x64* b, bool beta, int32_t* c) {
  // Block and thread indices
  int tx = threadIdx.x, ty = threadIdx.y;
  int bx = blockIdx.x, by = blockIdx.y;

  int gridx = bx*BLOCK_SIZE + tx;
  int gridy = by*BLOCK_SIZE + ty;

  //if(bx==0 && by==0 && tx==0 && ty==0) printf("In b1_comet_gemm_gpu_int_simple mnk=%d,%d,%d a=%dx%d * b=%dx%d = c=%dx%d gxy=%d,%d\n",m,n,k,m,k,k,n,m,n,gridx,gridy);

  if(gridx>=m || gridy>=n) return;

  enum { GM_TALLY1_MAX_VALUE_BITS = 26 };

  // Matrix block location
  int aBegin = k * BLOCK_SIZE * bx;
  int bBegin = k * BLOCK_SIZE * by;

  // Stores element of block sub-matrix computed by thread
  int32_t c0 = 0, c1 = 0, c2 = 0, c3 = 0;

  // Each thread computes one element of block sub-matrix
  for (int l=0; l<k; ++l) {

    // A row major and B col major
    int aInd = aBegin + k*tx + l;
    int bInd = bBegin + k*ty + l;

    //---Extract input values to process---
    const GMBits2x64 vi = a[aInd];
    const GMBits2x64 vj = b[bInd];

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
    const uint64_t vj0 = vj.data[0];
    const uint64_t vj1 = vj.data[1];

    //if(tx==0 && ty==0)
    //printf("b=%d,%d t=%d,%d g=%d,%d a=%d=%d b=%d=%d mnk=%d,%d,%d vi0=%lu vi1=%lu vj0=%lu vj1=%lu\n",
    //       bx,by,tx,ty,gridDim.x,gridDim.y,aBegin,aInd,bBegin,bInd,m,n,k,vi0,vi1,vj0,vj1);

    //---Get mask to ignore vi seminibbles with value of 1,0---
    //---NOTE: check that this handles pad properly---
    const uint64_t oddbits = 0x5555555555555555;
    const uint64_t vij0_10mask = (vi0 | ~(vi0 >> 1)) &
                                  (vj0 | ~(vj0 >> 1)) & oddbits;
    const uint64_t vij1_10mask = (vi1 | ~(vi1 >> 1)) &
                                  (vj1 | ~(vj1 >> 1)) & oddbits;

    //---Get even, odd bits for each semi-nibble, then mask---
    const uint64_t vi0_0 =    vi0 & vij0_10mask;
    const uint64_t vi1_0 =    vi1 & vij1_10mask;
    const uint64_t vj0_0 =    vj0 & vij0_10mask;
    const uint64_t vj1_0 =    vj1 & vij1_10mask;

    //---Get complements of the same bits, then mask---
    const uint64_t nvi0_0 = ~ vi0 & vij0_10mask;
    const uint64_t nvi1_0 = ~ vi1 & vij1_10mask;
    const uint64_t nvj0_0 = ~ vj0 & vij0_10mask;
    const uint64_t nvj1_0 = ~ vj1 & vij1_10mask;

    const uint64_t r00 = gm_popcount64((nvi0_0 & nvj0_0) |
                                      ( (nvi1_0 & nvj1_0) << 1 ));
    const uint64_t r01 = gm_popcount64((nvi0_0 &  vj0_0) |
                                      ( (nvi1_0 &  vj1_0) << 1 ));
    const uint64_t r10 = gm_popcount64(( vi0_0 & nvj0_0) |
                                      ( ( vi1_0 & nvj1_0) << 1 ));
    const uint64_t r11 = gm_popcount64(( vi0_0 &  vj0_0) |
                                      ( ( vi1_0 &  vj1_0) << 1 ));

    //---Accumulate---
    c0 += r00;
    c1 += r01;
    c2 += r10;
    c3 += r11;

    //printf("b=%d,%d t=%d,%d a=%d b=%d r00=%ld r01=%ld r10=%ld r11=%ld c0=%d c1=%d c2=%d c3=%d\n",
    //       bx,by,tx,ty,aInd,bInd,r00,r01,r10,r11,c0,c1,c2,c3);
  }

  // Each thread writes one element of block sub-matrix to memory
  // Assume c is row major
  int cBegin = n*bx*BLOCK_SIZE*4+by*BLOCK_SIZE*4;
  int cInd   = cBegin + tx*n*4 + ty*4;
  //printf("b=%d,%d t=%d,%d c=%d c0123=%d,%d,%d,%d\n",bx,by,tx,ty,cInd,c0,c1,c2,c3);
  c[cInd]   = c0;
  c[cInd+1] = c1;
  c[cInd+2] = c2;
  c[cInd+3] = c3;
}*/

//-----------------------------------------------------------------------------
/// \brief GPU kernel for simple 1-bit tensor core CoMet int GEMM

__global__ 
void b1_comet_xor_gemm_gpu_int_simple(int m, int n, int k,
  GMBits2x64* a, GMBits2x64* b, bool beta, int32_t* c) {

  // Block and thread indices
  int tx = threadIdx.x, ty = threadIdx.y;
  int bx = blockIdx.x, by = blockIdx.y;

  int gridx = bx*BLOCK_SIZE + tx;
  int gridy = by*BLOCK_SIZE + ty;

  //if(bx==0 && by==0 && tx==0 && ty==0) printf("In b1_comet_xor_gemm_gpu_int_simple mnk=%d,%d,%d a=%dx%d * b=%dx%d = c=%dx%d gxy=%d,%d\n",m,n,k,m,k,k,n,m,n,gridx,gridy);

  if(gridx>=m || gridy>=n) return;

  enum { GM_TALLY1_MAX_VALUE_BITS = 26 };

  // Matrix block location
  int aBegin = k * BLOCK_SIZE * bx;
  int bBegin = k * BLOCK_SIZE * by;

  // Stores element of block sub-matrix computed by thread
  int32_t c0 = 0, c1 = 0, c2 = 0, c3 = 0;

  // Each thread computes one element of block sub-matrix
  for (int l=0; l<k; ++l) {

    // A row major and B col major
    int aInd = aBegin + k*tx + l;
    int bInd = bBegin + k*ty + l;

    //---Extract input values to process---
    const GMBits2x64 vi = a[aInd];
    const GMBits2x64 vj = b[bInd];

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
    const uint64_t vj0 = vj.data[0];
    const uint64_t vj1 = vj.data[1];

    //if(tx==0 && ty==0)
    //printf("b=%d,%d t=%d,%d g=%d,%d a=%d=%d b=%d=%d mnk=%d,%d,%d vi0=%lu vi1=%lu vj0=%lu vj1=%lu\n",
    //       bx,by,tx,ty,gridDim.x,gridDim.y,aBegin,aInd,bBegin,bInd,m,n,k,vi0,vi1,vj0,vj1);

    // Compute masks to sample the single needed bit from each seminibble,
    // and to ignore undefined vector entries.

    const uint64_t oddbits = 0x5555555555555555;
    const uint64_t vi0mask = (vi0 | ~(vi0 >> 1)) & oddbits;
    const uint64_t vi1mask = (vi1 | ~(vi1 >> 1)) & oddbits;
    const uint64_t vj0mask = (vj0 | ~(vj0 >> 1)) & oddbits;
    const uint64_t vj1mask = (vj1 | ~(vj1 >> 1)) & oddbits;

    // Extract elts that are a "1" bit (=01).

    const uint64_t pvi0 =  vi0  & vi0mask;
    const uint64_t pvi1 =  vi1  & vi1mask;
    const uint64_t pvj0 =  vj0  & vj0mask;
    const uint64_t pvj1 =  vj1  & vj1mask;

    // Extract elts that are an "0" bit (=00).

    const uint64_t nvi0 = ~vi0  & vi0mask;
    const uint64_t nvi1 = ~vi1  & vi1mask;
    const uint64_t nvj0 = ~vj0  & vj0mask;
    const uint64_t nvj1 = ~vj1  & vj1mask;

    // Combine lower, upper words - each only uses odd bits - make packed.

    const uint64_t pvi = pvi0 | (pvi1 << 1);
    const uint64_t pvj = pvj0 | (pvj1 << 1);
    const uint64_t nvi = nvi0 | (nvi1 << 1);
    const uint64_t nvj = nvj0 | (nvj1 << 1);

    const uint64_t r00 = gm_popcount64(nvi ^ nvj);
    const uint64_t r01 = gm_popcount64(nvi ^ pvj);
    const uint64_t r10 = gm_popcount64(pvi ^ nvj);
    const uint64_t r11 = gm_popcount64(pvi ^ pvj);

    //---Accumulate---                
    c0 += r00; c1 += r01;
    c2 += r10; c3 += r11;
    
    //printf("b=%d,%d t=%d,%d a=%d b=%d r00=%ld r01=%ld r10=%ld r11=%ld c0=%d c1=%d c2=%d c3=%d\n",
    //       bx,by,tx,ty,aInd,bInd,r00,r01,r10,r11,c0,c1,c2,c3);
  }        

  // Different ordering from original 1-bit routines
  // Each thread writes one element of block sub-matrix to memory
  // Assume c is row major
  int cBegin = n*bx*BLOCK_SIZE*4+by*BLOCK_SIZE*4;
  int cInd   = cBegin + tx*n*4 + ty*4;
  //printf("b=%d,%d t=%d,%d c=%d c0123=%d,%d,%d,%d\n",bx,by,tx,ty,cInd,c0,c1,c2,c3);
  c[cInd]   = c0; c[cInd+1] = c1;
  c[cInd+2] = c2; c[cInd+3] = c3;
  
  // Ordering matching original 1-bit routines
  /*int cBegin = n*bx*BLOCK_SIZE*4+by*BLOCK_SIZE*4;
  int cInd   = cBegin + tx*n + ty;
  c[cInd]   = c0;
  c[cInd+1] = c1;*/
}

//-----------------------------------------------------------------------------
/// \brief Perform required GEMM.

template<int TC_METHOD>
static void tc_solve_comet_int_impl(bool is_first, int m, int n, int k,
  const void *matA, const void *matB, void* matC, TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST(matC);
  COMET_INSIST(m >= 0 && n >= 0 && k >= 0);

  if(env.print_details()) printf("In tc_solve_comet_int_impl mnk=%d,%d,%d\n",m,n,k);
  double tbegin = env.synced_time();

  const bool beta = 1;
  const int threadblockx = BLOCK_SIZE, threadblocky = BLOCK_SIZE;
  int gridblockx = (int)ceil((double)m/threadblockx);
  int gridblocky = (int)ceil((double)n/threadblocky);

  if(env.print_details())
    printf("Launching 1-bit GEMM kernel mnk=%d,%d,%d beta=%d "
           "gridDim=%d,%d threadDim=%d,%d\n",
           m,n,k,(int)beta,gridblockx,gridblocky,threadblockx,threadblocky);

  switch(env.num_kernel()) {
    // Basic GEMM
    case 20: {
      if(env.print_details()) printf("Calling b1_comet_xor_gemm_gpu_int_simple kernel\n");
      COMET_LAUNCH_KERNEL(b1_comet_xor_gemm_gpu_int_simple,
        dim3(gridblockx, gridblocky, 1),
        dim3(threadblockx, threadblocky, 1), 0, env.stream_compute(),
        m, n, k, (GMBits2x64*)matA,
        (GMBits2x64*)matB, beta, (int32_t*)matC);
    } break;

    // Optimized Cutlass GEMM that outputs ints
    case 25: {
      if(env.print_details()) printf("Calling tc_solve_comet_impl_cutlass_int\n");     
      tc_solve_comet_impl_cutlass_int(m,n,k,(GMBits2x64*)matA,
        (GMBits2x64*)matB, (int32_t*)matC);
    } break;
    default: {
      printf("Failed to call appropriate 1-bit GEMM kernel for num_kernel=%d\n",
         env.num_kernel());
      COMET_INSIST(false && "Failure to call GEMM function.");
    }
  }
  System::accel_last_call_succeeded();
  env.ops_local_inc(2 * m * (double)n * (double)k);

  env.gemmtime_inc(env.synced_time() - tbegin);
}

//-----------------------------------------------------------------------------
/// \brief Call to perform required GEMM.

template<int TC_METHOD>
void tc_solve_comet_int_(bool is_first, int nvll, int nvl, int npvfl_thisstep,
                         const void *matA, const void *matB,void* matC, TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST(matC);
  COMET_INSIST(nvll >= 0 && nvl >= 0 && nvll <= nvl);
  COMET_INSIST(npvfl_thisstep >= 0);
  COMET_INSIST(env.tc_eff() != TC::NO);

  const int nfl_thisstep = npvfl_thisstep;

  const int m = nvll; // metrics array dim
  const int n = nvl; // metrics array dim
  const int k = nfl_thisstep; // vectors array (as GemmIn_t) dim

  if(env.print_details()) printf("Calling tc_solve_comet_int_impl with mnk=%d,%d,%d\n",m,n,k);
  tc_solve_comet_int_impl<TC_METHOD>(is_first, m, n, k, matA, matB, matC, tc_bufs, env);
}

//-----------------------------------------------------------------------------
/// \brief GPU kernel for simple 1-bit xor CoMet GEMM kernel

__global__ void b1_comet_xor_gemm_gpu_simple(int m, int n, int k,
  GMBits2x64* a, GMBits2x64* b, bool beta, GMTally2x2* c) {

  // Block and thread indices
  int tx = threadIdx.x, ty = threadIdx.y;
  int bx = blockIdx.x, by = blockIdx.y;

  int gridx = bx*BLOCK_SIZE + tx;
  int gridy = by*BLOCK_SIZE + ty;

  //if(bx==0 && by==0 && tx==0 && ty==0)
  //  printf("In b1_comet_xor_gemm_gpu_simple mnk=%d,%d,%d a=%dx%d * b=%dx%d = c=%dx%d gxy=%d,%d\n",
  //         m,n,k,m,k,k,n,m,n,gridx,gridy);

  if(gridx>=m || gridy>=n) return;

  enum { GM_TALLY1_MAX_VALUE_BITS = 26 };

  // Matrix block location
  int aBegin = k * BLOCK_SIZE * bx;
  int bBegin = k * BLOCK_SIZE * by;

  // Stores element of block sub-matrix computed by thread
  double c0 = 0, c1 = 0;
  int32_t ci0 = 0, ci1 = 0, ci2 = 0, ci3 = 0;

  // Each thread computes one element of block sub-matrix
  for (int l=0; l<k; ++l) {

    // A row major and B col major
    int aInd = aBegin + k*tx + l;
    int bInd = bBegin + k*ty + l;

    //---Extract input values to process---
    const GMBits2x64 vi = a[aInd];
    const GMBits2x64 vj = b[bInd];
    //printf("b=%d,%d t=%d,%d k=%d l=%d a=%d=%d b=%d=%d\n",bx,by,tx,ty,k,l,aBegin,aInd,bBegin,bInd);

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
    const uint64_t vj0 = vj.data[0];
    const uint64_t vj1 = vj.data[1];

    // Compute masks to sample the single needed bit from each seminibble,
    // and to ignore undefined vector entries.

    const uint64_t oddbits = 0x5555555555555555;
    const uint64_t vi0mask = (vi0 | ~(vi0 >> 1)) & oddbits;
    const uint64_t vi1mask = (vi1 | ~(vi1 >> 1)) & oddbits;
    const uint64_t vj0mask = (vj0 | ~(vj0 >> 1)) & oddbits;
    const uint64_t vj1mask = (vj1 | ~(vj1 >> 1)) & oddbits;

    // Extract elts that are a "1" bit (=01).

    const uint64_t pvi0 =  vi0  & vi0mask;
    const uint64_t pvi1 =  vi1  & vi1mask;
    const uint64_t pvj0 =  vj0  & vj0mask;
    const uint64_t pvj1 =  vj1  & vj1mask;

    // Extract elts that are an "0" bit (=00).

    const uint64_t nvi0 = ~vi0  & vi0mask;
    const uint64_t nvi1 = ~vi1  & vi1mask;
    const uint64_t nvj0 = ~vj0  & vj0mask;
    const uint64_t nvj1 = ~vj1  & vj1mask;

    // Combine lower, upper words - each only uses odd bits - make packed.

    const uint64_t pvi = pvi0 | (pvi1 << 1);
    const uint64_t pvj = pvj0 | (pvj1 << 1);
    const uint64_t nvi = nvi0 | (nvi1 << 1);
    const uint64_t nvj = nvj0 | (nvj1 << 1);

    const uint64_t r00 = gm_popcount64(nvi ^ nvj);
    const uint64_t r01 = gm_popcount64(nvi ^ pvj);
    const uint64_t r10 = gm_popcount64(pvi ^ nvj);
    const uint64_t r11 = gm_popcount64(pvi ^ pvj);

    //---Accumulate---              
    c0 += r00 | (r01 << GM_TALLY1_MAX_VALUE_BITS);
    c1 += r10 | (r11 << GM_TALLY1_MAX_VALUE_BITS);  
    ci0 += r00; ci1 += r01; ci2 += r10; ci3 += r11;
 
    //printf("b=%d,%d t=%d,%d r=%ld %ld %ld %ld c=%d %d %d %d\n",
    //       bx,by,tx,ty,r00,r01,r10,r11,ci0,ci1,ci2,ci3);
    //printf("b=%d,%d t=%d,%d a=%d b=%d r00=%ld r01=%ld r10=%ld r11=%ld c0=%lf c1=%lf\n",
    //       bx,by,tx,ty,aInd,bInd,r00,r01,r10,r11,c0,c1);
  }

  // Each thread writes one element of block sub-matrix to memory
  int cBegin = n*bx*BLOCK_SIZE+by*BLOCK_SIZE;
  int cInd   = cBegin + tx*n + ty;
  c[cInd].data[0] = c0;
  c[cInd].data[1] = c1;

  printf("b=%d,%d t=%d,%d c=%d %d %d %d cInd=%d c0=%f c1=%f\n",
         bx,by,tx,ty,ci0,ci1,ci2,ci3,cInd,c0,c1);


  // Each thread writes one element of block sub-matrix to memory
  // Assume c is row major
  /*int cBegin = n*bx*BLOCK_SIZE*4+by*BLOCK_SIZE*4;
  int cInd   = cBegin + tx*n*4 + ty*4;
  printf("b=%d,%d t=%d,%d c=%d c0123=%d,%d,%d,%d\n",bx,by,tx,ty,cInd,c0,c1,c2,c3);
  c[cInd]   = c0;
  c[cInd+1] = c1;
  c[cInd+2] = c2;
  c[cInd+3] = c3;*/
}

//-----------------------------------------------------------------------------
/// \brief GPU kernel for custom 1-bit tensor core WMMA GEMM

__global__ void b1_comet_xor_gemm_gpu_tc_simple(int m, int n, int k,
  GMBits2x64* a, GMBits2x64* b, bool beta, GMTally2x2* c) {
  using namespace nvcuda;

  // Block and thread indices
  int tx = threadIdx.x, ty = threadIdx.y;
  int bx = blockIdx.x, by = blockIdx.y;

  int COMET_BLOCK_SIZE = 8;
  int CROWS = 2, CCOLS = 2;
  int KSIZE = 4;

  //if(bx==0 && by==0 && tx==0 && ty==0) printf("In b1_comet_xor_gemm_gpu_simple mnk=%d,%d,%d a=%dx%d * b=%dx%d = c=%dx%d\n",m,n,k,m,k,k,n,m,n);

  enum { GM_TALLY1_MAX_VALUE_BITS = 26 };

  // Matrix block location
  int aBegin = k * COMET_BLOCK_SIZE * bx;
  int bBegin = k * COMET_BLOCK_SIZE * by;

  // Declare fragments
  wmma::fragment<wmma::matrix_a, WMMA1B_M, WMMA1B_N, WMMA1B_K, wmma::experimental::precision::b1,
                 wmma::row_major> a_frag;
  wmma::fragment<wmma::matrix_b, WMMA1B_M, WMMA1B_N, WMMA1B_K, wmma::experimental::precision::b1,
                 wmma::col_major> b_frag;
  wmma::fragment<wmma::accumulator, WMMA1B_M, WMMA1B_N, WMMA1B_K, int> c_frag[2][2];
  for(int ii=0; ii<CROWS; ii++) {
    for(int jj=0; jj<CCOLS; jj++) {
      wmma::fill_fragment(c_frag[ii][jj], 0);
    }
  }
  __syncthreads();

  __shared__ uint64_t As[16][4], Bs[16][4];
  __shared__ int      Cs[16][16];

  int nrow  = tx;
  int kval  = ty;

  // Each thread computes one element of block sub-matrix
  for (int l=0; l<k; l+=KSIZE) {

    // A row major and B col major
    int aInd = aBegin + l + k*nrow + kval;
    int bInd = bBegin + l + k*nrow + kval;

    /*printf("b=%d,%d t=%d,%d l=%d/%d nrow=%d kval=%d KSIZE=%d a=%d=%d b=%d=%d a2=%d=%d=%d b2=%d=%d=%d\n",
           bx,by,tx,ty,l,k,nrow,kval,KSIZE,
           aBegin,aBegin+k*nrow+kval,aBegin,aBegin+l,aBegin+l+k*nrow+kval,
           bBegin,bBegin+k*nrow+kval,bBegin,bBegin+l,bBegin+l+k*nrow+kval);*/

    //---Extract input values to process---
    const GMBits2x64 vi = a[aInd];
    const GMBits2x64 vj = b[bInd];

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
    const uint64_t vj0 = vj.data[0];
    const uint64_t vj1 = vj.data[1];

    // Compute masks to sample the single needed bit from each seminibble,
    // and to ignore undefined vector entries.
    const uint64_t oddbits = 0x5555555555555555;
    const uint64_t vi0mask = (vi0 | ~(vi0 >> 1)) & oddbits;
    const uint64_t vi1mask = (vi1 | ~(vi1 >> 1)) & oddbits;
    const uint64_t vj0mask = (vj0 | ~(vj0 >> 1)) & oddbits;
    const uint64_t vj1mask = (vj1 | ~(vj1 >> 1)) & oddbits;

    // Extract elts that are a "1" bit (=01).
    const uint64_t pvi0 =  vi0  & vi0mask;
    const uint64_t pvi1 =  vi1  & vi1mask;
    const uint64_t pvj0 =  vj0  & vj0mask;
    const uint64_t pvj1 =  vj1  & vj1mask;

    // Extract elts that are an "0" bit (=00).
    const uint64_t nvi0 = ~vi0  & vi0mask;
    const uint64_t nvi1 = ~vi1  & vi1mask;
    const uint64_t nvj0 = ~vj0  & vj0mask;
    const uint64_t nvj1 = ~vj1  & vj1mask;

    // Combine lower, upper words - each only uses odd bits - make packed.
    const uint64_t pvi = pvi0 | (pvi1 << 1);
    const uint64_t pvj = pvj0 | (pvj1 << 1);
    const uint64_t nvi = nvi0 | (nvi1 << 1);
    const uint64_t nvj = nvj0 | (nvj1 << 1);

    int nbrow  = nrow * 2;
    int nbcol  = kval;

    As[nbrow  ][nbcol] = nvi;
    As[nbrow+1][nbcol] = pvi;

    Bs[nbrow  ][nbcol] = nvj;
    Bs[nbrow+1][nbcol] = pvj;
    __syncthreads();

    // Compute each tensor core operation
    for(int ii=0; ii<CROWS; ii++) {
      for(int jj=0; jj<CCOLS; jj++) {
        for(int kk=0; kk<2; kk++) {

          // Load the inputs
          wmma::load_matrix_sync(a_frag, &(As[ii*8][kk*2]), 4*64);
          wmma::load_matrix_sync(b_frag, &(Bs[jj*8][kk*2]), 4*64);

          // Perform the matrix multiplication
          wmma::bmma_sync(c_frag[ii][jj], a_frag, b_frag, c_frag[ii][jj]);
        }
      }
    }
  }

  // Store the output
  for(int ii=0; ii<CROWS; ii++) {
    for(int jj=0; jj<CCOLS; jj++) {
      wmma::store_matrix_sync(&(Cs[ii*8][jj*8]), c_frag[ii][jj], 16, wmma::mem_row_major);
    }
  }
  __syncthreads();

  // Store results in correct format
  for(int kk=0; kk<2; kk++) {

    //---Accumulate---
    int crow = tx*2;
    int ccol = ty*4 + kk*2;
    double c0 = 0, c1 = 0; 
    const uint64_t r00 = Cs[crow][ccol];
    const uint64_t r01 = Cs[crow][ccol+1];
    const uint64_t r10 = Cs[crow+1][ccol];
    const uint64_t r11 = Cs[crow+1][ccol+1];

    c0 = r00 | (r01 << GM_TALLY1_MAX_VALUE_BITS);
    c1 = r10 | (r11 << GM_TALLY1_MAX_VALUE_BITS);

    // Each thread writes one element of block sub-matrix to memory
    int cBegin = n*bx*COMET_BLOCK_SIZE+by*COMET_BLOCK_SIZE;
    int cInd   = cBegin + tx * 8 * gridDim.x + ty * 2 + kk;
    c[cInd].data[0] = c0;
    c[cInd].data[1] = c1;

    /*printf("b=%d,%d t=%d,%d crow=%d ccol=%d Cs[%d,%d]=%ld %ld %ld %ld cInd=%d c0=%f c1=%f\n",
           bx,by,tx,ty,crow,ccol,crow/2,ccol/2,
           r00,r01,r10,r11,cInd,c0,c1);*/
  }
}

//-----------------------------------------------------------------------------
/// \brief GPU kernel for custom 1-bit tensor core WMMA GEMM

#define BLOCK_DIM_X 8
#define BLOCK_DIM_Y 4

#define CROWS BLOCK_DIM_X*2/8
#define CCOLS BLOCK_DIM_Y/2

__global__ void b1_comet_xor_gemm_gpu_tc_opt(int m, int n, int k,
  GMBits2x64* a, GMBits2x64* b, bool beta, GMTally2x2* c) {
  using namespace nvcuda;

  // Block and thread indices
  int tx = threadIdx.x, ty = threadIdx.y;
  int bx = blockIdx.x, by = blockIdx.y;

  //if(bx==0 && by==0 && tx==0 && ty==0)
  //  printf("In b1_comet_xor_gemm_gpu_opt mnk=%d,%d,%d a=%dx%d * b=%dx%d = c=%dx%d\n",
  //         m,n,k,m,k,k,n,m,n);

  enum { GM_TALLY1_MAX_VALUE_BITS = 26 };

  // Matrix block location
  int aBegin = k * BLOCK_DIM_X * bx;
  int bBegin = k * BLOCK_DIM_X * by;

  // Declare fragments
  wmma::fragment<wmma::matrix_a, WMMA1B_M, WMMA1B_N, WMMA1B_K, wmma::experimental::precision::b1,
                 wmma::row_major> a_frag;
  wmma::fragment<wmma::matrix_b, WMMA1B_M, WMMA1B_N, WMMA1B_K, wmma::experimental::precision::b1,
                 wmma::col_major> b_frag;
  wmma::fragment<wmma::accumulator, WMMA1B_M, WMMA1B_N, WMMA1B_K, int> c_frag[CROWS][CCOLS];
  for(int ii=0; ii<CROWS; ii++) {
    for(int jj=0; jj<CCOLS; jj++) {
      wmma::fill_fragment(c_frag[ii][jj], 0);
    }
  }
  __syncthreads();

  __shared__ uint64_t As[BLOCK_DIM_X*2][BLOCK_DIM_Y], Bs[BLOCK_DIM_X*2][BLOCK_DIM_Y];
  __shared__ int      Cs[BLOCK_DIM_X*2][BLOCK_DIM_X*2];

  int nrow  = tx;
  int kval  = ty;

  // Each thread computes one element of block sub-matrix
  for (int l=0; l<k; l+=BLOCK_DIM_Y) {

    // A row major and B col major
    int aInd = aBegin + l + k*nrow + kval;
    int bInd = bBegin + l + k*nrow + kval;
    //printf("b=%d,%d t=%d,%d Starting loop l=%d/%d ngrid=%d nrow=%d kval=%d a=%d=%d b=%d=%d\n",
    //       bx,by,tx,ty,l,k/4,ngrid,nrow,kval,aBegin,aInd,bBegin,bInd);

    //---Extract input values to process---
    const GMBits2x64 vi = a[aInd];
    const GMBits2x64 vj = b[bInd];

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
    const uint64_t vj0 = vj.data[0];
    const uint64_t vj1 = vj.data[1];

    // Compute masks to sample the single needed bit from each seminibble,
    // and to ignore undefined vector entries.
    const uint64_t oddbits = 0x5555555555555555;
    const uint64_t vi0mask = (vi0 | ~(vi0 >> 1)) & oddbits;
    const uint64_t vi1mask = (vi1 | ~(vi1 >> 1)) & oddbits;
    const uint64_t vj0mask = (vj0 | ~(vj0 >> 1)) & oddbits;
    const uint64_t vj1mask = (vj1 | ~(vj1 >> 1)) & oddbits;

    // Extract elts that are a "1" bit (=01).
    const uint64_t pvi0 =  vi0  & vi0mask;
    const uint64_t pvi1 =  vi1  & vi1mask;
    const uint64_t pvj0 =  vj0  & vj0mask;
    const uint64_t pvj1 =  vj1  & vj1mask;

    // Extract elts that are an "0" bit (=00).
    const uint64_t nvi0 = ~vi0  & vi0mask;
    const uint64_t nvi1 = ~vi1  & vi1mask;
    const uint64_t nvj0 = ~vj0  & vj0mask;
    const uint64_t nvj1 = ~vj1  & vj1mask;

    // Combine lower, upper words - each only uses odd bits - make packed.
    const uint64_t pvi = pvi0 | (pvi1 << 1);
    const uint64_t pvj = pvj0 | (pvj1 << 1);
    const uint64_t nvi = nvi0 | (nvi1 << 1);
    const uint64_t nvj = nvj0 | (nvj1 << 1);

    int nbrow  = nrow * 2;
    int nbcol  = kval;

    //if(tx==0 && ty==0) {
    As[nbrow  ][nbcol] = nvi;
    As[nbrow+1][nbcol] = pvi;

    Bs[nbrow  ][nbcol] = nvj;
    Bs[nbrow+1][nbcol] = pvj;
    /*} else {
    As[nbrow  ][nbcol] = 0;
    As[nbrow+1][nbcol] = 0;

    Bs[nbrow  ][nbcol] = 0;
    Bs[nbrow+1][nbcol] = 0;
    }*/
    __syncthreads();

    //printf("b=%d,%d t=%d,%d ngrid=%d nrow=%d nbrow=%d nbcol=%d A=%ld %ld B=%ld %ld\n",
    //       bx,by,tx,ty,ngrid,nrow,nbrow,nbcol,
    //       As[nbrow][nbcol],As[nbrow+1][nbcol],Bs[nbrow][nbcol],Bs[nbrow+1][nbcol]);

    // Compute each tensor core operation
    for(int ii=0; ii<CROWS; ii++) {
      for(int jj=0; jj<CCOLS; jj++) {
        for(int kk=0; kk<2; kk++) {
          //if(tx==0 && ty==0) printf("Calling bmma with ntc=%d As[%d]=%ld %ld %ld %ld Bs[%d]=%ld %ld %ld %ld\n",
          //                          ntc,ntc*2,As[0][ntc*2],As[1][ntc*2+1],As[2][ntc*2+8],As[3][ntc*2+9],
          //                          ntc*2,Bs[0][ntc*2],Bs[1][ntc*2+1],Bs[2][ntc*2+8],Bs[3][ntc*2+9]);

          // Load the inputs
          wmma::load_matrix_sync(a_frag, &(As[ii*8][kk*2]), BLOCK_DIM_Y*64);
          wmma::load_matrix_sync(b_frag, &(Bs[jj*8][kk*2]), BLOCK_DIM_Y*64);

          // Perform the matrix multiplication
          wmma::bmma_sync(c_frag[ii][jj], a_frag, b_frag, c_frag[ii][jj]);
          //__syncthreads();
        }
      }
    }
  }

  // Store the output
  for(int ii=0; ii<CROWS; ii++) {
    for(int jj=0; jj<CCOLS; jj++) {
      wmma::store_matrix_sync(&(Cs[ii*8][jj*8]), c_frag[ii][jj], BLOCK_DIM_X*2, wmma::mem_row_major);
    }
  }
  __syncthreads();

  /*for(int ii=0; ii<c_frag[0][0].num_elements; ii++) {
    printf("b=%d,%d t=%d,%d c_frag[%d/%d]=%d\n",
           bx,by,tx,ty,ii,c_frag[0][0].num_elements,c_frag[0][0].x[ii]);
  }
  __syncthreads();

  for(int kk=0; kk<2; kk++) {
    int crow = tx;
    int ccol = ty*2+kk;
    printf("b=%d,%d t=%d,%d Cs[%d][%d]=%d\n",
           bx,by,tx,ty,crow,ccol,Cs[crow][ccol]);
  }
  __syncthreads();*/

  /*if(kval==0) {
    printf("b=%d,%d t=%d,%d nrow=%d Cs]=%d %d %d %d\n",
           bx,by,tx,ty,nrow,
           Cs[nrow*2][ncol*2],Cs[nrow*2][ncol*2+1],Cs[nrow*2+1][ncol*2],Cs[nrow*2+1][ncol*2+1]);
  }*/

  // Store results in correct format
  for(int kk=0; kk<2; kk++) {

    // Set result to double complex format
    int crow = tx*2;
    int ccol = ty*4 + kk*2;
    double c0 = 0, c1 = 0; 
    const uint64_t r00 = Cs[crow][ccol];
    const uint64_t r01 = Cs[crow][ccol+1];
    const uint64_t r10 = Cs[crow+1][ccol];
    const uint64_t r11 = Cs[crow+1][ccol+1];

    //
    /*int cInd = 8 * 2 * ty + tx * 2;
    double c0 = 0, c1 = 0;
    const uint64_t r00 = c_frag[0][0].x[cInd];
    const uint64_t r01 = Cs[crow][ccol+1];
    const uint64_t r10 = Cs[crow+1][ccol];
    const uint64_t r11 = Cs[crow+1][ccol+1];*/

    c0 = r00 | (r01 << GM_TALLY1_MAX_VALUE_BITS);
    c1 = r10 | (r11 << GM_TALLY1_MAX_VALUE_BITS);

    // Each thread writes one element of block sub-matrix to memory
    int cBegin = n*bx*BLOCK_DIM_X+by*BLOCK_DIM_X;
    int cInd   = cBegin + tx * 8 * gridDim.x + ty * 2 + kk;
    c[cInd].data[0] = c0;
    c[cInd].data[1] = c1;

    //printf("b=%d,%d t=%d,%d crow=%d ccol=%d Cs[%d,%d]=%d %d %d %d cInd=%d c0=%f c1=%f\n",
    //       bx,by,tx,ty,crow,ccol,crow/2,ccol/2,
    //       Cs[crow][ccol],Cs[crow][ccol+1],Cs[crow+1][ccol],Cs[crow+1][ccol+1],cInd,c0,c1);
  }
}

//-----------------------------------------------------------------------------
/// \brief Perform required GEMM.

template<int TC_METHOD>
static void tc_solve_comet_impl(bool is_first, int m, int n, int k,
  const void *matA, const void *matB, void* matC, TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST(matC);
  COMET_INSIST(m >= 0 && n >= 0 && k >= 0);

  if(env.print_details()) printf("In tc_solve_comet_impl mnk=%d,%d,%d\n",m,n,k);
  double tbegin = env.synced_time();

  const bool beta = 1;
  int threadblockx, threadblocky, gridblockx, gridblocky;

  if(env.print_details())
    printf("Launching 1-bit GEMM kernel mnk=%d,%d,%d beta=%d\n",m,n,k,(int)beta);

  switch(env.num_kernel()) {
    // Basic GEMM
    case 21: {
      threadblockx = BLOCK_SIZE; threadblocky = BLOCK_SIZE;
      gridblockx = (int)ceil((double)m/threadblockx);
      gridblocky = (int)ceil((double)n/threadblocky);
      if(env.print_details()) printf("Calling b1_comet_xor_gemm_gpu_simple kernel gridDim=%d,%d threadDim=%d,%d\n",gridblockx,gridblocky,threadblockx,threadblocky);
      COMET_LAUNCH_KERNEL(b1_comet_xor_gemm_gpu_simple,
        dim3(gridblockx, gridblocky, 1),
        dim3(threadblockx, threadblocky, 1), 0, env.stream_compute(),
        m, n, k, (GMBits2x64*)matA,
        (GMBits2x64*)matB, beta, (GMTally2x2*)matC);
    } break;

    // Simple tensor core GEMM
    case 22: {
      threadblockx = 8; threadblocky = 4;
      gridblockx = (int)ceil((double)m/threadblockx);
      gridblocky = (int)ceil((double)n/threadblockx);
      if(env.print_details()) printf("Calling b1_comet_xor_gemm_gpu_tc_simple kernel gridDim=%d,%d threadDim=%d,%d\n",gridblockx,gridblocky,threadblockx,threadblocky);
      COMET_LAUNCH_KERNEL(b1_comet_xor_gemm_gpu_tc_simple,
        dim3(gridblockx, gridblocky, 1),
        dim3(threadblockx, threadblocky, 1), 0, env.stream_compute(),
        m, n, k, (GMBits2x64*)matA,
        (GMBits2x64*)matB, beta, (GMTally2x2*)matC);
    } break;

    // Optimized tensor core GEMM
    case 23: {
      threadblockx = 8; threadblocky = 4;
      gridblockx = (int)ceil((double)m/threadblockx);
      gridblocky = (int)ceil((double)n/threadblockx);
      if(env.print_details()) printf("Calling b1_comet_xor_gemm_gpu_tc_opt kernel gridDim=%d,%d threadDim=%d,%d\n",gridblockx,gridblocky,threadblockx,threadblocky);
      COMET_LAUNCH_KERNEL(b1_comet_xor_gemm_gpu_tc_opt,
        dim3(gridblockx, gridblocky, 1),
        dim3(threadblockx, threadblocky, 1), 0, env.stream_compute(),
        m, n, k, (GMBits2x64*)matA,
        (GMBits2x64*)matB, beta, (GMTally2x2*)matC);
    } break;

    // Nvidia optimized Cutlass GEMM
    case 24: {
      if(env.print_details()) printf("Calling tc_solve_comet_impl_cutlass\n");
      tc_solve_comet_impl_cutlass(m,n,k,(GMBits2x64*)matA,
        (GMBits2x64*)matB, (GMTally2x2*)matC);
    } break;

    /*case 30: {
      
    } break;*/

    // Output error for invalid choice
    default: {
      printf("Failed to call appropriate 1-bit GEMM kernel for num_kernel=%d\n",
         env.num_kernel());
      COMET_INSIST(false && "Failure to call GEMM function.");
    }
  }
  System::accel_last_call_succeeded();
  env.ops_local_inc(2 * m * (double)n * (double)k);

  env.gemmtime_inc(env.synced_time() - tbegin);
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

  const int m = nvll; // metrics array dim
  const int n = nvl; // metrics array dim
  const int k = nfl_thisstep; // vectors array (as GemmIn_t) dim

  if(env.print_details()) printf("Calling tc_solve_comet_impl with mnk=%d,%d,%d\n",m,n,k);
  tc_solve_comet_impl<TC_METHOD>(is_first, m, n, k, matA, matB, matC, tc_bufs, env);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TC_SOLVE_XOR_I_HH_

//-----------------------------------------------------------------------------
