//-----------------------------------------------------------------------------
/*!
 * \file   tc_solve_comet_xor_int.i.hh
 * \author Paul Eller
 * \date   Tue April 27 17:44:00 EST 2021
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

#ifndef _COMET_TC_SOLVE_XOR_INT_I_HH_
#define _COMET_TC_SOLVE_XOR_INT_I_HH_

// Includes
#include "cstdlib"
#include <cuda_runtime.h>
#include "cuda_fp16.h"
#include "types.hh"

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
/// Assume A col major, B and C row major
__global__
void b1_comet_xor_gemm_gpu_int_simple(int m, int n, int k,
  GMBits2x64* a, GMBits2x64* b, bool beta, int32_t* c) {

  // Block and thread indices
  int tx = threadIdx.x, ty = threadIdx.y;
  int bx = blockIdx.x, by = blockIdx.y;

  int gridx = bx*BLOCK_SIZE + tx;
  int gridy = by*BLOCK_SIZE + ty;

  //if(bx==0 && by==0 && tx==0 && ty==0) printf("In b1_comet_xor_gemm_gpu_int_simple mnk=%d,%d,%d a=%dx%d * b=%dx%d = c=%dx%d gxy=%d,%d\n",m,n,k,m,k,k,n,m,n,gridx,gridy);

  if(gridy>=m || gridx>=n) return;

  // Matrix block location
  int aBegin = k * BLOCK_SIZE * by;
  int bBegin = k * BLOCK_SIZE * bx;

  // Stores element of block sub-matrix computed by thread
  int32_t c0 = 0, c1 = 0, c2 = 0, c3 = 0;

  // Each thread computes one element of block sub-matrix
  for (int l=0; l<k; ++l) {

    // A col major and B row major
    int aInd = aBegin + k*ty + l;
    int bInd = bBegin + k*tx + l;

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
    //  printf("b=%d,%d t=%d,%d g=%d,%d a=%d=%d b=%d=%d mnk=%d,%d,%d vi0=%lu vi1=%lu vj0=%lu vj1=%lu\n",
    //         bx,by,tx,ty,gridDim.x,gridDim.y,aBegin,aInd,bBegin,bInd,m,n,k,vi0,vi1,vj0,vj1);

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

    //if(tx==0 && ty==0)
    //  printf("b=%d,%d t=%d,%d a=%d b=%d r00=%ld r01=%ld r10=%ld r11=%ld c0123=%d,%d,%d,%d\n",
    //         bx,by,tx,ty,aInd,bInd,r00,r01,r10,r11,c0,c1,c2,c3);
  }

  // Different ordering from original 1-bit routines
  // Each thread writes one element of block sub-matrix to memory
  // Assume c is row major
  /*int cBegin = n*bx*BLOCK_SIZE*2+by*BLOCK_SIZE*2;
  int rind1 = ty * n*2;
  int rind2 = ty * n*2 + n;
  int cind = (tx % 4)*4 + (tx/4)*2;
  int cInd1 = cBegin + rind1 + cind;
  int cInd2 = cBegin + rind2 + cind;*/
  //int cBegin = n*bx*BLOCK_SIZE*2+by*BLOCK_SIZE*2;
  int cBegin = by*n*2*BLOCK_SIZE*2;
  int rind1 = ty * n*4;
  int rind2 = ty * n*4 + n*2;
  //int cind = (tx % (n/2))*4 + (tx/(n/2))*2;
  int cind = ((tx + bx*BLOCK_SIZE) % (n/2))*4 + ((tx + bx*BLOCK_SIZE) / (n/2))*2;
  int cInd1 = cBegin + rind1 + cind;
  int cInd2 = cBegin + rind2 + cind;

  //if(tx==0 && ty==0)
  //  printf("b=%d,%d t=%d,%d cb=%d ci1=%d=%d,%d ci2=%d=%d,%d c0123=%d,%d,%d,%d\n",
  //         bx,by,tx,ty,cBegin,cInd1,rind1,cind,cInd2,rind2,cind,c0,c1,c2,c3);

  c[cInd1] = c0; c[cInd1+1] = c1;
  c[cInd2] = c2; c[cInd2+1] = c3;
}

//-----------------------------------------------------------------------------
/// \brief Perform required GEMM.

template<int TC_METHOD>
static void tc_solve_comet_int_impl(bool is_first, int m, int n, int k,
  const void *matA, const void *matB, void* matC, TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST(matC);
  COMET_INSIST(m >= 0 && n >= 0 && k >= 0);

  if(env.print_details()) printf("In tc_solve_comet_int_impl mnk=%d,%d,%d\n",m,n,k);
  //double tbegin = env.get_cpu_time();

  const bool beta = is_first ? 0 : 1;
  const int threadblockx = BLOCK_SIZE, threadblocky = BLOCK_SIZE;
  int gridblockx = (int)ceil((double)m/threadblockx);
  int gridblocky = (int)ceil((double)n/threadblocky);

  if(env.print_details())
    printf("Launching 1-bit GEMM kernel mnk=%d,%d,%d beta=%d "
           "gridDim=%d,%d threadDim=%d,%d\n",
           m,n,k,(int)beta,gridblockx,gridblocky,threadblockx,threadblocky);

  switch(env.num_kernel()) {
    // Basic GEMM
    case 125: {
      if(env.print_details()) printf("Calling b1_comet_xor_gemm_gpu_int_simple kernel\n");
      COMET_LAUNCH_KERNEL(b1_comet_xor_gemm_gpu_int_simple,
        dim3(gridblockx, gridblocky, 1),
        dim3(threadblockx, threadblocky, 1), 0, env.stream_compute(),
        m, n, k, (GMBits2x64*)matA,
        (GMBits2x64*)matB, beta, (int32_t*)matC);
    } break;

    // Optimized Cutlass GEMM that outputs ints
    case 126: {
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

  //env.gemmtime_inc(env.get_cpu_time() - tbegin);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TC_SOLVE_XOR_INT_I_HH_

//-----------------------------------------------------------------------------
