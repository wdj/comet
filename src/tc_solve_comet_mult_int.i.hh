//-----------------------------------------------------------------------------
/*!
 * \file   tc_solve_comet_mult_int.i.hh
 * \author Paul Eller
 * \date   Tue April 27 17:43:00 EST 2021
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

#ifndef _COMET_TC_SOLVE_MULT_INT_I_HH_
#define _COMET_TC_SOLVE_MULT_INT_I_HH_

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
/// \brief GPU kernel for simple 1-bit tensor core CoMet int GEMM
/// Assume A col major, B and C row major
__global__
void b1_comet_mult_gemm_gpu_int_simple(int m, int n, int k,
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
  uint32_t c0=0, c1=0, c2=0, c3=0;

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
      printf("b=%d,%d t=%d,%d g=%d,%d a=%d=%d b=%d=%d mnk=%d,%d,%d vi0=%lu vi1=%lu vj0=%lu vj1=%lu\n",
             bx,by,tx,ty,gridDim.x,gridDim.y,aBegin,aInd,bBegin,bInd,m,n,k,vi0,vi1,vj0,vj1);

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

    const uint64_t r00 = gm_popcount64(nvi & nvj);
    const uint64_t r01 = gm_popcount64(nvi & pvj);
    const uint64_t r10 = gm_popcount64(pvi & nvj);
    const uint64_t r11 = gm_popcount64(pvi & pvj);

    //---Accumulate---
    c0 += r00; c1 += r01;
    c2 += r10; c3 += r11;

    //if(tx==0 && ty==0)
      printf("b=%d,%d t=%d,%d a=%d b=%d r00=%ld r01=%ld r10=%ld r11=%ld c0123=%d,%d,%d,%d\n",
             bx,by,tx,ty,aInd,bInd,r00,r01,r10,r11,c0,c1,c2,c3);
  }

  // Different ordering from original 1-bit routines
  // Each thread writes one element of block sub-matrix to memory
  // Assume c is row major
  int cBegin = n*2*bx*BLOCK_SIZE*2+by*BLOCK_SIZE*2;
  int cInd1  = cBegin + tx*2*n*2 + ty*2;
  int cInd2  = cBegin + (tx*2+1)*n*2 + ty*2;
  // Partial change
  //int cBegin = n*2*bx*BLOCK_SIZE*2+by*BLOCK_SIZE*2;
  //int cInd1  = cBegin + ty*2*m*2 + tx*2;
  //int cInd2  = cBegin + (ty*2+1)*m*2 + tx*2;

  // Assume c is col major
  //int cBegin = m*2*by*BLOCK_SIZE*2+bx*BLOCK_SIZE*2;
  //int cInd1  = cBegin + ty*2*m*2 + tx*2;
  //int cInd2  = cBegin + (ty*2+1)*m*2 + tx*2;

  // Assume all four entries are sequential
  //int cBegin = n*4*bx*BLOCK_SIZE + by*BLOCK_SIZE*4;
  //int cInd   = cBegin + tx*n*4 + ty*4;
  //if(tx==0 && ty==0)
  //  printf("b=%d,%d t=%d,%d cb=%d=%d,%d ci=%d=%d,%d c0123=%d,%d,%d,%d\n",
  //         bx,by,tx,ty,cBegin,n*4*bx*BLOCK_SIZE,by*BLOCK_SIZE*4,cInd,tx*n*4,ty*4,c0,c1,c2,c3);
  //c[cInd]=c0; c[cInd+1]=c1; c[cInd+2]=c2; c[cInd+3]=c3;

  //if(tx==0 && ty==0)
  //  printf("b=%d,%d t=%d,%d cb=%d=%d,%d ci1=%d=%d ci2=%d=%d c0123=%d,%d,%d,%d\n",
  //         bx,by,tx,ty,cBegin,n*2*bx*BLOCK_SIZE*2,by*BLOCK_SIZE*2,cInd1,tx*2*n*2+ty*2,cInd2,(tx*2+1)*n*2+ty*2,c0,c1,c2,c3);
  c[cInd1] = c0; c[cInd1+1] = c1;
  c[cInd2] = c2; c[cInd2+1] = c3;
}


//-----------------------------------------------------------------------------
/// \brief Perform required GEMM.

template<int TC_METHOD>
static void tc_solve_comet_mult_int_impl(bool is_first, int m, int n, int k,
  const void *matA, const void *matB, void* matC, TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST(matC);
  COMET_INSIST(m >= 0 && n >= 0 && k >= 0);

  if(env.print_details()) printf("In tc_solve_comet_mult_int_impl mnk=%d,%d,%d\n",m,n,k);
  //double tbegin = env.get_cpu_time();

  const bool beta = 1;
  const int threadblockx = BLOCK_SIZE, threadblocky = BLOCK_SIZE;
  int gridblockx = (int)ceil((double)m/threadblockx);
  int gridblocky = (int)ceil((double)n/threadblocky);

  if(env.print_details())
    printf("Launching 1-bit Mult Int GEMM kernel mnk=%d,%d,%d beta=%d "
           "gridDim=%d,%d threadDim=%d,%d\n",
           m,n,k,(int)beta,gridblockx,gridblocky,threadblockx,threadblocky);

  switch(env.num_kernel()) {
    // Basic GEMM
    case 175: {
      if(env.print_details()) printf("Calling b1_comet_mult_gemm_gpu_int_simple kernel\n");
      COMET_LAUNCH_KERNEL(b1_comet_mult_gemm_gpu_int_simple,
        dim3(gridblockx, gridblocky, 1),
        dim3(threadblockx, threadblocky, 1), 0, env.stream_compute(),
        m, n, k, (GMBits2x64*)matA,
        (GMBits2x64*)matB, beta, (int32_t*)matC);
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

#endif // _COMET_TC_SOLVE_MULT_INT_I_HH_

//-----------------------------------------------------------------------------
