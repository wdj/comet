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

#ifndef _COMET_TC_SOLVE_MULT_I_HH_
#define _COMET_TC_SOLVE_MULT_I_HH_

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
/// \brief GPU kernel for simple 1-bit CoMet GEMM kernel

__global__ void b1_comet_mult_gemm_gpu_simple(int m, int n, int k,
  GMBits2x64* a, GMBits2x64* b, bool beta, GMTally2x2* c) {

  // Block and thread indices
  int tx = threadIdx.x, ty = threadIdx.y;
  int bx = blockIdx.x, by = blockIdx.y;

  //if(tx==0 && ty==0 && bx==0 && by==0) printf("In b1_comet_mult_gemm_gpu_simple\n");

  int gridx = bx*BLOCK_SIZE + tx;
  int gridy = by*BLOCK_SIZE + ty;

  //if(bx==0 && by==0 && tx==0 && ty==0)
  //  printf("In b1_comet_mult_gemm_gpu_simple mnk=%d,%d,%d a=%dx%d * b=%dx%d = c=%dx%d gxy=%d,%d\n",
  //         m,n,k,m,k,k,n,m,n,gridx,gridy);

  if(gridy>=m || gridx>=n) return;

  enum { GM_TALLY1_MAX_VALUE_BITS = 26 };

  // Matrix block location
  int aBegin = k * BLOCK_SIZE * by;
  int bBegin = k * BLOCK_SIZE * bx;

  // Stores element of block sub-matrix computed by thread
  double c0 = 0, c1 = 0;
  int32_t ci0 = 0, ci1 = 0, ci2 = 0, ci3 = 0;

  // Each thread computes one element of block sub-matrix
  for (int l=0; l<k; ++l) {

    // A row major and B col major
    int aInd = aBegin + k*ty + l;
    int bInd = bBegin + k*tx + l;

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

    const uint64_t r00 = gm_popcount64(nvi & nvj);
    const uint64_t r01 = gm_popcount64(nvi & pvj);
    const uint64_t r10 = gm_popcount64(pvi & nvj);
    const uint64_t r11 = gm_popcount64(pvi & pvj);

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
  //if(beta) {
  //  c[cInd].data[0] += c0;
  //  c[cInd].data[1] += c1;
  //} else {
  c[cInd].data[0] = c0;
  c[cInd].data[1] = c1;
  //}
  //printf("b=%d,%d t=%d,%d cb=%d=%d,%d ci=%d c01=%lf,%lf ci0123=%d,%d,%d,%d\n",
  //       bx,by,tx,ty,cBegin,n*bx*BLOCK_SIZE,by*BLOCK_SIZE,cInd,c0,c1,ci0,ci1,ci2,ci3);
}

//-----------------------------------------------------------------------------
/// \brief Perform required GEMM.

template<int TC_METHOD>
static void tc_solve_comet_mult_impl(bool is_first, int m, int n, int k,
  const void *matA, const void *matB, void* matC, TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST(matC);
  COMET_INSIST(m >= 0 && n >= 0 && k >= 0);

  if(env.print_details()) printf("\nIn tc_solve_comet_mult_impl mnk=%d,%d,%d num_kernel=%d\n",
    m,n,k,env.num_kernel());
  //double tbegin = env.get_cpu_time();

  const bool beta = is_first ? 0 : 1;
  int threadblockx, threadblocky, gridblockx, gridblocky;

  if(env.print_details())
    printf("Launching 1-bit GEMM kernel mnk=%d,%d,%d beta=%d num_kernel=%d\n",
           m,n,k,(int)beta,env.num_kernel());

  env.gemm_timer.record();

  env.gemm_timer.start();

  switch(env.num_kernel()) {
    // Basic GEMM
    case 150: {
      threadblockx = BLOCK_SIZE; threadblocky = BLOCK_SIZE;
      gridblockx = (int)ceil((double)m/threadblockx);
      gridblocky = (int)ceil((double)n/threadblocky);
      if(env.print_details()) printf("Calling b1_comet_mult_gemm_gpu_simple kernel gridDim=%d,%d threadDim=%d,%d\n",gridblockx,gridblocky,threadblockx,threadblocky);
      COMET_LAUNCH_KERNEL(b1_comet_mult_gemm_gpu_simple,
        dim3(gridblockx, gridblocky, 1),
        dim3(threadblockx, threadblocky, 1), 0, env.stream_compute(),
        m, n, k, (GMBits2x64*)matA,
        (GMBits2x64*)matB, beta, (GMTally2x2*)matC);
    } break;
    // Output error for invalid choice
    default: {
      printf("Failed to call appropriate 1-bit GEMM kernel for num_kernel=%d\n",
         env.num_kernel());
      COMET_INSIST(false && "Failure to call GEMM function.");
    }
  }

  env.gemm_timer.end();

  int err = cudaGetLastError();
  if(env.print_details()) printf("tc_solve_comet_mult_impl computed 1-bit GEMM with 2x%dx%dx%d=%lf operations\n",
    2*m,2*n,k*64,2.0*(double)m*2.0*(double)n*2.0*(double)k*64.0);
  //env.stream_synchronize(env.stream_compute());
  System::accel_last_call_succeeded();
  env.ops_gemm_local_inc(2.0 * (double)m*2.0 * (double)n*2.0 * (double)k*64.0);
  env.ops_local_inc(2.0 * (double)m*2.0 * (double)n*2.0 * (double)k*64.0);

  //env.gemmtime_inc(env.get_cpu_time() - tbegin);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TC_SOLVE_MULT_I_HH_

//-----------------------------------------------------------------------------
