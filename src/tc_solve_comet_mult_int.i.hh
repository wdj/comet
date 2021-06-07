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

__global__
void b1_print_matrix_int(int m, int n, int32_t* c) {

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

  int px=0, py=0;

  if(bx==px && by==py)
    printf("b=%d,%d t=%d,%d cb=%d ci1=%d=%d,%d ci2=%d=%d,%d c0123=%d,%d,%d,%d\n",
           bx,by,tx,ty,cBegin,cInd1,rind1,cind,cInd2,rind2,cind,c0,c1,c2,c3);
}

//-----------------------------------------------------------------------------
/// \brief GPU kernel for simple 1-bit tensor core CoMet int GEMM
/// Assume A col major, B and C row major
__global__
void b1_comet_mult_gemm_gpu_int_simple(int m, int n, int k, int pfl_min, int step_2way, int nfal,
  GMBits2x64* a, GMBits2x64* a2, GMBits2x64* b, bool beta, int32_t* c) {

  // Block and thread indices
  int tx = threadIdx.x, ty = threadIdx.y;
  int bx = blockIdx.x, by = blockIdx.y;

  const int nvleD2 = n/2;
  int gridx = step_2way*nvleD2 + bx*BLOCK_SIZE + tx;
  int gridy = by*BLOCK_SIZE + ty;

  int gridxt = gridx % nvleD2 + step_2way*nvleD2;
  int gridyt = gridy;

  int px=0, py=0;
  //int px=95, py=95;
  int ptx=2, pty=2;
  int pf=16;

  //if(bx==0 && by==0 && tx==0 && ty==0) printf("In b1_comet_xor_gemm_gpu_int_simple mnk=%d,%d,%d a=%dx%d * b=%dx%d = c=%dx%d gxy=%d,%d\n",m,n,k,m,k,k,n,m,n,gridx,gridy);

  if(gridyt>=m || gridxt>=n) return;

  // Matrix block location
  int aBegin = step_2way*nvleD2*k + k * BLOCK_SIZE * by;
  int bBegin = k * BLOCK_SIZE * bx;
  int aThread = BLOCK_SIZE * by + ty;
  //int bThread = BLOCK_SIZE * bx + tx; 

  // Stores element of block sub-matrix computed by thread
  int32_t c0=0, c1=0, c2=0, c3=0;

  // Variables for 3-way
  enum {NUM_FIELD_PER_PACKEDVAL_FIELD = 64};
  enum {NUM_FL_PER_PFL = NUM_FIELD_PER_PACKEDVAL_FIELD};

  enum {NGIPT = TCTraits<TC::B1>::NGIPT};
  enum {NFPGI = TCTraits<TC::B1>::NFPGI};

  enum {NUM_FIELD_PER_THREAD = NFPGI * NGIPT};
  enum {NFPT = NUM_FIELD_PER_THREAD};

  enum {BITS_PER_BYTE = 8};
  enum {BPSN = 2}; // bits per seminibble (field)
  enum {SNPW = sizeof(TCWord_t) * BITS_PER_BYTE / BPSN}; // seminibbles (fields) / word

  //const int fl_min = pfl_min * NUM_FL_PER_PFL;
  //const int flT_min = fl_min / NFPT;

  // Each thread computes one element of block sub-matrix
  int nsteps = k;
  for (int l=0; l<nsteps; ++l) {

    // A col major and B row major
    int aInd = aBegin + k*ty + l;
    int bInd = bBegin + k*tx + l;
    int vInd = l;

    //---Extract input values to process---
    const GMBits2x64 vi  = a2[aInd];
    const GMBits2x64 vj  = b[bInd];
    const GMBits2x64 vv  = a[vInd];

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
    const uint64_t vv0 = vv.data[0];
    const uint64_t vv1 = vv.data[1];

    if(bx==px && by==py && tx<ptx && ty<pty && l<pf)
      printf("b=%d,%d t=%d,%d g=%d,%d a=%d=%d b=%d=%d vind=%d mnk=%d,%d,%d l=%d vi01=%lu,%lu vj01=%lu,%lu vv01=%lu,%lu ~vv01=%lu,%lu\n",
             bx,by,tx,ty,gridDim.x,gridDim.y,aBegin,aInd,bBegin,bInd,vInd,m,n,k,l,vi0,vi1,vj0,vj1,vv0,vv1,~vv0,~vv1);

    // Compute masks to sample the single needed bit from each seminibble,
    // and to ignore undefined vector entries.

    const uint64_t oddbits = 0x5555555555555555;
    const uint64_t vi0mask = (vi0 | ~(vi0 >> 1));// & oddbits;
    const uint64_t vi1mask = (vi1 | ~(vi1 >> 1));// & oddbits;
    const uint64_t vj0mask = (vj0 | ~(vj0 >> 1));// & oddbits;
    const uint64_t vj1mask = (vj1 | ~(vj1 >> 1));// & oddbits;
    const uint64_t vv0mask = (vv0 | ~(vv0 >> 1));// & oddbits;
    const uint64_t vv1mask = (vv1 | ~(vv1 >> 1));// & oddbits;

    // Upper bound on fields computed by this thread.
    const int fl_max = (l+1) * SNPW;

    const int fl_inactive = utils::min(32, utils::max(0, fl_max - nfal));

    const TCWord_t allbits = ~(TCWord_t)0;

    const uint64_t field_active_mask = fl_inactive >= 32 ?
      ((uint64_t)0) : allbits >> (2 * fl_inactive);

    const int kE0 = aThread / nvleD2;
    const int kE1 = step_2way;

    // Extract elts that are a "1" bit (=01).
    const uint64_t pvi0 =  (kE0 ? vi0 : ~vi0) & vv0 & oddbits & field_active_mask & vi0mask & vv0mask;
    const uint64_t pvi1 =  (kE1 ? vi1 : ~vi1) & vv1 & oddbits & field_active_mask & vi1mask & vv1mask;
    //const uint64_t pvi0 =  vv0 & vi0 & vi0mask & vv0mask & oddbits & field_active_mask;
    //const uint64_t pvi1 =  vv1 & vi1 & vi1mask & vv1mask & oddbits & field_active_mask;
    //const uint64_t pvi0 =  vi0 & vi0mask;
    //const uint64_t pvi1 =  vi1 & vi1mask;
    const uint64_t pvj0 =  vj0 & oddbits & field_active_mask & vj0mask;
    const uint64_t pvj1 =  vj1 & oddbits & field_active_mask & vj1mask;

    // Extract elts that are an "0" bit (=00).
    const uint64_t nvi0 = (kE0 ? vi0 : ~vi0) & ~vv0 & oddbits & field_active_mask & vi0mask & vv0mask;
    const uint64_t nvi1 = (kE1 ? vi1 : ~vi1) & ~vv1 & oddbits & field_active_mask & vi1mask & vv1mask;
    //const uint64_t nvi0 = ~vv0 & vi0 & vi0mask & vv0mask & field_active_mask;
    //const uint64_t nvi1 = ~vv1 & vi1 & vi1mask & vv1mask & field_active_mask;
    //const uint64_t nvi0 = ~vi0 & vi0mask;
    //const uint64_t nvi1 = ~vi1 & vi1mask;
    const uint64_t nvj0 = ~vj0  & oddbits & field_active_mask & vj0mask;
    const uint64_t nvj1 = ~vj1  & oddbits & field_active_mask & vj1mask;

    // Combine lower, upper words - each only uses odd bits - make packed.
    const uint64_t pvi = pvi0 | (pvi1 << 1);
    const uint64_t pvj = pvj0 | (pvj1 << 1);
    const uint64_t nvi = nvi0 | (nvi1 << 1);
    const uint64_t nvj = nvj0 | (nvj1 << 1);

    if(bx==px && by==py && tx<ptx && ty<pty && l<pf)
      printf("b=%d,%d t=%d,%d fmask=%lu v0m=%lu v1m=%lu vj=%lu,%lu ~vj=%lu,%lu vimask=%lu,%lu vjmask=%lu,%lu vvmask=%lu,%lu\n",
        bx,by,tx,ty,field_active_mask,(kE0 ? vi0 : ~vi0),(kE1 ? vi1 : ~vi1),vj0,vj1,~vj0,~vj1,vi0mask,vi1mask,vj0mask,vj1mask,vv0mask,vv1mask);

    if(bx==px && by==py && tx<ptx && ty<pty && l<pf)    
      printf("b=%d,%d t=%d,%d abThread=%d,%d nvleD2=%d kE01=%d,%d pvi01=%lu,%lu nvi01=%lu,%lu pvj01=%lu,%lu nvj01=%lu,%lu\n",
        bx,by,tx,ty,aThread,step_2way,nvleD2,kE0,kE1,pvi0,pvi1,nvi0,nvi1,pvj0,pvj1,nvj0,nvj1);

    if(bx==px && by==py && tx<ptx && ty<pty && l<pf)
      printf("b=%d,%d t=%d,%d l=%d pnvi=%lu,%lu pnvj=%lu,%lu pvi0=%d,%d pvi1=%d,%d nvi0=%d,%d nvi1=%d,%d pvi0=%d,%d pvi1=%d,%d nvi0=%d,%d nvi1=%d,%d\n",
        bx,by,tx,ty,l,pvi,nvi,pvj,nvj,((int32_t*)&pvi0)[0],((int32_t*)&pvi0)[1],((int32_t*)&pvi1)[0],((int32_t*)&pvi1)[1],
	((int32_t*)&nvi0)[0],((int32_t*)&nvi0)[1],((int32_t*)&nvi1)[0],((int32_t*)&nvi1)[1],
	((int32_t*)&pvj0)[0],((int32_t*)&pvj0)[1],((int32_t*)&pvj1)[0],((int32_t*)&pvj1)[1],
        ((int32_t*)&nvj0)[0],((int32_t*)&nvj0)[1],((int32_t*)&nvj1)[0],((int32_t*)&nvj1)[1]);

    const uint64_t r00 = gm_popcount64(nvi & nvj);
    const uint64_t r01 = gm_popcount64(nvi & pvj);
    const uint64_t r10 = gm_popcount64(pvi & nvj);
    const uint64_t r11 = gm_popcount64(pvi & pvj);

    //---Accumulate---
    c0 += r00; c1 += r01;
    c2 += r10; c3 += r11;

    if(bx==px && by==py && tx<ptx && ty<pty && l<pf)
      printf("b=%d,%d t=%d,%d a=%d b=%d r00=%lu r01=%lu r10=%lu r11=%lu c0123=%d,%d,%d,%d\n",
             bx,by,tx,ty,aInd,bInd,r00,r01,r10,r11,c0,c1,c2,c3);
  }

  // Different ordering from original 1-bit routines
  // Each thread writes one element of block sub-matrix to memory
  // Assume c is row major

  // Old attempts
  /*int cBegin = n*bx*BLOCK_SIZE*2+by*BLOCK_SIZE*2;
  int rind1 = ty * BLOCK_SIZE*4;
  int rind2 = ty * BLOCK_SIZE*4 + BLOCK_SIZE*2;
  int cind = (tx % 4)*4 + (tx/4)*2;
  int cInd1 = rind1 + cind;
  int cInd2 = rind2 + cind;*/ 
  //int cBegin = n*2*bx*BLOCK_SIZE*2+by*BLOCK_SIZE*2;
  //int cind = (tx % (n/2))*4 + (tx/(n/2))*2;

  // 2-way working
  /*int cBegin = by*n*2*BLOCK_SIZE*2;
  int rind1 = ty * n*4;
  int rind2 = ty * n*4 + n*2;
  int cind = ((tx + bx*BLOCK_SIZE) % (n/2))*4 + ((tx + bx*BLOCK_SIZE) / (n/2))*2;
  int cInd1 = cBegin + rind1 + cind;
  int cInd2 = cBegin + rind2 + cind;*/

  // 3-way attempt
  int cBegin = bx*m*2*BLOCK_SIZE*2;
  int rind1 = tx * m*4;
  int rind2 = tx * m*4 + m*2;
  int cind = ((ty + by*BLOCK_SIZE) % (m/2))*4 + ((ty + by*BLOCK_SIZE) / (m/2))*2;
  int cInd1 = cBegin + rind1 + cind;
  int cInd2 = cBegin + rind2 + cind;

  if(bx==px && by==py && tx<ptx && ty<pty)
    printf("b=%d,%d t=%d,%d cb=%d ci1=%d=%d,%d ci2=%d=%d,%d c0123=%d,%d,%d,%d\n",
           bx,by,tx,ty,cBegin,cInd1,rind1,cind,cInd2,rind2,cind,c0,c1,c2,c3);

  c[cInd1] = c0 - 512; c[cInd1+1] = c2;
  c[cInd2] = c1; c[cInd2+1] = c3;
}


//-----------------------------------------------------------------------------
/// \brief Perform required GEMM.

template<int TC_METHOD>
static void tc_solve_comet_mult_int_impl(bool is_first, int m, int n, int k, int pfl_min, int step_2way, int nfal,
  const void *matA1, const void *matA2, const void *matB, void* matC, TCBufs& tc_bufs, CEnv& env) {
  COMET_INSIST(matC);
  COMET_INSIST(m >= 0 && n >= 0 && k >= 0);

  if(env.print_details()) printf("In tc_solve_comet_mult_int_impl mnk=%d,%d,%d\n",m,n,k);
  //double tbegin = env.get_cpu_time();

  const bool beta = is_first ? 0 : 1;
  const int threadblockx = BLOCK_SIZE, threadblocky = BLOCK_SIZE;
  int gridblockx = (int)ceil((double)m/threadblockx);
  int gridblocky = (int)ceil((double)n/threadblocky);

  if(env.print_details())
    printf("Launching 1-bit Mult Int GEMM kernel with A1 and A2 mnk=%d,%d,%d beta=%d "
           "gridDim=%d,%d threadDim=%d,%d\n",
           m,n,k,(int)beta,gridblockx,gridblocky,threadblockx,threadblocky);

  switch(env.num_kernel()) {
    // Basic GEMM
    case 175: {
      if(env.print_details()) printf("Calling b1_comet_mult_gemm_gpu_int_simple kernel\n");
      COMET_LAUNCH_KERNEL(b1_comet_mult_gemm_gpu_int_simple,
        dim3(gridblockx, gridblocky, 1),
        dim3(threadblockx, threadblocky, 1), 0, env.stream_compute(),
        m, n, k, pfl_min, step_2way, nfal, (GMBits2x64*)matA1, (GMBits2x64*)matA2,
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

  cudaStreamSynchronize(env.stream_compute());
  printf("Printing matrix info\n");
  COMET_LAUNCH_KERNEL(b1_print_matrix_int,
        dim3(gridblockx, gridblocky, 1),
        dim3(threadblockx, threadblocky, 1), 0, env.stream_compute(),
        m, n, (int32_t*)matC);
  System::accel_last_call_succeeded();

  //env.gemmtime_inc(env.get_cpu_time() - tbegin);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TC_SOLVE_MULT_INT_I_HH_

//-----------------------------------------------------------------------------
