#ifndef _COMET_TC_SOLVE_FULL_I_HH_
#define _COMET_TC_SOLVE_FULL_I_HH_

// Includes
#include "cstdlib"
#include <stdlib.h>
#include <cuda_runtime.h>
#include "cuda_fp16.h"

//#include "env.hh"
//#include "tc.hh"
//#include "tc_helpers.i.hh"
//#include "tc_in.i.hh"
//#include "tc_solve.i.hh"
//#include "tc_out.i.hh"
#include <complex>

// Defines
#define gm_popcount64(x) __popcll(x)
#define BLOCK_SIZE 8

// Typedefs
typedef std::complex<double> Float_t;
typedef struct { double data[2]; } WGMTally2x2;
typedef unsigned long long int WGMUInt64;
typedef WGMUInt64 WGMBits1_2x64;
typedef struct { WGMBits1_2x64 data[2]; } WGMBits2x64;

namespace comet {

//-----------------------------------------------------------------------------
/////// \brief GPU kernel for custom 1-bit tensor core WMMA GEMM
__global__ void b1_xor_gemm_gpu_wmma_simple(int m, int n, int k, WGMBits2x64* a,
                                            WGMBits2x64* b, bool beta, WGMTally2x2* c) {
  // Block and thread indices
  int tx = threadIdx.x, ty = threadIdx.y;
  int bx = blockIdx.x, by = blockIdx.y;

  int gridx = bx*BLOCK_SIZE + tx;
  int gridy = by*BLOCK_SIZE + ty;

  //if(bx==0 && by==0 && tx==0 && ty==0) printf("In b1_xor_gemm_gpu_wmma_simple mnk=%d,%d,%d a=%dx%d * b=%dx%d = c=%dx%d gxy=%d,%d\n",m,n,k,m,k,k,n,m,n,gridx,gridy);

  if(gridx>=m || gridy>=n) return;

  enum { GM_TALLY1_MAX_VALUE_BITS = 26 };

  // Matrix block location
  int aBegin = k * BLOCK_SIZE * bx;
  int bBegin = k * BLOCK_SIZE * by;

  // Stores element of block sub-matrix computed by thread
  double c0 = 0, c1 = 0;

  // Each thread computes one element of block sub-matrix
  for (int l=0; l<k; ++l) {

    // A row major and B col major
    int aInd = aBegin + k*tx + l;
    int bInd = bBegin + k*ty + l;

    /*---Extract input values to process---*/
    const WGMBits2x64 vi = a[aInd];
    const WGMBits2x64 vj = b[bInd];

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
    const WGMUInt64 vi0 = vi.data[0];
    const WGMUInt64 vi1 = vi.data[1];
    const WGMUInt64 vj0 = vj.data[0];
    const WGMUInt64 vj1 = vj.data[1];

    //if(tx==0 && ty==0)
    //printf("b=%d,%d t=%d,%d g=%d,%d a=%d=%d b=%d=%d mnk=%d,%d,%d vi0=%lu vi1=%lu vj0=%lu vj1=%lu\n",
    //       bx,by,tx,ty,gridDim.x,gridDim.y,aBegin,aInd,bBegin,bInd,m,n,k,vi0,vi1,vj0,vj1);

    /*---Get mask to ignore vi seminibbles with value of 1,0---*/
    /*---NOTE: check that this handles pad properly---*/
    const WGMUInt64 oddbits = 0x5555555555555555;

    const WGMUInt64 vij0_10mask = (vi0 | ~(vi0 >> 1)) &
                                  (vj0 | ~(vj0 >> 1)) & oddbits;
    const WGMUInt64 vij1_10mask = (vi1 | ~(vi1 >> 1)) &
                                  (vj1 | ~(vj1 >> 1)) & oddbits;

    /*---Get even, odd bits for each semi-nibble, then mask---*/
    const WGMUInt64 vi0_0 =    vi0 & vij0_10mask;
    const WGMUInt64 vi1_0 =    vi1 & vij1_10mask;
    const WGMUInt64 vj0_0 =    vj0 & vij0_10mask;
    const WGMUInt64 vj1_0 =    vj1 & vij1_10mask;

    /*---Get complements of the same bits, then mask---*/
    const WGMUInt64 nvi0_0 = ~ vi0 & vij0_10mask;
    const WGMUInt64 nvi1_0 = ~ vi1 & vij1_10mask;
    const WGMUInt64 nvj0_0 = ~ vj0 & vij0_10mask;
    const WGMUInt64 nvj1_0 = ~ vj1 & vij1_10mask;

    const WGMUInt64 r00 = gm_popcount64((nvi0_0 & nvj0_0) |
                                      ( (nvi1_0 & nvj1_0) << 1 ));
    const WGMUInt64 r01 = gm_popcount64((nvi0_0 &  vj0_0) |
                                      ( (nvi1_0 &  vj1_0) << 1 ));
    const WGMUInt64 r10 = gm_popcount64(( vi0_0 & nvj0_0) |
                                      ( ( vi1_0 & nvj1_0) << 1 ));
    const WGMUInt64 r11 = gm_popcount64(( vi0_0 &  vj0_0) |
                                      ( ( vi1_0 &  vj1_0) << 1 ));

    /*---Accumulate---*/
    c0 += r00 | (r01 << GM_TALLY1_MAX_VALUE_BITS);
    c1 += r10 | (r11 << GM_TALLY1_MAX_VALUE_BITS);

    //printf("b=%d,%d t=%d,%d a=%d b=%d r00=%ld r01=%ld r10=%ld r11=%ld sum0=%lf sum1=%lf\n",
    //       bx,by,tx,ty,aInd,bInd,r00,r01,r10,r11,c0,c1);
  }

  // Each thread writes one element of block sub-matrix to memory
  // Assume c is row major
  int cBegin = n*bx*BLOCK_SIZE+by*BLOCK_SIZE;
  int cInd   = cBegin + tx*n + ty;
  //printf("b=%d,%d t=%d,%d c=%d\n",bx,by,tx,ty,cInd);
  c[cInd].data[0] = c0;
  c[cInd].data[1] = c1;
}

//-----------------------------------------------------------------------------
///// \brief Run custom 1-bit tensor core WMMA GEMM
void tc_solve_wmma_(size_t m, size_t n, size_t k,
  const void* matA, size_t ldda, const void* matB, size_t lddb,
  void* matC, size_t lddc, CEnv& env) {

  if(env.print_details()) printf("In tc_solve_wmma_\n");

  const bool beta = 1;
  const int threadblockx = BLOCK_SIZE, threadblocky = BLOCK_SIZE;
  int gridblockx = (int)ceil((double)m/threadblockx);
  int gridblocky = (int)ceil((double)n/threadblocky);

  if(env.print_details()) printf("Launching kernel with mnk=%zu,%zu,%zu threads=%d,%d grid=%d,%d\n",m,n,k,
         threadblockx,threadblocky,gridblockx,gridblocky);

  double tbegin = env.synced_time();

  // Assuming A and B are col oriented, C is row oriented
  COMET_LAUNCH_KERNEL(b1_xor_gemm_gpu_wmma_simple,
    dim3(gridblockx, gridblocky, 1),
    dim3(threadblockx, threadblocky, 1), 0, env.stream_compute(),
    m, n, k, (WGMBits2x64*)matA,
    (WGMBits2x64*)matB, beta, (WGMTally2x2*)matC);
  //cudaDeviceSynchronize();

  env.ops_local_inc(2 * m * (double)n * (double)k * 16); // 16 bytes per doublecomplex
  env.gemmtime_inc(env.synced_time() - tbegin);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif

