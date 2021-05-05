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
    printf("Launching 1-bit GEMM kernel mnk=%d,%d,%d beta=%d "
           "gridDim=%d,%d threadDim=%d,%d\n",
           m,n,k,(int)beta,gridblockx,gridblocky,threadblockx,threadblocky);

  switch(env.num_kernel()) {
    // Basic GEMM
    case 175: {
      /*if(env.print_details()) printf("Calling b1_comet_mult_gemm_gpu_int_simple kernel\n");
      COMET_LAUNCH_KERNEL(b1_comet_mult_gemm_gpu_int_simple,
        dim3(gridblockx, gridblocky, 1),
        dim3(threadblockx, threadblocky, 1), 0, env.stream_compute(),
        m, n, k, (GMBits2x64*)matA,
        (GMBits2x64*)matB, beta, (int32_t*)matC);*/
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
