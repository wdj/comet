//-----------------------------------------------------------------------------
/*!
 * \file   linalg.hh
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Perform (actual or modified) GEMM operations.
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

#ifndef _COMET_LINALG_HH_
#define _COMET_LINALG_HH_

#include "env.hh"
#include "decomp_mgr.hh"
#include "mirrored_buf.hh"
#include "magma_wrapper.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Class with static members to implement linalg package.

struct LinAlg {

  static void gemm_start(
    size_t m, size_t n, size_t k,
    const MirroredBuf* matA1, const MirroredBuf* matA2,
    const MirroredBuf* matB, MirroredBuf* matC,
    MirroredBuf* sums_I, MirroredBuf* sums_J, MirroredBuf* sums_K,
    MirroredBuf* counts_I, MirroredBuf* counts_J, MirroredBuf* counts_K, int J,
    int step_2way, GMDecompMgr& dm, MagmaWrapper& magma_wrapper,
    GemmShapes& gemm_shapes, CEnv& env, TCDebug tc_debug = {});

  static void gemm_wait(
    size_t m, size_t n, size_t k,
    const MirroredBuf* matA1, const MirroredBuf* matA2,
    const MirroredBuf* matB, MirroredBuf* matC,
    MirroredBuf* sums_I, MirroredBuf* sums_J, MirroredBuf* sums_K,
    MirroredBuf* counts_I, MirroredBuf* counts_J, MirroredBuf* counts_K, int J,
    int step_2way, GMDecompMgr& dm, GemmShapes& gemm_shapes, CEnv& env,
    TCDebug tc_debug = {});

  static void gemm_start(
    size_t m, size_t n, size_t k,
    const MirroredBuf* matA, const MirroredBuf* matB, MirroredBuf* matC,
    MirroredBuf* sums_I, MirroredBuf* sums_J,
    MirroredBuf* counts_I, MirroredBuf* counts_J,
    GMDecompMgr& dm, MagmaWrapper& magma_wrapper,
    GemmShapes& gemm_shapes, CEnv& env, TCDebug tc_debug = {});

  static void gemm_wait(
    size_t m, size_t n, size_t k,
    const MirroredBuf* matA, const MirroredBuf* matB, MirroredBuf* matC,
    MirroredBuf* sums_I, MirroredBuf* sums_J,
    MirroredBuf* counts_I, MirroredBuf* counts_J,
    GMDecompMgr& dm, GemmShapes& gemm_shapes, CEnv& env, TCDebug tc_debug = {});

  static void gemm(
    size_t m, size_t n, size_t k,
    const MirroredBuf* matA, const MirroredBuf* matB, MirroredBuf* matC,
    MirroredBuf* sums_I, MirroredBuf* sums_J,
    MirroredBuf* counts_I, MirroredBuf* counts_J,
    GMDecompMgr& dm, MagmaWrapper& magma_wrapper,
    GemmShapes& gemm_shapes, CEnv& env, TCDebug tc_debug = {});

};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

//#include "linalg.i.hh"

//-----------------------------------------------------------------------------

#endif // _COMET_LINALG_HH_

//-----------------------------------------------------------------------------
