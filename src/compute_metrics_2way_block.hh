//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_2way_block.hh
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Calculate metrics, 2-way, for a single block.
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

#ifndef _COMET_COMPUTE_METRICS_2WAY_BLOCK_HH_
#define _COMET_COMPUTE_METRICS_2WAY_BLOCK_HH_

#include "env.hh"
#include "mirrored_buf.hh"
#include "compressed_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"

//=============================================================================

namespace comet {

// Forward declaration.
class MagmaWrapper;

//-----------------------------------------------------------------------------

class ComputeMetrics2WayBlock {

public:

  static void compute_nums_start(
    GMVectors* vectors_left,
    GMVectors* vectors_right,
    GMMetrics* metrics,
    MirroredBuf* vectors_left_buf,
    MirroredBuf* vectors_right_buf,
    MirroredBuf* metrics_buf,
    VectorSums* vector_sums_left,
    VectorSums* vector_sums_right,
    int j_proc,
    bool compute_triang_only,
    MagmaWrapper& magma_wrapper,
    CEnv* env);

  static void compute_nums_wait(
    GMVectors* vectors_left,
    GMVectors* vectors_right,
    GMMetrics* metrics,
    MirroredBuf* vectors_left_buf,
    MirroredBuf* vectors_right_buf,
    MirroredBuf* metrics_buf,
    VectorSums* vector_sums_left,
    VectorSums* vector_sums_right,
    int j_proc,
    bool compute_triang_only,
    CEnv* env);

  static void finalize(
    GMMetrics* metrics,
    //MirroredBuf* metrics_buf,
    CompressedBuf* metrics_buf,
    const VectorSums* const vector_sums_left,
    const VectorSums* const vector_sums_right,
    int j_proc,
    bool do_compute_triang_only,
    CEnv* env);

};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_COMPUTE_METRICS_2WAY_BLOCK_HH_

//-----------------------------------------------------------------------------
