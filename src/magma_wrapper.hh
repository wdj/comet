//-----------------------------------------------------------------------------
/*!
 * \file   magma_wrapper.hh
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Interface to generalized linear algebra functions, e.g. MAGMA.
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

#ifndef _COMET_MAGMA_WRAPPER_HH_
#define _COMET_MAGMA_WRAPPER_HH_

#include "env.hh"
#include "mirrored_buf.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Helpers

static bool use_minproduct(CEnv* env) {
  return env->metric_type() == MetricType::CZEK;
}

static bool use_mgemm2(CEnv* env) {
  return env->metric_type() == MetricType::CCC &&
         env->num_way() == NumWay::_2 && ! env->sparse();
}

static bool use_mgemm3(CEnv* env) {
  return env->metric_type() == MetricType::CCC &&
         env->num_way() == NumWay::_3 && ! env->sparse();
}

static bool use_mgemm4(CEnv* env) {
  return env->metric_type() == MetricType::CCC && env->sparse();
}

static bool use_mgemm5(CEnv* env) {
  return env->metric_type() == MetricType::DUO;
}

//-----------------------------------------------------------------------------

void gm_linalg_initialize(CEnv* env);

void gm_linalg_finalize(CEnv* env);

//----------

void gm_linalg_malloc(MirroredBuf* p, size_t dim0, size_t dim1, CEnv* env);

void gm_linalg_free(MirroredBuf* p, CEnv* env);

//----------

void gm_linalg_set_matrix_zero_start_(MirroredBuf* matrix_buf, CEnv* env);

void gm_linalg_gemm_magma_start(size_t m, size_t n, size_t k,
  const void* matA, size_t ldda, const void* matB, size_t lddb,
  void* matC, size_t lddc, CEnv* env);

//----------

void gm_linalg_set_matrix_start(MirroredBuf* matrix_buf, CEnv* env);

void gm_linalg_set_matrix_wait(CEnv* env);

void gm_linalg_get_matrix_start(MirroredBuf* matrix_buf, CEnv* env);

void gm_linalg_get_matrix_wait(CEnv* env);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_MAGMA_WRAPPER_HH_

//-----------------------------------------------------------------------------
