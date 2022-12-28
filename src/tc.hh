//-----------------------------------------------------------------------------
/*!
 * \file   tc.hh
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  CUDA code, primarily for using tensor cores.
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

#ifndef _COMET_TC_HH_
#define _COMET_TC_HH_

#if defined COMET_USE_CUDA
#include "cublas_v2.h"
#elif defined COMET_USE_HIP
#include "rocblas.h"
#endif

#include "env.hh"
#include "histograms.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Blas handle typedefs.

#if defined COMET_USE_CUDA
typedef cublasHandle_t accelblasHandle_t;
#elif defined COMET_USE_HIP
typedef rocblas_handle accelblasHandle_t;
#else
typedef int accelblasHandle_t;
#endif

//-----------------------------------------------------------------------------
/// \brief Helper to track actual values needed from GEMM for 2-way metrics.

class GemmShape2Way {
public:

  GemmShape2Way() {}

  GemmShape2Way(bool do_compute_triang_only, int nvl, NV_t nva,
                int proc_num_I, int proc_num_J)
    : do_compute_triang_only_(do_compute_triang_only)
    , nvl_(nvl)
    , nva_(nva)
    , proc_num_I_(proc_num_I)
    , proc_num_J_(proc_num_J) {}

  void set(bool do_compute_triang_only, int nvl, NV_t nva,
           int proc_num_I, int proc_num_J) {
    do_compute_triang_only_ = do_compute_triang_only;
    nvl_ = nvl;
    nva_ = nva;
    proc_num_I_ = proc_num_I;
    proc_num_J_ = proc_num_J; 
  }

  __host__ __device__
  bool is_inside(int I, int J) const {
    return (I < J || !do_compute_triang_only_) &&
           I + nvl_*static_cast<NV_t>(proc_num_I_) < nva_ &&
           J + nvl_*static_cast<NV_t>(proc_num_J_) < nva_;
  }

private:

  bool do_compute_triang_only_;
  int nvl_;
  NV_t nva_;
  int proc_num_I_;
  int proc_num_J_;
};

//-----------------------------------------------------------------------------
/// \brief Helper to track actual values needed from GEMM for 3-way metrics.

class GemmShape3Way {
public:

  GemmShape3Way() {}

  GemmShape3Way(int I_max, int K_min, int K_max, int nvl, NV_t nva,
                int proc_num_I, int proc_num_J, int proc_num_K)
    : I_max_(I_max)
    , K_min_(K_min)
    , K_max_(K_max)
    , nvl_(nvl)
    , nva_(nva)
    , proc_num_I_(proc_num_I)
    , proc_num_J_(proc_num_J)
    , proc_num_K_(proc_num_K) {}

  void set(int I_max, int K_min, int K_max, int nvl, NV_t nva,
           int proc_num_I, int proc_num_J, int proc_num_K) {
    I_max_ = I_max;
    K_min_ = K_min;
    K_max_ = K_max;
    nvl_ = nvl;
    nva_ = nva;
    proc_num_I_ = proc_num_I;
    proc_num_J_ = proc_num_J;
    proc_num_K_ = proc_num_K;
  }

  __host__ __device__
  bool is_inside(int I, int J, int K) const {
    return I < I_max_ && K >= K_min_ && K < K_max_ &&
           I + nvl_*static_cast<NV_t>(proc_num_I_) < nva_ &&
           J + nvl_*static_cast<NV_t>(proc_num_J_) < nva_ &&
           K + nvl_*static_cast<NV_t>(proc_num_K_) < nva_;
  }

private:

  int I_max_;
  int K_min_;
  int K_max_;
  int nvl_;
  NV_t nva_;
  int proc_num_I_;
  int proc_num_J_;
  int proc_num_K_;
};

//-----------------------------------------------------------------------------
/// \brief GemmShape version for unused/uncollected case.

class GemmShapeNone {
public:

  __host__ __device__ bool is_inside(int I, int J, int K = 0) const {
    return false;
  }
};

//-----------------------------------------------------------------------------
/// \brief Aggregator for all possible kinds of gemm shape.

struct GemmShapes {

  GemmShapes(GemmShape2Way* p)
    : gemm_shape_2way(p)
    , gemm_shape_3way(NULL)
    , gemm_shape_none(NULL) {}

  GemmShapes(GemmShape3Way* p)
    : gemm_shape_2way(NULL)
    , gemm_shape_3way(p)
    , gemm_shape_none(NULL) {}

  GemmShapes(GemmShapeNone* p = 0)
    : gemm_shape_2way(NULL)
    , gemm_shape_3way(NULL)
    , gemm_shape_none(p) {}

  void attach(GemmShape2Way* gemm_shape_2way_) {
    gemm_shape_2way = gemm_shape_2way_;
    gemm_shape_3way = NULL;
    gemm_shape_none = NULL;
  }

  void attach(GemmShape3Way* gemm_shape_3way_) {
    gemm_shape_2way = NULL;
    gemm_shape_3way = gemm_shape_3way_;
    gemm_shape_none = NULL;
  }

  void attach(GemmShapeNone* gemm_shape_none_) {
    gemm_shape_2way = NULL;
    gemm_shape_3way = NULL;
    gemm_shape_none = gemm_shape_none_;
  }

  GemmShape2Way* gemm_shape_2way;
  GemmShape3Way* gemm_shape_3way;
  GemmShapeNone* gemm_shape_none;
};

//-----------------------------------------------------------------------------
/// \brief Struct to hold info for tc buffers.

struct TCBufs {

  TCBufs(CEnv& env);
  ~TCBufs();

#if 0
  static void malloc(
    int num_vector_local,
    int num_packedfield_local,
    TCBufs& tc_bufs,
    CEnv& env);

  static void free(TCBufs& tc_bufs, CEnv& env);
#endif

  void* tc_buf_left;
  void* tc_buf_right;
  size_t tc_buf_size;
  uint32_t* matX_counts;
  size_t matX_counts_size;
  accelblasHandle_t accelblas_handle;

  void allocate(int num_vector_local, int num_packedfield_local);
  void deallocate();

private:

  CEnv& env_;
  bool is_allocated_;
};

//-----------------------------------------------------------------------------

struct TCDebug {
//FIX
// size_t I_base;
// size_t J_base;
// size_t K_base;
};

//-----------------------------------------------------------------------------

void tc_gemm_start(
  int m, int n, int k,
  const void* matA1, int ldda1, const void* matA2, int ldda2,
  const void* matB, int lddb, void* matC, int lddc,
  void* sums_I, void* sums_J, void* sums_K,
  void* counts_I, void* counts_J, void* counts_K,
  int J, int nfal, int step_2way, TCBufs& tc_bufs, Histograms& histograms,
  GemmShapes& gemm_shapes, CEnv& env, TCDebug tc_debug);

size_t tc_gemm_vaxis_divisibility_required(const CEnv& env);
size_t tc_gemm_faxis_divisibility_required(const CEnv& env);
size_t tc_nvl_divisibility_required_for_gemm(const CEnv& env);
size_t tc_gemm_vaxis_size_required(size_t size_requested, const CEnv& env);
size_t tc_gemm_faxis_size_required(size_t size_requested, const CEnv& env);
size_t tc_nvl_size_required_for_gemm(size_t size_requested, const CEnv& env);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TC_HH_

//-----------------------------------------------------------------------------
