//-----------------------------------------------------------------------------
/*!
 * \file   linalg.cc
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Interface to generalized linear algebra functions, e.g. MAGMA.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifdef COMET_USE_MAGMA
#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"
#include "magma_mgemm2.h"
#include "magma_mgemm2_lapack.h"
#include "magma_mgemm3.h"
#include "magma_mgemm3_lapack.h"
#include "magma_mgemm4.h"
#include "magma_mgemm4_lapack.h"
#include "magma_mgemm5.h"
#include "magma_mgemm5_lapack.h"
//#elif defined COMET_USE_CUDA
//  #include "cublas_v2.h"
//#elif defined COMET_USE_HIP
//  #include "hip/hip_runtime_api.h"
//  #include "hip/hip_runtime.h"
#endif

#include "env.hh"
#include "assertions.hh"
#include "tc.hh"
#include "decomp_mgr.hh"
#include "linalg.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Helpers

static bool use_minproduct(CEnv* env) {
  return env->metric_type() == MetricType::CZEK;
}

static bool use_mgemm2(CEnv* env) {
  return env->metric_type() == MetricType::CCC &&
         env->num_way() == NUM_WAY::_2 && ! env->sparse();
}

static bool use_mgemm3(CEnv* env) {
  return env->metric_type() == MetricType::CCC &&
         env->num_way() == NUM_WAY::_3 && ! env->sparse();
}

static bool use_mgemm4(CEnv* env) {
  return env->metric_type() == MetricType::CCC && env->sparse();
}

static bool use_mgemm5(CEnv* env) {
  return env->metric_type() == MetricType::DUO;
}

//=============================================================================
/*---Magma setup, teardown---*/

void gm_linalg_initialize(CEnv* env) {
  COMET_INSIST(env);

  // need magma blasSetKernelStream -- see
  // http://on-demand.gputechconf.com/gtc/2014/presentations/S4158-cuda-streams-best-practices-common-pitfalls.pdf
  // page 14

// TODO: non-GPU calls should nt need to init magma, should use
// regular malloc instead of magma malloc.

#ifdef COMET_USE_MAGMA
  if (use_minproduct(env)) { //--------------------

    magma_minproduct_int_t magma_code = magma_minproduct_init();
    COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in call to magma_minproduct_init.");
    magma_code = magma_minproductblasSetKernelStream(env->stream_compute());
    COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in call to magma_minproductblasSetKernelStream.");

  } else if (use_mgemm4(env)) { //--------------------

    magma_mgemm4_int_t magma_code = magma_mgemm4_init();
    COMET_INSIST(magma_code == MAGMA_mgemm4_SUCCESS &&
                   "Error in call to magma_mgemm4_init.");
    magma_code = magma_mgemm4blasSetKernelStream(env->stream_compute());
    COMET_INSIST(magma_code == MAGMA_mgemm4_SUCCESS &&
                   "Error in call to magma_mgemm4blasSetKernelStream.");

  } else if (use_mgemm2(env)) { //--------------------

    magma_mgemm2_int_t magma_code = magma_mgemm2_init();
    COMET_INSIST(magma_code == MAGMA_mgemm2_SUCCESS &&
                   "Error in call to magma_mgemm2_init.");
    magma_code = magma_mgemm2blasSetKernelStream(env->stream_compute());
    COMET_INSIST(magma_code == MAGMA_mgemm2_SUCCESS &&
                   "Error in call to magma_mgemm2blasSetKernelStream.");

  } else if (use_mgemm3(env)) { //--------------------

    magma_mgemm3_int_t magma_code = magma_mgemm3_init();
    COMET_INSIST(magma_code == MAGMA_mgemm3_SUCCESS &&
                   "Error in call to magma_mgemm3_init.");
    magma_code = magma_mgemm3blasSetKernelStream(env->stream_compute());
    COMET_INSIST(magma_code == MAGMA_mgemm3_SUCCESS &&
                   "Error in call to magma_mgemm3blasSetKernelStream.");

  } else if (use_mgemm5(env)) { //--------------------

    magma_mgemm5_int_t magma_code = magma_mgemm5_init();
    COMET_INSIST(magma_code == MAGMA_mgemm5_SUCCESS &&
                   "Error in call to magma_mgemm5_init.");
    magma_code = magma_mgemm5blasSetKernelStream(env->stream_compute());
    COMET_INSIST(magma_code == MAGMA_mgemm5_SUCCESS &&
                   "Error in call to magma_mgemm5blasSetKernelStream.");

  } else { //--------------------

      COMET_INSIST_INTERFACE(env, false && "Unimplemented modified gemm method.");

  } // if //--------------------
#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------

void gm_linalg_finalize(CEnv* env) {
  COMET_INSIST(env);

  // TODO: (maybe) reset kernel stream (probably not really needed)

#ifdef COMET_USE_MAGMA
  if (use_minproduct(env)) { //--------------------

    magma_minproduct_int_t magma_code = magma_minproduct_finalize();
    COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in call to magma_minproduct_finalize.");

  } else if (use_mgemm4(env)) { //--------------------

    magma_mgemm4_int_t magma_code = magma_mgemm4_finalize();
    COMET_INSIST(magma_code == MAGMA_mgemm4_SUCCESS &&
                   "Error in call to magma_mgemm4_finalize.");

  } else if (use_mgemm2(env)) { //--------------------

    magma_mgemm2_int_t magma_code = magma_mgemm2_finalize();
    COMET_INSIST(magma_code == MAGMA_mgemm2_SUCCESS &&
                   "Error in call to magma_mgemm2_finalize.");

  } else if (use_mgemm3(env)) { //--------------------

    magma_mgemm3_int_t magma_code = magma_mgemm3_finalize();
    COMET_INSIST(magma_code == MAGMA_mgemm3_SUCCESS &&
                   "Error in call to magma_mgemm3_finalize.");

  } else if (use_mgemm5(env)) { //--------------------

    magma_mgemm5_int_t magma_code = magma_mgemm5_finalize();
    COMET_INSIST(magma_code == MAGMA_mgemm5_SUCCESS &&
                   "Error in call to magma_mgemm5_finalize.");

  } else { //--------------------

      COMET_INSIST_INTERFACE(env, false && "Unimplemented modified gemm method.");

  } // if //--------------------
#endif // COMET_USE_MAGMA
}

//=============================================================================
/*---Allocate/free host and device memory---*/

void gm_linalg_malloc(MirroredBuf* p, size_t dim0, size_t dim1, CEnv* env) {
  COMET_INSIST(p && env);
  COMET_INSIST(dim0 + 1 >= 1 && dim1 + 1 >= 1);

#ifdef COMET_USE_MAGMA

  p->dim0 = dim0;
  p->dim1 = dim1;

  const size_t n = dim0 * dim1;

  if (use_minproduct(env)) { //--------------------

    typedef GMFloat Float_t;

    magma_minproduct_int_t magma_code = 0;

    if (FP_PRECISION_DOUBLE) {
      magma_code = magma_minproduct_dmalloc_pinned((double**)&p->h, n);
      COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in magma_minproduct_dmalloc_pinned,"
                   " possibly insufficient memory.");
    } else {
      magma_code = magma_minproduct_smalloc_pinned((float**)&p->h, n);
      COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in magma_minproduct_smalloc_pinned,"
                   " possibly insufficient memory.");
    }
    GMFloat_fill_nan((GMFloat*)p->h, n);

    if (env->is_compute_method_gpu()) {
      if (FP_PRECISION_DOUBLE) {
        magma_code = magma_minproduct_dmalloc((double**)&p->d, n);
        COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in magma_minproduct_dmalloc,"
                     " possibly insufficient memory.");
      } else {
        magma_code = magma_minproduct_smalloc((float**)&p->d, n);
        COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS &&
                     "Error in call to magma_minproduct_smalloc,"
                     " possibly insufficient memory.");
      }
    }
    // TODO: ? fill GPU memory with NaNs
    p->size = n * sizeof(Float_t);

  } else if (use_mgemm4(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    magma_mgemm4_int_t magma_code = 0;

    magma_code = magma_mgemm4_zmalloc_pinned((Float_t**)&p->h, n);
    COMET_INSIST(magma_code == MAGMA_mgemm4_SUCCESS &&
      "Error in magma_mgemm4_zmalloc_pinned, possibly insufficient memory.");

    if (env->is_compute_method_gpu()) {
      magma_code = magma_mgemm4_zmalloc((Float_t**)&p->d, n);
      COMET_INSIST(magma_code == MAGMA_mgemm4_SUCCESS &&
        "Error in magma_mgemm4_zmalloc, possibly insufficient memory.");
    }
    p->size = n * sizeof(Float_t);

  } else if (use_mgemm2(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    magma_mgemm2_int_t magma_code = 0;

    magma_code = magma_mgemm2_zmalloc_pinned((Float_t**)&p->h, n);
    COMET_INSIST(magma_code == MAGMA_mgemm2_SUCCESS &&
      "Error in magma_mgemm2_zmalloc_pinned, possibly insufficient memory.");

    if (env->is_compute_method_gpu()) {
      magma_code = magma_mgemm2_zmalloc((Float_t**)&p->d, n);
      COMET_INSIST(magma_code == MAGMA_mgemm2_SUCCESS &&
        "Error in magma_mgemm2_zmalloc, possibly insufficient memory.");
    }
    p->size = n * sizeof(Float_t);

  } else if (use_mgemm3(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    magma_mgemm3_int_t magma_code = 0;

    magma_code = magma_mgemm3_zmalloc_pinned((Float_t**)&p->h, n);
    COMET_INSIST(magma_code == MAGMA_mgemm3_SUCCESS &&
      "Error in magma_mgemm3_zmalloc_pinned, possibly insufficient memory.");

    if (env->is_compute_method_gpu()) {
      magma_code = magma_mgemm3_zmalloc((Float_t**)&p->d, n);
      COMET_INSIST(magma_code == MAGMA_mgemm3_SUCCESS &&
        "Error in magma_mgemm3_zmalloc, possibly insufficient memory.");
    }
    p->size = n * sizeof(Float_t);

  } else if (use_mgemm5(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    magma_mgemm5_int_t magma_code = 0;

    magma_code = magma_mgemm5_zmalloc_pinned((Float_t**)&p->h, n);
    COMET_INSIST(magma_code == MAGMA_mgemm5_SUCCESS &&
      "Error in magma_mgemm5_zmalloc_pinned, possibly insufficient memory.");

    if (env->is_compute_method_gpu()) {
      magma_code = magma_mgemm5_zmalloc((Float_t**)&p->d, n);
      COMET_INSIST(magma_code == MAGMA_mgemm5_SUCCESS &&
        "Error in magma_mgemm5_zmalloc, possibly insufficient memory.");
    }
    p->size = n * sizeof(Float_t);

  } else { //--------------------

      COMET_INSIST_INTERFACE(env, false && "Unimplemented modified gemm method.");

  } // if //--------------------

  env->cpu_mem_local_inc(p->size);
  if (env->is_compute_method_gpu())
    env->gpu_mem_local_inc(p->size);

  COMET_INSIST(p->h &&
    "Invalid host pointer created, possibly due to insufficient memory.");
  COMET_INSIST((p->d || !env->is_compute_method_gpu()) &&
    "Invalid device pointer created, possibly due to insufficient memory.");

#else

  COMET_INSIST(false && "Magma library unavailable for this build.");

#endif
}

//-----------------------------------------------------------------------------

void gm_linalg_free(MirroredBuf* p, CEnv* env) {
  COMET_INSIST(p && env);
  COMET_INSIST(! p->is_alias);

#ifdef COMET_USE_MAGMA
  if (use_minproduct(env)) { //--------------------

    magma_minproduct_int_t magma_code = magma_minproduct_free_pinned(p->h);
    COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS &&
             "Error in magma_minproduct_free_pinned.");
    if (env->is_compute_method_gpu()) {
      magma_code = magma_minproduct_free(p->d);
      COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS &&
             "Error in magma_minproduct_free.");
    }

  } else if (use_mgemm4(env)) { //--------------------

    magma_mgemm4_int_t magma_code = magma_mgemm4_free_pinned(p->h);
    COMET_INSIST(magma_code == MAGMA_mgemm4_SUCCESS &&
             "Error in magma_mgemm4_free_pinned.");
    if (env->is_compute_method_gpu()) {
      magma_code = magma_mgemm4_free(p->d);
      COMET_INSIST(magma_code == MAGMA_mgemm4_SUCCESS &&
               "Error in magma_mgemm4_free.");
    }

  } else if (use_mgemm2(env)) { //--------------------

    magma_mgemm2_int_t magma_code = magma_mgemm2_free_pinned(p->h);
    COMET_INSIST(magma_code == MAGMA_mgemm2_SUCCESS &&
             "Error in magma_mgemm2_free_pinned.");
    if (env->is_compute_method_gpu()) {
      magma_code = magma_mgemm2_free(p->d);
      COMET_INSIST(magma_code == MAGMA_mgemm2_SUCCESS &&
               "Error in magma_mgemm2_free.");
    }

  } else if (use_mgemm3(env)) { //--------------------

    magma_mgemm3_int_t magma_code = magma_mgemm3_free_pinned(p->h);
    COMET_INSIST(magma_code == MAGMA_mgemm3_SUCCESS &&
             "Error in magma_mgemm3_free_pinned.");
    if (env->is_compute_method_gpu()) {
      magma_code = magma_mgemm3_free(p->d);
      COMET_INSIST(magma_code == MAGMA_mgemm3_SUCCESS &&
               "Error in magma_mgemm3_free.");
    }

  } else if (use_mgemm5(env)) { //--------------------

    magma_mgemm5_int_t magma_code = magma_mgemm5_free_pinned(p->h);
    COMET_INSIST(magma_code == MAGMA_mgemm5_SUCCESS &&
             "Error in magma_mgemm5_free_pinned.");
    if (env->is_compute_method_gpu()) {
      magma_code = magma_mgemm5_free(p->d);
      COMET_INSIST(magma_code == MAGMA_mgemm5_SUCCESS &&
               "Error in magma_mgemm5_free.");
    }

  } else { //--------------------

      COMET_INSIST_INTERFACE(env, false && "Unimplemented modified gemm method.");

  } // if //--------------------

  p->h = NULL;
  p->d = NULL;

  env->cpu_mem_local_dec(p->size);
  if (env->is_compute_method_gpu())
    env->gpu_mem_local_dec(p->size);

#else

  COMET_INSIST(false && "Magma library unavailable for this build.");

#endif
}

//-----------------------------------------------------------------------------

void gm_linalg_set_matrix_zero_start_(MirroredBuf* matrix_buf, CEnv* env) {
  COMET_INSIST(matrix_buf && env);

  if (!env->is_compute_method_gpu()) {
//    memset(matrix_buf->h, 0, matrix_buf->size);
    return;
  }

#ifdef COMET_USE_MAGMA

  const size_t mat_dim1 = matrix_buf->dim0;
  const size_t mat_dim2 = matrix_buf->dim1;

  // ISSUE: these MAGMA routines don't return an error code.

  if (use_minproduct(env)) { //--------------------

    if (FP_PRECISION_DOUBLE) {
      magma_minproductblas_dlaset
        (Magma_minproductFull, mat_dim1, mat_dim2, (double)0, (double)0,
         (double*)matrix_buf->d, mat_dim1);
    } else {
      magma_minproductblas_slaset
        (Magma_minproductFull, mat_dim1, mat_dim2, (float)0, (float)0,
         (float*)matrix_buf->d, mat_dim1);
    }

  } else if (use_mgemm4(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    Float_t zero = {0, 0};

    magma_mgemm4blas_zlaset(Magma_mgemm4Full, mat_dim1, mat_dim2, zero, zero,
                            (Float_t*)matrix_buf->d, mat_dim1);

  } else if (use_mgemm2(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    Float_t zero = {0, 0};

    magma_mgemm2blas_zlaset(Magma_mgemm2Full, mat_dim1, mat_dim2, zero, zero,
                            (Float_t*)matrix_buf->d, mat_dim1);

  } else if (use_mgemm3(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    Float_t zero = {0, 0};

    magma_mgemm3blas_zlaset(Magma_mgemm3Full, mat_dim1, mat_dim2, zero, zero,
                            (Float_t*)matrix_buf->d, mat_dim1);

  } else if (use_mgemm5(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    Float_t zero = {0, 0};

    magma_mgemm5blas_zlaset(Magma_mgemm5Full, mat_dim1, mat_dim2, zero, zero,
                            (Float_t*)matrix_buf->d, mat_dim1);

  } else { //--------------------

      COMET_INSIST_INTERFACE(env, false && "Unimplemented modified gemm method.");

  } // if //--------------------

#else

  COMET_INSIST(false && "Magma library unavailable for this build.");

//#elif defined COMET_USE_CUDA
// NOTE: shouldn't be needed.
//  cudaMemset(matrix_buf->d, 0, matrix_buf->size);
//#elif defined COMET_USE_HIP
// NOTE: shouldn't be needed.
//  hipMemset(matrix_buf->d, 0, matrix_buf->size);
#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------

void gm_linalg_gemm_magma_block_start(size_t m,
                                      size_t n,
                                      size_t k,
                                      const void* matA,
                                      size_t ldda,
                                      const void* matB,
                                      size_t lddb,
                                      void* matC,
                                      size_t lddc,
                                      bool is_beta_one,
                                      CEnv* env) {
  COMET_INSIST(matA && matB && matC && env);
  COMET_INSIST(env->is_compute_method_gpu());

#ifdef COMET_USE_MAGMA
  {
    // Ensure Magmablas function doesn't internally failover to CUBLAS call.
    int TransA = 1;
    int TransB = 0;

    size_t Am = ! TransA ? m : k;
    size_t An = ! TransA ? k : m;
    size_t Bm = ! TransB ? k : n;
    size_t Bn = ! TransB ? n : k;
    size_t sizeA = ldda * (An - 1) + Am;
    size_t sizeB = lddb * (Bn - 1) + Bm;

    size_t CUBLAS_MAX_1DBUF_SIZE = ((1 << 27) - 512);
    COMET_INSIST((! (sizeA >= CUBLAS_MAX_1DBUF_SIZE ||
                 sizeB >= CUBLAS_MAX_1DBUF_SIZE )) &&
             "Error in MAGMA block sizes.");
  }

  // ISSUE: these MAGMA routines don't return an error code.

  if (use_minproduct(env)) { //--------------------

    const GMFloat alpha = 1;
    const GMFloat beta = is_beta_one ? 1 : 0;

    typedef magma_minproduct_int_t Int_t;

    const Int_t m_ = m;
    const Int_t n_ = n;
    const Int_t k_ = k;
    const Int_t ldda_ = ldda;
    const Int_t lddb_ = lddb;
    const Int_t lddc_ = lddc;
    COMET_INSIST((size_t)m_ == m && "Integer overflow.");
    COMET_INSIST((size_t)n_ == n && "Integer overflow.");
    COMET_INSIST((size_t)k_ == k && "Integer overflow.");
    COMET_INSIST((size_t)ldda_ == ldda && "Integer overflow.");
    COMET_INSIST((size_t)lddb_ == lddb && "Integer overflow.");
    COMET_INSIST((size_t)lddc_ == lddc && "Integer overflow.");

    if (FP_PRECISION_DOUBLE) {
      magma_minproductblas_dgemm(
        Magma_minproductTrans,
        Magma_minproductNoTrans,
        m_,
        n_,
        k_,
        alpha,
        (double*)matA,
        ldda_,
        (double*)matB,
        lddb_,
        beta,
        (double*)matC,
        lddc_);
      COMET_INSIST(System::accel_last_call_succeeded() &&
               "Failure in call to magma_minproductblas_dgemm.");
    } else {
      magma_minproductblas_sgemm(
        Magma_minproductTrans,
        Magma_minproductNoTrans,
        m_,
        n_,
        k_,
        alpha,
        (float*)matA,
        ldda_,
        (float*)matB,
        lddb_,
        beta,
        (float*)matC,
        lddc_);
      COMET_INSIST(System::accel_last_call_succeeded() &&
               "Failure in call to magma_minproductblas_sgemm.");
    }

    env->ops_local_inc(2 * m * (double)n * (double)k);

  } else if (use_mgemm4(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const Float_t alpha = {1, 0};
    const Float_t beta = is_beta_one ? one : zero;

    typedef magma_mgemm4_int_t Int_t;

    const Int_t m_ = m;
    const Int_t n_ = n;
    const Int_t k_ = k;
    const Int_t ldda_ = ldda;
    const Int_t lddb_ = lddb;
    const Int_t lddc_ = lddc;
    COMET_INSIST((size_t)m_ == m && "Integer overflow.");
    COMET_INSIST((size_t)n_ == n && "Integer overflow.");
    COMET_INSIST((size_t)k_ == k && "Integer overflow.");
    COMET_INSIST((size_t)ldda_ == ldda && "Integer overflow.");
    COMET_INSIST((size_t)lddb_ == lddb && "Integer overflow.");
    COMET_INSIST((size_t)lddc_ == lddc && "Integer overflow.");

    magma_mgemm4blas_zgemm(
      Magma_mgemm4Trans,
      Magma_mgemm4NoTrans,
      m_,
      n_,
      k_,
      alpha,
      (Float_t*)matA,
      ldda_,
      (Float_t*)matB,
      lddb_,
      beta,
      (Float_t*)matC,
      lddc_);
    COMET_INSIST(System::accel_last_call_succeeded() &&
             "Failure in call to magma_mgemm4blas_zgemm.");

  } else if (use_mgemm2(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const Float_t alpha = {1, 0};
    const Float_t beta = is_beta_one ? one : zero;

    typedef magma_mgemm2_int_t Int_t;

    const Int_t m_ = m;
    const Int_t n_ = n;
    const Int_t k_ = k;
    const Int_t ldda_ = ldda;
    const Int_t lddb_ = lddb;
    const Int_t lddc_ = lddc;
    COMET_INSIST((size_t)m_ == m && "Integer overflow.");
    COMET_INSIST((size_t)n_ == n && "Integer overflow.");
    COMET_INSIST((size_t)k_ == k && "Integer overflow.");
    COMET_INSIST((size_t)ldda_ == ldda && "Integer overflow.");
    COMET_INSIST((size_t)lddb_ == lddb && "Integer overflow.");
    COMET_INSIST((size_t)lddc_ == lddc && "Integer overflow.");

    magma_mgemm2blas_zgemm(
      Magma_mgemm2Trans,
      Magma_mgemm2NoTrans,
      m_,
      n_,
      k_,
      alpha,
      (Float_t*)matA,
      ldda_,
      (Float_t*)matB,
      lddb_,
      beta,
      (Float_t*)matC,
      lddc_);
    COMET_INSIST(System::accel_last_call_succeeded() &&
             "Failure in call to magma_mgemm2blas_zgemm.");

  } else if (use_mgemm3(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const Float_t alpha = {1, 0};
    const Float_t beta = is_beta_one ? one : zero;

    typedef magma_mgemm3_int_t Int_t;

    const Int_t m_ = m;
    const Int_t n_ = n;
    const Int_t k_ = k;
    const Int_t ldda_ = ldda;
    const Int_t lddb_ = lddb;
    const Int_t lddc_ = lddc;
    COMET_INSIST((size_t)m_ == m && "Integer overflow.");
    COMET_INSIST((size_t)n_ == n && "Integer overflow.");
    COMET_INSIST((size_t)k_ == k && "Integer overflow.");
    COMET_INSIST((size_t)ldda_ == ldda && "Integer overflow.");
    COMET_INSIST((size_t)lddb_ == lddb && "Integer overflow.");
    COMET_INSIST((size_t)lddc_ == lddc && "Integer overflow.");

    magma_mgemm3blas_zgemm(
      Magma_mgemm3Trans,
      Magma_mgemm3NoTrans,
      m_,
      n_,
      k_,
      alpha,
      (Float_t*)matA,
      ldda_,
      (Float_t*)matB,
      lddb_,
      beta,
      (Float_t*)matC,
      lddc_);
    COMET_INSIST(System::accel_last_call_succeeded() &&
             "Failure in call to magma_mgemm3blas_zgemm.");

  } else if (use_mgemm5(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    typedef magma_mgemm5_int_t Int_t;

    const Int_t m_ = m;
    const Int_t n_ = n;
    const Int_t k_ = k;
    const Int_t ldda_ = ldda;
    const Int_t lddb_ = lddb;
    const Int_t lddc_ = lddc;
    COMET_INSIST((size_t)m_ == m && "Integer overflow.");
    COMET_INSIST((size_t)n_ == n && "Integer overflow.");
    COMET_INSIST((size_t)k_ == k && "Integer overflow.");
    COMET_INSIST((size_t)ldda_ == ldda && "Integer overflow.");
    COMET_INSIST((size_t)lddb_ == lddb && "Integer overflow.");
    COMET_INSIST((size_t)lddc_ == lddc && "Integer overflow.");

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const Float_t alpha = {1, 0};
    const Float_t beta = is_beta_one ? one : zero;

    magma_mgemm5blas_zgemm(
      Magma_mgemm5Trans,
      Magma_mgemm5NoTrans,
      m_,
      n_,
      k_,
      alpha,
      (Float_t*)matA,
      ldda_,
      (Float_t*)matB,
      lddb_,
      beta,
      (Float_t*)matC,
      lddc_);
    COMET_INSIST(System::accel_last_call_succeeded() &&
             "Failure in call to magma_mgemm5blas_zgemm.");

  } else { //--------------------

      COMET_INSIST_INTERFACE(env, false && "Unimplemented modified gemm method.");

  } // if //--------------------
#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------

void gm_linalg_gemm_magma_start(size_t m,
                                size_t n,
                                size_t k,
                                const void* matA,
                                size_t ldda,
                                const void* matB,
                                size_t lddb,
                                void* matC,
                                size_t lddc,
                                GMDecompMgr* dm,
                                CEnv* env) {
  COMET_INSIST(matA && matB && matC && env);
  COMET_INSIST(env->is_compute_method_gpu());

  // The purpose of this code is to workaround the magma size
  // limitation (for non CUBLAS failover) by doing gemm in blocks.

#ifdef COMET_USE_MAGMA
  const size_t rows = k;
  const size_t cols_A = m;
  const size_t cols_B = n;

  const size_t elt_size =
    env->metric_type() == MetricType::CZEK ? sizeof(GMFloat) :
   (env->metric_type() == MetricType::CCC && env->sparse()) ?
                                         sizeof(magma_mgemm4DoubleComplex) :
   (env->metric_type() == MetricType::CCC &&
    env->num_way() == NUM_WAY::_2) ? sizeof(magma_mgemm2DoubleComplex) :
   (env->metric_type() == MetricType::CCC &&
    env->num_way() == NUM_WAY::_3) ? sizeof(magma_mgemm3DoubleComplex) :
   (env->metric_type() == MetricType::DUO) ?
                                         sizeof(magma_mgemm5DoubleComplex) : 0;
  COMET_INSIST(elt_size > 0 && "Error in gemm block calculation.");

  const size_t align_factor = 128 / elt_size;
  const size_t max_elts = (1 << 27) - 512;

  /*---TODO: can we improve aspect ratios of submatrices---*/
//  const size_t max_rows_per_block_raw = (1 << 14);
//  const size_t max_cols_per_block_raw = max_elts / max_rows_per_block_raw;

  const size_t max_rows_per_block_raw = rows + align_factor;
  const size_t max_cols_per_block_raw = max_elts / rows;

  const size_t max_rows_per_block = (max_rows_per_block_raw / align_factor)
                                                            * align_factor;
  const size_t max_cols_per_block = (max_cols_per_block_raw / align_factor)
                                                            * align_factor;

  COMET_INSIST(max_rows_per_block != 0 && "Error in gemm block calculation.");
  COMET_INSIST(max_cols_per_block != 0 && "Error in gemm block calculation.");

  const size_t cols_per_block_A = utils::min(cols_A, max_cols_per_block);
  const size_t cols_per_block_B = utils::min(cols_B, max_cols_per_block);

  const size_t rows_per_block = utils::min(rows, max_rows_per_block);

  for (size_t row_base=0; row_base<rows; row_base+=rows_per_block) {
    const size_t rows_remaining = rows - row_base;
    const size_t rows_this = utils::min(rows_remaining, rows_per_block);

    for (size_t col_A_base=0; col_A_base<cols_A; col_A_base+=cols_per_block_A) {
      const size_t cols_A_remaining = cols_A - col_A_base;
      const size_t cols_A_this = utils::min(cols_A_remaining, cols_per_block_A);

      void* A_this = (char*)matA + (row_base + ldda*col_A_base)*elt_size;

      for (size_t col_B_base=0; col_B_base<cols_B;
           col_B_base+=cols_per_block_B) {

        const size_t cols_B_remaining = cols_B - col_B_base;
        const size_t cols_B_this = utils::min(cols_B_remaining,
                                             cols_per_block_B);

        void* B_this = (char*)matB + (row_base + ldda*col_B_base)*elt_size;

        void* C_this = (char*)matC + (col_A_base + lddc*col_B_base)*elt_size;

        gm_linalg_gemm_magma_block_start(cols_A_this, cols_B_this, rows_this,
          A_this, ldda, B_this, lddb, C_this, lddc, row_base > 0,  env);
      }
    }
  }
#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------

void gm_linalg_gemm_start(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA1, const MirroredBuf* matA2,
  const MirroredBuf* matB, MirroredBuf* matC,
  int step_2way, GMDecompMgr* dm, CEnv* env) {
  COMET_INSIST(matA1 && matA2 && matB && matC && env);

  if (m==0 || n==0 || k==0)
    return;

  if (env->is_compute_method_gpu()) {
    matA1->lock_d();
    if (matB != matA1) {
      matB->lock_d();
    }
    matC->lock_d();
  }

  if (env->is_using_tc()) {
    if (env->is_compute_method_gpu()) {
      tc_gemm_start(m, n, k,
        matA1->active, matA1->dim0, matA2->active, matA2->dim0,
        matB->active, matB->dim0, matC->active, matC->dim0,
        step_2way, dm->tc_bufs, *env);
    }
  } else {
    gm_linalg_set_matrix_zero_start_(matC, env); // apparently needed by magma.
    gm_linalg_gemm_magma_start(m, n, k, matA1->active, matA1->dim0,
      matB->active, matB->dim0, matC->active, matC->dim0, dm,  env);
  }
}

//-----------------------------------------------------------------------------

void gm_linalg_gemm_wait(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA1, const MirroredBuf* matA2,
  const MirroredBuf* matB, MirroredBuf* matC,
  int step_2way, GMDecompMgr* dm, CEnv* env) {
  COMET_INSIST(matA1 && matA2 && matB && matC && env);

  if (m==0 || n==0 || k==0)
    return;

  if (env->is_using_tc()) {
    if (!env->is_compute_method_gpu()) {
      matA1->lock_h();
      if (matA2 != matA1 && matA2 != matB) {
        matA2->lock_h();
      }
      if (matB != matA1) {
        matB->lock_h();
      }
      matC->lock_h();
      tc_gemm_start(m, n, k,
        matA1->active, matA1->dim0, matA2->active, matA2->dim0,
        matB->active, matB->dim0, matC->active, matC->dim0,
        step_2way, dm->tc_bufs, *env);
      matA1->unlock_h();
      if (matA2 != matA1 && matA2 != matB) {
        matA2->unlock_h();
      }
      if (matB != matA1) {
        matB->unlock_h();
      }
      matC->unlock_h();
    }
  }

  env->stream_synchronize(env->stream_compute());

  if (env->is_compute_method_gpu()) {
    matA1->unlock_d();
    if (matB != matA1) {
      matB->unlock_d();
    }
    matC->unlock_d();
  }
}

//-----------------------------------------------------------------------------

void gm_linalg_gemm_start(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA, const MirroredBuf* matB, MirroredBuf* matC,
  GMDecompMgr* dm, CEnv* env) {

  gm_linalg_gemm_start(m, n, k, matA, matA, matB, matC, 0, dm, env);
}

//-----------------------------------------------------------------------------

void gm_linalg_gemm_wait(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA, const MirroredBuf* matB, MirroredBuf* matC,
  GMDecompMgr* dm, CEnv* env) {

  gm_linalg_gemm_wait(m, n, k, matA, matA, matB, matC, 0, dm, env);
}

//-----------------------------------------------------------------------------

void gm_linalg_gemm(
  size_t m, size_t n, size_t k,
  const MirroredBuf* matA, const MirroredBuf* matB, MirroredBuf* matC,
  GMDecompMgr* dm, CEnv* env) {
  COMET_INSIST(matA && matB && matC && env);

  gm_linalg_gemm_start(m, n, k, matA, matB, matC, dm, env);
  gm_linalg_gemm_wait(m, n, k, matA, matB, matC, dm, env);
}

//=============================================================================
/*---Start/end transfer of generic matrix to GPU---*/

void gm_linalg_set_matrix_start(MirroredBuf* p, CEnv* env) {
  COMET_INSIST(p && env);

  if (!env->is_compute_method_gpu()) {
    return;
  }

  /*---Send vectors to GPU---*/

  // ISSUE: these MAGMA routines don't return an error code.

#ifdef COMET_USE_MAGMA

  if (use_minproduct(env)) { //--------------------

    if (FP_PRECISION_DOUBLE) {
      magma_minproduct_dsetmatrix_async(
        p->dim0, p->dim1, (double*)p->h, p->dim0,
        (double*)p->d, p->dim0, env->stream_togpu());
    } else {
      magma_minproduct_ssetmatrix_async(
        p->dim0, p->dim1, (float*)p->h, p->dim0,
        (float*)p->d, p->dim0, env->stream_togpu());
    }

  } else if (use_mgemm4(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    magma_mgemm4_zsetmatrix_async(p->dim0, p->dim1, (Float_t*)p->h,
                                  p->dim0, (Float_t*)p->d, p->dim0,
                                  env->stream_togpu());

  } else if (use_mgemm2(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    magma_mgemm2_zsetmatrix_async(p->dim0, p->dim1, (Float_t*)p->h,
                                  p->dim0, (Float_t*)p->d, p->dim0,
                                  env->stream_togpu());

  } else if (use_mgemm3(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    magma_mgemm3_zsetmatrix_async(p->dim0, p->dim1, (Float_t*)p->h,
                                  p->dim0, (Float_t*)p->d, p->dim0,
                                  env->stream_togpu());

  } else if (use_mgemm5(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    magma_mgemm5_zsetmatrix_async(p->dim0, p->dim1, (Float_t*)p->h,
                                  p->dim0, (Float_t*)p->d, p->dim0,
                                  env->stream_togpu());

  } else { //--------------------

      COMET_INSIST_INTERFACE(env, false && "Unimplemented modified gemm method.");

  } // if //--------------------

#else

  COMET_INSIST(false && "Magma library unavailable for this build.");

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------

void gm_linalg_set_matrix_wait(CEnv* env) {
  COMET_INSIST(env);

  env->stream_synchronize(env->stream_togpu());
}

//=============================================================================
/*---Start/end transfer of generic matrix from GPU---*/

void gm_linalg_get_matrix_start(MirroredBuf* p, CEnv* env) {
  COMET_INSIST(p && env);

  if (!env->is_compute_method_gpu()) {
    return;
  }

  /*---Get vectors from GPU---*/

  // ISSUE: these MAGMA routines don't return an error code.

#ifdef COMET_USE_MAGMA

  if (use_minproduct(env)) { //--------------------

    if (FP_PRECISION_DOUBLE) {
      magma_minproduct_dgetmatrix_async(
        p->dim0, p->dim1, (double*)p->d, p->dim0,
        (double*)p->h, p->dim0, env->stream_fromgpu());
    } else {
      magma_minproduct_sgetmatrix_async(
        p->dim0, p->dim1, (float*)p->d, p->dim0,
        (float*)p->h, p->dim0, env->stream_fromgpu());
    }

  } else if (use_mgemm4(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    magma_mgemm4_zgetmatrix_async(p->dim0, p->dim1, (Float_t*)p->d,
                                  p->dim0, (Float_t*)p->h, p->dim0,
                                  env->stream_fromgpu());

  } else if (use_mgemm2(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    magma_mgemm2_zgetmatrix_async(p->dim0, p->dim1, (Float_t*)p->d,
                                  p->dim0, (Float_t*)p->h, p->dim0,
                                  env->stream_fromgpu());

  } else if (use_mgemm3(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    magma_mgemm3_zgetmatrix_async(p->dim0, p->dim1, (Float_t*)p->d,
                                  p->dim0, (Float_t*)p->h, p->dim0,
                                  env->stream_fromgpu());

  } else if (use_mgemm5(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    magma_mgemm5_zgetmatrix_async(p->dim0, p->dim1, (Float_t*)p->d,
                                  p->dim0, (Float_t*)p->h, p->dim0,
                                  env->stream_fromgpu());

  } else { //--------------------

      COMET_INSIST_INTERFACE(env, false && "Unimplemented modified gemm method.");

  } // if //--------------------

#else

  COMET_INSIST(false && "Magma library unavailable for this build.");

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------

void gm_linalg_get_matrix_wait(CEnv* env) {
  COMET_INSIST(env);

  env->stream_synchronize(env->stream_fromgpu());
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
