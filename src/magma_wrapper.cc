//-----------------------------------------------------------------------------
/*!
 * \file   magma_wrapper.cc
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

#if defined COMET_USE_HIP
# define COMET_USE_MAGMA_V2
#endif

#ifdef COMET_USE_MAGMA
#ifdef COMET_USE_MAGMA_V2
#include "magma_minproduct_v2.h"
#include "magma_minproduct_auxiliary.h"
#include "magma_mgemm2_v2.h"
#include "magma_mgemm2_auxiliary.h"
#include "magma_mgemm3_v2.h"
#include "magma_mgemm3_auxiliary.h"
#include "magma_mgemm4_v2.h"
#include "magma_mgemm4_auxiliary.h"
#include "magma_mgemm5_v2.h"
#include "magma_mgemm5_auxiliary.h"
#else
#include "magma_minproduct.h"
#include "magma_minproductblas.h"
#include "magma_mgemm2.h"
#include "magma_mgemm2blas.h"
#include "magma_mgemm3.h"
#include "magma_mgemm3blas.h"
#include "magma_mgemm4.h"
#include "magma_mgemm4blas.h"
#include "magma_mgemm5.h"
#include "magma_mgemm5blas.h"
#endif
//#elif defined COMET_USE_CUDA
//  #include "cublas_v2.h"
#elif defined COMET_USE_HIP
  //#include "hip/hip_runtime_api.h"
//  #include "hip/hip_runtime.h"
#endif

#if defined COMET_USE_HIP
# define COMET_USE_MAGMA_V2
#endif

#if defined COMET_USE_MAGMA_V2
# define PREPEND_COMMA(s) , s
#else
# define PREPEND_COMMA(s)
#endif

#if defined COMET_USE_SEMIRING
#include "semiring.h"
#endif

#include "env.hh"
#include "assertions.hh"
#include "mirrored_buf.hh"
#include "magma_wrapper.hh"

//-----------------------------------------------------------------------------

#define MAGMA_SAFE_CALL(s, msg) \
{ \
  typedef magma_minproduct_int_t magma_generic_int_t; \
  const magma_generic_int_t MAGMA_generic_SUCCESS = MAGMA_minproduct_SUCCESS; \
  const magma_generic_int_t error_code = (s); \
  COMET_INSIST(MAGMA_generic_SUCCESS == error_code && msg); \
  COMET_INSIST(System::accel_last_call_succeeded()); \
}

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/*!
 * \brief Constructor for MagmaWrapper class.
 *
 */
MagmaWrapper::MagmaWrapper(CEnv& env)
  : env_(env) {
  initialize_();
}

//-----------------------------------------------------------------------------
/*!
 * \brief Destructor for MagmaWrapper class.
 *
 */
MagmaWrapper::~MagmaWrapper() {
  finalize_();
}

//-----------------------------------------------------------------------------
/*!
 * \brief Magma setup.
 *
 */
void MagmaWrapper::initialize_() {

  // TODO: non-GPU calls should not need to init magma, should use
  // regular malloc instead of magma malloc.

  if (!env_.is_compute_method_gpu())
    return;

  // TODO: skip this if not actually using magma.

#ifdef COMET_USE_MAGMA

  // may need magma blasSetKernelStream -- see
  // http://on-demand.gputechconf.com/gtc/2014/presentations/S4158-cuda-streams-best-practices-common-pitfalls.pdf
  // page 14

  if (use_minproduct_(env_)) { //--------------------

    MAGMA_SAFE_CALL(magma_minproduct_init(), "Init error.");
#   if ! defined COMET_USE_MAGMA_V2
      MAGMA_SAFE_CALL(magma_minproductblasSetKernelStream(
        env_.stream_compute()), "SetKernelStream.");
#   endif
    env_.queues_initialize<MagmaQueue<MagmaCloneId::MINPRODUCT>>();

  } else if (use_mgemm4_(env_)) { //--------------------

    MAGMA_SAFE_CALL(magma_mgemm4_init(), "Init error.");
#   if ! defined COMET_USE_MAGMA_V2
      MAGMA_SAFE_CALL(magma_mgemm4blasSetKernelStream(env_.stream_compute()),
                      "SetKernelStream.");
#   endif
    env_.queues_initialize<MagmaQueue<MagmaCloneId::MGEMM4>>();

  } else if (use_mgemm2_(env_)) { //--------------------

    MAGMA_SAFE_CALL(magma_mgemm2_init(), "Init error.");
#   if ! defined COMET_USE_MAGMA_V2
      MAGMA_SAFE_CALL(magma_mgemm2blasSetKernelStream(env_.stream_compute()),
                      "SetKernelStream.");
#   endif
    env_.queues_initialize<MagmaQueue<MagmaCloneId::MGEMM2>>();

  } else if (use_mgemm3_(env_)) { //--------------------

    MAGMA_SAFE_CALL(magma_mgemm3_init(), "Init error.");
#   if ! defined COMET_USE_MAGMA_V2
      MAGMA_SAFE_CALL(magma_mgemm3blasSetKernelStream(env_.stream_compute()),
                      "SetKernelStream.");
#   endif
    env_.queues_initialize<MagmaQueue<MagmaCloneId::MGEMM3>>();

  } else if (use_mgemm5_(env_)) { //--------------------

    MAGMA_SAFE_CALL(magma_mgemm5_init(), "Init error.");
#   if ! defined COMET_USE_MAGMA_V2
      MAGMA_SAFE_CALL(magma_mgemm5blasSetKernelStream(env_.stream_compute()),
                      "SetKernelStream.");
#   endif
    env_.queues_initialize<MagmaQueue<MagmaCloneId::MGEMM5>>();

  } // if (use_minproduct_(env_)) //--------------------

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------
/*!
 * \brief Magma teardown.
 *
 */
void MagmaWrapper::finalize_() {

  if (!env_.is_compute_method_gpu())
    return;

  // TODO: skip this if not actually using magma.

#ifdef COMET_USE_MAGMA

  // TODO: (maybe) reset kernel stream (probably not needed)

  if (use_minproduct_(env_)) { //--------------------

    MAGMA_SAFE_CALL(magma_minproduct_finalize(), "Finalize error.");
    env_.queues_terminate<MagmaQueue<MagmaCloneId::MINPRODUCT>>();

  } else if (use_mgemm4_(env_)) { //--------------------

    MAGMA_SAFE_CALL(magma_mgemm4_finalize(), "Finalize error.");
    env_.queues_terminate<MagmaQueue<MagmaCloneId::MGEMM4>>();

  } else if (use_mgemm2_(env_)) { //--------------------

    MAGMA_SAFE_CALL(magma_mgemm2_finalize(), "Finalize error.");
    env_.queues_terminate<MagmaQueue<MagmaCloneId::MGEMM3>>();

  } else if (use_mgemm3_(env_)) { //--------------------

    MAGMA_SAFE_CALL(magma_mgemm3_finalize(), "Finalize error.");
    env_.queues_terminate<MagmaQueue<MagmaCloneId::MGEMM3>>();

  } else if (use_mgemm5_(env_)) { //--------------------

    MAGMA_SAFE_CALL(magma_mgemm5_finalize(), "Finalize error.");
    env_.queues_terminate<MagmaQueue<MagmaCloneId::MGEMM5>>();

  } // if (use_minproduct_(env_)) //--------------------

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------
/*!
 * \brief Allocate host and device memory.
 *
 */
void MagmaWrapper::malloc(MirroredBuf* buf, size_t dim0, size_t dim1,
   CEnv& env) {
  COMET_INSIST(buf);
  COMET_INSIST(dim0 + 1 >= 1 && dim1 + 1 >= 1);
  COMET_INSIST(BuildHas::MAGMA);

#ifdef COMET_USE_MAGMA

  const size_t n = dim0 * dim1;

  if (use_minproduct_(env)) { //--------------------

    if (System::compute_capability() >= 800 || !buf->is_compute_method_gpu()) {
      buf->h = utils::malloc(n * env.sizeof_float(), env); // WORKAROUND
    } else {
      if (env.is_double_prec()) {
        MAGMA_SAFE_CALL(magma_minproduct_dmalloc_pinned(
          reinterpret_cast<double**>(&buf->h), n),
          "Memory allocation error, possibly due to insufficient CPU memory.");
        utils::fill_nan<double>(reinterpret_cast<double*>(buf->h), n);
      } else {
        MAGMA_SAFE_CALL(magma_minproduct_smalloc_pinned(
          reinterpret_cast<float**>(&buf->h), n),
          "Memory allocation error, possibly due to insufficient CPU memory.");
        utils::fill_nan<float>(reinterpret_cast<float*>(buf->h), n);
     }
    } // if (System::compute_capability() >= 800 || !buf->is_compute_method_gpu())

    if (env.is_compute_method_gpu()) {
      if (env.is_double_prec()) {
        MAGMA_SAFE_CALL(magma_minproduct_dmalloc(
          reinterpret_cast<double**>(&buf->d), n),
          "Memory allocation error, possibly due to insufficient GPU memory.");
        // TODO: ? fill GPU memory with NaNs
      } else {
        MAGMA_SAFE_CALL(magma_minproduct_smalloc(
          reinterpret_cast<float**>(&buf->d), n),
          "Memory allocation error, possibly due to insufficient GPU memory.");
        // TODO: ? fill GPU memory with NaNs
     }
    } // if (env.is_compute_method_gpu())

  } else if (use_mgemm4_(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    if (System::compute_capability() >= 800 || !buf->is_compute_method_gpu())
      buf->h = utils::malloc(n * sizeof(Float_t), env); // WORKAROUND
    else
      MAGMA_SAFE_CALL(magma_mgemm4_zmalloc_pinned(
        reinterpret_cast<Float_t**>(&buf->h), n),
        "Memory allocation error, possibly due to insufficient CPU memory.");

    if (env.is_compute_method_gpu())
      MAGMA_SAFE_CALL(magma_mgemm4_zmalloc(
        reinterpret_cast<Float_t**>(&buf->d), n),
        "Memory allocation error, possibly due to insufficient CPU memory.");

  } else if (use_mgemm2_(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    if (System::compute_capability() >= 800 || !buf->is_compute_method_gpu())
      buf->h = utils::malloc(n * sizeof(Float_t), env); // WORKAROUND
    else
      MAGMA_SAFE_CALL(magma_mgemm2_zmalloc_pinned(
        reinterpret_cast<Float_t**>(&buf->h), n),
        "Memory allocation error, possibly due to insufficient CPU memory.");

    if (env.is_compute_method_gpu())
      MAGMA_SAFE_CALL(magma_mgemm2_zmalloc(
        reinterpret_cast<Float_t**>(&buf->d), n),
        "Memory allocation error, possibly due to insufficient CPU memory.");

  } else if (use_mgemm3_(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    if (System::compute_capability() >= 800 || !buf->is_compute_method_gpu())
      buf->h = utils::malloc(n * sizeof(Float_t), env); // WORKAROUND
    else
      MAGMA_SAFE_CALL(magma_mgemm3_zmalloc_pinned(
        reinterpret_cast<Float_t**>(&buf->h), n),
        "Memory allocation error, possibly due to insufficient CPU memory.");

    if (env.is_compute_method_gpu())
      MAGMA_SAFE_CALL(magma_mgemm3_zmalloc(
        reinterpret_cast<Float_t**>(&buf->d), n),
        "Memory allocation error, possibly due to insufficient CPU memory.");

  } else if (use_mgemm5_(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    if (System::compute_capability() >= 800 || !buf->is_compute_method_gpu())
      buf->h = utils::malloc(n * sizeof(Float_t), env); // WORKAROUND
    else
      MAGMA_SAFE_CALL(magma_mgemm5_zmalloc_pinned(
        reinterpret_cast<Float_t**>(&buf->h), n),
        "Memory allocation error, possibly due to insufficient CPU memory.");

    if (env.is_compute_method_gpu())
      MAGMA_SAFE_CALL(magma_mgemm5_zmalloc(
        reinterpret_cast<Float_t**>(&buf->d), n),
        "Memory allocation error, possibly due to insufficient CPU memory.");

  } // if (use_minproduct_(env_)) //--------------------

  COMET_INSIST(buf->h &&
    "Invalid host pointer created, possibly due to insufficient CPU memory.");
  COMET_INSIST((buf->d || !env.is_compute_method_gpu()) &&
    "Invalid device pointer created, possibly due to insufficient GPU memory.");

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------
/*!
 * \brief Free host and device memory.
 *
 */
void MagmaWrapper::free(MirroredBuf* buf, size_t dim0, size_t dim1,
  CEnv& env) {
  COMET_INSIST(buf);
  COMET_INSIST(! buf->is_alias);
  COMET_INSIST(BuildHas::MAGMA);

#ifdef COMET_USE_MAGMA

  const size_t n = dim0 * dim1;

  if (use_minproduct_(env)) { //--------------------

    if (System::compute_capability() >= 800 || !buf->is_compute_method_gpu()) {
      utils::free(buf->h, n * env.sizeof_float(), env); // WORKAROUND
    } else {
      MAGMA_SAFE_CALL(magma_minproduct_free_pinned(buf->h),
                      "Error in CPU memory free.");
    }
    if (env.is_compute_method_gpu())
      MAGMA_SAFE_CALL(magma_minproduct_free(buf->d),
                      "Error in GPU memory free.");

  } else if (use_mgemm4_(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    if (System::compute_capability() >= 800 || !buf->is_compute_method_gpu()) {
      utils::free(buf->h, n * sizeof(Float_t), env); // WORKAROUND
    } else {
      MAGMA_SAFE_CALL(magma_mgemm4_free_pinned(buf->h),
                      "Error in CPU memory free.");
    }
    if (env.is_compute_method_gpu())
      MAGMA_SAFE_CALL(magma_mgemm4_free(buf->d),
                      "Error in GPU memory free.");

  } else if (use_mgemm2_(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    if (System::compute_capability() >= 800 || !buf->is_compute_method_gpu()) {
      utils::free(buf->h, n * sizeof(Float_t), env); // WORKAROUND
    } else {
      MAGMA_SAFE_CALL(magma_mgemm2_free_pinned(buf->h),
                      "Error in CPU memory free.");
    }
    if (env.is_compute_method_gpu())
      MAGMA_SAFE_CALL(magma_mgemm2_free(buf->d),
                      "Error in GPU memory free.");

  } else if (use_mgemm3_(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    if (System::compute_capability() >= 800 || !buf->is_compute_method_gpu()) {
      utils::free(buf->h, n * sizeof(Float_t), env); // WORKAROUND
    } else {
      MAGMA_SAFE_CALL(magma_mgemm3_free_pinned(buf->h),
                      "Error in CPU memory free.");
    }
    if (env.is_compute_method_gpu())
      MAGMA_SAFE_CALL(magma_mgemm3_free(buf->d),
                      "Error in GPU memory free.");

  } else { // if (use_mgemm5_(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    if (System::compute_capability() >= 800 || !buf->is_compute_method_gpu()) {
      utils::free(buf->h, n * sizeof(Float_t), env); // WORKAROUND
    } else {
      MAGMA_SAFE_CALL(magma_mgemm5_free_pinned(buf->h),
                      "Error in CPU memory free.");
    }
    if (env.is_compute_method_gpu())
      MAGMA_SAFE_CALL(magma_mgemm5_free(buf->d),
                      "Error in GPU memory free.");

  } // if (use_minproduct_(env_)) //--------------------

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------
/*!
 * \brief Set matrix to zero.
 *
 */
void MagmaWrapper::set_matrix_zero_start(MirroredBuf* buf, CEnv& env) {
  COMET_INSIST(buf);
  COMET_INSIST(BuildHas::MAGMA);

  if (!env.is_compute_method_gpu())
    return;

#ifdef COMET_USE_MAGMA

  const size_t mat_dim1 = buf->dim0;
  const size_t mat_dim2 = buf->dim1;

  // NOTE: these MAGMA routines don't return an error code.

  COMET_INSIST(env.queue_compute_.magma_queue);

  if (use_minproduct_(env) && env.is_double_prec()) { //--------------------

    typedef double Float_t;
    const Float_t zero = 0;

    magma_minproductblas_dlaset(Magma_minproductFull, mat_dim1, mat_dim2,
      zero, zero, reinterpret_cast<Float_t*>(buf->d), mat_dim1
      PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MINPRODUCT>>()));

  } else if (use_minproduct_(env)) { //--------------------

    typedef float Float_t;
    const Float_t zero = 0;

    magma_minproductblas_slaset(Magma_minproductFull, mat_dim1, mat_dim2,
      zero, zero, reinterpret_cast<Float_t*>(buf->d), mat_dim1
      PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MINPRODUCT>>()));

  } else if (use_mgemm4_(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;
    const Float_t zero = {0, 0};

    magma_mgemm4blas_zlaset(Magma_mgemm4Full, mat_dim1, mat_dim2,
      zero, zero, reinterpret_cast<Float_t*>(buf->d), mat_dim1
      PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MGEMM4>>()));

  } else if (use_mgemm2_(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;
    Float_t zero = {0, 0};

    magma_mgemm2blas_zlaset(Magma_mgemm2Full, mat_dim1, mat_dim2,
      zero, zero, reinterpret_cast<Float_t*>(buf->d), mat_dim1
      PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MGEMM2>>()));

  } else if (use_mgemm3_(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;
    Float_t zero = {0, 0};

    magma_mgemm3blas_zlaset(Magma_mgemm3Full, mat_dim1, mat_dim2,
      zero, zero, reinterpret_cast<Float_t*>(buf->d), mat_dim1
      PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MGEMM3>>()));

  } else if (use_mgemm5_(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;
    Float_t zero = {0, 0};

    magma_mgemm5blas_zlaset(Magma_mgemm5Full, mat_dim1, mat_dim2,
      zero, zero, reinterpret_cast<Float_t*>(buf->d), mat_dim1
      PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MGEMM5>>()));

  } // if (use_minproduct_(env_)) //--------------------

  COMET_INSIST(System::accel_last_call_succeeded());

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------
/*!
 * \brief GEMM on block of matrix, start.
 *
 */
void MagmaWrapper::gemm_block_start(size_t m, size_t n, size_t k,
  const void* matA, size_t ldda, const void* matB, size_t lddb,
  void* matC, size_t lddc, bool is_beta_one, CEnv& env) {
  COMET_INSIST(matA && matB && matC);
  COMET_INSIST(env.is_compute_method_gpu());
  COMET_INSIST(BuildHas::MAGMA);

#ifdef COMET_USE_MAGMA

  // Ensure Magmablas function doesn't internally failover to CUBLAS call.
  int TransA = 1;
  int TransB = 0;

  const size_t Am = ! TransA ? m : k;
  const size_t An = ! TransA ? k : m;
  const size_t Bm = ! TransB ? k : n;
  const size_t Bn = ! TransB ? n : k;
  const size_t sizeA = ldda * (An - 1) + Am;
  const size_t sizeB = lddb * (Bn - 1) + Bm;

  const size_t CUBLAS_MAX_1DBUF_SIZE = ((1 << 27) - 512);
  COMET_INSIST((! (sizeA >= CUBLAS_MAX_1DBUF_SIZE ||
               sizeB >= CUBLAS_MAX_1DBUF_SIZE )) &&
           "Error in MAGMA block sizes.");

  typedef magma_minproduct_int_t magma_generic_int_t;
  typedef magma_generic_int_t Int_t;

  const auto m_ = safe_cast_insist<Int_t>(m);
  const auto n_ = safe_cast_insist<Int_t>(n);
  const auto k_ = safe_cast_insist<Int_t>(k);
  const auto ldda_ = safe_cast_insist<Int_t>(ldda);
  const auto lddb_ = safe_cast_insist<Int_t>(lddb);
  const auto lddc_ = safe_cast_insist<Int_t>(lddc);

  // ISSUE: these MAGMA routines don't return an error code.

  if (use_minproduct_(env) && env.is_double_prec()) { //--------------------

    typedef double Float_t;

    const Float_t alpha = 1;
    const Float_t beta = is_beta_one ? 1 : 0;

    magma_minproductblas_dgemm(Magma_minproductTrans, Magma_minproductNoTrans,
      m_, n_, k_,
      alpha, reinterpret_cast<const Float_t*>(matA), ldda_,
      reinterpret_cast<const Float_t*>(matB), lddb_, beta,
      reinterpret_cast<Float_t*>(matC), lddc_
      PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MINPRODUCT>>()));
    COMET_INSIST(System::accel_last_call_succeeded() &&
                 "Failure in call to magma_minproductblas_dgemm.");

    env.ops_local_inc(static_cast<double>(2) * m * static_cast<double>(n) * k);

  } else if (use_minproduct_(env)) { // && !env.is_double_prec()) { //--------------------

    typedef float Float_t;

    const Float_t alpha = 1;
    const Float_t beta = is_beta_one ? 1 : 0;

    magma_minproductblas_sgemm(Magma_minproductTrans, Magma_minproductNoTrans,
      m_, n_, k_,
      alpha, reinterpret_cast<const Float_t*>(matA), ldda_,
      reinterpret_cast<const Float_t*>(matB), lddb_, beta,
      reinterpret_cast<Float_t*>(matC), lddc_
      PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MINPRODUCT>>()));
    COMET_INSIST(System::accel_last_call_succeeded() &&
                 "Failure in call to magma_minproductblas_sgemm.");

    env.ops_local_inc(static_cast<double>(2) * m * static_cast<double>(n) * k);

  } else if (use_mgemm4_(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const auto alpha = one;
    const auto beta = is_beta_one ? one : zero;

    magma_mgemm4blas_zgemm(Magma_mgemm4Trans, Magma_mgemm4NoTrans,
      m_, n_, k_,
      alpha, reinterpret_cast<const Float_t*>(matA), ldda_,
      reinterpret_cast<const Float_t*>(matB), lddb_, beta,
      reinterpret_cast<Float_t*>(matC), lddc_
      PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MGEMM4>>()));
    COMET_INSIST(System::accel_last_call_succeeded() &&
                 "Failure in call to magma_mgemm4blas_zgemm.");

  } else if (use_mgemm2_(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const auto alpha = one;
    const auto beta = is_beta_one ? one : zero;

    magma_mgemm2blas_zgemm(Magma_mgemm2Trans, Magma_mgemm2NoTrans,
      m_, n_, k_,
      alpha, reinterpret_cast<const Float_t*>(matA), ldda_,
      reinterpret_cast<const Float_t*>(matB), lddb_, beta,
      reinterpret_cast<Float_t*>(matC), lddc_
      PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MGEMM2>>()));
    COMET_INSIST(System::accel_last_call_succeeded() &&
                 "Failure in call to magma_mgemm2blas_zgemm.");

  } else if (use_mgemm3_(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const auto alpha = one;
    const auto beta = is_beta_one ? one : zero;

    magma_mgemm3blas_zgemm(Magma_mgemm3Trans, Magma_mgemm3NoTrans,
      m_, n_, k_,
      alpha, reinterpret_cast<const Float_t*>(matA), ldda_,
      reinterpret_cast<const Float_t*>(matB), lddb_, beta,
      reinterpret_cast<Float_t*>(matC), lddc_
      PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MGEMM3>>()));
    COMET_INSIST(System::accel_last_call_succeeded() &&
                 "Failure in call to magma_mgemm3blas_zgemm.");

  } else if (use_mgemm5_(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const auto alpha = one;
    const auto beta = is_beta_one ? one : zero;

    magma_mgemm5blas_zgemm(Magma_mgemm5Trans, Magma_mgemm5NoTrans,
      m_, n_, k_,
      alpha, reinterpret_cast<const Float_t*>(matA), ldda_,
      reinterpret_cast<const Float_t*>(matB), lddb_, beta,
      reinterpret_cast<Float_t*>(matC), lddc_
      PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MGEMM5>>()));
    COMET_INSIST(System::accel_last_call_succeeded() &&
                 "Failure in call to magma_mgemm5blas_zgemm.");

  } // if (use_minproduct_(env_)) //--------------------

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------
/*!
 * \brief GEMM on matrix, start.
 *
 */
void MagmaWrapper::gemm_start(size_t m, size_t n, size_t k,
  const void* matA, size_t ldda, const void* matB, size_t lddb,
  void* matC, size_t lddc, CEnv& env) {
  COMET_INSIST(matA && matB && matC);
  COMET_INSIST(env.is_compute_method_gpu());
  COMET_INSIST(BuildHas::MAGMA);

#if defined COMET_USE_SEMIRING

  COMET_STATIC_ASSERT(BuildHas::MAGMA &&
                      "seMIring build requires MAGMA enabled for support.");

  if (use_minproduct_(env)) { //--------------------

    const auto m_ = safe_cast_insist<int>(m);
    const auto n_ = safe_cast_insist<int>(n);
    const auto k_ = safe_cast_insist<int>(k);
    const auto ldda_ = safe_cast_insist<std::size_t>(ldda);
    const auto lddb_ = safe_cast_insist<std::size_t>(lddb);
    const auto lddc_ = safe_cast_insist<std::size_t>(lddc);

    if (env.is_double_prec())
      srPlusMinFP64TN(m_, n_, k_,
                      static_cast<const double*>(matA), ldda_,
                      static_cast<const double*>(matB), lddb_,
                      static_cast<double*>(matC), lddc_,
                      env.stream_compute());
    else
      srPlusMinFP32TN(m_, n_, k_,
                      static_cast<const float*>(matA), ldda_,
                      static_cast<const float*>(matB), lddb_,
                      static_cast<float*>(matC), lddc_,
                      env.stream_compute());

    COMET_INSIST(System::accel_last_call_succeeded() &&
             "Failure in call to seMIring.");

    env.ops_local_inc(static_cast<double>(2) * m * static_cast<double>(n) * k);

    return; // exit here, no need to call MAGMA.
  } // if (use_minproduct_(env))

#endif // COMET_USE_SEMIRING

#if defined COMET_USE_MAGMA

  // The purpose of this code is to workaround the magma size
  // limitation (for non CUBLAS failover) by doing gemm in blocks.

  const size_t rows = k;
  const size_t cols_A = m;
  const size_t cols_B = n;

  const size_t elt_size =
    use_minproduct_(env) ? env.sizeof_float() :
    use_mgemm4_(env) ? sizeof(magma_mgemm4DoubleComplex) :
    use_mgemm2_(env) ? sizeof(magma_mgemm2DoubleComplex) :
    use_mgemm3_(env) ? sizeof(magma_mgemm3DoubleComplex) :
    use_mgemm5_(env) ? sizeof(magma_mgemm5DoubleComplex) : 0;
  COMET_INSIST(elt_size > 0 && "Error in gemm block calculation.");

  const size_t align_factor = 128 / elt_size;
  const size_t max_elts = (1 << 27) - 512;

  // TODO: can we improve aspect ratios of submatrices.

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

      void* A_this = static_cast<char*>(const_cast<void*>(matA)) +
        (row_base + ldda*col_A_base)*elt_size;

      for (size_t col_B_base=0; col_B_base<cols_B;
           col_B_base+=cols_per_block_B) {

        const size_t cols_B_remaining = cols_B - col_B_base;
        const size_t cols_B_this = utils::min(cols_B_remaining,
                                             cols_per_block_B);

        void* B_this = static_cast<char*>(const_cast<void*>(matB)) +
          (row_base + ldda*col_B_base)*elt_size;

        void* C_this = static_cast<char*>(matC) +
          (col_A_base + lddc*col_B_base)*elt_size;

        MagmaWrapper::gemm_block_start(cols_A_this, cols_B_this, rows_this,
          A_this, ldda, B_this, lddb, C_this, lddc, row_base > 0,  env);
      } // for col_B_base
    } // for col_A_base
  } // for row_base

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------
/*!
 * \brief Start transfer of matrix to GPU.
 *
 */
void MagmaWrapper::set_matrix_start(MirroredBuf* buf, CEnv& env) {
  COMET_INSIST(buf);
  COMET_INSIST(BuildHas::MAGMA);

  if (!env.is_compute_method_gpu())
    return;

#ifdef COMET_USE_MAGMA

  // NOTE: these MAGMA routines don't return an error code.

  if (use_minproduct_(env) && env.is_double_prec()) { //--------------------

    typedef double Float_t;

    magma_minproduct_dsetmatrix_async(buf->dim0, buf->dim1,
      reinterpret_cast<Float_t*>(buf->h), buf->dim0,
      reinterpret_cast<Float_t*>(buf->d), buf->dim0,
      env.queue_togpu<MagmaQueue<MagmaCloneId::MINPRODUCT>>());

  } else if (use_minproduct_(env)) { //--------------------

    typedef float Float_t;

    magma_minproduct_ssetmatrix_async(buf->dim0, buf->dim1,
      reinterpret_cast<Float_t*>(buf->h), buf->dim0,
      reinterpret_cast<Float_t*>(buf->d), buf->dim0,
      env.queue_togpu<MagmaQueue<MagmaCloneId::MINPRODUCT>>());

  } else if (use_mgemm4_(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    magma_mgemm4_zsetmatrix_async(buf->dim0, buf->dim1,
      reinterpret_cast<Float_t*>(buf->h), buf->dim0,
      reinterpret_cast<Float_t*>(buf->d), buf->dim0,
      env.queue_togpu<MagmaQueue<MagmaCloneId::MGEMM4>>());

  } else if (use_mgemm2_(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    magma_mgemm2_zsetmatrix_async(buf->dim0, buf->dim1,
      reinterpret_cast<Float_t*>(buf->h), buf->dim0,
      reinterpret_cast<Float_t*>(buf->d), buf->dim0,
      env.queue_togpu<MagmaQueue<MagmaCloneId::MGEMM2>>());

  } else if (use_mgemm3_(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    magma_mgemm3_zsetmatrix_async(buf->dim0, buf->dim1,
      reinterpret_cast<Float_t*>(buf->h), buf->dim0,
      reinterpret_cast<Float_t*>(buf->d), buf->dim0,
      env.queue_togpu<MagmaQueue<MagmaCloneId::MGEMM3>>());

  } else if (use_mgemm5_(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    magma_mgemm5_zsetmatrix_async(buf->dim0, buf->dim1,
      reinterpret_cast<Float_t*>(buf->h), buf->dim0,
      reinterpret_cast<Float_t*>(buf->d), buf->dim0,
      env.queue_togpu<MagmaQueue<MagmaCloneId::MGEMM5>>());

  } // if (use_minproduct_(env_) && env.is_double_prec()) //--------------------

  COMET_INSIST(System::accel_last_call_succeeded());

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------
/*!
 * \brief Wait transfer of matrix to GPU.
 *
 */
void MagmaWrapper::set_matrix_wait(CEnv& env) {
  COMET_INSIST(BuildHas::MAGMA);
  env.stream_synchronize(env.stream_togpu());
}

//-----------------------------------------------------------------------------
/*!
 * \brief Start transfer of matrix from GPU.
 *
 */
void MagmaWrapper::get_matrix_start(MirroredBuf* buf, CEnv& env) {
  COMET_INSIST(buf);
  COMET_INSIST(BuildHas::MAGMA);

  if (!env.is_compute_method_gpu())
    return;

#ifdef COMET_USE_MAGMA

  // NOTE: these MAGMA routines don't return an error code.

  if (use_minproduct_(env) && env.is_double_prec()) { //--------------------

    typedef double Float_t;

    magma_minproduct_dgetmatrix_async(buf->dim0, buf->dim1,
      reinterpret_cast<Float_t*>(buf->d), buf->dim0,
      reinterpret_cast<Float_t*>(buf->h), buf->dim0,
      env.queue_fromgpu<MagmaQueue<MagmaCloneId::MINPRODUCT>>());

  } else if (use_minproduct_(env)) { //--------------------

    typedef float Float_t;

    magma_minproduct_sgetmatrix_async(buf->dim0, buf->dim1,
      reinterpret_cast<Float_t*>(buf->d), buf->dim0,
      reinterpret_cast<Float_t*>(buf->h), buf->dim0,
      env.queue_fromgpu<MagmaQueue<MagmaCloneId::MINPRODUCT>>());

  } else if (use_mgemm4_(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    magma_mgemm4_zgetmatrix_async(buf->dim0, buf->dim1,
      reinterpret_cast<Float_t*>(buf->d), buf->dim0,
      reinterpret_cast<Float_t*>(buf->h), buf->dim0,
      env.queue_fromgpu<MagmaQueue<MagmaCloneId::MGEMM4>>());

  } else if (use_mgemm2_(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    magma_mgemm2_zgetmatrix_async(buf->dim0, buf->dim1,
      reinterpret_cast<Float_t*>(buf->d), buf->dim0,
      reinterpret_cast<Float_t*>(buf->h), buf->dim0,
      env.queue_fromgpu<MagmaQueue<MagmaCloneId::MGEMM2>>());

  } else if (use_mgemm3_(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    magma_mgemm3_zgetmatrix_async(buf->dim0, buf->dim1,
      reinterpret_cast<Float_t*>(buf->d), buf->dim0,
      reinterpret_cast<Float_t*>(buf->h), buf->dim0,
      env.queue_fromgpu<MagmaQueue<MagmaCloneId::MGEMM3>>());

  } else if (use_mgemm5_(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    magma_mgemm5_zgetmatrix_async(buf->dim0, buf->dim1,
      reinterpret_cast<Float_t*>(buf->d), buf->dim0,
      reinterpret_cast<Float_t*>(buf->h), buf->dim0,
      env.queue_fromgpu<MagmaQueue<MagmaCloneId::MGEMM5>>());

  } // if (use_minproduct_(env_) && env.is_double_prec()) //--------------------

  COMET_INSIST(System::accel_last_call_succeeded());

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------
/*!
 * \brief Wait transfer of matrix from GPU.
 *
 */
void MagmaWrapper::get_matrix_wait(CEnv& env) {
  COMET_INSIST(BuildHas::MAGMA);
  env.stream_synchronize(env.stream_fromgpu());
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
