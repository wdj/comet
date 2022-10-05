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
# define PREPEND_COMMA(q) , q
#else
# define PREPEND_COMMA(q)
#endif

#if defined COMET_USE_SEMIRING
#include "semiring.h"
#endif

#include "env.hh"
#include "assertions.hh"
#include "mirrored_buf.hh"
#include "magma_wrapper.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Constructor for MagmaWrapper class.

MagmaWrapper::MagmaWrapper(CEnv& env)
  : env_(env) {
  initialize_(env_);
}

//-----------------------------------------------------------------------------
/// \brief Destructor for MagmaWrapper class.

MagmaWrapper::~MagmaWrapper() {
  finalize_(env_);
}

//-----------------------------------------------------------------------------
/// \brief Magma setup.

void MagmaWrapper::initialize_(CEnv& env) {

  if (!env.is_compute_method_gpu())
    return;

  // need magma blasSetKernelStream -- see
  // http://on-demand.gputechconf.com/gtc/2014/presentations/S4158-cuda-streams-best-practices-common-pitfalls.pdf
  // page 14

// TODO: non-GPU calls should not need to init magma, should use
// regular malloc instead of magma malloc.

#ifdef COMET_USE_MAGMA

  if (use_minproduct_(env)) { //--------------------

    magma_minproduct_int_t magma_code = magma_minproduct_init();
    COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS && "Init error.");
#   if ! defined COMET_USE_MAGMA_V2
      magma_code = magma_minproductblasSetKernelStream(env.stream_compute());
      COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS && "SetKernelStream.");
#   endif

//FIX
    env.queues_initialize<MagmaQueue<MagmaCloneId::MINPRODUCT>>();

    //env.queue_compute_.magma_queue = (void*)new MagmaQueue<MagmaCloneId::MINPRODUCT>(env.stream_compute(), env);
    //env.queue_togpu_.magma_queue = (void*)new MagmaQueue<MagmaCloneId::MINPRODUCT>(env.stream_togpu(), env);
    //env.queue_fromgpu_.magma_queue = (void*)new MagmaQueue<MagmaCloneId::MINPRODUCT>(env.stream_fromgpu(), env);
    //env.queue_compute_.is_initialized = true;
    //env.queue_togpu_.is_initialized = true;
    //env.queue_fromgpu_.is_initialized = true;

  } else if (use_mgemm4_(env)) { //--------------------

    magma_mgemm4_int_t magma_code = magma_mgemm4_init();
    COMET_INSIST(magma_code == MAGMA_mgemm4_SUCCESS && "Init error.");
#   if ! defined COMET_USE_MAGMA_V2
    magma_code = magma_mgemm4blasSetKernelStream(env.stream_compute());
    COMET_INSIST(magma_code == MAGMA_mgemm4_SUCCESS && "SetKernelStream.");
#   endif

//FIX
    env.queues_initialize<MagmaQueue<MagmaCloneId::MGEMM4>>();
    //env.queue_compute_.magma_queue = (void*)new MagmaQueue<MagmaCloneId::MGEMM4>(env.stream_compute(), env);
    //env.queue_togpu_.magma_queue = (void*)new MagmaQueue<MagmaCloneId::MGEMM4>(env.stream_togpu(), env);
    //env.queue_fromgpu_.magma_queue = (void*)new MagmaQueue<MagmaCloneId::MGEMM4>(env.stream_fromgpu(), env);
    //env.queue_compute_.is_initialized = true;
    //env.queue_togpu_.is_initialized = true;
    //env.queue_fromgpu_.is_initialized = true;

  } else if (use_mgemm2_(env)) { //--------------------

    magma_mgemm2_int_t magma_code = magma_mgemm2_init();
    COMET_INSIST(magma_code == MAGMA_mgemm2_SUCCESS && "Init error.");
#   if ! defined COMET_USE_MAGMA_V2
    magma_code = magma_mgemm2blasSetKernelStream(env.stream_compute());
    COMET_INSIST(magma_code == MAGMA_mgemm2_SUCCESS && "SetKernelStream.");
#   endif

//FIX
    env.queues_initialize<MagmaQueue<MagmaCloneId::MGEMM2>>();
    //env.queue_compute_.magma_queue = (void*)new MagmaQueue<MagmaCloneId::MGEMM2>(env.stream_compute(), env);
    //env.queue_togpu_.magma_queue = (void*)new MagmaQueue<MagmaCloneId::MGEMM2>(env.stream_togpu(), env);
    //env.queue_fromgpu_.magma_queue = (void*)new MagmaQueue<MagmaCloneId::MGEMM2>(env.stream_fromgpu(), env);
    //env.queue_compute_.is_initialized = true;
    //env.queue_togpu_.is_initialized = true;
    //env.queue_fromgpu_.is_initialized = true;

  } else if (use_mgemm3_(env)) { //--------------------

    magma_mgemm3_int_t magma_code = magma_mgemm3_init();
    COMET_INSIST(magma_code == MAGMA_mgemm3_SUCCESS && "Init error.");
#   if ! defined COMET_USE_MAGMA_V2
    magma_code = magma_mgemm3blasSetKernelStream(env.stream_compute());
    COMET_INSIST(magma_code == MAGMA_mgemm3_SUCCESS && "SetKernelStream.");
#   endif

//FIX
    env.queues_initialize<MagmaQueue<MagmaCloneId::MGEMM3>>();
    //env.queue_compute_.magma_queue = (void*)new MagmaQueue<MagmaCloneId::MGEMM3>(env.stream_compute(), env);
    //env.queue_togpu_.magma_queue = (void*)new MagmaQueue<MagmaCloneId::MGEMM3>(env.stream_togpu(), env);
    //env.queue_fromgpu_.magma_queue = (void*)new MagmaQueue<MagmaCloneId::MGEMM3>(env.stream_fromgpu(), env);
    //env.queue_compute_.is_initialized = true;
    //env.queue_togpu_.is_initialized = true;
    //env.queue_fromgpu_.is_initialized = true;

  } else if (use_mgemm5_(env)) { //--------------------

    magma_mgemm5_int_t magma_code = magma_mgemm5_init();
    COMET_INSIST(magma_code == MAGMA_mgemm5_SUCCESS && "Init error.");
#   if ! defined COMET_USE_MAGMA_V2
    magma_code = magma_mgemm5blasSetKernelStream(env.stream_compute());
    COMET_INSIST(magma_code == MAGMA_mgemm5_SUCCESS && "SetKernelStream.");
#   endif

//FIX
    env.queues_initialize<MagmaQueue<MagmaCloneId::MGEMM5>>();
    //env.queue_compute_.magma_queue = (void*)new MagmaQueue<MagmaCloneId::MGEMM5>(env.stream_compute(), env);
    //env.queue_togpu_.magma_queue = (void*)new MagmaQueue<MagmaCloneId::MGEMM5>(env.stream_togpu(), env);
    //env.queue_fromgpu_.magma_queue = (void*)new MagmaQueue<MagmaCloneId::MGEMM5>(env.stream_fromgpu(), env);
    //env.queue_compute_.is_initialized = true;
    //env.queue_togpu_.is_initialized = true;
    //env.queue_fromgpu_.is_initialized = true;

  } else { //--------------------

    COMET_INSIST_INTERFACE(&env, false && "Unimplemented method.");

  } // if //--------------------

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------
/// \brief Magma teardown.

void MagmaWrapper::finalize_(CEnv& env) {

  if (!env.is_compute_method_gpu())
    return;

  // TODO: (maybe) reset kernel stream (probably not really needed)

#ifdef COMET_USE_MAGMA

  if (use_minproduct_(env)) { //--------------------

    magma_minproduct_int_t magma_code = magma_minproduct_finalize();
    COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS && "Finalize error.");

//FIX
    env.queues_terminate<MagmaQueue<MagmaCloneId::MINPRODUCT>>();
    //delete (MagmaQueue<MagmaCloneId::MINPRODUCT>*)env.queue_compute_.magma_queue;
    //delete (MagmaQueue<MagmaCloneId::MINPRODUCT>*)env.queue_togpu_.magma_queue;
    //delete (MagmaQueue<MagmaCloneId::MINPRODUCT>*)env.queue_fromgpu_.magma_queue;
    //env.queue_compute_.magma_queue = NULL;
    //env.queue_togpu_.magma_queue = NULL;
    //env.queue_fromgpu_.magma_queue = NULL;
    //env.queue_compute_.is_initialized = false;
    //env.queue_togpu_.is_initialized = false;
    //env.queue_fromgpu_.is_initialized = false;

  } else if (use_mgemm4_(env)) { //--------------------

    magma_mgemm4_int_t magma_code = magma_mgemm4_finalize();
    COMET_INSIST(magma_code == MAGMA_mgemm4_SUCCESS && "Finalize error.");

//FIX
    env.queues_terminate<MagmaQueue<MagmaCloneId::MGEMM4>>();
    //delete (MagmaQueue<MagmaCloneId::MGEMM4>*)env.queue_compute_.magma_queue;
    //delete (MagmaQueue<MagmaCloneId::MGEMM4>*)env.queue_togpu_.magma_queue;
    //delete (MagmaQueue<MagmaCloneId::MGEMM4>*)env.queue_fromgpu_.magma_queue;
    //env.queue_compute_.magma_queue = NULL;
    //env.queue_togpu_.magma_queue = NULL;
    //env.queue_fromgpu_.magma_queue = NULL;
    //env.queue_compute_.is_initialized = false;
    //env.queue_togpu_.is_initialized = false;
    //env.queue_fromgpu_.is_initialized = false;

  } else if (use_mgemm2_(env)) { //--------------------

    magma_mgemm2_int_t magma_code = magma_mgemm2_finalize();
    COMET_INSIST(magma_code == MAGMA_mgemm2_SUCCESS && "Finalize error.");

//FIX
    env.queues_terminate<MagmaQueue<MagmaCloneId::MGEMM2>>();
    //delete (MagmaQueue<MagmaCloneId::MGEMM2>*)env.queue_compute_.magma_queue;
    //delete (MagmaQueue<MagmaCloneId::MGEMM2>*)env.queue_togpu_.magma_queue;
    //delete (MagmaQueue<MagmaCloneId::MGEMM2>*)env.queue_fromgpu_.magma_queue;
    //env.queue_compute_.magma_queue = NULL;
    //env.queue_togpu_.magma_queue = NULL;
    //env.queue_fromgpu_.magma_queue = NULL;
    //env.queue_compute_.is_initialized = false;
    //env.queue_togpu_.is_initialized = false;
    //env.queue_fromgpu_.is_initialized = false;

  } else if (use_mgemm3_(env)) { //--------------------

    magma_mgemm3_int_t magma_code = magma_mgemm3_finalize();
    COMET_INSIST(magma_code == MAGMA_mgemm3_SUCCESS && "Finalize error.");

//FIX
    env.queues_terminate<MagmaQueue<MagmaCloneId::MGEMM3>>();
    //delete (MagmaQueue<MagmaCloneId::MGEMM3>*)env.queue_compute_.magma_queue;
    //delete (MagmaQueue<MagmaCloneId::MGEMM3>*)env.queue_togpu_.magma_queue;
    //delete (MagmaQueue<MagmaCloneId::MGEMM3>*)env.queue_fromgpu_.magma_queue;
    //env.queue_compute_.magma_queue = NULL;
    //env.queue_togpu_.magma_queue = NULL;
    //env.queue_fromgpu_.magma_queue = NULL;
    //env.queue_compute_.is_initialized = false;
    //env.queue_togpu_.is_initialized = false;
    //env.queue_fromgpu_.is_initialized = false;

  } else if (use_mgemm5_(env)) { //--------------------

    magma_mgemm5_int_t magma_code = magma_mgemm5_finalize();
    COMET_INSIST(magma_code == MAGMA_mgemm5_SUCCESS && "Finalize error.");

//FIX
    env.queues_terminate<MagmaQueue<MagmaCloneId::MGEMM5>>();
    //delete (MagmaQueue<MagmaCloneId::MGEMM5>*)env.queue_compute_.magma_queue;
    //delete (MagmaQueue<MagmaCloneId::MGEMM5>*)env.queue_togpu_.magma_queue;
    //delete (MagmaQueue<MagmaCloneId::MGEMM5>*)env.queue_fromgpu_.magma_queue;
    //env.queue_compute_.magma_queue = NULL;
    //env.queue_togpu_.magma_queue = NULL;
    //env.queue_fromgpu_.magma_queue = NULL;
    //env.queue_compute_.is_initialized = false;
    //env.queue_togpu_.is_initialized = false;
    //env.queue_fromgpu_.is_initialized = false;

  } else { //--------------------

    COMET_INSIST_INTERFACE(&env, false && "Unimplemented method.");

  } // if //--------------------

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------
/// \brief Allocate host and device memory.

void MagmaWrapper::malloc(MirroredBuf* buf, size_t dim0, size_t dim1,
   CEnv& env) {
  COMET_INSIST(buf);
  COMET_INSIST(dim0 + 1 >= 1 && dim1 + 1 >= 1);

#ifdef COMET_USE_MAGMA

  const size_t n = dim0 * dim1;

  if (use_minproduct_(env)) { //--------------------

    //typedef GMFloat Float_t;

    magma_minproduct_int_t magma_code = 0;

//#ifdef COMET_PLATFORM_CORI_GPU
    if (System::compute_capability() >= 800 ||
        !buf->is_compute_method_gpu()) {
      // WORKAROUND
      buf->h = (GMFloat*)::malloc(n*sizeof(GMFloat));
//#else
    } else {
      if (env.is_double_prec()) {
        magma_code = magma_minproduct_dmalloc_pinned((double**)&buf->h, n);
        COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS &&
          "Memory allocation error, possibly due to insufficient CPU memory.");
        COMET_INSIST(System::accel_last_call_succeeded());
      } else {
        magma_code = magma_minproduct_smalloc_pinned((float**)&buf->h, n);
        COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS &&
          "Memory allocation error, possibly due to insufficient CPU memory.");
        COMET_INSIST(System::accel_last_call_succeeded());
      }
    }
//#endif
    GMFloat_fill_nan((GMFloat*)buf->h, n);

    if (env.is_compute_method_gpu()) {
      if (env.is_double_prec()) {
        magma_code = magma_minproduct_dmalloc((double**)&buf->d, n);
        COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS &&
          "Memory allocation error, possibly due to insufficient GPU memory.");
        COMET_INSIST(System::accel_last_call_succeeded());
      } else {
        magma_code = magma_minproduct_smalloc((float**)&buf->d, n);
        COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS &&
          "Memory allocation error, possibly due to insufficient GPU memory.");
        COMET_INSIST(System::accel_last_call_succeeded());
      }
    }
    // TODO: ? fill GPU memory with NaNs

  } else if (use_mgemm4_(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    magma_mgemm4_int_t magma_code = 0;

//#ifdef COMET_PLATFORM_CORI_GPU
    if (System::compute_capability() >= 800 ||
        !buf->is_compute_method_gpu()) {
      // WORKAROUND
      buf->h = (Float_t*)::malloc(n*sizeof(Float_t));
//#else
    } else {
      magma_code = magma_mgemm4_zmalloc_pinned((Float_t**)&buf->h, n);
      COMET_INSIST(magma_code == MAGMA_mgemm4_SUCCESS &&
        "Memory allocation error, possibly due to insufficient CPU memory.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }
//#endif

    if (env.is_compute_method_gpu()) {
      magma_code = magma_mgemm4_zmalloc((Float_t**)&buf->d, n);
      COMET_INSIST(magma_code == MAGMA_mgemm4_SUCCESS &&
        "Memory allocation error, possibly due to insufficient GPU memory.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }

  } else if (use_mgemm2_(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    magma_mgemm2_int_t magma_code = 0;

//#ifdef COMET_PLATFORM_CORI_GPU
    if (System::compute_capability() >= 800 ||
        !buf->is_compute_method_gpu()) {
      // WORKAROUND
      buf->h = (Float_t*)::malloc(n*sizeof(Float_t));
//#else
    } else {
      magma_code = magma_mgemm2_zmalloc_pinned((Float_t**)&buf->h, n);
      COMET_INSIST(magma_code == MAGMA_mgemm2_SUCCESS &&
        "Memory allocation error, possibly due to insufficient CPU memory.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }
//#endif

    if (env.is_compute_method_gpu()) {
      magma_code = magma_mgemm2_zmalloc((Float_t**)&buf->d, n);
      COMET_INSIST(magma_code == MAGMA_mgemm2_SUCCESS &&
        "Memory allocation error, possibly due to insufficient GPU memory.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }

  } else if (use_mgemm3_(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    magma_mgemm3_int_t magma_code = 0;

//#ifdef COMET_PLATFORM_CORI_GPU
    if (System::compute_capability() >= 800 ||
        !buf->is_compute_method_gpu()) {
      // WORKAROUND
      buf->h = (Float_t*)::malloc(n*sizeof(Float_t));
//#else
    } else {
      magma_code = magma_mgemm3_zmalloc_pinned((Float_t**)&buf->h, n);
      COMET_INSIST(magma_code == MAGMA_mgemm3_SUCCESS &&
        "Memory allocation error, possibly due to insufficient CPU memory.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }
//#endif

    if (env.is_compute_method_gpu()) {
      magma_code = magma_mgemm3_zmalloc((Float_t**)&buf->d, n);
      COMET_INSIST(magma_code == MAGMA_mgemm3_SUCCESS &&
        "Memory allocation error, possibly due to insufficient GPU memory.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }

  } else if (use_mgemm5_(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    magma_mgemm5_int_t magma_code = 0;

//#ifdef COMET_PLATFORM_CORI_GPU
    if (System::compute_capability() >= 800 ||
        !buf->is_compute_method_gpu()) {
      // WORKAROUND
      buf->h = (Float_t*)::malloc(n*sizeof(Float_t));
//#else
    } else {
      magma_code = magma_mgemm5_zmalloc_pinned((Float_t**)&buf->h, n);
      COMET_INSIST(magma_code == MAGMA_mgemm5_SUCCESS &&
        "Memory allocation error, possibly due to insufficient CPU memory.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }
//#endif

    if (env.is_compute_method_gpu()) {
      magma_code = magma_mgemm5_zmalloc((Float_t**)&buf->d, n);
      COMET_INSIST(magma_code == MAGMA_mgemm5_SUCCESS &&
        "Memory allocation error, possibly due to insufficient GPU memory.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }

  } else { //--------------------

    COMET_INSIST_INTERFACE(&env, false && "Unimplemented method.");

  } // if //--------------------

  COMET_INSIST(buf->h &&
    "Invalid host pointer created, possibly due to insufficient CPU memory.");
  COMET_INSIST((buf->d || !env.is_compute_method_gpu()) &&
    "Invalid device pointer created, possibly due to insufficient GPU memory.");

#else

  COMET_INSIST(false && "Magma library unavailable for this build.");

#endif
}

//-----------------------------------------------------------------------------
/// \brief Free host and device memory.

void MagmaWrapper::free(MirroredBuf* buf, CEnv& env) {
  COMET_INSIST(buf);
  COMET_INSIST(! buf->is_alias);

#ifdef COMET_USE_MAGMA

  if (use_minproduct_(env)) { //--------------------

    magma_minproduct_int_t magma_code = 0;

//#ifdef COMET_PLATFORM_CORI_GPU
    if (System::compute_capability() >= 800 ||
        !buf->is_compute_method_gpu()) {
      // WORKAROUND
      ::free(buf->h);
//#else
    } else {
      magma_code = magma_minproduct_free_pinned(buf->h);
      COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS &&
        "Error in CPU memory free.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }
//#endif
    if (env.is_compute_method_gpu()) {
      magma_code = magma_minproduct_free(buf->d);
      COMET_INSIST(magma_code == MAGMA_minproduct_SUCCESS &&
        "Error in GPU memory free.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }

  } else if (use_mgemm4_(env)) { //--------------------

    magma_mgemm4_int_t magma_code = 0;

//#ifdef COMET_PLATFORM_CORI_GPU
    if (System::compute_capability() >= 800 ||
        !buf->is_compute_method_gpu()) {
      // WORKAROUND
      ::free(buf->h);
//#else
    } else {
      magma_code = magma_mgemm4_free_pinned(buf->h);
      COMET_INSIST(magma_code == MAGMA_mgemm4_SUCCESS &&
               "Error in magma_mgemm4_free_pinned.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }
//#endif
    if (env.is_compute_method_gpu()) {
      magma_code = magma_mgemm4_free(buf->d);
      COMET_INSIST(magma_code == MAGMA_mgemm4_SUCCESS &&
        "Error in GPU memory free.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }

  } else if (use_mgemm2_(env)) { //--------------------

    magma_mgemm2_int_t magma_code = 0;

//#ifdef COMET_PLATFORM_CORI_GPU
    if (System::compute_capability() >= 800 ||
        !buf->is_compute_method_gpu()) {
      // WORKAROUND
      ::free(buf->h);
//#else
    } else {
      magma_code = magma_mgemm2_free_pinned(buf->h);
      COMET_INSIST(magma_code == MAGMA_mgemm2_SUCCESS &&
        "Error in CPU memory free.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }
//#endif
    if (env.is_compute_method_gpu()) {
      magma_code = magma_mgemm2_free(buf->d);
      COMET_INSIST(magma_code == MAGMA_mgemm2_SUCCESS &&
        "Error in GPU memory free.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }

  } else if (use_mgemm3_(env)) { //--------------------

    magma_mgemm3_int_t magma_code = 0;

//#ifdef COMET_PLATFORM_CORI_GPU
    if (System::compute_capability() >= 800 ||
        !buf->is_compute_method_gpu()) {
      // WORKAROUND
      ::free(buf->h);
//#else
    } else {
      magma_code = magma_mgemm3_free_pinned(buf->h);
      COMET_INSIST(magma_code == MAGMA_mgemm3_SUCCESS &&
        "Error in CPU memory free.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }
//#endif
    if (env.is_compute_method_gpu()) {
      magma_code = magma_mgemm3_free(buf->d);
      COMET_INSIST(magma_code == MAGMA_mgemm3_SUCCESS &&
        "Error in GPU memory free.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }

  } else if (use_mgemm5_(env)) { //--------------------

    magma_mgemm5_int_t magma_code = 0;

//#ifdef COMET_PLATFORM_CORI_GPU
    if (System::compute_capability() >= 800 ||
        !buf->is_compute_method_gpu()) {
      // WORKAROUND
      ::free(buf->h);
//#else
    } else {
      magma_code = magma_mgemm5_free_pinned(buf->h);
      COMET_INSIST(magma_code == MAGMA_mgemm5_SUCCESS &&
        "Error in CPU memory free.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }
//#endif
    if (env.is_compute_method_gpu()) {
      magma_code = magma_mgemm5_free(buf->d);
      COMET_INSIST(magma_code == MAGMA_mgemm5_SUCCESS &&
        "Error in GPU memory free.");
      COMET_INSIST(System::accel_last_call_succeeded());
    }

  } else { //--------------------

    COMET_INSIST_INTERFACE(&env, false && "Unimplemented method.");

  } // if //--------------------

#else

  COMET_INSIST(false && "Magma library unavailable for this build.");

#endif
}

//-----------------------------------------------------------------------------
/// \brief Set matrix to zero.

void MagmaWrapper::set_matrix_zero_start(MirroredBuf* buf, CEnv& env) {
  COMET_INSIST(buf);

  if (!env.is_compute_method_gpu())
    return;

#ifdef COMET_USE_MAGMA

  const size_t mat_dim1 = buf->dim0;
  const size_t mat_dim2 = buf->dim1;

  // ISSUE: these MAGMA routines don't return an error code.

  if (use_minproduct_(env) && env.is_double_prec()) { //--------------------

    COMET_INSIST(env.queue_compute_.magma_queue);
    magma_minproductblas_dlaset
      (Magma_minproductFull, mat_dim1, mat_dim2, (double)0, (double)0,
       //(double*)buf->d, mat_dim1 PREPEND_COMMA(MagmaQueue<MagmaCloneId::MINPRODUCT>().compute(env)));
       (double*)buf->d, mat_dim1 PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MINPRODUCT>>()));
    COMET_INSIST(System::accel_last_call_succeeded());

  } else if (use_minproduct_(env)) { //--------------------

    COMET_INSIST(env.queue_compute_.magma_queue);
    magma_minproductblas_slaset
      (Magma_minproductFull, mat_dim1, mat_dim2, (float)0, (float)0,
       //(float*)buf->d, mat_dim1 PREPEND_COMMA(MagmaQueue<MagmaCloneId::MINPRODUCT>().compute(env)));
       (float*)buf->d, mat_dim1 PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MINPRODUCT>>()));
    COMET_INSIST(System::accel_last_call_succeeded());

  } else if (use_mgemm4_(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    Float_t zero = {0, 0};

    COMET_INSIST(env.queue_compute_.magma_queue);
    magma_mgemm4blas_zlaset(Magma_mgemm4Full, mat_dim1, mat_dim2, zero, zero,
      //(Float_t*)buf->d, mat_dim1 PREPEND_COMMA(MagmaQueue<MagmaCloneId::MGEMM4>().compute(env)));
      (Float_t*)buf->d, mat_dim1 PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MGEMM4>>()));
    COMET_INSIST(System::accel_last_call_succeeded());

  } else if (use_mgemm2_(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    Float_t zero = {0, 0};

    COMET_INSIST(env.queue_compute_.magma_queue);
    magma_mgemm2blas_zlaset(Magma_mgemm2Full, mat_dim1, mat_dim2, zero, zero,
      //(Float_t*)buf->d, mat_dim1 PREPEND_COMMA(MagmaQueue<MagmaCloneId::MGEMM2>().compute(env)));
      (Float_t*)buf->d, mat_dim1 PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MGEMM2>>()));
    COMET_INSIST(System::accel_last_call_succeeded());

  } else if (use_mgemm3_(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    Float_t zero = {0, 0};

    COMET_INSIST(env.queue_compute_.magma_queue);
    magma_mgemm3blas_zlaset(Magma_mgemm3Full, mat_dim1, mat_dim2, zero, zero,
      //(Float_t*)buf->d, mat_dim1 PREPEND_COMMA(MagmaQueue<MagmaCloneId::MGEMM3>().compute(env)));
      (Float_t*)buf->d, mat_dim1 PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MGEMM3>>()));
    COMET_INSIST(System::accel_last_call_succeeded());

  } else if (use_mgemm5_(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    Float_t zero = {0, 0};

    COMET_INSIST(env.queue_compute_.magma_queue);
    magma_mgemm5blas_zlaset(Magma_mgemm5Full, mat_dim1, mat_dim2, zero, zero,
      //(Float_t*)buf->d, mat_dim1 PREPEND_COMMA(MagmaQueue<MagmaCloneId::MGEMM5>().compute(env)));
      (Float_t*)buf->d, mat_dim1 PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MGEMM5>>()));
    COMET_INSIST(System::accel_last_call_succeeded());

  } else { //--------------------

      COMET_INSIST_INTERFACE(&env, false && "Unimplemented method.");

  } // if //--------------------

#else

  COMET_INSIST(false && "Magma library unavailable for this build.");

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------
/// \brief GEMM on block of matrix, start.

void MagmaWrapper::gemm_block_start(size_t m, size_t n, size_t k,
  const void* matA, size_t ldda, const void* matB, size_t lddb,
  void* matC, size_t lddc, bool is_beta_one, CEnv& env) {
  COMET_INSIST(matA && matB && matC);
  COMET_INSIST(env.is_compute_method_gpu());

#ifdef COMET_USE_MAGMA

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

  // ISSUE: these MAGMA routines don't return an error code.

  if (use_minproduct_(env)) { //--------------------

    const GMFloat alpha = 1;
    const GMFloat beta = is_beta_one ? 1 : 0;

    typedef magma_minproduct_int_t Int_t;

    const Int_t m_ = safe_cast<Int_t>(m);
    const Int_t n_ = safe_cast<Int_t>(n);
    const Int_t k_ = safe_cast<Int_t>(k);
    const Int_t ldda_ = safe_cast<Int_t>(ldda);
    const Int_t lddb_ = safe_cast<Int_t>(lddb);
    const Int_t lddc_ = safe_cast<Int_t>(lddc);

    if (env.is_double_prec()) {
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
        //lddc_ PREPEND_COMMA(MagmaQueue<MagmaCloneId::MINPRODUCT>().compute(env)));
        lddc_ PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MINPRODUCT>>()));
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
        //lddc_ PREPEND_COMMA(MagmaQueue<MagmaCloneId::MINPRODUCT>().compute(env)));
        lddc_ PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MINPRODUCT>>()));
      COMET_INSIST(System::accel_last_call_succeeded() &&
               "Failure in call to magma_minproductblas_sgemm.");
    }

    env.ops_local_inc(2 * m * (double)n * (double)k);

  } else if (use_mgemm4_(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const Float_t alpha = {1, 0};
    const Float_t beta = is_beta_one ? one : zero;

    typedef magma_mgemm4_int_t Int_t;

    const Int_t m_ = safe_cast<Int_t>(m);
    const Int_t n_ = safe_cast<Int_t>(n);
    const Int_t k_ = safe_cast<Int_t>(k);
    const Int_t ldda_ = safe_cast<Int_t>(ldda);
    const Int_t lddb_ = safe_cast<Int_t>(lddb);
    const Int_t lddc_ = safe_cast<Int_t>(lddc);

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
      //lddc_ PREPEND_COMMA(MagmaQueue<MagmaCloneId::MGEMM4>().compute(env)));
      lddc_ PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MGEMM4>>()));
    COMET_INSIST(System::accel_last_call_succeeded() &&
             "Failure in call to magma_mgemm4blas_zgemm.");

  } else if (use_mgemm2_(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const Float_t alpha = {1, 0};
    const Float_t beta = is_beta_one ? one : zero;

    typedef magma_mgemm2_int_t Int_t;

    const Int_t m_ = safe_cast<Int_t>(m);
    const Int_t n_ = safe_cast<Int_t>(n);
    const Int_t k_ = safe_cast<Int_t>(k);
    const Int_t ldda_ = safe_cast<Int_t>(ldda);
    const Int_t lddb_ = safe_cast<Int_t>(lddb);
    const Int_t lddc_ = safe_cast<Int_t>(lddc);

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
      //lddc_ PREPEND_COMMA(MagmaQueue<MagmaCloneId::MGEMM2>().compute(env)));
      lddc_ PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MGEMM2>>()));
    COMET_INSIST(System::accel_last_call_succeeded() &&
             "Failure in call to magma_mgemm2blas_zgemm.");

  } else if (use_mgemm3_(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const Float_t alpha = {1, 0};
    const Float_t beta = is_beta_one ? one : zero;

    typedef magma_mgemm3_int_t Int_t;

    const Int_t m_ = safe_cast<Int_t>(m);
    const Int_t n_ = safe_cast<Int_t>(n);
    const Int_t k_ = safe_cast<Int_t>(k);
    const Int_t ldda_ = safe_cast<Int_t>(ldda);
    const Int_t lddb_ = safe_cast<Int_t>(lddb);
    const Int_t lddc_ = safe_cast<Int_t>(lddc);

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
      //lddc_ PREPEND_COMMA(MagmaQueue<MagmaCloneId::MGEMM3>().compute(env)));
      lddc_ PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MGEMM3>>()));
    COMET_INSIST(System::accel_last_call_succeeded() &&
             "Failure in call to magma_mgemm3blas_zgemm.");

  } else if (use_mgemm5_(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    typedef magma_mgemm5_int_t Int_t;

    const Int_t m_ = safe_cast<Int_t>(m);
    const Int_t n_ = safe_cast<Int_t>(n);
    const Int_t k_ = safe_cast<Int_t>(k);
    const Int_t ldda_ = safe_cast<Int_t>(ldda);
    const Int_t lddb_ = safe_cast<Int_t>(lddb);
    const Int_t lddc_ = safe_cast<Int_t>(lddc);

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
      //lddc_ PREPEND_COMMA(MagmaQueue<MagmaCloneId::MGEMM5>().compute(env)));
      lddc_ PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MGEMM5>>()));
    COMET_INSIST(System::accel_last_call_succeeded() &&
             "Failure in call to magma_mgemm5blas_zgemm.");

  } else { //--------------------

      COMET_INSIST_INTERFACE(&env, false && "Unimplemented method.");

  } // if //--------------------

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------
/// \brief GEMM on matrix, start.

void MagmaWrapper::gemm_start(size_t m, size_t n, size_t k,
  const void* matA, size_t ldda, const void* matB, size_t lddb,
  void* matC, size_t lddc, CEnv& env) {
  COMET_INSIST(matA && matB && matC);
  COMET_INSIST(env.is_compute_method_gpu());

#if defined COMET_USE_SEMIRING

  if (use_minproduct_(env)) { //--------------------
  //if (use_minproduct_(env) && !env.is_double_prec()) { //--------------------

    const auto m_ = safe_cast<int>(m);
    const auto n_ = safe_cast<int>(n);
    const auto k_ = safe_cast<int>(k);
    const auto ldda_ = safe_cast<std::size_t>(ldda);
    const auto lddb_ = safe_cast<std::size_t>(lddb);
    const auto lddc_ = safe_cast<std::size_t>(lddc);

    if (env.is_double_prec()) {

      //COMET_INSIST(m_ == n_);
      //COMET_INSIST(m_ == lddc_);
      //COMET_INSIST(k_ == lddb_);
      //COMET_INSIST(k_ == ldda_);

      srPlusMinFP64TN(m_, n_, k_,
                     (double*)static_cast<const double*>(matA), ldda_,
                     (double*)static_cast<const double*>(matB), lddb_,
                     static_cast<double*>(matC), lddc_,
                     env.stream_compute());

//printf("m %i n %i k %i lda %zu ldb %zu ldc %zu\n", m_, n_, k_, ldda_, lddb_, lddc_);

#if 0
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
        //lddc_ PREPEND_COMMA(MagmaQueue<MagmaCloneId::MINPRODUCT>().compute(env)));
        lddc_ PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MINPRODUCT>>()));
#endif

      COMET_INSIST(System::accel_last_call_succeeded() &&
               "Failure in call to srPlusMinFP64TN.");

    } else {

      srPlusMinFP32TN(m_, n_, k_,
                     (float*)static_cast<const float*>(matA), ldda_,
                     (float*)static_cast<const float*>(matB), lddb_,
                     static_cast<float*>(matC), lddc_,
                     env.stream_compute());

#if 0
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
        //lddc_ PREPEND_COMMA(MagmaQueue<MagmaCloneId::MINPRODUCT>().compute(env)));
        lddc_ PREPEND_COMMA(env.queue_compute<MagmaQueue<MagmaCloneId::MINPRODUCT>>()));
#endif

      COMET_INSIST(System::accel_last_call_succeeded() &&
               "Failure in call to srPlusMinFP32TN.");

    }

    env.ops_local_inc(2 * m * (double)n * (double)k);

    return; // exit here, no need to call MAGMA
  }

#endif // COMET_USE_SEMIRING

#if defined COMET_USE_MAGMA

  // The purpose of this code is to workaround the magma size
  // limitation (for non CUBLAS failover) by doing gemm in blocks.

  const size_t rows = k;
  const size_t cols_A = m;
  const size_t cols_B = n;

  const size_t elt_size =
    env.metric_type() == MetricType::CZEK ? sizeof(GMFloat) :
   (env.metric_type() == MetricType::CCC && env.sparse()) ?
                                         sizeof(magma_mgemm4DoubleComplex) :
   (env.metric_type() == MetricType::CCC &&
    env.num_way() == NumWay::_2) ? sizeof(magma_mgemm2DoubleComplex) :
   (env.metric_type() == MetricType::CCC &&
    env.num_way() == NumWay::_3) ? sizeof(magma_mgemm3DoubleComplex) :
   (env.metric_type() == MetricType::DUO) ?
                                         sizeof(magma_mgemm5DoubleComplex) : 0;
  COMET_INSIST(elt_size > 0 && "Error in gemm block calculation.");

  const size_t align_factor = 128 / elt_size;
  const size_t max_elts = (1 << 27) - 512;

  // TODO: can we improve aspect ratios of submatrices.
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

        MagmaWrapper::gemm_block_start(cols_A_this, cols_B_this, rows_this,
          A_this, ldda, B_this, lddb, C_this, lddc, row_base > 0,  env);
      }
    }
  }

#endif // COMET_USE_MAGMA

}

//-----------------------------------------------------------------------------
/// \brief Start transfer of matrix to GPU.

void MagmaWrapper::set_matrix_start(MirroredBuf* buf, CEnv& env) {
  COMET_INSIST(buf);

  if (!env.is_compute_method_gpu())
    return;

  // NOTE: these MAGMA routines don't return an error code.

#ifdef COMET_USE_MAGMA

  //typedef MagmaCloneId MW;

  if (use_minproduct_(env) && env.is_double_prec()) { //--------------------

    typedef double Float_t;

    magma_minproduct_dsetmatrix_async(
      buf->dim0, buf->dim1, (Float_t*)buf->h, buf->dim0,
      (Float_t*)buf->d, buf->dim0,
      env.queue_togpu<MagmaQueue<MagmaCloneId::MINPRODUCT>>());
      //MagmaQueue<MW::MINPRODUCT>().togpu(env));
      //MagmaQueue<MW::MINPRODUCT>().togpu(env));
      //(Float_t*)buf->d, buf->dim0, env.stream_togpu());
    COMET_INSIST(System::accel_last_call_succeeded());

  } else if (use_minproduct_(env)) { //--------------------

    typedef float Float_t;

    magma_minproduct_ssetmatrix_async(
      buf->dim0, buf->dim1, (Float_t*)buf->h, buf->dim0,
      (Float_t*)buf->d, buf->dim0,
      env.queue_togpu<MagmaQueue<MagmaCloneId::MINPRODUCT>>());
      //MagmaQueue<MW::MINPRODUCT>().togpu(env));
      //(Float_t*)buf->d, buf->dim0, env.stream_togpu());
    COMET_INSIST(System::accel_last_call_succeeded());

//#   if ! defined COMET_USE_MAGMA_V2
//#   endif

  } else if (use_mgemm4_(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    magma_mgemm4_zsetmatrix_async(buf->dim0, buf->dim1, (Float_t*)buf->h,
                                  buf->dim0, (Float_t*)buf->d, buf->dim0,
                                  env.queue_togpu<MagmaQueue<MagmaCloneId::MGEMM4>>());
                                  //MagmaQueue<MW::MGEMM4>().togpu(env));
                                  //env.stream_togpu());
    COMET_INSIST(System::accel_last_call_succeeded());

  } else if (use_mgemm2_(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    magma_mgemm2_zsetmatrix_async(buf->dim0, buf->dim1, (Float_t*)buf->h,
                                  buf->dim0, (Float_t*)buf->d, buf->dim0,
                                  env.queue_togpu<MagmaQueue<MagmaCloneId::MGEMM2>>());
                                  //MagmaQueue<MW::MGEMM2>().togpu(env));
                                  //env.stream_togpu());
    COMET_INSIST(System::accel_last_call_succeeded());

  } else if (use_mgemm3_(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    magma_mgemm3_zsetmatrix_async(buf->dim0, buf->dim1, (Float_t*)buf->h,
                                  buf->dim0, (Float_t*)buf->d, buf->dim0,
                                  env.queue_togpu<MagmaQueue<MagmaCloneId::MGEMM3>>());
                                  //MagmaQueue<MW::MGEMM3>().togpu(env));
                                  //env.stream_togpu());

  } else if (use_mgemm5_(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    magma_mgemm5_zsetmatrix_async(buf->dim0, buf->dim1, (Float_t*)buf->h,
                                  buf->dim0, (Float_t*)buf->d, buf->dim0,
                                  env.queue_togpu<MagmaQueue<MagmaCloneId::MGEMM5>>());
                                  //MagmaQueue<MW::MGEMM5>().togpu(env));
                                  //env.stream_togpu());
    COMET_INSIST(System::accel_last_call_succeeded());

  } else { //--------------------

      COMET_INSIST_INTERFACE(&env, false && "Unimplemented method.");

  } // if //--------------------

#else

  COMET_INSIST(false && "Magma library unavailable for this build.");

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------
/// \brief Wait transfer of matrix to GPU.

void MagmaWrapper::set_matrix_wait(CEnv& env) {
  env.stream_synchronize(env.stream_togpu());
}

//-----------------------------------------------------------------------------
/// \brief Start transfer of matrix from GPU.

void MagmaWrapper::get_matrix_start(MirroredBuf* buf, CEnv& env) {
  COMET_INSIST(buf);

  if (!env.is_compute_method_gpu())
    return;

  // NOTE: these MAGMA routines don't return an error code.

#ifdef COMET_USE_MAGMA

  //typedef MagmaCloneId MW;

  if (use_minproduct_(env) && env.is_double_prec()) { //--------------------

    typedef double Float_t;

    magma_minproduct_dgetmatrix_async(
      buf->dim0, buf->dim1, (Float_t*)buf->d, buf->dim0,
      (Float_t*)buf->h, buf->dim0,
      env.queue_fromgpu<MagmaQueue<MagmaCloneId::MINPRODUCT>>());
      //MagmaQueue<MW::MINPRODUCT>().fromgpu(env));
      //env.stream_fromgpu());
    COMET_INSIST(System::accel_last_call_succeeded());

  } else if (use_minproduct_(env)) { //--------------------

    typedef float Float_t;

    magma_minproduct_sgetmatrix_async(
      buf->dim0, buf->dim1, (Float_t*)buf->d, buf->dim0,
      (Float_t*)buf->h, buf->dim0,
      env.queue_fromgpu<MagmaQueue<MagmaCloneId::MINPRODUCT>>());
      //MagmaQueue<MW::MINPRODUCT>().fromgpu(env));
      //env.stream_fromgpu());
    COMET_INSIST(System::accel_last_call_succeeded());

  } else if (use_mgemm4_(env)) { //--------------------

    typedef magma_mgemm4DoubleComplex Float_t;

    magma_mgemm4_zgetmatrix_async(buf->dim0, buf->dim1, (Float_t*)buf->d,
                                  buf->dim0, (Float_t*)buf->h, buf->dim0,
                                  env.queue_fromgpu<MagmaQueue<MagmaCloneId::MGEMM4>>());
                                  //MagmaQueue<MW::MGEMM4>().fromgpu(env));
                                  //env.stream_fromgpu());
    COMET_INSIST(System::accel_last_call_succeeded());

  } else if (use_mgemm2_(env)) { //--------------------

    typedef magma_mgemm2DoubleComplex Float_t;

    magma_mgemm2_zgetmatrix_async(buf->dim0, buf->dim1, (Float_t*)buf->d,
                                  buf->dim0, (Float_t*)buf->h, buf->dim0,
                                  env.queue_fromgpu<MagmaQueue<MagmaCloneId::MGEMM2>>());
                                  //MagmaQueue<MW::MGEMM2>().fromgpu(env));
                                  //env.stream_fromgpu());
    COMET_INSIST(System::accel_last_call_succeeded());

  } else if (use_mgemm3_(env)) { //--------------------

    typedef magma_mgemm3DoubleComplex Float_t;

    magma_mgemm3_zgetmatrix_async(buf->dim0, buf->dim1, (Float_t*)buf->d,
                                  buf->dim0, (Float_t*)buf->h, buf->dim0,
                                  env.queue_fromgpu<MagmaQueue<MagmaCloneId::MGEMM3>>());
                                  //MagmaQueue<MW::MGEMM3>().fromgpu(env));
                                  //env.stream_fromgpu());
    COMET_INSIST(System::accel_last_call_succeeded());

  } else if (use_mgemm5_(env)) { //--------------------

    typedef magma_mgemm5DoubleComplex Float_t;

    magma_mgemm5_zgetmatrix_async(buf->dim0, buf->dim1, (Float_t*)buf->d,
                                  buf->dim0, (Float_t*)buf->h, buf->dim0,
                                  env.queue_fromgpu<MagmaQueue<MagmaCloneId::MGEMM5>>());
                                  //MagmaQueue<MW::MGEMM5>().fromgpu(env));
                                  //env.stream_fromgpu());
    COMET_INSIST(System::accel_last_call_succeeded());

  } else { //--------------------

      COMET_INSIST_INTERFACE(&env, false && "Unimplemented method.");

  } // if //--------------------

#else

  COMET_INSIST(false && "Magma library unavailable for this build.");

#endif // COMET_USE_MAGMA
}

//-----------------------------------------------------------------------------
/// \brief Wait transfer of matrix from GPU.

void MagmaWrapper::get_matrix_wait(CEnv& env) {
  env.stream_synchronize(env.stream_fromgpu());
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
