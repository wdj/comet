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
/// \brief Static-only wrapper class to access MAGMA library variants.

struct MagmaWrapper {

  enum {NONE = 0,
        MINPRODUCT = 1,
        MGEMM2 = 2,
        MGEMM3 = 3,
        MGEMM4 = 4,
        MGEMM5 = 5};

  static void initialize(CEnv& env);
  static void finalize(CEnv& env);
  static void malloc(MirroredBuf* buf, size_t dim0, size_t dim1, CEnv& env);
  static void free(MirroredBuf* buf, CEnv& env);
  static void set_matrix_start(MirroredBuf* buf, CEnv& env);
  static void set_matrix_wait(CEnv& env);
  static void get_matrix_start(MirroredBuf* buf, CEnv& env);
  static void get_matrix_wait(CEnv& env);
  static void gemm_start(size_t m, size_t n, size_t k,
    const void* matA, size_t ldda, const void* matB, size_t lddb,
    void* matC, size_t lddc, CEnv& env);
  static void set_matrix_zero_start(MirroredBuf* buf, CEnv& env);

private:

  static void gemm_block_start(size_t m, size_t n, size_t k,
    const void* matA, size_t ldda, const void* matB, size_t lddb,
    void* matC, size_t lddc, bool is_beta_one, CEnv& env);

  static bool use_minproduct_(CEnv& env) {
    return env.metric_type() == MetricType::CZEK;
  }

  static bool use_mgemm2_(CEnv& env) {
    return env.metric_type() == MetricType::CCC &&
           env.num_way() == NumWay::_2 && ! env.sparse();
  }

  static bool use_mgemm3_(CEnv& env) {
    return env.metric_type() == MetricType::CCC &&
           env.num_way() == NumWay::_3 && ! env.sparse();
  }

  static bool use_mgemm4_(CEnv& env) {
    return env.metric_type() == MetricType::CCC && env.sparse();
  }

  static bool use_mgemm5_(CEnv& env) {
    return env.metric_type() == MetricType::DUO;
  }

}; // MagmaWrapper

//=============================================================================

template<int>
struct MagmaQueueHelper {
  typedef CEnv::Stream_t queue_t;
  static void create(CEnv::Stream_t stream, queue_t& queue) {
    queue = stream;
  }
  static void destroy(queue_t& queue) {}
};

template<>
struct MagmaQueueHelper<MagmaWrapper::MINPRODUCT>
   : public MagmaQueueHelper<MagmaWrapper::NONE> {
# ifdef COMET_USE_MAGMA
# if defined COMET_USE_MAGMA_V2
    typedef magma_minproduct_queue_t queue_t;
    static void create(CEnv::Stream_t stream, queue_t& queue) {
      magma_minproduct_queue_create_from_hip(0, stream, NULL, NULL, &queue);
    }
    static void destroy(queue_t& queue) {
      magma_minproduct_queue_destroy(queue);
    }
# endif
# endif
};

template<>
struct MagmaQueueHelper<MagmaWrapper::MGEMM2>
   : public MagmaQueueHelper<MagmaWrapper::NONE> {
# ifdef COMET_USE_MAGMA
# if defined COMET_USE_MAGMA_V2
    typedef magma_mgemm2_queue_t queue_t;
    static void create(CEnv::Stream_t stream, queue_t& queue) {
      magma_mgemm2_queue_create_from_hip(0, stream, NULL, NULL, &queue);
    }
    static void destroy(queue_t& queue) {
      magma_mgemm2_queue_destroy(queue);
    }
# endif
# endif
};

template<>
struct MagmaQueueHelper<MagmaWrapper::MGEMM3>
   : public MagmaQueueHelper<MagmaWrapper::NONE> {
# ifdef COMET_USE_MAGMA
# if defined COMET_USE_MAGMA_V2
    typedef magma_mgemm3_queue_t queue_t;
    static void create(CEnv::Stream_t stream, queue_t& queue) {
      magma_mgemm3_queue_create_from_hip(0, stream, NULL, NULL, &queue);
    }
    static void destroy(queue_t& queue) {
      magma_mgemm3_queue_destroy(queue);
    }
# endif
# endif
};

template<>
struct MagmaQueueHelper<MagmaWrapper::MGEMM4>
   : public MagmaQueueHelper<MagmaWrapper::NONE> {
# ifdef COMET_USE_MAGMA
# if defined COMET_USE_MAGMA_V2
    typedef magma_mgemm4_queue_t queue_t;
    static void create(CEnv::Stream_t stream, queue_t& queue) {
      magma_mgemm4_queue_create_from_hip(0, stream, NULL, NULL, &queue);
    }
    static void destroy(queue_t& queue) {
      magma_mgemm4_queue_destroy(queue);
    }
# endif
# endif
};

template<>
struct MagmaQueueHelper<MagmaWrapper::MGEMM5>
   : public MagmaQueueHelper<MagmaWrapper::NONE> {
# ifdef COMET_USE_MAGMA
# if defined COMET_USE_MAGMA_V2
    typedef magma_mgemm5_queue_t queue_t;
    static void create(CEnv::Stream_t stream, queue_t& queue) {
      magma_mgemm5_queue_create_from_hip(0, stream, NULL, NULL, &queue);
    }
    static void destroy(queue_t& queue) {
      magma_mgemm5_queue_destroy(queue);
    }
# endif
# endif
};

//-----------------------------------------------------------------------------

template<int CLONE_ID>
class MagmaQueue {
public:

  MagmaQueue()
    : is_queue_initialized_(false) {
  }

  ~MagmaQueue() {
    COMET_INSIST(is_queue_initialized_);
    MagmaQueueHelper<CLONE_ID>::destroy(queue_);
  }

  typename MagmaQueueHelper<CLONE_ID>::queue_t togpu(CEnv& env) {
    COMET_INSIST(!is_queue_initialized_);
    MagmaQueueHelper<CLONE_ID>::create(env.stream_togpu(), queue_);
    is_queue_initialized_ = true;
    return queue_;
  }

  typename MagmaQueueHelper<CLONE_ID>::queue_t fromgpu(CEnv& env) {
    COMET_INSIST(!is_queue_initialized_);
    MagmaQueueHelper<CLONE_ID>::create(env.stream_fromgpu(), queue_);
    is_queue_initialized_ = true;
    return queue_;
  }

  typename MagmaQueueHelper<CLONE_ID>::queue_t compute(CEnv& env) {
    COMET_INSIST(!is_queue_initialized_);
    MagmaQueueHelper<CLONE_ID>::create(env.stream_compute(), queue_);
    is_queue_initialized_ = true;
    return queue_;
  }

private:

  bool is_queue_initialized_;
  typename MagmaQueueHelper<CLONE_ID>::queue_t queue_;

};

#if 0
template<int>
class MagmaQueue {
public:

  MagmaQueue()
    : is_queue_initialized_(false) {
  }

  ~MagmaQueue();

  CEnv::Stream_t togpu(CEnv& env) {
    return env.stream_togpu();
  }

protected:

  bool is_queue_initialized_;
};

//-----------------------------------------------------------------------------

template<>
class MagmaQueue<MagmaWrapper::MINPRODUCT> : public MagmaQueue<MagmaWrapper::NONE> {
# if defined COMET_USE_MAGMA_V2
public:

  magma_minproduct_queue_t togpu(CEnv& env) {
    COMET_INSIST(!is_queue_initialized_);
    magma_minproduct_queue_create_from_hip(0, env.stream_togpu(), NULL, NULL, &queue_);
    is_queue_initialized_ = true;
    return queue_;
  }

  ~MagmaQueue() {
    if (is_queue_initialized_)
      magma_minproduct_queue_destroy(queue_);
  }

private:

  magma_minproduct_queue_t queue_;
# endif
};
#endif

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_MAGMA_WRAPPER_HH_

//-----------------------------------------------------------------------------
