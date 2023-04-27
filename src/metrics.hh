//-----------------------------------------------------------------------------
/*!
 * \file   metrics.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Class to manage the calculated metrics output from the methods.
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

#ifndef _COMET_METRICS_HH_
#define _COMET_METRICS_HH_

#include "cstddef"
#include "cstdint"
#include "math.h"

#include "env.hh"
#include "mirrored_buf.hh"
#include "histograms.hh"
#include "decomp_mgr.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/*!
 * \class MetricsArrayId
 * \brief Helper class to designate which array (metrics proper, sums, counts.
 *
 */
//-----------------------------------------------------------------------------
class MetricsArrayId {

//  template<int> class Enum2Type {};

public:

  enum {M = 0, // Matrics
        S = 1, // Sums
        C = 2, // Counts
        O = 3, // coOrdinates
        NUM = 4};

  struct M_t {
    typedef void* __restrict__ Data_t;
  };

  struct S_t {
    typedef void* __restrict__ Data_t;
  };

  struct C_t {
    typedef void* __restrict__ Data_t;
  };

  struct O_t {
    typedef MetricItemCoords_t* __restrict__ Data_t;
  };


//  typedef Enum2Type<M> M_t;
//  typedef Enum2Type<S> S_t;
//  typedef Enum2Type<C> C_t;
//  typedef Enum2Type<O> O_t;
};

//-----------------------------------------------------------------------------
/// \brief Helper class for metrics memory.

// Forward declartion.
struct GMMetrics;

class MetricsMem {

  typedef MetricsArrayId MAI;

public:

  struct Array {
    static MAI::M_t::Data_t& data(MAI::M_t, MetricsMem* mm) {return mm->data_;}
    static MAI::S_t::Data_t& data(MAI::S_t, MetricsMem* mm) {return mm->data_S_;}
    static MAI::C_t::Data_t& data(MAI::C_t, MetricsMem* mm) {return mm->data_C_;}
    static MAI::O_t::Data_t& data(MAI::O_t, MetricsMem* mm) {return mm->coords_;}

    static size_t& size(MAI::M_t, MetricsMem* mm) {return mm->data_size_;}
    static size_t& size(MAI::S_t, MetricsMem* mm) {return mm->data_S_size_;}
    static size_t& size(MAI::C_t, MetricsMem* mm) {return mm->data_C_size_;}
    static size_t& size(MAI::O_t, MetricsMem* mm) {return mm->coords_size_;}
  };

#if 0
  template<int MAI> struct Array;

  template<> struct Array<MAI::M> {
    typedef void* __restrict__ Data_t;
    static Data_t& data(MetricsMem* mm) {return mm->data_;}
    static size_t& size(MetricsMem* mm) {return mm->data_size_;}
  };

  template<> struct Array<MAI::S> {
    typedef void* __restrict__ Data_t;
    static Data_t& data(MetricsMem* mm) {return mm->data_S_;}
    static size_t& size(MetricsMem* mm) {return mm->data_S_size_;}
  };

  template<> struct Array<MAI::C> {
    typedef void* __restrict__ Data_t;
    static Data_t& data(MetricsMem* mm) {return mm->data_C_;}
    static size_t& size(MetricsMem* mm) {return mm->data_C_size_;}
  };

  template<> struct Array<MAI::COORDS> {
    typedef MetricItemCoords_t* __restrict__ Data_t;
    static Data_t& data(MetricsMem* mm) {return mm->coords_;}
    static size_t& size(MetricsMem* mm) {return mm->coords_size_;}
  };
#endif

public:

  MetricsMem(CEnv& env);
  ~MetricsMem();
  void deallocate();

  MAI::M_t m_selector;
  MAI::S_t s_selector;
  MAI::C_t c_selector;
  MAI::O_t o_selector;

  // \brief Re/allocate one of the four arrays, as needed.

  template<class MAI_t>
  typename MAI_t::Data_t
  array_malloc(MAI_t mai_selector, size_t data_size) {

    if (!env_.is_proc_active())
      return NULL;

    COMET_INSIST(is_allocated_);

    typedef typename MAI_t::Data_t Data_t;

    Data_t& data_ref = Array::data(mai_selector, this);
    size_t& data_size_ref = Array::size(mai_selector, this);

    if (!data_ref || data_size > data_size_ref) {
      if (data_ref)
        utils::free(data_ref, data_size_ref, env_);
      data_ref = static_cast<Data_t>(utils::malloc(data_size, env_));
      data_size_ref = data_size;
    }

    return data_ref;
  }

#if 0
  template<int MAI>
  typename Array<MAI>::Data_t
  array_malloc(size_t data_size) {

    if (!env_.is_proc_active())
      return NULL;

    COMET_INSIST(is_allocated_);

    typename Array<MAI>::Data_t& data_ref = Array<MAI>::data(this);
    size_t& data_size_ref = Array<MAI>::size(this);

    if (!data_ref || data_size > data_size_ref) {
      if (data_ref)
        utils::free(data_ref, data_size_ref, env_);
      data_ref = static_cast<typename Array<MAI>::Data_t>(
        utils::malloc(data_size, env_));
      data_size_ref = data_size;
    }

    return data_ref;
  }
#endif

  // \brief Deallocate one of the four arrays.

  template<class MAI_t>
  void array_free(MAI_t mai_selector) {

    typedef typename MAI_t::Data_t Data_t;

    Data_t& data_ref = Array::data(mai_selector, this);
    size_t data_size = Array::size(mai_selector, this);

    if (data_ref)
      utils::free(data_ref, data_size, env_);

    data_ref = NULL;
  }

#if 0
  template<int MAI>
  void array_free() {

    typename Array<MAI>::Data_t& data_ref = Array<MAI>::data(this);
    size_t data_size = Array<MAI>::size(this);

    if (data_ref)
      utils::free(data_ref, data_size, env_);

    data_ref = NULL;
  }
#endif

private:

  CEnv& env_;
  bool is_allocated_;

  void* __restrict__ data_;
  void* __restrict__ data_S_;
  void* __restrict__ data_C_;
  MetricItemCoords_t* __restrict__ coords_;

  size_t data_size_;
  size_t data_S_size_;
  size_t data_C_size_;
  size_t coords_size_;

  friend GMMetrics;

  //---Disallowed methods.

  MetricsMem(  const MetricsMem&);
  void operator=(const MetricsMem&);
};

//=============================================================================
/// \brief Metrics struct declaration.

struct GMMetrics {

  enum {NUM_SECTION_MAX = 6};

  //-----
  // Sizes.
  //-----

  int num_field_local;
  int num_field;
  size_t num_field_active;
  int num_vector_local;
  NV_t num_vector;
  NV_t num_vector_active;
  int J_lo_part3_[NUM_SECTION_MAX];
  int J_wi_part3_[NUM_SECTION_MAX];

  // Number of metrics to compute for this compute_metrics call, this rank.
  NML_t num_metrics_local;

  // Same as num_metrics_local except excludes padding in number of vectors.
  NML_t num_metrics_active_local;

  // Same as num_metrics_local except measured in units of "metric item".
  NML_t num_metric_items_local;

  // Amout of storage space allocated to store metric items.
  NML_t num_metric_items_local_allocated;

//  // Number of metrics items currently buffered in memory.
//  NML_t num_metric_items_local_buffered;

  // Running total of number computed, this compute_metrics call, this rank.
  NML_t num_metric_items_local_computed;

  // Helper values.

  double recip_m;
  NML_t index_offset_part2_;
  NML_t index_offset_part3_;
  NML_t index_offset_section_part1_[NUM_SECTION_MAX];
  NML_t index_offset_section_part2_[NUM_SECTION_MAX];
  bool is_section_num_valid_part1_[NUM_SECTION_MAX];
  bool is_section_num_valid_part2_[NUM_SECTION_MAX];
  size_t section_size_part2_[NUM_SECTION_MAX];
  int phase_block_start_part2_[NUM_SECTION_MAX];
  int phase_block_start_part3_;
  int block_min_part2_;
  int num_steps_2way;

  // Data arrays.

  size_t data_elt_size;
  size_t data_S_elt_size;
  size_t data_C_elt_size;
  void* __restrict__ data;
  void* __restrict__ data_S;
  void* __restrict__ data_C;
  // Map of (contig) index to linearized Cartesian coords.
  MetricItemCoords_t* __restrict__ coords_;
  
  // Accessors.

  MetricItemCoords_t coords_value(NML_t index) const {
    COMET_ASSERT(index+1 >= 1 && index < num_metric_items_local_allocated);
    return coords_[index];
  }

  // Counters.

  void num_metric_items_local_computed_inc(size_t n) {
    num_metric_items_local_computed += n;
  }

  double shrink_achieved_local() const {
    double fuzz = 1;
    return num_metric_items_local == num_metric_items_local_computed ? 1 :
      num_metric_items_local / (num_metric_items_local_computed + fuzz);
  }

  // Other.
  int data_type_id;
  GMDecompMgr* dm;

  bool is_computing_histograms() const {
    return dm->histograms()->is_computing_histograms();
  }
};

//=============================================================================
/// \brief Class for manipulating MetricItemCoords_t values.

struct CoordsInfo {

// TODO: evaluate whether need to make these more performant.

  //---------- Get i/j/k G

  static NV_t getiG(MetricItemCoords_t coords, GMMetrics& metrics,
    CEnv& env) {
    const auto nv_ = static_cast<MetricItemCoords_t>(metrics.num_vector);
    if (env.coords_type_by_metric()) {
      const auto result = safe_cast_assert<NV_t>(
        coords % nv_);
      COMET_ASSERT(result+1 >= 1 && result < metrics.num_vector);
      return result;
    } else {
      const auto result = safe_cast_assert<NV_t>(
        (coords / 2) % nv_);
      COMET_ASSERT(result+1 >= 1 && result < metrics.num_vector);
      return result;
    }
  }

  static NV_t getjG(MetricItemCoords_t coords, GMMetrics& metrics,
    CEnv& env) {
    const auto nv_ = static_cast<MetricItemCoords_t>(metrics.num_vector);
    if (env.coords_type_by_metric()) {
      const auto result = safe_cast_assert<NV_t>(
        (coords / nv_) % nv_);
      COMET_ASSERT(result+1 >= 1 && result < metrics.num_vector);
      return result;
    } else {
      const auto result = safe_cast_assert<NV_t>(
        (coords / (2 * nv_ * 2)) % nv_);
      COMET_ASSERT(result+1 >= 1 && result < metrics.num_vector);
      return result;
    }
  }

  static NV_t getkG(MetricItemCoords_t coords, GMMetrics& metrics,
    CEnv& env) {
    COMET_ASSERT(env.num_way() >= NumWay::_3);
    const auto nv_ = static_cast<MetricItemCoords_t>(metrics.num_vector);
    if (env.coords_type_by_metric()) {
      const auto result = safe_cast_assert<NV_t>(
        coords / (nv_ * nv_));
      COMET_ASSERT(result+1 >= 1 && result < metrics.num_vector);
      return result;
    } else {
      const auto result = safe_cast_assert<NV_t>(
        coords / (2 * nv_ * 2 * nv_ * 2));
      COMET_ASSERT(result+1 >= 1 && result < metrics.num_vector);
      return result;
    }
  }

  static NV_t getG(MetricItemCoords_t coords, int ijk, GMMetrics& metrics,
    CEnv& env) {
    COMET_ASSERT(ijk >= 0 && ijk < env.num_way());
    return ijk==0 ? getiG(coords, metrics, env) :
           ijk==1 ? getjG(coords, metrics, env) :
                    getkG(coords, metrics, env);
  }

  static bool is_active(MetricItemCoords_t coords, GMMetrics& metrics,
    CEnv& env) {
    return getiG(coords, metrics, env) < metrics.num_vector_active &&
           getjG(coords, metrics, env) < metrics.num_vector_active &&
          (env.num_way() == NumWay::_2 ||
           getkG(coords, metrics, env) < metrics.num_vector_active);
  }

  //---------- Get i/j/k E

  static int getiE(MetricItemCoords_t coords, int entry_num, GMMetrics& metrics,
    CEnv& env) {
    COMET_ASSERT(entry_num >= 0 &&
        entry_num < env.coords_type_by_metric() ? env.pow2_num_way() : 1);
    if (env.coords_type_by_metric()) {
      const int result = entry_num / (env.pow2_num_way()/2);
      COMET_ASSERT(result >= 0 && result < 2);
      return result;
    } else {
      const int result = safe_cast_assert<int>(coords % 2);
      COMET_ASSERT(result >= 0 && result < 2);
      return result;
    }
  }

  static int getjE(MetricItemCoords_t coords, int entry_num, GMMetrics& metrics,
    CEnv& env) {
    COMET_ASSERT(entry_num >= 0 &&
        entry_num < env.coords_type_by_metric() ? env.pow2_num_way() : 1);
    if (env.coords_type_by_metric()) {
      const int result = (entry_num / (env.pow2_num_way()/4)) % 2;
      COMET_ASSERT(result >= 0 && result < 2);
      return result;
    } else {
      const auto nv_ = static_cast<MetricItemCoords_t>(metrics.num_vector);
      const int result = safe_cast_assert<int>((coords / (2 * nv_)) % 2);
      COMET_ASSERT(result >= 0 && result < 2);
      return result;
    }
  }

  static int getkE(MetricItemCoords_t coords, int entry_num, GMMetrics& metrics,
    CEnv& env) {
    COMET_ASSERT(entry_num >= 0 &&
        entry_num < env.coords_type_by_metric() ? env.pow2_num_way() : 1);
    COMET_ASSERT(env.num_way() >= NumWay::_3);
    if (env.coords_type_by_metric()) {
      const int result = entry_num % 2;
      COMET_ASSERT(result >= 0 && result < 2);
      return result;
    } else {
      const auto nv_ = static_cast<MetricItemCoords_t>(metrics.num_vector);
      const int result = safe_cast_assert<int>((coords / (2 * nv_ * 2 * nv_)) % 2);
      COMET_ASSERT(result >= 0 && result < 2);
      return result;
    }
  }

  static int getE(MetricItemCoords_t coords, int ijk, int entry_num,
    GMMetrics& metrics, CEnv& env) {
    COMET_ASSERT(ijk >= 0 && ijk < env.num_way());
    return ijk==0 ? getiE(coords, entry_num, metrics, env) :
           ijk==1 ? getjE(coords, entry_num, metrics, env) :
                    getkE(coords, entry_num, metrics, env);
  }

  //---------- Set

  static MetricItemCoords_t set(NV_t iG, NV_t jG, GMMetrics& metrics,
    CEnv& env) {
    COMET_ASSERT(env.coords_type() == CoordsType::BY_METRIC);
    COMET_ASSERT(env.num_way() == NumWay::_2);
    COMET_ASSERT(iG+1 >= 1 && iG < metrics.num_vector);
    COMET_ASSERT(jG+1 >= 1 && jG < metrics.num_vector);
    const auto nv_ = static_cast<MetricItemCoords_t>(metrics.num_vector);
    const auto result = safe_cast_assert<MetricItemCoords_t>(iG + nv_ * (
                                                      jG));
    COMET_ASSERT(getiG(result, metrics, env) == iG);
    COMET_ASSERT(getjG(result, metrics, env) == jG);
    return result;
  }

  static MetricItemCoords_t set(NV_t iG, NV_t jG,
    int iE, int jE, GMMetrics& metrics, CEnv& env) {
    COMET_ASSERT(env.coords_type() == CoordsType::BY_ENTRY);
    COMET_ASSERT(env.num_way() == NumWay::_2);
    COMET_ASSERT(iG+1 >= 1 && iG < metrics.num_vector);
    COMET_ASSERT(jG+1 >= 1 && jG < metrics.num_vector);
    COMET_ASSERT(iE >= 0 && iE < 2);
    COMET_ASSERT(jE >= 0 && jE < 2);
    const auto nv_ = static_cast<MetricItemCoords_t>(metrics.num_vector);
    const auto result = safe_cast_assert<MetricItemCoords_t>(iE + 2 * (
                                                      iG + nv_ * (
                                                      jE + 2 * (
                                                      jG))));
    COMET_ASSERT(getiG(result, metrics, env) == iG);
    COMET_ASSERT(getjG(result, metrics, env) == jG);
    COMET_ASSERT(getiE(result, 0, metrics, env) == iE);
    COMET_ASSERT(getjE(result, 0, metrics, env) == jE);
    return result;
  }

  static MetricItemCoords_t set(NV_t iG, NV_t jG, NV_t kG,
    GMMetrics& metrics, CEnv& env) {
    COMET_ASSERT(env.coords_type() == CoordsType::BY_METRIC);
    COMET_ASSERT(env.num_way() == NumWay::_3);
    COMET_ASSERT(iG+1 >= 1 && iG < metrics.num_vector);
    COMET_ASSERT(jG+1 >= 1 && jG < metrics.num_vector);
    COMET_ASSERT(kG+1 >= 1 && kG < metrics.num_vector);
    const auto nv_ = static_cast<MetricItemCoords_t>(metrics.num_vector);
    const auto result = safe_cast_assert<MetricItemCoords_t>(iG + nv_ * (
                                                      jG + nv_ * (
                                                      kG)));
    COMET_ASSERT(getiG(result, metrics, env) == iG);
    COMET_ASSERT(getjG(result, metrics, env) == jG);
    COMET_ASSERT(getkG(result, metrics, env) == kG);
    return result;
  }

  static MetricItemCoords_t set(NV_t iG, NV_t jG, NV_t kG,
    int iE, int jE, int kE, GMMetrics& metrics, CEnv& env) {
    COMET_ASSERT(env.coords_type() == CoordsType::BY_ENTRY);
    COMET_ASSERT(env.num_way() == NumWay::_3);
    COMET_ASSERT(iG+1 >= 1 && iG < metrics.num_vector);
    COMET_ASSERT(jG+1 >= 1 && jG < metrics.num_vector);
    COMET_ASSERT(kG+1 >= 1 && kG < metrics.num_vector);
    COMET_ASSERT(iE >= 0 && iE < 2);
    COMET_ASSERT(jE >= 0 && jE < 2);
    COMET_ASSERT(kE >= 0 && kE < 2);
    const auto nv_ = static_cast<MetricItemCoords_t>(metrics.num_vector);
    const auto result = safe_cast_assert<MetricItemCoords_t>(iE + 2 * (
                                                      iG + nv_ * (
                                                      jE + 2 * (
                                                      jG + nv_ * (
                                                      kE + 2 * (
                                                      kG))))));
    COMET_ASSERT(getiG(result, metrics, env) == iG);
    COMET_ASSERT(getjG(result, metrics, env) == jG);
    COMET_ASSERT(getkG(result, metrics, env) == kG);
    COMET_ASSERT(getiE(result, 0, metrics, env) == iE);
    COMET_ASSERT(getjE(result, 0, metrics, env) == jE);
    COMET_ASSERT(getkE(result, 0, metrics, env) == kE);
    return result;
  }
};

//=============================================================================
// Null object.

GMMetrics GMMetrics_null(void);

//=============================================================================
// Metrics pseudo-constructor.

void GMMetrics_create(GMMetrics* metrics, int data_type_id,
                      GMDecompMgr* dm, MetricsMem* metrics_mem, CEnv* env);

//-----------------------------------------------------------------------------

void GMMetrics_set_num_metrics(GMMetrics& metrics, int nvl, const CEnv& env);

//=============================================================================
// Metrics pseudo-destructor.

void GMMetrics_destroy(GMMetrics* metrics, CEnv* env);

//=============================================================================
// Accessors: indexing: global coord from (contig) index: generic.

static NV_t Metrics_coords_getG(GMMetrics& metrics, NML_t index, int ijk,
  CEnv& env) {
  COMET_ASSERT(index+1 >= 1 &&
               index < metrics.num_metric_items_local_allocated);
  COMET_ASSERT(ijk >= 0 && ijk < env.num_way());

  const NV_t result = CoordsInfo::getG(metrics.coords_value(index), ijk,
    metrics, env);

  return result;
}

//=============================================================================
// Adjustment required to compensate for padding.

void gm_metrics_pad_adjust(GMMetrics* metrics,
                           MirroredBuf* metrics_buf,
                           CEnv* env,
                           int weight = 1);

//=============================================================================
/// \brief  Is this (section_)block_num to be processed by current proc_repl.

//static bool gm_proc_r_active(int section_block_num, const CEnv* const env) {
//  COMET_ASSERT(env);
//  COMET_ASSERT(section_block_num >= 0);
//  return section_block_num % env->num_proc_repl() == env->proc_num_repl();
//}

static bool metrics_is_proc_repl_active(const GMMetrics& metrics,
  int section_block_num, const CEnv& env) {
  COMET_ASSERT(section_block_num >= 0);
  COMET_ASSERT((!(env.is_comm_ring() && env.num_way() != NumWay::_2)) &&
               "Unimplemented.");
  // NOTE section_block_num is based on a numbering within current phase.

  return env.is_comm_ring() ?
    // SRP ordering of axes for block diags.
    section_block_num / metrics.num_steps_2way == env.proc_num_repl() :
    // RSP ordering of axes for block diags.
    section_block_num % env.num_proc_repl() == env.proc_num_repl();
}

//=============================================================================

//-----------------------------------------------------------------------------
/*!
 * \class MetricsArrayData
 * \brief Helper class to access selected metrics array.
 *
 */
//-----------------------------------------------------------------------------
template<int MA> struct MetricsArrayData;

// Metrics proper.
template<> struct MetricsArrayData<MetricsArrayId::M> {
  static size_t elt_size(const GMMetrics& metrics) {
   return metrics.data_elt_size;
  }
  static void* __restrict__ p(const GMMetrics& metrics) {return metrics.data;}
};

// Sums.
template<> struct MetricsArrayData<MetricsArrayId::S> {
  static size_t elt_size(const GMMetrics& metrics) {
    return metrics.data_S_elt_size;
  }
  static void* __restrict__ p(const GMMetrics& metrics) {return metrics.data_S;}
};

// Counts.
template<> struct MetricsArrayData<MetricsArrayId::C> {
  static size_t elt_size(const GMMetrics& metrics) {
    return metrics.data_C_elt_size;
  }
  static void* __restrict__ p(const GMMetrics& metrics) {return metrics.data_C;}
};

//-----------------------------------------------------------------------------
/*!
 * \brief Const accessor for direct access to metrics array.
 *
 */
template<typename T, int MA = MetricsArrayId::M>
static T Metrics_elt_const(const GMMetrics& metrics, NML_t index, CEnv& env) {
  COMET_ASSERT(sizeof(T) == MetricsArrayData<MA>::elt_size(metrics));
  COMET_ASSERT(MetricsArrayData<MA>::p(metrics));
  COMET_ASSERT(index+1>=1 && index < metrics.num_metric_items_local_allocated);
  return reinterpret_cast<T*>(MetricsArrayData<MA>::p(metrics))[index];
}

//-----------------------------------------------------------------------------
/*!
 * \brief Accessor for direct access to metrics array.
 *
 */
template<typename T, int MA = MetricsArrayId::M>
static T& Metrics_elt(GMMetrics& metrics, NML_t index, CEnv& env) {
  COMET_ASSERT(sizeof(T) == MetricsArrayData<MA>::elt_size(metrics));
  COMET_ASSERT(MetricsArrayData<MA>::p(metrics));
  COMET_ASSERT(index+1>=1 && index < metrics.num_metric_items_local_allocated);
  return reinterpret_cast<T*>(MetricsArrayData<MA>::p(metrics))[index];
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
// Implementation include files.

#include "metrics_2way_indexing.i.hh"
#include "metrics_2way_accessors.i.hh"
#include "metrics_3way_indexing.i.hh"
#include "metrics_3way_accessors.i.hh"

#endif // _COMET_METRICS_HH_

//-----------------------------------------------------------------------------
