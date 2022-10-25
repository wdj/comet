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
/// \brief Helper class for metrics memory.

// Forward declartion.
struct GMMetrics;

class MetricsMem {

public:
  MetricsMem(CEnv* env);
  ~MetricsMem();
  void deallocate();

  void* malloc_data(size_t data_size);
  void* malloc_data_S(size_t data_size_S);
  void* malloc_data_C(size_t data_size_C);
  MetricItemCoords_t* malloc_coords(size_t coords_size);

private:

  CEnv* env_;
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

  // Logical sizes.
  int num_field_local;
  int num_field;
  size_t num_field_active;
  int num_vector_local;
  NV_t num_vector;
  NV_t num_vector_active;
  int J_lo_part3_[NUM_SECTION_MAX];
  int J_wi_part3_[NUM_SECTION_MAX];
  NML_t num_metrics_local;
  NML_t num_metric_items_local;
  NML_t num_metric_items_local_allocated;
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
  void* __restrict__ data;
  void* __restrict__ data_S;
  void* __restrict__ data_C;
  size_t data_elt_size;
  size_t data_S_elt_size;
  size_t data_C_elt_size;
  // Map of (contig) index to linearized Cartesian coords.
  MetricItemCoords_t* __restrict__ coords_;
  // Accessor.
  MetricItemCoords_t coords_value(NML_t index) const {
    COMET_ASSERT(index+1 >= 1 && index < num_metric_items_local_allocated);
    return coords_[index];
  }
  // Counters.
  NML_t num_metric_items_local_computed;
  void num_metric_items_local_computed_inc(size_t n) {
    num_metric_items_local_computed += n;
  }
  double shrink_achieved_local() const {
    double fuzz = 1;
    return num_metric_items_local == num_metric_items_local_computed ? 1 :
      num_metric_items_local / (num_metric_items_local_computed + fuzz);
  }
  NML_t num_metrics_active_local;
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
      const auto result = safe_cast<NV_t>(
        coords % nv_);
      COMET_ASSERT(result+1 >= 1 && result < metrics.num_vector);
      return result;
    } else {
      const auto result = safe_cast<NV_t>(
        (coords / 2) % nv_);
      COMET_ASSERT(result+1 >= 1 && result < metrics.num_vector);
      return result;
    }
  }

  static NV_t getjG(MetricItemCoords_t coords, GMMetrics& metrics,
    CEnv& env) {
    const auto nv_ = static_cast<MetricItemCoords_t>(metrics.num_vector);
    if (env.coords_type_by_metric()) {
      const auto result = safe_cast<NV_t>(
        (coords / nv_) % nv_);
      COMET_ASSERT(result+1 >= 1 && result < metrics.num_vector);
      return result;
    } else {
      const auto result = safe_cast<NV_t>(
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
      const auto result = safe_cast<NV_t>(
        coords / (nv_ * nv_));
      COMET_ASSERT(result+1 >= 1 && result < metrics.num_vector);
      return result;
    } else {
      const auto result = safe_cast<NV_t>(
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
      const int result = safe_cast<int>(coords % 2);
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
      const int result = safe_cast<int>((coords / (2 * nv_)) % 2);
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
      const int result = safe_cast<int>((coords / (2 * nv_ * 2 * nv_)) % 2);
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
    const auto result = safe_cast<MetricItemCoords_t>(iG + nv_ * (
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
    const auto result = safe_cast<MetricItemCoords_t>(iE + 2 * (
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
    const auto result = safe_cast<MetricItemCoords_t>(iG + nv_ * (
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
    const auto result = safe_cast<MetricItemCoords_t>(iE + 2 * (
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

struct MetricsArray {
  enum {_ = 0,
        S = 1,
        C = 2};
};

template<int MA> struct MetricsArrayData;

template<> struct MetricsArrayData<MetricsArray::_> {
  static size_t elt_size(const GMMetrics& metrics) {
   return metrics.data_elt_size;
  }
  static void* __restrict__ p(const GMMetrics& metrics) {return metrics.data;}
};

template<> struct MetricsArrayData<MetricsArray::S> {
  static size_t elt_size(const GMMetrics& metrics) {
    return metrics.data_S_elt_size;
  }
  static void* __restrict__ p(const GMMetrics& metrics) {return metrics.data_S;}
};

template<> struct MetricsArrayData<MetricsArray::C> {
  static size_t elt_size(const GMMetrics& metrics) {
    return metrics.data_C_elt_size;
  }
  static void* __restrict__ p(const GMMetrics& metrics) {return metrics.data_C;}
};

template<typename T, int MA = MetricsArray::_>
static T Metrics_elt_const(const GMMetrics& metrics, NML_t index, CEnv& env) {
  COMET_ASSERT(sizeof(T) == MetricsArrayData<MA>::elt_size(metrics));
  COMET_ASSERT(MetricsArrayData<MA>::p(metrics));
  //COMET_ASSERT(index+1 >= 1 && index < metrics.num_metrics_local);
  COMET_ASSERT(index+1>=1 && index < metrics.num_metric_items_local_allocated);
  return ((T*)MetricsArrayData<MA>::p(metrics))[index];
}

template<typename T, int MA = MetricsArray::_>
static T& Metrics_elt(GMMetrics& metrics, NML_t index, CEnv& env) {
  COMET_ASSERT(sizeof(T) == MetricsArrayData<MA>::elt_size(metrics));
  COMET_ASSERT(MetricsArrayData<MA>::p(metrics));
  //COMET_ASSERT(index+1 >= 1 && index < metrics.num_metrics_local);
  COMET_ASSERT(index+1>=1 && index < metrics.num_metric_items_local_allocated);
  return ((T*)MetricsArrayData<MA>::p(metrics))[index];
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
