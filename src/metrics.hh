//-----------------------------------------------------------------------------
/*!
 * \file   metrics.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Class to manage the calculated metrics output from the methods.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_metrics_hh_
#define _comet_metrics_hh_

#include "cstddef"
#include "cstdint"
#include "math.h"

#include "env.hh"
#include "decomp_mgr.hh"
#include "mirrored_buf.hh"

//=============================================================================

namespace comet {

typedef size_t Coords_t;

//-----------------------------------------------------------------------------
/// \brief Helper class for metrics memory.

// Forward declartion.
struct GMMetrics;

class MetricsMem {

public:
  MetricsMem(CEnv* env);
  ~MetricsMem();

  void* malloc_data(size_t data_size);
  void* malloc_data_S(size_t data_size_S);
  void* malloc_data_C(size_t data_size_C);
  Coords_t* malloc_coords_values(size_t coords_values_size);

private:

  CEnv* env_;
  void* __restrict__ data_;
  size_t data_size_;
  void* __restrict__ data_S_;
  size_t data_S_size_;
  void* __restrict__ data_C_;
  size_t data_C_size_;
  Coords_t* coords_values_;
  size_t coords_values_size_;

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
  int num_field;
  int num_field_local;
  size_t num_field_active;
  size_t num_vector;
  int num_vector_local;
  size_t num_vector_active;
  int J_lo_part3_[NUM_SECTION_MAX];
  int J_wi_part3_[NUM_SECTION_MAX];
  size_t num_metrics_local;
  // Helper values.
  double recip_m;
  int64_t index_offset_part2_;
  int64_t index_offset_part3_;
  int64_t index_offset_section_part1_[NUM_SECTION_MAX];
  int64_t index_offset_section_part2_[NUM_SECTION_MAX];
  bool is_section_num_valid_part1_[NUM_SECTION_MAX];
  bool is_section_num_valid_part2_[NUM_SECTION_MAX];
  size_t section_size_part2_[NUM_SECTION_MAX];
  int phase_block_start_part2_[NUM_SECTION_MAX];
  int phase_block_start_part3_;
  int block_min_part2_;
  // Data arrays.
  void* __restrict__ data;
  void* __restrict__ data_S;
  void* __restrict__ data_C;
  size_t data_size;
  size_t data_S_size;
  size_t data_C_size;
  size_t data_elt_size;
  size_t data_S_elt_size;
  size_t data_C_elt_size;
  // Map of (contig) index to linearized Cartesian coords.
  Coords_t* coords_values_;
  Coords_t coords_value(size_t index) const {
    COMET_ASSERT(index+1 >= 1 && index < num_metrics_local);
    return coords_values_[index];
  }
  // Counters.
  size_t num_metrics_local_computed;
  // Other.
  int data_type_id;
  int num_entries_per_metric;
  GMDecompMgr* dm;

private:

};

//=============================================================================

struct CoordsInfo {

  static size_t getiG(Coords_t coords, GMMetrics& metrics, CEnv& env) {
    const size_t result = coords % metrics.num_vector;
    return result;
  }

  static size_t getjG(Coords_t coords, GMMetrics& metrics, CEnv& env) {
    const size_t result = (coords / metrics.num_vector) % metrics.num_vector;
    return result;
  }

  static size_t getkG(Coords_t coords, GMMetrics& metrics, CEnv& env) {
    COMET_ASSERT(env.num_way() >= NUM_WAY::_3);
    const size_t result = coords / (metrics.num_vector * metrics.num_vector);
    COMET_ASSERT(result+1 >= 1 && result < metrics.num_vector);
    return result;
  }

  static size_t getG(Coords_t coords, int ijk, GMMetrics& metrics, CEnv& env) {
    COMET_ASSERT(ijk >= 0 && ijk < env.num_way());
    return ijk==0 ? getiG(coords, metrics, env) :
           ijk==1 ? getjG(coords, metrics, env) :
                    getkG(coords, metrics, env);
  }

  static int getiE(Coords_t coords, int entry_num, GMMetrics& metrics,
    CEnv& env) {
    COMET_ASSERT(entry_num >= 0 && entry_num < (1 << env.num_way()));
    // TODO: make this and the sequels more performant.
    const size_t result = entry_num / (1 << (env.num_way()-1));
    COMET_ASSERT(result >= 0 && result < 2);
    return result;
  }

  static int getjE(Coords_t coords, int entry_num, GMMetrics& metrics,
    CEnv& env) {
    COMET_ASSERT(entry_num >= 0 && entry_num < (1 << env.num_way()));
    const size_t result = (entry_num / (1 << (env.num_way()-2))) % 2;
    COMET_ASSERT(result >= 0 && result < 2);
    return result;
  }

  static int getkE(Coords_t coords, int entry_num, GMMetrics& metrics,
    CEnv& env) {
    COMET_ASSERT(entry_num >= 0 && entry_num < (1 << env.num_way()));
    COMET_ASSERT(env.num_way() >= NUM_WAY::_3);
    const size_t result = entry_num % 2;
    COMET_ASSERT(result >= 0 && result < 2);
    return result;
  }

  static int getE(Coords_t coords, int ijk, int entry_num, GMMetrics& metrics,
    CEnv& env) {
    COMET_ASSERT(ijk >= 0 && ijk < env.num_way());
    return ijk==0 ? getiE(coords, entry_num, metrics, env) :
           ijk==1 ? getjE(coords, entry_num, metrics, env) :
                    getkE(coords, entry_num, metrics, env);
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

void GMMetrics_3way_num_metrics_local(GMMetrics* metrics, int nvl, CEnv* env);

//=============================================================================
// Metrics pseudo-destructor.

void GMMetrics_destroy(GMMetrics* metrics, CEnv* env);

//=============================================================================
// Accessors: indexing: global coord from (contig) index: generic.

size_t Metrics_coords_getG(GMMetrics& metrics, size_t index, int ijk,
 CEnv& env);

//=============================================================================
// Adjustment required to compensate for padding.

void gm_metrics_pad_adjust(GMMetrics* metrics,
                           MirroredBuf* metrics_buf,
                           CEnv* env,
                           int weight = 1);

//=============================================================================
// Helper: is this (section_)block_num to be processed by this proc_r.

static bool gm_proc_r_active(int section_block_num, const CEnv* const env) {
  COMET_ASSERT(env);
  COMET_ASSERT(section_block_num >= 0);
  return section_block_num % env->num_proc_repl()
         == env->proc_num_repl();
}

//=============================================================================

struct MetricsArray {
  enum {_ = 0,
        S = 1,
        C = 2};
};

template<int MA> struct MetricsArrayData;

template<> struct MetricsArrayData<MetricsArray::_> {
  static size_t elt_size(const GMMetrics& metrics) {return metrics.data_elt_size;}
  static void* __restrict__ p(const GMMetrics& metrics) {return metrics.data;}
};

template<> struct MetricsArrayData<MetricsArray::S> {
  static size_t elt_size(const GMMetrics& metrics) {return metrics.data_S_elt_size;}
  static void* __restrict__ p(const GMMetrics& metrics) {return metrics.data_S;}
};

template<> struct MetricsArrayData<MetricsArray::C> {
  static size_t elt_size(const GMMetrics& metrics) {return metrics.data_C_elt_size;}
  static void* __restrict__ p(const GMMetrics& metrics) {return metrics.data_C;}
};

template<typename T, int MA = MetricsArray::_>
static T Metrics_elt_const(const GMMetrics& metrics, size_t index, CEnv& env) {
  COMET_ASSERT(sizeof(T) == MetricsArrayData<MA>::elt_size(metrics));
  COMET_ASSERT(MetricsArrayData<MA>::p(metrics));
  COMET_ASSERT(index+1 >= 1 && index < metrics.num_metrics_local);
  return ((T*)MetricsArrayData<MA>::p(metrics))[index];
}

template<typename T, int MA = MetricsArray::_>
static T& Metrics_elt(GMMetrics& metrics, size_t index, CEnv& env) {
  COMET_ASSERT(sizeof(T) == MetricsArrayData<MA>::elt_size(metrics));
  COMET_ASSERT(MetricsArrayData<MA>::p(metrics));
  COMET_ASSERT(index+1 >= 1 && index < metrics.num_metrics_local);
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

#endif // _comet_metrics_hh_

//-----------------------------------------------------------------------------
