//-----------------------------------------------------------------------------
/*!
 * \file   metrics.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_metrics_hh_
#define _gm_metrics_hh_

#include "stddef.h"

#include "env.hh"
#include "decomp_mgr.hh"
#include "mirrored_buf.hh"

//=============================================================================
/*---Helper class for memory---*/

class GMMetricsMem {

public:
  GMMetricsMem(GMEnv* env);
  ~GMMetricsMem();

  void* malloc_data(size_t data_size);
  void* malloc_data_S(size_t data_size_S);
  void* malloc_data_C(size_t data_size_C);
  void* malloc_coords_global_from_index(size_t coords_global_from_index_size);

private:

  GMEnv* env_;
  void* __restrict__ data_;
  size_t data_size_;
  void* __restrict__ data_S_;
  size_t data_S_size_;
  void* __restrict__ data_C_;
  size_t data_C_size_;
  void* coords_global_from_index_;
  size_t coords_global_from_index_size_;

  //---Disallowed methods.

  GMMetricsMem(  const GMMetricsMem&);
  void operator=(const GMMetricsMem&);
};

//=============================================================================
/*---Struct declaration---*/

typedef struct {
  /*---Logical sizes---*/
  int num_field;
  int num_field_local;
  size_t num_field_active;
  int num_vector;
  int num_vector_local;
  size_t num_vector_active;
  int nvl6;
  int J_lo_part3_[6];
  int J_wi_part3_[6];
  int pad1;
  size_t num_elts_local;
  /*---Helper values---*/
  GMInt64 index_offset_0_;
  GMInt64 index_offset_01_;
  GMInt64 index_offset_section_part1_[6];
  GMInt64 index_offset_section_part2_[6];
  bool section_num_valid_part1_[6];
  bool section_num_valid_part2_[6];
  size_t section_size_part2[6];
  size_t phase_block_start_2_[6];
  size_t phase_block_start_3_;
  GMFloat m;
  GMFloat recip_m;
  int block_min;
  /*---map of (contig) index to linearized Cartesian coords---*/
  size_t* coords_global_from_index;
  /*---Other---*/
  int data_type_id;
  int data_type_num_values;
  void* __restrict__ data;
  void* __restrict__ data_S;
  void* __restrict__ data_C;
  size_t data_size;
  size_t data_S_size;
  size_t data_C_size;
  size_t num_elts_local_computed;
  GMDecompMgr* dm;
} GMMetrics;

//=============================================================================
/*---Null object---*/

GMMetrics GMMetrics_null(void);

//=============================================================================
/*---Metrics pseudo-constructor---*/

void GMMetrics_create(GMMetrics* metrics, int data_type_id,
                      GMDecompMgr* dm, GMMetricsMem* metrics_mem, GMEnv* env);

//-----------------------------------------------------------------------------

void GMMetrics_3way_num_elts_local(GMMetrics* metrics, int nvl,
                                   GMEnv* env);

//=============================================================================
/*---Metrics pseudo-destructor---*/

void GMMetrics_destroy(GMMetrics* metrics, GMEnv* env);

//=============================================================================
/*---Accessors: indexing: global coord from (contig) index: generic---*/

int GMMetrics_coord_global_from_index(GMMetrics* metrics,
                                      size_t index,
                                      int coord_num,
                                      GMEnv* env);

//=============================================================================
// Adjustment required to compensate for padding.

void gm_metrics_pad_adjust(GMMetrics* metrics,
                           GMMirroredBuf* metrics_buf,
                           GMEnv* env);

//=============================================================================
/*---Helper: is this (section_)block_num to be processed by this proc_r---*/

static bool gm_proc_r_active(int section_block_num, const GMEnv* const env) {
  GMAssert(env);
  GMAssert(section_block_num >= 0);
  return section_block_num % GMEnv_num_proc_repl(env)
         == GMEnv_proc_num_repl(env);
}

//=============================================================================
/*---Companion include files---*/

#include "metrics_2way.hh"
#include "metrics_3way.hh"

#endif // _gm_metrics_hh_

//-----------------------------------------------------------------------------
