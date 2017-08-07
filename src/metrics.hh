/*---------------------------------------------------------------------------*/
/*!
 * \file   metrics.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _gm_metrics_hh_
#define _gm_metrics_hh_

#include <stddef.h>

#include "env.hh"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
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
} GMMetrics;

/*===========================================================================*/
/*---Null object---*/

GMMetrics GMMetrics_null(void);

/*===========================================================================*/
/*---Metrics pseudo-constructor---*/

void GMMetrics_create(GMMetrics* metrics,
                      int data_type_id,
                      int num_field,
                      size_t num_field_active,
                      int num_vector_local,
                      size_t num_vector_active,
                      GMEnv* env);

/*---------------------------------------------------------------------------*/

void GMMetrics_3way_num_elts_local(GMMetrics* metrics, int nvl,
                                   GMEnv* env);

/*===========================================================================*/
/*---Metrics pseudo-destructor---*/

void GMMetrics_destroy(GMMetrics* metrics, GMEnv* env);

/*===========================================================================*/
/*---Accessors: indexing: global coord from (contig) index: generic---*/

int GMMetrics_coord_global_from_index(GMMetrics* metrics,
                                      size_t index,
                                      int coord_num,
                                      GMEnv* env);

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/
/*---Companion include files---*/

#include "metrics_2way.hh"
#include "metrics_3way.hh"

/*===========================================================================*/

#endif /*---_gm_metrics_hh_---*/

/*---------------------------------------------------------------------------*/
