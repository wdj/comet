/*---------------------------------------------------------------------------*/
/*!
 * \file   metrics.h
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

/*=============================================================================


=============================================================================*/

#ifndef _metrics_h_
#define _metrics_h_

#include <stddef.h>

#include "env.h"

/*===========================================================================*/
/*---Struct declaration---*/

typedef struct {
  int num_vector;
  int num_vector_local;
  int num_vector_local_max;
  size_t num_elts_local;
  int data_type_id;
  size_t* index_map; /*---map contig index to linearized Cartesian index---*/
  void* __restrict__ data;
} GMMetrics;

/*===========================================================================*/
/*---Null object---*/

GMMetrics GMMetrics_null(void);

/*===========================================================================*/
/*---Metrics pseudo-constructor---*/

void GMMetrics_create(GMMetrics* metrics,
                    int data_type_id,
                    int num_vector_local,
                    GMEnv* env);

/*===========================================================================*/
/*---Metrics pseudo-destructor---*/

void GMMetrics_destroy(GMMetrics* metrics, GMEnv* env);

/*===========================================================================*/
/*---Metrics checksum---*/

double GMMetrics_checksum(GMMetrics* metrics, GMEnv* env);

/*===========================================================================*/
/*---Accessors: indexing: contig index from coord---*/

static size_t GMMetrics_index_from_coord_2(GMMetrics* metrics,
                                         int i,
                                         int j,
                                         GMEnv* env) {
  GMAssert(metrics);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(env);

  size_t index = ((j * (size_t)(j - 1)) >> 1) + i;
  return index;
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_index_from_coord_3(GMMetrics* metrics,
                                         int i,
                                         int j,
                                         int k,
                                         GMEnv* env) {
  GMAssert(metrics);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(k >= 0);
  GMAssert(k < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(j < k);
  GMAssert(env);

  size_t index = (k * (size_t)(k - 1) * (size_t)(k - 2)) / 6 +
                 (j * (size_t)(j - 1)) / 2 + i;
  return index;
}

/*===========================================================================*/
/*---Accessors: value from coord---*/

static void GMMetrics_float_set_2(GMMetrics* metrics,
                                int i,
                                int j,
                                GMFloat value,
                                GMEnv* env) {
  GMAssert(metrics);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(env);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_float_set_3(GMMetrics* metrics,
                                int i,
                                int j,
                                int k,
                                GMFloat value,
                                GMEnv* env) {
  GMAssert(metrics);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(k >= 0);
  GMAssert(k < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(j < k);
  GMAssert(env);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

/*---------------------------------------------------------------------------*/

static GMFloat GMMetrics_float_get_from_index(GMMetrics* metrics, int index,
                                            GMEnv* env) {
  GMAssert(metrics);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env);

  return ((GMFloat*)(metrics->data))[index];
}

/*---------------------------------------------------------------------------*/

static GMFloat GMMetrics_float_get_2(GMMetrics* metrics, int i, int j, GMEnv* env) {
  GMAssert(metrics);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(env);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  return GMMetrics_float_get_from_index( metrics, index, env );
}

/*---------------------------------------------------------------------------*/

static GMFloat GMMetrics_float_get_3(GMMetrics* metrics,
                                   int i,
                                   int j,
                                   int k,
                                   GMEnv* env) {
  GMAssert(metrics);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(k >= 0);
  GMAssert(k < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(j < k);
  GMAssert(env);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  return GMMetrics_float_get_from_index( metrics, index, env );
}

/*===========================================================================*/
/*---Accessors: indexing: coord from contig index---*/

static int GMMetrics_coord0_from_index_2(GMMetrics* metrics,
                                       size_t index,
                                       GMEnv* env) {
  GMAssert(metrics);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env);
  GMAssert(env->num_way == 2);

  int i = metrics->index_map[index] % metrics->num_vector_local;
  return i;
}

/*---------------------------------------------------------------------------*/

static int GMMetrics_coord1_from_index_2(GMMetrics* metrics,
                                       size_t index,
                                       GMEnv* env) {
  GMAssert(metrics);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env);
  GMAssert(env->num_way == 2);

  int j = metrics->index_map[index] / metrics->num_vector_local;
  return j;
}

/*---------------------------------------------------------------------------*/

static int GMMetrics_coord0_from_index_3(GMMetrics* metrics,
                                       size_t index,
                                       GMEnv* env) {
  GMAssert(metrics);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env);
  GMAssert(env->num_way == 3);

  int i = metrics->index_map[index] % metrics->num_vector_local;
  return i;
}

/*---------------------------------------------------------------------------*/

static int GMMetrics_coord1_from_index_3(GMMetrics* metrics,
                                       size_t index,
                                       GMEnv* env) {
  GMAssert(metrics);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env);
  GMAssert(env->num_way == 3);

  int i = (metrics->index_map[index] / metrics->num_vector_local) %
          metrics->num_vector_local;
  return i;
}

/*---------------------------------------------------------------------------*/

static int GMMetrics_coord2_from_index_3(GMMetrics* metrics,
                                       size_t index,
                                       GMEnv* env) {
  GMAssert(metrics);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env);
  GMAssert(env->num_way == 3);

  int i = metrics->index_map[index] /
          (metrics->num_vector_local * metrics->num_vector_local);
  return i;
}

/*---------------------------------------------------------------------------*/
static int GMMetrics_coord_from_index(GMMetrics* metrics,
                                    size_t index,
                                    int coord_num,
                                    GMEnv* env) {
  GMAssert(metrics);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env);
  GMAssert(coord_num >= 0);
  GMAssert(coord_num < env->num_way);

  int result = 0;

  switch (env->num_way + 4 * coord_num) {
    case 2 + 4 * 0: /* 2-way, coord 0 */
      result = GMMetrics_coord0_from_index_2(metrics, index, env);
      break;
    case 2 + 4 * 1: /* 2-way, coord 1 */
      result = GMMetrics_coord1_from_index_2(metrics, index, env);
      break;
    case 3 + 4 * 0: /* 3-way, coord 0 */
      result = GMMetrics_coord0_from_index_3(metrics, index, env);
      break;
    case 3 + 4 * 1: /* 3-way, coord 1 */
      result = GMMetrics_coord1_from_index_3(metrics, index, env);
      break;
    case 3 + 4 * 2: /* 3-way, coord 2 */
      result = GMMetrics_coord2_from_index_3(metrics, index, env);
      break;
    default:
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/

  return result;
}

/*===========================================================================*/

#endif /*---_metrics_h_---*/

/*---------------------------------------------------------------------------*/
