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

#include <stdio.h> /*FIX*/

#include <stddef.h>

#include "env.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Struct declaration---*/

typedef struct {
  int num_vector;
  int num_vector_local;
  size_t num_elts_local;
  int data_type_id;
  /*---map (contig) index to linearized Cartesian coords---*/
  size_t* coords_global_from_index;
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
/*---Accessors: indexing: (contig) index from coord---*/

static size_t GMMetrics_index_from_coord_2(GMMetrics* metrics,
                                           int i,
                                           int j,
                                           GMEnv* env) {
  GMAssert(metrics);
  GMAssert(!env->all2all);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(env);

  size_t index = ((j * (size_t)(j - 1)) >> 1) + i;
  GMAssert(i + metrics->num_vector_local * env->proc_num ==
           metrics->coords_global_from_index[index] % metrics->num_vector);
  GMAssert(j + metrics->num_vector_local * env->proc_num ==
           metrics->coords_global_from_index[index] / metrics->num_vector);
  return index;
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_index_from_coord_all2all_2(GMMetrics* metrics,
                                                   int i,
                                                   int j,
                                                   int j_proc,
                                                   GMEnv* env) {
  GMAssert(metrics);
  GMAssert(env->all2all);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_proc >= 0);
  GMAssert(j_proc < env->num_proc);
  GMAssert(i < j || j_proc != env->proc_num);
  /*---WARNING: these conditions on j_proc are not exhaustive---*/
  GMAssert(env);

  const int nvl = metrics->num_vector_local;

  size_t index = j_proc == env->proc_num
                     ? ((j * (size_t)(j - 1)) >> 1) + i
                     : ((nvl * (size_t)(nvl - 1)) >> 1) +
                           (((j_proc - env->proc_num - 1 + 2 * env->num_proc) %
                             env->num_proc) *
                                nvl +
                            j) *
                               nvl +
                           i;
  GMAssert(index >= 0 && index < metrics->num_elts_local);
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
  GMAssert(i + metrics->num_vector_local * env->proc_num ==
           metrics->coords_global_from_index[index] % metrics->num_vector);
  GMAssert(j + metrics->num_vector_local * env->proc_num ==
           (metrics->coords_global_from_index[index] / metrics->num_vector) %
               metrics->num_vector);
  GMAssert(k + metrics->num_vector_local * env->proc_num ==
           metrics->coords_global_from_index[index] /
               (metrics->num_vector * metrics->num_vector));
  return index;
}

/*===========================================================================*/
/*---Accessors: value from (contig) index---*/

static GMFloat GMMetrics_float_get_from_index(GMMetrics* metrics,
                                              int index,
                                              GMEnv* env) {
  GMAssert(metrics);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env);

  return ((GMFloat*)(metrics->data))[index];
}

/*===========================================================================*/
/*---Accessors: value from coord---*/

static void GMMetrics_float_set_2(GMMetrics* metrics,
                                  int i,
                                  int j,
                                  GMFloat value,
                                  GMEnv* env) {
  GMAssert(metrics);
  GMAssert(!env->all2all);
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

static void GMMetrics_float_set_all2all_2(GMMetrics* metrics,
                                          int i,
                                          int j,
                                          int j_proc,
                                          GMFloat value,
                                          GMEnv* env) {
  GMAssert(metrics);
  GMAssert(env->all2all);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_proc >= 0);
  GMAssert(j_proc < env->num_proc);
  GMAssert(i < j || j_proc != env->proc_num);
  /*---WARNING: these conditions on j_proc are not exhaustive---*/
  GMAssert(env);

  size_t index =
      GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_proc, env);
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

static GMFloat GMMetrics_float_get_2(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     GMEnv* env) {
  GMAssert(metrics);
  GMAssert(!env->all2all);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(env);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
}

/*---------------------------------------------------------------------------*/

static GMFloat GMMetrics_float_get_all2all_2(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int j_proc,
                                             GMEnv* env) {
  GMAssert(metrics);
  GMAssert(env->all2all);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_proc >= 0);
  GMAssert(j_proc < env->num_proc);
  GMAssert(i < j || j_proc != env->proc_num);
  /*---WARNING: these conditions on j_proc are not exhaustive---*/
  GMAssert(env);

  size_t index =
      GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_proc, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
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
  return GMMetrics_float_get_from_index(metrics, index, env);
}

/*===========================================================================*/
/*---Accessors: indexing: global coord from (contig) index---*/

static int GMMetrics_coord0_global_from_index_2(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env);
  GMAssert(env->num_way == 2);

  const int i = metrics->coords_global_from_index[index] % metrics->num_vector;
  return i;
}

/*---------------------------------------------------------------------------*/

static int GMMetrics_coord1_global_from_index_2(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env);
  GMAssert(env->num_way == 2);

  const int j = metrics->coords_global_from_index[index] / metrics->num_vector;
  return j;
}

/*---------------------------------------------------------------------------*/

static int GMMetrics_coord0_global_from_index_3(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env);
  GMAssert(env->num_way == 3);

  const int i = metrics->coords_global_from_index[index] % metrics->num_vector;
  return i;
}

/*---------------------------------------------------------------------------*/

static int GMMetrics_coord1_global_from_index_3(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env);
  GMAssert(env->num_way == 3);

  const int j =
      (metrics->coords_global_from_index[index] / metrics->num_vector) %
      metrics->num_vector;
  return j;
}

/*---------------------------------------------------------------------------*/

static int GMMetrics_coord2_global_from_index_3(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env);
  GMAssert(env->num_way == 3);

  const int k = metrics->coords_global_from_index[index] /
                (metrics->num_vector * metrics->num_vector);
  return k;
}

/*---------------------------------------------------------------------------*/
static int GMMetrics_coord_global_from_index(GMMetrics* metrics,
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

  GMAssert(env->num_way < 4
               ? "num_way > 3 not currently functional; please modify code"
               : 0);

  switch (env->num_way + 4 * coord_num) {
    case 2 + 4 * 0: /* 2-way, coord 0 */
      result = GMMetrics_coord0_global_from_index_2(metrics, index, env);
      break;
    case 2 + 4 * 1: /* 2-way, coord 1 */
      result = GMMetrics_coord1_global_from_index_2(metrics, index, env);
      break;
    case 3 + 4 * 0: /* 3-way, coord 0 */
      result = GMMetrics_coord0_global_from_index_3(metrics, index, env);
      break;
    case 3 + 4 * 1: /* 3-way, coord 1 */
      result = GMMetrics_coord1_global_from_index_3(metrics, index, env);
      break;
    case 3 + 4 * 2: /* 3-way, coord 2 */
      result = GMMetrics_coord2_global_from_index_3(metrics, index, env);
      break;
    default:
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/

  return result;
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_metrics_h_---*/

/*---------------------------------------------------------------------------*/
