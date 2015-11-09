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
  size_t num_elts_0;
  size_t num_elts_01;
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
/*---Helper functions for 3-way all2all case---*/

static int gm_metrics_3way_section_axis(GMMetrics* metrics,
                                        int i_proc,
                                        int j_proc,
                                        int k_proc,
                                        GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert( i_proc >= 0 );
  GMAssert( j_proc >= 0 );
  GMAssert( k_proc >= 0 );
  GMAssert( i_proc < env->num_proc );
  GMAssert( j_proc < env->num_proc );
  GMAssert( k_proc < env->num_proc );
  GMAssert( i_proc != j_proc );
  GMAssert( i_proc != k_proc );
  GMAssert( k_proc != j_proc );

  /*---NOTE: this could possibly be implemented somewhat more efficiently---*/

  return i_proc < j_proc && i_proc < k_proc  ? 0 : /*---i axis---*/
         j_proc < i_proc && j_proc < k_proc  ? 1 : /*---j axis---*/
                                               2;  /*---k axis---*/
}

/*---------------------------------------------------------------------------*/

static int gm_metrics_3way_section_num(GMMetrics* metrics,
                                       int i_proc,
                                       int j_proc,
                                       int k_proc,
                                       GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert( i_proc >= 0 );
  GMAssert( j_proc >= 0 );
  GMAssert( k_proc >= 0 );
  GMAssert( i_proc < env->num_proc );
  GMAssert( j_proc < env->num_proc );
  GMAssert( k_proc < env->num_proc );
  GMAssert( i_proc != j_proc );
  GMAssert( i_proc != k_proc );
  GMAssert( k_proc != j_proc );

  /*---NOTE: this could possibly be implemented somewhat more efficiently---*/

  return i_proc < k_proc && k_proc < j_proc ?    0 :
         i_proc < j_proc && j_proc < k_proc ?    1 :
         j_proc < i_proc && i_proc < k_proc ?    2 :
         j_proc < k_proc && k_proc < i_proc ?    3 :
         k_proc < j_proc && j_proc < i_proc ?    4 :
      /* k_proc < i_proc && i_proc < j_proc ? */ 5;
}

/*===========================================================================*/
/*---Accessors: indexing: (contig) index from coord---*/

static size_t GMMetrics_index_from_coord_2(GMMetrics* metrics,
                                           int i,
                                           int j,
                                           GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(!env->all2all);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);

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
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(env->all2all);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_proc >= 0);
  GMAssert(j_proc < env->num_proc);
  GMAssert(i < j || j_proc != env->proc_num);
  /*---WARNING: these conditions on j_proc are not exhaustive---*/

  const int i_proc = env->proc_num;

  size_t index = j_proc == i_proc
               ? GMMetrics_index_from_coord_2(metrics, i, j, env)
               : metrics->num_elts_0 +
                 i + metrics->num_vector_local * (
                 j + metrics->num_vector_local * (
                 (j_proc - i_proc - 1 + 2*env->num_proc) % env->num_proc ));
  GMAssert(index >= 0 && index < metrics->num_elts_local);
  return index;
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_index_from_coord_3(GMMetrics* metrics,
                                           int i,
                                           int j,
                                           int k,
                                           GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(k >= 0);
  GMAssert(k < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(j < k);

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

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_index_from_coord_all2all_3(GMMetrics* metrics,
                                                   int i,
                                                   int j,
                                                   int k,
                                                   int j_proc,
                                                   int k_proc,
                                                   GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(env->all2all);
  GMAssert(i >= 0);
  GMAssert(j >= 0);
  GMAssert(k >= 0);
  GMAssert(j_proc >= 0);
  GMAssert(j_proc < env->num_proc);
  GMAssert(k_proc >= 0);
  GMAssert(k_proc < env->num_proc);
  /*---WARNING: these conditions are not exhaustive---*/

  const int i_proc = env->proc_num;

  const int nvl = metrics->num_vector_local;

  const int section_axis = gm_metrics_3way_section_axis(
                                metrics, i_proc, j_proc, k_proc, env);
  const int section_num = gm_metrics_3way_section_num(
                                metrics, i_proc, j_proc, k_proc, env);

  size_t index = j_proc == i_proc && k_proc == i_proc
               ? GMMetrics_index_from_coord_3(metrics, i, j, k, env)
               : j_proc == k_proc
               ? metrics->num_elts_0 +
                 i + nvl * (
                 GMMetrics_index_from_coord_2(metrics, j, k, env) + nvl * (
                 i_proc - ( i_proc > j_proc ) ))
               : metrics->num_elts_0 +
                 i + ( section_axis == 0 ? section_num * nvl : 0 ) +
                     ( section_axis == 0 ? nvl / 6 : nvl ) * (
                 j + ( section_axis == 1 ? section_num * nvl : 0 ) +
                     ( section_axis == 1 ? nvl / 6 : nvl ) * (
                 k + ( section_axis == 2 ? section_num * nvl : 0 ) +
                     ( section_axis == 2 ? nvl / 6 : nvl ) * (
                 i_proc - ( i_proc > j_proc ) - ( i_proc > k_proc ) )));
  GMAssert(index >= 0 && index < metrics->num_elts_local);
  return index;
}

/*===========================================================================*/
/*---Accessors: value from (contig) index---*/

static GMFloat GMMetrics_float_get_from_index(GMMetrics* metrics,
                                              int index,
                                              GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);

  return ((GMFloat*)(metrics->data))[index];
}

/*===========================================================================*/
/*---Accessors: value from (local) coord: set---*/

static void GMMetrics_float_set_2(GMMetrics* metrics,
                                  int i,
                                  int j,
                                  GMFloat value,
                                  GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(!env->all2all);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);

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
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(env->all2all);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_proc >= 0);
  GMAssert(j_proc < env->num_proc);
  GMAssert(i < j || j_proc != env->proc_num);
  /*---WARNING: these conditions on j_proc are not exhaustive---*/

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
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(k >= 0);
  GMAssert(k < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(j < k);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_float_set_all2all_3(GMMetrics* metrics,
                                          int i,
                                          int j,
                                          int k,
                                          int j_proc,
                                          int k_proc,
                                          GMFloat value,
                                          GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(env->all2all);
  GMAssert(i >= 0);
  GMAssert(j >= 0);
  GMAssert(k >= 0);
  GMAssert(j_proc >= 0);
  GMAssert(j_proc < env->num_proc);
  GMAssert(k_proc >= 0);
  GMAssert(k_proc < env->num_proc);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index =
    GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_proc, k_proc, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

/*===========================================================================*/
/*---Accessors: value from coord: (local) get---*/

static GMFloat GMMetrics_float_get_2(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(!env->all2all);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
}

/*---------------------------------------------------------------------------*/

static GMFloat GMMetrics_float_get_all2all_2(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int j_proc,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(env->all2all);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_proc >= 0);
  GMAssert(j_proc < env->num_proc);
  GMAssert(i < j || j_proc != env->proc_num);
  /*---WARNING: these conditions on j_proc are not exhaustive---*/

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
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
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

/*---------------------------------------------------------------------------*/

static GMFloat GMMetrics_float_get_all2all_3(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int k,
                                             int j_proc,
                                             int k_proc,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(env->all2all);
  GMAssert(i >= 0);
  GMAssert(j >= 0);
  GMAssert(k >= 0);
  GMAssert(j_proc >= 0);
  GMAssert(j_proc < env->num_proc);
  GMAssert(k_proc >= 0);
  GMAssert(k_proc < env->num_proc);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index =
    GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_proc, k_proc, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
}

/*===========================================================================*/
/*---Accessors: indexing: global coord from (contig) index: 2-way---*/

static int GMMetrics_coord0_global_from_index_2(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env->num_way == 2);

  const int i = metrics->coords_global_from_index[index] % metrics->num_vector;
  return i;
}

/*---------------------------------------------------------------------------*/

static int GMMetrics_coord1_global_from_index_2(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env->num_way == 2);

  const int j = metrics->coords_global_from_index[index] / metrics->num_vector;
  return j;
}

/*===========================================================================*/
/*---Accessors: indexing: global coord from (contig) index: 3-way---*/

static int GMMetrics_coord0_global_from_index_3(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env->num_way == 3);

  const int i = metrics->coords_global_from_index[index] % metrics->num_vector;
  return i;
}

/*---------------------------------------------------------------------------*/

static int GMMetrics_coord1_global_from_index_3(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
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
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(env->num_way == 3);

  const int k = metrics->coords_global_from_index[index] /
                (metrics->num_vector * metrics->num_vector);
  return k;
}

/*===========================================================================*/
/*---Accessors: indexing: global coord from (contig) index: generic---*/

static int GMMetrics_coord_global_from_index(GMMetrics* metrics,
                                             size_t index,
                                             int coord_num,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
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
