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
  size_t* index_map;  /*---map contig index to linearized Cartesian index---*/
  void* __restrict__ data;
} Metrics;

/*===========================================================================*/
/*---Null object---*/

Metrics Metrics_null(void);

/*===========================================================================*/
/*---Vectors pseudo-constructor---*/

void Metrics_create(Metrics* metrics,
                    int data_type_id,
                    int num_vector_local,
                    Env* env);

/*===========================================================================*/
/*---Vectors pseudo-destructor---*/

void Metrics_destroy(Metrics* metrics, Env* env);

/*===========================================================================*/
/*---Metrics checksum---*/

double Metrics_checksum(Metrics * metrics, Env * env);

/*===========================================================================*/
/*---Accessors: index from coord---*/

static size_t Metrics_index_from_coord_2(Metrics* metrics,
                                         int i,
                                         int j,
                                         Env* env) {
  Assert(metrics);
  Assert(i >= 0);
  Assert(i < metrics->num_vector_local);
  Assert(j >= 0);
  Assert(j < metrics->num_vector_local);
  Assert(env);

  Assert(i != j);

  size_t index = i > j ? ((i * (size_t)(i - 1)) >> 1) + j
                       : ((j * (size_t)(j - 1)) >> 1) + i;
  return index;
}

/*---------------------------------------------------------------------------*/

static size_t Metrics_index_from_coord_3(Metrics* metrics,
                                         int i,
                                         int j,
                                         int k,
                                         Env* env) {
  Assert(metrics);
  Assert(i >= 0);
  Assert(i < metrics->num_vector_local);
  Assert(j >= 0);
  Assert(j < metrics->num_vector_local);
  Assert(k >= 0);
  Assert(k < metrics->num_vector_local);
  Assert(env);

  Assert(i != j);
  Assert(i != k);

  size_t index = (i*(size_t)(i-1)*(size_t)(i-2))/6 + (j*(size_t)(j-1))/2 + k;



/*
FIX
      int index = i>j ? ( (i*(i-1))>>1 ) + j : ( (j*(j-1))>>1 ) + i;
  */
  return index;
}

/*===========================================================================*/
/*---Accessors: coord from index---*/

static int Metrics_coord0_from_index_2(Metrics* metrics,
                                       size_t index,
                                       Env* env) {
  Assert(metrics);
  Assert(index >= 0);
  Assert(index < metrics->num_elts_local);
  Assert(env);
  Assert(env->num_way == 2);

  int i = metrics->index_map[index] % metrics->num_elts_local;
  return i;
}

/*---------------------------------------------------------------------------*/

static int Metrics_coord1_from_index_2(Metrics* metrics,
                                       size_t index,
                                       Env* env) {
  Assert(metrics);
  Assert(index >= 0);
  Assert(index < metrics->num_elts_local);
  Assert(env);
  Assert(env->num_way == 2);

  int j = metrics->index_map[index] / metrics->num_elts_local;
  return j;
}

/*---------------------------------------------------------------------------*/

static int Metrics_coord0_from_index_3(Metrics* metrics,
                                       size_t index,
                                       Env* env) {
  Assert(metrics);
  Assert(index >= 0);
  Assert(index < metrics->num_elts_local);
  Assert(env);
  Assert(env->num_way == 3);

  int i = metrics->index_map[index] % metrics->num_elts_local;
  return i;
}

/*---------------------------------------------------------------------------*/

static int Metrics_coord1_from_index_3(Metrics* metrics,
                                       size_t index,
                                       Env* env) {
  Assert(metrics);
  Assert(index >= 0);
  Assert(index < metrics->num_elts_local);
  Assert(env);
  Assert(env->num_way == 3);

  int i = ( metrics->index_map[index] / metrics->num_elts_local )
                                      % metrics->num_elts_local;
  return i;
}

/*---------------------------------------------------------------------------*/

static int Metrics_coord2_from_index_3(Metrics* metrics,
                                       size_t index,
                                       Env* env) {
  Assert(metrics);
  Assert(index >= 0);
  Assert(index < metrics->num_elts_local);
  Assert(env);
  Assert(env->num_way == 3);

  int i = metrics->index_map[index] / ( metrics->num_elts_local *
                                        metrics->num_elts_local );
  return i;
}

/*---------------------------------------------------------------------------*/
static int Metrics_coord_from_index(Metrics* metrics,
                                    size_t index,
                                    int coord_num,
                                    Env* env) {
  Assert(metrics);
  Assert(index >= 0);
  Assert(index < metrics->num_elts_local);
  Assert(env);
  Assert(coord_num >= 0);
  Assert(coord_num < env->num_way);

  int result = 0;

  switch (env->num_way + 4 * coord_num) {
    case 2:   /* 2-way, coord 0 */
      result = Metrics_coord0_from_index_2(metrics, index, env);
      break;
    case 6:   /* 2-way, coord 1 */
      result = Metrics_coord1_from_index_2(metrics, index, env);
      break;
    case 3:   /* 3-way, coord 0 */
      result = Metrics_coord0_from_index_3(metrics, index, env);
      break;
    case 7:   /* 3-way, coord 1 */
      result = Metrics_coord1_from_index_3(metrics, index, env);
      break;
    case 11:   /* 3-way, coord 2 */
      result = Metrics_coord2_from_index_3(metrics, index, env);
      break;
    default:
      Insist(env, Bool_false ? "Unimplemented." : 0);
  } /*---case---*/

  return result;
}

/*===========================================================================*/

#endif /*---_metrics_h_---*/

/*---------------------------------------------------------------------------*/
