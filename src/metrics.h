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

#include "env.h"

/*===========================================================================*/
/*---Struct declaration---*/

typedef struct {
  int num_vector;
  int num_vector_local;
  int num_vector_local_max;
  int num_elts_local;
  int data_type_id;
  int* index_map;
  void* __restrict__ data;
} Metrics;

/*===========================================================================*/
/*---Functions---*/

void Metrics_create(Metrics* metrics,
                    int data_type_id,
                    int num_vector_local,
                    Env* env);

void Metrics_destroy(Metrics* metrics, Env* env);

/*---------------------------------------------------------------------------*/
/*---Accessors---*/

int Metrics_index_from_coord_2(Metrics* metrics, int i, int j, Env* env) {
  Assert(metrics);
  Assert(i >= 0);
  Assert(i < metrics->num_vector_local);
  Assert(j >= 0);
  Assert(j < metrics->num_vector_local);
  Assert(env);

  Assert(i != j);

  int index = i > j ? ((i * (i - 1)) >> 1) + j : ((j * (j - 1)) >> 1) + i;
  return index;
}

/*---------------------------------------------------------------------------*/

int Metrics_coord0_from_index_2(Metrics* metrics, int index, Env* env) {
  Assert(metrics);
  Assert(index >= 0);
  Assert(index < metrics->num_elts_local);
  Assert(env);

  int i = metrics->index_map[index] / metrics->num_elts_local;
  return i;
}

/*---------------------------------------------------------------------------*/

int Metrics_coord1_from_index_2(Metrics* metrics, int index, Env* env) {
  Assert(metrics);
  Assert(index >= 0);
  Assert(index < metrics->num_elts_local);
  Assert(env);

  int j = metrics->index_map[index] % metrics->num_elts_local;
  return j;
}

/*---------------------------------------------------------------------------*/

int Metrics_index_from_coord_3(Metrics* metrics,
                               int i,
                               int j,
                               int k,
                               Env* env) {
  Assert(metrics);
  Assert(i >= 0);
  Assert(i < metrics->num_vector_local);
  Assert(j >= 0);
  Assert(j < metrics->num_vector_local);
  Assert(env);

  Assert(i != j);

  /*
      int index = (i*(i-1)*(i-2))/6 + (j*(j-1))/2 + k;



      int index = i>j ? ( (i*(i-1))>>1 ) + j : ( (j*(j-1))>>1 ) + i;
      return index;
  */
}

/*===========================================================================*/

#endif /*---_metrics_h_---*/

/*---------------------------------------------------------------------------*/
