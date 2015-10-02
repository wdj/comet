/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Functions for computing metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h>


#include <stdlib.h>

#include "env.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_metrics.h"

/*===========================================================================*/

void compute_metrics(Metrics* metrics, Vectors* vectors, Env* env) {
  Assert(metrics != NULL);
  Assert(vectors != NULL);
  Assert(env != NULL);

  switch (env->metric_type + 10 * ( env->num_way + 10 * env->compute_method )) {
    case METRIC_TYPE_SORENSON + 10 * ( 2 + 10 * COMPUTE_METHOD_CPU ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;
    case METRIC_TYPE_SORENSON + 10 * ( 2 + 10 * COMPUTE_METHOD_GPU ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;
    case METRIC_TYPE_SORENSON + 10 * ( 3 + 10 * COMPUTE_METHOD_CPU ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;
    case METRIC_TYPE_SORENSON + 10 * ( 3 + 10 * COMPUTE_METHOD_GPU ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;
    case METRIC_TYPE_CZEKANOWSKI + 10 * ( 2 + 10 * COMPUTE_METHOD_CPU ):
      compute_metrics_czek_2way_cpu(metrics, vectors, env);
      break;
    case METRIC_TYPE_CZEKANOWSKI + 10 * ( 2 + 10 * COMPUTE_METHOD_GPU ):
      compute_metrics_czek_2way_gpu(metrics, vectors, env);
      break;
    case METRIC_TYPE_CZEKANOWSKI + 10 * ( 3 + 10 * COMPUTE_METHOD_CPU ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;
    case METRIC_TYPE_CZEKANOWSKI + 10 * ( 3 + 10 * COMPUTE_METHOD_GPU ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;
    case METRIC_TYPE_CCC + 10 * ( 2 + 10 * COMPUTE_METHOD_CPU ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;
    case METRIC_TYPE_CCC + 10 * ( 2 + 10 * COMPUTE_METHOD_GPU ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;
    case METRIC_TYPE_CCC + 10 * ( 3 + 10 * COMPUTE_METHOD_CPU ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;
    case METRIC_TYPE_CCC + 10 * ( 3 + 10 * COMPUTE_METHOD_GPU ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;
    default:
      Insist(env, Bool_false ? "Unimplemented." : 0);
  } /*---switch---*/
}

/*===========================================================================*/

void compute_metrics_czek_2way_cpu(Metrics* metrics,
                                   Vectors* vectors,
                                   Env* env) {
  Assert(metrics != NULL);
  Assert(vectors != NULL);
  Assert(env != NULL);

  Insist(env, ( ! env->global_all2all ) ? "Unimplemented." : 0);

  Float_t* vector_sums = malloc(metrics->num_vector_local*sizeof(Float_t));

  /*---Denominator---*/

  int i = 0;
  for ( i = 0; i < metrics->num_vector_local; ++i ) {
    Float_t sum = 0;
    int field = 0;
    for ( field = 0; field < vectors->num_field; ++field ) {
      Float_t value = Vectors_float_get(vectors, field, i, env);
      sum += value;
    }
    vector_sums[i] = sum;
  }

  /*---Numerator---*/

  int j = 0;
  for ( i = 0; i < metrics->num_vector_local; ++i ) {
    for ( j = i+1; j < metrics->num_vector_local; ++j ) {
      Float_t sum = 0;
      int field = 0;
      for ( field = 0; field < vectors->num_field; ++field ) {
        const Float_t value1 = Vectors_float_get(vectors, field, i, env);
        const Float_t value2 = Vectors_float_get(vectors, field, j, env);
        sum += value1 < value2 ? value1 : value2;
      } /*---for k---*/
      Metrics_Float_set_2(metrics, i, j, sum, env);
    } /*---for j---*/
  } /*---for i---*/

  /*---Combine---*/

  for ( i = 0; i < metrics->num_vector_local; ++i ) {
    for ( j = i+1; j < metrics->num_vector_local; ++j ) {
      const Float_t numerator = Metrics_Float_get_2(metrics, i, j, env);
      const Float_t denominator = vector_sums[i] + vector_sums[j];
      Metrics_Float_set_2(metrics, i, j, 2 * numerator / denominator, env);
    } /*---for j---*/
  } /*---for i---*/

  free( vector_sums );
}

/*===========================================================================*/

void compute_metrics_czek_2way_gpu(Metrics* metrics,
                                   Vectors* vectors,
                                   Env* env) {
  Assert(metrics != NULL);
  Assert(vectors != NULL);
  Assert(env != NULL);

  Insist(env, ( ! env->global_all2all ) ? "Unimplemented." : 0);




  Insist(env, Bool_false ? "Unimplemented." : 0);


}

/*===========================================================================*/
/* Compute simple sum of elements of a vector. */

Float_t vector_sum(int len, Float_t * const __restrict__ v1) {
  int i = 0;
  Float_t result = 0;

  for(i = 0; i < len; ++i) {
    result += ( v1[i] );
  }

  return result;
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
