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
#include "compute_metrics_czekanowski_2way.h"
#include "compute_metrics_czekanowski_3way.h"

/*===========================================================================*/

void compute_metrics(Metrics* metrics, Vectors* vectors, Env* env) {
  Assert(metrics != NULL);
  Assert(vectors != NULL);
  Assert(env != NULL);

  switch (env->metric_type + NUM_METRIC_TYPE * (
          env->compute_method + NUM_COMPUTE_METHOD * ( env->num_way ))) {

    case METRIC_TYPE_SORENSON + NUM_METRIC_TYPE * (
         COMPUTE_METHOD_CPU + NUM_COMPUTE_METHOD * ( 2 ) ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;

    case METRIC_TYPE_SORENSON + NUM_METRIC_TYPE * (
         COMPUTE_METHOD_GPU + NUM_COMPUTE_METHOD * ( 2 ) ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;

    case METRIC_TYPE_SORENSON + NUM_METRIC_TYPE * (
         COMPUTE_METHOD_CPU + NUM_COMPUTE_METHOD * ( 3 ) ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;

    case METRIC_TYPE_SORENSON + NUM_METRIC_TYPE * (
         COMPUTE_METHOD_GPU + NUM_COMPUTE_METHOD * ( 3 ) ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;

    case METRIC_TYPE_CZEKANOWSKI + NUM_METRIC_TYPE * (
         COMPUTE_METHOD_CPU + NUM_COMPUTE_METHOD * ( 2 ) ):
      compute_metrics_czekanowski_2way_cpu(metrics, vectors, env);
      break;

    case METRIC_TYPE_CZEKANOWSKI + NUM_METRIC_TYPE * (
         COMPUTE_METHOD_GPU + NUM_COMPUTE_METHOD * ( 2 ) ):
      compute_metrics_czekanowski_2way_gpu(metrics, vectors, env);
      break;

    case METRIC_TYPE_CZEKANOWSKI + NUM_METRIC_TYPE * (
         COMPUTE_METHOD_CPU + NUM_COMPUTE_METHOD * ( 3 ) ):
      compute_metrics_czekanowski_3way_cpu(metrics, vectors, env);
      break;

    case METRIC_TYPE_CZEKANOWSKI + NUM_METRIC_TYPE * (
         COMPUTE_METHOD_GPU + NUM_COMPUTE_METHOD * ( 3 ) ):
      compute_metrics_czekanowski_3way_gpu(metrics, vectors, env);
      break;

    case METRIC_TYPE_CCC + NUM_METRIC_TYPE * (
         COMPUTE_METHOD_CPU + NUM_COMPUTE_METHOD * ( 2 ) ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;

    case METRIC_TYPE_CCC + NUM_METRIC_TYPE * (
         COMPUTE_METHOD_GPU + NUM_COMPUTE_METHOD * ( 2 ) ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;

    case METRIC_TYPE_CCC + NUM_METRIC_TYPE * (
         COMPUTE_METHOD_CPU + NUM_COMPUTE_METHOD * ( 3 ) ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;

    case METRIC_TYPE_CCC + NUM_METRIC_TYPE * (
         COMPUTE_METHOD_GPU + NUM_COMPUTE_METHOD * ( 3 ) ):
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;

    default:
      Insist(env, Bool_false ? "Unimplemented." : 0);
  } /*---switch---*/
}

/*===========================================================================*/
/* Compute simple sum of elements of a vector. */

/*
Float_t vector_sum(int len, Float_t * const __restrict__ v1) {
  int i = 0;
  Float_t result = 0;

  for(i = 0; i < len; ++i) {
    result += ( v1[i] );
  }

  return result;
}
*/

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
