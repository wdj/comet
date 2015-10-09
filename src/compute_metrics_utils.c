/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_utils.c
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing metrics, utilities.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include "env.h"
#include "vectors.h"
#include "metrics.h"

/*===========================================================================*/

void compute_vector_sums(Vectors* vectors,
                         Float_t* __restrict__ vector_sums,
                         Env* env) {
  Assert(vectors != NULL);
  Assert(vector_sums != NULL);
  Assert(env != NULL);

  int i = 0;
  for (i = 0; i < vectors->num_vector_local; ++i) {
    Float_t sum = 0;
    int field = 0;
    for (field = 0; field < vectors->num_field; ++field) {
      Float_t value = Vectors_float_get(vectors, field, i, env);
      sum += value;
    }
    vector_sums[i] = sum;
  }
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
