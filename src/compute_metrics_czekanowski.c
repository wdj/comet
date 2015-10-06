/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_czekanowski.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Functions for computing Czekanowski metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

/*FIX*/
#include <stdio.h>

#include <stdlib.h>

#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"

#include "env.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_metrics_czekanowski.h"

/*===========================================================================*/

void compute_metrics_czekanowski_2way_cpu(Metrics* metrics,
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

void compute_metrics_czekanowski_2way_gpu(Metrics* metrics,
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

  /*---------------*/
  /*---Numerator---*/
  /*---------------*/

  magma_minproduct_init();

#ifdef FP_PRECISION_DOUBLE

  const int numvec = metrics->num_vector_local;
  const int numfield = vectors->num_field;

  Float_t* h_vectors = NULL;
  Float_t* h_numer = NULL;

  /*---Allocate magma CPU memory for vectors and for result */
  magma_minproduct_dmalloc_pinned(&h_vectors,numvec*numfield);
  magma_minproduct_dmalloc_pinned(&h_numer,numvec*numvec);

  /*---Copy in vectors---*/

  for (i = 0; i < numvec; ++i) {
    int k = 0;
    for (k = 0; k < numfield; ++k) {
      h_vectors[k+numfield*i] = Vectors_float_get( vectors, k, i, env );
    }
  }

  /*---Allocate GPU mirrors for CPU arrays---*/

  Float_t* d_vectors = NULL;
  Float_t* d_numer = NULL;

  magma_minproduct_dmalloc(&d_vectors, numvec*numfield);
  magma_minproduct_dmalloc(&d_numer, numvec*numvec);

  /*---Send vectors to GPU---*/

  magma_minproduct_dsetmatrix(numfield, numvec, h_vectors, numfield,
                                   d_vectors, numfield);

  /*---Initialize result matrix to zero (apparently magma requires)---*/

  for (i = 0; i < numvec; ++i) {
    int j = 0;
    for (j = 0; j < numvec; ++j) {
      h_numer[j+numvec*i] = 0;
    }
  }

  /*---Send initialized result matrix to GPU---*/

  magma_minproduct_dsetmatrix(numvec, numvec, h_numer, numvec,
                              d_numer, numvec);

  /*---Perform pseudo matrix-matrix product---*/

  magma_minproductblas_dgemm(
    Magma_minproductTrans,
    Magma_minproductNoTrans,
    numvec,
    numvec,
    numfield,
    1.0,
    d_vectors,
    numfield,
    d_vectors,
    numfield,
    0.0,
    d_numer,
    numvec );

  /*---Copy result from GPU---*/

  magma_minproduct_dgetmatrix(numvec, numvec, d_numer, numvec,
                              h_numer, numvec);

  /*---Combine---*/

  for ( i = 0; i < metrics->num_vector_local; ++i ) {
    int j = 0;
    for ( j = i+1; j < metrics->num_vector_local; ++j ) {
      const Float_t numerator = h_numer[j+numvec*i];
      const Float_t denominator = vector_sums[i] + vector_sums[j];
      Metrics_Float_set_2(metrics, i, j, 2 * numerator / denominator, env);
    } /*---for j---*/
  } /*---for i---*/

  /*---Free memory---*/

  magma_minproduct_free(d_vectors);
  magma_minproduct_free(d_numer);

  magma_minproduct_free_pinned(h_vectors);
  magma_minproduct_free_pinned(h_numer);
#endif

#ifdef FP_PRECISION_SINGLE
  Insist(env, Bool_false ? "Unimplemented." : 0);
#endif

  magma_minproduct_finalize();

  free( vector_sums );
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
