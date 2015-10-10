/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_czekanowski_2way.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Functions for computing 2-way Czekanowski metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h> /*FIX*/
#include <stdlib.h>

#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"

#include "env.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_metrics_utils.h"
#include "compute_metrics_czekanowski_2way.h"

/*===========================================================================*/

void compute_metrics_czekanowski_2way_cpu(Metrics* metrics,
                                          Vectors* vectors,
                                          Env* env) {
  Assert(metrics != NULL);
  Assert(vectors != NULL);
  Assert(env != NULL);

  Insist(env, (!env->all2all) ? "Unimplemented." : 0);

  /*---Denominator---*/

  Float_t* vector_sums = malloc(metrics->num_vector_local * sizeof(Float_t));

  compute_vector_sums(vectors, vector_sums, env);

  /*---Numerator---*/

  int i = 0;
  int j = 0;
  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      Float_t sum = 0;
      int field = 0;
      for (field = 0; field < vectors->num_field; ++field) {
        const Float_t value1 = Vectors_float_get(vectors, field, i, env);
        const Float_t value2 = Vectors_float_get(vectors, field, j, env);
        sum += value1 < value2 ? value1 : value2;
      } /*---for k---*/
      Metrics_Float_set_2(metrics, i, j, sum, env);
    } /*---for j---*/
  }   /*---for i---*/

  /*---Combine---*/

  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      const Float_t numerator = Metrics_Float_get_2(metrics, i, j, env);
      const Float_t denominator = vector_sums[i] + vector_sums[j];
      Metrics_Float_set_2(metrics, i, j, 2 * numerator / denominator, env);
    } /*---for j---*/
  }   /*---for i---*/

  free(vector_sums);
}

/*===========================================================================*/

void compute_metrics_czekanowski_2way_gpu(Metrics* metrics,
                                          Vectors* vectors,
                                          Env* env) {
  Assert(metrics != NULL);
  Assert(vectors != NULL);
  Assert(env != NULL);

  Insist(env, (!env->all2all) ? "Unimplemented." : 0);

  /*---Denominator---*/

  Float_t* vector_sums = malloc(metrics->num_vector_local * sizeof(Float_t));

  /* .02 / 1.56 */
  compute_vector_sums(vectors, vector_sums, env);

  /*---------------*/
  /*---Numerator---*/
  /*---------------*/

  int i = 0;

  magma_minproduct_int_t magma_code = 0;
  if ( magma_code ) {} /*---Avoid unused variable warning---*/

  magma_code = magma_minproduct_init();
  Assert(magma_code == MAGMA_minproduct_SUCCESS);

  const int numvec = metrics->num_vector_local;
  const int numfield = vectors->num_field;

  Float_t* h_vectors = NULL;
  Float_t* h_numer = NULL;

/*---Allocate magma CPU memory for vectors and for result */

#ifdef FP_PRECISION_DOUBLE
  magma_code = magma_minproduct_dmalloc_pinned(&h_vectors, numvec * numfield);
  Assert(magma_code == MAGMA_minproduct_SUCCESS);
  magma_code = magma_minproduct_dmalloc_pinned(&h_numer, numvec * numvec);
  Assert(magma_code == MAGMA_minproduct_SUCCESS);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_code = magma_minproduct_smalloc_pinned(&h_vectors, numvec * numfield);
  Assert(magma_code == MAGMA_minproduct_SUCCESS);
  magma_code = magma_minproduct_smalloc_pinned(&h_numer, numvec * numvec);
  Assert(magma_code == MAGMA_minproduct_SUCCESS);
#endif

  /*---Allocate GPU mirrors for CPU arrays---*/

  Float_t* d_vectors = NULL;
  Float_t* d_numer = NULL;

#ifdef FP_PRECISION_DOUBLE
  magma_code = magma_minproduct_dmalloc(&d_vectors, numvec * numfield);
  Assert(magma_code == MAGMA_minproduct_SUCCESS);
  magma_code = magma_minproduct_dmalloc(&d_numer, numvec * numvec);
  Assert(magma_code == MAGMA_minproduct_SUCCESS);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_code = magma_minproduct_smalloc(&d_vectors, numvec * numfield);
  Assert(magma_code == MAGMA_minproduct_SUCCESS);
  magma_code = magma_minproduct_smalloc(&d_numer, numvec * numvec);
  Assert(magma_code == MAGMA_minproduct_SUCCESS);
#endif

  /*---Initialize result matrix to zero (apparently magma requires)---*/

#ifdef FP_PRECISION_DOUBLE
  magma_minproductblas_dlaset(Magma_minproductFull, numvec, numvec, 0.0, 0.0, d_numer, numvec);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproductblas_slaset(Magma_minproductFull, numvec, numvec, 0.0, 0.0, d_numer, numvec);
#endif

  /*---Copy in vectors---*/

  /* .08 / 1.56 */
  for (i = 0; i < numvec; ++i) {
    int k = 0;
    for (k = 0; k < numfield; ++k) {
      h_vectors[k + numfield * i] = Vectors_float_get(vectors, k, i, env);
    }
  }

/*---Send vectors to GPU---*/

#ifdef FP_PRECISION_DOUBLE
  magma_minproduct_dsetmatrix(numfield, numvec, h_vectors,
                              numfield, d_vectors, numfield);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproduct_ssetmatrix(numfield, numvec, h_vectors,
                              numfield, d_vectors, numfield);
#endif

/*---Perform pseudo matrix-matrix product---*/

/* .63 / 1.56 */
#ifdef FP_PRECISION_DOUBLE
  magma_minproductblas_dgemm
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproductblas_sgemm
#endif
      (Magma_minproductTrans, Magma_minproductNoTrans, numvec, numvec, numfield,
       1.0, d_vectors, numfield, d_vectors, numfield, 0.0, d_numer, numvec);

/*---Copy result from GPU---*/

#ifdef FP_PRECISION_DOUBLE
  magma_minproduct_dgetmatrix(numvec, numvec, d_numer, numvec,
                              h_numer, numvec);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproduct_sgetmatrix(numvec, numvec, d_numer, numvec,
                              h_numer, numvec);
#endif

  /*---Combine---*/

  /* .22 / 1.56 */
  for (i = 0; i < metrics->num_vector_local; ++i) {
    int j = 0;
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      const Float_t numerator = h_numer[j + numvec * i];
      const Float_t denominator = vector_sums[i] + vector_sums[j];
      Metrics_Float_set_2(metrics, i, j, 2 * numerator / denominator, env);
    } /*---for j---*/
  }   /*---for i---*/

  /*---Free memory---*/

  magma_minproduct_free(d_vectors);
  magma_minproduct_free(d_numer);

  magma_minproduct_free_pinned(h_vectors);
  magma_minproduct_free_pinned(h_numer);

  magma_minproduct_finalize();

  free(vector_sums);
}

/*===========================================================================*/

void compute_metrics_czekanowski_3way_cpu(Metrics* metrics,
                                          Vectors* vectors,
                                          Env* env) {
  Assert(metrics != NULL);
  Assert(vectors != NULL);
  Assert(env != NULL);

  Insist(env, (!env->all2all) ? "Unimplemented." : 0);

  Insist(env, Bool_false ? "Unimplemented." : 0);
}

/*===========================================================================*/

void compute_metrics_czekanowski_3way_gpu(Metrics* metrics,
                                          Vectors* vectors,
                                          Env* env) {
  Assert(metrics != NULL);
  Assert(vectors != NULL);
  Assert(env != NULL);

  Insist(env, (!env->all2all) ? "Unimplemented." : 0);

  Insist(env, Bool_false ? "Unimplemented." : 0);
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
