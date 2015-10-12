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

void gm_compute_metrics_czekanowski_2way_cpu(GMMetrics* metrics,
                                          GMVectors* vectors,
                                          GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, (!env->all2all) ? "Unimplemented." : 0);

  /*---Denominator---*/

  GMFloat* vector_sums = malloc(metrics->num_vector_local * sizeof(GMFloat));

  gm_compute_float_vector_sums(vectors, vector_sums, env);

  /*---Numerator---*/

  int i = 0;
  int j = 0;
  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      GMFloat sum = 0;
      int field = 0;
      for (field = 0; field < vectors->num_field; ++field) {
        const GMFloat value1 = GMVectors_float_get(vectors, field, i, env);
        const GMFloat value2 = GMVectors_float_get(vectors, field, j, env);
        sum += value1 < value2 ? value1 : value2;
      } /*---for k---*/
      GMMetrics_float_set_2(metrics, i, j, sum, env);
    } /*---for j---*/
  }   /*---for i---*/

  /*---Combine---*/

  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      const GMFloat numerator = GMMetrics_float_get_2(metrics, i, j, env);
      const GMFloat denominator = vector_sums[i] + vector_sums[j];
      GMMetrics_float_set_2(metrics, i, j, 2 * numerator / denominator, env);
    } /*---for j---*/
  }   /*---for i---*/

  free(vector_sums);
}

/*===========================================================================*/

void gm_compute_metrics_czekanowski_2way_gpu(GMMetrics* metrics,
                                          GMVectors* vectors,
                                          GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, (!env->all2all) ? "Unimplemented." : 0);

  /*---Denominator---*/

  GMFloat* vector_sums = malloc(metrics->num_vector_local * sizeof(GMFloat));

  /* .02 / 1.56 */
  gm_compute_float_vector_sums(vectors, vector_sums, env);

  /*---------------*/
  /*---Numerator---*/
  /*---------------*/

  int i = 0;

  magma_minproduct_int_t magma_code = 0;
  if ( magma_code ) {} /*---Avoid unused variable warning---*/

  magma_code = magma_minproduct_init();
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);

  const int numvec = metrics->num_vector_local;
  const int numfield = vectors->num_field;

  GMFloat* h_vectors = NULL;
  GMFloat* h_numer = NULL;

/*---Allocate magma CPU memory for vectors and for result */

#ifdef FP_PRECISION_DOUBLE
  magma_code = magma_minproduct_dmalloc_pinned(&h_vectors, numvec * numfield);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
  magma_code = magma_minproduct_dmalloc_pinned(&h_numer, numvec * numvec);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_code = magma_minproduct_smalloc_pinned(&h_vectors, numvec * numfield);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
  magma_code = magma_minproduct_smalloc_pinned(&h_numer, numvec * numvec);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
#endif

  /*---Allocate GPU mirrors for CPU arrays---*/

  GMFloat* d_vectors = NULL;
  GMFloat* d_numer = NULL;

#ifdef FP_PRECISION_DOUBLE
  magma_code = magma_minproduct_dmalloc(&d_vectors, numvec * numfield);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
  magma_code = magma_minproduct_dmalloc(&d_numer, numvec * numvec);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_code = magma_minproduct_smalloc(&d_vectors, numvec * numfield);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
  magma_code = magma_minproduct_smalloc(&d_numer, numvec * numvec);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
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
      h_vectors[k + numfield * i] = GMVectors_float_get(vectors, k, i, env);
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
      const GMFloat numerator = h_numer[j + numvec * i];
      const GMFloat denominator = vector_sums[i] + vector_sums[j];
      GMMetrics_float_set_2(metrics, i, j, 2 * numerator / denominator, env);
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

void gm_compute_metrics_czekanowski_3way_cpu(GMMetrics* metrics,
                                          GMVectors* vectors,
                                          GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, (!env->all2all) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

void gm_compute_metrics_czekanowski_3way_gpu(GMMetrics* metrics,
                                          GMVectors* vectors,
                                          GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, (!env->all2all) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
