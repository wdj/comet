/****************************************************************/
/* 3-Way Czekanowski Similarity Metric                          */
/*                                                              */
/* Code Author: James Nance                                     */
/* Created: October, 2015                                       */
/****************************************************************/

#include <stdio.h> /*FIX*/
#include <stdlib.h>

#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"

#include "env.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_metrics_utils.h"
#include "compute_metrics_czekanowski_3way.h"

/*===========================================================================*/

void gm_compute_metrics_czekanowski_3way_cpu(GMMetrics* metrics,
                                             GMVectors* vectors,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, (!env->all2all) ? "Unimplemented." : 0);

  /*---Denominator---*/

  GMFloat* vector_sums = GMFLoat_malloc((metrics->num_vector_local);

  gm_compute_float_vector_sums(vectors, vector_sums, env);

  /*---Numerator---*/

  int i = 0;
  int j = 0;
  int k = 0;
  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      for (k = j + 1; k < metrics->num_vector_local; ++k) {
        GMFloat sum = 0;
        int field = 0;
        for (field = 0; field < vectors->num_field; ++field) {
          const GMFloat value1 = GMVectors_float_get(vectors, field, i, env);
          const GMFloat value2 = GMVectors_float_get(vectors, field, j, env);
          const GMFloat value3 = GMVectors_float_get(vectors, field, k, env);
          GMFloat min12 = value1 < value2 ? value1 : value2;
          sum += min12;
          sum += value1 < value3 ? value1 : value3;
          sum += value2 < value3 ? value2 : value3;
          sum -= min12 < value3 ? min12 : value3;
        } /*---for field---*/
        GMMetrics_float_set_3(metrics, i, j, k, sum, env);
      } /*---for k---*/
    }   /*---for j---*/
  }     /*---for i---*/

  /*---Combine---*/

  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      for (k = j + 1; k < metrics->num_vector_local; ++k) {
        const GMFloat numerator = GMMetrics_float_get_3(metrics, i, j, k, env);
        const GMFloat denominator =
            vector_sums[i] + vector_sums[j] + vector_sums[k];
        GMMetrics_float_set_3(metrics, i, j, k,
                              3 * numerator / (2 * denominator), env);
      } /*---for k---*/
    }   /*---for j---*/
  }     /*---for i---*/

  free(vector_sums);
}

/*===========================================================================*/

void gm_compute_metrics_czekanowski_3way_gpu(GMMetrics* metrics,
                                             GMVectors* vectors,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, (!env->all2all) ? "Unimplemented." : 0);

  /*---Denominator---*/

  GMVectorSums vector_sums = GMVectorSums_null();
  GMVectorSums_create(&vector_sums, vectors, env);
  gm_compute_vector_sums(vectors, &vector_sums, env);

  /*---Numerator---*/

  int i = 0;
  int j = 0;
  int k = 0;

  const int numvec = metrics->num_vector_local;
  const int numfield = vectors->num_field;

  /*---Initialize MAGMA library---*/
  gm_magma_initialize(env);

  /*---Allocate magma CPU memory for vectors and for result---*/
  GMMirroredPointer h_vectors = 
    gm_malloc_magma ( numvec * (size_t)numfield, env); //Data matrix X
  GMMirroredPointer h_matM = 
    gm_malloc_magma ( numvec * (size_t)numvec, env); // M = X^T minprod X
  GMMirroredPointer h_matV = 
    gm_malloc_magma ( numvec * (size_t)numfield, env); // V = elementwise min of X and one column of X
  GMMirroredPointer h_matB =
    gm_malloc_magma ( numvec * (size_t)numvec, env); // B = X^T minprod V

  /*---Copy in vectors---*/
  gm_vectors_to_buf(vectors, &h_vectors, env);

  /*---Send matrix vectors to GPU---*/
  gm_set_vectors_start(vectors, &h_vectors, env);
  gm_set_vectors_wait(env);


  /*---Perform pseudo matrix-matrix product matM=vectors^T*vectors---*/
  /*
    My matM is the same as the 2way "numerators" variable, so I should 
    be able to use the same function from the 2-way implementation.
  */
  gm_compute_numerators_start(vectors, vectors, 
       metrics, &h_vectors, &h_vectors, &h_matM,
       env->proc_num, GM_BOOL_TRUE, env);

  gm_compute_numerators_wait(env);

/*---Copy matM from GPU---*/

  gm_get_metrics_start(metrics, &h_matM, env);
  gm_get_metrics_wait(env);

  /*---Process all combinations starting with j, i, k---*/
  for (j = 1; j < numvec - 1; ++j) {
    /*---Populate first j-1 columns of matV---*/
    for (i = 0; i < numvec; ++i) {
      // Compare columns x_i and x_j element-wise
      for (k = 0; k < numfield; ++k) {
        // row = k, column = i
        h_matV[k + i * numfield] =
            h_vectors[k + numfield * i] < h_vectors[k + numfield * j]
                ? h_vectors[k + numfield * i]
                : h_vectors[k + numfield * j];
      }
    }
    /*---Send matrix matV to GPU---*/
    gm_set_matrix_start(h_matV, numfield, numvec, env);
    gm_set_matrix_wait(env);

    /*---Initialize result matrix to zero (apparently magma requires)---*/

#ifdef FP_PRECISION_DOUBLE
  magma_minproductblas_dlaset(Magma_minproductFull, numvec, numvec, 0.0, 0.0,
                              d_matB, numvec);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproductblas_slaset(Magma_minproductFull, numvec, numvec, 0.0, 0.0,
                              d_matB, numvec);
#endif

/*---Perform matrix-matrix product matB = matV^T*vectors---*/
#ifdef FP_PRECISION_DOUBLE
    magma_minproductblas_dgemm
#endif
#ifdef FP_PRECISION_SINGLE
        magma_minproductblas_sgemm
#endif
        (Magma_minproductTrans, Magma_minproductNoTrans, numvec, numvec,
         numfield, 1.0, d_matV, numfield, d_vectors, numfield, 0.0, d_matB,
         numvec);
  gm_compute_numerators_wait(env); /*---added - WJ---*/

/*---Copy matB from GPU---*/
  /*
  gm_get_matrix_start(h_matB, numvec, numvec, env);
  gm_get_matrix_wait(env);
  */

#ifdef FP_PRECISION_DOUBLE
    magma_minproduct_dgetmatrix(numvec, numvec, d_matB, numvec, h_matB, numvec);
#endif
#ifdef FP_PRECISION_SINGLE
    magma_minproduct_sgetmatrix(numvec, numvec, d_matB, numvec, h_matB, numvec);
#endif

    /*---Combine results---*/
    for (i = 0; i < j; ++i) {
      GMFloat min_ij = h_matM[j + numvec * i];
      for (k = j + 1; k < numvec; ++k) {
        GMFloat min_ik = h_matM[k + numvec * i];
        GMFloat min_jk = h_matM[k + numvec * j];
        // sum of mins vectors i, j, and k is matB(i,k)
        GMFloat min_ijk = h_matB[k + numvec * i];
        const GMFloat numerator = min_ij + min_ik + min_jk - min_ijk;
        const GMFloat denominator =
            vector_sums[i] + vector_sums[j] + vector_sums[k];
        GMMetrics_float_set_3(metrics, i, j, k,
                              3 * numerator / (2 * denominator), env);
      } /*---for k---*/
    }   /*---for i---*/
  }     /*---for j---*/

/*
  gm_compute_combine(metrics, &h_numerators,
                     &vector_sums, &vector_sums, env->proc_num,
                     GM_BOOL_TRUE, env);
*/

  /*---Free memory and finalize---*/
  magma_minproduct_free(d_matM);
  magma_minproduct_free(d_vectors);
  magma_minproduct_free(d_matV);
  magma_minproduct_free(d_matB);
  magma_minproduct_free_pinned(h_matM);
  magma_minproduct_free_pinned(h_vectors);
  magma_minproduct_free_pinned(h_matB);
  magma_minproduct_free_pinned(h_matV);

  magma_minproduct_finalize();

  free(vector_sums);
}
