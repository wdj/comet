/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_czekanowski_3way.c
 * \author James Nance, Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing 3-way Czekanowski metrics.
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
#include "compute_metrics_czekanowski_3way.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/

/*---NOTE: This routine currently only handles Czekanowski, but the
    intent is that it will later be adapted to handle all the metrics---*/

void gm_compute_metrics_czekanowski_3way_all2all(GMMetrics* metrics,
                                                 GMVectors* vectors,
                                                 GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);


  /*---Initializations---*/

#ifdef i_am_undefined


  int i = 0;
  int i_proc = env->proc_num;















#endif



  /*---TODO: remove when ready---*/
  GMInsist(env, (!env->all2all) ? "Unimplemented." : 0);


}

/*===========================================================================*/

void gm_compute_metrics_czekanowski_3way_cpu(GMMetrics* metrics,
                                             GMVectors* vectors,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  if (env->all2all) {
    gm_compute_metrics_czekanowski_3way_all2all(metrics, vectors, env);
    return;
  }

  /*---Denominator---*/

  GMFloat* vector_sums = GMFloat_malloc(metrics->num_vector_local);

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
  //printf("---CPU Implementation---\n");
  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      for (k = j + 1; k < metrics->num_vector_local; ++k) {
        const GMFloat numerator = GMMetrics_float_get_3(metrics, i, j, k, env);
        const GMFloat denominator =
            vector_sums[i] + vector_sums[j] + vector_sums[k];
        //printf("%i,%i,%i . . . numerator = %f . . . denominator = %f\n",i,j,k,numerator,denominator);
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

  if (env->all2all) {
    gm_compute_metrics_czekanowski_3way_all2all(metrics, vectors, env);
    return;
  }

  //printf("---GPU Implementation---\n");

  /*---Denominator---*/

  GMVectorSums vector_sums = GMVectorSums_null();
  GMVectorSums_create(&vector_sums, vectors, env);
  gm_compute_vector_sums(vectors, &vector_sums, env);

  /*---Numerator---*/
  const int numvec = metrics->num_vector_local;
  const int numfield = vectors->num_field;

  /*---Initialize MAGMA library---*/
  gm_magma_initialize(env);

  /*---Allocate magma CPU memory for vectors and for result---*/
  GMMirroredPointer vectors_buf = 
    gm_malloc_magma ( numvec * (size_t)numfield, env); //Data matrix X

  /*---Copy in vectors---*/
  gm_vectors_to_buf(vectors, &vectors_buf, env);

  /*---Send matrix vectors to GPU---*/
  
  gm_set_vectors_start(vectors, &vectors_buf, env);
  gm_set_vectors_wait(env);
  
  /*---Compute numerators---*/
  gm_compute_czekanowski_numerators_3way_start(
       vectors, vectors, vectors, metrics,
       &vectors_buf, &vectors_buf, &vectors_buf, 
       env->proc_num, GM_BOOL_TRUE, env);
  gm_compute_wait(env);

  /*---Combine results---*/
  gm_compute_czekanowski_3way_combine(metrics,
                     (GMFloat*)(&vector_sums)->data,
                     (GMFloat*)(&vector_sums)->data,
                     (GMFloat*)(&vector_sums)->data,
                     env->proc_num,
                     GM_BOOL_TRUE, env);

  /*---Free memory and finalize---*/

  GMVectorSums_destroy( &vector_sums, env);  

  gm_free_magma( &vectors_buf, env);
  
  gm_magma_finalize(env);

}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
