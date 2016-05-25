/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_czekanowski_2way.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Functions for computing 2-way Czekanowski metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>

#include "env.h"
#include "vector_sums.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_utils_magma.h"
#include "compute_metrics_utils.h"
#include "compute_metrics_czekanowski_2way.h"

#ifdef __cplusplus
extern "C" {
#endif

#if 0
/*===========================================================================*/

void gm_compute_metrics_czekanowski_2way_cpu(GMMetrics* metrics,
                                             GMVectors* vectors,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, Env_num_proc_field(env) == 1
    ? "num_proc_field>1 for CPU case not supported" : 0);

  GMAssert(!Env_all2all(env));

  /*---Compute sums for denominators---*/

  GMFloat* vector_sums = GMFloat_malloc(metrics->num_vector_local);
  GMFloat* vector_sums_tmp = GMFloat_malloc(metrics->num_vector_local);

  gm_compute_float_vector_sums(vectors, vector_sums, vector_sums_tmp, env);

  /*---Compute numerators---*/

  int i = 0;
  int j = 0;
  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      GMFloat sum = 0;
      int field_local = 0;
      for (field_local = 0; field_local < vectors->num_field_local;
           ++field_local) {
        const GMFloat value_i = GMVectors_float_get(vectors,
                                                    field_local, i, env);
        const GMFloat value_j = GMVectors_float_get(vectors,
                                                    field_local, j, env);
        sum += value_i < value_j ? value_i : value_j;
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

  /*---Deallocations---*/

  free(vector_sums);
  free(vector_sums_tmp);
}

/*===========================================================================*/

void gm_compute_metrics_czekanowski_2way_gpu(GMMetrics* metrics,
                                             GMVectors* vectors,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMAssert(!Env_all2all(env));

  /*---------------*/
  /*---Denominator---*/
  /*---------------*/

  GMVectorSums vector_sums = GMVectorSums_null();
  GMVectorSums_create(&vector_sums, vectors, env);

  /* .02 / 1.56 */
  GMVectorSums_compute(&vector_sums, vectors, env);

  /*---------------*/
  /*---Numerator---*/
  /*---------------*/

  gm_magma_initialize(env);

  const int numvecl = metrics->num_vector_local;
  const int numpfieldl = vectors->num_field_local;

  /*---Allocate magma CPU memory for vectors and for result---*/

  GMMirroredPointer vectors_buf =
      gm_malloc_magma(numvecl * (size_t)numpfieldl, env);

  GMMirroredPointer metrics_buf =
      gm_malloc_magma(numvecl * (size_t)numvecl, env);

  GMMirroredPointer metrics_buf_tmp =
      gm_malloc_magma(numvecl * (size_t)numvecl, env);

  GMMirroredPointer* metrics_buf_local = Env_num_proc_field(env) == 1 ?
    &metrics_buf : &metrics_buf_tmp;

  /*---Copy in vectors---*/

  /* .08 / 1.56 */
  gm_vectors_to_buf(vectors, &vectors_buf, env);

  /*---Send vectors to GPU---*/

  gm_set_vectors_start(vectors, &vectors_buf, env);
  gm_set_vectors_wait(env);

  /*---Compute numerators---*/

  gm_compute_numerators_2way_start(vectors, vectors, metrics, &vectors_buf,
                                   &vectors_buf, metrics_buf_local,
                                   Env_proc_num_vector_i(env),
                                   GM_BOOL_TRUE, env);
  gm_compute_wait(env);

  /*---Copy result from GPU---*/

  gm_get_metrics_start(metrics, metrics_buf_local, env);
  gm_get_metrics_wait(metrics, metrics_buf_local, env);

  if (Env_num_proc_field(env) > 1) {
    int mpi_code = 0;
    mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
    mpi_code = MPI_Allreduce(metrics_buf_local->h, metrics_buf.h,
                 numvecl*(size_t)numvecl, GM_MPI_FLOAT, MPI_SUM,
                 Env_mpi_comm_field(env));
    GMAssert(mpi_code == MPI_SUCCESS);
  }

  /*---Combine---*/

  /* .22 / 1.56 */
  gm_compute_2way_combine(metrics, &metrics_buf, &vector_sums, &vector_sums,
                          Env_proc_num_vector_i(env), GM_BOOL_TRUE, env);

  /*---------------*/
  /*---Free memory and finalize---*/
  /*---------------*/

  GMVectorSums_destroy(&vector_sums, env);

  gm_free_magma(&vectors_buf, env);
  gm_free_magma(&metrics_buf, env);
  gm_free_magma(&metrics_buf_tmp, env);

  gm_magma_finalize(env);
}

/*---------------------------------------------------------------------------*/
#endif

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
