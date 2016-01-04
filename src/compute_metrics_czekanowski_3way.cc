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
extern "C" {
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

  const int num_proc = Env_num_proc_vector(env);

  /*---Initialize MAGMA library---*/

  gm_magma_initialize(env);

  /*---Initializations---*/

  const int numvec = metrics->num_vector_local;
  const int numfield = vectors->num_field_local;

  // int i = 0;
  const int i_proc = Env_proc_num_vector(env);

  /*------------------------*/
  /*---Part 1 Computation: tetrahedron---*/
  /*------------------------*/

  /*---Denominator---*/

  GMVectors* vectors_i = vectors;

  GMVectorSums vector_sums_i = GMVectorSums_null();
  GMVectorSums_create(&vector_sums_i, vectors, env);
  GMVectorSums_compute(&vector_sums_i, vectors_i, env);

  /*---Numerator---*/

  /*---Allocate magma CPU memory for vectors and for result---*/

  GMMirroredPointer vectors_i_buf =
      gm_malloc_magma(numvec * (size_t)numfield, env);

  /*---Copy in vectors---*/

  gm_vectors_to_buf(vectors_i, &vectors_i_buf, env);

  /*---Send matrix vectors to GPU---*/

  gm_set_vectors_start(vectors_i, &vectors_i_buf, env);
  gm_set_vectors_wait(env);

  /*---Compute numerators---*/

  gm_compute_czekanowski_numerators_3way_start(
      vectors_i, vectors_i, vectors_i, metrics, &vectors_i_buf, &vectors_i_buf,
      &vectors_i_buf, i_proc, i_proc, env);
  gm_compute_wait(env);

  /*---Combine results---*/

  gm_compute_czekanowski_3way_combine(metrics, (GMFloat*)(&vector_sums_i)->data,
                                      (GMFloat*)(&vector_sums_i)->data,
                                      (GMFloat*)(&vector_sums_i)->data, i_proc,
                                      i_proc, env);

  /*------------------------*/
  /*---Part 2 Computation: triangular prisms---*/
  /*------------------------*/

  const int data_type = Env_data_type_vectors(env);

  GMVectors vectors_j_value = GMVectors_null();
  GMVectors* vectors_j = &vectors_j_value;
  GMVectors_create(vectors_j, data_type, vectors->num_field, numvec, env);

  GMMirroredPointer vectors_j_buf =
      gm_malloc_magma(numvec * (size_t)numfield, env);

  GMVectorSums vector_sums_j = GMVectorSums_null();
  GMVectorSums_create(&vector_sums_j, vectors, env);

  int proc_diff_j = 0;
  for (proc_diff_j = 1; proc_diff_j < num_proc; ++proc_diff_j) {
    MPI_Request mpi_requests_j[2];

    const int proc_send_j = (i_proc - proc_diff_j + num_proc) % num_proc;
    const int proc_recv_j = (i_proc + proc_diff_j) % num_proc;
    const int j_proc = proc_recv_j;

    mpi_requests_j[0] = gm_send_vectors_start(vectors_i, proc_send_j, env);
    mpi_requests_j[1] = gm_recv_vectors_start(vectors_j, proc_recv_j, env);

    gm_send_vectors_wait(&(mpi_requests_j[0]), env);
    gm_recv_vectors_wait(&(mpi_requests_j[1]), env);

    gm_vectors_to_buf(vectors_j, &vectors_j_buf, env);

    gm_set_vectors_start(vectors_j, &vectors_j_buf, env);
    gm_set_vectors_wait(env);

    GMVectorSums_compute(&vector_sums_j, vectors_j, env);

    /*---Compute numerators---*/

    gm_compute_czekanowski_numerators_3way_start(
        vectors_i, vectors_j, vectors_j, metrics, &vectors_i_buf,
        &vectors_j_buf, &vectors_j_buf, j_proc, j_proc, env);
    gm_compute_wait(env);

    /*---Combine results---*/

    gm_compute_czekanowski_3way_combine(
        metrics, (GMFloat*)(&vector_sums_i)->data,
        (GMFloat*)(&vector_sums_j)->data, (GMFloat*)(&vector_sums_j)->data,
        j_proc, j_proc, env);

  } /*---proc_diff_j---*/

  /*------------------------*/
  /*---Part 3 Computation: block sections---*/
  /*------------------------*/

  GMVectors vectors_k_value = GMVectors_null();
  GMVectors* vectors_k = &vectors_k_value;
  GMVectors_create(vectors_k, data_type, vectors->num_field, numvec, env);

  GMMirroredPointer vectors_k_buf =
      gm_malloc_magma(numvec * (size_t)numfield, env);

  GMVectorSums vector_sums_k = GMVectorSums_null();
  GMVectorSums_create(&vector_sums_k, vectors, env);

  for (proc_diff_j = 1; proc_diff_j < num_proc; ++proc_diff_j) {
    const int proc_send_j = (i_proc - proc_diff_j + num_proc) % num_proc;
    const int proc_recv_j = (i_proc + proc_diff_j) % num_proc;
    const int j_proc = proc_recv_j;

    MPI_Request mpi_requests_j[2];

    mpi_requests_j[0] = gm_send_vectors_start(vectors_i, proc_send_j, env);
    mpi_requests_j[1] = gm_recv_vectors_start(vectors_j, proc_recv_j, env);

    gm_send_vectors_wait(&(mpi_requests_j[0]), env);
    gm_recv_vectors_wait(&(mpi_requests_j[1]), env);

    gm_vectors_to_buf(vectors_j, &vectors_j_buf, env);

    gm_set_vectors_start(vectors_j, &vectors_j_buf, env);
    gm_set_vectors_wait(env);

    GMVectorSums_compute(&vector_sums_j, vectors_j, env);

    int proc_diff_k = 0;

    for (proc_diff_k = 1; proc_diff_k < num_proc; ++proc_diff_k) {
      const int proc_send_k = (i_proc - proc_diff_k + num_proc) % num_proc;
      const int proc_recv_k = (i_proc + proc_diff_k) % num_proc;
      const int k_proc = proc_recv_k;

      MPI_Request mpi_requests_k[2];

      mpi_requests_k[0] = gm_send_vectors_start(vectors_i, proc_send_k, env);
      mpi_requests_k[1] = gm_recv_vectors_start(vectors_k, proc_recv_k, env);

      gm_send_vectors_wait(&(mpi_requests_k[0]), env);
      gm_recv_vectors_wait(&(mpi_requests_k[1]), env);

      //---TODO: skip the above communication if not needed.

      if (k_proc == j_proc) {
        continue;
      }

      gm_vectors_to_buf(vectors_k, &vectors_k_buf, env);

      gm_set_vectors_start(vectors_k, &vectors_k_buf, env);
      gm_set_vectors_wait(env);

      GMVectorSums_compute(&vector_sums_k, vectors_k, env);

      /*---Compute numerators---*/

      gm_compute_czekanowski_numerators_3way_start(
          vectors_i, vectors_j, vectors_k, metrics, &vectors_i_buf,
          &vectors_j_buf, &vectors_k_buf, j_proc, k_proc, env);
      gm_compute_wait(env);

      /*---Combine results---*/

      gm_compute_czekanowski_3way_combine(
          metrics, (GMFloat*)(&vector_sums_i)->data,
          (GMFloat*)(&vector_sums_j)->data, (GMFloat*)(&vector_sums_k)->data,
          j_proc, k_proc, env);

    } /*---proc_diff_k---*/
  }   /*---proc_diff_j---*/

  /*---Free memory and finalize---*/

  GMVectors_destroy(vectors_k, env);
  GMVectors_destroy(vectors_j, env);

  GMVectorSums_destroy(&vector_sums_k, env);
  GMVectorSums_destroy(&vector_sums_j, env);
  GMVectorSums_destroy(&vector_sums_i, env);

  gm_free_magma(&vectors_k_buf, env);
  gm_free_magma(&vectors_j_buf, env);
  gm_free_magma(&vectors_i_buf, env);

  gm_magma_finalize(env);
}

/*===========================================================================*/

void gm_compute_metrics_czekanowski_3way_cpu(GMMetrics* metrics,
                                             GMVectors* vectors,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  if (Env_all2all(env)) {
    gm_compute_metrics_czekanowski_3way_all2all(metrics, vectors, env);
    return;
  }

  GMInsist(env, Env_num_proc_field(env) == 1
    ? "num_proc_field>1 for CPU case not supported" : 0);

  /*---Denominator---*/

  GMFloat* vector_sums = GMFloat_malloc(metrics->num_vector_local);
  GMFloat* vector_sums_tmp = GMFloat_malloc(metrics->num_vector_local);

  gm_compute_float_vector_sums(vectors, vector_sums, vector_sums_tmp, env);

  /*---Numerator---*/

  int i = 0;
  int j = 0;
  int k = 0;
  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      for (k = j + 1; k < metrics->num_vector_local; ++k) {
        GMFloat sum = 0;
        int field_local = 0;
        for (field_local = 0; field_local < vectors->num_field_local;
             ++field_local) {
          const GMFloat value1 = GMVectors_float_get(vectors,
                                                     field_local, i, env);
          const GMFloat value2 = GMVectors_float_get(vectors,
                                                     field_local, j, env);
          const GMFloat value3 = GMVectors_float_get(vectors,
                                                     field_local, k, env);
          GMFloat min12 = value1 < value2 ? value1 : value2;
          sum += min12;
          sum += value1 < value3 ? value1 : value3;
          sum += value2 < value3 ? value2 : value3;
          sum -= min12 < value3 ? min12 : value3;
        } /*---for field_local---*/
        GMMetrics_float_set_3(metrics, i, j, k, sum, env);
      } /*---for k---*/
    }   /*---for j---*/
  }     /*---for i---*/

  /*---Combine---*/
  // printf("---CPU Implementation---\n");
  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      for (k = j + 1; k < metrics->num_vector_local; ++k) {
        const GMFloat numerator = GMMetrics_float_get_3(metrics, i, j, k, env);
        const GMFloat denominator =
            vector_sums[i] + vector_sums[j] + vector_sums[k];
        // printf("%i,%i,%i . . . numerator = %f . . . denominator =
        // %f\n",i,j,k,numerator,denominator);
        GMMetrics_float_set_3(metrics, i, j, k,
                              3 * numerator / (2 * denominator), env);
      } /*---for k---*/
    }   /*---for j---*/
  }     /*---for i---*/

  free(vector_sums);
  free(vector_sums_tmp);
}

/*===========================================================================*/

void gm_compute_metrics_czekanowski_3way_gpu(GMMetrics* metrics,
                                             GMVectors* vectors,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  if (Env_all2all(env)) {
    gm_compute_metrics_czekanowski_3way_all2all(metrics, vectors, env);
    return;
  }

  /*---Denominator---*/

  GMVectorSums vector_sums = GMVectorSums_null();
  GMVectorSums_create(&vector_sums, vectors, env);
  GMVectorSums_compute(&vector_sums, vectors, env);

  /*---Numerator---*/
  const int numvec = metrics->num_vector_local;
  const int numfield = vectors->num_field_local;

  /*---Initialize MAGMA library---*/
  gm_magma_initialize(env);

  /*---Allocate magma CPU memory for vectors and for result---*/
  GMMirroredPointer vectors_buf =
      gm_malloc_magma(numvec * (size_t)numfield, env);  // Data matrix X

  /*---Copy in vectors---*/
  gm_vectors_to_buf(vectors, &vectors_buf, env);

  /*---Send matrix vectors to GPU---*/

  gm_set_vectors_start(vectors, &vectors_buf, env);
  gm_set_vectors_wait(env);

  /*---Compute numerators---*/
  gm_compute_czekanowski_numerators_3way_start(
      vectors, vectors, vectors, metrics, &vectors_buf, &vectors_buf,
      &vectors_buf, Env_proc_num_vector(env), Env_proc_num_vector(env), env);
  gm_compute_wait(env);

  /*---Combine results---*/
  gm_compute_czekanowski_3way_combine(metrics,
      (GMFloat*)(&vector_sums)->data, (GMFloat*)(&vector_sums)->data,
      (GMFloat*)(&vector_sums)->data, Env_proc_num_vector(env),
      Env_proc_num_vector(env), env);

  /*---Free memory and finalize---*/

  GMVectorSums_destroy(&vector_sums, env);

  gm_free_magma(&vectors_buf, env);

  gm_magma_finalize(env);
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
