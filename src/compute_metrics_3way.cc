/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_3way.c
 * \author Wayne Joubert
 * \date   Thu Jan 21 19:07:47 EST 2016
 * \brief  Functions for computing 3-way metrics.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include "env.h"
#include "vector_sums.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_utils_magma.h"
#include "compute_metrics_utils.h"
#include "compute_metrics_3way.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

void gm_compute_metrics_3way_notall2all(GMMetrics* metrics,
                                        GMVectors* vectors,
                                        GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMAssert(!Env_all2all(env));

  /*---Denominator---*/

  GMVectorSums vector_sums = GMVectorSums_null();
  GMVectorSums_create(&vector_sums, vectors, env);
  GMVectorSums_compute(&vector_sums, vectors, env);

  /*---------------*/
  /*---Numerator---*/
  /*---------------*/

  gm_magma_initialize(env);

  const int numvecl = vectors->num_vector_local;
  const int numpfieldl = vectors->num_packedval_field_local;

  /*---Allocate magma CPU memory for vectors and for result */

  GMMirroredPointer vectors_buf =
      gm_malloc_magma(numvecl * (size_t)numpfieldl, env);

  GMMirroredPointer metrics_buf =
      gm_malloc_magma(numvecl * (size_t)numvecl, env);

  GMMirroredPointer metrics_buf_tmp =
      gm_malloc_magma(numvecl * (size_t)numvecl, env);

  GMMirroredPointer* metrics_buf_local =
      Env_num_proc_field(env) == 1 ? &metrics_buf : &metrics_buf_tmp;

  /*---Copy in vectors---*/

  gm_vectors_to_buf(vectors, &vectors_buf, env);

  /*---Send vectors to GPU---*/

  gm_set_vectors_start(vectors, &vectors_buf, env);
  gm_set_vectors_wait(env);

  /*---Compute numerators---*/

  gm_compute_numerators_3way_start(
      vectors, vectors, vectors, metrics, &vectors_buf, &vectors_buf,
      &vectors_buf, Env_proc_num_vector_i(env), Env_proc_num_vector_i(env),
      env);
  gm_compute_wait(env);

  /*---Copy result from GPU---*/

  gm_get_metrics_start(metrics, metrics_buf_local, env);
  gm_get_metrics_wait(metrics, metrics_buf_local, env);

  /*---Do reduction across field procs if needed---*/

  if (Env_num_proc_field(env) > 1) {
    gm_allreduce_metrics(metrics, &metrics_buf, metrics_buf_local, env);
  }

  /*---------------*/
  /*---Combine---*/
  /*---------------*/

  gm_compute_3way_combine(metrics, &vector_sums, &vector_sums, &vector_sums,
                          Env_proc_num_vector_i(env),
                          Env_proc_num_vector_i(env), env);

  /*---------------*/
  /*---Free memory---*/
  /*---------------*/

  GMVectorSums_destroy(&vector_sums, env);

  gm_free_magma(&vectors_buf, env);
  gm_free_magma(&metrics_buf, env);
  gm_free_magma(&metrics_buf_tmp, env);

  gm_magma_finalize(env);
}

/*===========================================================================*/

void gm_compute_metrics_3way_all2all(GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMAssert(Env_all2all(env));

  const int num_proc = Env_num_block_vector(env);

  /*---Initialize MAGMA library---*/

  gm_magma_initialize(env);

  /*---Initializations---*/

  const int numvecl = metrics->num_vector_local;
  const int numpfieldl = vectors->num_packedval_field_local;

  const int i_proc = Env_proc_num_vector_i(env);

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
      gm_malloc_magma(numvecl * (size_t)numpfieldl, env);

  /*---Copy in vectors---*/

  gm_vectors_to_buf(vectors_i, &vectors_i_buf, env);

  /*---Send matrix vectors to GPU---*/

  gm_set_vectors_start(vectors_i, &vectors_i_buf, env);
  gm_set_vectors_wait(env);

  /*---Compute numerators---*/

  gm_compute_numerators_3way_start(vectors_i, vectors_i, vectors_i, metrics,
                                   &vectors_i_buf, &vectors_i_buf,
                                   &vectors_i_buf, i_proc, i_proc, env);
  gm_compute_wait(env);

  /*---Combine results---*/

  gm_compute_3way_combine(metrics, &vector_sums_i, &vector_sums_i,
                          &vector_sums_i, i_proc, i_proc, env);

  /*------------------------*/
  /*---Part 2 Computation: triangular prisms---*/
  /*------------------------*/

  const int data_type = Env_data_type_vectors(env);

  GMVectors vectors_j_value = GMVectors_null();
  GMVectors* vectors_j = &vectors_j_value;
  GMVectors_create(vectors_j, data_type, vectors->num_field, numvecl, env);

  GMMirroredPointer vectors_j_buf =
      gm_malloc_magma(numvecl * (size_t)numpfieldl, env);

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

    gm_compute_numerators_3way_start(vectors_i, vectors_j, vectors_j, metrics,
                                     &vectors_i_buf, &vectors_j_buf,
                                     &vectors_j_buf, j_proc, j_proc, env);
    gm_compute_wait(env);

    /*---Combine results---*/

    gm_compute_3way_combine(metrics, &vector_sums_i, &vector_sums_j,
                            &vector_sums_j, j_proc, j_proc, env);

  } /*---proc_diff_j---*/

  /*------------------------*/
  /*---Part 3 Computation: block sections---*/
  /*------------------------*/

  GMVectors vectors_k_value = GMVectors_null();
  GMVectors* vectors_k = &vectors_k_value;
  GMVectors_create(vectors_k, data_type, vectors->num_field, numvecl, env);

  GMMirroredPointer vectors_k_buf =
      gm_malloc_magma(numvecl * (size_t)numpfieldl, env);

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

      gm_compute_numerators_3way_start(vectors_i, vectors_j, vectors_k, metrics,
                                       &vectors_i_buf, &vectors_j_buf,
                                       &vectors_k_buf, j_proc, k_proc, env);
      gm_compute_wait(env);

      /*---Combine results---*/

      gm_compute_3way_combine(metrics, &vector_sums_i, &vector_sums_j,
                              &vector_sums_k, j_proc, k_proc, env);

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

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
