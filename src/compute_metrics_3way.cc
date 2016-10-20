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
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);

  GMAssertAlways(!Env_all2all(env));

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
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(Env_all2all(env));

  /*---Initializations---*/

  gm_magma_initialize(env);

  const int num_block = Env_num_block_vector(env);

  const int numvecl = metrics->num_vector_local;
  const int numpfieldl = vectors->num_packedval_field_local;

  const int i_block = Env_proc_num_vector_i(env);

  const int data_type = Env_data_type_vectors(env);

  const int proc_num_r = Env_proc_num_repl(env);
  const int num_proc_r = Env_num_proc_repl(env);

  const int proc_num_ir = proc_num_r + num_proc_r * i_block;
  const int num_proc_ir = num_block * num_proc_r;

  /*------------------------*/
  /*---Allocations: Part 1---*/
  /*------------------------*/

  GMVectors* vectors_i = vectors;

  GMVectorSums vector_sums_i = GMVectorSums_null();
  GMVectorSums_create(&vector_sums_i, vectors_i, env);

  GMMirroredPointer vectors_i_buf =
      gm_malloc_magma(numvecl * (size_t)numpfieldl, env);

  /*------------------------*/
  /*---Allocations: Part 2---*/
  /*------------------------*/

  GMVectors vectors_j_value = GMVectors_null();
  GMVectors* vectors_j = &vectors_j_value;
  GMVectors_create(vectors_j, data_type, vectors->num_field, numvecl, env);

  GMMirroredPointer vectors_j_buf =
      gm_malloc_magma(numvecl * (size_t)numpfieldl, env);

  GMVectorSums vector_sums_j = GMVectorSums_null();
  GMVectorSums_create(&vector_sums_j, vectors, env);

  /*------------------------*/
  /*---Allocations: Part 3---*/
  /*------------------------*/

  GMVectors vectors_k_value = GMVectors_null();
  GMVectors* vectors_k = &vectors_k_value;
  GMVectors_create(vectors_k, data_type, vectors->num_field, numvecl, env);

  GMMirroredPointer vectors_k_buf =
      gm_malloc_magma(numvecl * (size_t)numpfieldl, env);

  GMVectorSums vector_sums_k = GMVectorSums_null();
  GMVectorSums_create(&vector_sums_k, vectors, env);

  int block_num = 0;

  /*------------------------*/
  /*---Part 1 Computation: tetrahedron---*/
  /*------------------------*/

  /*---Denominator---*/

  GMVectorSums_compute(&vector_sums_i, vectors_i, env);

  /*---Copy in vectors---*/

  gm_vectors_to_buf(vectors_i, &vectors_i_buf, env);

  /*---Send vectors to GPU---*/

  gm_set_vectors_start(vectors_i, &vectors_i_buf, env);
  gm_set_vectors_wait(env);

  if (block_num % num_proc_r == proc_num_r) {

    /*---Compute numerators---*/

    gm_compute_numerators_3way_start(vectors_i, vectors_i, vectors_i, metrics,
                                     &vectors_i_buf, &vectors_i_buf,
                                     &vectors_i_buf, i_block, i_block, env);
    gm_compute_wait(env);

    /*---Combine results---*/

    gm_compute_3way_combine(metrics, &vector_sums_i, &vector_sums_i,
                            &vector_sums_i, i_block, i_block, env);

  } /*---if (block_num ...)---*/

  block_num++;

  /*------------------------*/
  /*---Part 2 Computation: triangular prisms---*/
  /*------------------------*/

  int j_i_block_delta = 0;
  for (j_i_block_delta = 1; j_i_block_delta < num_block; ++j_i_block_delta) {

    const int j_block = gm_mod_i(i_block + j_i_block_delta, num_block);

    const int proc_send_j = gm_mod_i(proc_num_ir - j_i_block_delta*num_proc_r,
                                     num_proc_ir);
    const int proc_recv_j = gm_mod_i(proc_num_ir + j_i_block_delta*num_proc_r,
                                     num_proc_ir);

    if (block_num % num_proc_r == proc_num_r) {

      /*---Communicate vectors---*/

      MPI_Request mpi_requests_j[2];

      mpi_requests_j[0] = gm_send_vectors_start(vectors_i, proc_send_j, env);
      mpi_requests_j[1] = gm_recv_vectors_start(vectors_j, proc_recv_j, env);

      gm_send_vectors_wait(&(mpi_requests_j[0]), env);
      gm_recv_vectors_wait(&(mpi_requests_j[1]), env);

      /*---Copy in vectors---*/

      gm_vectors_to_buf(vectors_j, &vectors_j_buf, env);

      /*---Send vectors to GPU---*/

      gm_set_vectors_start(vectors_j, &vectors_j_buf, env);
      gm_set_vectors_wait(env);

      /*---Compute numerators---*/

      gm_compute_numerators_3way_start(vectors_i, vectors_j, vectors_j, metrics,
                                       &vectors_i_buf, &vectors_j_buf,
                                       &vectors_j_buf, j_block, j_block, env);
      gm_compute_wait(env);

      /*---Denominator---*/

      GMVectorSums_compute(&vector_sums_j, vectors_j, env);

      /*---Combine results---*/

      gm_compute_3way_combine(metrics, &vector_sums_i, &vector_sums_j,
                              &vector_sums_j, j_block, j_block, env);

    } /*---if (block_num ...)---*/

    block_num++;
  } /*---j_i_block_delta---*/

  /*------------------------*/
  /*---Part 3 Computation: block sections---*/
  /*------------------------*/

  int k_i_block_delta = 0;
  for (k_i_block_delta = 1; k_i_block_delta < num_block; ++k_i_block_delta) {

    const int k_block = gm_mod_i(i_block + k_i_block_delta, num_block);

    const int proc_send_k = gm_mod_i(proc_num_ir - k_i_block_delta*num_proc_r,
                                     num_proc_ir);
    const int proc_recv_k = gm_mod_i(proc_num_ir + k_i_block_delta*num_proc_r,
                                     num_proc_ir);
    /*---Communicate vectors---*/

    MPI_Request mpi_requests_k[2];

    mpi_requests_k[0] = gm_send_vectors_start(vectors_i, proc_send_k, env);
    mpi_requests_k[1] = gm_recv_vectors_start(vectors_k, proc_recv_k, env);

    gm_send_vectors_wait(&(mpi_requests_k[0]), env);
    gm_recv_vectors_wait(&(mpi_requests_k[1]), env);

    /*---Copy in vectors---*/

    gm_vectors_to_buf(vectors_k, &vectors_k_buf, env);

    /*---Send vectors to GPU---*/

    gm_set_vectors_start(vectors_k, &vectors_k_buf, env);
    gm_set_vectors_wait(env);

    /*---Denominator---*/

    GMVectorSums_compute(&vector_sums_k, vectors_k, env);

    for (j_i_block_delta = 1; j_i_block_delta < num_block; ++j_i_block_delta) {

      const int j_block = gm_mod_i(i_block + j_i_block_delta, num_block);

      const int proc_send_j = gm_mod_i(proc_num_ir - j_i_block_delta*num_proc_r,
                                       num_proc_ir);
      const int proc_recv_j = gm_mod_i(proc_num_ir + j_i_block_delta*num_proc_r,
                                       num_proc_ir);
      if (j_block == k_block) {
        /*---NOTE: this condition occurs on all procs at exactly the same
             j/k iteration in lockstep, so there is no chance the immediately
             following communication will deadlock/mispair---*/
        continue;
      }
      GMAssertAlways((j_block == k_block) == (j_i_block_delta == k_i_block_delta));

      if (block_num % num_proc_r == proc_num_r) {

#ifdef GM_ASSERTIONS_ON
        const int block_num_calculated =
          (num_block) +
          ((num_block-2) * (k_i_block_delta - 1)) +
          (j_i_block_delta - 1 - (j_i_block_delta > k_i_block_delta));
        GMAssert(block_num_calculated == block_num);
#endif

        /*---Communicate vectors---*/

        MPI_Request mpi_requests_j[2];

        mpi_requests_j[0] = gm_send_vectors_start(vectors_i, proc_send_j, env);
        mpi_requests_j[1] = gm_recv_vectors_start(vectors_j, proc_recv_j, env);

        gm_send_vectors_wait(&(mpi_requests_j[0]), env);
        gm_recv_vectors_wait(&(mpi_requests_j[1]), env);

        /*---Copy in vectors---*/

        gm_vectors_to_buf(vectors_j, &vectors_j_buf, env);

        /*---Send vectors to GPU---*/

        gm_set_vectors_start(vectors_j, &vectors_j_buf, env);
        gm_set_vectors_wait(env);

        /*---Compute numerators---*/

        gm_compute_numerators_3way_start(vectors_i, vectors_j, vectors_k,
                                         metrics,
                                         &vectors_i_buf, &vectors_j_buf,
                                         &vectors_k_buf, j_block, k_block, env);
        gm_compute_wait(env);

        /*---Denominator---*/

        GMVectorSums_compute(&vector_sums_j, vectors_j, env);

        /*---Combine results---*/

        gm_compute_3way_combine(metrics, &vector_sums_i, &vector_sums_j,
                                &vector_sums_k, j_block, k_block, env);

      } /*---if (block_num ...)---*/

      block_num++;
    } /*---k_i_block_delta---*/
  }   /*---j_i_block_delta---*/

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
