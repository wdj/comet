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

  const int section_step = 0;
  GMAssertAlways(GMEnv_num_section_steps(env, 1) == 1);

  gm_compute_numerators_3way_start(
      vectors, vectors, vectors, metrics, &vectors_buf, &vectors_buf,
      &vectors_buf, Env_proc_num_vector_i(env), Env_proc_num_vector_i(env),
      &vector_sums, &vector_sums, &vector_sums,
      section_step, env);
  gm_compute_wait(env);

  /*---Copy result from GPU---*/

  gm_get_metrics_start(metrics, metrics_buf_local, env);
  gm_get_metrics_wait(metrics, metrics_buf_local, env);

  /*---Do reduction across field procs if needed---*/

  if (Env_num_proc_field(env) > 1) {
    gm_allreduce_metrics(metrics, &metrics_buf, metrics_buf_local, env);
  }

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

  const int numvecl = metrics->num_vector_local;
  const int numpfieldl = vectors->num_packedval_field_local;

  const int data_type = Env_data_type_vectors(env);

  const int num_block = Env_num_block_vector(env);

  const int i_block = Env_proc_num_vector_i(env);

  const int proc_num_r = Env_proc_num_repl(env);
  const int num_proc_r = Env_num_proc_repl(env);

  /*---Create flattened index within space of procs assigned to
       vectors (non-field procs) - i.e., vector_i (=block) X repl ---*/

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

  int section_block_num = 0;

  /*------------------------*/
  /*---Part 1 Computation: tetrahedron---*/
  /*------------------------*/

//vvv
  /*---Denominator---*/
  GMVectorSums_compute(&vector_sums_i, vectors_i, env);

  /*---Copy in vectors---*/
  gm_vectors_to_buf(vectors_i, &vectors_i_buf, env);

  /*---Send vectors to GPU---*/
  gm_set_vectors_start(vectors_i, &vectors_i_buf, env);
  gm_set_vectors_wait(env);
//^^^

  int section_step = 0;

  for (section_step=0; section_step<GMEnv_num_section_steps(env, 1);
       ++section_step) {

    if (section_block_num % num_proc_r == proc_num_r) {

//vvv
      /*---Compute numerators---*/
      gm_compute_numerators_3way_start(vectors_i, vectors_i, vectors_i, metrics,
                                       &vectors_i_buf, &vectors_i_buf,
                                       &vectors_i_buf, i_block, i_block,
                                       &vector_sums_i, &vector_sums_i,
                                       &vector_sums_i,
                                       section_step, env);
      gm_compute_wait(env);
//^^^

    } /*---if (section_block_num ...)---*/

    ++section_block_num;
  } /*---section_step---*/

  /*------------------------*/
  /*---Part 2 Computation: triangular prisms---*/
  /*------------------------*/

#if 0

Possible plan for setting up async comm

- each loop kept as-is
- each loop post communications for its step
- after post of send/recv, call compute numerators start for
  any pending prev step, if any.
- after this, post wait for completion of this step send/recv

- after final of the three loops, final cleanup call to do any pending
  compute numerators

- will need to set up double buffers and indexing for various things:
  - mpi requests
  - vectors_j/k
  - _Bool is_pending_compute_numerators
  - ind, ind_prev, j_block_prev, k_block_prev, section_step_prev

- do we need to push some k-loop calculations down into the j-loop

#endif

  MPI_Request mpi_request_send_j;
  MPI_Request mpi_request_recv_j;

  for (section_step=0; section_step<GMEnv_num_section_steps(env, 2);
       ++section_step) {

    int j_i_block_delta = 0;
    for (j_i_block_delta = 1; j_i_block_delta < num_block; ++j_i_block_delta) {

      const int j_block = gm_mod_i(i_block + j_i_block_delta, num_block);

      const int proc_send_j = gm_mod_i(proc_num_ir - j_i_block_delta*num_proc_r,
                                       num_proc_ir);
      const int proc_recv_j = gm_mod_i(proc_num_ir + j_i_block_delta*num_proc_r,
                                       num_proc_ir);

      if (section_block_num % num_proc_r == proc_num_r) {

        /*---Communicate vectors start---*/
        mpi_request_send_j = gm_send_vectors_start(vectors_i, proc_send_j, env);
        mpi_request_recv_j = gm_recv_vectors_start(vectors_j, proc_recv_j, env);

        /*---Communicate vectors wait---*/
        gm_send_vectors_wait(&mpi_request_send_j, env);
        gm_recv_vectors_wait(&mpi_request_recv_j, env);

//vvv
        /*---Copy in vectors---*/
        gm_vectors_to_buf(vectors_j, &vectors_j_buf, env);

        /*---Send vectors to GPU---*/
        gm_set_vectors_start(vectors_j, &vectors_j_buf, env);
        gm_set_vectors_wait(env);

        /*---Denominator---*/
        GMVectorSums_compute(&vector_sums_j, vectors_j, env);

        /*---Compute numerators---*/
        gm_compute_numerators_3way_start(vectors_i, vectors_j, vectors_j,
                                         metrics, &vectors_i_buf,
                                         &vectors_j_buf, &vectors_j_buf,
                                         j_block, j_block,
                                         &vector_sums_i, &vector_sums_j,
                                         &vector_sums_j,
                                         section_step, env);
        gm_compute_wait(env);
//^^^

      } /*---if (section_block_num ...)---*/

      ++section_block_num;
    } /*---j_i_block_delta---*/
  } /*---section_step---*/

  /*------------------------*/
  /*---Part 3 Computation: block sections---*/
  /*------------------------*/

  MPI_Request mpi_request_send_k;
  MPI_Request mpi_request_recv_k;

  for (section_step=0; section_step<GMEnv_num_section_steps(env, 3);
       ++section_step) {

    int k_i_block_delta = 0;
    for (k_i_block_delta = 1; k_i_block_delta < num_block; ++k_i_block_delta) {

      const int k_block = gm_mod_i(i_block + k_i_block_delta, num_block);

      const int proc_send_k = gm_mod_i(proc_num_ir - k_i_block_delta*num_proc_r,
                                       num_proc_ir);
      const int proc_recv_k = gm_mod_i(proc_num_ir + k_i_block_delta*num_proc_r,
                                       num_proc_ir);

      /*---Communicate vectors start---*/
      mpi_request_send_k = gm_send_vectors_start(vectors_i, proc_send_k, env);
      mpi_request_recv_k = gm_recv_vectors_start(vectors_k, proc_recv_k, env);

      /*---Communicate vectors wait---*/
      gm_send_vectors_wait(&mpi_request_send_k, env);
      gm_recv_vectors_wait(&mpi_request_recv_k, env);

//vvv
      /*---Copy in vectors---*/
      gm_vectors_to_buf(vectors_k, &vectors_k_buf, env);

      /*---Send vectors to GPU---*/
      gm_set_vectors_start(vectors_k, &vectors_k_buf, env);
      gm_set_vectors_wait(env);

      /*---Denominator---*/
      GMVectorSums_compute(&vector_sums_k, vectors_k, env);
//^^^

      int j_i_block_delta = 0;
      for (j_i_block_delta = 1; j_i_block_delta < num_block; ++j_i_block_delta){

        const int j_block = gm_mod_i(i_block + j_i_block_delta, num_block);

        const int proc_send_j = gm_mod_i(proc_num_ir-j_i_block_delta*num_proc_r,
                                         num_proc_ir);
        const int proc_recv_j = gm_mod_i(proc_num_ir+j_i_block_delta*num_proc_r,
                                         num_proc_ir);
        if (j_block == k_block) {
          /*---NOTE: this condition occurs on all procs at exactly the same
               j/k iteration in lockstep, so there is no chance the immediately
               following communication will deadlock/mispair---*/
          continue;
        }
        GMAssertAlways((j_block == k_block) ==
                       (j_i_block_delta == k_i_block_delta));

        if (section_block_num % num_proc_r == proc_num_r) {

          /*---Communicate vectors start---*/
          mpi_request_send_j = gm_send_vectors_start(vectors_i, proc_send_j,
                                                     env);
          mpi_request_recv_j = gm_recv_vectors_start(vectors_j, proc_recv_j,
                                                     env);

          /*---Communicate vectors wait---*/
          gm_send_vectors_wait(&mpi_request_send_j, env);
          gm_recv_vectors_wait(&mpi_request_recv_j, env);

//vvv
          /*---Copy in vectors---*/
          gm_vectors_to_buf(vectors_j, &vectors_j_buf, env);

          /*---Send vectors to GPU---*/
          gm_set_vectors_start(vectors_j, &vectors_j_buf, env);
          gm_set_vectors_wait(env);

          /*---Denominator---*/
          GMVectorSums_compute(&vector_sums_j, vectors_j, env);

          /*---Compute numerators---*/
          gm_compute_numerators_3way_start(vectors_i, vectors_j, vectors_k,
                                           metrics, &vectors_i_buf,
                                           &vectors_j_buf, &vectors_k_buf,
                                           j_block, k_block,
                                           &vector_sums_i, &vector_sums_j,
                                           &vector_sums_k,
                                           section_step, env);
          gm_compute_wait(env);
//^^^

        } /*---if (section_block_num ...)---*/

        ++section_block_num;
      } /*---k_i_block_delta---*/
    }   /*---j_i_block_delta---*/
  } /*---section_step---*/

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
