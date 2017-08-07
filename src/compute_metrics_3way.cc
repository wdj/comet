/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_3way.cc
 * \author Wayne Joubert
 * \date   Thu Jan 21 19:07:47 EST 2016
 * \brief  Functions for computing 3-way metrics.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include "env.hh"
#include "mirrored_buf.hh"
#include "vector_sums.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "linalg.hh"
#include "compute_metrics_utils.hh"
#include "compute_metrics_utils_3way.hh"
#include "compute_metrics_3way.hh"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

void gm_compute_metrics_3way_notall2all(GMMetrics* metrics,
                                        GMVectors* vectors,
                                        GMEnv* env) {
  GMAssertAlways(metrics && vectors && env);
  GMAssertAlways(!GMEnv_all2all(env));

  /*---Denominator---*/

  GMVectorSums vector_sums = GMVectorSums_null();
  GMVectorSums_create(&vector_sums, vectors, env);
  GMVectorSums_compute(&vector_sums, vectors, env);

  /*---------------*/
  /*---Numerator---*/
  /*---------------*/

  gm_linalg_initialize(env);

  const int nvl = vectors->num_vector_local;
  const int npvfl = vectors->num_packedval_field_local;

  /*---Allocate magma CPU memory for vectors and for result */

  GMMirroredBuf vectors_buf = GMMirroredBuf_null();
  GMMirroredBuf_create(&vectors_buf, npvfl, nvl, env);

  /*---Copy in vectors---*/

  gm_vectors_to_buf(&vectors_buf, vectors, env);

  /*---Send vectors to GPU---*/

  gm_set_vectors_start(vectors, &vectors_buf, env);
  gm_set_vectors_wait(env);

  /*---Compute numerators---*/

  const int section_step = 0;
  GMAssertAlways(gm_num_section_steps(env, 1) == 1);

  GMComputeNumerators3Way gm_compute_numerators_3way = {0};
  GMComputeNumerators3Way_create(&gm_compute_numerators_3way, nvl, npvfl, env);

  GMComputeNumerators3Way_start(
      &gm_compute_numerators_3way,
      vectors, vectors, vectors, metrics, &vectors_buf, &vectors_buf,
      &vectors_buf, GMEnv_proc_num_vector_i(env), GMEnv_proc_num_vector_i(env),
      &vector_sums, &vector_sums, &vector_sums,
      section_step, env);
  gm_compute_wait(env);

  /*---------------*/
  /*---Free memory---*/
  /*---------------*/

  GMComputeNumerators3Way_destroy(&gm_compute_numerators_3way, env);

  GMVectorSums_destroy(&vector_sums, env);

  GMMirroredBuf_destroy(&vectors_buf, env);

  gm_linalg_finalize(env);
}

/*===========================================================================*/

void gm_compute_metrics_3way_all2all(GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env) {
  GMAssertAlways(metrics && vectors && env);
  GMAssertAlways(GMEnv_all2all(env));

  /*---Initializations---*/

  gm_linalg_initialize(env);

  const int nvl = metrics->num_vector_local;
  const int npvfl = vectors->num_packedval_field_local;

  const int data_type = GMEnv_data_type_vectors(env);

  const int num_block = GMEnv_num_block_vector(env);

  const int i_block = GMEnv_proc_num_vector_i(env);

  const int proc_num_r = GMEnv_proc_num_repl(env);
  const int num_proc_r = GMEnv_num_proc_repl(env);

  /*---Create flattened index within space of procs assigned to
       vectors (non-field procs) - i.e., vector_i (=block) X repl ---*/

  const int proc_num_ir = proc_num_r + num_proc_r * i_block;
  const int num_proc_ir = num_block * num_proc_r;

  /*------------------------*/
  /*---Allocations: Part 1---*/
  /*------------------------*/

  GMVectorSums vector_sums_i_value = GMVectorSums_null();
  GMVectorSums* const vector_sums_i = &vector_sums_i_value;
  GMVectorSums_create(vector_sums_i, vectors, env);

  GMVectors* vectors_i = vectors;

  GMMirroredBuf vectors_i_buf_value = GMMirroredBuf_null();
  GMMirroredBuf* const vectors_i_buf = &vectors_i_buf_value;
  GMMirroredBuf_create(vectors_i_buf, npvfl, nvl, env);

  /*------------------------*/
  /*---Allocations: Part 2---*/
  /*------------------------*/

  GMVectorSums vector_sums_j_value = GMVectorSums_null();
  GMVectorSums* const vector_sums_j = &vector_sums_j_value;
  GMVectorSums_create(vector_sums_j, vectors, env);

  GMVectors vectors_j_value_0 = GMVectors_null();
  GMVectors vectors_j_value_1 = GMVectors_null();
  GMVectors* vectors_j[2] = {&vectors_j_value_0, &vectors_j_value_1};
  GMVectors_create(vectors_j[0], data_type, vectors->num_field,
                   vectors->num_field_active, nvl, env);
  GMVectors_create(vectors_j[1], data_type, vectors->num_field,
                   vectors->num_field_active, nvl, env);

  GMMirroredBuf vectors_j_buf_value = GMMirroredBuf_null();
  GMMirroredBuf* const vectors_j_buf = &vectors_j_buf_value;
  GMMirroredBuf_create(vectors_j_buf, npvfl, nvl, env);

  /*------------------------*/
  /*---Allocations: Part 3---*/
  /*------------------------*/

  GMVectorSums vector_sums_k_value = GMVectorSums_null();
  GMVectorSums* const vector_sums_k = &vector_sums_k_value;
  GMVectorSums_create(vector_sums_k, vectors, env);

  GMVectors vectors_k_value_0 = GMVectors_null();
  GMVectors vectors_k_value_1 = GMVectors_null();
  GMVectors* vectors_k[2] = {&vectors_k_value_0, &vectors_k_value_1};
  GMVectors_create(vectors_k[0], data_type, vectors->num_field,
                   vectors->num_field_active, nvl, env);
  GMVectors_create(vectors_k[1], data_type, vectors->num_field,
                   vectors->num_field_active, nvl, env);

  GMMirroredBuf vectors_k_buf_value = GMMirroredBuf_null();
  GMMirroredBuf* const vectors_k_buf = &vectors_k_buf_value;
  GMMirroredBuf_create(vectors_k_buf, npvfl, nvl, env);

  /*------------------------*/
  /*---Prepare to compute---*/
  /*------------------------*/

  GMVectors* vectors_j_prev = NULL;
  GMVectors* vectors_k_prev = NULL;

  GMMirroredBuf* vectors_j_buf_prev = NULL;
  GMMirroredBuf* vectors_k_buf_prev = NULL;

  GMVectorSums* vector_sums_j_prev = NULL;
  GMVectorSums* vector_sums_k_prev = NULL;

  int j_block_prev = -1;
  int k_block_prev = -1;

  int section_step_prev = -1;

  bool have_unprocessed_section_block = false;

  GMComputeNumerators3Way gm_compute_numerators_3way = {0};
  GMComputeNumerators3Way_create(&gm_compute_numerators_3way, nvl, npvfl, env);

  // Counter for quantum of work:
  //   for part 1 or part 2: 1/6 section of work needed for block
  //   for part 3: all work needed for block
  int section_block_num = 0;

  int index_j_comm = 0;
  int index_k_comm = 0;

  // The following three sections process part 1, 2 and 3 blocks in
  // sequence. Conceptually they can be thought of as a single looping
  // process over section blocks.
  // To achieve overlap of communication with computation, there is a
  // "lag" in the computation: the communication is performed on the relevant
  // active step, but computation is lagged one loop cycle so as to
  // compute on data that is already communicated.

  /*------------------------*/
  /*---Part 1 Computation: tetrahedron---*/
  /*------------------------*/

  /*---Denominator---*/
  GMVectorSums_compute(vector_sums_i, vectors_i, env);

  /*---Copy in vectors---*/
  gm_vectors_to_buf(vectors_i_buf, vectors_i, env);

  /*---Send vectors to GPU---*/
  gm_set_vectors_start(vectors_i, vectors_i_buf, env);
  gm_set_vectors_wait(env);

  for (int section_step=0; section_step<gm_num_section_steps(env, 1);
       ++section_step) {
    if (gm_proc_r_active(section_block_num, env)) {

      if (have_unprocessed_section_block) {
        /*---Compute numerators---*/
        GMComputeNumerators3Way_start(
          &gm_compute_numerators_3way,
          vectors_i, vectors_j_prev, vectors_k_prev, metrics,
          vectors_i_buf, vectors_j_buf_prev, vectors_k_buf_prev,
          j_block_prev, k_block_prev,
          vector_sums_i, vector_sums_j_prev, vector_sums_k_prev,
          section_step_prev, env);
        gm_compute_wait(env);
        have_unprocessed_section_block = false;
      }

      /*---Remember processing to do next time---*/
      vectors_j_prev = vectors_i;
      vectors_k_prev = vectors_i;
      vectors_j_buf_prev = vectors_i_buf;
      vectors_k_buf_prev = vectors_i_buf;
      vector_sums_j_prev = vector_sums_i;
      vector_sums_k_prev = vector_sums_i;
      j_block_prev = i_block;
      k_block_prev = i_block;
      section_step_prev = section_step;
      have_unprocessed_section_block = true;

    } /*---if (section_block_num ...)---*/
    ++section_block_num;
  } /*---section_step---*/

  /*------------------------*/
  /*---Part 2 Computation: triangular prisms---*/
  /*------------------------*/

  for (int section_step=0; section_step<gm_num_section_steps(env, 2);
       ++section_step) {
    for (int j_i_offset = 1; j_i_offset < num_block; ++j_i_offset) {

      const int j_block = gm_mod_i(i_block + j_i_offset, num_block);

      //TODO: can possibly simplify this - mod by num_proc_i instead

      const int proc_send_j = gm_mod_i(proc_num_ir - j_i_offset*num_proc_r,
                                       num_proc_ir);
      const int proc_recv_j = gm_mod_i(proc_num_ir + j_i_offset*num_proc_r,
                                       num_proc_ir);

      if (gm_proc_r_active(section_block_num, env)) {

        /*---Communicate vectors start---*/
        GMVectors* const vectors_j_this = vectors_j[index_j_comm];
        index_j_comm = 1 - index_j_comm;
        MPI_Request req_recv_j = gm_recv_vectors_start(vectors_j_this,
                                       proc_recv_j, 0+3*section_block_num, env);
        MPI_Request req_send_j = gm_send_vectors_start(vectors_i,
                                       proc_send_j, 0+3*section_block_num, env);

        if (have_unprocessed_section_block) {
          /*---Compute numerators---*/
          GMComputeNumerators3Way_start(
            &gm_compute_numerators_3way,
            vectors_i, vectors_j_prev, vectors_k_prev, metrics,
            vectors_i_buf, vectors_j_buf_prev, vectors_k_buf_prev,
            j_block_prev, k_block_prev,
            vector_sums_i, vector_sums_j_prev, vector_sums_k_prev,
            section_step_prev, env);
          gm_compute_wait(env);
          have_unprocessed_section_block = false;
        }

        /*---Communicate vectors wait---*/
        gm_send_vectors_wait(&req_send_j, env);
        gm_recv_vectors_wait(&req_recv_j, env);

        /*---Copy in vectors---*/
        gm_vectors_to_buf(vectors_j_buf, vectors_j_this, env);

        /*---Send vectors to GPU start---*/
        gm_set_vectors_start(vectors_j_this, vectors_j_buf, env);

        /*---Denominator---*/
        GMVectorSums_compute(vector_sums_j, vectors_j_this, env);

        /*---Send vectors to GPU wait---*/
        gm_set_vectors_wait(env);

        /*---Remember processing to do next time---*/
        vectors_j_prev = vectors_j_this;
        vectors_k_prev = vectors_j_prev;
        vectors_j_buf_prev = vectors_j_buf;
        vectors_k_buf_prev = vectors_j_buf;
        vector_sums_j_prev = vector_sums_j;
        vector_sums_k_prev = vector_sums_j;
        j_block_prev = j_block;
        k_block_prev = j_block;
        section_step_prev = section_step;
        have_unprocessed_section_block = true;

      } /*---if (section_block_num ...)---*/
      ++section_block_num;
    } /*---j_i_offset---*/
  } /*---section_step---*/

  /*------------------------*/
  /*---Part 3 Computation: block sections---*/
  /*------------------------*/

  int k_block_currently_resident = -1;

  for (int section_step=0; section_step<gm_num_section_steps(env, 3);
       ++section_step) {
    for (int k_i_offset = 1; k_i_offset < num_block; ++k_i_offset) {
      const int k_block = gm_mod_i(i_block + k_i_offset, num_block);

      const int proc_send_k = gm_mod_i(proc_num_ir - k_i_offset*num_proc_r,
                                       num_proc_ir);
      const int proc_recv_k = gm_mod_i(proc_num_ir + k_i_offset*num_proc_r,
                                       num_proc_ir);

      for (int j_i_offset = 1; j_i_offset < num_block; ++j_i_offset){

        const int j_block = gm_mod_i(i_block + j_i_offset, num_block);

        const int proc_send_j = gm_mod_i(proc_num_ir-j_i_offset*num_proc_r,
                                         num_proc_ir);
        const int proc_recv_j = gm_mod_i(proc_num_ir+j_i_offset*num_proc_r,
                                         num_proc_ir);
        if (j_block == k_block) {
          /*---NOTE: this condition occurs on all procs at exactly the same
               j/k iteration in lockstep, so there is no chance the immediately
               following communication will deadlock/mispair---*/
          continue;
        }
        GMAssertAlways((j_block == k_block) ==
                       (j_i_offset == k_i_offset));
        if (gm_proc_r_active(section_block_num, env)) {

          const bool do_k_comm = k_block != k_block_currently_resident;

          GMVectors* vectors_k_this = 0;
          MPI_Request req_send_k;
          MPI_Request req_recv_k;
          if (do_k_comm) {
            /*---Communicate vectors start---*/
            vectors_k_this = vectors_k[index_k_comm];
            index_k_comm = 1 - index_k_comm;
            // NOTE: in some cases may not need double buffer, one may be enough
            req_recv_k = gm_recv_vectors_start(vectors_k_this,
                                       proc_recv_k, 1+3*section_block_num, env);
            req_send_k = gm_send_vectors_start(vectors_i,
                                       proc_send_k, 1+3*section_block_num, env);
          }

          /*---Communicate vectors start---*/
          GMVectors* const vectors_j_this = vectors_j[index_j_comm];
          index_j_comm = 1 - index_j_comm;
          MPI_Request req_recv_j = gm_recv_vectors_start(vectors_j_this,
                                      proc_recv_j, 2+3*section_block_num, env);
          MPI_Request req_send_j = gm_send_vectors_start(vectors_i,
                                      proc_send_j, 2+3*section_block_num, env);

          if (have_unprocessed_section_block) {
            /*---Compute numerators---*/
            GMComputeNumerators3Way_start(
              &gm_compute_numerators_3way,
              vectors_i, vectors_j_prev, vectors_k_prev, metrics,
              vectors_i_buf, vectors_j_buf_prev, vectors_k_buf_prev,
              j_block_prev, k_block_prev,
              vector_sums_i, vector_sums_j_prev, vector_sums_k_prev,
              section_step_prev, env);
            gm_compute_wait(env);
            have_unprocessed_section_block = false;
          }

          if (do_k_comm) {
            /*---Communicate vectors wait---*/
            gm_send_vectors_wait(&req_send_k, env);
            gm_recv_vectors_wait(&req_recv_k, env);
            k_block_currently_resident = k_block;

            /*---Copy in vectors---*/
            gm_vectors_to_buf(vectors_k_buf, vectors_k_this, env);

            /*---Send vectors to GPU start---*/
            gm_set_vectors_start(vectors_k_this, vectors_k_buf, env);

            /*---Denominator---*/
            GMVectorSums_compute(vector_sums_k, vectors_k_this, env);

            /*---Send vectors to GPU wait---*/
            gm_set_vectors_wait(env);

            /*---Remember processing to do next time---*/
            vectors_k_prev = vectors_k_this;
            vectors_k_buf_prev = vectors_k_buf;
            vector_sums_k_prev = vector_sums_k;
            k_block_prev = k_block;
          }

          /*---Communicate vectors wait---*/
          gm_send_vectors_wait(&req_send_j, env);
          gm_recv_vectors_wait(&req_recv_j, env);

          /*---Copy in vectors---*/
          gm_vectors_to_buf(vectors_j_buf, vectors_j_this, env);

          /*---Send vectors to GPU start---*/
          gm_set_vectors_start(vectors_j_this, vectors_j_buf, env);

          /*---Denominator---*/
          GMVectorSums_compute(vector_sums_j, vectors_j_this, env);

          /*---Send vectors to GPU wait---*/
          gm_set_vectors_wait(env);

          /*---Remember processing to do next time---*/
          vectors_j_prev = vectors_j_this;
          vectors_j_buf_prev = vectors_j_buf;
          vector_sums_j_prev = vector_sums_j;
          j_block_prev = j_block;
          section_step_prev = section_step;
          have_unprocessed_section_block = true;

        } /*---if (section_block_num ...)---*/
        ++section_block_num;
      } /*---k_i_offset---*/
    }   /*---j_i_offset---*/
  } /*---section_step---*/

  /*------------------------*/
  /*---Cleanup---*/
  /*------------------------*/

  if (have_unprocessed_section_block) {
    /*---Compute numerators---*/
    GMComputeNumerators3Way_start(
      &gm_compute_numerators_3way,
      vectors_i, vectors_j_prev, vectors_k_prev, metrics,
      vectors_i_buf, vectors_j_buf_prev, vectors_k_buf_prev,
      j_block_prev, k_block_prev,
      vector_sums_i, vector_sums_j_prev, vector_sums_k_prev,
      section_step_prev, env);
    gm_compute_wait(env);
    have_unprocessed_section_block = false;
  }

  /*------------------------*/
  /*---Free memory and finalize---*/
  /*------------------------*/

  GMComputeNumerators3Way_destroy(&gm_compute_numerators_3way, env);

  GMVectors_destroy(vectors_k[0], env);
  GMVectors_destroy(vectors_k[1], env);
  GMVectors_destroy(vectors_j[0], env);
  GMVectors_destroy(vectors_j[1], env);

  GMVectorSums_destroy(vector_sums_k, env);
  GMVectorSums_destroy(vector_sums_j, env);
  GMVectorSums_destroy(vector_sums_i, env);

  GMMirroredBuf_destroy(vectors_k_buf, env);
  GMMirroredBuf_destroy(vectors_j_buf, env);
  GMMirroredBuf_destroy(vectors_i_buf, env);

  gm_linalg_finalize(env);
}

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
