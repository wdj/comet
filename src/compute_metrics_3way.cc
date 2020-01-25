//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way.cc
 * \author Wayne Joubert
 * \date   Thu Jan 21 19:07:47 EST 2016
 * \brief  Calculate metrics, 3-way.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "env.hh"
#include "linalg.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "comm_xfer_utils.hh"
#include "compute_metrics_3way_block.hh"
#include "compute_metrics_3way.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Constructor for ComputeMetrics3Way class.

ComputeMetrics3Way::ComputeMetrics3Way(GMDecompMgr& dm, GMEnv& env)
  : env_(env) {
  COMET_INSIST(env_.is_proc_active());

  COMET_INSIST_INTERFACE(&env_, env_.compute_method() != MetricType::DUO &&
                                "3-way DUO method not supported.");
}

//-----------------------------------------------------------------------------
/// \brief Destructor for ComputeMetrics3Way class.

ComputeMetrics3Way::~ComputeMetrics3Way() {
  COMET_INSIST(env_.is_proc_active());
}

//-----------------------------------------------------------------------------
/// \brief Perform the 3-way metrics computation.

void ComputeMetrics3Way::compute(GMMetrics& metrics, GMVectors& vectors) {
  COMET_INSIST(env_.is_proc_active());

  if (!env_.all2all()) {
    compute_notall2all_(metrics, vectors);
  } else {
    compute_all2all_(metrics, vectors);
  }
}

//-----------------------------------------------------------------------------
/// \brief Perform the 3-way metrics computation, non-all2all case.

void ComputeMetrics3Way::compute_notall2all_(GMMetrics& metrics,
                                             GMVectors& vectors) {
  Env* const env = &env_;

  //---------------
  // Denominator
  //---------------

  VectorSums vector_sums(vectors, env_);

  //GMVectorSums vector_sums = GMVectorSums_null();
  //GMVectorSums_create(&vector_sums, vectors.num_vector_local, env);
  //GMVectorSums_compute(&vector_sums, &vectors, env);

  //---------------
  // Numerator
  //---------------

  gm_linalg_initialize(env);

  {

  const int nvl = vectors.num_vector_local;
  const int npvfl = vectors.num_packedval_field_local;

  // Allocate magma CPU memory for vectors and for result

  GMMirroredBuf vectors_buf(npvfl, nvl, env_);

  // Copy in vectors

  gm_vectors_to_buf(&vectors_buf, &vectors, env);

  // Send vectors to GPU

  //gm_set_vectors_start(&vectors, &vectors_buf, env);
  //gm_set_vectors_wait(env);
  vectors_buf.to_accel();

  // Compute numerators

  const int section_step = 0;
  COMET_INSIST(gm_num_section_steps(env, 1) == 1 &&
           "not all2all case always has 1 section step.");

  {

  ComputeMetrics3WayBlock compute_metrics_3way_block(nvl, npvfl, env_);

  compute_metrics_3way_block.compute(
    VData(&vectors, &vectors_buf, &vector_sums),
    VData(&vectors, &vectors_buf, &vector_sums),
    VData(&vectors, &vectors_buf, &vector_sums),
    metrics, env_.proc_num_vector(), env_.proc_num_vector(), section_step);

  //---------------
  // Terminations
  //---------------

  }

  //GMVectorSums_destroy(&vector_sums, env);

  }

  gm_linalg_finalize(env);
}

//-----------------------------------------------------------------------------
/// \brief Perform the 3-way metrics computation,  all2all case.

void ComputeMetrics3Way::compute_all2all_(GMMetrics& metrics,
                                          GMVectors& vectors) {
  Env* const env = &env_;

  /*---Initializations---*/

  gm_linalg_initialize(env);

  {

  const int nvl = metrics.num_vector_local;
  const int npvfl = vectors.num_packedval_field_local;

  const int data_type = env->data_type_vectors();

  const int num_block = env->num_block_vector();

  const int i_block = env->proc_num_vector();

  const int proc_num_r = env->proc_num_repl();
  const int num_proc_r = env->num_proc_repl();

  /*---Create flattened index within space of procs assigned to
       vectors (non-field procs) - i.e., vector_i (=block) X repl ---*/

  const int proc_num_rv = proc_num_r + num_proc_r * i_block;
  const int num_proc_rv = num_block * num_proc_r;

  /*------------------------*/
  /*---Allocations: Part 1---*/
  /*------------------------*/

  VectorSums vector_sums_i_value(vectors.num_vector_local, env_);
  VectorSums* const vector_sums_i = &vector_sums_i_value;

  //GMVectorSums vector_sums_i_value = GMVectorSums_null();
  //GMVectorSums* const vector_sums_i = &vector_sums_i_value;
  //GMVectorSums_create(vector_sums_i, vectors.num_vector_local, env);

  GMVectors* const vectors_i = &vectors;

  GMMirroredBuf vectors_i_buf_value(npvfl, nvl, env_);
  GMMirroredBuf* const vectors_i_buf = &vectors_i_buf_value;

  /*------------------------*/
  /*---Allocations: Part 2---*/
  /*------------------------*/

  VectorSums vector_sums_j_value(vectors.num_vector_local, env_);
  VectorSums* const vector_sums_j = &vector_sums_j_value;

  //GMVectorSums vector_sums_j_value = GMVectorSums_null();
  //GMVectorSums* const vector_sums_j = &vector_sums_j_value;
  //GMVectorSums_create(vector_sums_j, vectors.num_vector_local, env);

  GMVectors vectors_j_value_0 = GMVectors_null();
  GMVectors vectors_j_value_1 = GMVectors_null();
  GMVectors* const vectors_j[2] = {&vectors_j_value_0, &vectors_j_value_1};
  GMVectors_create(vectors_j[0], data_type, vectors.dm, env);
  GMVectors_create(vectors_j[1], data_type, vectors.dm, env);

  GMMirroredBuf vectors_j_buf_value(npvfl, nvl,env_);
  GMMirroredBuf* const vectors_j_buf = &vectors_j_buf_value;

  /*------------------------*/
  /*---Allocations: Part 3---*/
  /*------------------------*/

  VectorSums vector_sums_k_value(vectors.num_vector_local, env_);
  VectorSums* const vector_sums_k = &vector_sums_k_value;

  //GMVectorSums vector_sums_k_value = GMVectorSums_null();
  //GMVectorSums* const vector_sums_k = &vector_sums_k_value;
  //GMVectorSums_create(vector_sums_k, vectors.num_vector_local, env);

  GMVectors vectors_k_value_0 = GMVectors_null();
  GMVectors vectors_k_value_1 = GMVectors_null();
  GMVectors* const vectors_k[2] = {&vectors_k_value_0, &vectors_k_value_1};
  GMVectors_create(vectors_k[0], data_type, vectors.dm, env);
  GMVectors_create(vectors_k[1], data_type, vectors.dm, env);

  GMMirroredBuf vectors_k_buf_value(npvfl, nvl,env_);
  GMMirroredBuf* const vectors_k_buf = &vectors_k_buf_value;

  /*------------------------*/
  /*---Prepare to compute---*/
  /*------------------------*/

  GMVectors* vectors_j_prev = NULL;
  GMVectors* vectors_k_prev = NULL;

  GMMirroredBuf* vectors_j_buf_prev = NULL;
  GMMirroredBuf* vectors_k_buf_prev = NULL;

  VectorSums* vector_sums_j_prev = NULL;
  VectorSums* vector_sums_k_prev = NULL;

  int j_block_prev = -1;
  int k_block_prev = -1;

  int section_step_prev = -1;

  bool have_unprocessed_section_block = false;

  ComputeMetrics3WayBlock compute_metrics_3way_block(nvl, npvfl, env_);

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
  vector_sums_i->compute(*vectors_i);
  //GMVectorSums_compute(vector_sums_i, vectors_i, env);

  /*---Copy in vectors---*/
  gm_vectors_to_buf(vectors_i_buf, vectors_i, env);

  /*---Send vectors to GPU---*/
  //gm_set_vectors_start(vectors_i, vectors_i_buf, env);
  //gm_set_vectors_wait(env);
  vectors_i_buf->to_accel();

  const int num_section_steps_1 = gm_num_section_steps(env, 1);
  for (int section_step=0; section_step<num_section_steps_1; ++section_step) {
    if (gm_proc_r_active(section_block_num, env)) {

      if (have_unprocessed_section_block) {
        /*---Compute numerators---*/
        compute_metrics_3way_block.compute(
          VData(vectors_i, vectors_i_buf, vector_sums_i),
          VData(vectors_j_prev, vectors_j_buf_prev, vector_sums_j_prev),
          VData(vectors_k_prev, vectors_k_buf_prev, vector_sums_k_prev),
          metrics, j_block_prev, k_block_prev, section_step_prev);
        have_unprocessed_section_block = false;
      }

      if (gm_is_section_block_in_phase(env, section_block_num)) {
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
      } // if

    } /*---if (section_block_num ...)---*/
    ++section_block_num;
  } /*---section_step---*/

  /*------------------------*/
  /*---Part 2 Computation: triangular prisms---*/
  /*------------------------*/

  GMVectors* vectors_j_this = 0;
  MPI_Request req_send_j;
  MPI_Request req_recv_j;

  const int num_section_steps_2 = gm_num_section_steps(env, 2);
  for (int section_step=0; section_step<num_section_steps_2; ++section_step) {
    for (int j_i_offset = 1; j_i_offset < num_block; ++j_i_offset) {

      const int j_block = utils::mod_i(i_block + j_i_offset, num_block);

      //TODO: can possibly simplify this - mod by num_proc_i instead

      const int proc_send_j = utils::mod_i(proc_num_rv - j_i_offset*num_proc_r,
                                       num_proc_rv);
      const int proc_recv_j = utils::mod_i(proc_num_rv + j_i_offset*num_proc_r,
                                       num_proc_rv);

      if (gm_proc_r_active(section_block_num, env)) {

        if (gm_is_section_block_in_phase(env, section_block_num)) {
          /*---Communicate vectors start---*/
          vectors_j_this = vectors_j[index_j_comm];
          index_j_comm = 1 - index_j_comm; // toggle buffer num
          req_recv_j = gm_recv_vectors_start(vectors_j_this, proc_recv_j,
                                             0+3*section_block_num, env);
          req_send_j = gm_send_vectors_start(vectors_i, proc_send_j,
                                             0+3*section_block_num, env);
        }
  
        if (have_unprocessed_section_block) {
          /*---Compute numerators---*/
          compute_metrics_3way_block.compute(
            VData(vectors_i, vectors_i_buf, vector_sums_i),
            VData(vectors_j_prev, vectors_j_buf_prev, vector_sums_j_prev),
            VData(vectors_k_prev, vectors_k_buf_prev, vector_sums_k_prev),
            metrics, j_block_prev, k_block_prev, section_step_prev);
          have_unprocessed_section_block = false;
        }

        if (gm_is_section_block_in_phase(env, section_block_num)) {
          /*---Communicate vectors wait---*/
          gm_send_vectors_wait(&req_send_j, env);
          gm_recv_vectors_wait(&req_recv_j, env);

          /*---Copy in vectors---*/
          gm_vectors_to_buf(vectors_j_buf, vectors_j_this, env);

          /*---Send vectors to GPU start---*/
          //gm_set_vectors_start(vectors_j_this, vectors_j_buf, env);
          vectors_j_buf->to_accel_start();

          /*---Denominator---*/
          vector_sums_j->compute(*vectors_j_this);
          //GMVectorSums_compute(vector_sums_j, vectors_j_this, env);

          /*---Send vectors to GPU wait---*/
          //gm_set_vectors_wait(env);
          vectors_j_buf->to_accel_wait();

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
        } // if

      } /*---if (section_block_num ...)---*/
      ++section_block_num;
    } /*---j_i_offset---*/
  } /*---section_step---*/

  /*------------------------*/
  /*---Part 3 Computation: block sections---*/
  /*------------------------*/

  GMVectors* vectors_k_this = 0;
  MPI_Request req_send_k;
  MPI_Request req_recv_k;

  int k_block_currently_resident = -1;

  const int num_section_steps_3 = gm_num_section_steps(env, 3); // = 1
  for (int section_step=0; section_step<num_section_steps_3; ++section_step) {
    for (int k_i_offset = 1; k_i_offset < num_block; ++k_i_offset) {
      const int k_block = utils::mod_i(i_block + k_i_offset, num_block);

      const int proc_send_k = utils::mod_i(proc_num_rv - k_i_offset*num_proc_r,
                                       num_proc_rv);
      const int proc_recv_k = utils::mod_i(proc_num_rv + k_i_offset*num_proc_r,
                                       num_proc_rv);

      for (int j_i_offset = 1; j_i_offset < num_block; ++j_i_offset){

        const int j_block = utils::mod_i(i_block + j_i_offset, num_block);

        const int proc_send_j = utils::mod_i(proc_num_rv-j_i_offset*num_proc_r,
                                         num_proc_rv);
        const int proc_recv_j = utils::mod_i(proc_num_rv+j_i_offset*num_proc_r,
                                         num_proc_rv);
        if (j_block == k_block) {
          /*---NOTE: this condition occurs on all procs at exactly the same
               j/k iteration in lockstep, so there is no chance the immediately
               following communication will deadlock/mispair---*/
          continue;
        }
        COMET_INSIST((j_block == k_block) == (j_i_offset == k_i_offset) &&
                  "Error in block indexing for communication.");
        if (gm_proc_r_active(section_block_num, env)) {

          const bool do_k_comm = k_block != k_block_currently_resident;

          if (gm_is_section_block_in_phase(env, section_block_num)) {
            if (do_k_comm) {
              /*---Communicate vectors start---*/
              vectors_k_this = vectors_k[index_k_comm];
              index_k_comm = 1 - index_k_comm; // toggle buffer num
              // NOTE: in some cases may not need double buf, one may be enough
              req_recv_k = gm_recv_vectors_start(vectors_k_this,
                                       proc_recv_k, 1+3*section_block_num, env);
              req_send_k = gm_send_vectors_start(vectors_i,
                                       proc_send_k, 1+3*section_block_num, env);
            } // if do_k_comm

            /*---Communicate vectors start---*/
            vectors_j_this = vectors_j[index_j_comm];
            index_j_comm = 1 - index_j_comm; // toggle buffer num
            req_recv_j = gm_recv_vectors_start(vectors_j_this, proc_recv_j,
                                               2+3*section_block_num, env);
            req_send_j = gm_send_vectors_start(vectors_i, proc_send_j,
                                                2+3*section_block_num, env);
          } // if

          if (have_unprocessed_section_block) {
            /*---Compute numerators---*/
            compute_metrics_3way_block.compute(
              VData(vectors_i, vectors_i_buf, vector_sums_i),
              VData(vectors_j_prev, vectors_j_buf_prev, vector_sums_j_prev),
              VData(vectors_k_prev, vectors_k_buf_prev, vector_sums_k_prev),
              metrics, j_block_prev, k_block_prev, section_step_prev);
            have_unprocessed_section_block = false;
          }

          if (gm_is_section_block_in_phase(env, section_block_num)) {
            if (do_k_comm) {
              /*---Communicate vectors wait---*/
              gm_send_vectors_wait(&req_send_k, env);
              gm_recv_vectors_wait(&req_recv_k, env);
              k_block_currently_resident = k_block;

              /*---Copy in vectors---*/
              gm_vectors_to_buf(vectors_k_buf, vectors_k_this, env);

              /*---Send vectors to GPU start---*/
              //gm_set_vectors_start(vectors_k_this, vectors_k_buf, env);
              vectors_k_buf->to_accel_start();

              /*---Denominator---*/
              vector_sums_k->compute(*vectors_k_this);
              //GMVectorSums_compute(vector_sums_k, vectors_k_this, env);

              /*---Send vectors to GPU wait---*/
              //gm_set_vectors_wait(env);
              vectors_k_buf->to_accel_wait();

              /*---Remember processing to do next time---*/
              vectors_k_prev = vectors_k_this;
              vectors_k_buf_prev = vectors_k_buf;
              vector_sums_k_prev = vector_sums_k;
              k_block_prev = k_block;
            } // if do_k_comm

            /*---Communicate vectors wait---*/
            gm_send_vectors_wait(&req_send_j, env);
            gm_recv_vectors_wait(&req_recv_j, env);

            /*---Copy in vectors---*/
            gm_vectors_to_buf(vectors_j_buf, vectors_j_this, env);

            /*---Send vectors to GPU start---*/
            //gm_set_vectors_start(vectors_j_this, vectors_j_buf, env);
            vectors_j_buf->to_accel_start();

            /*---Denominator---*/
            vector_sums_j->compute(*vectors_j_this);
            //GMVectorSums_compute(vector_sums_j, vectors_j_this, env);

            /*---Send vectors to GPU wait---*/
            //gm_set_vectors_wait(env);
            vectors_j_buf->to_accel_wait();

            /*---Remember processing to do next time---*/
            vectors_j_prev = vectors_j_this;
            vectors_j_buf_prev = vectors_j_buf;
            vector_sums_j_prev = vector_sums_j;
            j_block_prev = j_block;
            section_step_prev = section_step;
            have_unprocessed_section_block = true;
          } // if

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
    compute_metrics_3way_block.compute(
      VData(vectors_i, vectors_i_buf, vector_sums_i),
      VData(vectors_j_prev, vectors_j_buf_prev, vector_sums_j_prev),
      VData(vectors_k_prev, vectors_k_buf_prev, vector_sums_k_prev),
      metrics, j_block_prev, k_block_prev, section_step_prev);
    have_unprocessed_section_block = false;
  }

  /*------------------------*/
  /*---Free memory and finalize---*/
  /*------------------------*/

  GMVectors_destroy(vectors_k[0], env);
  GMVectors_destroy(vectors_k[1], env);
  GMVectors_destroy(vectors_j[0], env);
  GMVectors_destroy(vectors_j[1], env);

  //GMVectorSums_destroy(vector_sums_k, env);
  //GMVectorSums_destroy(vector_sums_j, env);
  //GMVectorSums_destroy(vector_sums_i, env);

  }

  gm_linalg_finalize(env);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
