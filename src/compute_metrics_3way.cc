//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way.cc
 * \author Wayne Joubert
 * \date   Thu Jan 21 19:07:47 EST 2016
 * \brief  Calculate metrics, 3-way.
 */
//-----------------------------------------------------------------------------
/*-----------------------------------------------------------------------------

Copyright 2020, UT-Battelle, LLC

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-----------------------------------------------------------------------------*/

#include "env.hh"
#include "linalg.hh"
#include "mirrored_buf.hh"
#include "magma_wrapper.hh"
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

ComputeMetrics3Way::ComputeMetrics3Way(GMDecompMgr& dm, CEnv& env)
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
  CEnv* const env = &env_;

  //---------------
  // Denominator
  //---------------

  VectorSums vector_sums(vectors, env_);

  //---------------
  // Numerator
  //---------------

  //MagmaWrapper::initialize(env_);
  MagmaWrapper magma_wrapper(env_);

  {

  const int nvl = vectors.num_vector_local;
  const int npfl = vectors.num_packedfield_local;

  // Allocate magma CPU memory for vectors and for result

  MirroredBuf vectors_buf(npfl, nvl, env_);

  // Copy in vectors

  gm_vectors_to_buf(&vectors_buf, &vectors, env);

  // Send vectors to GPU

  vectors_buf.to_accel();

  // Compute numerators

  const int section_step = 0;
  COMET_INSIST(gm_num_section_steps(env, 1) == 1 &&
           "not all2all case always has 1 section step.");

  {

  ComputeMetrics3WayBlock compute_metrics_3way_block(nvl, npfl, env_);

  compute_metrics_3way_block.compute(
    VData(&vectors, &vectors_buf, &vector_sums),
    VData(&vectors, &vectors_buf, &vector_sums),
    VData(&vectors, &vectors_buf, &vector_sums),
    metrics, env_.proc_num_vector(), env_.proc_num_vector(), section_step,
    magma_wrapper);

  //---------------
  // Terminations
  //---------------

  }

  }

  //MagmaWrapper::finalize(env_);
}

//-----------------------------------------------------------------------------
/// \brief Perform the 3-way metrics computation,  all2all case.

void ComputeMetrics3Way::compute_all2all_(GMMetrics& metrics,
                                          GMVectors& vectors) {
  CEnv* const env = &env_;

  // Initializations.

  MagmaWrapper magma_wrapper(env_);

  {

  const int nvl = metrics.num_vector_local;
  const int npfl = vectors.num_packedfield_local;

  const int data_type = env->data_type_vectors();

  const int num_block = env->num_block_vector();

  const int i_block = env->proc_num_vector();

  const int proc_num_repl = env->proc_num_repl();

  // ------------------
  // Allocations: Part 1.
  // ------------------

  VectorSums vector_sums_i_value(vectors.num_vector_local, env_);
  VectorSums* const vector_sums_i = &vector_sums_i_value;

  GMVectors* const vectors_i = &vectors;

  MirroredBuf vectors_i_buf_value(npfl, nvl, env_);
  MirroredBuf* const vectors_i_buf = &vectors_i_buf_value;

  // ------------------
  // Allocations: Part 2.
  // ------------------

  VectorSums vector_sums_j_value(vectors.num_vector_local, env_);
  VectorSums* const vector_sums_j = &vector_sums_j_value;

  GMVectors vectors_j_value_0;
  GMVectors vectors_j_value_1;
  GMVectors* const vectors_j[2] = {&vectors_j_value_0, &vectors_j_value_1};
  vectors_j[0]->create(data_type, *vectors.dm, env_);
  vectors_j[1]->create(data_type, *vectors.dm, env_);

  MirroredBuf vectors_j_buf_value(npfl, nvl,env_);
  MirroredBuf* const vectors_j_buf = &vectors_j_buf_value;

  // ------------------
  // Allocations: Part 3.
  // ------------------

  VectorSums vector_sums_k_value(vectors.num_vector_local, env_);
  VectorSums* const vector_sums_k = &vector_sums_k_value;

  GMVectors vectors_k_value_0;
  GMVectors vectors_k_value_1;
  GMVectors* const vectors_k[2] = {&vectors_k_value_0, &vectors_k_value_1};
  vectors_k[0]->create(data_type, *vectors.dm, env_);
  vectors_k[1]->create(data_type, *vectors.dm, env_);

  MirroredBuf vectors_k_buf_value(npfl, nvl,env_);
  MirroredBuf* const vectors_k_buf = &vectors_k_buf_value;

  // ------------------
  // Prepare to compute.
  // ------------------

  GMVectors* vectors_j_prev = NULL;
  GMVectors* vectors_k_prev = NULL;

  MirroredBuf* vectors_j_buf_prev = NULL;
  MirroredBuf* vectors_k_buf_prev = NULL;

  VectorSums* vector_sums_j_prev = NULL;
  VectorSums* vector_sums_k_prev = NULL;

  int j_block_prev = -1;
  int k_block_prev = -1;

  int section_step_prev = -1;

  bool have_unprocessed_section_block = false;

  ComputeMetrics3WayBlock compute_metrics_3way_block(nvl, npfl, env_);

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

  // ------------------
  // Part 1 Computation: tetrahedron.
  // ------------------

  // Denominator.
  vector_sums_i->compute(*vectors_i);
  if (env_.is_threshold_tc())
    vector_sums_i->to_accel();

  // Copy in vectors.
  gm_vectors_to_buf(vectors_i_buf, vectors_i, env);

  // Send vectors to GPU.
  vectors_i_buf->to_accel();

  const int num_section_steps_1 = gm_num_section_steps(env, 1);
  for (int section_step=0; section_step<num_section_steps_1; ++section_step) {
    if (metrics_is_proc_repl_active(metrics, section_block_num, env_)) {

      if (have_unprocessed_section_block) {
        // Compute numerators.
        compute_metrics_3way_block.compute(
          VData(vectors_i, vectors_i_buf, vector_sums_i),
          VData(vectors_j_prev, vectors_j_buf_prev, vector_sums_j_prev),
          VData(vectors_k_prev, vectors_k_buf_prev, vector_sums_k_prev),
          metrics, j_block_prev, k_block_prev, section_step_prev,
          magma_wrapper);
        have_unprocessed_section_block = false;
      }

      if (gm_is_section_block_in_phase(env, section_block_num)) {
        // Remember processing to do next time.
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

    } // if (section_block_num ...)
    ++section_block_num;
  } // section_step

  // ------------------
  // Part 2 Computation: triangular prisms.
  // ------------------

  GMVectors* vectors_j_recv = NULL;
  GMVectors* vectors_k_recv = NULL;
  if(0) printf("%zu\n",(size_t)vectors_j_recv);
  if(0) printf("%zu\n",(size_t)vectors_k_recv);

  GMVectors* vectors_j_this = 0;
  CommRequest comm_request_send_j;
  CommRequest comm_request_recv_j;
  //MPI_Request req_send_j;
  //MPI_Request req_recv_j;

  const int num_section_steps_2 = gm_num_section_steps(env, 2);
  for (int section_step=0; section_step<num_section_steps_2; ++section_step) {
    for (int j_i_offset = 1; j_i_offset < num_block; ++j_i_offset) {

      const int j_block = utils::mod_i(i_block + j_i_offset, num_block);

      const int proc_send_j = env_.proc_num_repl_vector(proc_num_repl,
        utils::mod_i(i_block - j_i_offset, num_block));

      const int proc_recv_j = env_.proc_num_repl_vector(proc_num_repl,
        utils::mod_i(i_block + j_i_offset, num_block));

      if (metrics_is_proc_repl_active(metrics, section_block_num, env_)) {

        if (gm_is_section_block_in_phase(env, section_block_num)) {
          // Communicate vectors start.
          vectors_j_this = vectors_j[index_j_comm];
          index_j_comm = 1 - index_j_comm; // toggle buffer num

          comm_recv_vectors_start(*vectors_j_this, proc_recv_j,
              0+3*section_block_num, comm_request_recv_j, *env);
          comm_send_vectors_start(*vectors_i, proc_send_j,
              0+3*section_block_num, comm_request_send_j, *env);
          vectors_j_recv = vectors_j_this;
          //req_recv_j = gm_recv_vectors_start(vectors_j_this, proc_recv_j,
          //                                   0+3*section_block_num, env);
          //req_send_j = gm_send_vectors_start(vectors_i, proc_send_j,
          //                                   0+3*section_block_num, env);
        } // if (gm_is_section_block_in_phase ...)
  
        if (have_unprocessed_section_block) {
          // Compute numerators.
          compute_metrics_3way_block.compute(
            VData(vectors_i, vectors_i_buf, vector_sums_i),
            VData(vectors_j_prev, vectors_j_buf_prev, vector_sums_j_prev),
            VData(vectors_k_prev, vectors_k_buf_prev, vector_sums_k_prev),
            metrics, j_block_prev, k_block_prev, section_step_prev,
            magma_wrapper);
          have_unprocessed_section_block = false;
        }

        if (gm_is_section_block_in_phase(env, section_block_num)) {
          // Communicate vectors wait.
          comm_request_send_j.wait();
          comm_request_recv_j.wait();
          //fprintf(stderr, "%zu %zu\n", comm_request_recv_j.cksum_, GMVectors_cksum(vectors_j_recv, env));
          COMET_ASSERT(comm_request_recv_j.cksum_ == GMVectors_cksum(vectors_j_recv, env));
          //gm_send_vectors_wait(&req_send_j, env);
          //gm_recv_vectors_wait(&req_recv_j, env);

          // Copy in vectors.
          gm_vectors_to_buf(vectors_j_buf, vectors_j_this, env);

          // Send vectors to GPU start.
          vectors_j_buf->to_accel_start();

          // Denominator.
          vector_sums_j->compute(*vectors_j_this);
          if (env_.is_threshold_tc())
            vector_sums_j->to_accel();

          // Send vectors to GPU wait.
          vectors_j_buf->to_accel_wait();

          // Remember processing to do next time.
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
        } // if (gm_is_section_block_in_phase ...)

      } // if (metrics_is_proc_repl_active ...)
      ++section_block_num;
    } // j_i_offset
  } // section_step

  // ------------------
  // Part 3 Computation: block sections.
  // ------------------

  GMVectors* vectors_k_this = 0;
  CommRequest comm_request_send_k;
  CommRequest comm_request_recv_k;
  //MPI_Request req_send_k;
  //MPI_Request req_recv_k;

  int k_block_currently_resident = -1;

  const int num_section_steps_3 = gm_num_section_steps(env, 3); // = 1
  for (int section_step=0; section_step<num_section_steps_3; ++section_step) {
    for (int k_i_offset = 1; k_i_offset < num_block; ++k_i_offset) {

      const int k_block = utils::mod_i(i_block + k_i_offset, num_block);

      const int proc_send_k = env_.proc_num_repl_vector(proc_num_repl,
        utils::mod_i(i_block - k_i_offset, num_block));

      const int proc_recv_k = env_.proc_num_repl_vector(proc_num_repl,
        utils::mod_i(i_block + k_i_offset, num_block));

      for (int j_i_offset = 1; j_i_offset < num_block; ++j_i_offset) {

        const int j_block = utils::mod_i(i_block + j_i_offset, num_block);

        const int proc_send_j = env_.proc_num_repl_vector(proc_num_repl,
          utils::mod_i(i_block - j_i_offset, num_block));

        const int proc_recv_j = env_.proc_num_repl_vector(proc_num_repl,
          utils::mod_i(i_block + j_i_offset, num_block));

        if (j_block == k_block) {
          /*---NOTE: this condition occurs on all procs at exactly the same
               j/k iteration in lockstep, so there is no chance the immediately
               following communication will deadlock/mispair---*/
          continue;
        }
        COMET_INSIST((j_block == k_block) == (j_i_offset == k_i_offset) &&
                  "Error in block indexing for communication.");
        if (metrics_is_proc_repl_active(metrics, section_block_num, env_)) {

          const bool do_k_comm = k_block != k_block_currently_resident;

          if (gm_is_section_block_in_phase(env, section_block_num)) {
            if (do_k_comm) {
              // Communicate vectors start.
              vectors_k_this = vectors_k[index_k_comm];
              index_k_comm = 1 - index_k_comm; // toggle buffer num
              // NOTE: in some cases may not need double buf, one may be enough
              comm_recv_vectors_start(*vectors_k_this, proc_recv_k,
                  1+3*section_block_num, comm_request_recv_k, *env);
              comm_send_vectors_start(*vectors_i, proc_send_k,
                  1+3*section_block_num, comm_request_send_k, *env);
              vectors_k_recv = vectors_k_this;
              //req_recv_k = gm_recv_vectors_start(vectors_k_this,
              //                         proc_recv_k, 1+3*section_block_num, env);
              //req_send_k = gm_send_vectors_start(vectors_i,
              //                         proc_send_k, 1+3*section_block_num, env);
            } // if do_k_comm

            // Communicate vectors start.
            vectors_j_this = vectors_j[index_j_comm];
            index_j_comm = 1 - index_j_comm; // toggle buffer num
            comm_recv_vectors_start(*vectors_j_this, proc_recv_j,
                2+3*section_block_num, comm_request_recv_j, *env);
            comm_send_vectors_start(*vectors_i, proc_send_j,
                2+3*section_block_num, comm_request_send_j, *env);
            vectors_j_recv = vectors_j_this;
            //req_recv_j = gm_recv_vectors_start(vectors_j_this, proc_recv_j,
            //                                   2+3*section_block_num, env);
            //req_send_j = gm_send_vectors_start(vectors_i, proc_send_j,
            //                                    2+3*section_block_num, env);
          } // if (gm_is_section_block_in_phase ...)

          if (have_unprocessed_section_block) {
            // Compute numerators.
            compute_metrics_3way_block.compute(
              VData(vectors_i, vectors_i_buf, vector_sums_i),
              VData(vectors_j_prev, vectors_j_buf_prev, vector_sums_j_prev),
              VData(vectors_k_prev, vectors_k_buf_prev, vector_sums_k_prev),
              metrics, j_block_prev, k_block_prev, section_step_prev,
              magma_wrapper);
            have_unprocessed_section_block = false;
          }

          if (gm_is_section_block_in_phase(env, section_block_num)) {
            if (do_k_comm) {
              // Communicate vectors wait.
              comm_request_send_k.wait();
              comm_request_recv_k.wait();
              
              COMET_ASSERT(comm_request_recv_k.cksum_ == GMVectors_cksum(vectors_k_recv, env));
              //gm_send_vectors_wait(&req_send_k, env);
              //gm_recv_vectors_wait(&req_recv_k, env);
              k_block_currently_resident = k_block;

              // Copy in vectors.
              gm_vectors_to_buf(vectors_k_buf, vectors_k_this, env);

              // Send vectors to GPU start.
              vectors_k_buf->to_accel_start();

              // Denominator.
              vector_sums_k->compute(*vectors_k_this);
              if (env_.is_threshold_tc())
                vector_sums_k->to_accel();

              // Send vectors to GPU wait.
              vectors_k_buf->to_accel_wait();

              // Remember processing to do next time.
              vectors_k_prev = vectors_k_this;
              vectors_k_buf_prev = vectors_k_buf;
              vector_sums_k_prev = vector_sums_k;
              k_block_prev = k_block;
            } // if do_k_comm

            // Communicate vectors wait.
            comm_request_send_j.wait();
            comm_request_recv_j.wait();
            COMET_ASSERT(comm_request_recv_j.cksum_ == GMVectors_cksum(vectors_j_recv, env));
            //gm_send_vectors_wait(&req_send_j, env);
            //gm_recv_vectors_wait(&req_recv_j, env);

            // Copy in vectors.
            gm_vectors_to_buf(vectors_j_buf, vectors_j_this, env);

            // Send vectors to GPU start.
            vectors_j_buf->to_accel_start();

            // Denominator.
            vector_sums_j->compute(*vectors_j_this);
            if (env_.is_threshold_tc())
              vector_sums_j->to_accel();

            // Send vectors to GPU wait.
            vectors_j_buf->to_accel_wait();

            // Remember processing to do next time.
            vectors_j_prev = vectors_j_this;
            vectors_j_buf_prev = vectors_j_buf;
            vector_sums_j_prev = vector_sums_j;
            j_block_prev = j_block;
            section_step_prev = section_step;
            have_unprocessed_section_block = true;
          } // if (gm_is_section_block_in_phase ...)

        } // if (metrics_is_proc_repl_active ...)
        ++section_block_num;
      } // k_i_offset
    }   // j_i_offset
  } // section_step

  // ------------------
  // Cleanup.
  // ------------------

  if (have_unprocessed_section_block) {
    // Compute numerators.
    compute_metrics_3way_block.compute(
      VData(vectors_i, vectors_i_buf, vector_sums_i),
      VData(vectors_j_prev, vectors_j_buf_prev, vector_sums_j_prev),
      VData(vectors_k_prev, vectors_k_buf_prev, vector_sums_k_prev),
      metrics, j_block_prev, k_block_prev, section_step_prev,
      magma_wrapper);
    have_unprocessed_section_block = false;
  }

  // ------------------
  // Free memory and finalize.
  // ------------------

  GMVectors_destroy(vectors_k[0], env);
  GMVectors_destroy(vectors_k[1], env);
  GMVectors_destroy(vectors_j[0], env);
  GMVectors_destroy(vectors_j[1], env);

  }
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
