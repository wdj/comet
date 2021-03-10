//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_2way.cc
 * \author Wayne Joubert
 * \date   Thu Jan  7 10:21:09 EST 2016
 * \brief  Calculate metrics, 2-way.
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

#include "string.h"

#include "env.hh"
#include "linalg.hh"
#include "magma_wrapper.hh"
#include "mirrored_buf.hh"
#include "compressed_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "comm_xfer_utils.hh"
#include "compute_metrics_2way_block.hh"
#include "compute_metrics_2way.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Constructor for ComputeMetrics2Way class.

ComputeMetrics2Way::ComputeMetrics2Way(GMDecompMgr& dm, CEnv& env)
  : env_(env) 
  , vectors_01_{}
  , metrics_buf_0_(env)
  , metrics_buf_1_(env)
  , metrics_buf_01_{&metrics_buf_0_, &metrics_buf_1_}
  , vectors_buf_(env)
  , metrics_tmp_buf_(env)
  , vector_sums_onproc_(dm.num_vector_local, env)
  , vector_sums_offproc_0_(env.all2all() ? dm.num_vector_local : 0, env)
  , vector_sums_offproc_1_(env.all2all() ? dm.num_vector_local : 0, env)
  , vector_sums_offproc_01_{&vector_sums_offproc_0_, &vector_sums_offproc_1_} {
  COMET_INSIST(env_.is_proc_active());

  if (!env_.all2all())
    return;

  if(env.print_details()) printf("Creating ComputeMetrics2Way\n");

  for (int i = 0; i < NUM_BUF; ++i) {
    GMVectors_create_with_buf(&vectors_01_[i], env_.data_type_vectors(),
      &dm, &env_);
    metrics_buf_01_[i]->allocate(dm.num_vector_local, dm.num_vector_local);
  }

  vectors_buf_.allocate(dm.num_packedfield_local, dm.num_vector_local);

  if (env_.do_reduce())
    metrics_tmp_buf_.allocate(dm.num_vector_local, dm.num_vector_local);
}

//-----------------------------------------------------------------------------
/// \brief Destructor for ComputeMetrics2Way class.

ComputeMetrics2Way::~ComputeMetrics2Way() {
  COMET_INSIST(env_.is_proc_active());

  if (!env_.all2all())
    return;

  for (int i = 0; i < NUM_BUF; ++i) {
    GMVectors_destroy(&vectors_01_[i], &env_);
  }
}

//-----------------------------------------------------------------------------
/// \brief Perform the 2-way metrics computation.

void ComputeMetrics2Way::compute(GMMetrics& metrics, GMVectors& vectors) {
  COMET_INSIST(env_.is_proc_active());

  if (!env_.all2all()) {
    compute_notall2all_(metrics, vectors);
  } else {
    compute_all2all_(metrics, vectors);
  }
}

//-----------------------------------------------------------------------------
/// \brief Perform the 2-way metrics computation, non-all2all case.

void ComputeMetrics2Way::compute_notall2all_(GMMetrics& metrics,
                                             GMVectors& vectors) {
  COMET_INSIST(!env_.all2all());

  if(env_.print_details()) printf("In compute_notall2all\n");

  //---------------
  // Denominator
  //---------------

  vector_sums_onproc_.compute(vectors);
  if (env_.is_threshold_tc())
    vector_sums_onproc_.to_accel();

  //---------------
  // Numerator
  //---------------

  MagmaWrapper magma_wrapper(env_);

  {

  const int nvl = vectors.num_vector_local;
  const int npfl = vectors.num_packedfield_local;

  // Allocate memory for vectors and for result 

  MirroredBuf vectors_buf(npfl, nvl, env_);

  MirroredBuf metrics_buf(nvl, nvl, env_);

  MirroredBuf metrics_tmp_buf(env_);
  if (env_.do_reduce())
    metrics_tmp_buf.allocate(nvl, nvl);

  MirroredBuf* metrics_buf_ptr =
      env_.do_reduce() ?  &metrics_tmp_buf : &metrics_buf;

  // Copy in vectors

  gm_vectors_to_buf(&vectors_buf, &vectors, &env_);

  // Send vectors to GPU

  vectors_buf.to_accel();

  env_.gemm_start_timer.start();
  ComputeMetrics2WayBlock::compute_nums_start(
    &vectors, &vectors, &metrics, &vectors_buf,
     &vectors_buf, metrics_buf_ptr,
     &vector_sums_onproc_,
     &vector_sums_onproc_,
     env_.proc_num_vector(),
     true, magma_wrapper, &env_);
  env_.gemm_start_timer.end();

  env_.gemm_wait_timer.start();
  ComputeMetrics2WayBlock::compute_nums_wait(
    &vectors, &vectors, &metrics, &vectors_buf,
    &vectors_buf, metrics_buf_ptr,
    &vector_sums_onproc_,
    &vector_sums_onproc_,
    env_.proc_num_vector(),
    true, &env_);
  env_.gemm_wait_timer.end();

  // Copy result from GPU

  metrics_buf_ptr->from_accel();
  gm_metrics_pad_adjust(&metrics, metrics_buf_ptr, &env_);

  // Do reduction across field procs if needed

  gm_reduce_metrics(&metrics, &metrics_buf, metrics_buf_ptr, &env_);

  // Combine

  CompressedBuf matB_buf_compressed(metrics_buf, env_);
  env_.finalize_timer.start();
  ComputeMetrics2WayBlock::finalize(&metrics, &matB_buf_compressed,
                                    &vector_sums_onproc_, &vector_sums_onproc_,
                                    env_.proc_num_vector(), true, &env_);
  env_.finalize_timer.end();

  //---------------
  // Terminations
  //---------------

  }

  if(env_.print_details()) printf("Done in compute_notall2all\n");
}

//=============================================================================
/// \brief Perform the 2-way metrics computation, all2all case.

void ComputeMetrics2Way::compute_all2all_(GMMetrics& metrics,
                                          GMVectors& vectors) {
  COMET_INSIST(env_.all2all());

  int rank = 0;
  COMET_MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  if(env_.print_details()) printf("rank=%d In compute_all2all\n",rank);

  // Initializations

  const int num_block = env_.num_block_vector();
  const int i_block = env_.proc_num_vector();

  MagmaWrapper magma_wrapper(env_);

  // Create double buffer of vectors objects for send/recv

  // Allocate GPU buffers
  // To overlap transfers with compute, set up double buffers for the
  // vectors sent to the GPU and the metrics received from the GPU.

  // Result matrix is diagonal block and half the blocks to the right
  // (including wraparound to left side of matrix when appropriate).
  // For even number of vector blocks, block rows of lower half of matrix
  //  have one less block to make correct count.

  const int num_proc_repl = env_.num_proc_repl();
  const int proc_num_repl = env_.proc_num_repl();

//  const int num_proc_vector = env_.num_proc_vector();
//  const int proc_num_vector = env_.proc_num_vector();

  // Flatten the proc_vector and proc_repl indices into a single index.

  //const int num_proc_rv = num_block * num_proc_repl;
  //const int proc_num_rv = proc_num_repl + num_proc_repl * i_block;

  MPI_Request mpi_requests[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

  // Prepare for loop over blocks of result.

  /* Summary of the opertions in this loop:

    send VECTORS next step start
    recv VECTORS next step start
    set VECTORS this step wait
    compute numerators start
    get METRICS prev step wait
    combine prev step wait (GPU case)
    recv VECTORS next step wait
    set VECTORS next step start
    compute numerators wait
    get METRICS this step start
    compute denominators
    combine this step (CPU case)
    send VECTORS next step wait

  */

  // Lowest/highest (block) diag to be computed for this phase,
  // measured from (block) main diag.
  // For all repl procs.

  const int j_i_offset_min = metrics_bdiag_thisphase_min(env_);
  const int j_i_offset_this_row_max = metrics_bdiag_thisphase_thisbrow_max(env_);

  // Convenience struct to remember loop state across cycles.

  struct LoopVars {
    GMVectors* vectors_right;
    MirroredBuf* vectors_right_buf;
    MirroredBuf* metrics_buf;
    VectorSums* vector_sums_right;
    bool is_compute_step;
    bool is_first_compute_step;
    bool do_compute_block;
    bool needs_comm;
    bool is_main_diag;
    bool is_right_aliased;
    int step_num;
    int index_01;
    int j_i_offset;
    int j_block;
  };

  LoopVars vars = {};
  LoopVars vars_prev = {};
  LoopVars vars_next = {};

  // Optionally compress the result data on the GPU.

  CompressedBuf matB_buf_compressed(*metrics_buf_01_[0], env_);

  // Num steps to take to compute the blocks
  // Note: at each step, proc_repl procs each compute a block
  // (except for a few cases when it could be less)
  // NOTE: num_step is consistent across procs.
  // This is to ensure all communication is the same, even if some
  // procs don't need to compute on a step.

  const int num_step = metrics_num_steps_2way(env_);

  // Add extra step at begin/end to fill/drain pipeline.

  const int extra_step = 1;
  const int first_step = 0 - extra_step;

  //if(env_.print_details()) printf("rank=%d Starting for loop %d-%d\n",rank,0-extra_step,num_step+extra_step);

  //========================================
  for (int step_num = first_step; step_num < num_step+extra_step; ++step_num) {
  //========================================

    //if(env_.print_details()) printf("rank=%d In loop step_num=%d\n",rank,step_num);

    // Set per-step variables

    vars_prev = vars;
    vars = vars_next;

    vars_next.step_num = step_num + 1;
    vars_next.is_compute_step = vars_next.step_num >= 0 &&
                                vars_next.step_num < num_step;
    vars_next.is_first_compute_step = vars_next.step_num == 0;
    vars_next.index_01 = utils::mod_i(vars_next.step_num, 2);

    // Find bdiag number being computed.

    vars_next.j_i_offset = j_i_offset_min + (
       env_.is_comm_ring() ?
       vars_next.step_num + num_step * proc_num_repl :
       proc_num_repl + num_proc_repl * vars_next.step_num);

    vars_next.j_block = utils::mod_i(i_block + vars_next.j_i_offset, num_block);
    vars_next.is_main_diag = 0 == vars_next.j_i_offset;

    // Ensure that the respective block is in the range of what is
    // to be computed in this phase, this block row.
    // This is to handle the issue when num_step * num_proc_repl doesn't equal
    // num computed this phase - due to divisibility.

    vars_next.do_compute_block = vars_next.is_compute_step &&
      vars_next.j_i_offset < j_i_offset_this_row_max;

    // Pointers to left/right-side vecs.
    // Here we are computing V^T W, for V, W containing column vectors.

    vars_next.is_right_aliased = vars_next.is_main_diag;

    vars_next.needs_comm = vars_next.j_i_offset < num_block &&
      ! vars_next.is_main_diag;

    GMVectors* vectors_left = &vectors;
    vars_next.vectors_right = vars_next.is_right_aliased ?
      vectors_left : &vectors_01_[vars_next.index_01];

    MirroredBuf* vectors_left_buf = &vectors_buf_;
    vars_next.vectors_right_buf = vars_next.is_right_aliased ?
      vectors_left_buf : vectors_01_[vars_next.index_01].buf;

    VectorSums* vector_sums_left = &vector_sums_onproc_;
    vars_next.vector_sums_right = vars_next.is_right_aliased ?
      vector_sums_left : vector_sums_offproc_01_[vars_next.index_01];

    // Pointer to metrics buffer

    vars_next.metrics_buf = metrics_buf_01_[vars_next.index_01];

// TODO: move "Send left matrix to GPU" here if gpu direct.

    //========== MPI sends, receives - START

    if (vars_next.is_compute_step && vars_next.needs_comm) {

      const bool is_first_comm_step = ! vars.needs_comm; // && vars_next.needs_comm

      // Select send/recv procs to use.
      // if ring comm step, shift by 1; else shift vectors_left by offset.

      const int proc_offset = env_.is_comm_ring() && ! is_first_comm_step ?
        1 : vars_next.j_i_offset;
    
      const int proc_send = env_.proc_num_repl_vector(proc_num_repl,
        utils::mod_i(i_block - proc_offset, num_block));

      const int proc_recv = env_.proc_num_repl_vector(proc_num_repl,
        utils::mod_i(i_block + proc_offset, num_block));

      const int mpi_tag = step_num + 1;

      COMET_INSIST((!vars_next.is_right_aliased) &&
               "Next step should always compute off-diag block.");

      GMVectors* vectors_send = env_.is_comm_ring() && ! is_first_comm_step ?
        &vectors_01_[1-vars_next.index_01] :
        vectors_left;

      // Initiate sends/recvs for vecs needed on next step
      //if(env_.print_details()) printf("rank=%d Initiating sends/recvs\n",rank);
      GMVectors* vectors_recv = vars_next.vectors_right;

      // NOTE: the following order seems to help performance.
      mpi_requests[1] = gm_recv_vectors_start(vectors_recv,
                                              proc_recv, mpi_tag, &env_);
      mpi_requests[0] = gm_send_vectors_start(vectors_send,
                                              proc_send, mpi_tag, &env_);
    }

    //========== Send right matrix to GPU - WAIT.

    //if(env_.print_details()) printf("rank=%d Sending right vectors to GPU end\n",rank);

    if (vars.is_compute_step && vars.do_compute_block &&
        ! vars.is_right_aliased) {
      env_.vec2_wait_timer.record();
      env_.vec2_wait_timer.start();
      vars.vectors_right_buf->to_accel_wait();
      env_.vec2_wait_timer.end();
    }

    //========== Send left matrix to GPU on first step.
    //if(env_.print_details()) printf("rank=%d Sending left vecs to GPU\n",rank);

    if (vars_next.is_first_compute_step) {
      gm_vectors_to_buf(vectors_left_buf, vectors_left, &env_);
      env_.vec1_to_gpu_timer.record();
      env_.vec1_to_gpu_timer.start();
      vectors_left_buf->to_accel_start();
      env_.vec1_to_gpu_timer.end();
      // TODO: examine whether overlap possible.
      // May not be possible for general repl and phase (??).
      env_.vec1_wait_timer.record();
      env_.vec1_wait_timer.start();
      vectors_left_buf->to_accel_wait();
      env_.vec1_wait_timer.end();
    }

    //========== Compute sums for denominators
    //if(env_.print_details()) printf("rank=%d Computing sums for denominators\n",rank);

    //const int compute_sums_this = CEnv_is_ppc64() ? 1 : 2;
    const int compute_sums_this = 0; //FIX env_.is_threshold_tc() ? 0 : 1;

    if (0 == compute_sums_this) { // needed here for is_threshold_tc
      if (vars.is_compute_step && vars.do_compute_block) {
        //TODO: possibly move this
        if (vars.is_first_compute_step) {
          vector_sums_left->compute(*vectors_left);
          if (env_.is_threshold_tc())
            vector_sums_left->to_accel();;
        }
        if (! vars.is_main_diag) {
          vars.vector_sums_right->compute(*vars.vectors_right);
          if (env_.is_threshold_tc())
            vars.vector_sums_right->to_accel();
        }
      }
    }

    //========== Perform pseudo GEMM - START

    if (vars.is_compute_step && vars.do_compute_block) {
      if(env_.print_details()) printf("rank=%d Calling ComputeMetrics2WayBlock::compute_nums_start\n",rank);
      env_.gemm_start_timer.start();
      ComputeMetrics2WayBlock::compute_nums_start(
        vectors_left, vars.vectors_right, &metrics,
        vectors_left_buf, vars.vectors_right_buf, vars.metrics_buf,
        vector_sums_left, vars.vector_sums_right,
        vars.j_block, vars.is_main_diag, magma_wrapper, &env_);
      env_.gemm_start_timer.end();
      if(env_.print_details()) printf("rank=%d Done calling ComputeMetrics2WayBlock::compute_nums_start\n",rank);
    }

    //========== Copy result matrix from GPU - WAIT
    if(env_.print_details()) printf("rank=%d Copying result matrix from GPU\n",rank);

    if (env_.is_using_linalg()) {
      if (vars_prev.is_compute_step && vars_prev.do_compute_block) {
        matB_buf_compressed.from_accel_wait();
        gm_metrics_pad_adjust(&metrics, vars_prev.metrics_buf, &env_);

        //TODO: remove need to allocate metrics_tmp_buf device array
        MirroredBuf* metrics_buf_prev_ptr =
            env_.do_reduce() ?  &metrics_tmp_buf_ : vars_prev.metrics_buf;

        //========== Reduce along field procs

        if (env_.do_reduce()) {
          gm_reduce_metrics(&metrics, metrics_buf_prev_ptr,
                            vars_prev.metrics_buf, &env_);
          matB_buf_compressed.attach(*metrics_buf_prev_ptr);
        }

        //========== Combine numerators, denominators: CPU case

	if(env_.print_details()) printf("rank=%d Calling ComputeMetrics2WayBlock::finalize\n",rank);
        env_.finalize_timer.start();
	ComputeMetrics2WayBlock::finalize(
          &metrics,
          &matB_buf_compressed, 
          vector_sums_left, vars_prev.vector_sums_right,
          vars_prev.j_block,
          vars_prev.is_main_diag, &env_);
	env_.finalize_timer.end();
        if(env_.print_details()) printf("Done calling ComputeMetrics2WayBlock::finalize\n");
      }
    }

    // ISSUE: it may be possible to increase performance by swapping the
    // some code blocks below and the one code block above.  It depends
    // on the relative speeds.  If these would be put in two different
    // CPU threads, then it wouldn't matter.

    //========== Compute sums for denominators: case 1
    //if(env_.print_details()) printf("rank=%d computing sums for denominators\n",rank);

    if (1 == compute_sums_this) { // put it here for speed on this arch
      if (vars.is_compute_step && vars.do_compute_block) {
        //TODO: possibly move this
        if (vars.is_first_compute_step) {
          vector_sums_left->compute(*vectors_left);
          if (env_.is_threshold_tc())
            vector_sums_left->to_accel();
        }
        if (! vars.is_main_diag) {
          vars.vector_sums_right->compute(*vars.vectors_right);
          if (env_.is_threshold_tc())
            vars.vector_sums_right->to_accel();
        }
      }
    }

    //========== MPI receives - WAIT
    //if(env_.print_details()) printf("rank=%d Waiting for recvs to complete\n",rank);

    if (vars_next.is_compute_step && vars_next.needs_comm) {
      gm_recv_vectors_wait(&(mpi_requests[1]), &env_);
      COMET_INSIST((!vars_next.is_right_aliased) &&
               "Next step should always compute off-diag block.");
    }

    //========== Send right matrix to GPU - START.

    if (vars_next.is_compute_step && vars_next.do_compute_block &&
        ! vars_next.is_right_aliased) {
      // ISSUE: make sure not necessary if vars_next.is_right_aliased
      env_.vec2_to_gpu_timer.record();
      env_.vec2_to_gpu_timer.start();
      vars_next.vectors_right_buf->to_accel_start();
      env_.vec2_to_gpu_timer.end();
    }

    //========== Perform pseudo GEMM - WAIT

    if (vars.is_compute_step && vars.do_compute_block) {
      if(env_.print_details()) printf("rank=%d Calling compute_nums_wait\n",rank);
      env_.gemm_wait_timer.start();
      ComputeMetrics2WayBlock::compute_nums_wait(
        vectors_left, vars.vectors_right, &metrics,
        vectors_left_buf, vars.vectors_right_buf, vars.metrics_buf,
        vector_sums_left, vars.vector_sums_right,
        vars.j_block, vars.is_main_diag, &env_);
      env_.gemm_wait_timer.end();
      if(env_.print_details()) printf("rank=%d Done calling compute_nums_wait\n",rank);
      matB_buf_compressed.attach(*vars.metrics_buf);
      matB_buf_compressed.compress();
    }

    //========== Copy result matrix from GPU - START

    if (vars.is_compute_step && vars.do_compute_block)
      matB_buf_compressed.from_accel_start();

    //========== Compute sums for denominators: case 2

    if (2 == compute_sums_this) { // put it here for speed on this arch
      if (vars.is_compute_step && vars.do_compute_block) {
        //TODO: possibly move this
        if (vars.is_first_compute_step) {
          vector_sums_left->compute(*vectors_left);
          if (env_.is_threshold_tc())
            vector_sums_left->to_accel();
        }
        if (! vars.is_main_diag) {
          vars.vector_sums_right->compute(*vars.vectors_right);
          if (env_.is_threshold_tc())
            vars.vector_sums_right->to_accel();
        }
      }
    }

    //========== Combine numerators, denominators: CPU case
    //if(env_.print_details()) printf("rank=%d Combining numerators and denominators\n",rank);

    if (!env_.is_using_linalg()) {
      if (vars.is_compute_step && vars.do_compute_block) {
        matB_buf_compressed.from_accel_wait();
        if(env_.print_details()) printf("rank=%d Calling 2nd Block::finalize\n",rank);
        env_.finalize_timer.start();
	ComputeMetrics2WayBlock::finalize(
          &metrics,
          &matB_buf_compressed,
          vector_sums_left,
          vars.vector_sums_right, vars.j_block,
          vars.is_main_diag, &env_);
        env_.finalize_timer.end();
	if(env_.print_details()) printf("rank=%d Done calling 2nd Block::finalize\n",rank);
      }
    }

    //========== MPI sends - WAIT

    if (vars_next.is_compute_step && vars_next.needs_comm)
      gm_send_vectors_wait(&(mpi_requests[0]), &env_);

    //if(env_.print_details()) printf("rank=%d Done with loop step_num=%d\n",rank,step_num);

  //========================================
  } // step_num
  //========================================

  if(env_.print_details()) printf("rank=%d Done in compute_all2all\n",rank);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

