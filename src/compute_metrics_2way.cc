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

  //---------------
  // Denominator
  //---------------

  vector_sums_onproc_.compute(vectors);

  //---------------
  // Numerator
  //---------------

  //MagmaWrapper::initialize(env_);
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

  ComputeMetrics2WayBlock::compute_nums_start(
    &vectors, &vectors, &metrics, &vectors_buf,
     &vectors_buf, metrics_buf_ptr,
     &vector_sums_onproc_,
     &vector_sums_onproc_,
     env_.proc_num_vector(),
     true, magma_wrapper, &env_);

  ComputeMetrics2WayBlock::compute_nums_wait(
    &vectors, &vectors, &metrics, &vectors_buf,
    &vectors_buf, metrics_buf_ptr,
    &vector_sums_onproc_,
    &vector_sums_onproc_,
    env_.proc_num_vector(),
    true, &env_);

  // Copy result from GPU

  metrics_buf_ptr->from_accel();
  gm_metrics_pad_adjust(&metrics, metrics_buf_ptr, &env_);

  // Do reduction across field procs if needed

  gm_reduce_metrics(&metrics, &metrics_buf, metrics_buf_ptr, &env_);

  // Combine

//>>>
  CompressedBuf matB_buf_compressed(metrics_buf, env_);
  ComputeMetrics2WayBlock::finalize(&metrics, &matB_buf_compressed,
                                    &vector_sums_onproc_, &vector_sums_onproc_,
                                    env_.proc_num_vector(), true, &env_);

  //---------------
  // Terminations
  //---------------

  }

  //MagmaWrapper::finalize(env_);
}

//=============================================================================
/// \brief Perform the 2-way metrics computation, all2all case.

void ComputeMetrics2Way::compute_all2all_(GMMetrics& metrics,
                                          GMVectors& vectors) {
  COMET_INSIST(env_.all2all());

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

  //const int j_i_offset_min = gm_bdiag_computed_min(&env_);
  //const int j_i_offset_max = gm_bdiag_computed_max(&env_);
  const int j_i_offset_min = metrics_bdiag_thisphase_min(env_);
  //const int j_i_offset_max = metrics_bdiag_thisphase_max(env_);
  //const int j_i_offset_this_row_max = gm_block_computed_this_row_max(&env_);
  const int j_i_offset_this_row_max = metrics_bdiag_thisphase_thisbrow_max(env_);

  //const int num_bdiag_computed = j_i_offset_max - j_i_offset_min;

  // Convenience struct to remember loop state across cycles.

  struct LoopVars {
    GMVectors* vectors_right;
    MirroredBuf* vectors_right_buf;
    MirroredBuf* metrics_buf;
    VectorSums* vector_sums_right;
    bool is_compute_step;
    bool is_first_compute_step;
    bool do_compute_block;
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
  // (note: at each step, num_proc_repl processors each compute a block)
  // NOTE: num_step should be consistent within same proc_r.

  //const int num_step = utils::ceil(num_bdiag_computed, num_proc_repl);
  const int num_step = metrics_2way_num_steps(env_);

  // Add extra step/ at begin/end to fill/drain pipeline.

  const int extra_step = 1;
  const int first_step = 0 - extra_step;

  //========================================
  for (int step_num = first_step; step_num < num_step+extra_step; ++step_num) {
  //========================================

    // Set per-step variables

    vars_prev = vars;
    vars = vars_next;

    vars_next.step_num = step_num + 1;
    vars_next.is_compute_step = vars_next.step_num >= 0 &&
                                vars_next.step_num < num_step;
    vars_next.is_first_compute_step = vars_next.step_num == 0;
    vars_next.index_01 = utils::mod_i(vars_next.step_num, 2);

    // Find bdiag number being computed.

//    vars_next.j_i_offset = j_i_offset_min +
//       proc_num_repl + num_proc_repl * vars_next.step_num;
    vars_next.j_i_offset = j_i_offset_min + (
       env_.is_comm_ring() ?
       vars_next.step_num + num_step * proc_num_repl :
       proc_num_repl + num_proc_repl * vars_next.step_num);

    vars_next.j_block = utils::mod_i(i_block + vars_next.j_i_offset, num_block);
    vars_next.is_main_diag = vars_next.j_i_offset == 0;

    // Ensure that the respective block is in the range of what is
    // to be computed in this phase.
    // This may be an issue if num_step * num_proc_repl doesn't equal
    // num computed this phase - due to divisibility.

    vars_next.do_compute_block = vars_next.is_compute_step &&
                   vars_next.j_i_offset < j_i_offset_this_row_max;

    // Pointers to left/right-side vecs.
    // Here we are computing V^T W, for V, W containing column vectors.

    vars_next.is_right_aliased = vars_next.is_main_diag;

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

    //========== MPI sends, receives - START

//    const int proc_send = utils::mod_i(proc_num_rv
//      - vars_next.j_i_offset * num_proc_repl,
//       num_proc_rv);

//    const int proc_recv = utils::mod_i(proc_num_rv
//      + vars_next.j_i_offset * num_proc_repl,
//      num_proc_rv);

//FIXRING
//    const int proc_send = env_.proc_num_repl_vector(proc_num_repl,
//      utils::mod_i(i_block - env_.is_comm_ring() ? 1 : vars_next.j_i_offset, num_block));
//    const int proc_recv = env_.proc_num_repl_vector(proc_num_repl,
//      utils::mod_i(i_block + env_.is_comm_ring() ? 1 : vars_next.j_i_offset, num_block));
// ISSUE: does metrics storage depend on this.
// probably need to revise gm_proc_r_active
// also Metrics_index_2_part2 ...

    const int proc_send = env_.proc_num_repl_vector(proc_num_repl,
      utils::mod_i(i_block - vars_next.j_i_offset, num_block));

    const int proc_recv = env_.proc_num_repl_vector(proc_num_repl,
      utils::mod_i(i_block + vars_next.j_i_offset, num_block));

//if (env_.num_proc() > 1)
//if (env_.proc_num_vector() == 0)
//printf("pn %i pnr %i step %i step_next %i ji %i ps %i pr %i  \n", env_.proc_num(), proc_num_repl, step_num, vars_next.step_num, vars_next.j_i_offset, proc_send, proc_recv );

    const bool do_comm = vars_next.j_i_offset < num_block &&
      ! vars_next.is_main_diag;

    if (vars_next.is_compute_step && do_comm) {
      const int mpi_tag = step_num + 1;
      COMET_INSIST((!vars_next.is_right_aliased) &&
               "Next step should always compute off-diag block.");
      // NOTE: the following order seems to help performance.
      mpi_requests[1] = gm_recv_vectors_start(vars_next.vectors_right,
                                              proc_recv, mpi_tag, &env_);
      mpi_requests[0] = gm_send_vectors_start(vectors_left,
                                              proc_send, mpi_tag, &env_);
//FIXRING
//      mpi_requests[0] = gm_send_vectors_start(
//        vars_next.step_num == 0 || !env_.is_comm_ring() ?
//          vectors_left : vars.vectors_right,
//                                              proc_send, mpi_tag, &env_);

    }

    //========== Send right matrix to GPU - WAIT.

    if (vars.is_compute_step && vars.do_compute_block &&
        ! vars.is_right_aliased)
      vars.vectors_right_buf->to_accel_wait();

    //========== Send left matrix to GPU on first step.

    if (vars_next.is_first_compute_step) {
      gm_vectors_to_buf(vectors_left_buf, vectors_left, &env_);
      vectors_left_buf->to_accel_start();
      // TODO: examine whether overlap possible.
      // May not be possible for general repl and phase (??).
      vectors_left_buf->to_accel_wait();
    }

    //========== Compute sums for denominators

    //const int compute_sums_this = CEnv_is_ppc64() ? 1 : 2;
    const int compute_sums_this = 0; //FIX env_.is_threshold_tc() ? 0 : 1;

    if (0 == compute_sums_this) { // needed here for is_threshold_tc
      if (vars.is_compute_step && vars.do_compute_block) {
        //TODO: possibly move this
        if (vars.is_first_compute_step)
          vector_sums_left->compute(*vectors_left);
        if (! vars.is_main_diag)
          vars.vector_sums_right->compute(*vars.vectors_right);
      }
    }

    //========== Perform pseudo GEMM - START

    if (vars.is_compute_step && vars.do_compute_block) {
      ComputeMetrics2WayBlock::compute_nums_start(
        vectors_left, vars.vectors_right, &metrics,
        vectors_left_buf, vars.vectors_right_buf, vars.metrics_buf,
        vector_sums_left, vars.vector_sums_right,
        vars.j_block, vars.is_main_diag, magma_wrapper, &env_);
    }

    //========== Copy result matrix from GPU - WAIT

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

        ComputeMetrics2WayBlock::finalize(
          &metrics,
          &matB_buf_compressed, 
          vector_sums_left, vars_prev.vector_sums_right,
          vars_prev.j_block,
          vars_prev.is_main_diag, &env_);
      }
    }

    // ISSUE: it may be possible to increase performance by swapping the
    // some code blocks below and the one code block above.  It depends
    // on the relative speeds.  If these would be put in two different
    // CPU threads, then it wouldn't matter.

    //========== Compute sums for denominators: case 1

    if (1 == compute_sums_this) { // put it here for speed on this arch
      if (vars.is_compute_step && vars.do_compute_block) {
        //TODO: possibly move this
        if (vars.is_first_compute_step) {
          vector_sums_left->compute(*vectors_left);
        }
        if (! vars.is_main_diag) {
          vars.vector_sums_right->compute(*vars.vectors_right);
        }
      }
    }

    //========== MPI receives - WAIT

    if (vars_next.is_compute_step && do_comm) {
      gm_recv_vectors_wait(&(mpi_requests[1]), &env_);
      COMET_INSIST((!vars_next.is_right_aliased) &&
               "Next step should always compute off-diag block.");
    }

    //========== Send right matrix to GPU - START.

    if (vars_next.is_compute_step && vars_next.do_compute_block &&
        ! vars_next.is_right_aliased) {
      // ISSUE: make sure not necessary if vars_next.is_right_aliased
      vars_next.vectors_right_buf->to_accel_start();
    }

    //========== Perform pseudo GEMM - WAIT

    if (vars.is_compute_step && vars.do_compute_block) {
      ComputeMetrics2WayBlock::compute_nums_wait(
        vectors_left, vars.vectors_right, &metrics,
        vectors_left_buf, vars.vectors_right_buf, vars.metrics_buf,
        vector_sums_left, vars.vector_sums_right,
        vars.j_block, vars.is_main_diag, &env_);
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
        }
        if (! vars.is_main_diag) {
          vars.vector_sums_right->compute(*vars.vectors_right);
        }
      }
    }

    //========== Combine numerators, denominators: CPU case

    if (!env_.is_using_linalg()) {
      if (vars.is_compute_step && vars.do_compute_block) {
        matB_buf_compressed.from_accel_wait();
        ComputeMetrics2WayBlock::finalize(
          &metrics,
          &matB_buf_compressed,
          vector_sums_left,
          vars.vector_sums_right, vars.j_block,
          vars.is_main_diag, &env_);
      }
    }

    //========== MPI sends - WAIT

    if (vars_next.is_compute_step && do_comm)
      gm_send_vectors_wait(&(mpi_requests[0]), &env_);

  //========================================
  } // step_num
  //========================================

  //MagmaWrapper::finalize(env_);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

