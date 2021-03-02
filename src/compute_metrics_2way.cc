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

  //MagmaWrapper::initialize(env_);
  MagmaWrapper magma_wrapper(env_);

  // Create double buffer of vectors objects for send/recv

  // Allocate GPU buffers
  // To overlap transfers with compute, set up double buffers for the
  // vectors sent to the GPU and the metrics received from the GPU.

  // Result matrix is diagonal block and half the blocks to the right
  // (including wraparound to left side of matrix when appropriate).
  // For even number of vector blocks, block rows of lower half of matrix
  //  have one less block to make correct count.

  const int num_proc_r = env_.num_proc_repl();
  const int proc_num_r = env_.proc_num_repl();

  // Flatten the proc_vector and proc_repl indices into a single index.

  const int num_proc_rv = num_block * num_proc_r;
  const int proc_num_rv = proc_num_r + num_proc_r * i_block;

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

  // Add extra step at begin/end to fill/drain pipeline.

  const int extra_step = 1;

  // Lowest/highest (block) diag to be computed for this phase,
  // measured from (block) main diag.
  // For all repl procs.

  const int j_i_offset_min = gm_bdiag_computed_min(&env_);
  const int j_i_offset_max = gm_bdiag_computed_max(&env_);
  const int j_i_offset_this_row_max = gm_block_computed_this_row_max(&env_);

  const int num_bdiag_computed = j_i_offset_max - j_i_offset_min;

  // Num steps to take to compute blocks
  // (note: at each step, num_proc_r processors each compute a block)
  // NOTE: num_step should be consistent within same proc_r.

  const int num_step = utils::ceil(num_bdiag_computed, num_proc_r);

  typedef struct {
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
  } LoopVars;

  LoopVars vars = {};
  LoopVars vars_prev = {};
  LoopVars vars_next = {};

  // Use locks to verify no race condition on a buffer.
  // Lock buffer when in use for read or write, unlock when done.

  bool lock_vectors_01_buf_h[2] = {false, false};
  bool lock_vectors_01_buf_d[2] = {false, false};
  bool lock_metrics_buf_01_h[2] = {false, false};
  bool lock_metrics_buf_01_d[2] = {false, false};
  bool lock_vectors_buf_h = false;
  bool lock_vectors_buf_d = false;
  //bool lock_metrics_tmp_buf_d = false; // Not needed
  bool lock_metrics_tmp_buf_h = false;

//>>>
  CompressedBuf matB_buf_compressed(*metrics_buf_01_[0], env_);

  //if(env_.print_details()) printf("rank=%d Starting for loop %d-%d\n",rank,0-extra_step,num_step+extra_step);

  //========================================
  for (int step_num = 0-extra_step; step_num < num_step+extra_step; ++step_num){
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
    vars_next.j_i_offset = j_i_offset_min + vars_next.step_num * num_proc_r
                           + proc_num_r;
    vars_next.is_main_diag = vars_next.j_i_offset == 0;
    vars_next.j_block = utils::mod_i(i_block + vars_next.j_i_offset, num_block);
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

    // Set up lock aliases

    bool& lock_metrics_buf_ptr_h_prev
                                  = lock_metrics_buf_01_h[vars_prev.index_01];
    bool& lock_metrics_buf_ptr_d_prev
                                  = lock_metrics_buf_01_d[vars_prev.index_01];

    bool& lock_metrics_buf_ptr_h = lock_metrics_buf_01_h[vars.index_01];
    bool& lock_metrics_buf_ptr_d = lock_metrics_buf_01_d[vars.index_01];

    bool& lock_vectors_left_buf_h = lock_vectors_buf_h;
    bool& lock_vectors_left_buf_d = lock_vectors_buf_d;

    bool& lock_vectors_right_buf_h_next = vars_next.is_right_aliased ?
      lock_vectors_left_buf_h : lock_vectors_01_buf_h[vars_next.index_01];

    bool& lock_vectors_right_buf_d_next = vars_next.is_right_aliased ?
      lock_vectors_left_buf_d : lock_vectors_01_buf_d[vars_next.index_01];

    bool& lock_vectors_right_buf_h = vars.is_right_aliased ?
      lock_vectors_left_buf_h : lock_vectors_01_buf_h[vars.index_01];

    bool& lock_vectors_right_buf_d = vars.is_right_aliased ?
      lock_vectors_left_buf_d : lock_vectors_01_buf_d[vars.index_01];

    // Prepare for sends/recvs: procs for communication

    const int proc_send = utils::mod_i(proc_num_rv
        - vars_next.j_i_offset*num_proc_r, num_proc_rv);

    const int proc_recv = utils::mod_i(proc_num_rv
        + vars_next.j_i_offset*num_proc_r, num_proc_rv);

    const bool comm_with_self = vars_next.is_main_diag;

    // Initiate sends/recvs for vecs needed on next step
    //if(env_.print_details()) printf("rank=%d Initiating sends/recvs\n",rank);
    if (vars_next.is_compute_step && ! comm_with_self) {
      const int mpi_tag = step_num + 1;
      // NOTE: the following order helps performance
      COMET_INSIST((!vars_next.is_right_aliased) &&
               "Next step should always compute off-diag block.");
      lock(lock_vectors_right_buf_h_next);
      mpi_requests[1] = gm_recv_vectors_start(vars_next.vectors_right,
                                              proc_recv, mpi_tag, &env_);
      mpi_requests[0] = gm_send_vectors_start(vectors_left,
                                              proc_send, mpi_tag, &env_);
    }

    // Send right vectors to GPU end
    //if(env_.print_details()) printf("rank=%d Sending right vectors to GPU end\n",rank);
    if (vars.is_compute_step && vars.do_compute_block &&
        ! vars.is_right_aliased) {
      vars.vectors_right_buf->to_accel_wait();
      unlock(lock_vectors_right_buf_h);
      unlock(lock_vectors_right_buf_d);
    }

    // First step (for any repl or phase): send (left) vecs to GPU
    //if(env_.print_details()) printf("rank=%d Sending left vecs to GPU\n",rank);
    if (vars_next.is_first_compute_step) {
      lock(lock_vectors_left_buf_h);
      gm_vectors_to_buf(vectors_left_buf, vectors_left, &env_);
      lock(lock_vectors_left_buf_d);
      vectors_left_buf->to_accel_start();
      // TODO: examine whether overlap possible.
      // May not be possible for general repl and phase (??).
      vectors_left_buf->to_accel_wait();
      unlock(lock_vectors_left_buf_h);
      unlock(lock_vectors_left_buf_d);
    }

    // Compute sums for denominators
    //if(env_.print_details()) printf("rank=%d Computing sums for denominators\n",rank);
    //const int compute_sums_this = CEnv_is_ppc64() ? 1 : 2;
    const int compute_sums_this = 0; //FIX env_.is_threshold_tc() ? 0 : 1;

    if (0 == compute_sums_this) { // needed here for is_threshold_tc
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

    // Commence numerators computation
    //if(env_.print_details()) printf("rank=%d Computing numerators\n",rank);
    if (vars.is_compute_step && vars.do_compute_block) {
      lock(lock_vectors_left_buf_d);
      if (! vars.is_right_aliased) {
        lock(lock_vectors_right_buf_d);
      }
      lock(lock_metrics_buf_ptr_d);
      if(env_.print_details()) printf("rank=%d Calling ComputeMetrics2WayBlock::compute_nums_start\n",rank);
      ComputeMetrics2WayBlock::compute_nums_start(
        vectors_left, vars.vectors_right, &metrics,
        vectors_left_buf, vars.vectors_right_buf, vars.metrics_buf,
        vector_sums_left, vars.vector_sums_right,
        vars.j_block, vars.is_main_diag, magma_wrapper, &env_);
      if(env_.print_details()) printf("rank=%d Done calling ComputeMetrics2WayBlock::compute_nums_start\n",rank);
    }

    // GPU case: wait for prev step get metrics to complete, then combine.
    // Note this is hidden under GPU computation
    //if(env_.print_details()) printf("rank=%d waiting for prev step get metrics to complete\n",rank);
    if (env_.is_using_linalg()) {
      if (vars_prev.is_compute_step && vars_prev.do_compute_block) {
        //vars_prev.metrics_buf->from_accel_wait();
//>>>
        matB_buf_compressed.from_accel_wait();
        unlock(lock_metrics_buf_ptr_d_prev);
        unlock(lock_metrics_buf_ptr_h_prev);
        lock(lock_metrics_buf_ptr_h_prev);
        gm_metrics_pad_adjust(&metrics, vars_prev.metrics_buf, &env_);
//>>>
        unlock(lock_metrics_buf_ptr_h_prev);

        //TODO: remove need to allocate metrics_tmp_buf device array
        MirroredBuf* metrics_buf_prev_ptr =
            env_.do_reduce() ?  &metrics_tmp_buf_ : vars_prev.metrics_buf;
//>>>

        lock(lock_metrics_buf_ptr_h_prev); // semantics not perfect but ok

        if (env_.do_reduce()) {
          lock(lock_metrics_tmp_buf_h);
          gm_reduce_metrics(&metrics, metrics_buf_prev_ptr,
                            vars_prev.metrics_buf, &env_);
//>>>
          matB_buf_compressed.attach(*metrics_buf_prev_ptr);
        }

        if(env_.print_details()) printf("rank=%d Calling ComputeMetrics2WayBlock::finalize\n",rank);
        ComputeMetrics2WayBlock::finalize(
          &metrics,
          //metrics_buf_prev_ptr,
//>>>
        &matB_buf_compressed, 
          vector_sums_left, vars_prev.vector_sums_right,
          vars_prev.j_block,
          vars_prev.is_main_diag, &env_);
        if(env_.print_details()) printf("Done calling ComputeMetrics2WayBlock::finalize\n");

        unlock(lock_metrics_buf_ptr_h_prev); // semantics not perfect but ok

        if (env_.do_reduce()) {
          unlock(lock_metrics_tmp_buf_h);
        }
      }
    }

    // ISSUE: it may be possible to increase performance by swapping the
    // some code blocks below and the one code block above.  It depends
    // on the relative speeds.  If these would be put in two different
    // CPU threads, then it wouldn't matter.

    // Compute sums for denominators
    //if(env_.print_details()) printf("rank=%d computing sums for denominators\n",rank);
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

    // Wait for recvs to complete
    //if(env_.print_details()) printf("rank=%d Waiting for recvs to complete\n",rank);
    if (vars_next.is_compute_step && ! comm_with_self) {
      gm_recv_vectors_wait(&(mpi_requests[1]), &env_);
      COMET_INSIST((!vars_next.is_right_aliased) &&
               "Next step should always compute off-diag block.");
      unlock(lock_vectors_right_buf_h_next);
    }

    // Send right vectors for next step to GPU start

    if (vars_next.is_compute_step && vars_next.do_compute_block &&
        ! vars_next.is_right_aliased) {
      // ISSUE: make sure not necessary if vars_next.is_right_aliased
      lock(lock_vectors_right_buf_h_next);
      lock(lock_vectors_right_buf_d_next);
      vars_next.vectors_right_buf->to_accel_start();
    }

    // Wait for numerators computation to complete

    if (vars.is_compute_step && vars.do_compute_block) {
      if(env_.print_details()) printf("rank=%d Calling compute_nums_wait\n",rank);
      ComputeMetrics2WayBlock::compute_nums_wait(
        vectors_left, vars.vectors_right, &metrics,
        vectors_left_buf, vars.vectors_right_buf, vars.metrics_buf,
//>>>
        vector_sums_left, vars.vector_sums_right,
        vars.j_block, vars.is_main_diag, &env_);

      if(env_.print_details()) printf("rank=%d Done calling compute_nums_wait\n",rank);
        matB_buf_compressed.attach(*vars.metrics_buf);
        matB_buf_compressed.compress();
      unlock(lock_vectors_left_buf_d);
      if (! vars.is_right_aliased) {
        unlock(lock_vectors_right_buf_d);
      }
      unlock(lock_metrics_buf_ptr_d);
    }

    // Commence copy of completed numerators back from GPU

    if (vars.is_compute_step && vars.do_compute_block) {
      lock(lock_metrics_buf_ptr_h);
      lock(lock_metrics_buf_ptr_d);
      //vars.metrics_buf->from_accel_start();
//>>>
      matB_buf_compressed.from_accel_start();
    }

    // Compute sums for denominators

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

    // CPU case: combine numerators, denominators to obtain final result
    //if(env_.print_details()) printf("rank=%d Combining numerators and denominators\n",rank);
    if (!env_.is_using_linalg()) {
      if (vars.is_compute_step && vars.do_compute_block) {
        //vars.metrics_buf->from_accel_wait(); // NO-OP
//>>>
        matB_buf_compressed.from_accel_wait();
        unlock(lock_metrics_buf_ptr_d);
        unlock(lock_metrics_buf_ptr_h);
        lock(lock_metrics_buf_ptr_h);
        if(env_.print_details()) printf("rank=%d Calling 2nd Block::finalize\n",rank);
        ComputeMetrics2WayBlock::finalize(
          &metrics,
          //vars.metrics_buf,
//>>>
          &matB_buf_compressed,
          vector_sums_left,
          vars.vector_sums_right, vars.j_block,
          vars.is_main_diag, &env_);
        if(env_.print_details()) printf("rank=%d Done calling 2nd Block::finalize\n",rank);
        unlock(lock_metrics_buf_ptr_h);
      }
    }

    // Wait for sends to complete

    if (vars_next.is_compute_step && ! comm_with_self) {
      gm_send_vectors_wait(&(mpi_requests[0]), &env_);
    }

    //if(env_.print_details()) printf("rank=%d Done with loop step_num=%d\n",rank,step_num);

  //========================================
  } // step_num
  //========================================

  //---------------
  // Terminations
  //---------------

  for (int i=0; i<2; ++i) {
    COMET_INSIST(!lock_vectors_01_buf_h[i]);
    COMET_INSIST(!lock_vectors_01_buf_d[i]);
    COMET_INSIST(!lock_metrics_buf_01_h[i]);
    COMET_INSIST(!lock_metrics_buf_01_d[i]);
  }
  COMET_INSIST(!lock_vectors_buf_h);
  COMET_INSIST(!lock_vectors_buf_d);
  COMET_INSIST(!lock_metrics_tmp_buf_h);

  //MagmaWrapper::finalize(env_);
  if(env_.print_details()) printf("rank=%d Done in compute_all2all\n",rank);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

