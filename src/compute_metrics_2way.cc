//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_2way.cc
 * \author Wayne Joubert
 * \date   Thu Jan  7 10:21:09 EST 2016
 * \brief  Calculate metrics, 2-way.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "string.h"

#include "env.hh"
#include "linalg.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "comm_xfer_utils.hh"
#include "compute_metrics_2way_block_nums.hh"
#include "compute_metrics_2way_block_combine.hh"
#include "compute_metrics_2way.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

void ComputeMetrics2Way::create(
    GMDecompMgr* dm,
    GMEnv* env) {
  COMET_INSIST(dm && env);

  if (!(env->num_way() == 2 && env->all2all())) {
    return;
  }

  if (! env->is_proc_active()) {
    return;
  }

  memset((void*)this, 0, sizeof(*this));

  GMVectorSums_create(&this->vector_sums_onproc, dm->num_vector_local, env);
  GMVectorSums_create(&this->vector_sums_offproc, dm->num_vector_local, env);

  for (int i = 0; i < 2; ++i) {
    GMVectors_create_with_buf(&this->vectors_01[i],
                              env->data_type_vectors(), dm, env);
  }

  for (int i = 0; i < 2; ++i) {
    GMMirroredBuf_create(&this->metrics_buf_01[i],
                         dm->num_vector_local, dm->num_vector_local, env);
  }

  GMMirroredBuf_create(&this->vectors_buf, dm->num_packedfield_local,
                       dm->num_vector_local, env);

  if (env->do_reduce()) {
    GMMirroredBuf_create(&this->metrics_tmp_buf,
                         dm->num_vector_local, dm->num_vector_local, env);
  }
}

//-----------------------------------------------------------------------------

void ComputeMetrics2Way::destroy(
    GMEnv* env) {
  COMET_INSIST(env);

  if (!(env->num_way() == 2 && env->all2all())) {
    return;
  }

  if (! env->is_proc_active()) {
    return;
  }

  GMVectorSums_destroy(&this->vector_sums_onproc, env);
  GMVectorSums_destroy(&this->vector_sums_offproc, env);

  for (int i = 0; i < 2; ++i) {
    GMVectors_destroy(&this->vectors_01[i], env);
  }

  for (int i = 0; i < 2; ++i) {
    GMMirroredBuf_destroy(&this->metrics_buf_01[i], env);
  }

  GMMirroredBuf_destroy(&this->vectors_buf, env);

  if (env->do_reduce()) {
    GMMirroredBuf_destroy(&this->metrics_tmp_buf, env);
  }
}


//=============================================================================

void ComputeMetrics2Way::compute_notall2all(
  GMMetrics* metrics,
  GMVectors* vectors,
  GMEnv* env) {

  COMET_INSIST(metrics && vectors && env);
  COMET_INSIST(! env->all2all());

  // Denominator

  GMVectorSums vector_sums = GMVectorSums_null();
  GMVectorSums_create(&vector_sums, vectors->num_vector_local, env);

  GMVectorSums_compute(&vector_sums, vectors, env);

  // Numerator

  gm_linalg_initialize(env);

  const int nvl = vectors->num_vector_local;
  const int npvfl = vectors->num_packedval_field_local;

  // Allocate memory for vectors and for result 

  GMMirroredBuf vectors_buf = GMMirroredBuf_null();
  GMMirroredBuf_create(&vectors_buf, npvfl, nvl, env);

  GMMirroredBuf metrics_buf = GMMirroredBuf_null();
  GMMirroredBuf_create(&metrics_buf, nvl, nvl, env);

  GMMirroredBuf metrics_tmp_buf = GMMirroredBuf_null();
  if (env->do_reduce()) {
    GMMirroredBuf_create(&metrics_tmp_buf, nvl, nvl, env);
  }

  GMMirroredBuf* metrics_buf_ptr =
      env->do_reduce() ?  &metrics_tmp_buf : &metrics_buf;

  // Copy in vectors

  gm_vectors_to_buf(&vectors_buf, vectors, env);

  // Send vectors to GPU

  gm_set_vectors_start(vectors, &vectors_buf, env);
  gm_set_vectors_wait(env);

  gm_compute_2way_proc_nums_start(vectors, vectors, metrics, &vectors_buf,
                                  &vectors_buf, metrics_buf_ptr,
                                  env->proc_num_vector(),
                                  true, env);

  //gm_compute_wait(env);
  gm_compute_2way_proc_nums_wait(vectors, vectors, metrics, &vectors_buf,
                                 &vectors_buf, metrics_buf_ptr,
                                 env->proc_num_vector(),
                                 true, env);

  // Copy result from GPU

  gm_get_metrics_start(metrics, metrics_buf_ptr, env);
  gm_get_metrics_wait(metrics, metrics_buf_ptr, env);
  gm_metrics_pad_adjust(metrics, metrics_buf_ptr, env);

  // Do reduction across field procs if needed

  if (env->do_reduce()) {
    gm_reduce_metrics(metrics, &metrics_buf, metrics_buf_ptr, env);
  }

  // Combine

  gm_compute_2way_proc_combine(metrics, &metrics_buf,
                               &vector_sums, &vector_sums,
                               env->proc_num_vector(), true, env);

  // Terminations

  GMVectorSums_destroy(&vector_sums, env);

  GMMirroredBuf_destroy(&vectors_buf, env);
  GMMirroredBuf_destroy(&metrics_buf, env);

  if (env->do_reduce()) {
    GMMirroredBuf_destroy(&metrics_tmp_buf, env);
  }

  gm_linalg_finalize(env);
}

//=============================================================================

void ComputeMetrics2Way::compute_all2all(
  GMMetrics* metrics,
  GMVectors* vectors,
  GMEnv* env) {

  COMET_INSIST(metrics && vectors && env);
  COMET_INSIST(env->all2all());

  // Initializations

  const int num_block = env->num_block_vector();
  const int i_block = env->proc_num_vector();

  GMVectorSums vector_sums_onproc = this->vector_sums_onproc;
  GMVectorSums vector_sums_offproc = this->vector_sums_offproc;

  gm_linalg_initialize(env);

  // Create double buffer of vectors objects for send/recv

  GMVectors* const & vectors_01 = this->vectors_01;

  // Allocate GPU buffers
  // To overlap transfers with compute, set up double buffers for the
  // vectors sent to the GPU and the metrics received from the GPU.

  GMMirroredBuf* const & metrics_buf_01 = this->metrics_buf_01;
  GMMirroredBuf& vectors_buf = this->vectors_buf;
  GMMirroredBuf& metrics_tmp_buf = this->metrics_tmp_buf;

  // Result matrix is diagonal block and half the blocks to the right
  // (including wraparound to left side of matrix when appropriate).
  // For even number of vector blocks, block rows of lower half of matrix
  //  have one less block to make correct count.

  const int num_proc_r = env->num_proc_repl();
  const int proc_num_r = env->proc_num_repl();

  // Flatten the proc_vector and proc_repl indices into a single index.

  const int num_proc_rv = num_block * num_proc_r;
  const int proc_num_rv = proc_num_r + num_proc_r * i_block;

  MPI_Request mpi_requests[2];

  // Prepare for loop over blocks of result.

  /*---Summary of the opertions in this loop:

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

  ---*/

  // Add extra step at begin/end to fill/drain pipeline.

  const int extra_step = 1;

  // Lowest/highest (block) diag to be computed for this phase,
  // measured from (block) main diag.
  // For all repl procs.

  const int j_i_offset_min = gm_bdiag_computed_min(env);
  const int j_i_offset_max = gm_bdiag_computed_max(env);
  const int j_i_offset_this_row_max = gm_block_computed_this_row_max(env);

  const int num_bdiag_computed = j_i_offset_max - j_i_offset_min;

  // Num steps to take to compute blocks
  // (note: at each step, num_proc_r processors each compute a block)
  // NOTE: num_step should be consistent within same proc_r.

  const int num_step = utils::ceil(num_bdiag_computed, num_proc_r);

  typedef struct {
    GMVectors* vectors_right;
    GMMirroredBuf* vectors_right_buf;
    GMMirroredBuf* metrics_buf;
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

  LoopVars vars = {0};
  LoopVars vars_prev = {0};
  LoopVars vars_next = {0};

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

  //========================================
  for (int step_num = 0-extra_step; step_num < num_step+extra_step; ++step_num){
  //========================================

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

    GMVectors* vectors_left = vectors;
    vars_next.vectors_right = vars_next.is_right_aliased ?
      vectors_left : &vectors_01[vars_next.index_01];

    GMMirroredBuf* vectors_left_buf = &vectors_buf;
    vars_next.vectors_right_buf = vars_next.is_right_aliased ?
      vectors_left_buf : &vectors_01[vars_next.index_01].buf;

    // Pointer to metrics buffer

    vars_next.metrics_buf = &metrics_buf_01[vars_next.index_01];

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

    if (vars_next.is_compute_step && ! comm_with_self) {
      const int mpi_tag = step_num + 1;
      // NOTE: the following order helps performance
      COMET_INSIST((!vars_next.is_right_aliased) &&
               "Next step should always compute off-diag block.");
      lock(lock_vectors_right_buf_h_next);
      mpi_requests[1] = gm_recv_vectors_start(vars_next.vectors_right,
                                              proc_recv, mpi_tag, env);
      mpi_requests[0] = gm_send_vectors_start(vectors_left,
                                              proc_send, mpi_tag, env);
    }

    // Send right vectors to GPU end

    if (vars.is_compute_step && vars.do_compute_block &&
        ! vars.is_right_aliased) {
      gm_set_vectors_wait(env);
      unlock(lock_vectors_right_buf_h);
      unlock(lock_vectors_right_buf_d);
    }

    // First step (for any repl or phase): send (left) vecs to GPU

    if (vars_next.is_first_compute_step) {
      lock(lock_vectors_left_buf_h);
      gm_vectors_to_buf(vectors_left_buf, vectors_left, env);
      lock(lock_vectors_left_buf_d);
      gm_set_vectors_start(vectors_left, vectors_left_buf, env);
      // TODO: examine whether overlap possible.
      // May not be possible for general repl and phase (??).
      gm_set_vectors_wait(env);
      unlock(lock_vectors_left_buf_h);
      unlock(lock_vectors_left_buf_d);
    }

    // Commence numerators computation

    if (vars.is_compute_step && vars.do_compute_block) {
      lock(lock_vectors_left_buf_d);
      if (! vars.is_right_aliased) {
        lock(lock_vectors_right_buf_d);
      }
      lock(lock_metrics_buf_ptr_d);
      gm_compute_2way_proc_nums_start(
        vectors_left, vars.vectors_right, metrics,
        vectors_left_buf, vars.vectors_right_buf, vars.metrics_buf,
        vars.j_block, vars.is_main_diag, env);
    }

    // GPU case: wait for prev step get metrics to complete, then combine.
    // Note this is hidden under GPU computation

    if (env->is_using_linalg()) {
      if (vars_prev.is_compute_step && vars_prev.do_compute_block) {
        gm_get_metrics_wait(metrics, vars_prev.metrics_buf, env);
        unlock(lock_metrics_buf_ptr_d_prev);
        unlock(lock_metrics_buf_ptr_h_prev);
        lock(lock_metrics_buf_ptr_h_prev);
        gm_metrics_pad_adjust(metrics, vars_prev.metrics_buf, env);
        unlock(lock_metrics_buf_ptr_h_prev);

        GMVectorSums* vector_sums_left = &vector_sums_onproc;
        GMVectorSums* vector_sums_right =
          vars_prev.is_main_diag
          ? &vector_sums_onproc : &vector_sums_offproc;

        //TODO: remove need to allocate metrics_tmp_buf device array
        GMMirroredBuf* metrics_buf_prev_ptr =
            env->do_reduce() ?  &metrics_tmp_buf : vars_prev.metrics_buf;

        lock(lock_metrics_buf_ptr_h_prev); // semantics not perfect but ok

        if (env->do_reduce()) {
          lock(lock_metrics_tmp_buf_h);
          gm_reduce_metrics(metrics, metrics_buf_prev_ptr,
                            vars_prev.metrics_buf, env);
        }

        gm_compute_2way_proc_combine(
          metrics, metrics_buf_prev_ptr,
          vector_sums_left, vector_sums_right,
          vars_prev.j_block,
          vars_prev.is_main_diag, env);

        unlock(lock_metrics_buf_ptr_h_prev); // semantics not perfect but ok

        if (env->do_reduce()) {
          unlock(lock_metrics_tmp_buf_h);
        }
      }
    }

    // ISSUE: it may be possible to increase performance by swapping the
    // some code blocks below and the one code block above.  It depends
    // on the relative speeds.  If these would be put in two different
    // CPU threads, then it wouldn't matter.

    // Compute sums for denominators

    //const bool compute_sums_early = GMEnv_is_ppc64();
    const bool compute_sums_early = true;

    if (compute_sums_early) { // put it here for speed on this arch
      if (vars.is_compute_step && vars.do_compute_block) {
        //TODO: possibly move this
        if (vars.is_first_compute_step) {
          GMVectorSums_compute(&vector_sums_onproc, vectors_left, env);
        }
        if (! vars.is_main_diag) {
          GMVectorSums_compute(&vector_sums_offproc, vars.vectors_right, env);
        }
      }
    }

    // Wait for recvs to complete

    if (vars_next.is_compute_step && ! comm_with_self) {
      gm_recv_vectors_wait(&(mpi_requests[1]), env);
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
      gm_set_vectors_start(vars_next.vectors_right,
                           vars_next.vectors_right_buf, env);
    }

    // Wait for numerators computation to complete

    if (vars.is_compute_step && vars.do_compute_block) {
      //gm_compute_wait(env);
      gm_compute_2way_proc_nums_wait(
        vectors_left, vars.vectors_right, metrics,
        vectors_left_buf, vars.vectors_right_buf, vars.metrics_buf,
        vars.j_block, vars.is_main_diag, env);
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
      gm_get_metrics_start(metrics, vars.metrics_buf, env);
    }

    // Compute sums for denominators

    if (! compute_sums_early) { // put it here for speed on this arch
      if (vars.is_compute_step && vars.do_compute_block) {
        //TODO: possibly move this
        if (vars.is_first_compute_step) {
          GMVectorSums_compute(&vector_sums_onproc, vectors_left, env);
        }
        if (! vars.is_main_diag) {
          GMVectorSums_compute(&vector_sums_offproc, vars.vectors_right, env);
        }
      }
    }

    // CPU case: combine numerators, denominators to obtain final result

    if (!env->is_using_linalg()) {
      if (vars.is_compute_step && vars.do_compute_block) {
        GMVectorSums* vector_sums_left = &vector_sums_onproc;
        GMVectorSums* vector_sums_right =
            vars.is_main_diag
            ? &vector_sums_onproc : &vector_sums_offproc;
        gm_get_metrics_wait(metrics, vars.metrics_buf, env); // NO-OP
        unlock(lock_metrics_buf_ptr_d);
        unlock(lock_metrics_buf_ptr_h);
        lock(lock_metrics_buf_ptr_h);
        gm_compute_2way_proc_combine(
          metrics, vars.metrics_buf, vector_sums_left,
          vector_sums_right, vars.j_block,
          vars.is_main_diag, env);
        unlock(lock_metrics_buf_ptr_h);
      }
    }

    // Wait for sends to complete

    if (vars_next.is_compute_step && ! comm_with_self) {
      gm_send_vectors_wait(&(mpi_requests[0]), env);
    }

  //========================================
  } // step_num
  //========================================

  // Terminations

  for (int i=0; i<2; ++i) {
    COMET_INSIST(!lock_vectors_01_buf_h[i]);
    COMET_INSIST(!lock_vectors_01_buf_d[i]);
    COMET_INSIST(!lock_metrics_buf_01_h[i]);
    COMET_INSIST(!lock_metrics_buf_01_d[i]);
  }
  COMET_INSIST(!lock_vectors_buf_h);
  COMET_INSIST(!lock_vectors_buf_d);
  COMET_INSIST(!lock_metrics_tmp_buf_h);

  gm_linalg_finalize(env);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

