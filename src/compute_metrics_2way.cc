/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_2way.cc
 * \author Wayne Joubert
 * \date   Thu Jan  7 10:21:09 EST 2016
 * \brief  Functions for computing 2-way metrics.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include "env.hh"
#include "vector_sums.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "compute_utils_linalg.hh"
#include "compute_metrics_utils.hh"
#include "compute_metrics_utils_2way.hh"
#include "compute_metrics_2way.hh"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

void gm_compute_metrics_2way_notall2all(GMMetrics* metrics,
                                        GMVectors* vectors,
                                        GMEnv* env) {
  GMAssertAlways(metrics && vectors && env);
  GMAssertAlways(!GMEnv_all2all(env));

  /*---------------*/
  /*---Denominator---*/
  /*---------------*/

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

  GMMirroredPointer vectors_buf = gm_linalg_malloc(npvfl, nvl, env);

  GMMirroredPointer metrics_buf = gm_linalg_malloc(nvl, nvl, env);

  GMMirroredPointer metrics_buf_tmp;
  if (env->do_reduce) {
    metrics_buf_tmp = gm_linalg_malloc(nvl, nvl, env);
  }

  GMMirroredPointer* metrics_buf_ptr =
      env->do_reduce ?  &metrics_buf_tmp : &metrics_buf;

  /*---Copy in vectors---*/

  gm_vectors_to_buf(&vectors_buf, vectors, env);

  /*---Send vectors to GPU---*/

  gm_set_vectors_start(vectors, &vectors_buf, env);
  gm_set_vectors_wait(env);

  gm_compute_numerators_2way_start(vectors, vectors, metrics, &vectors_buf,
                                   &vectors_buf, metrics_buf_ptr,
                                   GMEnv_proc_num_vector_i(env),
                                   true, env);
  gm_compute_wait(env);

  /*---Copy result from GPU---*/

  gm_get_metrics_start(metrics, metrics_buf_ptr, env);
  gm_get_metrics_wait(metrics, metrics_buf_ptr, env);
  gm_metrics_gpu_adjust(metrics, metrics_buf_ptr, env);

  /*---Do reduction across field procs if needed---*/

  if (env->do_reduce) {
    gm_reduce_metrics(metrics, &metrics_buf, metrics_buf_ptr, env);
  }

  /*---------------*/
  /*---Combine---*/
  /*---------------*/

  gm_compute_2way_combine(metrics, &metrics_buf, &vector_sums, &vector_sums,
                          GMEnv_proc_num_vector_i(env), true, env);

  /*---------------*/
  /*---Free memory---*/
  /*---------------*/

  GMVectorSums_destroy(&vector_sums, env);

  gm_linalg_free(&vectors_buf, env);
  gm_linalg_free(&metrics_buf, env);
  if (env->do_reduce) {
    gm_linalg_free(&metrics_buf_tmp, env);
  }

  gm_linalg_finalize(env);
}

/*===========================================================================*/

void gm_compute_metrics_2way_all2all(GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env) {
  GMAssertAlways(metrics && vectors && env);
  GMAssertAlways(GMEnv_all2all(env));

  /*---Initializations---*/

  const int num_block = GMEnv_num_block_vector(env);

  const int nvl = vectors->num_vector_local;
  const int npvfl = vectors->num_packedval_field_local;

  const int i_block = GMEnv_proc_num_vector_i(env);

  GMVectorSums vector_sums_onproc = GMVectorSums_null();
  GMVectorSums vector_sums_offproc = GMVectorSums_null();
  GMVectorSums_create(&vector_sums_onproc, vectors, env);
  GMVectorSums_create(&vector_sums_offproc, vectors, env);

  /*---Magma initializations---*/

  gm_linalg_initialize(env);

  /*---Create double buffer of vectors objects for send/recv---*/

  GMVectors vectors_01[2];
  for (int i = 0; i < 2; ++i) {
    GMVectors_create_with_buf(&vectors_01[i], GMEnv_data_type_vectors(env),
                     vectors->num_field, vectors->num_field_active,
                     nvl, env);
  }

  /*---Allocate GPU buffers---*/
  /*---To overlap transfers with compute, set up double buffers for the
       vectors sent to the GPU and the metrics received from the GPU.---*/

  GMMirroredPointer metrics_buf_01[2];
  GMMirroredPointer vectors_buf = GMMirroredPointer_null();
  GMMirroredPointer metrics_buf_tmp = GMMirroredPointer_null();
  for (int i = 0; i < 2; ++i) {
    metrics_buf_01[i] = gm_linalg_malloc(nvl, nvl, env);
  }
  vectors_buf = gm_linalg_malloc(npvfl, nvl, env);
  metrics_buf_tmp = gm_linalg_malloc(nvl, nvl, env);

  /*---Result matrix is diagonal block and half the blocks to the right
       (including wraparound to left side of matrix when appropriate).
       For even number of vector blocks, block rows of lower half of matrix
       have one less block to make correct count---*/

  const int num_proc_r = GMEnv_num_proc_repl(env);
  const int proc_num_r = GMEnv_proc_num_repl(env);

  /*---Flatten the proc_vector and proc_repl indices into a single index---*/

  const int num_proc_ir = num_block * num_proc_r;
  const int proc_num_ir = proc_num_r + num_proc_r * i_block;

  MPI_Request mpi_requests[2];

  /*----------------------------------------*/
  /*---Begin loop over computing blocks of the result---*/
  /*----------------------------------------*/

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

  /*---Add extra step at begin/end to fill/drain pipeline---*/

  const int extra_step = 1;

  /*---Lowest/highest (block) diag to be computed for this phase,
       measured from (block) main diag---*/

  const int j_i_offset_min = gm_diag_computed_min(env);
  const int j_i_offset_max = gm_diag_computed_max(env);
  const int j_i_offset_this_row_max = gm_diag_computed_this_row_max(env);

  const int num_diag_computed = j_i_offset_max - j_i_offset_min;

  /*---Num steps to take to compute blocks---*/
  /*---(note: at each step, num_proc_r processors each compute a block)---*/

  const int num_step = gm_ceil_i(num_diag_computed, num_proc_r);

  typedef struct {
    GMVectors* vectors_right;
    GMMirroredPointer* vectors_right_buf;
    GMMirroredPointer* metrics_buf;
    bool is_compute_step;
    bool is_first_compute_step;
    bool do_compute_block;
    bool is_main_diag;
    int step_num;
    int index_01;
    int j_i_offset;
    int j_block;
  } LoopVars;

  LoopVars vars = {0};
  LoopVars vars_prev = {0};
  LoopVars vars_next = {0};

  /*========================================*/
  for (int step_num = 0-extra_step; step_num < num_step+extra_step; ++step_num){
  /*========================================*/

    vars_prev = vars;
    vars = vars_next;
    vars_next.step_num = step_num + 1;

    /*---Determine what kind of step this is---*/

    vars_next.is_compute_step = vars_next.step_num >= 0 &&
                                vars_next.step_num < num_step;
    vars_next.is_first_compute_step = vars_next.step_num == 0;

    /*---Which entry of double-buffered data items to use---*/

    vars_next.index_01 = gm_mod_i(vars_next.step_num, 2);

    /*---Offset of the block diagonal to be computed this step and proc_r---*/

    vars_next.j_i_offset = j_i_offset_min + vars_next.step_num * num_proc_r
                           + proc_num_r;

    vars_next.is_main_diag = vars_next.j_i_offset == 0;

    /*---Block num for "right-side" vecs - wrap above offset to num_blocks---*/

    vars_next.j_block = gm_mod_i(i_block + vars_next.j_i_offset, num_block);

    /*---Only compute blocks in rectangle/phase, maybe different per row---*/

    vars_next.do_compute_block = vars_next.is_compute_step &&
                   vars_next.j_i_offset < j_i_offset_this_row_max;

    /*---Pointers to left/right-side vecs.
         Here we are computing V^T W, for V, W containing column vectors---*/

    vars_next.vectors_right = vars_next.is_main_diag ?
      vectors : &vectors_01[vars_next.index_01];
    GMVectors* vectors_left = vectors;

    vars_next.vectors_right_buf = vars_next.is_main_diag ?
      &vectors_buf : &vectors_01[vars_next.index_01].buf;
    GMMirroredPointer* vectors_left_buf = &vectors_buf;

    /*---Pointer to metrics buffer---*/

    vars_next.metrics_buf = &metrics_buf_01[vars_next.index_01];

    /*---Prepare for sends/recvs: procs for communication---*/

    const int proc_send = gm_mod_i(proc_num_ir
        - vars_next.j_i_offset*num_proc_r, num_proc_ir);

    const int proc_recv = gm_mod_i(proc_num_ir
        + vars_next.j_i_offset*num_proc_r, num_proc_ir);

    const bool comm_with_self = vars_next.is_main_diag;

    /*---Initiate sends/recvs for vecs needed on next step---*/

    if (vars_next.is_compute_step && !comm_with_self) {
      const int mpi_tag = step_num + 1;

      /*---NOTE: this order helps performance---*/
      mpi_requests[1] = gm_recv_vectors_start(vars_next.vectors_right,
                                              proc_recv, mpi_tag, env);
      mpi_requests[0] = gm_send_vectors_start(vectors_left,
                                              proc_send, mpi_tag, env);
    }

    /*---First step: send (left) vecs to GPU---*/

    if (vars.is_first_compute_step) {
      gm_vectors_to_buf(vectors_left_buf, vectors_left, env);
      gm_set_vectors_start(vectors_left, vectors_left_buf, env);
      gm_set_vectors_wait(env);
    }

    /*---Send right vectors to GPU end---*/

    if (vars.is_compute_step && vars.do_compute_block) {
      gm_set_vectors_wait(env);
    }

    /*--------------------*/
    /*---Commence numerators computation---*/
    /*--------------------*/

    if (vars.is_compute_step && vars.do_compute_block) {
      gm_compute_numerators_2way_start(
          vectors_left, vars.vectors_right, metrics,
          vectors_left_buf, vars.vectors_right_buf, vars.metrics_buf,
          vars.j_block, vars.is_main_diag, env);
    }

    /*---GPU case: wait for prev step get metrics to complete, then combine.
         Note this is hidden under GPU computation---*/

    if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
      if (vars_prev.is_compute_step && vars_prev.do_compute_block) {
        gm_get_metrics_wait(metrics, vars_prev.metrics_buf, env);
        gm_metrics_gpu_adjust(metrics, vars_prev.metrics_buf, env);

        GMMirroredPointer* metrics_buf_prev_ptr =
            env->do_reduce ?  &metrics_buf_tmp : vars_prev.metrics_buf;

        if (env->do_reduce) {
          gm_reduce_metrics(metrics, metrics_buf_prev_ptr,
                            vars_prev.metrics_buf, env);
        }

        GMVectorSums* vector_sums_left = &vector_sums_onproc;
        GMVectorSums* vector_sums_right =
          vars_prev.is_main_diag
          ? &vector_sums_onproc : &vector_sums_offproc;
        gm_compute_2way_combine(metrics, metrics_buf_prev_ptr,
                                vector_sums_left, vector_sums_right,
                                vars_prev.j_block,
                                vars_prev.is_main_diag, env);
      }
    }

    /*---ISSUE: it may be possible to increase performance by swapping the
         some code blocks below and the one code block above.  It depends
         on the relative speeds.  If these would be put in two different
         CPU threads, then it wouldn't matter---*/

#if 1
    /*---Wait for recvs to complete---*/

    if (vars_next.is_compute_step && !comm_with_self) {
      gm_recv_vectors_wait(&(mpi_requests[1]), env);
    }

    /*---Send right vectors for next step to GPU start---*/

    if (vars_next.is_compute_step && vars_next.do_compute_block) {
      gm_set_vectors_start(vars_next.vectors_right,
                           vars_next.vectors_right_buf, env);
    }

    /*--------------------*/
    /*---Wait for numerators computation to complete---*/
    /*--------------------*/

    if (vars.is_compute_step && vars.do_compute_block) {
      gm_compute_wait(env);
    }

    /*---Commence copy of completed numerators back from GPU---*/

    if (vars.is_compute_step && vars.do_compute_block) {
      gm_get_metrics_start(metrics, vars.metrics_buf, env);
    }

    /*---Compute sums for denominators---*/

    if (vars.is_compute_step && vars.do_compute_block) {
//TODO: possibly move this
      if (vars.is_first_compute_step) {
        GMVectorSums_compute(&vector_sums_onproc, vectors_left, env);
      }
      if (!vars.is_main_diag) {
        GMVectorSums_compute(&vector_sums_offproc, vars.vectors_right, env);
      }
    }

    /*---CPU case: combine numerators, denominators to obtain final result---*/

    if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
      if (vars.is_compute_step && vars.do_compute_block) {
        GMVectorSums* vector_sums_left = &vector_sums_onproc;
        GMVectorSums* vector_sums_right =
            vars.is_main_diag
            ? &vector_sums_onproc : &vector_sums_offproc;
        gm_compute_2way_combine(metrics, vars.metrics_buf, vector_sums_left,
                                vector_sums_right, vars.j_block,
                                vars.is_main_diag, env);
      }
    }

    /*---Wait for sends to complete---*/

    if (vars_next.is_compute_step && !comm_with_self) {
      gm_send_vectors_wait(&(mpi_requests[0]), env);
    }
#endif

#if 0
    /*---Wait for sends to complete---*/
    /*---NOTE: putting this here instead of end of loop seems faster---*/

    if (vars_next.is_compute_step && !comm_with_self) {
      gm_send_vectors_wait(&(mpi_requests[0]), env);
    }

    /*---Wait for recvs to complete---*/

    if (vars_next.is_compute_step && !comm_with_self) {
      gm_recv_vectors_wait(&(mpi_requests[1]), env);
    }

    /*---Compute sums for denominators---*/

    if (vars.is_compute_step && vars.do_compute_block) {
      if (is_first_compute_step) {
        GMVectorSums_compute(&vector_sums_onproc, vectors_left, env);
      }
      if (!vars.is_main_diag) {
        GMVectorSums_compute(&vector_sums_offproc, vars.vectors_right, env);
      }
    }

    /*---Send right vectors for next step to GPU start---*/

    if (vars_next.is_compute_step && vars_next.do_compute_block) {
      //gm_vectors_to_buf(vars_next.vectors_right_buf, vars_next.vectors_right, env);
      gm_set_vectors_start(vars_next.vectors_right, vars_next.vectors_right_buf, env);
    }

    /*--------------------*/
    /*---Wait for numerators computation to complete---*/
    /*--------------------*/

    if (vars.is_compute_step && vars.do_compute_block) {
      gm_compute_wait(env);
    }

    /*---Commence copy of completed numerators back from GPU---*/

    if (vars.is_compute_step && vars.do_compute_block) {
      gm_get_metrics_start(metrics, vars.metrics_buf, env);
    }

    /*---CPU case: combine numerators, denominators to obtain final result---*/

    if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
      if (vars.is_compute_step && vars.do_compute_block) {
        GMVectorSums* vector_sums_left = &vector_sums_onproc;
        GMVectorSums* vector_sums_right =
            vars.is_main_diag
            ? &vector_sums_onproc : &vector_sums_offproc;
        gm_compute_2way_combine(metrics, vars.metrics_buf, vector_sums_left,
                                vector_sums_right, vars.j_block,
                                vars.is_main_diag, env);
      }
    }
#endif

  /*========================================*/
  } /*---step_num---*/
  /*========================================*/

  /*----------------------------------------*/
  /*---End loop over steps of circular shift of vectors objects---*/
  /*----------------------------------------*/

  GMVectorSums_destroy(&vector_sums_onproc, env);
  GMVectorSums_destroy(&vector_sums_offproc, env);

  for (int i = 0; i < 2; ++i) {
    GMVectors_destroy(&vectors_01[i], env);
  }

  /*---Magma terminations---*/

  for (int i = 0; i < 2; ++i) {
    gm_linalg_free(&metrics_buf_01[i], env);
  }
  gm_linalg_free(&vectors_buf, env);
  gm_linalg_free(&metrics_buf_tmp, env);

  gm_linalg_finalize(env);
}

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
