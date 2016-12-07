/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_2way.c
 * \author Wayne Joubert
 * \date   Thu Jan  7 10:21:09 EST 2016
 * \brief  Functions for computing 2-way metrics.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include "env.h"
#include "vector_sums.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_utils_magma.h"
#include "compute_metrics_utils.h"
#include "compute_metrics_2way.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

void gm_compute_metrics_2way_notall2all(GMMetrics* metrics,
                                        GMVectors* vectors,
                                        GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);

  GMAssertAlways(!Env_all2all(env));

  /*---------------*/
  /*---Denominator---*/
  /*---------------*/

  GMVectorSums vector_sums = GMVectorSums_null();
  GMVectorSums_create(&vector_sums, vectors, env);

  /* .02 / 1.56 */
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

  GMMirroredPointer metrics_buf_tmp;
  if (Env_num_proc_field(env) > 1) {
    metrics_buf_tmp = gm_malloc_magma(numvecl * (size_t)numvecl, env);
  }

  GMMirroredPointer* metrics_buf_local =
      Env_num_proc_field(env) > 1 ? &metrics_buf_tmp : &metrics_buf;

  /*---Copy in vectors---*/

  /* .08 / 1.56 */
  gm_vectors_to_buf(&vectors_buf, vectors, env);

  /*---Send vectors to GPU---*/

  gm_set_vectors_start(vectors, &vectors_buf, env);
  gm_set_vectors_wait(env);

  gm_compute_numerators_2way_start(vectors, vectors, metrics, &vectors_buf,
                                   &vectors_buf, metrics_buf_local,
                                   Env_proc_num_vector_i(env),
                                   GM_BOOL_TRUE, env);
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

  /* .22 / 1.56 */
  gm_compute_2way_combine(metrics, &metrics_buf, &vector_sums, &vector_sums,
                          Env_proc_num_vector_i(env), GM_BOOL_TRUE, env);

  /*---------------*/
  /*---Free memory---*/
  /*---------------*/

  GMVectorSums_destroy(&vector_sums, env);

  gm_free_magma(&vectors_buf, env);
  gm_free_magma(&metrics_buf, env);
  if (Env_num_proc_field(env) > 1) {
    gm_free_magma(&metrics_buf_tmp, env);
  }

  gm_magma_finalize(env);
}

/*===========================================================================*/

void gm_compute_metrics_2way_all2all(GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);

  GMAssertAlways(Env_all2all(env));

  /*---Initializations---*/

  const int num_block = Env_num_block_vector(env);

  const int numvecl = vectors->num_vector_local;
  const int numpfieldl = vectors->num_packedval_field_local;

  int i = 0;
  int i_block = Env_proc_num_vector_i(env);

  GMVectorSums vector_sums_onproc = GMVectorSums_null();
  GMVectorSums vector_sums_offproc = GMVectorSums_null();
  GMVectorSums_create(&vector_sums_onproc, vectors, env);
  GMVectorSums_create(&vector_sums_offproc, vectors, env);

  /*---Create double buffer of vectors objects for send/recv---*/

  GMVectors vectors_01[2];
  for (i = 0; i < 2; ++i) {
    GMVectors_create(&vectors_01[i], Env_data_type_vectors(env),
                     vectors->num_field, numvecl, env);
  }

  /*---Magma initializations---*/

  gm_magma_initialize(env);

  /*---Allocate GPU buffers---*/
  /*---To overlap transfers with compute, set up double buffers for the
       vectors sent to the GPU and the metrics received from the GPU.---*/

  GMMirroredPointer metrics_buf_01[2];
  GMMirroredPointer vectors_buf_01[2];
  GMMirroredPointer vectors_buf = GMMirroredPointer_null();
  GMMirroredPointer metrics_buf_tmp = GMMirroredPointer_null();
  for (i = 0; i < 2; ++i) {
    vectors_buf_01[i] = gm_malloc_magma(numvecl * (size_t)numpfieldl, env);
    metrics_buf_01[i] = gm_malloc_magma(numvecl * (size_t)numvecl, env);
  }
  vectors_buf = gm_malloc_magma(numvecl * (size_t)numpfieldl, env);
  metrics_buf_tmp = gm_malloc_magma(numvecl * (size_t)numvecl, env);

  /*---Result matrix is diagonal block and half the blocks to the right
       (including wraparound to left side of matrix when appropriate).
       For even number of vector blocks, block rows of lower half of matrix
       have one less block to make correct count---*/

  const int num_proc_r = Env_num_proc_repl(env);
  const int proc_num_r = Env_proc_num_repl(env);

  const int proc_num_ir = proc_num_r + num_proc_r * i_block;
  const int num_proc_ir = num_block * num_proc_r;

  MPI_Request mpi_requests[2];

  /*----------------------------------------*/
  /*---Begin loop over steps of circular shift of vectors objects---*/
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

//XXXCHANGE
  const int max_rectangle_width = 1 + (num_block / 2);

  const int rectangle_width =
    (num_block % 2 == 0) && (2 * i_block >= num_block) ?
    max_rectangle_width - 1 : max_rectangle_width;

  const int num_step = gm_ceil_i(max_rectangle_width, num_proc_r);

  int step_num = 0;

  /*========================================*/
  for (step_num = 0-extra_step; step_num < num_step+extra_step; ++step_num) {
  /*========================================*/
    /*---Determine what kind of step this is---*/

    const _Bool is_compute_step = step_num >= 0 && step_num < num_step;
    const _Bool is_compute_step_prev =
        step_num - 1 >= 0 && step_num - 1 < num_step;
    const _Bool is_compute_step_next =
        step_num + 1 >= 0 && step_num + 1 < num_step;

    const _Bool is_first_compute_step = step_num == 0;
    const _Bool is_first_compute_step_prev = step_num - 1 == 0;
    //const _Bool is_first_compute_step_next = step_num + 1 == 0;

    /*---Which entry of double buffered data items to use---*/

    const int index_01 = (step_num + 2) % 2;
    const int index_01_prev = (step_num - 1 + 2) % 2;
    const int index_01_next = (step_num + 1 + 2) % 2;

    /*---Main diagonal block only computes strict upper triangular part---*/

    const _Bool do_compute_triang_only = is_first_compute_step &&
        proc_num_r == 0;

    const _Bool do_compute_triang_only_prev = is_first_compute_step_prev &&
        proc_num_r == 0;

    //const _Bool do_compute_triang_only_next = is_first_compute_step_next &&
    //    proc_num_r == 0;

    /*---Possibly skip the block computation on the final compute step---*/

//XXXCHANGE
    const int proc_num_offset = num_proc_r * ( proc_num_r + num_proc_r *
                                             ( step_num ) );

    const int proc_num_offset_prev = num_proc_r * ( proc_num_r + num_proc_r *
                                             ( step_num - 1 ) );

    const int proc_num_offset_next = num_proc_r * ( proc_num_r + num_proc_r *
                                             ( step_num + 1 ) );

    const _Bool do_compute_block = is_compute_step &&
      proc_num_offset/num_proc_r < rectangle_width;

    const _Bool do_compute_block_prev = is_compute_step_prev &&
      proc_num_offset_prev/num_proc_r < rectangle_width;

    const _Bool do_compute_block_next = is_compute_step_next &&
      proc_num_offset_next/num_proc_r < rectangle_width;

    /*---Point to left/right-side vecs, also right-side vecs for next step.
         Here we are computing V^T W, for V, W containing column vectors---*/

    GMVectors* vectors_left = vectors;
    GMVectors* vectors_right =
        do_compute_triang_only ? vectors : &vectors_01[index_01];
    GMVectors* vectors_right_next = &vectors_01[index_01_next];

    GMMirroredPointer* vectors_left_buf = &vectors_buf;
    GMMirroredPointer* vectors_right_buf =
        do_compute_triang_only ? &vectors_buf : &vectors_buf_01[index_01];
    GMMirroredPointer* vectors_right_buf_next = &vectors_buf_01[index_01_next];

    /*---Point to metrics buffers---*/

    GMMirroredPointer* metrics_buf = &metrics_buf_01[index_01];
    GMMirroredPointer* metrics_buf_prev = &metrics_buf_01[index_01_prev];

    /*---Prepare for sends/recvs: procs for communication---*/

//XXXCHANGE
    const int proc_recv = gm_mod_i(proc_num_ir + proc_num_offset_next,
                                   num_proc_ir);
    const int proc_send = gm_mod_i(proc_num_ir - proc_num_offset_next,
                                   num_proc_ir);

    const _Bool comm_with_self = proc_num_offset_next % num_proc_ir == 0;

    /*---Initiate sends/recvs for vecs needed on next step---*/

    if (is_compute_step_next && !comm_with_self) {
      const int mpi_tag = step_num + 1;
      mpi_requests[0] = gm_send_vectors_start(vectors_left, proc_send,
                                              mpi_tag, env);
      mpi_requests[1] = gm_recv_vectors_start(vectors_right_next, proc_recv,
                                              mpi_tag, env);
    }

    /*---First step: send (left) vecs to GPU---*/

    if (is_first_compute_step) {
      gm_vectors_to_buf(vectors_left_buf, vectors_left, env);
      gm_set_vectors_start(vectors_left, vectors_left_buf, env);
      gm_set_vectors_wait(env);
    }

    /*---The block num for the "right-side" vecs for the pseudo-product---*/

//XXXCHANGE

    const int j_block = gm_mod_i(proc_num_ir + proc_num_offset,
                                 num_proc_ir)/num_proc_r;

    const int j_block_prev = gm_mod_i(proc_num_ir + proc_num_offset_prev,
                                 num_proc_ir)/num_proc_r;

    /*---Send right vectors to GPU end---*/

    if (is_compute_step && do_compute_block) {
      gm_set_vectors_wait(env);
    }

    /*--------------------*/
    /*---Commence numerators computation---*/
    /*--------------------*/

    if (is_compute_step && do_compute_block) {
      gm_compute_numerators_2way_start(
          vectors_left, vectors_right, metrics, vectors_left_buf,
          vectors_right_buf, metrics_buf, j_block, do_compute_triang_only, env);
    }

    /*---GPU case: wait for prev step get metrics to complete, then combine.
         Note this is hidden under GPU computation---*/

    if (Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
      if (is_compute_step_prev && do_compute_block_prev) {
        gm_get_metrics_wait(metrics, metrics_buf_prev, env);

        GMMirroredPointer* metrics_buf_prev_global =
            Env_num_proc_field(env) == 1 ? metrics_buf_prev : &metrics_buf_tmp;

        if (Env_num_proc_field(env) > 1) {
          gm_allreduce_metrics(metrics, metrics_buf_prev_global,
                               metrics_buf_prev, env);
        }

        GMVectorSums* vector_sums_left = &vector_sums_onproc;
//XXXCHANGE
        GMVectorSums* vector_sums_right = proc_num_offset_prev % num_proc_ir == 0
                                              ? &vector_sums_onproc
                                              : &vector_sums_offproc;
        gm_compute_2way_combine(metrics, metrics_buf_prev_global,
                                vector_sums_left, vector_sums_right,
                                j_block_prev, do_compute_triang_only_prev, env);
      }
    }

    /*---ISSUE: it may be possible to increase performance by swapping the
         two code blocks below and the one code block above.  It depends
         on the relative speeds.  If these would be put in two different
         CPU threads, then it wouldn't matter---*/

    /*---Wait for recvs to complete---*/

    if (is_compute_step_next && !comm_with_self) {
      gm_recv_vectors_wait(&(mpi_requests[1]), env);
    }

    /*---Send right vectors for next step to GPU start---*/

    if (is_compute_step_next && do_compute_block_next) {
      gm_vectors_to_buf(vectors_right_buf_next, vectors_right_next, env);
      gm_set_vectors_start(vectors_right_next, vectors_right_buf_next, env);
    }

    /*--------------------*/
    /*---Wait for numerators computation to complete---*/
    /*--------------------*/

    if (is_compute_step && do_compute_block) {
      gm_compute_wait(env);
    }

    /*---Commence copy of completed numerators back from GPU---*/

    if (is_compute_step && do_compute_block) {
      gm_get_metrics_start(metrics, metrics_buf, env);
    }

    /*---Compute sums for denominators---*/

    if (is_compute_step && do_compute_block) {
//XXXCHANGE
      if (is_first_compute_step) {
        GMVectorSums_compute(&vector_sums_onproc, vectors_left, env);
      }
      if (proc_num_offset % num_proc_ir != 0) {
        GMVectorSums_compute(&vector_sums_offproc, vectors_right, env);
      }
    }

    /*---CPU case: combine numerators, denominators to obtain final result---*/

    if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
      if (is_compute_step && do_compute_block) {
        GMVectorSums* vector_sums_left = &vector_sums_onproc;
//XXXCHANGE
        GMVectorSums* vector_sums_right =
            proc_num_offset % num_proc_ir == 0
            ? &vector_sums_onproc : &vector_sums_offproc;
        gm_compute_2way_combine(metrics, metrics_buf, vector_sums_left,
                                vector_sums_right, j_block,
                                do_compute_triang_only, env);
      }
    }

    /*---Wait for sends to complete---*/

    if (is_compute_step_next && !comm_with_self) {
      gm_send_vectors_wait(&(mpi_requests[0]), env);
    }

  /*========================================*/
  } /*---step_num---*/
  /*========================================*/

  /*----------------------------------------*/
  /*---End loop over steps of circular shift of vectors objects---*/
  /*----------------------------------------*/

  GMVectorSums_destroy(&vector_sums_onproc, env);
  GMVectorSums_destroy(&vector_sums_offproc, env);

  for (i = 0; i < 2; ++i) {
    GMVectors_destroy(&vectors_01[i], env);
  }

  /*---Magma terminations---*/

  for (i = 0; i < 2; ++i) {
    gm_free_magma(&metrics_buf_01[i], env);
    gm_free_magma(&vectors_buf_01[i], env);
  }
  gm_free_magma(&vectors_buf, env);
  gm_free_magma(&metrics_buf_tmp, env);

  gm_magma_finalize(env);
}

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
