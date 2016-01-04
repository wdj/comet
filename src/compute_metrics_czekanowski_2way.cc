/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_czekanowski_2way.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Functions for computing 2-way Czekanowski metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h> /*FIX*/
#include <stdlib.h>

#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"

#include "env.h"
#include "vector_sums.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_metrics_utils.h"
#include "compute_metrics_czekanowski_2way.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

/*---NOTE: This routine currently only handles Czekanowski, but the
    intent is that it will later be adapted to handle all the metrics---*/

void gm_compute_metrics_czekanowski_2way_all2all(GMMetrics* metrics,
                                                 GMVectors* vectors,
                                                 GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  /*---Initializations---*/

  const int num_proc = Env_num_proc_vector(env);

  const int numvec = metrics->num_vector_local;
  const int numfield = vectors->num_field_local;

  int i = 0;
  int i_proc = Env_proc_num_vector(env);

  GMVectorSums vector_sums_onproc = GMVectorSums_null();
  GMVectorSums vector_sums_offproc = GMVectorSums_null();
  GMVectorSums_create(&vector_sums_onproc, vectors, env);
  GMVectorSums_create(&vector_sums_offproc, vectors, env);

  /*---Create double buffer of vectors objects for send/recv---*/

  GMVectors vectors_01[2];
  for (i = 0; i < 2; ++i) {
    const int data_type = Env_data_type_vectors(env);
    GMVectors_create(&vectors_01[i], data_type, vectors->num_field,
      numvec, env);
  }

  /*---Magma initializations---*/

  gm_magma_initialize(env);

  /*---Allocate GPU buffers---*/
  /*---To overlap transfers with compute, set up double buffers for the
       vectors sent to the GPU and the metrics received from the GPU.---*/

  GMMirroredPointer metrics_buf_01[2];
  GMMirroredPointer vectors_buf_01[2];
  GMMirroredPointer vectors_buf = GMMirroredPointer_null();
  for (i = 0; i < 2; ++i) {
    vectors_buf_01[i] = gm_malloc_magma(
        numvec * (size_t)numfield, env);
    metrics_buf_01[i] = gm_malloc_magma(
        numvec * (size_t)numvec, env);
  }
  vectors_buf = gm_malloc_magma(
      numvec * (size_t)numfield, env);
  GMMirroredPointer metrics_buf_tmp = gm_malloc_magma(
        numvec * (size_t)numvec, env);

  /*---Result matrix is diagonal block and half the blocks to the right
       (including wraparound to left side of matrix when appropriate).
       For even number of procs, block rows of lower half of matrix
       have one less block to make correct count---*/

  const int num_step = 1 + (num_proc / 2);

  /*----------------------------------------*/
  /*---Begin loop over steps of circular shift of vectors objects---*/
  /*----------------------------------------*/

  /*---Add extra step at end to drain pipeline---*/
  const int extra_step = 1;

  int step_num = 0;
  for (step_num = 0; step_num < num_step + extra_step; ++step_num) {
    /*---Determine what kind of step this is---*/

    const _Bool is_compute_step = step_num >= 0 && step_num < num_step;
    const _Bool is_compute_step_prev =
        step_num - 1 >= 0 && step_num - 1 < num_step;
    const _Bool is_compute_step_next =
        step_num + 1 >= 0 && step_num + 1 < num_step;

    const _Bool is_first_compute_step = step_num == 0;
    const _Bool is_first_compute_step_prev = step_num - 1 == 0;

    const _Bool is_last_compute_step = step_num == num_step - 1;
    const _Bool is_last_compute_step_prev = step_num - 1 == num_step - 1;
    const _Bool is_last_compute_step_next = step_num + 1 == num_step - 1;

    /*---Which entry of double buffered data items to use---*/

    const int index_01 = (step_num + 2) % 2;
    const int index_01_prev = (step_num - 1 + 2) % 2;
    const int index_01_next = (step_num + 1 + 2) % 2;

    /*---Point to left/right-side vecs, also right-side vecs for next step.
         Here we are computing V^T W, for V, W containing column vectors---*/

    GMVectors* vectors_left = vectors;
    GMVectors* vectors_right =
        is_first_compute_step ? vectors : &vectors_01[index_01];
    GMVectors* vectors_right_next = &vectors_01[index_01_next];

    GMMirroredPointer* vectors_left_buf = &vectors_buf;
    GMMirroredPointer* vectors_right_buf =
        is_first_compute_step ? &vectors_buf : &vectors_buf_01[index_01];
    GMMirroredPointer* vectors_right_buf_next = &vectors_buf_01[index_01_next];

    /*---Point to metrics buffers---*/

    GMMirroredPointer* metrics_buf = &metrics_buf_01[index_01];
    GMMirroredPointer* metrics_buf_prev = &metrics_buf_01[index_01_prev];

    /*---Prepare for sends/recvs: procs for communication---*/

    const int proc_up = (i_proc + 1) % num_proc;
    const int proc_dn = (i_proc - 1 + num_proc) % num_proc;

    MPI_Request mpi_requests[2];

    /*---Initiate sends/recvs for vecs needed on next step---*/

    if (is_compute_step_next) {
      mpi_requests[0] = gm_send_vectors_start(vectors_right, proc_dn, env);
      mpi_requests[1] = gm_recv_vectors_start(vectors_right_next, proc_up, env);
    }

    /*---First step: send (left) vecs to GPU---*/

    if (is_first_compute_step) {
      gm_vectors_to_buf(vectors_left, vectors_left_buf, env);
      gm_set_vectors_start(vectors_left, vectors_left_buf, env);
      gm_set_vectors_wait(env);
    }

    /*---The proc that owns the "right-side" vecs for the minproduct---*/

    const int j_proc = (i_proc + step_num) % num_proc;
    const int j_proc_prev = (i_proc + step_num - 1) % num_proc;

    /*---To remove redundancies from symmetry, skip some blocks---*/

    const _Bool skipping_active = (num_proc % 2 == 0) &&
                                  (2 * i_proc >= num_proc);

    const _Bool skipped_last_block_lower_half =
        skipping_active && is_last_compute_step;

    const _Bool skipped_last_block_lower_half_prev =
        skipping_active && is_last_compute_step_prev;

    const _Bool skipped_last_block_lower_half_next =
        skipping_active && is_last_compute_step_next;

    const _Bool do_compute_block =
        is_compute_step && !skipped_last_block_lower_half;

    const _Bool do_compute_block_prev =
        is_compute_step_prev && !skipped_last_block_lower_half_prev;

    const _Bool do_compute_block_next =
        is_compute_step_next && !skipped_last_block_lower_half_next;

    /*---Main diagonal block only computes strict upper triangular part---*/

    const _Bool do_compute_triang_only = is_first_compute_step;
    const _Bool do_compute_triang_only_prev = is_first_compute_step_prev;

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
          vectors_right_buf, metrics_buf, j_proc, do_compute_triang_only, env);
    }

    /*---GPU case: wait for prev step get metrics to complete, then combine.
         Note this is hidden under GPU computation---*/

    if (Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
      if (is_compute_step_prev && do_compute_block_prev) {
        gm_get_metrics_wait(env);

        GMMirroredPointer* metrics_buf_prev_global =
          Env_num_proc_field(env) == 1 ? metrics_buf_prev : &metrics_buf_tmp;

        if (Env_num_proc_field(env) > 1) {
          int mpi_code = 0;
          mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
          mpi_code = MPI_Allreduce(metrics_buf_prev->h,
                       metrics_buf_prev_global->h,
                       numvec*(size_t)numvec, GM_MPI_FLOAT, MPI_SUM,
                       Env_mpi_comm_field(env));
          GMAssert(mpi_code == MPI_SUCCESS);
        }


        GMVectorSums* vector_sums_left = &vector_sums_onproc;
        GMVectorSums* vector_sums_right = is_first_compute_step_prev
                                              ? &vector_sums_onproc
                                              : &vector_sums_offproc;
        gm_compute_2way_combine(metrics, metrics_buf_prev_global,
                                vector_sums_left, vector_sums_right,
                                j_proc_prev, do_compute_triang_only_prev, env);
      }
    }

    /*---ISSUE: it may be possible to increase performance by swapping the
         two code blocks below and the one code block above.  It depends
         on the relative speeds.  If these would be put in two different
         CPU threads, then it wouldn't matter---*/

    /*---Wait for recvs to complete---*/

    if (is_compute_step_next) {
      gm_recv_vectors_wait(&(mpi_requests[1]), env);
    }

    /*---Send right vectors for next step to GPU start---*/

    if (is_compute_step_next && do_compute_block_next) {
      gm_vectors_to_buf(vectors_right_next, vectors_right_buf_next, env);
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
      if (is_first_compute_step) {
        GMVectorSums_compute(&vector_sums_onproc, vectors_left, env);
      } else {
        GMVectorSums_compute(&vector_sums_offproc, vectors_right, env);
      }
    }

    /*---CPU case: combine numerators, denominators to obtain final result---*/

    if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
      if (is_compute_step && do_compute_block) {
        GMVectorSums* vector_sums_left = &vector_sums_onproc;
        GMVectorSums* vector_sums_right =
            is_first_compute_step ? &vector_sums_onproc : &vector_sums_offproc;
        gm_compute_2way_combine(metrics, metrics_buf, vector_sums_left,
                                vector_sums_right, j_proc,
                                do_compute_triang_only, env);
      }
    }

    /*---Wait for sends to complete---*/

    if (is_compute_step_next) {
      gm_send_vectors_wait(&(mpi_requests[0]), env);
    }

  } /*---step_num---*/

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

/*===========================================================================*/

void gm_compute_metrics_czekanowski_2way_cpu(GMMetrics* metrics,
                                             GMVectors* vectors,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, Env_num_proc_field(env) == 1
    ? "num_proc_field>1 for CPU case not supported" : 0);

  if (Env_all2all(env)) {
    gm_compute_metrics_czekanowski_2way_all2all(metrics, vectors, env);
    return;
  }

  /*---Compute sums for denominators---*/

  GMFloat* vector_sums = GMFloat_malloc(metrics->num_vector_local);
  GMFloat* vector_sums_tmp = GMFloat_malloc(metrics->num_vector_local);

  gm_compute_float_vector_sums(vectors, vector_sums, vector_sums_tmp, env);

  /*---Compute numerators---*/

  int i = 0;
  int j = 0;
  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      GMFloat sum = 0;
      int field_local = 0;
      for (field_local = 0; field_local < vectors->num_field_local;
           ++field_local) {
        const GMFloat value_i = GMVectors_float_get(vectors,
                                                    field_local, i, env);
        const GMFloat value_j = GMVectors_float_get(vectors,
                                                    field_local, j, env);
        sum += value_i < value_j ? value_i : value_j;
      } /*---for k---*/
      GMMetrics_float_set_2(metrics, i, j, sum, env);
    } /*---for j---*/
  }   /*---for i---*/

  /*---Combine---*/

  for (i = 0; i < metrics->num_vector_local; ++i) {
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      const GMFloat numerator = GMMetrics_float_get_2(metrics, i, j, env);
      const GMFloat denominator = vector_sums[i] + vector_sums[j];
      GMMetrics_float_set_2(metrics, i, j, 2 * numerator / denominator, env);
    } /*---for j---*/
  }   /*---for i---*/

  /*---Deallocations---*/

  free(vector_sums);
  free(vector_sums_tmp);
}

/*===========================================================================*/

void gm_compute_metrics_czekanowski_2way_gpu(GMMetrics* metrics,
                                             GMVectors* vectors,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  if (Env_all2all(env)) {
    gm_compute_metrics_czekanowski_2way_all2all(metrics, vectors, env);
    return;
  }

  /*---Denominator---*/

  GMVectorSums vector_sums = GMVectorSums_null();
  GMVectorSums_create(&vector_sums, vectors, env);

  /* .02 / 1.56 */
  GMVectorSums_compute(&vector_sums, vectors, env);

  /*---------------*/
  /*---Numerator---*/
  /*---------------*/

  gm_magma_initialize(env);

  const int numvec = metrics->num_vector_local;
  const int numfield = vectors->num_field_local;

  /*---Allocate magma CPU memory for vectors and for result */

  GMMirroredPointer vectors_buf =
      gm_malloc_magma(numvec * (size_t)numfield, env);

  GMMirroredPointer metrics_buf =
      gm_malloc_magma(numvec * (size_t)numvec, env);

  GMMirroredPointer metrics_buf_tmp =
      gm_malloc_magma(numvec * (size_t)numvec, env);

  GMMirroredPointer* metrics_buf_local = Env_num_proc_field(env) == 1 ?
    &metrics_buf : &metrics_buf_tmp;

  /*---Copy in vectors---*/

  /* .08 / 1.56 */
  gm_vectors_to_buf(vectors, &vectors_buf, env);

  /*---Send vectors to GPU---*/

  gm_set_vectors_start(vectors, &vectors_buf, env);
  gm_set_vectors_wait(env);

  gm_compute_numerators_2way_start(vectors, vectors, metrics, &vectors_buf,
                                   &vectors_buf, metrics_buf_local,
                                   Env_proc_num_vector(env),
                                   GM_BOOL_TRUE, env);
  gm_compute_wait(env);

  /*---Copy result from GPU---*/

  gm_get_metrics_start(metrics, metrics_buf_local, env);
  gm_get_metrics_wait(env);

  if (Env_num_proc_field(env) > 1) {
    int mpi_code = 0;
    mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
    mpi_code = MPI_Allreduce(metrics_buf_local->h, metrics_buf.h,
                 numvec*(size_t)numvec, GM_MPI_FLOAT, MPI_SUM,
                 Env_mpi_comm_field(env));
    GMAssert(mpi_code == MPI_SUCCESS);
  }

  /*---Combine---*/

  /* .22 / 1.56 */
  gm_compute_2way_combine(metrics, &metrics_buf, &vector_sums, &vector_sums,
                          Env_proc_num_vector(env), GM_BOOL_TRUE, env);

  /*---Free memory---*/

  GMVectorSums_destroy(&vector_sums, env);

  gm_free_magma(&vectors_buf, env);
  gm_free_magma(&metrics_buf, env);
  gm_free_magma(&metrics_buf_tmp, env);

  gm_magma_finalize(env);
}

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
