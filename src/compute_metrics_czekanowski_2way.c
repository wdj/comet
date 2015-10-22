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
#include "vectors.h"
#include "metrics.h"
#include "compute_metrics_utils.h"
#include "compute_metrics_czekanowski_2way.h"

/*===========================================================================*/

void gm_compute_metrics_czekanowski_2way_all2all_gpu(GMMetrics* metrics,
                                                     GMVectors* vectors,
                                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  /*---Initializations---*/

  const int data_type = gm_data_type_from_metric_type(env->metric_type, env);
  int i = 0;

  GMFloat* vector_sums_onproc = GMFloat_malloc(metrics->num_vector_local);
  GMFloat* vector_sums_offproc = GMFloat_malloc(metrics->num_vector_local);

  /*---Create double buffer of vectors objects for send/recv---*/

  GMVectors vectors_01[2];
  for ( i = 0 ; i < 2 ; ++i ) {
    GMVectors_create(&vectors_01[i], data_type, vectors->num_field,
                     vectors->num_vector_local, env);
  }

  /*---Magma initializations---*/

  magma_minproduct_int_t magma_code = 0;
  magma_code = magma_code*1; /*---Avoid unused variable warning---*/

  if (env->compute_method == GM_COMPUTE_METHOD_GPU ) {
    magma_code = magma_minproduct_init();
    GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
  }

  /*---Allocate GPU buffers---*/
  /*---To overlap transfers with compute, set up double buffers for the
       vectors sent to the GPU and the metrics received from the GPU.---*/

  GMFloatMirroredPointer metrics_buf_01[2];
  GMFloatMirroredPointer vectors_buf_01[2];
  GMFloatMirroredPointer vectors_buf = GMFloatMirroredPointer_null();
  for ( i = 0 ; i < 2 ; ++i ) {
    vectors_buf_01[i] =
      GMFloat_malloc_magma_minproduct( vectors->num_vector_local * (size_t)
                                       vectors->num_field, env);
    metrics_buf_01[i] =
      GMFloat_malloc_magma_minproduct( metrics->num_vector_local * (size_t)
                                       metrics->num_vector_local, env);
  }
  vectors_buf =
    GMFloat_malloc_magma_minproduct( vectors->num_vector_local * (size_t)
                                     vectors->num_field, env);

  /*---Result matrix is diagonal block and half the blocks to the right
       (including wraparound to left side of matrix when appropriate).
       For even number of procs, block rows of lower half of matrix
       have one less block to make correct count---*/

  const int num_step = 1 + (env->num_proc / 2);

  /*----------------------------------------*/
  /*---Begin loop over steps of circular shift of vectors objects---*/
  /*----------------------------------------*/

  /*---Add extra step at begin and at end to prime/drain pipeline---*/
  const int extra_step = 1;

  int step_num = 0;
  for (step_num = 0; step_num < num_step+extra_step; ++step_num) {

    /*---Determine what kind of step this is---*/

    const GMBool is_compute_step = step_num >= 0 && step_num < num_step;
    const GMBool is_compute_step_prev = step_num-1 >=0 && step_num-1 < num_step;
    const GMBool is_compute_step_next = step_num+1 >=0 && step_num+1 < num_step;

    const GMBool is_first_compute_step = step_num == 0;
    const GMBool is_first_compute_step_prev = step_num-1 == 0;

    const GMBool is_last_compute_step = step_num == num_step-1;
    const GMBool is_last_compute_step_prev = step_num-1 == num_step-1;

    /*---Which entry of double buffered data items to use---*/

    const int index_01 = ( step_num + 2 ) % 2;
    const int index_01_prev = ( step_num-1 + 2 ) % 2;
    const int index_01_next = ( step_num+1 + 2 ) % 2;

    /*---Denote left/right-side vecs, also right-side vecs for next step.
         Here we are computing V^T W, for V, W containing column vectors---*/

    GMVectors* vectors_left = vectors;
    GMVectors* vectors_right = is_first_compute_step ? vectors :
                                                       &vectors_01[index_01];
    GMVectors* vectors_right_next = &vectors_01[index_01_next];

    GMFloatMirroredPointer* vectors_left_buf = &vectors_buf;
    GMFloatMirroredPointer* vectors_right_buf = is_first_compute_step ?
                                     &vectors_buf : &vectors_buf_01[index_01];
/*FIX
    GMFloatMirroredPointer* vectors_right_buf_next 
                                             = &vectors_buf_01[index_01_next];
*/
    /*---Denote metrics buffers---*/

    GMFloatMirroredPointer* metrics_buf = &metrics_buf_01[index_01];
    GMFloatMirroredPointer* metrics_buf_prev = &metrics_buf_01[index_01_prev];

    /*---Prepare for sends/recvs procs for communication---*/

    const int proc_up = (env->proc_num + 1) % env->num_proc;
    const int proc_dn = (env->proc_num - 1 + env->num_proc) % env->num_proc;

    MPI_Request mpi_requests[2];

    /*---Initiate sends/recvs for vecs needed on next step---*/

    if ( is_compute_step_next ) {
      mpi_requests[0] = gm_send_vectors_start(vectors_right, proc_dn, env );
      mpi_requests[1] = gm_recv_vectors_start(vectors_right_next,
                                              proc_up, env );
    }

    /*---The proc that owns the "right-side" vecs for the minproduct---*/

    const int j_proc = (env->proc_num + step_num) % env->num_proc;
    const int j_proc_prev = (env->proc_num + step_num-1) % env->num_proc;

    /*---To remove redundancies from symmetry, skip some blocks---*/

    const GMBool skipping_active =
                ( env->num_proc % 2 == 0 ) &&
                ( 2 * env->proc_num >= env->num_proc );

    const GMBool skipped_last_block_lower_half = skipping_active &&
                is_last_compute_step;

    const GMBool skipped_last_block_lower_half_prev = skipping_active &&
                is_last_compute_step_prev;

    const GMBool do_compute_block = is_compute_step &&
                                 ! skipped_last_block_lower_half;

    const GMBool do_compute_block_prev = is_compute_step_prev &&
                                 ! skipped_last_block_lower_half_prev;

    /*---Main diagonal block only computes strict upper triangular part---*/

    const GMBool do_compute_triang_only = is_first_compute_step;
    const GMBool do_compute_triang_only_prev = is_first_compute_step_prev;

    /*---Commence numerators computation---*/

    if ( is_compute_step && do_compute_block ) {

      if ( env->compute_method == GM_COMPUTE_METHOD_GPU ) {
        /*---Copy vectors into GPU buffers if needed---*/
        for (i = 0; i < vectors->num_vector_local; ++i) {
          int k = 0;
          for (k = 0; k < vectors->num_field; ++k) {
            vectors_left_buf->h[k + vectors->num_field*i] =
              GMVectors_float_get(vectors_left, k, i, env);
          }
        }
        for (i = 0; i < vectors->num_vector_local; ++i) {
          int k = 0;
          for (k = 0; k < vectors->num_field; ++k) {
            vectors_right_buf->h[k + vectors->num_field*i] =
              GMVectors_float_get(vectors_right, k, i, env);
          }
        }
        /*---Send vectors to GPU---*/
        gm_set_float_vectors_start(vectors_left, vectors_left_buf, env);
        gm_set_float_vectors_wait(env);
        gm_set_float_vectors_start(vectors_right, vectors_right_buf, env);
        gm_set_float_vectors_wait(env);
      }

      /*---Compute numerators---*/
      gm_compute_czekanowski_numerators_start(vectors_left, vectors_right,
           metrics, vectors_left_buf, vectors_right_buf, metrics_buf,
           j_proc, do_compute_triang_only, env);

    } /*---is_compute_step ...---*/

    /*---GPU case: wait for prev step get metrics to complete, then combine.
         Note this is hidden under GPU computation---*/

    if ( env->compute_method == GM_COMPUTE_METHOD_GPU ) {
      if ( is_compute_step_prev && do_compute_block_prev ) {
        gm_get_float_metrics_wait(env);
        GMFloat* vector_sums_left = vector_sums_onproc;
        GMFloat* vector_sums_right = is_first_compute_step_prev ?
                                     vector_sums_onproc : vector_sums_offproc;
        gm_compute_czekanowski_combine(metrics, metrics_buf_prev,
                                       vector_sums_left,
                                       vector_sums_right, j_proc_prev,
                                       do_compute_triang_only_prev, env);
      }
    }

    /*---Wait for numerators computation to complete---*/

    if ( is_compute_step && do_compute_block ) {
        gm_compute_czekanowski_numerators_wait(env);
    }

    /*---Commence copy of completed numerators back from GPU---*/

    if ( is_compute_step && do_compute_block ) {
        gm_get_float_metrics_start(metrics, metrics_buf, env);
    }

    /*---Compute sums for denominators---*/

    if ( is_compute_step && do_compute_block ) {
      if (is_first_compute_step) {
        gm_compute_float_vector_sums(vectors_left, vector_sums_onproc, env);
      } else {
        gm_compute_float_vector_sums(vectors_right, vector_sums_offproc, env);
      }
    }

    /*---CPU case: combine numerators, denominators to obtain final result---*/

    if ( env->compute_method != GM_COMPUTE_METHOD_GPU ) {
      if ( is_compute_step && do_compute_block ) {
        GMFloat* vector_sums_left = vector_sums_onproc;
        GMFloat* vector_sums_right =
            is_first_compute_step ? vector_sums_onproc : vector_sums_offproc;
        gm_compute_czekanowski_combine(metrics, metrics_buf,
                                       vector_sums_left,
                                       vector_sums_right, j_proc,
                                       do_compute_triang_only, env);
      }
    }

    /*---Wait for sends/recvs to complete---*/

    if ( is_compute_step_next ) {
      gm_send_vectors_wait(&(mpi_requests[0]), env);
      gm_recv_vectors_wait(&(mpi_requests[1]), env);
    }

  } /*---step_num---*/

  /*----------------------------------------*/
  /*---End loop over steps of circular shift of vectors objects---*/
  /*----------------------------------------*/

  free(vector_sums_onproc);
  free(vector_sums_offproc);

  for ( i = 0 ; i < 2 ; ++i ) {
    GMVectors_destroy(&vectors_01[i], env);
  }

  /*---Magma terminations---*/

  for ( i = 0 ; i < 2 ; ++i ) {
    GMFloat_free_magma_minproduct( &metrics_buf_01[i], env );
    GMFloat_free_magma_minproduct( &vectors_buf_01[i], env );
  }
  GMFloat_free_magma_minproduct( &vectors_buf, env );

  if (env->compute_method == GM_COMPUTE_METHOD_GPU ) {
    magma_code = magma_minproduct_finalize();
    GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
  }
}

/*===========================================================================*/

void gm_compute_metrics_czekanowski_2way_cpu(GMMetrics* metrics,
                                             GMVectors* vectors,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  if (env->all2all) {

    /*---Initializations---*/

    int data_type = gm_data_type_from_metric_type(env->metric_type, env);

    GMFloat* vector_sums_onproc = GMFloat_malloc(metrics->num_vector_local);
    GMFloat* vector_sums_offproc = GMFloat_malloc(metrics->num_vector_local);

    /*---Create size-2 circular buffer of vectors objects for send/recv---*/

    GMVectors vectors_0 = GMVectors_null();
    GMVectors vectors_1 = GMVectors_null();

    GMVectors_create(&vectors_0, data_type, vectors->num_field,
                     vectors->num_vector_local, env);
    GMVectors_create(&vectors_1, data_type, vectors->num_field,
                     vectors->num_vector_local, env);

    GMVectors* vectors_01[2] = {&vectors_0, &vectors_1};

    /*---Result matrix is diagonal block and half the blocks to the right
         (including wraparound to left side of matrix when appropriate).
         For even number of procs, block rows of lower half of matrix
         have one less block to make correct count---*/

    const int num_step = 1 + (env->num_proc / 2);

    /*---Loop over steps of circular shift of vectors objects---*/

    int step_num = 0;
    for (step_num = 0; step_num < num_step; ++step_num) {

      const GMBool is_first_step = step_num == 0;
      const GMBool is_last_step = step_num == num_step - 1;

      /*---The proc that owns the "right-side" vectors for the minproduct---*/

      const int j_proc = (env->proc_num + step_num) % env->num_proc;

      /*---Left- and right-side vecs, also right-side vecs for next step---*/

      GMVectors* vectors_left = vectors;
      GMVectors* vectors_right =
          is_first_step ? vectors : vectors_01[step_num % 2];
      GMVectors* vectors_right_next = vectors_01[(step_num + 1) % 2];

      /*---Initialize sends/recvs---*/

      MPI_Request mpi_requests[2];

      int mpi_code = 0;
      mpi_code = mpi_code*1; /*---Avoid unused variable warning---*/

      const int proc_up = (env->proc_num + 1) % env->num_proc;
      const int proc_dn = (env->proc_num - 1 + env->num_proc) % env->num_proc;

      const int mpi_tag = 0;

      /*---Initiate sends/recvs for vecs needed on next step---*/

      if (! is_last_step) {
        mpi_code = MPI_Isend(
            (void*)vectors_right->data, vectors_right->num_dataval_local,
            GM_MPI_FLOAT, proc_dn, mpi_tag, env->mpi_comm, &(mpi_requests[0]));
        GMAssert(mpi_code == MPI_SUCCESS);

        mpi_code =
            MPI_Irecv((void*)vectors_right_next->data,
                      vectors_right_next->num_dataval_local, GM_MPI_FLOAT,
                      proc_up, mpi_tag, env->mpi_comm, &(mpi_requests[1]));
        GMAssert(mpi_code == MPI_SUCCESS);

      } /*---if step_num---*/

      /*---Compute sums for denominators---*/

      if (is_first_step) {
        gm_compute_float_vector_sums(vectors_left, vector_sums_onproc, env);
      } else {
        gm_compute_float_vector_sums(vectors_right, vector_sums_offproc, env);
      }

      GMFloat* vector_sums_left = vector_sums_onproc;
      GMFloat* vector_sums_right =
          is_first_step ? vector_sums_onproc : vector_sums_offproc;

      /*---Compute numerators---*/

      const GMBool compute_triang_only = is_first_step;

      const GMBool skipped_last_block_lower_half = env->num_proc % 2 == 0 &&
                        ( 2 * env->proc_num >= env->num_proc ) && is_last_step;

      if ( ! skipped_last_block_lower_half ) {
        int j = 0;
        for (j = 0; j < metrics->num_vector_local; ++j) {
          const int i_max = compute_triang_only ? j : metrics->num_vector_local;
          int i = 0;
          for (i = 0; i < i_max; ++i) {
          GMFloat sum = 0;
            int field = 0;
            for (field = 0; field < vectors->num_field; ++field) {
              const GMFloat value1 =
                  GMVectors_float_get(vectors_left, field, i, env);
              const GMFloat value2 =
                  GMVectors_float_get(vectors_right, field, j, env);
              sum += value1 < value2 ? value1 : value2;
            } /*---for k---*/
            GMMetrics_float_set_all2all_2(metrics, i, j, j_proc, sum, env);
          } /*---for i---*/
        }   /*---for j---*/
      } /*---if---*/

      /*---Combine---*/

      if ( ! skipped_last_block_lower_half ) {
        int j = 0;
        for (j = 0; j < metrics->num_vector_local; ++j) {
          const int i_max = compute_triang_only ? j : metrics->num_vector_local;
          int i = 0;
          for (i = 0; i < i_max; ++i) {
            const GMFloat numerator =
                GMMetrics_float_get_all2all_2(metrics, i, j, j_proc, env);
            const GMFloat denominator =
                vector_sums_left[i] + vector_sums_right[j];
            GMMetrics_float_set_all2all_2(metrics, i, j, j_proc,
                                          2 * numerator / denominator, env);
          } /*---for i---*/
        }   /*---for j---*/
      } /*---if---*/

      /*---Wait for sends/recvs to complete---*/

      if (! is_last_step) {
        MPI_Status mpi_status;
        mpi_code = MPI_Waitall(2, &(mpi_requests[0]), &mpi_status);
        GMAssert(mpi_code == MPI_SUCCESS);
      }

    } /*---step_num---*/

    /*---Deallocations---*/

    free(vector_sums_onproc);
    free(vector_sums_offproc);

    GMVectors_destroy(&vectors_0, env);
    GMVectors_destroy(&vectors_1, env);

  } else /*---if ( ! env->all2all )---*/ {

    /*---Compute sums for denominators---*/

    GMFloat* vector_sums = GMFloat_malloc(metrics->num_vector_local);

    gm_compute_float_vector_sums(vectors, vector_sums, env);

    /*---Compute numerators---*/

    int i = 0;
    int j = 0;
    for (i = 0; i < metrics->num_vector_local; ++i) {
      for (j = i + 1; j < metrics->num_vector_local; ++j) {
        GMFloat sum = 0;
        int field = 0;
        for (field = 0; field < vectors->num_field; ++field) {
          const GMFloat value1 = GMVectors_float_get(vectors, field, i, env);
          const GMFloat value2 = GMVectors_float_get(vectors, field, j, env);
          sum += value1 < value2 ? value1 : value2;
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

  } /*---if ( env->all2all )---*/
}

/*===========================================================================*/

void gm_compute_metrics_czekanowski_2way_gpu(GMMetrics* metrics,
                                             GMVectors* vectors,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  if (env->all2all) {
    gm_compute_metrics_czekanowski_2way_all2all_gpu(metrics, vectors, env);
    return;
  }

  /*---Denominator---*/

  GMFloat* vector_sums = GMFloat_malloc(metrics->num_vector_local);

  /* .02 / 1.56 */
  gm_compute_float_vector_sums(vectors, vector_sums, env);

  /*---------------*/
  /*---Numerator---*/
  /*---------------*/

  int i = 0;

  magma_minproduct_int_t magma_code = 0;
  magma_code = magma_code*1; /*---Avoid unused variable warning---*/

  magma_code = magma_minproduct_init();
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);

  const int num_vec = metrics->num_vector_local;
  const int num_field = vectors->num_field;

  /*---Allocate magma CPU memory for vectors and for result */

  GMFloatMirroredPointer vectors_buf =
    GMFloat_malloc_magma_minproduct( num_vec * (size_t)num_field, env);

  GMFloatMirroredPointer numerators_buf =
    GMFloat_malloc_magma_minproduct( num_vec * (size_t)num_vec, env);

  /*---Copy in vectors---*/

  /* .08 / 1.56 */
  for (i = 0; i < num_vec; ++i) {
    int k = 0;
    for (k = 0; k < num_field; ++k) {
      vectors_buf.h[k + num_field*i] = GMVectors_float_get(vectors, k, i, env);
    }
  }

  /*---Send vectors to GPU---*/

  gm_set_float_vectors_start(vectors, &vectors_buf, env);
  gm_set_float_vectors_wait(env);

  /*---Initialize result matrix to zero (apparently magma requires)---*/

#ifdef FP_PRECISION_DOUBLE
  magma_minproductblas_dlaset(Magma_minproductFull, num_vec, num_vec,
                              0.0, 0.0, numerators_buf.d, num_vec);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproductblas_slaset(Magma_minproductFull, num_vec, num_vec,
                              0.0, 0.0, numerators_buf.d, num_vec);
#endif

  /*---Perform pseudo matrix-matrix product---*/

/* .63 / 1.56 */
#ifdef FP_PRECISION_DOUBLE
  magma_minproductblas_dgemm
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproductblas_sgemm
#endif
      (Magma_minproductTrans, Magma_minproductNoTrans, num_vec, num_vec,
       num_field, 1.0, vectors_buf.d, num_field, vectors_buf.d, num_field, 0.0,
       numerators_buf.d, num_vec);

  /*---Copy result from GPU---*/

  gm_get_float_metrics_start(metrics, &numerators_buf, env);
  gm_get_float_metrics_wait(env);

  /*---Combine---*/

  /* .22 / 1.56 */
  for (i = 0; i < metrics->num_vector_local; ++i) {
    int j = 0;
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      const GMFloat numerator = numerators_buf.h[j + num_vec * i];
      const GMFloat denominator = vector_sums[i] + vector_sums[j];
      GMMetrics_float_set_2(metrics, i, j, 2 * numerator / denominator, env);
    } /*---for j---*/
  }   /*---for i---*/

  /*---Free memory---*/

  GMFloat_free_magma_minproduct( &vectors_buf, env);
  GMFloat_free_magma_minproduct( &numerators_buf, env);

  magma_code = magma_minproduct_finalize();
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);

  free(vector_sums);
}

/*---------------------------------------------------------------------------*/
