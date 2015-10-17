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

    const GMBool is_comp_step = step_num >= 0 && step_num < num_step;
    const GMBool is_first_comp_step = step_num == 0;
    const GMBool is_last_comp_step = step_num == num_step - 1;

if(is_comp_step){}/*FIX*/

    /*---Left- and right-side vecs, also right-side vecs for next step---*/

    GMVectors* vectors_left = vectors;
    GMVectors* vectors_right =
        is_first_comp_step ? vectors : vectors_01[step_num % 2];
    GMVectors* vectors_right_next = vectors_01[(step_num + 1) % 2];

    /*---Prepare for sends/recvs---*/

    MPI_Request mpi_requests[2];

    int mpi_code = 0;
    mpi_code = mpi_code ? 0 : 0; /*---Avoid unused variable warning---*/

    const int proc_up = (env->proc_num + 1) % env->num_proc;
    const int proc_dn = (env->proc_num - 1 + env->num_proc) % env->num_proc;

    /*---Initiate sends/recvs for vecs needed on next step---*/

    if (! is_last_comp_step) {
      mpi_requests[0] = gm_send_vectors_start(vectors_right, proc_dn, env );
      mpi_requests[1] = gm_recv_vectors_start(vectors_right_next,
                                              proc_up, env );

    }

    /*---The proc that owns the "right-side" vectors for the minproduct---*/

    const int j_proc = (env->proc_num + step_num) % env->num_proc;

    /*---Compute numerators---*/

    const GMBool compute_triang_only = is_first_comp_step;

    const GMBool skipped_last_block_lower_half = env->num_proc % 2 == 0 &&
                  ( 2 * env->proc_num >= env->num_proc ) && is_last_comp_step;

    if ( ! skipped_last_block_lower_half ) {
      gm_compute_czekanowski_numerators_start(vectors_left, vectors_right,
                                              metrics, j_proc,
                                              compute_triang_only, env);
      gm_compute_czekanowski_numerators_wait(env);
    }

    /*---Compute sums for denominators---*/

    if (is_first_comp_step) {
      gm_compute_float_vector_sums(vectors_left, vector_sums_onproc, env);
    } else {
      gm_compute_float_vector_sums(vectors_right, vector_sums_offproc, env);
    }

    GMFloat* vector_sums_left = vector_sums_onproc;
    GMFloat* vector_sums_right =
        is_first_comp_step ? vector_sums_onproc : vector_sums_offproc;

    /*---Combine numerators, denominators to obtain final result---*/

    if ( ! skipped_last_block_lower_half ) {
      gm_compute_czekanowski_combine(metrics, vector_sums_left,
                                     vector_sums_right, j_proc,
                                     compute_triang_only, env);
    }

    /*---Wait for sends/recvs to complete---*/

    if (! is_last_comp_step) {
      gm_send_vectors_wait(&(mpi_requests[0]), env);
      gm_recv_vectors_wait(&(mpi_requests[1]), env);
    }

  } /*---step_num---*/

  /*---Deallocations---*/

  free(vector_sums_onproc);
  free(vector_sums_offproc);

  GMVectors_destroy(&vectors_0, env);
  GMVectors_destroy(&vectors_1, env);
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
      mpi_code = mpi_code ? 0 : 0; /*---Avoid unused variable warning---*/

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
  magma_code = magma_code ? 0 : 0; /*---Avoid unused variable warning---*/

  magma_code = magma_minproduct_init();
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);

  const int numvec = metrics->num_vector_local;
  const int numfield = vectors->num_field;

  GMFloat* h_vectors = NULL;
  GMFloat* h_numer = NULL;

/*---Allocate magma CPU memory for vectors and for result */

#ifdef FP_PRECISION_DOUBLE
  magma_code = magma_minproduct_dmalloc_pinned(&h_vectors, numvec * numfield);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
  magma_code = magma_minproduct_dmalloc_pinned(&h_numer, numvec * numvec);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_code = magma_minproduct_smalloc_pinned(&h_vectors, numvec * numfield);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
  magma_code = magma_minproduct_smalloc_pinned(&h_numer, numvec * numvec);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
#endif

  /*---Allocate GPU mirrors for CPU arrays---*/

  GMFloat* d_vectors = NULL;
  GMFloat* d_numer = NULL;

#ifdef FP_PRECISION_DOUBLE
  magma_code = magma_minproduct_dmalloc(&d_vectors, numvec * numfield);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
  magma_code = magma_minproduct_dmalloc(&d_numer, numvec * numvec);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_code = magma_minproduct_smalloc(&d_vectors, numvec * numfield);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
  magma_code = magma_minproduct_smalloc(&d_numer, numvec * numvec);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
#endif

/*---Initialize result matrix to zero (apparently magma requires)---*/

#ifdef FP_PRECISION_DOUBLE
  magma_minproductblas_dlaset(Magma_minproductFull, numvec, numvec, 0.0, 0.0,
                              d_numer, numvec);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproductblas_slaset(Magma_minproductFull, numvec, numvec, 0.0, 0.0,
                              d_numer, numvec);
#endif

  /*---Copy in vectors---*/

  /* .08 / 1.56 */
  for (i = 0; i < numvec; ++i) {
    int k = 0;
    for (k = 0; k < numfield; ++k) {
      h_vectors[k + numfield * i] = GMVectors_float_get(vectors, k, i, env);
    }
  }

/*---Send vectors to GPU---*/

#ifdef FP_PRECISION_DOUBLE
  magma_minproduct_dsetmatrix(numfield, numvec, h_vectors, numfield, d_vectors,
                              numfield);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproduct_ssetmatrix(numfield, numvec, h_vectors, numfield, d_vectors,
                              numfield);
#endif

/*---Perform pseudo matrix-matrix product---*/

/* .63 / 1.56 */
#ifdef FP_PRECISION_DOUBLE
  magma_minproductblas_dgemm
#endif
#ifdef FP_PRECISION_SINGLE
      magma_minproductblas_sgemm
#endif
      (Magma_minproductTrans, Magma_minproductNoTrans, numvec, numvec, numfield,
       1.0, d_vectors, numfield, d_vectors, numfield, 0.0, d_numer, numvec);

/*---Copy result from GPU---*/

#ifdef FP_PRECISION_DOUBLE
  magma_minproduct_dgetmatrix(numvec, numvec, d_numer, numvec, h_numer, numvec);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproduct_sgetmatrix(numvec, numvec, d_numer, numvec, h_numer, numvec);
#endif

  /*---Combine---*/

  /* .22 / 1.56 */
  for (i = 0; i < metrics->num_vector_local; ++i) {
    int j = 0;
    for (j = i + 1; j < metrics->num_vector_local; ++j) {
      const GMFloat numerator = h_numer[j + numvec * i];
      const GMFloat denominator = vector_sums[i] + vector_sums[j];
      GMMetrics_float_set_2(metrics, i, j, 2 * numerator / denominator, env);
    } /*---for j---*/
  }   /*---for i---*/

  /*---Free memory---*/

  magma_minproduct_free(d_vectors);
  magma_minproduct_free(d_numer);

  magma_minproduct_free_pinned(h_vectors);
  magma_minproduct_free_pinned(h_numer);

  magma_minproduct_finalize();

  free(vector_sums);
}

/*---------------------------------------------------------------------------*/
