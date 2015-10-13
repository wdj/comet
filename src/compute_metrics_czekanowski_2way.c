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

void gm_compute_metrics_czekanowski_2way_cpu(GMMetrics* metrics,
                                          GMVectors* vectors,
                                          GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  if ( env->all2all ) {

    /*---Initializations---*/

    GMVectors vectors_0 = GMVectors_null();
    GMVectors vectors_1 = GMVectors_null();

    GMVectors_create(&vectors_0,
                     gm_data_type_from_metric_type(env->metric_type, env),
                     vectors->num_field, vectors->num_vector_local, env);
    GMVectors_create(&vectors_1,
                     gm_data_type_from_metric_type(env->metric_type, env),
                     vectors->num_field, vectors->num_vector_local, env);

    GMVectors* vectors_01[2] = { &vectors_0, &vectors_1 };

    const int num_step = 1 + ( env->num_proc / 2 );

    GMFloat* vector_sums_onproc = malloc(metrics->num_vector_local *
                                         sizeof(GMFloat));
    GMFloat* vector_sums_offproc = malloc(metrics->num_vector_local *
                                          sizeof(GMFloat));

    /*---Loop over steps of circular shift---*/

    int step_num = 0;
    for ( step_num = 0; step_num < num_step; ++step_num ) {

      const int j_proc = ( env->proc_num + step_num ) % env->num_proc;

      GMVectors* vectors_left = vectors;
      GMVectors* vectors_right = step_num == 0 ? vectors
                                               : vectors_01[ step_num % 2 ];
      GMVectors* vectors_right_next = vectors_01[ ( step_num + 1 ) % 2 ];

      /*---Initialize sends/recvs---*/

      MPI_Request send_req;
      MPI_Request recv_req;

      int mpi_code = 0;
      if ( mpi_code ) {
      } /*---Avoid unused variable warning---*/

      const int proc_up = ( env->proc_num + 1 ) % env->num_proc;
      const int proc_dn = ( env->proc_num - 1 + env->num_proc ) % env->num_proc;

      const int mpi_tag = 0;

      /*---Post sends/receives---*/

      if ( step_num != num_step - 1 ) {

        mpi_code = MPI_Isend( (void*)vectors_right->data,
                               vectors_right->num_vector_local*
                               vectors_right->num_field,
                               GM_MPI_FLOAT, proc_dn, mpi_tag,
                               env->mpi_comm, &send_req );
        GMAssert(mpi_code == MPI_SUCCESS);

        mpi_code = MPI_Irecv( (void*)vectors_right_next->data,
                               vectors_right_next->num_vector_local*
                               vectors_right_next->num_field,
                               GM_MPI_FLOAT, proc_up, mpi_tag,
                               env->mpi_comm, &recv_req );
        GMAssert(mpi_code == MPI_SUCCESS);

      } /*---if step_num---*/

      /*---Compute sums for denominators---*/

      if ( step_num == 0 ) {
        gm_compute_float_vector_sums(vectors_left, vector_sums_onproc, env);
      } else {
        gm_compute_float_vector_sums(vectors_right, vector_sums_offproc, env);
      }

      GMFloat* vector_sums_left = vector_sums_onproc;
      GMFloat* vector_sums_right = step_num == 0
                                 ? vector_sums_onproc : vector_sums_offproc;

      /*---Compute numerators---*/

      int j = 0;
      for (j = 0; j < metrics->num_vector_local; ++j) {
        const int i_max = step_num == 0
                                      ? j
                                      : env->num_proc % 2 == 0 && 
                                        2 * env->proc_num >= env->num_proc &&
                                        step_num == num_step - 1
                                      ? 0
                                      : metrics->num_vector_local;
        int i = 0;
        for (i = 0; i < i_max; ++i) {
          GMFloat sum = 0;
          int field = 0;
          for (field = 0; field < vectors->num_field; ++field) {
            const GMFloat value1 = GMVectors_float_get(vectors_left,
                                                       field, i, env);
            const GMFloat value2 = GMVectors_float_get(vectors_right,
                                                       field, j, env);
            sum += value1 < value2 ? value1 : value2;
          } /*---for k---*/
          GMMetrics_float_set_all2all_2(metrics, i, j, j_proc, sum, env);
        }   /*---for i---*/
      } /*---for j---*/

      /*---Combine---*/

      for (j = 0; j < metrics->num_vector_local; ++j) {
        const int i_max = step_num == 0
                                      ? j
                                      : env->num_proc % 2 == 0 && 
                                        2 * env->proc_num >= env->num_proc &&
                                        step_num == num_step - 1
                                      ? 0
                                      : metrics->num_vector_local;
        int i = 0;
        for (i = 0; i < i_max; ++i) {
          const GMFloat numerator = GMMetrics_float_get_all2all_2(metrics,
                                        i, j, j_proc, env);
          const GMFloat denominator = vector_sums_left[i] +
                                      vector_sums_right[j];
          GMMetrics_float_set_all2all_2(metrics, i, j, j_proc,
                                        2 * numerator / denominator, env);
        }   /*---for i---*/
      } /*---for j---*/

      /*---Wait for sends/recvs to complete---*/

      if ( step_num != num_step - 1 ) {
        MPI_Status mpi_status;
        mpi_code = MPI_Waitall( 1, &send_req, &mpi_status );
        GMAssert(mpi_code == MPI_SUCCESS);
        mpi_code = MPI_Waitall( 1, &recv_req, &mpi_status );
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

    GMFloat* vector_sums = malloc(metrics->num_vector_local * sizeof(GMFloat));

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

  GMInsist(env, (!env->all2all) ? "Unimplemented." : 0);

  /*---Denominator---*/

  GMFloat* vector_sums = malloc(metrics->num_vector_local * sizeof(GMFloat));

  /* .02 / 1.56 */
  gm_compute_float_vector_sums(vectors, vector_sums, env);

  /*---------------*/
  /*---Numerator---*/
  /*---------------*/

  int i = 0;

  magma_minproduct_int_t magma_code = 0;
  if ( magma_code ) {} /*---Avoid unused variable warning---*/

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
  magma_minproductblas_dlaset(Magma_minproductFull, numvec, numvec, 0.0, 0.0, d_numer, numvec);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproductblas_slaset(Magma_minproductFull, numvec, numvec, 0.0, 0.0, d_numer, numvec);
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
  magma_minproduct_dsetmatrix(numfield, numvec, h_vectors,
                              numfield, d_vectors, numfield);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproduct_ssetmatrix(numfield, numvec, h_vectors,
                              numfield, d_vectors, numfield);
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
  magma_minproduct_dgetmatrix(numvec, numvec, d_numer, numvec,
                              h_numer, numvec);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproduct_sgetmatrix(numvec, numvec, d_numer, numvec,
                              h_numer, numvec);
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

/*===========================================================================*/

void gm_compute_metrics_czekanowski_3way_cpu(GMMetrics* metrics,
                                          GMVectors* vectors,
                                          GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, (!env->all2all) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

void gm_compute_metrics_czekanowski_3way_gpu(GMMetrics* metrics,
                                          GMVectors* vectors,
                                          GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, (!env->all2all) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
