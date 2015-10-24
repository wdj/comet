/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_utils.c
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing metrics, utilities.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include "cuda.h"

#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"

#include "env.h"
#include "vectors.h"
#include "metrics.h"

/*===========================================================================*/
/*---Compute the sum of elements of each vector on CPU, for denom---*/

void gm_compute_float_vector_sums(GMVectors* vectors,
                                  GMFloat* __restrict__ vector_sums,
                                  GMEnv* env) {
  GMAssert(vectors != NULL);
  GMAssert(vector_sums != NULL);
  GMAssert(env != NULL);

  int i = 0;
  for (i = 0; i < vectors->num_vector_local; ++i) {
    GMFloat sum = 0;
    int field = 0;
    for (field = 0; field < vectors->num_field; ++field) {
      GMFloat value = GMVectors_float_get(vectors, field, i, env);
      sum += value;
    }
    vector_sums[i] = sum;
  }
}

/*===========================================================================*/
/*---Allocate/free host and device memory---*/

GMFloatMirroredPointer GMFloat_malloc_magma_minproduct(size_t n, GMEnv* env) {
  GMAssert(n+1 >= 1);
  GMAssert(env != NULL);

  magma_minproduct_int_t magma_code = 0;
  magma_code = magma_code*1; /*---Avoid unused variable warning---*/

  GMFloatMirroredPointer p = GMFloatMirroredPointer_null();

#ifdef FP_PRECISION_DOUBLE
  magma_code = magma_minproduct_dmalloc_pinned((GMFloat**)&p.h, n);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_code = magma_minproduct_smalloc_pinned((GMFloat**)&p.h, n);
#endif
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);

#ifdef FP_PRECISION_DOUBLE
  magma_code = magma_minproduct_dmalloc((GMFloat**)&p.d, n);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_code = magma_minproduct_smalloc((GMFloat**)&p.d, n);
#endif
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);

  return p;
}

/*---------------------------------------------------------------------------*/

void GMFloat_free_magma_minproduct(GMFloatMirroredPointer* p, GMEnv* env) {
  GMAssert(p != NULL);
  GMAssert(env != NULL);

  magma_minproduct_int_t magma_code = 0;
  magma_code = magma_code*1; /*---Avoid unused variable warning---*/

  magma_minproduct_free_pinned(p->h);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
  magma_minproduct_free(p->d);
  GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
}

/*===========================================================================*/
/*---Start/end MPI send/receive of vectors data---*/

MPI_Request gm_send_vectors_start(GMVectors* vectors,
                                  int proc_num,
                                  GMEnv* env) {
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);
  GMAssert(proc_num >= 0 && proc_num < env->num_proc);

  const int mpi_tag = 0;
  MPI_Request mpi_request;
  int mpi_code = 0;
  mpi_code = mpi_code*1; /*---Avoid unused variable warning---*/

  mpi_code = MPI_Isend(
      (void*)vectors->data, vectors->num_dataval_local,
      GM_MPI_FLOAT, proc_num, mpi_tag, env->mpi_comm, &mpi_request);
  GMAssert(mpi_code == MPI_SUCCESS);

  return mpi_request;
}

/*---------------------------------------------------------------------------*/

MPI_Request gm_recv_vectors_start(GMVectors* vectors,
                                  int proc_num,
                                  GMEnv* env) {
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);
  GMAssert(proc_num >= 0 && proc_num < env->num_proc);

  const int mpi_tag = 0;
  MPI_Request mpi_request;
  int mpi_code = 0;
  mpi_code = mpi_code*1; /*---Avoid unused variable warning---*/

  mpi_code = MPI_Irecv(
      (void*)vectors->data, vectors->num_dataval_local,
      GM_MPI_FLOAT, proc_num, mpi_tag, env->mpi_comm, &mpi_request);
  GMAssert(mpi_code == MPI_SUCCESS);

  return mpi_request;
}

/*---------------------------------------------------------------------------*/

void gm_send_vectors_wait(MPI_Request* mpi_request, GMEnv* env) {
  GMAssert(mpi_request != NULL);
  GMAssert(env != NULL);

  MPI_Status mpi_status;
  int mpi_code = 0;
  mpi_code = mpi_code*1; /*---Avoid unused variable warning---*/

  mpi_code = MPI_Wait(mpi_request, &mpi_status);
  GMAssert(mpi_code == MPI_SUCCESS);
}

/*---------------------------------------------------------------------------*/

void gm_recv_vectors_wait(MPI_Request* mpi_request, GMEnv* env) {
  GMAssert(mpi_request != NULL);
  GMAssert(env != NULL);

  MPI_Status mpi_status;
  int mpi_code = 0;
  mpi_code = mpi_code*1; /*---Avoid unused variable warning---*/

  mpi_code = MPI_Wait(mpi_request, &mpi_status);
  GMAssert(mpi_code == MPI_SUCCESS);
}

/*===========================================================================*/
/*---Start/end transfer of vectors data to GPU---*/

void gm_set_float_vectors_start(GMVectors* vectors,
                                GMFloatMirroredPointer* vectors_buf,
                                GMEnv* env) {
  GMAssert(vectors != NULL);
  GMAssert(vectors_buf != NULL);
  GMAssert(env != NULL);

  const int numvec = vectors->num_vector_local;
  const int numfield = vectors->num_field;

  /*---Send vectors to GPU---*/

  if (env->compute_method == GM_COMPUTE_METHOD_GPU) {
    GMEnv_initialize_streams(env);
#ifdef FP_PRECISION_DOUBLE
    magma_minproduct_dsetmatrix_async(numfield, numvec, vectors_buf->h,
                     numfield, vectors_buf->d, numfield, env->stream_vectors);
#endif
#ifdef FP_PRECISION_SINGLE
    magma_minproduct_ssetmatrix_async(numfield, numvec, vectors_buf->h,
                     numfield, vectors_buf->d, numfield, env->stream_vectors);
#endif
  }
}

/*---------------------------------------------------------------------------*/

void gm_set_float_vectors_wait(GMEnv* env) {
  GMAssert(env != NULL);

  if (env->compute_method == GM_COMPUTE_METHOD_GPU) {
    GMEnv_initialize_streams(env);
    cudaStreamSynchronize( env->stream_vectors );
    GMAssert(GMEnv_cuda_last_call_succeeded(env));
  }
}

/*===========================================================================*/
/*---Start/end transfer of metrics data from GPU---*/

void gm_get_float_metrics_start(GMMetrics* metrics,
                                GMFloatMirroredPointer* metrics_buf,
                                GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(metrics_buf != NULL);
  GMAssert(env != NULL);

  const int numvec = metrics->num_vector_local;

  /*---Send vectors to GPU---*/

  if (env->compute_method == GM_COMPUTE_METHOD_GPU) {
    GMEnv_initialize_streams(env);
#ifdef FP_PRECISION_DOUBLE
    magma_minproduct_dgetmatrix_async(numvec, numvec, metrics_buf->d,
                          numvec, metrics_buf->h, numvec, env->stream_metrics);
#endif
#ifdef FP_PRECISION_SINGLE
    magma_minproduct_sgetmatrix_async(numvec, numvec, metrics_buf->d,
                          numvec, metrics_buf->h, numvec, env->stream_metrics);
#endif
  }
}

/*---------------------------------------------------------------------------*/

void gm_get_float_metrics_wait(GMEnv* env) {
  GMAssert(env != NULL);

  if (env->compute_method == GM_COMPUTE_METHOD_GPU) {
    GMEnv_initialize_streams(env);
    cudaStreamSynchronize( env->stream_metrics );
    GMAssert(GMEnv_cuda_last_call_succeeded(env));
  }
}

/*===========================================================================*/
/*---CPU-GPU transfer nbuffer manipulation---*/

void gm_float_vectors_to_buf(GMVectors* vectors,
                             GMFloatMirroredPointer* vectors_buf,
                             GMEnv* env) {
  GMAssert(vectors != NULL);
  GMAssert(vectors_buf != NULL);
  GMAssert(env != NULL);

  if ( env->compute_method == GM_COMPUTE_METHOD_GPU ) {
    int i = 0;
    /*---Copy vectors into GPU buffers if needed---*/
    for (i = 0; i < vectors->num_vector_local; ++i) {
      int k = 0;
      for (k = 0; k < vectors->num_field; ++k) {
        vectors_buf->h[k + vectors->num_field*i] =
          GMVectors_float_get(vectors, k, i, env);
      }
    }
  }
}

/*===========================================================================*/
/*---Start/end calculation of numerators---*/

void gm_compute_czekanowski_numerators_start(
                                      GMVectors* vectors_left,
                                      GMVectors* vectors_right,
                                      GMMetrics* numerators,
                                      GMFloatMirroredPointer* vectors_left_buf,
                                      GMFloatMirroredPointer* vectors_right_buf,
                                      GMFloatMirroredPointer* numerators_buf,
                                      int j_proc,
                                      GMBool compute_triang_only,
                                      GMEnv* env) {
  GMAssert(vectors_left != NULL);
  GMAssert(vectors_right != NULL);
  GMAssert(numerators != NULL);
  GMAssert(env != NULL);
  GMAssert(j_proc >= 0 && j_proc < env->num_proc);

  /*----------------------------------------*/
  if ( env->compute_method != GM_COMPUTE_METHOD_GPU ) {
  /*----------------------------------------*/

    /*---Perform pseudo matrix-matrix product---*/

    int j = 0;
    for (j = 0; j < numerators->num_vector_local; ++j) {
      const int i_max = compute_triang_only ? j : numerators->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        GMFloat numerator = 0;
        int field = 0;
        for (field = 0; field < vectors_left->num_field; ++field) {
          const GMFloat value1 =
              GMVectors_float_get(vectors_left, field, i, env);
          const GMFloat value2 =
              GMVectors_float_get(vectors_right, field, j, env);
          numerator += value1 < value2 ? value1 : value2;
        } /*---for k---*/
        GMMetrics_float_set_all2all_2(numerators, i, j, j_proc, numerator, env);
      } /*---for i---*/
    }   /*---for j---*/

  /*----------------------------------------*/
  } else /* if (env->compute_method == GM_COMPUTE_METHOD_GPU) */ {
  /*----------------------------------------*/

    /*---Initialize result matrix to zero (apparently magma requires)---*/

#ifdef FP_PRECISION_DOUBLE
    magma_minproductblas_dlaset
#endif
#ifdef FP_PRECISION_SINGLE
    magma_minproductblas_slaset
#endif
        (Magma_minproductFull,
        numerators->num_vector_local, numerators->num_vector_local,
        0.0, 0.0, numerators_buf->d, numerators->num_vector_local);

    /*---Perform pseudo matrix-matrix product---*/

/* .63 / 1.56 */
#ifdef FP_PRECISION_DOUBLE
    magma_minproductblas_dgemm
#endif
#ifdef FP_PRECISION_SINGLE
    magma_minproductblas_sgemm
#endif
        (Magma_minproductTrans, Magma_minproductNoTrans,
         vectors_left->num_vector_local, vectors_left->num_vector_local,
         vectors_left->num_field, 1.0,
         vectors_left_buf->d, vectors_left->num_field,
         vectors_right_buf->d, vectors_left->num_field,
         0.0, numerators_buf->d, vectors_left->num_vector_local);

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*---------------------------------------------------------------------------*/

void gm_compute_czekanowski_numerators_wait(GMEnv* env) {
  GMAssert(env != NULL);

  if (env->compute_method == GM_COMPUTE_METHOD_GPU) {
    GMEnv_initialize_streams(env);
    cudaStreamSynchronize( env->stream_compute );
    GMAssert(GMEnv_cuda_last_call_succeeded(env));
  }
}

/*===========================================================================*/
/*---Combine numerators and denominators on CPU to get final result---*/

void gm_compute_czekanowski_combine(GMMetrics* metrics,
                                    GMFloatMirroredPointer* metrics_buf,
                                    GMFloat* __restrict__ vector_sums_left,
                                    GMFloat* __restrict__ vector_sums_right,
                                    int j_proc,
                                    GMBool compute_triang_only,
                                    GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vector_sums_left != NULL);
  GMAssert(vector_sums_right != NULL);
  GMAssert(env != NULL);
  GMAssert(j_proc >= 0 && j_proc < env->num_proc);

  /*---For CPU case, copy numerator out of metrics struct which is temporarily
       holding numerators.
       For GPU case, directly access the metrics_buf holding the numerators.
  ---*/


  /*----------------------------------------*/
  if ( env->compute_method != GM_COMPUTE_METHOD_GPU ) {
  /*----------------------------------------*/

    int j = 0;
    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        const GMFloat numerator = GMMetrics_float_get_all2all_2(
                                                   metrics, i, j, j_proc, env);
        const GMFloat denominator = vector_sums_left[i] + vector_sums_right[j];
        GMMetrics_float_set_all2all_2(metrics, i, j, j_proc,
                                      2 * numerator / denominator, env);
      } /*---for i---*/
    }   /*---for j---*/
/*---TODO: here and elsewhere check for unlikely case denom is/nearly zero---*/

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

    int j = 0;
    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        const GMFloat numerator = metrics_buf->h[i +
                                              metrics->num_vector_local * j];
        const GMFloat denominator = vector_sums_left[i] + vector_sums_right[j];
        GMMetrics_float_set_all2all_2(metrics, i, j, j_proc,
                                      2 * numerator / denominator, env);
      } /*---for i---*/
    }   /*---for j---*/

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
