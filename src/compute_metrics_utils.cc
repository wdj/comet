/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_utils.c
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing metrics, utilities.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include "cuda.h"

#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"

#include "env.h"
#include "vector_sums.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_metrics_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Compute the sum of elements of each vector on CPU, for denom---*/

void gm_compute_float_vector_sums(GMVectors* vectors,
                                  GMFloat* __restrict__ vector_sums,
                                  GMFloat* __restrict__ vector_sums_tmp,
                                  GMEnv* env) {
  GMAssert(vectors != NULL);
  GMAssert(vector_sums != NULL);
  GMAssert(env != NULL);

  const int num_proc = Env_num_proc_field(env);
  GMFloat* __restrict__ vector_sums_local = num_proc==1
                                          ? vector_sums : vector_sums_tmp;

  int i = 0;
  for (i = 0; i < vectors->num_vector_local; ++i) {
    GMFloat sum = 0;
    int field_local = 0;
    for (field_local = 0; field_local < vectors->num_field_local;
         ++field_local) {
      GMFloat value = GMVectors_float_get(vectors, field_local, i, env);
      sum += value;
    }
    vector_sums_local[i] = sum;
  }

  if (num_proc > 1) {
    int mpi_code = 0;
    mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
    mpi_code = MPI_Allreduce(vector_sums_local, vector_sums,
                 vectors->num_vector_local, GM_MPI_FLOAT, MPI_SUM,
                 Env_mpi_comm_field(env));
    GMAssert(mpi_code == MPI_SUCCESS);
  }
}

/*---------------------------------------------------------------------------*/

void gm_compute_vector_sums(GMVectors* vectors,
                            GMVectorSums* vector_sums,
                            GMVectorSums* vector_sums_tmp,
                            GMEnv* env) {
  GMAssert(vectors != NULL);
  GMAssert(vector_sums != NULL);
  GMAssert(env != NULL);

  switch (Env_metric_type(env)) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      /*----------------------------------------*/

      gm_compute_float_vector_sums(vectors, (GMFloat*)vector_sums->data,
                                        (GMFloat*)vector_sums_tmp->data,  env);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CCC: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    default:
      /*----------------------------------------*/
      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/
}

/*===========================================================================*/
/*---Magma setup, teardown---*/

void gm_magma_initialize(GMEnv* env) {
  GMAssert(env != NULL);

  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  switch (Env_metric_type(env)) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      /*----------------------------------------*/

      magma_minproduct_int_t magma_code = 0;
      magma_code = magma_code * 1; /*---Avoid unused variable warning---*/

      magma_code = magma_minproduct_init();
      GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
      /*---need this -- see
       * http://on-demand.gputechconf.com/gtc/2014/presentations/S4158-cuda-streams-best-practices-common-pitfalls.pdf
       * page 14 ---*/
      magma_code = magma_minproductblasSetKernelStream(Env_stream_compute(env));
      GMAssert(magma_code == MAGMA_minproduct_SUCCESS);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CCC: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    default:
      /*----------------------------------------*/
      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/
}

/*---------------------------------------------------------------------------*/

void gm_magma_finalize(GMEnv* env) {
  GMAssert(env != NULL);

  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  switch (Env_metric_type(env)) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      /*----------------------------------------*/

      magma_minproduct_int_t magma_code = 0;
      magma_code = magma_code * 1; /*---Avoid unused variable warning---*/

      magma_code = magma_minproduct_finalize();
      GMAssert(magma_code == MAGMA_minproduct_SUCCESS);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CCC: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    default:
      /*----------------------------------------*/
      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/
}

/*===========================================================================*/
/*---Allocate/free host and device memory---*/

GMMirroredPointer gm_malloc_magma(size_t n, GMEnv* env) {
  GMAssert(n + 1 >= 1);
  GMAssert(env != NULL);

  GMMirroredPointer p = GMMirroredPointer_null();

  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return p;
  }

  switch (Env_metric_type(env)) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      /*----------------------------------------*/

      magma_minproduct_int_t magma_code = 0;
      magma_code = magma_code * 1; /*---Avoid unused variable warning---*/

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

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CCC: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    default:
      /*----------------------------------------*/
      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/

  return p;
}

/*---------------------------------------------------------------------------*/

void gm_free_magma(GMMirroredPointer* p, GMEnv* env) {
  GMAssert(p != NULL);
  GMAssert(env != NULL);

  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  switch (Env_metric_type(env)) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      /*----------------------------------------*/

      magma_minproduct_int_t magma_code = 0;
      magma_code = magma_code * 1; /*---Avoid unused variable warning---*/

      magma_minproduct_free_pinned(p->h);
      GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
      magma_minproduct_free(p->d);
      GMAssert(magma_code == MAGMA_minproduct_SUCCESS);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CCC: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    default:
      /*----------------------------------------*/
      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/
}

/*---------------------------------------------------------------------------*/

void gm_magma_set_matrix_zero_start(GMMirroredPointer* matrix_buf,
                                    int mat_dim1,
                                    int mat_dim2,
                                    GMEnv* env) {
#ifdef FP_PRECISION_DOUBLE
  magma_minproductblas_dlaset
#endif
#ifdef FP_PRECISION_SINGLE
      magma_minproductblas_slaset
#endif
      (Magma_minproductFull, mat_dim1, mat_dim2, (GMFloat)0, (GMFloat)0,
       (GMFloat*)matrix_buf->d, mat_dim1);
}

/*---------------------------------------------------------------------------*/

void magma_gemm_start(magma_minproduct_int_t m,
                      magma_minproduct_int_t n,
                      magma_minproduct_int_t k,
                      void* dA,
                      magma_minproduct_int_t ldda,
                      void* dB,
                      magma_minproduct_int_t lddb,
                      void* dC,
                      magma_minproduct_int_t lddc) {
  const GMFloat alpha = 1;
  const GMFloat beta = 0;

#ifdef FP_PRECISION_DOUBLE
  magma_minproductblas_dgemm
#endif
#ifdef FP_PRECISION_SINGLE
      magma_minproductblas_sgemm
#endif
      (Magma_minproductTrans, Magma_minproductNoTrans, m, n, k, alpha,
       (GMFloat*)dA, ldda, (GMFloat*)dB, lddb, beta, (GMFloat*)dC, lddc);
}

/*---------------------------------------------------------------------------*/
/*---Wait for any computation on the GPU top complete---*/

void gm_compute_wait(GMEnv* env) {
  GMAssert(env != NULL);

  if (Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    cudaStreamSynchronize(Env_stream_compute(env));
    GMAssert(GMEnv_cuda_last_call_succeeded(env));
  }
}

/*===========================================================================*/
/*---Start/end MPI send/receive of vectors data---*/

MPI_Request gm_send_vectors_start(GMVectors* vectors,
                                  int proc_num,
                                  GMEnv* env) {
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);
  GMAssert(proc_num >= 0 && proc_num < Env_num_proc_vector(env));

  const int mpi_tag = 0;
  MPI_Request mpi_request;
  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/

  mpi_code =
      MPI_Isend((void*)vectors->data, vectors->num_packedval_local,
                GM_MPI_FLOAT, proc_num, mpi_tag, Env_mpi_comm_vector(env),
                &mpi_request);
  GMAssert(mpi_code == MPI_SUCCESS);

  return mpi_request;
}

/*---------------------------------------------------------------------------*/

MPI_Request gm_recv_vectors_start(GMVectors* vectors,
                                  int proc_num,
                                  GMEnv* env) {
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);
  GMAssert(proc_num >= 0 && proc_num < Env_num_proc_vector(env));

  const int mpi_tag = 0;
  MPI_Request mpi_request;
  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/

  mpi_code =
      MPI_Irecv((void*)vectors->data, vectors->num_packedval_local,
                GM_MPI_FLOAT, proc_num, mpi_tag, Env_mpi_comm_vector(env),
                &mpi_request);
  GMAssert(mpi_code == MPI_SUCCESS);

  return mpi_request;
}

/*---------------------------------------------------------------------------*/

void gm_send_vectors_wait(MPI_Request* mpi_request, GMEnv* env) {
  GMAssert(mpi_request != NULL);
  GMAssert(env != NULL);

  MPI_Status mpi_status;
  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/

  mpi_code = MPI_Wait(mpi_request, &mpi_status);
  GMAssert(mpi_code == MPI_SUCCESS);
}

/*---------------------------------------------------------------------------*/

void gm_recv_vectors_wait(MPI_Request* mpi_request, GMEnv* env) {
  GMAssert(mpi_request != NULL);
  GMAssert(env != NULL);

  MPI_Status mpi_status;
  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/

  mpi_code = MPI_Wait(mpi_request, &mpi_status);
  GMAssert(mpi_code == MPI_SUCCESS);
}

/*===========================================================================*/
/*---Start/end transfer of generic matrix to GPU---*/

void gm_set_matrix_start(GMMirroredPointer* matrix_buf,
                         int mat_dim1,
                         int mat_dim2,
                         GMEnv* env) {
  GMAssert(matrix_buf != NULL);
  GMAssert(env != NULL);

  if (!Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    return;
  }

  /*---Send vectors to GPU---*/

#ifdef FP_PRECISION_DOUBLE
  magma_minproduct_dsetmatrix_async(mat_dim1, mat_dim2, (GMFloat*)matrix_buf->h,
                                    mat_dim1, (GMFloat*)matrix_buf->d, mat_dim1,
                                    Env_stream_togpu(env));
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproduct_ssetmatrix_async(mat_dim1, mat_dim2, (GMFloat*)matrix_buf->h,
                                    mat_dim1, (GMFloat*)matrix_buf->d, mat_dim1,
                                    Env_stream_togpu(env));
#endif
}

/*---------------------------------------------------------------------------*/

void gm_set_matrix_wait(GMEnv* env) {
  GMAssert(env != NULL);

  if (!Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    return;
  }

  cudaStreamSynchronize(Env_stream_togpu(env));
  GMAssert(GMEnv_cuda_last_call_succeeded(env));
}

/*===========================================================================*/
/*---Start/end transfer of generic matrix from GPU---*/

void gm_get_matrix_start(GMMirroredPointer* matrix_buf,
                         int mat_dim1,
                         int mat_dim2,
                         GMEnv* env) {
  GMAssert(matrix_buf != NULL);
  GMAssert(env != NULL);

  if (!Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    return;
  }

  /*---Get vectors from GPU---*/

#ifdef FP_PRECISION_DOUBLE
  magma_minproduct_dgetmatrix_async(mat_dim1, mat_dim2, (GMFloat*)matrix_buf->d,
                                    mat_dim1, (GMFloat*)matrix_buf->h, mat_dim1,
                                    Env_stream_fromgpu(env));
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproduct_sgetmatrix_async(mat_dim1, mat_dim2, (GMFloat*)matrix_buf->d,
                                    mat_dim1, (GMFloat*)matrix_buf->h, mat_dim1,
                                    Env_stream_fromgpu(env));
#endif
}

/*---------------------------------------------------------------------------*/

void gm_get_matrix_wait(GMEnv* env) {
  GMAssert(env != NULL);

  if (!Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    return;
  }

  cudaStreamSynchronize(Env_stream_fromgpu(env));
  GMAssert(GMEnv_cuda_last_call_succeeded(env));
}

/*===========================================================================*/
/*---Start/end transfer of vectors data to GPU---*/

void gm_set_vectors_start(GMVectors* vectors,
                          GMMirroredPointer* vectors_buf,
                          GMEnv* env) {
  GMAssert(vectors != NULL);
  GMAssert(vectors_buf != NULL);
  GMAssert(env != NULL);

  gm_set_matrix_start(vectors_buf, vectors->num_field_local,
                      vectors->num_vector_local, env);
}

/*---------------------------------------------------------------------------*/

void gm_set_vectors_wait(GMEnv* env) {
  GMAssert(env != NULL);

  gm_set_matrix_wait(env);
}

/*===========================================================================*/
/*---Start/end transfer of metrics data from GPU---*/

void gm_get_metrics_start(GMMetrics* metrics,
                          GMMirroredPointer* metrics_buf,
                          GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(metrics_buf != NULL);
  GMAssert(env != NULL);

  gm_get_matrix_start(metrics_buf, metrics->num_vector_local,
                      metrics->num_vector_local, env);
}

/*---------------------------------------------------------------------------*/

void gm_get_metrics_wait(GMEnv* env) {
  GMAssert(env != NULL);

  gm_get_matrix_wait(env);
}

/*===========================================================================*/
/*---CPU-GPU transfer buffer manipulation---*/

void gm_vectors_to_buf(GMVectors* vectors,
                       GMMirroredPointer* vectors_buf,
                       GMEnv* env) {
  GMAssert(vectors != NULL);
  GMAssert(vectors_buf != NULL);
  GMAssert(env != NULL);

  if (!Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    return;
  }

  int i = 0;
  /*---Copy vectors into GPU buffers if needed---*/
  for (i = 0; i < vectors->num_vector_local; ++i) {
    int k = 0;
    for (k = 0; k < vectors->num_field_local; ++k) {
      ((GMFloat*)vectors_buf->h)[k + vectors->num_field_local * i] =
          GMVectors_float_get(vectors, k, i, env);
    }
  }
}

/*===========================================================================*/
/*---Start calculation of numerators, 2-way generic---*/

void gm_compute_numerators_2way_start(GMVectors* vectors_left,
                                      GMVectors* vectors_right,
                                      GMMetrics* numerators,
                                      GMMirroredPointer* vectors_left_buf,
                                      GMMirroredPointer* vectors_right_buf,
                                      GMMirroredPointer* numerators_buf,
                                      int j_proc,
                                      _Bool do_compute_triang_only,
                                      GMEnv* env) {
  GMAssert(vectors_left != NULL);
  GMAssert(vectors_right != NULL);
  GMAssert(numerators != NULL);
  GMAssert(env != NULL);
  GMAssert(j_proc >= 0 && j_proc < Env_num_proc_vector(env));

  switch (Env_metric_type(env)) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      /*----------------------------------------*/
      gm_compute_czekanowski_numerators_2way_start(
          vectors_left, vectors_right, numerators, vectors_left_buf,
          vectors_right_buf, numerators_buf, j_proc, do_compute_triang_only,
          env);
    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CCC: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    default:
      /*----------------------------------------*/
      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/
}

/*===========================================================================*/
/*---Combine nums and denoms on CPU to get final result, 2-way generic---*/

void gm_compute_2way_combine(GMMetrics* metrics,
                             GMMirroredPointer* metrics_buf,
                             GMVectorSums* vector_sums_left,
                             GMVectorSums* vector_sums_right,
                             int j_proc,
                             _Bool do_compute_triang_only,
                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vector_sums_left != NULL);
  GMAssert(vector_sums_right != NULL);
  GMAssert(env != NULL);
  GMAssert(j_proc >= 0 && j_proc < Env_num_proc_vector(env));

  switch (Env_metric_type(env)) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      /*----------------------------------------*/

      gm_compute_czekanowski_2way_combine(metrics, metrics_buf,
                                          (GMFloat*)vector_sums_left->data,
                                          (GMFloat*)vector_sums_right->data,
                                          j_proc, do_compute_triang_only, env);
    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CCC: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    default:
      /*----------------------------------------*/
      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/
}

/*===========================================================================*/
/*---Start calculation of numerators, 2-way Czekanowski---*/

void gm_compute_czekanowski_numerators_2way_start(
    GMVectors* vectors_left,
    GMVectors* vectors_right,
    GMMetrics* numerators,
    GMMirroredPointer* vectors_left_buf,
    GMMirroredPointer* vectors_right_buf,
    GMMirroredPointer* numerators_buf,
    int j_proc,
    _Bool do_compute_triang_only,
    GMEnv* env) {
  GMAssert(vectors_left != NULL);
  GMAssert(vectors_right != NULL);
  GMAssert(numerators != NULL);
  GMAssert(env != NULL);
  GMAssert(j_proc >= 0 && j_proc < Env_num_proc_vector(env));

  /*----------------------------------------*/
  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    /*----------------------------------------*/

    GMInsist(env, Env_num_proc_field(env) == 1
      ? "num_proc_field>1 for CPU case not supported" : 0);

    /*---Perform pseudo matrix-matrix product---*/

    int j = 0;
    for (j = 0; j < numerators->num_vector_local; ++j) {
      const int i_max =
          do_compute_triang_only ? j : numerators->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        GMFloat numerator = 0;
        int field_local = 0;
        for (field_local = 0; field_local < vectors_left->num_field_local;
             ++field_local) {
          const GMFloat value1 =
              GMVectors_float_get(vectors_left, field_local, i, env);
          const GMFloat value2 =
              GMVectors_float_get(vectors_right, field_local, j, env);
          numerator += value1 < value2 ? value1 : value2;
        } /*---for k---*/
        GMMetrics_float_set_all2all_2(numerators, i, j, j_proc, numerator, env);
      } /*---for i---*/
    }   /*---for j---*/

    /*----------------------------------------*/
  } else /* if (Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) */ {
    /*----------------------------------------*/

    /*---Initialize result matrix to zero (apparently magma requires)---*/

    gm_magma_set_matrix_zero_start(numerators_buf, numerators->num_vector_local,
                                   numerators->num_vector_local, env);

    /*---Perform pseudo matrix-matrix product---*/

    /* .63 / 1.56 */
    magma_gemm_start(vectors_left->num_vector_local,
                  vectors_left->num_vector_local, vectors_left->num_field_local,
                  vectors_left_buf->d, vectors_left->num_field_local,
                  vectors_right_buf->d, vectors_left->num_field_local,
                  numerators_buf->d, vectors_left->num_vector_local);

    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*===========================================================================*/
/*---Start calculation of numerators, 3-way Czekanowski---*/

void gm_compute_czekanowski_numerators_3way_start(
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* metrics,
    GMMirroredPointer* vectors_i_buf,
    GMMirroredPointer* vectors_j_buf,
    GMMirroredPointer* vectors_k_buf,
    int j_proc,
    int k_proc,
    GMEnv* env) {
  GMAssert(vectors_i != NULL);
  GMAssert(vectors_j != NULL);
  GMAssert(vectors_k != NULL);
  GMAssert(metrics != NULL);
  GMAssert(vectors_i_buf != NULL);
  GMAssert(vectors_j_buf != NULL);
  GMAssert(vectors_k_buf != NULL);
  GMAssert(env != NULL);
  GMAssert(j_proc >= 0 && j_proc < Env_num_proc_vector(env));
  GMAssert(k_proc >= 0 && k_proc < Env_num_proc_vector(env));
  GMAssert(!(Env_proc_num_vector(env) == j_proc &&
             Env_proc_num_vector(env) != k_proc));
  GMAssert(!(Env_proc_num_vector(env) == k_proc &&
             Env_proc_num_vector(env) != j_proc));

  /*---Initializations---*/

  const int numvec = metrics->num_vector_local;
  const int numfield = vectors_i->num_field_local;

  const int i_proc = Env_proc_num_vector(env);

  const _Bool is_part1 = i_proc == j_proc && j_proc == k_proc;
  const _Bool is_part3 =
      i_proc != j_proc && j_proc != k_proc && i_proc != k_proc;

  /*---Get specification of region to be computed for Part 3---*/

  const int section_axis =
      gm_metrics_3way_section_axis(metrics, i_proc, j_proc, k_proc, env);
  const int section_num =
      gm_metrics_3way_section_num(metrics, i_proc, j_proc, k_proc, env);

  /*---Define bounding box containing region to be computed---*/

  const int i_lb = is_part3 && section_axis == 0 ? (section_num * numvec)/6 : 0;

  const int j_lb = is_part3 && section_axis == 1 ? (section_num * numvec)/6 : 0;

  const int k_lb = is_part3 && section_axis == 2 ? (section_num * numvec)/6 : 0;

  const int i_ub = is_part3 && section_axis == 0
                       ? ((section_num + 1) * numvec) / 6
                       : numvec;

  const int j_ub = is_part3 && section_axis == 1
                       ? ((section_num + 1) * numvec) / 6
                       : numvec;

  const int k_ub = is_part3 && section_axis == 2
                       ? ((section_num + 1) * numvec) / 6
                       : numvec;

  /*----------------------------------------*/
  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU && ! Env_all2all(env)) {
  /*----------------------------------------*/

    GMInsist(env, Env_num_proc_field(env) == 1
      ? "num_proc_field>1 for CPU case not supported" : 0);

    /*---No off-proc all2all: compute tetrahedron of values---*/

    int k = 0;
    for (k = 0; k < numvec; ++k) {
      int j = 0;
      for (j = 0; j < k; ++j) {
        int i = 0;
        for (i = 0; i < j; ++i) {
          GMFloat sum = 0;
          int field_local = 0;
          for (field_local = 0; field_local < numfield; ++field_local) {
            const GMFloat val1 = GMVectors_float_get(vectors_i,
                                                     field_local, i, env);
            const GMFloat val2 = GMVectors_float_get(vectors_j,
                                                     field_local, j, env);
            const GMFloat val3 = GMVectors_float_get(vectors_k,
                                                     field_local, k, env);
            GMFloat min12 = val1 < val2 ? val1 : val2;
            sum += min12;
            sum += val1 < val3 ? val1 : val3;
            sum += val2 < val3 ? val2 : val3;
            sum -= min12 < val3 ? min12 : val3;
          } /*---for field_local---*/
          GMMetrics_float_set_3(metrics, i, j, k, sum, env);
        }
      }
    }

  /*----------------------------------------*/
  } else if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
  /*----------------------------------------*/

    GMInsist(env, Env_num_proc_field(env) == 1
      ? "num_proc_field>1 for CPU case not supported" : 0);

    /*---Compute tetrahedron, triang prism or block section---*/

    int k = 0;
    for (k = k_lb; k < k_ub; ++k) {
      const int j_max = is_part3 ? j_ub : k;
      int j = 0;
      for (j = j_lb; j < j_max; ++j) {
        const int i_max = is_part1 ? j : i_ub;
        int i = 0;
        for (i = i_lb; i < i_max; ++i) {
          GMFloat numerator = 0;
          int field_local = 0;
          for (field_local = 0; field_local < numfield; ++field_local) {
            const GMFloat val1 = GMVectors_float_get(vectors_i,
                                                     field_local, i, env);
            const GMFloat val2 = GMVectors_float_get(vectors_j,
                                                     field_local, j, env);
            const GMFloat val3 = GMVectors_float_get(vectors_k,
                                                     field_local, k, env);
            const GMFloat min_ij = val1 < val2 ? val1 : val2;
            const GMFloat min_ik = val1 < val3 ? val1 : val3;
            const GMFloat min_jk = val2 < val3 ? val2 : val3;
            const GMFloat min_ijk = min_ij < val3 ? min_ij : val3;
            numerator += min_ij + min_ik + min_jk - min_ijk;
          } /*---for field_local---*/
          GMMetrics_float_set_all2all_3(metrics, i, j, k, j_proc, k_proc,
                                        numerator, env);
        }
      }
    }

  /*----------------------------------------*/
  } else /* if (Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) */ {
  /*----------------------------------------*/

    /*----------------------------------------*/
    /*---First get the required 2-way ij, jk, ik metrics---*/
    /*----------------------------------------*/

    /*--------------------*/
    /*---Compute i_proc - j_proc minproduct---*/
    /*--------------------*/

    /*---Allocate magma CPU/GPU memory for M = X^T minprod X---*/

    GMMirroredPointer matM_ij_buf_value =
        gm_malloc_magma(numvec * (size_t)numvec, env);  // M = X^T minprod X
    GMMirroredPointer* const matM_ij_buf = &matM_ij_buf_value;

    GMMirroredPointer mat_buf_tmp =
        gm_malloc_magma(numvec * (size_t)numvec, env);

    if (GM_BOOL_TRUE) {
      GMMirroredPointer* matM_ij_buf_local = Env_num_proc_field(env) == 1 ?
        matM_ij_buf : &mat_buf_tmp;

      /*---Initialize result matrix to zero (apparently magma requires)---*/

      gm_magma_set_matrix_zero_start(matM_ij_buf_local, numvec, numvec, env);

      /*---Perform pseudo matrix-matrix min product for M = X^T minprod X---*/

      magma_gemm_start(numvec, numvec, numfield, vectors_i_buf->d, numfield,
                    vectors_j_buf->d, numfield, matM_ij_buf_local->d, numvec);
      gm_compute_wait(env);

      /*---Copy matM_ij from GPU---*/

      gm_get_matrix_start(matM_ij_buf_local, numvec, numvec, env);
      gm_get_matrix_wait(env);

      if (Env_num_proc_field(env) > 1) {
        int mpi_code = 0;
        mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
        mpi_code = MPI_Allreduce(matM_ij_buf_local->h, matM_ij_buf->h,
                     numvec*(size_t)numvec, GM_MPI_FLOAT, MPI_SUM,
                     Env_mpi_comm_field(env));
        GMAssert(mpi_code == MPI_SUCCESS);
      }
    }

    /*--------------------*/
    /*---Compute j_proc - k_proc minproduct---*/
    /*--------------------*/

    /*---Need to compute only if not identical to already computed values---*/

    GMMirroredPointer matM_jk_buf_value =
        !is_part1 ? gm_malloc_magma(numvec * (size_t)numvec, env)
                     : GMMirroredPointer_null();

    GMMirroredPointer* const matM_jk_buf =
        !is_part1 ? &matM_jk_buf_value : matM_ij_buf;

    if (!is_part1) {
      GMMirroredPointer* matM_jk_buf_local = Env_num_proc_field(env) == 1 ?
        matM_jk_buf : &mat_buf_tmp;

      gm_magma_set_matrix_zero_start(matM_jk_buf_local, numvec, numvec, env);

      magma_gemm_start(numvec, numvec, numfield, vectors_j_buf->d, numfield,
                     vectors_k_buf->d, numfield, matM_jk_buf_local->d, numvec);
      gm_compute_wait(env);

      gm_get_matrix_start(matM_jk_buf_local, numvec, numvec, env);
      gm_get_matrix_wait(env);

      if (Env_num_proc_field(env) > 1) {
        int mpi_code = 0;
        mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
        mpi_code = MPI_Allreduce(matM_jk_buf_local->h, matM_jk_buf->h,
                     numvec*(size_t)numvec, GM_MPI_FLOAT, MPI_SUM,
                     Env_mpi_comm_field(env));
        GMAssert(mpi_code == MPI_SUCCESS);
      }
    }

    /*--------------------*/
    /*---Compute k_proc - i_proc minproduct---*/
    /*--------------------*/

    /*---Need to compute only if not identical to already computed values---*/

    /*---NOTE: for Part 3, this is indexed directly as (k,i).
         Otherwise, it is indexed through an alias as (i,k)---*/

    GMMirroredPointer matM_kik_buf_value =
        is_part3 ? gm_malloc_magma(numvec * (size_t)numvec, env)
                          : GMMirroredPointer_null();

    GMMirroredPointer* const matM_kik_buf =
        is_part3 ? &matM_kik_buf_value : matM_ij_buf;

    if (is_part3) {
      GMMirroredPointer* matM_kik_buf_local = Env_num_proc_field(env) == 1 ?
        matM_kik_buf : &mat_buf_tmp;

      gm_magma_set_matrix_zero_start(matM_kik_buf_local, numvec, numvec, env);

      magma_gemm_start(numvec, numvec, numfield, vectors_k_buf->d, numfield,
                    vectors_i_buf->d, numfield, matM_kik_buf_local->d, numvec);
      gm_compute_wait(env);

      gm_get_matrix_start(matM_kik_buf_local, numvec, numvec, env);
      gm_get_matrix_wait(env);

      if (Env_num_proc_field(env) > 1) {
        int mpi_code = 0;
        mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
        mpi_code = MPI_Allreduce(matM_kik_buf_local->h, matM_kik_buf->h,
                     numvec*(size_t)numvec, GM_MPI_FLOAT, MPI_SUM,
                     Env_mpi_comm_field(env));
        GMAssert(mpi_code == MPI_SUCCESS);
      }
    }

    /*----------------------------------------*/
    /*---Now compute ijk piece, via an outer loop over j values---*/
    /*----------------------------------------*/

    /*---Allocate magma CPU/GPU memory for matrices V and B---*/
    /*
       V = elementwise min of one vector with the rest of the vectors.
       The for the jth iteration the ith column of V is the elementwise min
         of vectors i and j
       B = X^T minprod V = three way min product
    */
    GMMirroredPointer matV_buf = gm_malloc_magma(numvec*(size_t)numfield, env);
    GMMirroredPointer matB_buf = gm_malloc_magma(numvec*(size_t)numvec, env);

    /*---Set up pointers to permute the access of axes for Part 3---*/
    /*---We use capitals I, J, K here to denote the permuted axes---*/

    const _Bool sax0 = section_axis==0;
    const _Bool sax1 = section_axis==1;
    const _Bool sax2 = section_axis==2;

    /* clang-format off */
    GMMirroredPointer* const vectors_I_buf = !is_part3 ? vectors_i_buf :
                                                  sax0 ? vectors_k_buf :
                                                  sax1 ? vectors_i_buf :
                                                  sax2 ? vectors_j_buf : 0;

    GMMirroredPointer* const vectors_J_buf = !is_part3 ? vectors_j_buf :
                                                  sax0 ? vectors_i_buf :
                                                  sax1 ? vectors_j_buf :
                                                  sax2 ? vectors_k_buf : 0;

    GMMirroredPointer* const vectors_K_buf = !is_part3 ? vectors_k_buf :
                                                  sax0 ? vectors_j_buf :
                                                  sax1 ? vectors_k_buf :
                                                  sax2 ? vectors_i_buf : 0;

    /*---NOTE: must pay attention that these permuted matrices
         are indexed the right way by the permuted indices---*/

    GMMirroredPointer* const matM_IJ_buf  = !is_part3 ? matM_ij_buf  :
                                                 sax0 ? matM_kik_buf :
                                                 sax1 ? matM_ij_buf  :
                                                 sax2 ? matM_jk_buf  : 0;

    GMMirroredPointer* const matM_JK_buf  = !is_part3 ? matM_jk_buf  :
                                                 sax0 ? matM_ij_buf  :
                                                 sax1 ? matM_jk_buf  :
                                                 sax2 ? matM_kik_buf : 0;

    GMMirroredPointer* const matM_KIK_buf = !is_part3 ? matM_kik_buf :
                                                 sax0 ? matM_jk_buf  :
                                                 sax1 ? matM_kik_buf :
                                                 sax2 ? matM_ij_buf  : 0;
    /* clang-format on */

    /*---Process all combinations starting with j, i, k---*/

    const int J_min = is_part3 ? (section_num + 0) * numvec / 6 : 0;
    const int J_max = is_part3 ? (section_num + 1) * numvec / 6 : numvec;
    int J = 0;

    /*--------------------*/
    /*---J loop---*/
    /*--------------------*/

    for (J = J_min; J < J_max; ++J) {

      const int I_min = 0;
      const int I_max = is_part1 ? J : numvec;
      if (I_min >= I_max) {
        continue;
      }

      const int K_min = is_part3 ? 0 : J + 1;
      const int K_max = numvec;
      if (K_min >= K_max) {
        continue;
      }

      /*---Populate leading columns of matV---*/

      int I = 0;
      for (I = I_min; I < I_max; ++I) {
        // Compare columns x_i and x_j element-wise
        int field_local = 0;
        for (field_local = 0; field_local < numfield; ++field_local) {
          const GMFloat a =
              ((GMFloat*)(vectors_I_buf->h))[field_local + numfield * I];
          const GMFloat b =
              ((GMFloat*)(vectors_J_buf->h))[field_local + numfield * J];
          ((GMFloat*)(matV_buf.h))[field_local + numfield * I] = a < b ? a : b;
        }  //---for field_local---//
      }    //---for i---//

      /*---Send matrix matV to GPU---*/

      gm_set_matrix_start(&matV_buf, numfield, I_max, env);
      gm_set_matrix_wait(env);

      /*---Initialize result matrix to zero (apparently magma requires)---*/

      GMMirroredPointer* matB_buf_local = Env_num_proc_field(env) == 1 ?
        &matB_buf : &mat_buf_tmp;

      gm_magma_set_matrix_zero_start(matB_buf_local, numvec, I_max, env);

      /*---Perform matrix-matrix product matB = matV^T minprod X---*/

      magma_gemm_start(I_max, numvec, numfield, matV_buf.d, numfield,
                       vectors_K_buf->d, numfield, matB_buf_local->d, I_max);
      gm_compute_wait(env);

      /*---Copy result matrix matB from GPU---*/

      gm_get_matrix_start(matB_buf_local, I_max, numvec, env);
      gm_get_matrix_wait(env);

      if (Env_num_proc_field(env) > 1) {
        int mpi_code = 0;
        mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
        mpi_code = MPI_Allreduce(matB_buf_local->h, matB_buf.h,
                     numvec*(size_t)numvec, GM_MPI_FLOAT, MPI_SUM,
                     Env_mpi_comm_field(env));
        GMAssert(mpi_code == MPI_SUCCESS);
      }

      /*---Compute numerators using 2-way pieces and ijk piece---*/

      /*----------*/
      if (!Env_all2all(env)) {
      /*----------*/

        for (I = I_min; I < I_max; ++I) {
          const GMFloat min_IJ = ((GMFloat*)(matM_IJ_buf->h))[I + numvec*J];
          int K = 0;
          for (K = K_min; K < K_max; ++K) {
            const GMFloat min_JK  = ((GMFloat*)(matM_JK_buf->h))[J + numvec*K];
            const GMFloat min_KIK = ((GMFloat*)(matM_KIK_buf->h))[K + numvec*I];
            // sum of mins vectors i, j, and k is matB(k,i)
            const GMFloat min_IJK = ((GMFloat*)(matB_buf.h))[I + I_max*K];
            const GMFloat numerator = min_IJ + min_JK + min_KIK - min_IJK;
            const int i = I;
            const int j = J;
            const int k = K;
            GMMetrics_float_set_3(metrics, i, j, k, numerator, env);
          } /*---for K---*/
        }   /*---for I---*/

      /*----------*/
      } else /*---if (Env_all2all(env))---*/ {
      /*----------*/

        for (I = I_min; I < I_max; ++I) {
          const GMFloat min_IJ = ((GMFloat*)(matM_IJ_buf->h))[I + numvec*J];
          int K = 0;
          for (K = K_min; K < K_max; ++K) {
            const GMFloat min_JK = ((GMFloat*)(matM_JK_buf->h))[J + numvec*K];
            const GMFloat min_KIK = is_part3 ?
                                   ((GMFloat*)(matM_KIK_buf->h))[K + numvec*I] :
                                   ((GMFloat*)(matM_KIK_buf->h))[I + numvec*K];
            // sum of mins vectors i, j, and k is matB(k,i)
            const GMFloat min_IJK = ((GMFloat*)(matB_buf.h))[I + I_max*K];
            const GMFloat numerator = min_IJ + min_JK + min_KIK - min_IJK;
            /* clang-format off */
            const int i = !is_part3 ?   I :
                               sax0 ?   J :
                               sax1 ?   I :
                            /* sax2 ?*/ K;
            const int j = !is_part3 ?   J :
                               sax0 ?   K :
                               sax1 ?   J :
                            /* sax2 ?*/ I;
            const int k = !is_part3 ?   K :
                               sax0 ?   I :
                               sax1 ?   K :
                            /* sax2 ?*/ J;
            /* clang-format on */
            GMMetrics_float_set_all2all_3(metrics, i, j, k, j_proc, k_proc,
                                          numerator, env);
          } /*---for K---*/
        }   /*---for I---*/

      /*----------*/
      } /*---if (Env_all2all(env))---*/
      /*----------*/

    } /*---for J---*/

    /*--------------------*/
    /*---Free memory---*/
    /*--------------------*/

    gm_free_magma(&matM_ij_buf_value, env);
    if (!is_part1) {
      gm_free_magma(&matM_jk_buf_value, env);
    }
    if (is_part3) {
      gm_free_magma(&matM_kik_buf_value, env);
    }
    gm_free_magma(&matV_buf, env);
    gm_free_magma(&matB_buf, env);
    gm_free_magma(&mat_buf_tmp, env);

  } /*---if GPU---*/
}

/*===========================================================================*/
/*---Combine nums and denoms on CPU to get final result, 2-way Czek---*/

void gm_compute_czekanowski_2way_combine(
    GMMetrics* metrics,
    GMMirroredPointer* metrics_buf,
    GMFloat* __restrict__ vector_sums_left,
    GMFloat* __restrict__ vector_sums_right,
    int j_proc,
    _Bool do_compute_triang_only,
    GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vector_sums_left != NULL);
  GMAssert(vector_sums_right != NULL);
  GMAssert(env != NULL);
  GMAssert(j_proc >= 0 && j_proc < Env_num_proc_vector(env));

  /*---For CPU case, copy numerator out of metrics struct which is temporarily
       holding numerators.
       For GPU case, directly access the metrics_buf holding the numerators.
  ---*/

  const _Bool are_vector_sums_aliased = vector_sums_left == vector_sums_right;

  /*----------------------------------------*/
  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU && Env_all2all(env)) {
    /*----------------------------------------*/

    int j = 0;
    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        const GMFloat numerator =
            GMMetrics_float_get_all2all_2(metrics, i, j, j_proc, env);
        /*---Don't use two pointers pointing to the same thing---*/
        const GMFloat denominator =
            are_vector_sums_aliased
                ? vector_sums_left[i] + vector_sums_left[j]
                : vector_sums_left[i] + vector_sums_right[j];
        GMMetrics_float_set_all2all_2(metrics, i, j, j_proc,
                                      2 * numerator / denominator, env);
      } /*---for i---*/
    }   /*---for j---*/
        /*---TODO: here and elsewhere check for unlikely case denom is/nearly
         * zero---*/

    /*----------------------------------------*/
  } else if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    /*----------------------------------------*/

    int j = 0;
    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        const GMFloat numerator =
            GMMetrics_float_get_all2all_2(metrics, i, j, j_proc, env);
        /*---Don't use two different pointers pointing to the same thing---*/
        const GMFloat denominator = vector_sums_left[i] + vector_sums_left[j];
        GMMetrics_float_set_2(metrics, i, j, 2 * numerator / denominator, env);
      } /*---for i---*/
    }   /*---for j---*/

    /*----------------------------------------*/
  } else if (Env_all2all(env)) {
    /*----------------------------------------*/

    int j = 0;
    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        const GMFloat numerator =
            ((GMFloat*)metrics_buf->h)[i + metrics->num_vector_local * j];
        /*---Don't use two pointers pointing to the same thing---*/
        const GMFloat denominator =
            are_vector_sums_aliased
                ? vector_sums_left[i] + vector_sums_left[j]
                : vector_sums_left[i] + vector_sums_right[j];
        GMMetrics_float_set_all2all_2(metrics, i, j, j_proc,
                                      2 * numerator / denominator, env);
      } /*---for i---*/
    }   /*---for j---*/

    /*----------------------------------------*/
  } else {
    /*----------------------------------------*/

    int j = 0;
    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        const GMFloat numerator =
            ((GMFloat*)metrics_buf->h)[i + metrics->num_vector_local * j];
        /*---Don't use two different pointers pointing to the same thing---*/
        const GMFloat denominator = vector_sums_left[i] + vector_sums_left[j];
        GMMetrics_float_set_2(metrics, i, j, 2 * numerator / denominator, env);
      } /*---for i---*/
    }   /*---for j---*/

    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*===========================================================================*/
/*---Combine nums and denoms on CPU to get final result, 3-way Czek---*/

void gm_compute_czekanowski_3way_combine(GMMetrics* metrics,
                                         GMFloat* __restrict__ vector_sums_i,
                                         GMFloat* __restrict__ vector_sums_j,
                                         GMFloat* __restrict__ vector_sums_k,
                                         int j_proc,
                                         int k_proc,
                                         GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vector_sums_i != NULL);
  GMAssert(vector_sums_j != NULL);
  GMAssert(vector_sums_k != NULL);
  GMAssert(env != NULL);
  GMAssert(j_proc >= 0 && j_proc < Env_num_proc_vector(env));
  GMAssert(k_proc >= 0 && k_proc < Env_num_proc_vector(env));
  GMAssert(Env_proc_num_vector(env) != j_proc || j_proc == k_proc);
  GMAssert(Env_proc_num_vector(env) != k_proc || j_proc == k_proc);

  const int i_proc = Env_proc_num_vector(env);

  /*----------------------------------------*/
  if (Env_all2all(env)) {
  /*----------------------------------------*/

    /*----------------------------------------*/
    if (i_proc == j_proc && j_proc == k_proc) {
    /*----------------------------------------*/

      int k = 0;
      for (k = 0; k < metrics->num_vector_local; ++k) {
        int j = 0;
        for (j = 0; j < k; ++j) {
          int i = 0;
          for (i = 0; i < j; ++i) {
            const GMFloat numerator = GMMetrics_float_get_all2all_3(
                metrics, i, j, k, j_proc, k_proc, env);
            const GMFloat denominator =
                vector_sums_i[i] + vector_sums_i[j] + vector_sums_i[k];
            GMMetrics_float_set_all2all_3(metrics, i, j, k, j_proc, k_proc,
                                          3 * numerator / (2 * denominator),
                                          env);
          } /*---for i---*/
        }   /*---for j---*/
      }     /*---for k---*/

    /*----------------------------------------*/
    } else if (j_proc == k_proc) {
    /*----------------------------------------*/

      int k = 0;
      for (k = 0; k < metrics->num_vector_local; ++k) {
        int j = 0;
        for (j = 0; j < k; ++j) {
          int i = 0;
          for (i = 0; i < metrics->num_vector_local; ++i) {
            const GMFloat numerator = GMMetrics_float_get_all2all_3(
                metrics, i, j, k, j_proc, k_proc, env);
            const GMFloat denominator =
                vector_sums_i[i] + vector_sums_j[j] + vector_sums_j[k];
            GMMetrics_float_set_all2all_3(metrics, i, j, k, j_proc, k_proc,
                                          3 * numerator / (2 * denominator),
                                          env);
          } /*---for i---*/
        }   /*---for j---*/
      }     /*---for k---*/


    /*----------------------------------------*/
    } else /*---i_proc != j_proc && i_proc != k_proc && j_proc != k_proc---*/ {
    /*----------------------------------------*/

      const int numvec = metrics->num_vector_local;

      /*---Get specification of region to be computed for Part 3---*/

      const int section_axis =
          gm_metrics_3way_section_axis(metrics, i_proc, j_proc, k_proc, env);
      const int section_num =
          gm_metrics_3way_section_num(metrics, i_proc, j_proc, k_proc, env);

      /*---Define bounding box containing region to be computed---*/

      const int i_lb = section_axis == 0 ? (section_num * numvec) / 6 : 0;

      const int j_lb = section_axis == 1 ? (section_num * numvec) / 6 : 0;

      const int k_lb = section_axis == 2 ? (section_num * numvec) / 6 : 0;

      const int i_ub =
          section_axis == 0 ? ((section_num + 1) * numvec) / 6 : numvec;

      const int j_ub =
          section_axis == 1 ? ((section_num + 1) * numvec) / 6 : numvec;

      const int k_ub =
          section_axis == 2 ? ((section_num + 1) * numvec) / 6 : numvec;

      int k = 0;
      for (k = k_lb; k < k_ub; ++k) {
        int j = 0;
        for (j = j_lb; j < j_ub; ++j) {
          int i = 0;
          for (i = i_lb; i < i_ub; ++i) {
            const GMFloat numerator = GMMetrics_float_get_all2all_3(
                metrics, i, j, k, j_proc, k_proc, env);
            const GMFloat denominator =
                vector_sums_i[i] + vector_sums_j[j] + vector_sums_k[k];
            const GMFloat value =  ((GMFloat)3) * numerator /
                                  (((GMFloat)2) * denominator);
            GMMetrics_float_set_all2all_3(metrics, i, j, k, j_proc, k_proc,
                                          value, env);
          } /*---for i---*/
        }   /*---for j---*/
      }     /*---for k---*/
    }

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

    int i = 0;
    for (i = 0; i < metrics->num_vector_local; ++i) {
      int j = 0;
      for (j = i + 1; j < metrics->num_vector_local; ++j) {
        int k = 0;
        for (k = j + 1; k < metrics->num_vector_local; ++k) {
          const GMFloat numerator =
              GMMetrics_float_get_3(metrics, i, j, k, env);
          /*---Don't use two different pointers pointing to the same thing---*/
          const GMFloat denominator =
              vector_sums_i[i] + vector_sums_i[j] + vector_sums_i[k];
          GMMetrics_float_set_3(metrics, i, j, k,
                                3 * numerator / (2 * denominator), env);
        } /*---for k---*/
      }   /*---for j---*/
    }     /*---for i---*/

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
