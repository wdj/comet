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
#include "vector_sums.h"
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

/*---------------------------------------------------------------------------*/

void gm_compute_vector_sums(GMVectors* vectors,
                            GMVectorSums* vector_sums,
                            GMEnv* env) {
  GMAssert(vectors != NULL);
  GMAssert(vector_sums != NULL);
  GMAssert(env != NULL);

  switch (env->metric_type) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
    /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
    /*----------------------------------------*/

      gm_compute_float_vector_sums(vectors, (GMFloat*)vector_sums->data, env);

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

  if (env->compute_method != GM_COMPUTE_METHOD_GPU ) {
    return;
  }

  switch (env->metric_type) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
    /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
    /*----------------------------------------*/

      magma_minproduct_int_t magma_code = 0;
      magma_code = magma_code*1; /*---Avoid unused variable warning---*/

      magma_code = magma_minproduct_init();
      GMAssert(magma_code == MAGMA_minproduct_SUCCESS);
      /*---need this -- see http://on-demand.gputechconf.com/gtc/2014/presentations/S4158-cuda-streams-best-practices-common-pitfalls.pdf page 14 ---*/
      magma_code = magma_minproductblasSetKernelStream(env->stream_compute);
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

  if (env->compute_method != GM_COMPUTE_METHOD_GPU ) {
    return;
  }

  switch (env->metric_type) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
    /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
    /*----------------------------------------*/

      magma_minproduct_int_t magma_code = 0;
      magma_code = magma_code*1; /*---Avoid unused variable warning---*/

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
  GMAssert(n+1 >= 1);
  GMAssert(env != NULL);

  GMMirroredPointer p = GMMirroredPointer_null();

  if (env->compute_method != GM_COMPUTE_METHOD_GPU ) {
    return p;
  }

  switch (env->metric_type) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
    /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
    /*----------------------------------------*/

      magma_minproduct_int_t magma_code = 0;
      magma_code = magma_code*1; /*---Avoid unused variable warning---*/

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

  if (env->compute_method != GM_COMPUTE_METHOD_GPU ) {
    return;
  }

  switch (env->metric_type) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
    /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
    /*----------------------------------------*/

      magma_minproduct_int_t magma_code = 0;
      magma_code = magma_code*1; /*---Avoid unused variable warning---*/

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

void gm_set_vectors_start(GMVectors* vectors,
                          GMMirroredPointer* vectors_buf,
                          GMEnv* env) {
  GMAssert(vectors != NULL);
  GMAssert(vectors_buf != NULL);
  GMAssert(env != NULL);

  if (! env->compute_method == GM_COMPUTE_METHOD_GPU) {
    return;
  }

  const int numvec = vectors->num_vector_local;
  const int numfield = vectors->num_field;

  /*---Send vectors to GPU---*/

  GMEnv_initialize_streams(env);
#ifdef FP_PRECISION_DOUBLE
  magma_minproduct_dsetmatrix_async(numfield, numvec,
                   (GMFloat*)vectors_buf->h, numfield,
                   (GMFloat*)vectors_buf->d, numfield, env->stream_vectors);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproduct_ssetmatrix_async(numfield, numvec,
                   (GMFloat*)vectors_buf->h, numfield,
                   (GMFloat*)vectors_buf->d, numfield, env->stream_vectors);
#endif
}

/*---------------------------------------------------------------------------*/

void gm_set_vectors_wait(GMEnv* env) {
  GMAssert(env != NULL);

  if (! env->compute_method == GM_COMPUTE_METHOD_GPU) {
    return;
  }

  GMEnv_initialize_streams(env);
  cudaStreamSynchronize( env->stream_vectors );
  GMAssert(GMEnv_cuda_last_call_succeeded(env));
}

/*===========================================================================*/
/*---Start/end transfer of generic matrix to GPU---*/
//---Added by James Nance---//
void gm_set_matrix_start(GMMirroredPointer* matrix_buf,
                                 int mat_dim1, int mat_dim2, GMEnv* env) {
  GMAssert(matrix_buf != NULL);
  GMAssert(env != NULL);

  if (! env->compute_method == GM_COMPUTE_METHOD_GPU) {
    return;
  }

  /*---Send vectors to GPU---*/

  GMEnv_initialize_streams(env);
#ifdef FP_PRECISION_DOUBLE
  magma_minproduct_dsetmatrix_async(mat_dim1, mat_dim2,
                   (GMFloat*)matrix_buf->h, mat_dim1,
                   (GMFloat*)matrix_buf->d, mat_dim1, env->stream_vectors);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproduct_ssetmatrix_async(mat_dim1, mat_dim2,
                   (GMFloat*)matrix_buf->h, mat_dim1,
                   (GMFloat*)matrix_buf->d, mat_dim1, env->stream_vectors);
#endif

}

/*---------------------------------------------------------------------------*/

void gm_set_matrix_wait(GMEnv* env) {
  GMAssert(env != NULL);

  if (! env->compute_method == GM_COMPUTE_METHOD_GPU) {
    return;
  }

  GMEnv_initialize_streams(env);
  cudaStreamSynchronize( env->stream_vectors );//Does this line need to be changed?? -JN
  GMAssert(GMEnv_cuda_last_call_succeeded(env));
}

/*===========================================================================*/
/*---Start/end transfer of metrics data from GPU---*/

void gm_get_metrics_start(GMMetrics* metrics,
                          GMMirroredPointer* metrics_buf,
                          GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(metrics_buf != NULL);
  GMAssert(env != NULL);

  if (! env->compute_method == GM_COMPUTE_METHOD_GPU) {
    return;
  }

  const int numvec = metrics->num_vector_local;

  /*---Send vectors to GPU---*/

  GMEnv_initialize_streams(env);
#ifdef FP_PRECISION_DOUBLE
  magma_minproduct_dgetmatrix_async(numvec, numvec,
                        (GMFloat*)metrics_buf->d, numvec,
                        (GMFloat*)metrics_buf->h, numvec, env->stream_metrics);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproduct_sgetmatrix_async(numvec, numvec,
                        (GMFloat*)metrics_buf->d, numvec,
                        (GMFloat*)metrics_buf->h, numvec, env->stream_metrics);
#endif
}

/*---------------------------------------------------------------------------*/

void gm_get_metrics_wait(GMEnv* env) {
  GMAssert(env != NULL);

  if (! env->compute_method == GM_COMPUTE_METHOD_GPU) {
    return;
  }

  GMEnv_initialize_streams(env);
  cudaStreamSynchronize( env->stream_metrics );
  GMAssert(GMEnv_cuda_last_call_succeeded(env));
}

/*===========================================================================*/
/*---Start/end transfer of generic matrix from GPU---*/
//---Added by James Nance---//
void gm_get_matrix_start(GMMirroredPointer* matrix_buf,
                         int mat_dim1, int mat_dim2,
                         GMEnv* env) {
  
  GMAssert(matrix_buf != NULL);
  GMAssert(env != NULL);

  if (! env->compute_method == GM_COMPUTE_METHOD_GPU) {
    return;
  }

  /*---Send vectors to GPU---*/

  GMEnv_initialize_streams(env);
#ifdef FP_PRECISION_DOUBLE
  magma_minproduct_dgetmatrix_async(mat_dim1, mat_dim2,
                        (GMFloat*)matrix_buf->d, mat_dim1,
                        (GMFloat*)matrix_buf->h, mat_dim1, env->stream_metrics);
#endif
#ifdef FP_PRECISION_SINGLE
  magma_minproduct_sgetmatrix_async(mat_dim1, mat_dim2,
                        (GMFloat*)matrix_buf->d, mat_dim1,
                        (GMFloat*)matrix_buf->h, mat_dim1, env->stream_metrics);
#endif
}

/*---------------------------------------------------------------------------*/

void gm_get_matrix_wait(GMEnv* env) {
  GMAssert(env != NULL);

  if (! env->compute_method == GM_COMPUTE_METHOD_GPU) {
    return;
  }

  GMEnv_initialize_streams(env);
  cudaStreamSynchronize( env->stream_metrics );
  GMAssert(GMEnv_cuda_last_call_succeeded(env));
}

/*===========================================================================*/
/*---CPU-GPU transfer nbuffer manipulation---*/

void gm_vectors_to_buf(GMVectors* vectors,
                       GMMirroredPointer* vectors_buf,
                       GMEnv* env) {
  GMAssert(vectors != NULL);
  GMAssert(vectors_buf != NULL);
  GMAssert(env != NULL);

  if (! env->compute_method == GM_COMPUTE_METHOD_GPU) {
    return;
  }

  int i = 0;
  /*---Copy vectors into GPU buffers if needed---*/
  for (i = 0; i < vectors->num_vector_local; ++i) {
    int k = 0;
    for (k = 0; k < vectors->num_field; ++k) {
      ((GMFloat*)vectors_buf->h)[k + vectors->num_field*i] =
        GMVectors_float_get(vectors, k, i, env);
    }
  }
}

/*===========================================================================*/
/*---Start/end calculation of numerators---*/

void gm_compute_czekanowski_numerators_start(
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
  GMAssert(j_proc >= 0 && j_proc < env->num_proc);

  /*----------------------------------------*/
  if ( env->compute_method != GM_COMPUTE_METHOD_GPU ) {
  /*----------------------------------------*/

    /*---Perform pseudo matrix-matrix product---*/

    int j = 0;
    for (j = 0; j < numerators->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j
                                               : numerators->num_vector_local;
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
        0.0, 0.0, (GMFloat*)numerators_buf->d, numerators->num_vector_local);

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
         (GMFloat*)vectors_left_buf->d, vectors_left->num_field,
         (GMFloat*)vectors_right_buf->d, vectors_left->num_field,
         0.0, (GMFloat*)numerators_buf->d, vectors_left->num_vector_local);

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}


/*---------------------------------------------------------------------------*/
void gm_compute_czekanowski_numerators_3way_start(
                                      GMVectors* vectors_left,
                                      GMVectors* vectors_right,
                                      GMMetrics* metrics,
                                      GMMirroredPointer* vectors_left_buf,
                                      GMMirroredPointer* vectors_right_buf,
                                      int j_proc,
                                      _Bool do_compute_triang_only,
                                      GMEnv* env) {
  GMAssert(vectors_left != NULL);
  GMAssert(vectors_right != NULL);
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(j_proc >= 0 && j_proc < env->num_proc);

  const int numvec = metrics->num_vector_local;
  const int numfield = vectors_left->num_field;
 
  /*----------------------------------------*/
  if ( env->compute_method != GM_COMPUTE_METHOD_GPU ) {
  /*----------------------------------------*/

  /*----------------------------------------*/
  } else /* if (env->compute_method == GM_COMPUTE_METHOD_GPU) */ {
  /*----------------------------------------*/
   
    /*---Allocate magma CPU/GPU memory for M = X^T minprod X---*/ 
    GMMirroredPointer matM_buf =
      gm_malloc_magma( numvec * (size_t)numvec, env); // M = X^T minprod X    

    /*---Initialize result matrix to zero (apparently magma requires)---*/

#ifdef FP_PRECISION_DOUBLE
    magma_minproductblas_dlaset
#endif
#ifdef FP_PRECISION_SINGLE
    magma_minproductblas_slaset
#endif
        (Magma_minproductFull,
        numvec, numvec,
        0.0, 0.0, (GMFloat*)(matM_buf.d), numvec);

    /*---Perform pseudo matrix-matrix min product for M = X^T minprod X---*/
///*
    for (int ri = 0; ri < numvec; ++ri) {
      for (int rj = 0; rj < numfield; ++rj) {
         printf("X(%i,%i) = %f\n",rj,ri,((GMFloat*)(vectors_left_buf->d))[rj + ri * numfield]);
      }
    }
//*/
/* .63 / 1.56 */
#ifdef FP_PRECISION_DOUBLE
    magma_minproductblas_dgemm
#endif
#ifdef FP_PRECISION_SINGLE
    magma_minproductblas_sgemm
#endif
        (Magma_minproductTrans, Magma_minproductNoTrans,
         numvec, numvec, numfield, 1.0,
         (GMFloat*)vectors_left_buf->d, numfield,
         (GMFloat*)vectors_right_buf->d, numfield,
         0.0, (GMFloat*)matM_buf.d, numvec);

    /*---Copy matM from GPU---*/
    gm_get_matrix_start(&matM_buf, numvec, numvec, env);
    gm_get_matrix_wait(env);
     

    for (int ii = 0; ii < numvec; ++ii) {
      for (int jj = 0; jj < numvec; ++jj) {
         printf("M(%i,%i) = %f\n",ii,jj,((GMFloat*)(matM_buf.h))[ii + jj * numvec]);
      }
    }

    /*---Allocate magma CPU/GPU memory for matrices V and B---*/
    /* 
       V = elementwise min of one vector with the rest of the vectors. 
       The for the jth iteration the ith column of V is the elementwise min of vectors i and j
       B = X^T minprod V = three way min product
    */
    GMMirroredPointer matV_buf =
      gm_malloc_magma ( numvec * (size_t)numfield, env); 
    GMMirroredPointer matB_buf =
      gm_malloc_magma ( numvec * (size_t)numvec, env);
    
    /*---Process all combinations starting with j, i, k---*/
    int j = 0;
    int i = 0;
    int k = 0;
    for (j = 1; j < numvec-1; ++j) {
      /*---Populate first j-1 columns of matV---*/
      for (i = 0; i < numvec; ++i) {
        // Compare columns x_i and x_j element-wise
        for (k = 0; k < numfield; ++k) {
          const GMFloat a = ((GMFloat*)(vectors_left_buf->h))[k + numfield * i];
          const GMFloat b = ((GMFloat*)(vectors_right_buf->h))[k + numfield * j];
          ((GMFloat*)(matV_buf.h))[k + i * numfield] = a < b ? a : b;
          //printf("V(%i,%i) = %f\n",k,i,((GMFloat*)(matV_buf.h))[k + i * numfield]);
        }//---for k---//
      }//---for i---//

      /*---Send matrix matV to GPU---*/
      gm_set_matrix_start(&matV_buf, numfield, numvec, env);
      gm_set_matrix_wait(env); 

      /*---Initialize result matrix to zero (apparently magma requires)---*/
#ifdef FP_PRECISION_DOUBLE
      magma_minproductblas_dlaset(Magma_minproductFull, numvec, numvec, 0.0, 0.0,
                                matB_buf.d, numvec);
#endif
#ifdef FP_PRECISION_SINGLE
      magma_minproductblas_slaset(Magma_minproductFull, numvec, numvec, 0.0, 0.0,
                                matB_buf.d, numvec);
#endif

      /*---Perform matrix-matrix product matB = matV^T minprod X---*/
#ifdef FP_PRECISION_DOUBLE
      magma_minproductblas_dgemm
#endif
#ifdef FP_PRECISION_SINGLE
      magma_minproductblas_sgemm
#endif
          (Magma_minproductTrans, Magma_minproductNoTrans,
           numvec, numvec, numfield, 1.0,
           (GMFloat*)matV_buf.d, numfield, 
           (GMFloat*)vectors_left_buf->d, numfield, 
           0.0, (GMFloat*)matB_buf.d, numvec);

       /*---Copy matB from GPU---*/
       gm_get_matrix_start(&matB_buf, numvec, numvec, env);
       gm_get_matrix_wait(env);
  
     /*---Compute numerators---*/
     for (i = 0; i < j; ++i) {
       const GMFloat min_ij = ((GMFloat*)(matM_buf.h))[j + numvec * i];
       for (k = j + 1; k < numvec; ++k) {
         const GMFloat min_ik = ((GMFloat*)(matM_buf.h))[k + numvec * i];
         const GMFloat min_jk = ((GMFloat*)(matM_buf.h))[k + numvec * j];
         // sum of mins vectors i, j, and k is matB(i,k)
         const GMFloat min_ijk = ((GMFloat*)(matB_buf.h))[k + numvec * i];
         const GMFloat numerator = min_ij + min_ik + min_jk - min_ijk;
         //printf("%i,%i,%i . . . numerator = %f\n",i,j,k,numerator);
         GMMetrics_float_set_3(metrics, i, j, k,
                             numerator, env);
        } /*---for k---*/
      }   /*---for i---*/
    }     /*---for j---*/
  
    /*---Free memory---*/
    gm_free_magma( &matM_buf, env);
    gm_free_magma( &matV_buf, env);
    gm_free_magma( &matB_buf, env);
  
  }/*---if GPU---*/
}
/*---------------------------------------------------------------------------*/

void gm_compute_numerators_start( GMVectors* vectors_left,
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
  GMAssert(j_proc >= 0 && j_proc < env->num_proc);

  switch (env->metric_type) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
    /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
    /*----------------------------------------*/
      if ( env-> num_way == 2) {
        gm_compute_czekanowski_numerators_start(vectors_left, vectors_right,
             numerators, vectors_left_buf, vectors_right_buf, numerators_buf,
             j_proc, do_compute_triang_only, env);
      } else { //---3way added by James Nance---//
        gm_compute_czekanowski_numerators_3way_start(vectors_left, vectors_right,
             numerators, vectors_left_buf, vectors_right_buf,
             j_proc, do_compute_triang_only, env); 
      }
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

void gm_compute_numerators_wait(GMEnv* env) {
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
  GMAssert(j_proc >= 0 && j_proc < env->num_proc);

  /*---For CPU case, copy numerator out of metrics struct which is temporarily
       holding numerators.
       For GPU case, directly access the metrics_buf holding the numerators.
  ---*/

  const _Bool are_vector_sums_aliased = vector_sums_left == vector_sums_right;

  /*----------------------------------------*/
  if ( env->compute_method != GM_COMPUTE_METHOD_GPU && env->all2all ) {
  /*----------------------------------------*/

    int j = 0;
    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        const GMFloat numerator = GMMetrics_float_get_all2all_2(
                                                   metrics, i, j, j_proc, env);
        /*---Don't use two pointers pointing to the same thing---*/
        const GMFloat denominator = are_vector_sums_aliased
                       ? vector_sums_left[i] + vector_sums_left[j]
                       : vector_sums_left[i] + vector_sums_right[j];
        GMMetrics_float_set_all2all_2(metrics, i, j, j_proc,
                                      2 * numerator / denominator, env);
      } /*---for i---*/
    }   /*---for j---*/
/*---TODO: here and elsewhere check for unlikely case denom is/nearly zero---*/

  /*----------------------------------------*/
  } else if ( env->compute_method != GM_COMPUTE_METHOD_GPU ) {
  /*----------------------------------------*/
 

    int j = 0;
    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        const GMFloat numerator = GMMetrics_float_get_all2all_2(
                                                   metrics, i, j, j_proc, env);
        /*---Don't use two pointers pointing to the same thing---*/
        const GMFloat denominator = vector_sums_left[i] + vector_sums_left[j];
        GMMetrics_float_set_2(metrics, i, j, 2 * numerator / denominator, env);
      } /*---for i---*/
    }   /*---for j---*/  
  
  /*----------------------------------------*/
  } else if ( env->all2all ) {
  /*----------------------------------------*/

    int j = 0;
    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        const GMFloat numerator = ((GMFloat*)metrics_buf->h)[i +
                                              metrics->num_vector_local * j];
        /*---Don't use two pointers pointing to the same thing---*/
        const GMFloat denominator = are_vector_sums_aliased
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
        const GMFloat numerator = ((GMFloat*)metrics_buf->h)[i +
                                              metrics->num_vector_local * j];
        /*---Don't use two pointers pointing to the same thing---*/
        const GMFloat denominator = vector_sums_left[i] + vector_sums_left[j];
        GMMetrics_float_set_2(metrics, i, j, 2 * numerator / denominator, env);
      } /*---for i---*/
    }   /*---for j---*/


  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*===========================================================================*/
//--Added by James Nance---//
void gm_compute_czekanowski_3way_combine(GMMetrics* metrics,
                                    GMFloat* __restrict__ vector_sums_left,
                                    GMFloat* __restrict__ vector_sums_right,
                                    int j_proc,
                                    _Bool do_compute_triang_only,
                                    GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vector_sums_left != NULL);
  GMAssert(vector_sums_right != NULL);
  GMAssert(env != NULL);
  GMAssert(j_proc >= 0 && j_proc < env->num_proc);
 
  /*----------------------------------------*/
  if ( env->all2all ) {
  /*----------------------------------------*/

    GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/
    int i = 0;
    int j = 0;
    int k = 0;
    for (i = 0; i < metrics->num_vector_local-2; ++i) {
      for (j = i + 1; j < metrics->num_vector_local-1; ++j) {
        for (k = j + 1; k < metrics->num_vector_local; ++k) {
          const GMFloat numerator = GMMetrics_float_get_3(metrics, i, j, k, env);
          const GMFloat denominator = vector_sums_left[i] 
                          + vector_sums_left[j] + vector_sums_left[k];
          printf("%i,%i,%i . . . numerator = %f . . . denominator = %f\n",i,j,k,numerator,denominator);
          GMMetrics_float_set_3(metrics, i, j, k,
                              3 * numerator / (2 * denominator), env);
        } /*---for k---*/
      } /*---for j---*/
    } /*---for i---*/
  } /*--- if (env->all2all) ---*/

}

/*---------------------------------------------------------------------------*/

void gm_compute_combine(GMMetrics* metrics,
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
  GMAssert(j_proc >= 0 && j_proc < env->num_proc);

  switch (env->metric_type) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
    /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
    /*----------------------------------------*/
      
      if (env->num_way == 2) {
        gm_compute_czekanowski_combine(metrics,
                                      metrics_buf,
                                      (GMFloat*)vector_sums_left->data,
                                      (GMFloat*)vector_sums_right->data,
                                      j_proc,
                                      do_compute_triang_only,
                                      env);
      } else {
        gm_compute_czekanowski_3way_combine(metrics,
                                           (GMFloat*)vector_sums_left->data,
                                           (GMFloat*)vector_sums_right->data,
                                           j_proc,
                                           do_compute_triang_only,
                                           env);
      }
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

/*---------------------------------------------------------------------------*/
