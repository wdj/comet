/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_utils_magma.c
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Magma utilities.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include "cuda.h"

#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"

#include "magma_tally4.h"
#include "magma_tally4_lapack.h"

#include "magma_tally3.h"
#include "magma_tally3_lapack.h"

#include "env.h"
#include "vector_sums.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_metrics_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Magma setup, teardown---*/

void gm_magma_initialize(GMEnv* env) {
  GMAssert(env != NULL);

  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  /*----------------------------------------*/
  if (Env_metric_type(env) == GM_METRIC_TYPE_SORENSON) {
  /*----------------------------------------*/

    GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI) {
  /*----------------------------------------*/

    magma_minproduct_int_t magma_code = 0;
    magma_code = magma_code * 1; /*---Avoid unused variable warning---*/

    magma_code = magma_minproduct_init();
    GMAssertAlways(magma_code == MAGMA_minproduct_SUCCESS);
    /*---need this -- see
     * http://on-demand.gputechconf.com/gtc/2014/presentations/S4158-cuda-streams-best-practices-common-pitfalls.pdf
     * page 14 ---*/
    magma_code = magma_minproductblasSetKernelStream(Env_stream_compute(env));
    GMAssertAlways(magma_code == MAGMA_minproduct_SUCCESS);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_num_way(env) == GM_NUM_WAY_2) {
  /*----------------------------------------*/

    magma_tally4_int_t magma_code = 0;
    magma_code = magma_code * 1; /*---Avoid unused variable warning---*/

    magma_code = magma_tally4_init();
    GMAssertAlways(magma_code == MAGMA_tally4_SUCCESS);
    /*---need this -- see
     * http://on-demand.gputechconf.com/gtc/2014/presentations/S4158-cuda-streams-best-practices-common-pitfalls.pdf
     * page 14 ---*/
    magma_code = magma_tally4blasSetKernelStream(Env_stream_compute(env));
    GMAssertAlways(magma_code == MAGMA_tally4_SUCCESS);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_num_way(env) == GM_NUM_WAY_3) {
  /*----------------------------------------*/

      magma_tally3_int_t magma_code = 0;
      magma_code = magma_code * 1; /*---Avoid unused variable warning---*/

      magma_code = magma_tally3_init();
      GMAssertAlways(magma_code == MAGMA_tally3_SUCCESS);
      /*---need this -- see
       * http://on-demand.gputechconf.com/gtc/2014/presentations/S4158-cuda-streams-best-practices-common-pitfalls.pdf
       * page 14 ---*/
      magma_code = magma_tally3blasSetKernelStream(Env_stream_compute(env));
      GMAssertAlways(magma_code == MAGMA_tally3_SUCCESS);

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*---------------------------------------------------------------------------*/

void gm_magma_finalize(GMEnv* env) {
  GMAssert(env != NULL);

  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  /*----------------------------------------*/
  if (Env_metric_type(env) == GM_METRIC_TYPE_SORENSON) {
  /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI) {
  /*----------------------------------------*/

    magma_minproduct_int_t magma_code = 0;
    magma_code = magma_code * 1; /*---Avoid unused variable warning---*/

    // TODO: reset kernel stream (not really needed)
    magma_code = magma_minproduct_finalize();
    GMAssertAlways(magma_code == MAGMA_minproduct_SUCCESS);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_num_way(env) == GM_NUM_WAY_2) {
  /*----------------------------------------*/

    magma_tally4_int_t magma_code = 0;
    magma_code = magma_code * 1; /*---Avoid unused variable warning---*/

    // TODO: reset kernel stream (not really needed)
    magma_code = magma_tally4_finalize();
    GMAssertAlways(magma_code == MAGMA_tally4_SUCCESS);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_num_way(env) == GM_NUM_WAY_3) {
  /*----------------------------------------*/

    magma_tally3_int_t magma_code = 0;
    magma_code = magma_code * 1; /*---Avoid unused variable warning---*/

    // TODO: reset kernel stream (not really needed)
    magma_code = magma_tally3_finalize();
    GMAssertAlways(magma_code == MAGMA_tally3_SUCCESS);

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
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

  /*----------------------------------------*/
  if (Env_metric_type(env) == GM_METRIC_TYPE_SORENSON) {
  /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI) {
  /*----------------------------------------*/

    magma_minproduct_int_t magma_code = 0;
    magma_code = magma_code * 1; /*---Avoid unused variable warning---*/

#ifdef FP_PRECISION_DOUBLE
    magma_code = magma_minproduct_dmalloc_pinned((GMFloat**)&p.h, n);
#endif
#ifdef FP_PRECISION_SINGLE
    magma_code = magma_minproduct_smalloc_pinned((GMFloat**)&p.h, n);
#endif
    GMAssertAlways(magma_code == MAGMA_minproduct_SUCCESS);

#ifdef FP_PRECISION_DOUBLE
    magma_code = magma_minproduct_dmalloc((GMFloat**)&p.d, n);
#endif
#ifdef FP_PRECISION_SINGLE
    magma_code = magma_minproduct_smalloc((GMFloat**)&p.d, n);
#endif
    GMAssertAlways(magma_code == MAGMA_minproduct_SUCCESS);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_num_way(env) == GM_NUM_WAY_2) {
  /*----------------------------------------*/

    typedef magma_tally4DoubleComplex Float_t;

    magma_tally4_int_t magma_code = 0;
    magma_code = magma_code * 1; /*---Avoid unused variable warning---*/

    magma_code = magma_tally4_zmalloc_pinned((Float_t**)&p.h, n);
    GMAssertAlways(magma_code == MAGMA_tally4_SUCCESS);

    magma_code = magma_tally4_zmalloc((Float_t**)&p.d, n);
    GMAssertAlways(magma_code == MAGMA_tally4_SUCCESS);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_num_way(env) == GM_NUM_WAY_3) {
  /*----------------------------------------*/

    typedef magma_tally3DoubleComplex Float_t;

    magma_tally3_int_t magma_code = 0;
    magma_code = magma_code * 1; /*---Avoid unused variable warning---*/

    magma_code = magma_tally3_zmalloc_pinned((Float_t**)&p.h, n);
    GMAssertAlways(magma_code == MAGMA_tally3_SUCCESS);

    magma_code = magma_tally3_zmalloc((Float_t**)&p.d, n);
    GMAssertAlways(magma_code == MAGMA_tally3_SUCCESS);

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/

  return p;
}

/*---------------------------------------------------------------------------*/

void gm_free_magma(GMMirroredPointer* p, GMEnv* env) {
  GMAssert(p != NULL);
  GMAssert(env != NULL);

  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  /*----------------------------------------*/
  if (Env_metric_type(env) == GM_METRIC_TYPE_SORENSON) {
  /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI) {
  /*----------------------------------------*/

    magma_minproduct_int_t magma_code = 0;
    magma_code = magma_code * 1; /*---Avoid unused variable warning---*/

    magma_minproduct_free_pinned(p->h);
    GMAssertAlways(magma_code == MAGMA_minproduct_SUCCESS);
    magma_minproduct_free(p->d);
    GMAssertAlways(magma_code == MAGMA_minproduct_SUCCESS);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_num_way(env) == GM_NUM_WAY_2) {
  /*----------------------------------------*/

    magma_tally4_int_t magma_code = 0;
    magma_code = magma_code * 1; /*---Avoid unused variable warning---*/

    magma_tally4_free_pinned(p->h);
    GMAssertAlways(magma_code == MAGMA_tally4_SUCCESS);
    magma_tally4_free(p->d);
    GMAssertAlways(magma_code == MAGMA_tally4_SUCCESS);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_num_way(env) == GM_NUM_WAY_3) {
  /*----------------------------------------*/

    magma_tally3_int_t magma_code = 0;
    magma_code = magma_code * 1; /*---Avoid unused variable warning---*/

    magma_tally3_free_pinned(p->h);
    GMAssertAlways(magma_code == MAGMA_tally3_SUCCESS);
    magma_tally3_free(p->d);
    GMAssertAlways(magma_code == MAGMA_tally3_SUCCESS);

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*---------------------------------------------------------------------------*/

void gm_magma_set_matrix_zero_start(GMMirroredPointer* matrix_buf,
                                    int mat_dim1,
                                    int mat_dim2,
                                    GMEnv* env) {
  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  /*----------------------------------------*/
  if (Env_metric_type(env) == GM_METRIC_TYPE_SORENSON) {
  /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI) {
  /*----------------------------------------*/

#ifdef FP_PRECISION_DOUBLE
    magma_minproductblas_dlaset
#endif
#ifdef FP_PRECISION_SINGLE
    magma_minproductblas_slaset
#endif
      (Magma_minproductFull, mat_dim1, mat_dim2, (GMFloat)0, (GMFloat)0,
       (GMFloat*)matrix_buf->d, mat_dim1);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_num_way(env) == GM_NUM_WAY_2) {
  /*----------------------------------------*/

    typedef magma_tally4DoubleComplex Float_t;

    Float_t zero = {0, 0};

    magma_tally4blas_zlaset(Magma_tally4Full, mat_dim1, mat_dim2, zero, zero,
                            (Float_t*)matrix_buf->d, mat_dim1);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_num_way(env) == GM_NUM_WAY_3) {
  /*----------------------------------------*/

    typedef magma_tally3DoubleComplex Float_t;

    Float_t zero = {0, 0};

    magma_tally3blas_zlaset(Magma_tally3Full, mat_dim1, mat_dim2, zero, zero,
                            (Float_t*)matrix_buf->d, mat_dim1);

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*---------------------------------------------------------------------------*/

void gm_magma_gemm_start(magma_minproduct_int_t m,
                         magma_minproduct_int_t n,
                         magma_minproduct_int_t k,
                         void* dA,
                         magma_minproduct_int_t ldda,
                         void* dB,
                         magma_minproduct_int_t lddb,
                         void* dC,
                         magma_minproduct_int_t lddc,
                         GMEnv* env) {
  GMAssert(m >= 0);
  GMAssert(n >= 0);
  GMAssert(k >= 0);
  GMAssert(dA != NULL);
  GMAssert(dB != NULL);
  GMAssert(dC != NULL);
  GMAssert(ldda >= 0);
  GMAssert(lddb >= 0);
  GMAssert(env != NULL);

  GMAssert(Env_compute_method(env) == GM_COMPUTE_METHOD_GPU);

  /*----------------------------------------*/
  if (Env_metric_type(env) == GM_METRIC_TYPE_SORENSON) {
  /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI) {
  /*----------------------------------------*/

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

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_num_way(env) == GM_NUM_WAY_2) {
  /*----------------------------------------*/

    typedef magma_tally4DoubleComplex Float_t;

    const Float_t alpha = {1, 0};
    const Float_t beta = {0, 0};

    magma_tally4blas_zgemm(Magma_tally4Trans, Magma_tally4NoTrans, m, n, k,
                           alpha, (Float_t*)dA, ldda, (Float_t*)dB, lddb,
                           beta, (Float_t*)dC, lddc);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_num_way(env) == GM_NUM_WAY_3) {
  /*----------------------------------------*/

    typedef magma_tally3DoubleComplex Float_t;

    const Float_t alpha = {1, 0};
    const Float_t beta = {0, 0};

    magma_tally3blas_zgemm(Magma_tally3Trans, Magma_tally3NoTrans, m, n, k,
                           alpha, (Float_t*)dA, ldda, (Float_t*)dB, lddb,
                           beta, (Float_t*)dC, lddc);

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*---------------------------------------------------------------------------*/
/*---Wait for any computation on the GPU top complete---*/

void gm_compute_wait(GMEnv* env) {
  GMAssert(env != NULL);

  if (Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    cudaStreamSynchronize(Env_stream_compute(env));
    GMAssertAlways(GMEnv_cuda_last_call_succeeded(env));
  }
}

/*===========================================================================*/
/*---Start/end transfer of generic matrix to GPU---*/

void gm_set_matrix_start(GMMirroredPointer* matrix_buf,
                         int mat_dim1,
                         int mat_dim2,
                         GMEnv* env) {
  GMAssert(matrix_buf != NULL);
  GMAssert(env != NULL);

  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  /*---Send vectors to GPU---*/

  /*----------------------------------------*/
  if (Env_metric_type(env) == GM_METRIC_TYPE_SORENSON) {
  /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI) {
  /*----------------------------------------*/

#ifdef FP_PRECISION_DOUBLE
    magma_minproduct_dsetmatrix_async(
        mat_dim1, mat_dim2, (GMFloat*)matrix_buf->h, mat_dim1,
        (GMFloat*)matrix_buf->d, mat_dim1, Env_stream_togpu(env));
#endif
#ifdef FP_PRECISION_SINGLE
    magma_minproduct_ssetmatrix_async(
        mat_dim1, mat_dim2, (GMFloat*)matrix_buf->h, mat_dim1,
        (GMFloat*)matrix_buf->d, mat_dim1, Env_stream_togpu(env));
#endif

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_num_way(env) == GM_NUM_WAY_2) {
  /*----------------------------------------*/

    typedef magma_tally4DoubleComplex Float_t;

    magma_tally4_zsetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->h,
                                  mat_dim1, (Float_t*)matrix_buf->d, mat_dim1,
                                  Env_stream_togpu(env));

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_num_way(env) == GM_NUM_WAY_3) {
  /*----------------------------------------*/

    typedef magma_tally3DoubleComplex Float_t;

    magma_tally3_zsetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->h,
                                  mat_dim1, (Float_t*)matrix_buf->d, mat_dim1,
                                  Env_stream_togpu(env));

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*---------------------------------------------------------------------------*/

void gm_set_matrix_wait(GMEnv* env) {
  GMAssert(env != NULL);

  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  cudaStreamSynchronize(Env_stream_togpu(env));
  GMAssertAlways(GMEnv_cuda_last_call_succeeded(env));
}

/*===========================================================================*/
/*---Start/end transfer of generic matrix from GPU---*/

void gm_get_matrix_start(GMMirroredPointer* matrix_buf,
                         int mat_dim1,
                         int mat_dim2,
                         GMEnv* env) {
  GMAssert(matrix_buf != NULL);
  GMAssert(env != NULL);

  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  /*---Get vectors from GPU---*/

  /*----------------------------------------*/
  if (Env_metric_type(env) == GM_METRIC_TYPE_SORENSON) {
  /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI) {
  /*----------------------------------------*/

#ifdef FP_PRECISION_DOUBLE
    magma_minproduct_dgetmatrix_async(
        mat_dim1, mat_dim2, (GMFloat*)matrix_buf->d, mat_dim1,
        (GMFloat*)matrix_buf->h, mat_dim1, Env_stream_fromgpu(env));
#endif
#ifdef FP_PRECISION_SINGLE
    magma_minproduct_sgetmatrix_async(
        mat_dim1, mat_dim2, (GMFloat*)matrix_buf->d, mat_dim1,
        (GMFloat*)matrix_buf->h, mat_dim1, Env_stream_fromgpu(env));
#endif

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_num_way(env) == GM_NUM_WAY_2) {
  /*----------------------------------------*/

    typedef magma_tally4DoubleComplex Float_t;

    magma_tally4_zgetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->d,
                                  mat_dim1, (Float_t*)matrix_buf->h, mat_dim1,
                                  Env_stream_fromgpu(env));

  /*----------------------------------------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_num_way(env) == GM_NUM_WAY_3) {
  /*----------------------------------------*/

    typedef magma_tally3DoubleComplex Float_t;

    magma_tally3_zgetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->d,
                                  mat_dim1, (Float_t*)matrix_buf->h, mat_dim1,
                                  Env_stream_fromgpu(env));

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*---------------------------------------------------------------------------*/

void gm_get_matrix_wait(GMEnv* env) {
  GMAssert(env != NULL);

  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  cudaStreamSynchronize(Env_stream_fromgpu(env));
  GMAssertAlways(GMEnv_cuda_last_call_succeeded(env));
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
