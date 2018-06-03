//-----------------------------------------------------------------------------
/*!
 * \file   linalg.cc
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Magma interface.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "cuda.h"

#include "magma_minproduct.h"
#include "magma_minproduct_lapack.h"

#include "magma_tally4.h"
#include "magma_tally4_lapack.h"

#include "magma_tally3.h"
#include "magma_tally3_lapack.h"

#include "magma_tally2.h"
#include "magma_tally2_lapack.h"

#include "env.hh"
#include "vector_sums.hh"
#include "vectors.hh"
#include "metrics.hh"

#include "linalg.hh"
#include "linalg_cuda.hh"

//=============================================================================
/*---Magma setup, teardown---*/

void gm_linalg_initialize(GMEnv* env) {
  GMInsist(env);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  /*----------------------------------------*/
  if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK) {
  /*----------------------------------------*/

    magma_minproduct_int_t magma_code = magma_minproduct_init();
    GMInsist(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in call to magma_minproduct_init.");
    /*---need this -- see
     * http://on-demand.gputechconf.com/gtc/2014/presentations/S4158-cuda-streams-best-practices-common-pitfalls.pdf
     * page 14 ---*/
    magma_code = magma_minproductblasSetKernelStream(GMEnv_stream_compute(env));
    GMInsist(magma_code == MAGMA_minproduct_SUCCESS ?
                   "Error in call to magma_minproductblasSetKernelStream." : 0);

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC && env->sparse) {
  /*----------------------------------------*/

    magma_tally2_int_t magma_code = magma_tally2_init();
    GMInsist(magma_code == MAGMA_tally2_SUCCESS &&
                   "Error in call to magma_tally2_init.");
    /*---need this -- see
     * http://on-demand.gputechconf.com/gtc/2014/presentations/S4158-cuda-streams-best-practices-common-pitfalls.pdf
     * page 14 ---*/
    magma_code = magma_tally2blasSetKernelStream(GMEnv_stream_compute(env));
    GMInsist(magma_code == MAGMA_tally2_SUCCESS ?
                   "Error in call to magma_tally2blasSetKernelStream." : 0);

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_num_way(env) == GM_NUM_WAY_2) {
  /*----------------------------------------*/

    magma_tally4_int_t magma_code = magma_tally4_init();
    GMInsist(magma_code == MAGMA_tally4_SUCCESS &&
                   "Error in call to magma_tally4_init.");
    /*---need this -- see
     * http://on-demand.gputechconf.com/gtc/2014/presentations/S4158-cuda-streams-best-practices-common-pitfalls.pdf
     * page 14 ---*/
    magma_code = magma_tally4blasSetKernelStream(GMEnv_stream_compute(env));
    GMInsist(magma_code == MAGMA_tally4_SUCCESS &&
                   "Error in call to magma_tally4blasSetKernelStream.");

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_num_way(env) == GM_NUM_WAY_3) {
  /*----------------------------------------*/

    magma_tally3_int_t magma_code = magma_tally3_init();
    GMInsist(magma_code == MAGMA_tally3_SUCCESS &&
                   "Error in call to magma_tally3_init.");
    /*---need this -- see
     * http://on-demand.gputechconf.com/gtc/2014/presentations/S4158-cuda-streams-best-practices-common-pitfalls.pdf
     * page 14 ---*/
    magma_code = magma_tally3blasSetKernelStream(GMEnv_stream_compute(env));
    GMInsist(magma_code == MAGMA_tally3_SUCCESS &&
                   "Error in call to magma_tally3blasSetKernelStream.");

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

      GMInsistInterface(env, false && "Unimplemented.");

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//-----------------------------------------------------------------------------

void gm_linalg_finalize(GMEnv* env) {
  GMInsist(env);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  /*----------------------------------------*/
  if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK) {
  /*----------------------------------------*/


    // TODO: reset kernel stream (not really needed)
    magma_minproduct_int_t magma_code = magma_minproduct_finalize();
    GMInsist(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in call to magma_minproduct_finalize.");

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC && env->sparse) {
  /*----------------------------------------*/

    // TODO: reset kernel stream (not really needed)
    magma_tally2_int_t magma_code = magma_tally2_finalize();
    GMInsist(magma_code == MAGMA_tally2_SUCCESS &&
                   "Error in call to magma_tally2_finalize.");

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_num_way(env) == GM_NUM_WAY_2) {
  /*----------------------------------------*/

    // TODO: reset kernel stream (not really needed)
    magma_tally4_int_t magma_code = magma_tally4_finalize();
    GMInsist(magma_code == MAGMA_tally4_SUCCESS &&
                   "Error in call to magma_tally4_finalize.");

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_num_way(env) == GM_NUM_WAY_3) {
  /*----------------------------------------*/

    // TODO: reset kernel stream (not really needed)
    magma_tally3_int_t magma_code = magma_tally3_finalize();
    GMInsist(magma_code == MAGMA_tally3_SUCCESS &&
                   "Error in call to magma_tally3_finalize.");

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

      GMInsistInterface(env, false && "Unimplemented.");

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//=============================================================================
/*---Allocate/free host and device memory---*/

void gm_linalg_malloc(GMMirroredBuf* p, size_t dim0, size_t dim1, GMEnv* env) {
  GMInsist(p && env);
  GMInsist(dim0 + 1 >= 1 && dim1 + 1 >= 1);

  *p = GMMirroredBuf_null();

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  p->dim0 = dim0;
  p->dim1 = dim1;

  const size_t n = dim0 * dim1;

  /*----------------------------------------*/
  if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK) {
  /*----------------------------------------*/

    magma_minproduct_int_t magma_code = 0;

    if (GM_FP_PRECISION_DOUBLE) {
      magma_code = magma_minproduct_dmalloc_pinned((double**)&p->h, n);
      GMInsist(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in call to magma_minproduct_dmalloc_pinned.");
    } else {
      magma_code = magma_minproduct_smalloc_pinned((float**)&p->h, n);
      GMInsist(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in call to magma_minproduct_smalloc_pinned.");
    }
    GMFloat_fill_nan((GMFloat*)p->h, n);

    if (GM_FP_PRECISION_DOUBLE) {
      magma_code = magma_minproduct_dmalloc((double**)&p->d, n);
      GMInsist(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in call to magma_minproduct_dmalloc.");
    } else {
      magma_code = magma_minproduct_smalloc((float**)&p->d, n);
      GMInsist(magma_code == MAGMA_minproduct_SUCCESS &&
                   "Error in call to magma_minproduct_smalloc.");
    }

    p->size = n*sizeof(GMFloat);
    env->cpu_mem += p->size;
    env->cpu_mem_max = gm_max_i8(env->cpu_mem_max, env->cpu_mem);
    env->gpu_mem += p->size;
    env->gpu_mem_max = gm_max_i8(env->gpu_mem_max, env->gpu_mem);

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC && env->sparse) {
  /*----------------------------------------*/

    typedef magma_tally2DoubleComplex Float_t;

    magma_tally2_int_t magma_code = 0;

    magma_code = magma_tally2_zmalloc_pinned((Float_t**)&p->h, n);
    GMInsist(magma_code == MAGMA_tally2_SUCCESS &&
                   "Error in call to magma_tally2_zmalloc_pinned.");

    magma_code = magma_tally2_zmalloc((Float_t**)&p->d, n);
    GMInsist(magma_code == MAGMA_tally2_SUCCESS &&
                   "Error in call to magma_tally2_zmalloc.");

    p->size = n*sizeof(Float_t);
    env->cpu_mem += p->size;
    env->cpu_mem_max = gm_max_i8(env->cpu_mem_max, env->cpu_mem);
    env->gpu_mem += p->size;
    env->gpu_mem_max = gm_max_i8(env->gpu_mem_max, env->gpu_mem);

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_num_way(env) == GM_NUM_WAY_2) {
  /*----------------------------------------*/

    typedef magma_tally4DoubleComplex Float_t;

    magma_tally4_int_t magma_code = 0;

    magma_code = magma_tally4_zmalloc_pinned((Float_t**)&p->h, n);
    GMInsist(magma_code == MAGMA_tally4_SUCCESS &&
                   "Error in call to magma_tally4_zmalloc_pinned.");

    magma_code = magma_tally4_zmalloc((Float_t**)&p->d, n);
    GMInsist(magma_code == MAGMA_tally4_SUCCESS &&
                   "Error in call to magma_tally4_zmalloc.");

    p->size = n*sizeof(Float_t);
    env->cpu_mem += p->size;
    env->cpu_mem_max = gm_max_i8(env->cpu_mem_max, env->cpu_mem);
    env->gpu_mem += p->size;
    env->gpu_mem_max = gm_max_i8(env->gpu_mem_max, env->gpu_mem);

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_num_way(env) == GM_NUM_WAY_3) {
  /*----------------------------------------*/

    typedef magma_tally3DoubleComplex Float_t;

    magma_tally3_int_t magma_code = 0;

    magma_code = magma_tally3_zmalloc_pinned((Float_t**)&p->h, n);
    GMInsist(magma_code == MAGMA_tally3_SUCCESS &&
                   "Error in call to magma_tally3_zmalloc_pinned.");

    magma_code = magma_tally3_zmalloc((Float_t**)&p->d, n);
    GMInsist(magma_code == MAGMA_tally3_SUCCESS &&
                   "Error in call to magma_tally3_zmalloc.");

    p->size = n*sizeof(Float_t);
    env->cpu_mem += p->size;
    env->cpu_mem_max = gm_max_i8(env->cpu_mem_max, env->cpu_mem);
    env->gpu_mem += p->size;
    env->gpu_mem_max = gm_max_i8(env->gpu_mem_max, env->gpu_mem);

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

      GMInsistInterface(env, false && "Unimplemented.");

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/

  GMInsist(p->h && "Invalid host pointer created in gm_linalg_malloc.");
  GMInsist(p->d && "Invalid device pointer created in gm_linalg_malloc.");
  p->is_alias = false;
}

//-----------------------------------------------------------------------------

void gm_linalg_free(GMMirroredBuf* p, GMEnv* env) {
  GMInsist(p && env);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  GMInsist(! p->is_alias);

  const size_t size = p->size;

  /*----------------------------------------*/
  if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK) {
  /*----------------------------------------*/

    magma_minproduct_int_t magma_code = magma_minproduct_free_pinned(p->h);
    GMInsist(magma_code == MAGMA_minproduct_SUCCESS);
    magma_code = magma_minproduct_free(p->d);
    GMInsist(magma_code == MAGMA_minproduct_SUCCESS);

    env->cpu_mem -= size;
    env->gpu_mem -= size;

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC && env->sparse) {
  /*----------------------------------------*/

    magma_tally2_int_t magma_code = magma_tally2_free_pinned(p->h);
    GMInsist(magma_code == MAGMA_tally2_SUCCESS);
    magma_code = magma_tally2_free(p->d);
    GMInsist(magma_code == MAGMA_tally2_SUCCESS);

    env->cpu_mem -= size;
    env->gpu_mem -= size;

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_num_way(env) == GM_NUM_WAY_2) {
  /*----------------------------------------*/

    magma_tally4_int_t magma_code = magma_tally4_free_pinned(p->h);
    GMInsist(magma_code == MAGMA_tally4_SUCCESS);
    magma_code = magma_tally4_free(p->d);
    GMInsist(magma_code == MAGMA_tally4_SUCCESS);

    env->cpu_mem -= size;
    env->gpu_mem -= size;

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_num_way(env) == GM_NUM_WAY_3) {
  /*----------------------------------------*/

    magma_tally3_int_t magma_code = magma_tally3_free_pinned(p->h);
    GMInsist(magma_code == MAGMA_tally3_SUCCESS);
    magma_code = magma_tally3_free(p->d);
    GMInsist(magma_code == MAGMA_tally3_SUCCESS);

    env->cpu_mem -= size;
    env->gpu_mem -= size;

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

      GMInsistInterface(env, false && "Unimplemented.");

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//-----------------------------------------------------------------------------

void gm_linalg_set_matrix_zero_start(GMMirroredBuf* matrix_buf,
                                     GMEnv* env) {
  GMInsist(matrix_buf && env);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  const size_t mat_dim1 = matrix_buf->dim0;
  const size_t mat_dim2 = matrix_buf->dim1;

  // ISSUE: these MAGMA routines don't return an error code.

  /*----------------------------------------*/
  if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK) {
  /*----------------------------------------*/

    if (GM_FP_PRECISION_DOUBLE) {
      magma_minproductblas_dlaset
        (Magma_minproductFull, mat_dim1, mat_dim2, (double)0, (double)0,
         (double*)matrix_buf->d, mat_dim1);
    } else {
      magma_minproductblas_slaset
        (Magma_minproductFull, mat_dim1, mat_dim2, (float)0, (float)0,
         (float*)matrix_buf->d, mat_dim1);
    }

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC && env->sparse) {
  /*----------------------------------------*/

    typedef magma_tally2DoubleComplex Float_t;

    Float_t zero = {0, 0};

    magma_tally2blas_zlaset(Magma_tally2Full, mat_dim1, mat_dim2, zero, zero,
                            (Float_t*)matrix_buf->d, mat_dim1);

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_num_way(env) == GM_NUM_WAY_2) {
  /*----------------------------------------*/

    typedef magma_tally4DoubleComplex Float_t;

    Float_t zero = {0, 0};

    magma_tally4blas_zlaset(Magma_tally4Full, mat_dim1, mat_dim2, zero, zero,
                            (Float_t*)matrix_buf->d, mat_dim1);

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_num_way(env) == GM_NUM_WAY_3) {
  /*----------------------------------------*/

    typedef magma_tally3DoubleComplex Float_t;

    Float_t zero = {0, 0};

    magma_tally3blas_zlaset(Magma_tally3Full, mat_dim1, mat_dim2, zero, zero,
                            (Float_t*)matrix_buf->d, mat_dim1);

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

      GMInsistInterface(env, false && "Unimplemented.");

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//-----------------------------------------------------------------------------

void gm_linalg_gemm_block_start(magma_minproduct_int_t m,
                                magma_minproduct_int_t n,
                                magma_minproduct_int_t k,
                                void* dA,
                                magma_minproduct_int_t ldda,
                                void* dB,
                                magma_minproduct_int_t lddb,
                                void* dC,
                                magma_minproduct_int_t lddc,
                                bool is_beta_one,
                                GMEnv* env) {
  GMInsist(dA && dB && dC && env);
  GMInsist(m >= 0 && n >= 0 && k >= 0);
  GMInsist(ldda >= 0 && lddb >= 0);
  GMInsist(GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU);

  {
    int TransA = 1;
    int TransB = 0;

    magma_minproduct_int_t Am = ( ! TransA ? m : k);
    magma_minproduct_int_t An = ( ! TransA ? k : m);
    magma_minproduct_int_t Bm = ( ! TransB ? k : n);
    magma_minproduct_int_t Bn = ( ! TransB ? n : k);
    size_t sizeA = (size_t) ldda * (An - 1) + Am;
    size_t sizeB = (size_t) lddb * (Bn - 1) + Bm;

    size_t CUBLAS_MAX_1DBUF_SIZE = ((1 << 27) - 512);
    GMInsist( ! (sizeA >= CUBLAS_MAX_1DBUF_SIZE ||
                 sizeB >= CUBLAS_MAX_1DBUF_SIZE ));
  }

  // ISSUE: these MAGMA routines don't return an error code.

  /*----------------------------------------*/
  if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK) {
  /*----------------------------------------*/

    const GMFloat alpha = 1;
    const GMFloat beta = is_beta_one ? 1 : 0;

    if (GM_FP_PRECISION_DOUBLE) {
      magma_minproductblas_dgemm
        (Magma_minproductTrans, Magma_minproductNoTrans, m, n, k, alpha,
         (double*)dA, ldda, (double*)dB, lddb, beta, (double*)dC, lddc);
    } else {
      magma_minproductblas_sgemm
        (Magma_minproductTrans, Magma_minproductNoTrans, m, n, k, alpha,
         (float*)dA, ldda, (float*)dB, lddb, beta, (float*)dC, lddc);
    }
    GMInsist(GMEnv_cuda_last_call_succeeded(env));

    env->ops_local += 2 * m * (double)n * (double)k;

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC && env->sparse) {
  /*----------------------------------------*/

    typedef magma_tally2DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const Float_t alpha = {1, 0};
    const Float_t beta = is_beta_one ? one : zero;

    magma_tally2blas_zgemm(Magma_tally2Trans, Magma_tally2NoTrans, m, n, k,
                           alpha, (Float_t*)dA, ldda, (Float_t*)dB, lddb,
                           beta, (Float_t*)dC, lddc);
    GMInsist(GMEnv_cuda_last_call_succeeded(env));

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_num_way(env) == GM_NUM_WAY_2) {
  /*----------------------------------------*/

    typedef magma_tally4DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const Float_t alpha = {1, 0};
    const Float_t beta = is_beta_one ? one : zero;

    magma_tally4blas_zgemm(Magma_tally4Trans, Magma_tally4NoTrans, m, n, k,
                           alpha, (Float_t*)dA, ldda, (Float_t*)dB, lddb,
                           beta, (Float_t*)dC, lddc);
    GMInsist(GMEnv_cuda_last_call_succeeded(env));

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_num_way(env) == GM_NUM_WAY_3) {
  /*----------------------------------------*/

    typedef magma_tally3DoubleComplex Float_t;

    const Float_t zero = {0, 0};
    const Float_t one = {1, 0};
    const Float_t alpha = {1, 0};
    const Float_t beta = is_beta_one ? one : zero;

    magma_tally3blas_zgemm(Magma_tally3Trans, Magma_tally3NoTrans, m, n, k,
                           alpha, (Float_t*)dA, ldda, (Float_t*)dB, lddb,
                           beta, (Float_t*)dC, lddc);
    GMInsist(GMEnv_cuda_last_call_succeeded(env));

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

      GMInsistInterface(env, false && "Unimplemented.");

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//-----------------------------------------------------------------------------

void gm_linalg_gemm_start(magma_minproduct_int_t m,
                          magma_minproduct_int_t n,
                          magma_minproduct_int_t k,
                          void* dA,
                          magma_minproduct_int_t ldda,
                          void* dB,
                          magma_minproduct_int_t lddb,
                          void* dC,
                          magma_minproduct_int_t lddc,
                          GMEnv* env) {
  GMInsist(dA && dB && dC && env);
  GMInsist(m >= 0 && n >= 0 && k >= 0);
  GMInsist(ldda >= 0 && lddb >= 0);
  GMInsist(GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU);

  if (m==0 || n==0 || k==0) {
    return;
  }

  if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC && env->tc) {
    gm_tc_gemm_start(m, n, k, dA, ldda, dB, lddb, dC, lddc, env);
#if 0
    // NOTE: may be possible to use (smaller) m somehow here.
    const int I_max = m;
    const int nvl = n;
    gm_tc_buf_write(0, I_max, nvl, k, dA, env);
    //for (int i=0; i<2; ++i)
    gm_tc_buf_write(1, I_max, nvl, k, dB, env);
    //for (int i=0; i<10; ++i) //FIX
    gm_tc_solve(I_max, nvl, k, dA, ldda, dB, lddb, dC, lddc, env);
    gm_tc_fix_metrics(I_max, nvl, dC, env);
#endif
    return;
  }

  const size_t rows = k;
  const size_t cols_A = m;
  const size_t cols_B = n;

  const size_t elt_size =
    GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK ? sizeof(GMFloat) :
   (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC && env->sparse) ?
                                         sizeof(magma_tally2DoubleComplex) :
   (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
    GMEnv_num_way(env) == GM_NUM_WAY_2) ? sizeof(magma_tally4DoubleComplex) :
   (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
    GMEnv_num_way(env) == GM_NUM_WAY_3) ? sizeof(magma_tally3DoubleComplex) : 0;
  GMInsist(elt_size > 0);

  const size_t align_factor = 128 / elt_size;
  const size_t max_elts = (1 << 27) - 512;

  /*---TODO: can we improve aspect ratios of submatrices---*/
//  const size_t max_rows_per_block_raw = (1 << 14);
//  const size_t max_cols_per_block_raw = max_elts / max_rows_per_block_raw;

  const size_t max_rows_per_block_raw = rows + align_factor;
  const size_t max_cols_per_block_raw = max_elts / rows;

  const size_t max_rows_per_block = (max_rows_per_block_raw / align_factor)
                                                            * align_factor;
  const size_t max_cols_per_block = (max_cols_per_block_raw / align_factor)
                                                            * align_factor;

  GMInsist(max_rows_per_block != 0);
  GMInsist(max_cols_per_block != 0);

  const size_t cols_per_block_A = gm_min_i8(cols_A, max_cols_per_block);
  const size_t cols_per_block_B = gm_min_i8(cols_B, max_cols_per_block);

  const size_t rows_per_block = gm_min_i8(rows, max_rows_per_block);

  for (size_t row_base=0; row_base<rows; row_base+=rows_per_block) {
    const size_t rows_remaining = rows - row_base;
    const size_t rows_this = gm_min_i8(rows_remaining, rows_per_block);

    for (size_t col_A_base=0; col_A_base<cols_A; col_A_base+=cols_per_block_A) {
      const size_t cols_A_remaining = cols_A - col_A_base;
      const size_t cols_A_this = gm_min_i8(cols_A_remaining, cols_per_block_A);

      void* dA_this = (char*)dA + (row_base + ldda*col_A_base)*elt_size;

      for (size_t col_B_base=0; col_B_base<cols_B;
           col_B_base+=cols_per_block_B) {

        const size_t cols_B_remaining = cols_B - col_B_base;
        const size_t cols_B_this = gm_min_i8(cols_B_remaining,
                                             cols_per_block_B);

        void* dB_this = (char*)dB + (row_base + ldda*col_B_base)*elt_size;

        void* dC_this = (char*)dC + (col_A_base + lddc*col_B_base)*elt_size;

        gm_linalg_gemm_block_start(cols_A_this, cols_B_this, rows_this,
          dA_this, ldda, dB_this, lddb, dC_this, lddc, row_base > 0,  env);
      }
    }
  }
}

//-----------------------------------------------------------------------------
/*---Wait for any computation on the GPU to complete---*/

void gm_compute_wait(GMEnv* env) {
  GMInsist(env);

  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    cudaStreamSynchronize(GMEnv_stream_compute(env));
    GMInsist(GMEnv_cuda_last_call_succeeded(env));
  }
}

//=============================================================================
/*---Start/end transfer of generic matrix to GPU---*/

void gm_linalg_set_matrix_start(GMMirroredBuf* matrix_buf, GMEnv* env) {
  GMInsist(matrix_buf && env);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  const size_t mat_dim1 = matrix_buf->dim0;
  const size_t mat_dim2 = matrix_buf->dim1;

  /*---Send vectors to GPU---*/

  // ISSUE: these MAGMA routines don't return an error code.

  /*----------------------------------------*/
  if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK) {
  /*----------------------------------------*/

    if (GM_FP_PRECISION_DOUBLE) {
      magma_minproduct_dsetmatrix_async(
        mat_dim1, mat_dim2, (double*)matrix_buf->h, mat_dim1,
        (double*)matrix_buf->d, mat_dim1, GMEnv_stream_togpu(env));
    } else {
      magma_minproduct_ssetmatrix_async(
        mat_dim1, mat_dim2, (float*)matrix_buf->h, mat_dim1,
        (float*)matrix_buf->d, mat_dim1, GMEnv_stream_togpu(env));
    }

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC && env->sparse) {
  /*----------------------------------------*/

    typedef magma_tally2DoubleComplex Float_t;

    magma_tally2_zsetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->h,
                                  mat_dim1, (Float_t*)matrix_buf->d, mat_dim1,
                                  GMEnv_stream_togpu(env));

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_num_way(env) == GM_NUM_WAY_2) {
  /*----------------------------------------*/

    typedef magma_tally4DoubleComplex Float_t;

    magma_tally4_zsetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->h,
                                  mat_dim1, (Float_t*)matrix_buf->d, mat_dim1,
                                  GMEnv_stream_togpu(env));

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_num_way(env) == GM_NUM_WAY_3) {
  /*----------------------------------------*/

    typedef magma_tally3DoubleComplex Float_t;

    magma_tally3_zsetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->h,
                                  mat_dim1, (Float_t*)matrix_buf->d, mat_dim1,
                                  GMEnv_stream_togpu(env));

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

      GMInsistInterface(env, false && "Unimplemented.");

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//-----------------------------------------------------------------------------

void gm_linalg_set_matrix_wait(GMEnv* env) {
  GMInsist(env);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  cudaStreamSynchronize(GMEnv_stream_togpu(env));
  GMInsist(GMEnv_cuda_last_call_succeeded(env));
}

//=============================================================================
/*---Start/end transfer of generic matrix from GPU---*/

void gm_linalg_get_matrix_start(GMMirroredBuf* matrix_buf,
                                GMEnv* env) {
  GMInsist(matrix_buf && env);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  const size_t mat_dim1 = matrix_buf->dim0;
  const size_t mat_dim2 = matrix_buf->dim1;

  /*---Get vectors from GPU---*/

  // ISSUE: these MAGMA routines don't return an error code.

  /*----------------------------------------*/
  if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK) {
  /*----------------------------------------*/

    if (GM_FP_PRECISION_DOUBLE) {
      magma_minproduct_dgetmatrix_async(
        mat_dim1, mat_dim2, (double*)matrix_buf->d, mat_dim1,
        (double*)matrix_buf->h, mat_dim1, GMEnv_stream_fromgpu(env));
    } else {
      magma_minproduct_sgetmatrix_async(
        mat_dim1, mat_dim2, (float*)matrix_buf->d, mat_dim1,
        (float*)matrix_buf->h, mat_dim1, GMEnv_stream_fromgpu(env));
    }

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC && env->sparse) {
  /*----------------------------------------*/

    typedef magma_tally2DoubleComplex Float_t;

    magma_tally2_zgetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->d,
                                  mat_dim1, (Float_t*)matrix_buf->h, mat_dim1,
                                  GMEnv_stream_fromgpu(env));

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_num_way(env) == GM_NUM_WAY_2) {
  /*----------------------------------------*/

    typedef magma_tally4DoubleComplex Float_t;

    magma_tally4_zgetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->d,
                                  mat_dim1, (Float_t*)matrix_buf->h, mat_dim1,
                                  GMEnv_stream_fromgpu(env));

  /*----------------------------------------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_num_way(env) == GM_NUM_WAY_3) {
  /*----------------------------------------*/

    typedef magma_tally3DoubleComplex Float_t;

    magma_tally3_zgetmatrix_async(mat_dim1, mat_dim2, (Float_t*)matrix_buf->d,
                                  mat_dim1, (Float_t*)matrix_buf->h, mat_dim1,
                                  GMEnv_stream_fromgpu(env));

  /*----------------------------------------*/
  } else {
  /*----------------------------------------*/

      GMInsistInterface(env, false && "Unimplemented.");

  /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//-----------------------------------------------------------------------------

void gm_linalg_get_matrix_wait(GMEnv* env) {
  GMInsist(env);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  cudaStreamSynchronize(GMEnv_stream_fromgpu(env));
  GMInsist(GMEnv_cuda_last_call_succeeded(env));
}

//-----------------------------------------------------------------------------
