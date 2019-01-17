//-----------------------------------------------------------------------------
/*!
 * \file   linalg_cuda.cuh
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  CUDA code to support linear algebra operations, headers.
 * \note   Copyright (C) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_linalg_cuda_cuh_
#define _gm_linalg_cuda_cuh_

#include "cublas_v2.h"

#include "env.hh"

//-----------------------------------------------------------------------------

struct TCBufs {
  size_t tc_buf_size;
  void* tc_buf_left;
  void* tc_buf_right;
  cublasHandle_t cublas_handle;
};

//-----------------------------------------------------------------------------

void gm_tc_gemm_start(int m, int n, int k,
                      void* dA, int ldda,
                      void* dB, int lddb,
                      void* dC, int lddc,
                      TCBufs& tc_bufs,
                      GMEnv* env);

void gm_tc_bufs_malloc(int num_vector_local,
                       int num_field_local,
                       int num_packedval_field_local,
                       TCBufs& tc_bufs,
                       GMEnv* env);

void gm_tc_bufs_free(TCBufs& tc_bufs,
                     GMEnv* env);

//-----------------------------------------------------------------------------

#endif // _gm_linalg_cuda_cuh_

//-----------------------------------------------------------------------------
