//-----------------------------------------------------------------------------
/*!
 * \file   linalg_cuda.hh
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  Supporting CUDA functions, header..
 * \note   Copyright (C) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_linalg_cuda_hh_
#define _gm_linalg_cuda_hh_

#include "env.hh"
#include "vectors.hh"

//-----------------------------------------------------------------------------

#if 0
void gm_tc_buf_write(
  int left_right,
  int I_max,
  int num_vector_local,
  int num_packedval_field_local,
  void* bufd, GMEnv* env);

void gm_tc_solve(
  int I_max,
  int num_vector_local,
  int num_packedval_field_local,
  void* dA,
  int ldda,
  void* dB,
  int lddb,
  void* dC,
  int lddc,
  GMEnv* env);

void gm_tc_fix_metrics(
  int I_max,
  int num_vector_local,
  void* bufd,
  GMEnv* env);
#endif

void gm_tc_gemm_start(int m, int n, int k,
                      void* dA, int ldda,
                      void* dB, int lddb,
                      void* dC, int lddc,
                      GMEnv* env);

//-----------------------------------------------------------------------------

#endif // _gm_linalg_cuda_hh_

//-----------------------------------------------------------------------------
