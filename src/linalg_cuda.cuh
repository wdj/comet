//-----------------------------------------------------------------------------
/*!
 * \file   linalg_cuda.cuh
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  Supporting CUDA functions, headers.
 * \note   Copyright (C) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_linalg_cuda_cuh_
#define _gm_linalg_cuda_cuh_

#include "env.hh"

//-----------------------------------------------------------------------------

void gm_tc_gemm_start(int m, int n, int k,
                      void* dA, int ldda,
                      void* dB, int lddb,
                      void* dC, int lddc,
                      GMEnv* env);

//-----------------------------------------------------------------------------

#endif // _gm_linalg_cuda_cuh_

//-----------------------------------------------------------------------------
