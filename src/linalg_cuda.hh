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

//-----------------------------------------------------------------------------

void gm_tc_gemm_start(int m, int n, int k,
                      void* dA, int ldda,
                      void* dB, int lddb,
                      void* dC, int lddc,
                      GMEnv* env);

//-----------------------------------------------------------------------------

#endif // _gm_linalg_cuda_hh_

//-----------------------------------------------------------------------------
