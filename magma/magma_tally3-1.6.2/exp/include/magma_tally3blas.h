/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#ifndef _MAGMA_tally3BLAS_
#define _MAGMA_tally3BLAS_

typedef int magma_tally3_int_t;

#include <cublas.h>
#include <cuda.h>

#include "magma_tally3blas_z.h"
#include "magma_tally3blas_c.h"
#include "magma_tally3blas_d.h"
#include "magma_tally3blas_s.h"
#include "magma_tally3blas_zc.h"
#include "magma_tally3blas_ds.h"

#if (GPUSHMEM < 200)  
  #define magma_tally3blas_zgemm cublasZgemm
#endif
#define magma_tally3blas_cgemm cublasCgemm

#endif
