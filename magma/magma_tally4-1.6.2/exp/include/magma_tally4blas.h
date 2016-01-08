/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#ifndef _MAGMA_tally4BLAS_
#define _MAGMA_tally4BLAS_

typedef int magma_tally4_int_t;

#include <cublas.h>
#include <cuda.h>

#include "magma_tally4blas_z.h"
#include "magma_tally4blas_c.h"
#include "magma_tally4blas_d.h"
#include "magma_tally4blas_s.h"
#include "magma_tally4blas_zc.h"
#include "magma_tally4blas_ds.h"

#if (GPUSHMEM < 200)  
  #define magma_tally4blas_zgemm cublasZgemm
#endif
#define magma_tally4blas_cgemm cublasCgemm

#endif
