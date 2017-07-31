/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#ifndef _MAGMA_tally2BLAS_
#define _MAGMA_tally2BLAS_

typedef int magma_tally2_int_t;

#include <cublas.h>
#include <cuda.h>

#include "magma_tally2blas_z.h"
#include "magma_tally2blas_c.h"
#include "magma_tally2blas_d.h"
#include "magma_tally2blas_s.h"
#include "magma_tally2blas_zc.h"
#include "magma_tally2blas_ds.h"

#if (GPUSHMEM < 200)  
  #define magma_tally2blas_zgemm cublasZgemm
#endif
#define magma_tally2blas_cgemm cublasCgemm

#endif
