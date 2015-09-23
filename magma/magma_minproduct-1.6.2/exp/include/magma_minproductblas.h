/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#ifndef _MAGMA_minproductBLAS_
#define _MAGMA_minproductBLAS_

typedef int magma_minproduct_int_t;

#include <cublas.h>
#include <cuda.h>

#include "magma_minproductblas_z.h"
#include "magma_minproductblas_c.h"
#include "magma_minproductblas_d.h"
#include "magma_minproductblas_s.h"
#include "magma_minproductblas_zc.h"
#include "magma_minproductblas_ds.h"

#if (GPUSHMEM < 200)  
  #define magma_minproductblas_zgemm cublasZgemm
#endif
#define magma_minproductblas_cgemm cublasCgemm

#endif
