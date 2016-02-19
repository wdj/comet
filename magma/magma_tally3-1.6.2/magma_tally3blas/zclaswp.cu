/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions mixed zc -> ds

*/
#include "common_magma_tally3.h"

#define NB 64

// TODO check precision, as in zlag2c?

__global__ void
zclaswp_kernel(int n, magma_tally3DoubleComplex *A, int lda, magma_tally3FloatComplex *SA, int m, const magma_tally3_int_t *ipiv)
{
    int ind = blockIdx.x*NB + threadIdx.x;
    int newind;
    magma_tally3FloatComplex res;
    
    if (ind < m) {
        SA   += ind;
        ipiv += ind;
        
        newind = ipiv[0];
        
        for(int i=0; i < n; i++) {
            res = MAGMA_tally3_C_MAKE( (float)cuCreal(A[newind+i*lda]),
                                (float)cuCimag(A[newind+i*lda]) );
            SA[i*lda] = res; 
        }
    }
}

__global__ void
zclaswp_inv_kernel(int n, magma_tally3DoubleComplex *A, int lda, magma_tally3FloatComplex *SA, int m, const magma_tally3_int_t *ipiv)
{
    int ind = blockIdx.x*NB + threadIdx.x;
    int newind;
    magma_tally3DoubleComplex res;

    if (ind < m) {
        A    += ind;
        ipiv += ind;

        newind = ipiv[0];

        for(int i=0; i < n; i++) {
            res = MAGMA_tally3_Z_MAKE( (double)cuCrealf(SA[newind+i*lda]),
                                (double)cuCimagf(SA[newind+i*lda]) );
            A[i*lda] = res;
        }
    }
}


/**
    Purpose
    -------
    Row i of  A is cast to single precision in row ipiv[i] of SA (incx > 0), or
    row i of SA is cast to double precision in row ipiv[i] of  A (incx < 0),
    for 0 <= i < M.

    @param[in]
    n       INTEGER.
            On entry, N specifies the number of columns of the matrix A.

    @param[in,out]
    A       DOUBLE PRECISION array on the GPU, dimension (LDA,N)
            On entry, the M-by-N matrix to which the row interchanges will be applied.
            TODO update docs

    @param[in]
    lda     INTEGER.
            LDA specifies the leading dimension of A.

    @param[in,out]
    SA      REAL array on the GPU, dimension (LDA,N)
            On exit, the single precision, permuted matrix.
            TODO update docs
        
    @param[in]
    m       The number of rows to be interchanged.

    @param[in]
    ipiv    INTEGER array on the GPU, dimension (M)
            The vector of pivot indices. Row i of A is cast to single 
            precision in row ipiv[i] of SA, for 0 <= i < m. 

    @param[in]
    incx    INTEGER
            If INCX is negative, the pivots are applied in reverse order,
            otherwise in straight-forward order.
    
    @param[in]
    queue   magma_tally3_queue_t
            Queue to execute in.

    @ingroup magma_tally3_zaux2
    ********************************************************************/
extern "C" void
magma_tally3blas_zclaswp_q(
    magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr A, magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr SA, magma_tally3_int_t m,
    const magma_tally3_int_t *ipiv, magma_tally3_int_t incx,
    magma_tally3_queue_t queue )
{
    int blocks = (m - 1)/NB + 1;
    dim3 grid(blocks, 1, 1);
    dim3 threads(NB, 1, 1);

    if (incx >= 0)
        zclaswp_kernel<<< grid, threads, 0, queue >>>(n, A, lda, SA, m, ipiv);
    else
        zclaswp_inv_kernel<<< grid, threads, 0, queue >>>(n, A, lda, SA, m, ipiv);
}


/**
    @see magma_tally3blas_zclaswp_q
    @ingroup magma_tally3_zaux2
    ********************************************************************/
extern "C" void
magma_tally3blas_zclaswp(
    magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr A, magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr SA, magma_tally3_int_t m,
    const magma_tally3_int_t *ipiv, magma_tally3_int_t incx )
{
    magma_tally3blas_zclaswp_q( n, A, lda, SA, m, ipiv, incx, magma_tally3_stream );
}
