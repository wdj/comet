/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zlascl2.cu normal z -> c, Fri Jan 30 19:00:09 2015

       @author Theo Mary
*/
#include "common_magma_tally2.h"

#define NB 64


// each thread block does one NB x n block row of A.
// each thread does one row, starting from left edge and moving right.
__global__ void
clascl2_full(int m, int n, const float* D, magma_tally2FloatComplex* A, int lda)
{
    int ind = blockIdx.x * NB + threadIdx.x;

    float mul = D[ind];
    A += ind;
    if (ind < m) {
        for(int j=0; j < n; j++ )
            A[j*lda] *= mul;
    }
}


// each thread block does one NB x n block row of A.
// each thread does one row, starting from left edge and moving right to diagonal.
__global__ void
clascl2_lower(int m, int n, const float* D, magma_tally2FloatComplex* A, int lda)
{
    int ind = blockIdx.x * NB + threadIdx.x;

    int break_d = (ind < n) ? ind : n-1;

    float mul = D[ind];
    A += ind;
    if (ind < m) {
        for(int j=0; j <= break_d; j++ )
            A[j*lda] *= mul;
    }
}


// each thread block does one NB x n block row of A.
// each thread does one row, starting from right edge and moving left to diagonal.
__global__ void
clascl2_upper(int m, int n, const float *D, magma_tally2FloatComplex* A, int lda)
{
    int ind = blockIdx.x * NB + threadIdx.x;

    float mul = D[ind];
    A += ind;
    if (ind < m) {
        for(int j=n-1; j >= ind; j--)
            A[j*lda] *= mul;
    }
}


/**
    Purpose
    -------
    CLASCL2 scales the M by N complex matrix A by the real diagonal matrix dD.
    TYPE specifies that A may be full, upper triangular, lower triangular.

    Arguments
    ---------
    \param[in]
    type    magma_tally2_type_t
            TYPE indices the storage type of the input matrix A.
            = Magma_tally2Full:   full matrix.
            = Magma_tally2Lower:  lower triangular matrix.
            = Magma_tally2Upper:  upper triangular matrix.
            Other formats that LAPACK supports, MAGMA_tally2 does not currently support.

    \param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    \param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    \param[in]
    dD      REAL vector, dimension (M)
            The diagonal matrix containing the scalar factors. Stored as a vector.

    \param[in,out]
    dA      COMPLEX array, dimension (LDDA,N)
            The matrix to be scaled by dD.  See TYPE for the
            storage type.

    \param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

    \param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value.
    
    @param[in]
    queue   magma_tally2_queue_t
            Queue to execute in.

    @ingroup magma_tally2_caux2
    ********************************************************************/
extern "C" void
magma_tally2blas_clascl2_q(
    magma_tally2_type_t type, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_const_ptr dD,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_queue_t queue,
    magma_tally2_int_t *info )
{
    *info = 0;
    if ( type != Magma_tally2Lower && type != Magma_tally2Upper && type != Magma_tally2Full )
        *info = -1;
    else if ( m < 0 )
        *info = -2;
    else if ( n < 0 )
        *info = -3;
    else if ( ldda < max(1,m) )
        *info = -5;
    
    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return;  //info;
    }
    
    dim3 grid( (m + NB - 1)/NB );
    dim3 threads( NB );
    
    if (type == Magma_tally2Lower) {
        clascl2_lower <<< grid, threads, 0, queue >>> (m, n, dD, dA, ldda);
    }
    else if (type == Magma_tally2Upper) {
        clascl2_upper <<< grid, threads, 0, queue >>> (m, n, dD, dA, ldda);
    }
    else if (type == Magma_tally2Full) {
        clascl2_full  <<< grid, threads, 0, queue >>> (m, n, dD, dA, ldda);
    }
}


/**
    @see magma_tally2blas_clascl2_q
    @ingroup magma_tally2_caux2
    ********************************************************************/
extern "C" void
magma_tally2blas_clascl2(
    magma_tally2_type_t type, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_const_ptr dD,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *info )
{
    magma_tally2blas_clascl2_q( type, m, n, dD, dA, ldda, magma_tally2_stream, info );
}
