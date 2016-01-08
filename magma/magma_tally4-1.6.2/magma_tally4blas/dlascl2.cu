/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zlascl2.cu normal z -> d, Fri Jan 30 19:00:09 2015

       @author Theo Mary
*/
#include "common_magma_tally4.h"

#define NB 64


// each thread block does one NB x n block row of A.
// each thread does one row, starting from left edge and moving right.
__global__ void
dlascl2_full(int m, int n, const double* D, double* A, int lda)
{
    int ind = blockIdx.x * NB + threadIdx.x;

    double mul = D[ind];
    A += ind;
    if (ind < m) {
        for(int j=0; j < n; j++ )
            A[j*lda] *= mul;
    }
}


// each thread block does one NB x n block row of A.
// each thread does one row, starting from left edge and moving right to diagonal.
__global__ void
dlascl2_lower(int m, int n, const double* D, double* A, int lda)
{
    int ind = blockIdx.x * NB + threadIdx.x;

    int break_d = (ind < n) ? ind : n-1;

    double mul = D[ind];
    A += ind;
    if (ind < m) {
        for(int j=0; j <= break_d; j++ )
            A[j*lda] *= mul;
    }
}


// each thread block does one NB x n block row of A.
// each thread does one row, starting from right edge and moving left to diagonal.
__global__ void
dlascl2_upper(int m, int n, const double *D, double* A, int lda)
{
    int ind = blockIdx.x * NB + threadIdx.x;

    double mul = D[ind];
    A += ind;
    if (ind < m) {
        for(int j=n-1; j >= ind; j--)
            A[j*lda] *= mul;
    }
}


/**
    Purpose
    -------
    DLASCL2 scales the M by N real matrix A by the real diagonal matrix dD.
    TYPE specifies that A may be full, upper triangular, lower triangular.

    Arguments
    ---------
    \param[in]
    type    magma_tally4_type_t
            TYPE indices the storage type of the input matrix A.
            = Magma_tally4Full:   full matrix.
            = Magma_tally4Lower:  lower triangular matrix.
            = Magma_tally4Upper:  upper triangular matrix.
            Other formats that LAPACK supports, MAGMA_tally4 does not currently support.

    \param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    \param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    \param[in]
    dD      DOUBLE PRECISION vector, dimension (M)
            The diagonal matrix containing the scalar factors. Stored as a vector.

    \param[in,out]
    dA      DOUBLE PRECISION array, dimension (LDDA,N)
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
    queue   magma_tally4_queue_t
            Queue to execute in.

    @ingroup magma_tally4_daux2
    ********************************************************************/
extern "C" void
magma_tally4blas_dlascl2_q(
    magma_tally4_type_t type, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_const_ptr dD,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_queue_t queue,
    magma_tally4_int_t *info )
{
    *info = 0;
    if ( type != Magma_tally4Lower && type != Magma_tally4Upper && type != Magma_tally4Full )
        *info = -1;
    else if ( m < 0 )
        *info = -2;
    else if ( n < 0 )
        *info = -3;
    else if ( ldda < max(1,m) )
        *info = -5;
    
    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return;  //info;
    }
    
    dim3 grid( (m + NB - 1)/NB );
    dim3 threads( NB );
    
    if (type == Magma_tally4Lower) {
        dlascl2_lower <<< grid, threads, 0, queue >>> (m, n, dD, dA, ldda);
    }
    else if (type == Magma_tally4Upper) {
        dlascl2_upper <<< grid, threads, 0, queue >>> (m, n, dD, dA, ldda);
    }
    else if (type == Magma_tally4Full) {
        dlascl2_full  <<< grid, threads, 0, queue >>> (m, n, dD, dA, ldda);
    }
}


/**
    @see magma_tally4blas_dlascl2_q
    @ingroup magma_tally4_daux2
    ********************************************************************/
extern "C" void
magma_tally4blas_dlascl2(
    magma_tally4_type_t type, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_const_ptr dD,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *info )
{
    magma_tally4blas_dlascl2_q( type, m, n, dD, dA, ldda, magma_tally4_stream, info );
}
