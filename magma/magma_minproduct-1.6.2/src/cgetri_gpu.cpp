/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgetri_gpu.cpp normal z -> c, Fri Jan 30 19:00:14 2015

*/
#include "common_magma_minproduct.h"

#define PRECISION_c

// === Define what BLAS to use ============================================
#define magma_minproduct_ctrsm magma_minproductblas_ctrsm
// === End defining what BLAS to use ======================================

/**
    Purpose
    -------
    CGETRI computes the inverse of a matrix using the LU factorization
    computed by CGETRF. This method inverts U and then computes inv(A) by
    solving the system inv(A)*L = inv(U) for inv(A).
    
    Note that it is generally both faster and more accurate to use CGESV,
    or CGETRF and CGETRS, to solve the system AX = B, rather than inverting
    the matrix and multiplying to form X = inv(A)*B. Only in special
    instances should an explicit inverse be computed with this routine.

    Arguments
    ---------
    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    dA      COMPLEX array on the GPU, dimension (LDDA,N)
            On entry, the factors L and U from the factorization
            A = P*L*U as computed by CGETRF_GPU.
            On exit, if INFO = 0, the inverse of the original matrix A.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[in]
    ipiv    INTEGER array, dimension (N)
            The pivot indices from CGETRF; for 1 <= i <= N, row i of the
            matrix was interchanged with row IPIV(i).

    @param[out]
    dwork   (workspace) COMPLEX array on the GPU, dimension (MAX(1,LWORK))
  
    @param[in]
    lwork   INTEGER
            The dimension of the array DWORK.  LWORK >= N*NB, where NB is
            the optimal blocksize returned by magma_minproduct_get_cgetri_nb(n).
    \n
            Unlike LAPACK, this version does not currently support a
            workspace query, because the workspace is on the GPU.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
                  singular and its cannot be computed.

    @ingroup magma_minproduct_cgesv_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_cgetri_gpu(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductFloatComplex_ptr dwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info )
{
    #define dA(i, j)  (dA + (i) + (j)*ldda)
    #define dL(i, j)  (dL + (i) + (j)*lddl)
    
    /* Local variables */
    magma_minproductFloatComplex c_zero    = MAGMA_minproduct_C_ZERO;
    magma_minproductFloatComplex c_one     = MAGMA_minproduct_C_ONE;
    magma_minproductFloatComplex c_neg_one = MAGMA_minproduct_C_NEG_ONE;
    magma_minproductFloatComplex_ptr dL = dwork;
    magma_minproduct_int_t lddl = n;
    magma_minproduct_int_t nb   = magma_minproduct_get_cgetri_nb(n);
    magma_minproduct_int_t j, jmax, jb, jp;
    
    *info = 0;
    if (n < 0)
        *info = -1;
    else if (ldda < max(1,n))
        *info = -3;
    else if ( lwork < n*nb )
        *info = -6;

    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if ( n == 0 )
        return *info;
    
    /* Invert the triangular factor U */
    magma_minproduct_ctrtri_gpu( Magma_minproductUpper, Magma_minproductNonUnit, n, dA, ldda, info );
    if ( *info != 0 )
        return *info;
    
    jmax = ((n-1) / nb)*nb;
    for( j = jmax; j >= 0; j -= nb ) {
        jb = min( nb, n-j );
        
        // copy current block column of A to work space dL
        // (only needs lower trapezoid, but we also copy upper triangle),
        // then zero the strictly lower trapezoid block column of A.
        magma_minproductblas_clacpy( Magma_minproductFull, n-j, jb,
                          dA(j,j), ldda,
                          dL(j,0), lddl );
        magma_minproductblas_claset( Magma_minproductLower, n-j-1, jb, c_zero, c_zero, dA(j+1,j), ldda );
        
        // compute current block column of Ainv
        // Ainv(:, j:j+jb-1)
        //   = ( U(:, j:j+jb-1) - Ainv(:, j+jb:n) L(j+jb:n, j:j+jb-1) )
        //   * L(j:j+jb-1, j:j+jb-1)^{-1}
        // where L(:, j:j+jb-1) is stored in dL.
        if ( j+jb < n ) {
            magma_minproduct_cgemm( Magma_minproductNoTrans, Magma_minproductNoTrans, n, jb, n-j-jb,
                         c_neg_one, dA(0,j+jb), ldda,
                                    dL(j+jb,0), lddl,
                         c_one,     dA(0,j),    ldda );
        }
        // TODO use magma_minproductblas work interface
        magma_minproduct_ctrsm( Magma_minproductRight, Magma_minproductLower, Magma_minproductNoTrans, Magma_minproductUnit,
                     n, jb, c_one,
                     dL(j,0), lddl,
                     dA(0,j), ldda );
    }

    // Apply column interchanges
    for( j = n-2; j >= 0; --j ) {
        jp = ipiv[j] - 1;
        if ( jp != j ) {
            magma_minproductblas_cswap( n, dA(0,j), 1, dA(0,jp), 1 );
        }
    }
    
    return *info;
}
