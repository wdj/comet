/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       @author Adrien REMY

       @precisions normal z -> s d c

*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    Solves a system of linear equations A*X = B with a complex
    Hermitian matrix A using the factorization A = U*D*U**H or
    A = L*D*L**H computed by ZHETRF_NOPIV_GPU.
    
    Arguments
    ---------
    
    @param[in]
    uplo    magma_minproduct_uplo_t
      -     = Magma_minproductUpper:  Upper triangle of A is stored;
      -     = Magma_minproductLower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in]
    dA      COMPLEX_16 array on the GPU, dimension (LDA,N)
            The block diagonal matrix D and the multipliers used to
            obtain the factor U or L as computed by ZHETRF_NOPIV_GPU.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    param[in,out]
    dB      COMPLEX_16 array on the GPU, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_minproduct_zhesv_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_zhetrs_nopiv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info)
{
    magma_minproductDoubleComplex c_one = MAGMA_minproduct_Z_ONE;

    int                upper = (uplo == Magma_minproductUpper);
    *info = 0;
    if (! upper && uplo != Magma_minproductLower) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (nrhs < 0) {
        *info = -3;
    } else if (ldda < max(1,n)) {
        *info = -5;
    } else if (lddb < max(1,n)) {
        *info = -7;
    }
    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return *info;
    }



  if (upper) {
            magma_minproductblas_ztrsm( Magma_minproductLeft, Magma_minproductUpper, 
                           Magma_minproductConjTrans, Magma_minproductUnit, 
                           n, nrhs, c_one,
                           dA, ldda, dB, lddb );
            magma_minproductblas_zlascl_diag(Magma_minproductUpper, 1, n, dA, ldda, dB,1, info);
            magma_minproductblas_ztrsm( Magma_minproductLeft, Magma_minproductUpper, 
                           Magma_minproductNoTrans, Magma_minproductUnit, 
                           n, nrhs, c_one,
                           dA, ldda, dB, lddb );
        } else {

            magma_minproduct_ztrsm( Magma_minproductLeft, Magma_minproductLower, 
                           Magma_minproductNoTrans, Magma_minproductUnit, 
                           n, nrhs, c_one,
                           dA, ldda, dB, lddb );
            magma_minproductblas_zlascl_diag(Magma_minproductLower, 1, n, dA, ldda, dB,1, info);
            magma_minproductblas_ztrsm( Magma_minproductLeft, Magma_minproductLower, 
                           Magma_minproductConjTrans, Magma_minproductUnit, 
                           n, nrhs, c_one,
                           dA, ldda, dB, lddb );
        }

    
     return *info;
}
