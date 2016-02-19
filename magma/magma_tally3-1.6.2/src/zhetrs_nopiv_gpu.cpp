/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       @author Adrien REMY

       @precisions normal z -> s d c

*/
#include "common_magma_tally3.h"

/**
    Purpose
    -------
    Solves a system of linear equations A*X = B with a complex
    Hermitian matrix A using the factorization A = U*D*U**H or
    A = L*D*L**H computed by ZHETRF_NOPIV_GPU.
    
    Arguments
    ---------
    
    @param[in]
    uplo    magma_tally3_uplo_t
      -     = Magma_tally3Upper:  Upper triangle of A is stored;
      -     = Magma_tally3Lower:  Lower triangle of A is stored.

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

    @ingroup magma_tally3_zhesv_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_zhetrs_nopiv_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info)
{
    magma_tally3DoubleComplex c_one = MAGMA_tally3_Z_ONE;

    int                upper = (uplo == Magma_tally3Upper);
    *info = 0;
    if (! upper && uplo != Magma_tally3Lower) {
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
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return *info;
    }



  if (upper) {
            magma_tally3blas_ztrsm( Magma_tally3Left, Magma_tally3Upper, 
                           Magma_tally3ConjTrans, Magma_tally3Unit, 
                           n, nrhs, c_one,
                           dA, ldda, dB, lddb );
            magma_tally3blas_zlascl_diag(Magma_tally3Upper, 1, n, dA, ldda, dB,1, info);
            magma_tally3blas_ztrsm( Magma_tally3Left, Magma_tally3Upper, 
                           Magma_tally3NoTrans, Magma_tally3Unit, 
                           n, nrhs, c_one,
                           dA, ldda, dB, lddb );
        } else {

            magma_tally3_ztrsm( Magma_tally3Left, Magma_tally3Lower, 
                           Magma_tally3NoTrans, Magma_tally3Unit, 
                           n, nrhs, c_one,
                           dA, ldda, dB, lddb );
            magma_tally3blas_zlascl_diag(Magma_tally3Lower, 1, n, dA, ldda, dB,1, info);
            magma_tally3blas_ztrsm( Magma_tally3Left, Magma_tally3Lower, 
                           Magma_tally3ConjTrans, Magma_tally3Unit, 
                           n, nrhs, c_one,
                           dA, ldda, dB, lddb );
        }

    
     return *info;
}
