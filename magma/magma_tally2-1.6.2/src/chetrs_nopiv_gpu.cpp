/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       @author Adrien REMY

       @generated from zhetrs_nopiv_gpu.cpp normal z -> c, Fri Jan 30 19:00:16 2015

*/
#include "common_magma_tally2.h"

/**
    Purpose
    -------
    Solves a system of linear equations A*X = B with a complex
    Hermitian matrix A using the factorization A = U*D*U**H or
    A = L*D*L**H computed by CHETRF_NOPIV_GPU.
    
    Arguments
    ---------
    
    @param[in]
    uplo    magma_tally2_uplo_t
      -     = Magma_tally2Upper:  Upper triangle of A is stored;
      -     = Magma_tally2Lower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in]
    dA      COMPLEX array on the GPU, dimension (LDA,N)
            The block diagonal matrix D and the multipliers used to
            obtain the factor U or L as computed by CHETRF_NOPIV_GPU.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    param[in,out]
    dB      COMPLEX array on the GPU, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally2_chesv_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_chetrs_nopiv_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info)
{
    magma_tally2FloatComplex c_one = MAGMA_tally2_C_ONE;

    int                upper = (uplo == Magma_tally2Upper);
    *info = 0;
    if (! upper && uplo != Magma_tally2Lower) {
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
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return *info;
    }



  if (upper) {
            magma_tally2blas_ctrsm( Magma_tally2Left, Magma_tally2Upper, 
                           Magma_tally2ConjTrans, Magma_tally2Unit, 
                           n, nrhs, c_one,
                           dA, ldda, dB, lddb );
            magma_tally2blas_clascl_diag(Magma_tally2Upper, 1, n, dA, ldda, dB,1, info);
            magma_tally2blas_ctrsm( Magma_tally2Left, Magma_tally2Upper, 
                           Magma_tally2NoTrans, Magma_tally2Unit, 
                           n, nrhs, c_one,
                           dA, ldda, dB, lddb );
        } else {

            magma_tally2_ctrsm( Magma_tally2Left, Magma_tally2Lower, 
                           Magma_tally2NoTrans, Magma_tally2Unit, 
                           n, nrhs, c_one,
                           dA, ldda, dB, lddb );
            magma_tally2blas_clascl_diag(Magma_tally2Lower, 1, n, dA, ldda, dB,1, info);
            magma_tally2blas_ctrsm( Magma_tally2Left, Magma_tally2Lower, 
                           Magma_tally2ConjTrans, Magma_tally2Unit, 
                           n, nrhs, c_one,
                           dA, ldda, dB, lddb );
        }

    
     return *info;
}
