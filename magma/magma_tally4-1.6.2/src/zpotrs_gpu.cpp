/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    ZPOTRS solves a system of linear equations A*X = B with a Hermitian
    positive definite matrix A using the Cholesky factorization
    A = U**H*U or A = L*L**H computed by ZPOTRF.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally4_uplo_t
      -     = Magma_tally4Upper:  Upper triangle of A is stored;
      -     = Magma_tally4Lower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in]
    dA      COMPLEX_16 array on the GPU, dimension (LDDA,N)
            The triangular factor U or L from the Cholesky factorization
            A = U**H*U or A = L*L**H, as computed by ZPOTRF.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[in,out]
    dB      COMPLEX_16 array on the GPU, dimension (LDDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally4_zposv_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_zpotrs_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info)
{
    magma_tally4DoubleComplex c_one = MAGMA_tally4_Z_ONE;

    *info = 0;
    if ( uplo != Magma_tally4Upper && uplo != Magma_tally4Lower )
        *info = -1;
    if ( n < 0 )
        *info = -2;
    if ( nrhs < 0)
        *info = -3;
    if ( ldda < max(1, n) )
        *info = -5;
    if ( lddb < max(1, n) )
        *info = -7;
    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if ( (n == 0) || (nrhs == 0) ) {
        return *info;
    }

    if ( uplo == Magma_tally4Upper ) {
        if ( nrhs == 1) {
            magma_tally4_ztrsv(Magma_tally4Upper, Magma_tally4ConjTrans, Magma_tally4NonUnit, n, dA, ldda, dB, 1 );
            magma_tally4_ztrsv(Magma_tally4Upper, Magma_tally4NoTrans,   Magma_tally4NonUnit, n, dA, ldda, dB, 1 );
        } else {
            magma_tally4_ztrsm(Magma_tally4Left, Magma_tally4Upper, Magma_tally4ConjTrans, Magma_tally4NonUnit, n, nrhs, c_one, dA, ldda, dB, lddb);
            magma_tally4_ztrsm(Magma_tally4Left, Magma_tally4Upper, Magma_tally4NoTrans,   Magma_tally4NonUnit, n, nrhs, c_one, dA, ldda, dB, lddb);
        }
    }
    else {
        if ( nrhs == 1) {
            magma_tally4_ztrsv(Magma_tally4Lower, Magma_tally4NoTrans,   Magma_tally4NonUnit, n, dA, ldda, dB, 1 );
            magma_tally4_ztrsv(Magma_tally4Lower, Magma_tally4ConjTrans, Magma_tally4NonUnit, n, dA, ldda, dB, 1 );
        } else {
            magma_tally4_ztrsm(Magma_tally4Left, Magma_tally4Lower, Magma_tally4NoTrans,   Magma_tally4NonUnit, n, nrhs, c_one, dA, ldda, dB, lddb);
            magma_tally4_ztrsm(Magma_tally4Left, Magma_tally4Lower, Magma_tally4ConjTrans, Magma_tally4NonUnit, n, nrhs, c_one, dA, ldda, dB, lddb);
        }
    }

    return *info;
}
