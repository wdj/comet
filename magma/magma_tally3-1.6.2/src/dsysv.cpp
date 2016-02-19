/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zhesv.cpp normal z -> d, Fri Jan 30 19:00:16 2015
*/

#include "common_magma_tally3.h"

/**
    Purpose
    =======

    DSYSV computes the solution to a real system of linear equations
       A * X = B,
    where A is an n-by-n symmetric matrix and X and B are n-by-nrhs
    matrices.

    The diagonal pivoting method is used to factor A as
       A = U * D * U**H,  if uplo = 'U', or
       A = L * D * L**H,  if uplo = 'L',
    where U (or L) is a product of permutation and unit upper (lower)
    triangular matrices, and D is symmetric and block diagonal with
    1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
    used to solve the system of equations A * X = B.

    Arguments
    =========
    @param[in]
    uplo    CHARACTER*1
            = 'U':  Upper triangle of A is stored;
            = 'L':  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The number of linear equations, i.e., the order of the
            matrix A.  n >= 0.
 
    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  nrhs >= 0.

    @param[in,out]
    A       DOUBLE PRECISION array, dimension (lda,n)
            On entry, the symmetric matrix A.  If uplo = 'U', the leading
            n-by-n upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If uplo = 'L', the
            leading n-by-n lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.

            On exit, if info = 0, the block diagonal matrix D and the
            multipliers used to obtain the factor U or L from the
            factorization A = U*D*U**H or A = L*D*L**H as computed by
            DSYTRF.
 
    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  lda >= max(1,n).
 
    @param[out]
    ipiv    INTEGER array, dimension (n)
            Details of the interchanges and the block structure of D, as
            determined by DSYTRF.  If ipiv(k) > 0, then rows and columns
            k and ipiv(k) were interchanged, and D(k,k) is a 1-by-1
            diagonal block.  If uplo = 'U' and ipiv(k) = ipiv(k-1) < 0,
            then rows and columns k-1 and -ipiv(k) were interchanged and
            D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If uplo = 'L' and
            ipiv(k) = ipiv(k+1) < 0, then rows and columns k+1 and
            -ipiv(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2
            diagonal block.

    @param[in,out]
    B       (input/output) DOUBLE PRECISION array, dimension (ldb,nrhs)
            On entry, the n-by-nrhs right hand side matrix B.
            On exit, if info = 0, the n-by-nrhs solution matrix X.

    @param[in]
    ldb     INTEGER
            The leading dimension of the array B.  ldb >= max(1,n).

    @param[out]
    info    INTEGER
            = 0: successful exit
            < 0: if info = -i, the i-th argument had an illegal value
            > 0: if info = i, D(i,i) is exactly zero.  The factorization
                 has been completed, but the block diagonal matrix D is
                 exactly singular, so the solution could not be computed.

    @ingroup magma_tally3_dsysv_driver
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_dsysv(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs, 
    double *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv, 
    double *B, magma_tally3_int_t ldb, 
    magma_tally3_int_t *info ) 
{
    /* .. Local Scalars .. */
    magma_tally3_int_t upper = (uplo == Magma_tally3Upper);;

    /* Test the input parameters. */
    *info = 0;
    if( !upper && uplo != Magma_tally3Lower ) {
       *info = -1;
    } else if ( n < 0 ) {
       *info = -2;
    } else if ( nrhs < 0 ) {
       *info = -3;
    } else if ( lda < max( 1, n ) ) {
       *info = -5;
    } else if ( ldb < max( 1, n ) ) {
       *info = -8;
    }

    if( *info != 0 ) {
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Compute the factorization A = U*D*U' or A = L*D*L'. */
    magma_tally3_dsytrf( uplo, n, A, lda, ipiv, info );
    if( *info == 0 ) {

        /* Solve the system A*X = B, overwriting B with X. */
        lapackf77_dsytrs( (upper ? Magma_tally3UpperStr: Magma_tally3LowerStr), 
                           &n, &nrhs, A, &lda, ipiv, B, &ldb, info );
    }

    return *info;
    /* End of DSYSV */
}

