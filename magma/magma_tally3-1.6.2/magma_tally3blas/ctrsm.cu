/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from ztrsm.cu normal z -> c, Fri Jan 30 19:00:10 2015

       @author Peng Du
       @author Tingxing Dong
       @author Mark Gates
*/
#include "common_magma_tally3.h"
#include "ctrtri.h"  // get NB from ctrtri

/**
    Purpose
    -------
    ctrsm_outofplace solves one of the matrix equations on gpu

        op(A)*X = alpha*B,   or   X*op(A) = alpha*B,

    where alpha is a scalar, X and B are m by n matrices, A is a unit, or
    non-unit, upper or lower triangular matrix and op(A) is one of

        op(A) = A,   or   op(A) = A^T,  or  op(A) = A^H.

    The matrix X is output.

    This is an asynchronous version of magma_tally3blas_ctrsm with flag,
    d_dinvA and dX workspaces as arguments.

    Arguments
    ----------
    @param[in]
    side    magma_tally3_side_t.
            On entry, side specifies whether op(A) appears on the left
            or right of X as follows:
      -     = Magma_tally3Left:       op(A)*X = alpha*B.
      -     = Magma_tally3Right:      X*op(A) = alpha*B.

    @param[in]
    uplo    magma_tally3_uplo_t.
            On entry, uplo specifies whether the matrix A is an upper or
            lower triangular matrix as follows:
      -     = Magma_tally3Upper:  A is an upper triangular matrix.
      -     = Magma_tally3Lower:  A is a  lower triangular matrix.

    @param[in]
    transA  magma_tally3_trans_t.
            On entry, transA specifies the form of op(A) to be used in
            the matrix multiplication as follows:
      -     = Magma_tally3NoTrans:    op(A) = A.
      -     = Magma_tally3Trans:      op(A) = A^T.
      -     = Magma_tally3ConjTrans:  op(A) = A^H.

    @param[in]
    diag    magma_tally3_diag_t.
            On entry, diag specifies whether or not A is unit triangular
            as follows:
      -     = Magma_tally3Unit:     A is assumed to be unit triangular.
      -     = Magma_tally3NonUnit:  A is not assumed to be unit triangular.

    @param[in]
    m       INTEGER.
            On entry, m specifies the number of rows of B. m >= 0.

    @param[in]
    n       INTEGER.
            On entry, n specifies the number of columns of B. n >= 0.

    @param[in]
    alpha   COMPLEX.
            On entry, alpha specifies the scalar alpha. When alpha is
            zero then A is not referenced and B need not be set before
            entry.

    @param[in]
    dA      COMPLEX array of dimension ( ldda, k ), where k is m
            when side = Magma_tally3Left and is n when side = Magma_tally3Right.
            Before entry with uplo = Magma_tally3Upper, the leading k by k
            upper triangular part of the array A must contain the upper
            triangular matrix and the strictly lower triangular part of
            A is not referenced.
            Before entry with uplo = Magma_tally3Lower, the leading k by k
            lower triangular part of the array A must contain the lower
            triangular matrix and the strictly upper triangular part of
            A is not referenced.
            Note that when diag = Magma_tally3Unit, the diagonal elements of
            A are not referenced either, but are assumed to be unity.

    @param[in]
    ldda    INTEGER.
            On entry, ldda specifies the first dimension of A.
            When side = Magma_tally3Left,  ldda >= max( 1, m ),
            when side = Magma_tally3Right, ldda >= max( 1, n ).

    @param[in]
    dB      COMPLEX array of dimension ( lddb, n ).
            Before entry, the leading m by n part of the array B must
            contain the right-hand side matrix B.

    @param[in]
    lddb    INTEGER.
            On entry, lddb specifies the first dimension of B.
            lddb >= max( 1, m ).

    @param[in]
    flag    BOOLEAN.
            If flag is true, invert diagonal blocks.
            If flag is false, assume diagonal blocks (stored in d_dinvA) are already inverted.

    @param
    d_dinvA (workspace) on device.
            If side == Magma_tally3Left,  d_dinvA must be of size >= ((m+NB-1)/NB)*NB*NB,
            If side == Magma_tally3Right, d_dinvA must be of size >= ((n+NB-1)/NB)*NB*NB,
            where NB = 128.

    @param[out]
    dX      COMPLEX array of dimension ( m, n ).
            On exit it contain the solution matrix X.

    @ingroup magma_tally3_cblas3
    ********************************************************************/
extern "C"
void magma_tally3blas_ctrsm_outofplace(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transA, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr       dB, magma_tally3_int_t lddb,
    magma_tally3_int_t flag,
    magma_tally3FloatComplex_ptr d_dinvA, magma_tally3FloatComplex_ptr dX)
{
    #define dA(i_, j_) (dA + (i_) + (j_)*ldda)
    #define dB(i_, j_) (dB + (i_) + (j_)*lddb)
    #define dX(i_, j_) (dX + (i_) + (j_)*m)
    #define d_dinvA(i_) (d_dinvA + (i_)*NB)

    const magma_tally3FloatComplex c_neg_one = MAGMA_tally3_C_NEG_ONE;
    const magma_tally3FloatComplex c_one     = MAGMA_tally3_C_ONE;
    const magma_tally3FloatComplex c_zero    = MAGMA_tally3_C_ZERO;

    magma_tally3_int_t i, jb;
    magma_tally3_int_t nrowA = (side == Magma_tally3Left ? m : n);

    magma_tally3_int_t info = 0;
    if ( side != Magma_tally3Left && side != Magma_tally3Right ) {
        info = -1;
    } else if ( uplo != Magma_tally3Upper && uplo != Magma_tally3Lower ) {
        info = -2;
    } else if ( transA != Magma_tally3NoTrans && transA != Magma_tally3Trans && transA != Magma_tally3ConjTrans ) {
        info = -3;
    } else if ( diag != Magma_tally3Unit && diag != Magma_tally3NonUnit ) {
        info = -4;
    } else if (m < 0) {
        info = -5;
    } else if (n < 0) {
        info = -6;
    } else if (ldda < max(1,nrowA)) {
        info = -9;
    } else if (lddb < max(1,m)) {
        info = -11;
    }

    if (info != 0) {
        magma_tally3_xerbla( __func__, -(info) );
        return;
    }

    // quick return if possible.
    if (m == 0 || n == 0)
        return;

    if (side == Magma_tally3Left) {
        // invert diagonal blocks
        if (flag)
            magma_tally3blas_ctrtri_diag( uplo, diag, m, dA, ldda, d_dinvA );

        if (transA == Magma_tally3NoTrans) {
            if (uplo == Magma_tally3Lower) {
                // left, lower no-transpose
                // handle first block seperately with alpha
                jb = min(NB, m);
                magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, jb, n, jb, alpha, d_dinvA(0), NB, dB, lddb, c_zero, dX, m );
                if (NB < m) {
                    magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, m-NB, n, NB, c_neg_one, dA(NB,0), ldda, dX, m, alpha, dB(NB,0), lddb );

                    // remaining blocks
                    for( i=NB; i < m; i += NB ) {
                        jb = min(m-i, NB);
                        magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, jb, n, jb, c_one, d_dinvA(i), NB, dB(i,0), lddb, c_zero, dX(i,0), m );
                        if (i+NB >= m)
                            break;
                        magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, m-i-NB, n, NB, c_neg_one, dA(i+NB,i), ldda, dX(i,0), m, c_one, dB(i+NB,0), lddb );
                    }
                }
            }
            else {
                // left, upper no-transpose
                // handle first block seperately with alpha
                jb = (m % NB == 0) ? NB : (m % NB);
                i = m-jb;
                magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, jb, n, jb, alpha, d_dinvA(i), NB, dB(i,0), lddb, c_zero, dX(i,0), m );
                if (i-NB >= 0) {
                    magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, i, n, jb, c_neg_one, dA(0,i), ldda, dX(i,0), m, alpha, dB, lddb );

                    // remaining blocks
                    for( i=m-jb-NB; i >= 0; i -= NB ) {
                        magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, NB, n, NB, c_one, d_dinvA(i), NB, dB(i,0), lddb, c_zero, dX(i,0), m );
                        if (i-NB < 0)
                            break;
                        magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, i, n, NB, c_neg_one, dA(0,i), ldda, dX(i,0), m, c_one, dB, lddb );
                    }
                }
            }
        }
        else {  // transA == Magma_tally3Trans || transA == Magma_tally3ConjTrans
            if (uplo == Magma_tally3Lower) {
                // left, lower transpose
                // handle first block seperately with alpha
                jb = (m % NB == 0) ? NB : (m % NB);
                i = m-jb;
                magma_tally3_cgemm( transA, Magma_tally3NoTrans, jb, n, jb, alpha, d_dinvA(i), NB, dB(i,0), lddb, c_zero, dX(i,0), m );
                if (i-NB >= 0) {
                    magma_tally3_cgemm( transA, Magma_tally3NoTrans, i, n, jb, c_neg_one, dA(i,0), ldda, dX(i,0), m, alpha, dB, lddb );

                    // remaining blocks
                    for( i=m-jb-NB; i >= 0; i -= NB ) {
                        magma_tally3_cgemm( transA, Magma_tally3NoTrans, NB, n, NB, c_one, d_dinvA(i), NB, dB(i,0), lddb, c_zero, dX(i,0), m );
                        if (i-NB < 0)
                            break;
                        magma_tally3_cgemm( transA, Magma_tally3NoTrans, i, n, NB, c_neg_one, dA(i,0), ldda, dX(i,0), m, c_one, dB, lddb );
                    }
                }
            }
            else {
                // left, upper transpose
                // handle first block seperately with alpha
                jb = min(NB, m);
                magma_tally3_cgemm( transA, Magma_tally3NoTrans, jb, n, jb, alpha, d_dinvA(0), NB, dB, lddb, c_zero, dX, m );
                if (NB < m) {
                    magma_tally3_cgemm( transA, Magma_tally3NoTrans, m-NB, n, NB, c_neg_one, dA(0,NB), ldda, dX, m, alpha, dB(NB,0), lddb );

                    // remaining blocks
                    for( i=NB; i < m; i += NB ) {
                        jb = min(m-i, NB);
                        magma_tally3_cgemm( transA, Magma_tally3NoTrans, jb, n, jb, c_one, d_dinvA(i), NB, dB(i,0), lddb, c_zero, dX(i,0), m );
                        if (i+NB >= m)
                            break;
                        magma_tally3_cgemm( transA, Magma_tally3NoTrans, m-i-NB, n, NB, c_neg_one, dA(i,i+NB), ldda, dX(i,0), m, c_one, dB(i+NB,0), lddb );
                    }
                }
            }
        }
    }
    else {  // side == Magma_tally3Right
        // invert diagonal blocks
        if (flag)
            magma_tally3blas_ctrtri_diag( uplo, diag, n, dA, ldda, d_dinvA );

        if (transA == Magma_tally3NoTrans) {
            if (uplo == Magma_tally3Lower) {
                // right, lower no-transpose
                // handle first block seperately with alpha
                jb = (n % NB == 0) ? NB : (n % NB);
                i = n-jb;
                magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, m, jb, jb, alpha, dB(0,i), lddb, d_dinvA(i), NB, c_zero, dX(0,i), m );
                if (i-NB >= 0) {
                    magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, m, i, jb, c_neg_one, dX(0,i), m, dA(i,0), ldda, alpha, dB, lddb );

                    // remaining blocks
                    for( i=n-jb-NB; i >= 0; i -= NB ) {
                        magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, m, NB, NB, c_one, dB(0,i), lddb, d_dinvA(i), NB, c_zero, dX(0,i), m );
                        if (i-NB < 0)
                            break;
                        magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, m, i, NB, c_neg_one, dX(0,i), m, dA(i,0), ldda, c_one, dB, lddb );
                    }
                }
            }
            else {
                // right, upper no-transpose
                // handle first block seperately with alpha
                jb = min(NB, n);
                magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, m, jb, jb, alpha, dB, lddb, d_dinvA(0), NB, c_zero, dX, m );
                if (NB < n) {
                    magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, m, n-NB, NB, c_neg_one, dX, m, dA(0,NB), ldda, alpha, dB(0,NB), lddb );

                    // remaining blocks
                    for( i=NB; i < n; i += NB ) {
                        jb = min(NB, n-i);
                        magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, m, jb, jb, c_one, dB(0,i), lddb, d_dinvA(i), NB, c_zero, dX(0,i), m );
                        if (i+NB >= n)
                            break;
                        magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, m, n-i-NB, NB, c_neg_one, dX(0,i), m, dA(i,i+NB), ldda, c_one, dB(0,i+NB), lddb );
                    }
                }
            }
        }
        else { // transA == Magma_tally3Trans || transA == Magma_tally3ConjTrans
            if (uplo == Magma_tally3Lower) {
                // right, lower transpose
                // handle first block seperately with alpha
                jb = min(NB, n);
                magma_tally3_cgemm( Magma_tally3NoTrans, transA, m, jb, jb, alpha, dB, lddb, d_dinvA(0), NB, c_zero, dX, m );
                if (NB < n) {
                    magma_tally3_cgemm( Magma_tally3NoTrans, transA, m, n-NB, NB, c_neg_one, dX, m, dA(NB,0), ldda, alpha, dB(0,NB), lddb );

                    // remaining blocks
                    for( i=NB; i < n; i += NB ) {
                        jb = min(NB, n-i);
                        magma_tally3_cgemm( Magma_tally3NoTrans, transA, m, jb, jb, c_one, dB(0,i), lddb, d_dinvA(i), NB, c_zero, dX(0,i), m );
                        if (i+NB >= n)
                            break;
                        magma_tally3_cgemm( Magma_tally3NoTrans, transA, m, n-i-NB, NB, c_neg_one, dX(0,i), m, dA(NB+i,i), ldda, c_one, dB(0,i+NB), lddb );
                    }
                }
            }
            else {
                // right, upper transpose
                // handle first block seperately with alpha
                jb = (n % NB == 0) ? NB : (n % NB);
                i = n-jb;
                magma_tally3_cgemm( Magma_tally3NoTrans, transA, m, jb, jb, alpha, dB(0,i), lddb, d_dinvA(i), NB, c_zero, dX(0,i), m );
                if (i-NB >= 0) {
                    magma_tally3_cgemm( Magma_tally3NoTrans, transA, m, i, jb, c_neg_one, dX(0,i), m, dA(0,i), ldda, alpha, dB, lddb );

                    // remaining blocks
                    for( i=n-jb-NB; i >= 0; i -= NB ) {
                        magma_tally3_cgemm( Magma_tally3NoTrans, transA, m, NB, NB, c_one, dB(0,i), lddb, d_dinvA(i), NB, c_zero, dX(0,i), m );
                        if (i-NB < 0)
                            break;
                        magma_tally3_cgemm( Magma_tally3NoTrans, transA, m, i, NB, c_neg_one, dX(0,i), m, dA(0,i), ldda, c_one, dB, lddb );
                    }
                }
            }
        }
    }
}

/**
    @see magma_tally3blas_ctrsm_outofplace
    @ingroup magma_tally3_cblas3
    ********************************************************************/
extern "C"
void magma_tally3blas_ctrsm_work(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transA, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr       dB, magma_tally3_int_t lddb,
    magma_tally3_int_t flag,
    magma_tally3FloatComplex_ptr d_dinvA, magma_tally3FloatComplex_ptr dX)
{

    magma_tally3blas_ctrsm_outofplace( side, uplo, transA, diag, m, n, alpha,
                                dA, ldda, dB, lddb, 1, d_dinvA, dX );
    // copy X to B
    magma_tally3blas_clacpy( Magma_tally3Full, m, n, dX, m, dB, lddb );
}

/**
    @see magma_tally3blas_ctrsm_work
    @ingroup magma_tally3_cblas3
    ********************************************************************/
extern "C"
void magma_tally3blas_ctrsm(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transA, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr       dB, magma_tally3_int_t lddb )
{
    magma_tally3_int_t nrowA = (side == Magma_tally3Left ? m : n);

    magma_tally3_int_t info = 0;
    if ( side != Magma_tally3Left && side != Magma_tally3Right ) {
        info = -1;
    } else if ( uplo != Magma_tally3Upper && uplo != Magma_tally3Lower ) {
        info = -2;
    } else if ( transA != Magma_tally3NoTrans && transA != Magma_tally3Trans && transA != Magma_tally3ConjTrans ) {
        info = -3;
    } else if ( diag != Magma_tally3Unit && diag != Magma_tally3NonUnit ) {
        info = -4;
    } else if (m < 0) {
        info = -5;
    } else if (n < 0) {
        info = -6;
    } else if (ldda < max(1,nrowA)) {
        info = -9;
    } else if (lddb < max(1,m)) {
        info = -11;
    }

    if (info != 0) {
        magma_tally3_xerbla( __func__, -(info) );
        return;
    }

    magma_tally3FloatComplex_ptr d_dinvA, dX;
    magma_tally3_int_t size_dinvA;
    magma_tally3_int_t size_x = m*n;
    if ( side == Magma_tally3Left ) {
        size_dinvA = ((m+NB-1)/NB)*NB*NB;
    }
    else {
        size_dinvA = ((n+NB-1)/NB)*NB*NB;
    }

    magma_tally3_cmalloc( &d_dinvA, size_dinvA );
    magma_tally3_cmalloc( &dX, size_x );
    if ( d_dinvA == NULL || dX == NULL ) {
        info = MAGMA_tally3_ERR_DEVICE_ALLOC;
        magma_tally3_xerbla( __func__, -(info) );
    }
    else {
        magma_tally3blas_claset(Magma_tally3Full, size_dinvA, 1, MAGMA_tally3_C_ZERO, MAGMA_tally3_C_ZERO, d_dinvA, size_dinvA);
        magma_tally3blas_claset(Magma_tally3Full, m, n, MAGMA_tally3_C_ZERO, MAGMA_tally3_C_ZERO, dX, m);
        magma_tally3blas_ctrsm_work( side, uplo, transA, diag, m, n, alpha,
                              dA, ldda, dB, lddb, 1, d_dinvA, dX );
    }

    magma_tally3_free( d_dinvA );
    magma_tally3_free( dX );
}
