/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from ztrsm.cu normal z -> d, Fri Jan 30 19:00:10 2015

       @author Peng Du
       @author Tingxing Dong
       @author Mark Gates
*/
#include "common_magma_tally4.h"
#include "dtrtri.h"  // get NB from dtrtri

/**
    Purpose
    -------
    dtrsm_outofplace solves one of the matrix equations on gpu

        op(A)*X = alpha*B,   or   X*op(A) = alpha*B,

    where alpha is a scalar, X and B are m by n matrices, A is a unit, or
    non-unit, upper or lower triangular matrix and op(A) is one of

        op(A) = A,   or   op(A) = A^T,  or  op(A) = A^H.

    The matrix X is output.

    This is an asynchronous version of magma_tally4blas_dtrsm with flag,
    d_dinvA and dX workspaces as arguments.

    Arguments
    ----------
    @param[in]
    side    magma_tally4_side_t.
            On entry, side specifies whether op(A) appears on the left
            or right of X as follows:
      -     = Magma_tally4Left:       op(A)*X = alpha*B.
      -     = Magma_tally4Right:      X*op(A) = alpha*B.

    @param[in]
    uplo    magma_tally4_uplo_t.
            On entry, uplo specifies whether the matrix A is an upper or
            lower triangular matrix as follows:
      -     = Magma_tally4Upper:  A is an upper triangular matrix.
      -     = Magma_tally4Lower:  A is a  lower triangular matrix.

    @param[in]
    transA  magma_tally4_trans_t.
            On entry, transA specifies the form of op(A) to be used in
            the matrix multiplication as follows:
      -     = Magma_tally4NoTrans:    op(A) = A.
      -     = Magma_tally4Trans:      op(A) = A^T.
      -     = Magma_tally4ConjTrans:  op(A) = A^H.

    @param[in]
    diag    magma_tally4_diag_t.
            On entry, diag specifies whether or not A is unit triangular
            as follows:
      -     = Magma_tally4Unit:     A is assumed to be unit triangular.
      -     = Magma_tally4NonUnit:  A is not assumed to be unit triangular.

    @param[in]
    m       INTEGER.
            On entry, m specifies the number of rows of B. m >= 0.

    @param[in]
    n       INTEGER.
            On entry, n specifies the number of columns of B. n >= 0.

    @param[in]
    alpha   DOUBLE_PRECISION.
            On entry, alpha specifies the scalar alpha. When alpha is
            zero then A is not referenced and B need not be set before
            entry.

    @param[in]
    dA      DOUBLE_PRECISION array of dimension ( ldda, k ), where k is m
            when side = Magma_tally4Left and is n when side = Magma_tally4Right.
            Before entry with uplo = Magma_tally4Upper, the leading k by k
            upper triangular part of the array A must contain the upper
            triangular matrix and the strictly lower triangular part of
            A is not referenced.
            Before entry with uplo = Magma_tally4Lower, the leading k by k
            lower triangular part of the array A must contain the lower
            triangular matrix and the strictly upper triangular part of
            A is not referenced.
            Note that when diag = Magma_tally4Unit, the diagonal elements of
            A are not referenced either, but are assumed to be unity.

    @param[in]
    ldda    INTEGER.
            On entry, ldda specifies the first dimension of A.
            When side = Magma_tally4Left,  ldda >= max( 1, m ),
            when side = Magma_tally4Right, ldda >= max( 1, n ).

    @param[in]
    dB      DOUBLE_PRECISION array of dimension ( lddb, n ).
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
            If side == Magma_tally4Left,  d_dinvA must be of size >= ((m+NB-1)/NB)*NB*NB,
            If side == Magma_tally4Right, d_dinvA must be of size >= ((n+NB-1)/NB)*NB*NB,
            where NB = 128.

    @param[out]
    dX      DOUBLE_PRECISION array of dimension ( m, n ).
            On exit it contain the solution matrix X.

    @ingroup magma_tally4_dblas3
    ********************************************************************/
extern "C"
void magma_tally4blas_dtrsm_outofplace(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t transA, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n,
    double alpha,
    magma_tally4Double_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr       dB, magma_tally4_int_t lddb,
    magma_tally4_int_t flag,
    magma_tally4Double_ptr d_dinvA, magma_tally4Double_ptr dX)
{
    #define dA(i_, j_) (dA + (i_) + (j_)*ldda)
    #define dB(i_, j_) (dB + (i_) + (j_)*lddb)
    #define dX(i_, j_) (dX + (i_) + (j_)*m)
    #define d_dinvA(i_) (d_dinvA + (i_)*NB)

    const double c_neg_one = MAGMA_tally4_D_NEG_ONE;
    const double c_one     = MAGMA_tally4_D_ONE;
    const double c_zero    = MAGMA_tally4_D_ZERO;

    magma_tally4_int_t i, jb;
    magma_tally4_int_t nrowA = (side == Magma_tally4Left ? m : n);

    magma_tally4_int_t info = 0;
    if ( side != Magma_tally4Left && side != Magma_tally4Right ) {
        info = -1;
    } else if ( uplo != Magma_tally4Upper && uplo != Magma_tally4Lower ) {
        info = -2;
    } else if ( transA != Magma_tally4NoTrans && transA != Magma_tally4Trans && transA != Magma_tally4ConjTrans ) {
        info = -3;
    } else if ( diag != Magma_tally4Unit && diag != Magma_tally4NonUnit ) {
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
        magma_tally4_xerbla( __func__, -(info) );
        return;
    }

    // quick return if possible.
    if (m == 0 || n == 0)
        return;

    if (side == Magma_tally4Left) {
        // invert diagonal blocks
        if (flag)
            magma_tally4blas_dtrtri_diag( uplo, diag, m, dA, ldda, d_dinvA );

        if (transA == Magma_tally4NoTrans) {
            if (uplo == Magma_tally4Lower) {
                // left, lower no-transpose
                // handle first block seperately with alpha
                jb = min(NB, m);
                magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, jb, n, jb, alpha, d_dinvA(0), NB, dB, lddb, c_zero, dX, m );
                if (NB < m) {
                    magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, m-NB, n, NB, c_neg_one, dA(NB,0), ldda, dX, m, alpha, dB(NB,0), lddb );

                    // remaining blocks
                    for( i=NB; i < m; i += NB ) {
                        jb = min(m-i, NB);
                        magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, jb, n, jb, c_one, d_dinvA(i), NB, dB(i,0), lddb, c_zero, dX(i,0), m );
                        if (i+NB >= m)
                            break;
                        magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, m-i-NB, n, NB, c_neg_one, dA(i+NB,i), ldda, dX(i,0), m, c_one, dB(i+NB,0), lddb );
                    }
                }
            }
            else {
                // left, upper no-transpose
                // handle first block seperately with alpha
                jb = (m % NB == 0) ? NB : (m % NB);
                i = m-jb;
                magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, jb, n, jb, alpha, d_dinvA(i), NB, dB(i,0), lddb, c_zero, dX(i,0), m );
                if (i-NB >= 0) {
                    magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, i, n, jb, c_neg_one, dA(0,i), ldda, dX(i,0), m, alpha, dB, lddb );

                    // remaining blocks
                    for( i=m-jb-NB; i >= 0; i -= NB ) {
                        magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, NB, n, NB, c_one, d_dinvA(i), NB, dB(i,0), lddb, c_zero, dX(i,0), m );
                        if (i-NB < 0)
                            break;
                        magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, i, n, NB, c_neg_one, dA(0,i), ldda, dX(i,0), m, c_one, dB, lddb );
                    }
                }
            }
        }
        else {  // transA == Magma_tally4Trans || transA == Magma_tally4ConjTrans
            if (uplo == Magma_tally4Lower) {
                // left, lower transpose
                // handle first block seperately with alpha
                jb = (m % NB == 0) ? NB : (m % NB);
                i = m-jb;
                magma_tally4_dgemm( transA, Magma_tally4NoTrans, jb, n, jb, alpha, d_dinvA(i), NB, dB(i,0), lddb, c_zero, dX(i,0), m );
                if (i-NB >= 0) {
                    magma_tally4_dgemm( transA, Magma_tally4NoTrans, i, n, jb, c_neg_one, dA(i,0), ldda, dX(i,0), m, alpha, dB, lddb );

                    // remaining blocks
                    for( i=m-jb-NB; i >= 0; i -= NB ) {
                        magma_tally4_dgemm( transA, Magma_tally4NoTrans, NB, n, NB, c_one, d_dinvA(i), NB, dB(i,0), lddb, c_zero, dX(i,0), m );
                        if (i-NB < 0)
                            break;
                        magma_tally4_dgemm( transA, Magma_tally4NoTrans, i, n, NB, c_neg_one, dA(i,0), ldda, dX(i,0), m, c_one, dB, lddb );
                    }
                }
            }
            else {
                // left, upper transpose
                // handle first block seperately with alpha
                jb = min(NB, m);
                magma_tally4_dgemm( transA, Magma_tally4NoTrans, jb, n, jb, alpha, d_dinvA(0), NB, dB, lddb, c_zero, dX, m );
                if (NB < m) {
                    magma_tally4_dgemm( transA, Magma_tally4NoTrans, m-NB, n, NB, c_neg_one, dA(0,NB), ldda, dX, m, alpha, dB(NB,0), lddb );

                    // remaining blocks
                    for( i=NB; i < m; i += NB ) {
                        jb = min(m-i, NB);
                        magma_tally4_dgemm( transA, Magma_tally4NoTrans, jb, n, jb, c_one, d_dinvA(i), NB, dB(i,0), lddb, c_zero, dX(i,0), m );
                        if (i+NB >= m)
                            break;
                        magma_tally4_dgemm( transA, Magma_tally4NoTrans, m-i-NB, n, NB, c_neg_one, dA(i,i+NB), ldda, dX(i,0), m, c_one, dB(i+NB,0), lddb );
                    }
                }
            }
        }
    }
    else {  // side == Magma_tally4Right
        // invert diagonal blocks
        if (flag)
            magma_tally4blas_dtrtri_diag( uplo, diag, n, dA, ldda, d_dinvA );

        if (transA == Magma_tally4NoTrans) {
            if (uplo == Magma_tally4Lower) {
                // right, lower no-transpose
                // handle first block seperately with alpha
                jb = (n % NB == 0) ? NB : (n % NB);
                i = n-jb;
                magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, m, jb, jb, alpha, dB(0,i), lddb, d_dinvA(i), NB, c_zero, dX(0,i), m );
                if (i-NB >= 0) {
                    magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, m, i, jb, c_neg_one, dX(0,i), m, dA(i,0), ldda, alpha, dB, lddb );

                    // remaining blocks
                    for( i=n-jb-NB; i >= 0; i -= NB ) {
                        magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, m, NB, NB, c_one, dB(0,i), lddb, d_dinvA(i), NB, c_zero, dX(0,i), m );
                        if (i-NB < 0)
                            break;
                        magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, m, i, NB, c_neg_one, dX(0,i), m, dA(i,0), ldda, c_one, dB, lddb );
                    }
                }
            }
            else {
                // right, upper no-transpose
                // handle first block seperately with alpha
                jb = min(NB, n);
                magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, m, jb, jb, alpha, dB, lddb, d_dinvA(0), NB, c_zero, dX, m );
                if (NB < n) {
                    magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, m, n-NB, NB, c_neg_one, dX, m, dA(0,NB), ldda, alpha, dB(0,NB), lddb );

                    // remaining blocks
                    for( i=NB; i < n; i += NB ) {
                        jb = min(NB, n-i);
                        magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, m, jb, jb, c_one, dB(0,i), lddb, d_dinvA(i), NB, c_zero, dX(0,i), m );
                        if (i+NB >= n)
                            break;
                        magma_tally4_dgemm( Magma_tally4NoTrans, Magma_tally4NoTrans, m, n-i-NB, NB, c_neg_one, dX(0,i), m, dA(i,i+NB), ldda, c_one, dB(0,i+NB), lddb );
                    }
                }
            }
        }
        else { // transA == Magma_tally4Trans || transA == Magma_tally4ConjTrans
            if (uplo == Magma_tally4Lower) {
                // right, lower transpose
                // handle first block seperately with alpha
                jb = min(NB, n);
                magma_tally4_dgemm( Magma_tally4NoTrans, transA, m, jb, jb, alpha, dB, lddb, d_dinvA(0), NB, c_zero, dX, m );
                if (NB < n) {
                    magma_tally4_dgemm( Magma_tally4NoTrans, transA, m, n-NB, NB, c_neg_one, dX, m, dA(NB,0), ldda, alpha, dB(0,NB), lddb );

                    // remaining blocks
                    for( i=NB; i < n; i += NB ) {
                        jb = min(NB, n-i);
                        magma_tally4_dgemm( Magma_tally4NoTrans, transA, m, jb, jb, c_one, dB(0,i), lddb, d_dinvA(i), NB, c_zero, dX(0,i), m );
                        if (i+NB >= n)
                            break;
                        magma_tally4_dgemm( Magma_tally4NoTrans, transA, m, n-i-NB, NB, c_neg_one, dX(0,i), m, dA(NB+i,i), ldda, c_one, dB(0,i+NB), lddb );
                    }
                }
            }
            else {
                // right, upper transpose
                // handle first block seperately with alpha
                jb = (n % NB == 0) ? NB : (n % NB);
                i = n-jb;
                magma_tally4_dgemm( Magma_tally4NoTrans, transA, m, jb, jb, alpha, dB(0,i), lddb, d_dinvA(i), NB, c_zero, dX(0,i), m );
                if (i-NB >= 0) {
                    magma_tally4_dgemm( Magma_tally4NoTrans, transA, m, i, jb, c_neg_one, dX(0,i), m, dA(0,i), ldda, alpha, dB, lddb );

                    // remaining blocks
                    for( i=n-jb-NB; i >= 0; i -= NB ) {
                        magma_tally4_dgemm( Magma_tally4NoTrans, transA, m, NB, NB, c_one, dB(0,i), lddb, d_dinvA(i), NB, c_zero, dX(0,i), m );
                        if (i-NB < 0)
                            break;
                        magma_tally4_dgemm( Magma_tally4NoTrans, transA, m, i, NB, c_neg_one, dX(0,i), m, dA(0,i), ldda, c_one, dB, lddb );
                    }
                }
            }
        }
    }
}

/**
    @see magma_tally4blas_dtrsm_outofplace
    @ingroup magma_tally4_dblas3
    ********************************************************************/
extern "C"
void magma_tally4blas_dtrsm_work(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t transA, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n,
    double alpha,
    magma_tally4Double_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr       dB, magma_tally4_int_t lddb,
    magma_tally4_int_t flag,
    magma_tally4Double_ptr d_dinvA, magma_tally4Double_ptr dX)
{

    magma_tally4blas_dtrsm_outofplace( side, uplo, transA, diag, m, n, alpha,
                                dA, ldda, dB, lddb, 1, d_dinvA, dX );
    // copy X to B
    magma_tally4blas_dlacpy( Magma_tally4Full, m, n, dX, m, dB, lddb );
}

/**
    @see magma_tally4blas_dtrsm_work
    @ingroup magma_tally4_dblas3
    ********************************************************************/
extern "C"
void magma_tally4blas_dtrsm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t transA, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n,
    double alpha,
    magma_tally4Double_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr       dB, magma_tally4_int_t lddb )
{
    magma_tally4_int_t nrowA = (side == Magma_tally4Left ? m : n);

    magma_tally4_int_t info = 0;
    if ( side != Magma_tally4Left && side != Magma_tally4Right ) {
        info = -1;
    } else if ( uplo != Magma_tally4Upper && uplo != Magma_tally4Lower ) {
        info = -2;
    } else if ( transA != Magma_tally4NoTrans && transA != Magma_tally4Trans && transA != Magma_tally4ConjTrans ) {
        info = -3;
    } else if ( diag != Magma_tally4Unit && diag != Magma_tally4NonUnit ) {
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
        magma_tally4_xerbla( __func__, -(info) );
        return;
    }

    magma_tally4Double_ptr d_dinvA, dX;
    magma_tally4_int_t size_dinvA;
    magma_tally4_int_t size_x = m*n;
    if ( side == Magma_tally4Left ) {
        size_dinvA = ((m+NB-1)/NB)*NB*NB;
    }
    else {
        size_dinvA = ((n+NB-1)/NB)*NB*NB;
    }

    magma_tally4_dmalloc( &d_dinvA, size_dinvA );
    magma_tally4_dmalloc( &dX, size_x );
    if ( d_dinvA == NULL || dX == NULL ) {
        info = MAGMA_tally4_ERR_DEVICE_ALLOC;
        magma_tally4_xerbla( __func__, -(info) );
    }
    else {
        magma_tally4blas_dlaset(Magma_tally4Full, size_dinvA, 1, MAGMA_tally4_D_ZERO, MAGMA_tally4_D_ZERO, d_dinvA, size_dinvA);
        magma_tally4blas_dlaset(Magma_tally4Full, m, n, MAGMA_tally4_D_ZERO, MAGMA_tally4_D_ZERO, dX, m);
        magma_tally4blas_dtrsm_work( side, uplo, transA, diag, m, n, alpha,
                              dA, ldda, dB, lddb, 1, d_dinvA, dX );
    }

    magma_tally4_free( d_dinvA );
    magma_tally4_free( dX );
}
