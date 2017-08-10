/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from ztrsm_batched.cu normal z -> c, Fri Jan 30 19:00:10 2015

       @author Peng Du
       @author Tingxing Dong
       @author Mark Gates
       @author Azzam Haidar
*/
#include "common_magma_tally2.h"
#include "batched_kernel_param.h"
/**
    Purpose
    -------
    ctrsm_work solves one of the matrix equations on gpu

        op(A)*X = alpha*B,   or   X*op(A) = alpha*B,

    where alpha is a scalar, X and B are m by n matrices, A is a unit, or
    non-unit, upper or lower triangular matrix and op(A) is one of

        op(A) = A,   or   op(A) = A^T,  or  op(A) = A^H.

    The matrix X is overwritten on B.

    This is an asynchronous version of magma_tally2blas_ctrsm with flag,
    d_dinvA and dX workspaces as arguments.

    Arguments
    ----------
    @param[in]
    side    magma_tally2_side_t.
            On entry, side specifies whether op(A) appears on the left
            or right of X as follows:
      -     = Magma_tally2Left:       op(A)*X = alpha*B.
      -     = Magma_tally2Right:      X*op(A) = alpha*B.

    @param[in]
    uplo    magma_tally2_uplo_t.
            On entry, uplo specifies whether the matrix A is an upper or
            lower triangular matrix as follows:
      -     = Magma_tally2Upper:  A is an upper triangular matrix.
      -     = Magma_tally2Lower:  A is a  lower triangular matrix.

    @param[in]
    transA  magma_tally2_trans_t.
            On entry, transA specifies the form of op(A) to be used in
            the matrix multiplication as follows:
      -     = Magma_tally2NoTrans:    op(A) = A.
      -     = Magma_tally2Trans:      op(A) = A^T.
      -     = Magma_tally2ConjTrans:  op(A) = A^H.

    @param[in]
    diag    magma_tally2_diag_t.
            On entry, diag specifies whether or not A is unit triangular
            as follows:
      -     = Magma_tally2Unit:     A is assumed to be unit triangular.
      -     = Magma_tally2NonUnit:  A is not assumed to be unit triangular.

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
            when side = Magma_tally2Left and is n when side = Magma_tally2Right.
            Before entry with uplo = Magma_tally2Upper, the leading k by k
            upper triangular part of the array A must contain the upper
            triangular matrix and the strictly lower triangular part of
            A is not referenced.
            Before entry with uplo = Magma_tally2Lower, the leading k by k
            lower triangular part of the array A must contain the lower
            triangular matrix and the strictly upper triangular part of
            A is not referenced.
            Note that when diag = Magma_tally2Unit, the diagonal elements of
            A are not referenced either, but are assumed to be unity.

    @param[in]
    ldda    INTEGER.
            On entry, ldda specifies the first dimension of A.
            When side = Magma_tally2Left,  ldda >= max( 1, m ),
            when side = Magma_tally2Right, ldda >= max( 1, n ).

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
            If side == Magma_tally2Left,  d_dinvA must be of size >= ((m+TRI_NB-1)/TRI_NB)*TRI_NB*TRI_NB,
            If side == Magma_tally2Right, d_dinvA must be of size >= ((n+TRI_NB-1)/TRI_NB)*TRI_NB*TRI_NB,
            where TRI_NB = 128.

    @param[out]
    dX      COMPLEX array of dimension ( lddx, n ).
            On exit it contain the solution matrix X.

    @ingroup magma_tally2_cblas3
    ********************************************************************/
extern "C"
void magma_tally2blas_ctrsm_outofplace_batched(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t transA, magma_tally2_diag_t diag,
    magma_tally2_int_t flag, magma_tally2_int_t m, magma_tally2_int_t n, 
    magma_tally2FloatComplex alpha, 
    magma_tally2FloatComplex** dA_array,    magma_tally2_int_t ldda,
    magma_tally2FloatComplex** dB_array,    magma_tally2_int_t lddb,
    magma_tally2FloatComplex** dX_array,    magma_tally2_int_t lddx, 
    magma_tally2FloatComplex** dinvA_array, magma_tally2_int_t dinvA_length,
    magma_tally2FloatComplex** dA_displ, magma_tally2FloatComplex** dB_displ, 
    magma_tally2FloatComplex** dX_displ, magma_tally2FloatComplex** dinvA_displ,
    magma_tally2_int_t resetozero, magma_tally2_int_t batchCount, magma_tally2_queue_t queue)
{
/*
    #define dA(i_, j_) (dA + (i_) + (j_)*ldda)
    #define dB(i_, j_) (dB + (i_) + (j_)*lddb)
    #define dX(i_, j_) (dX + (i_) + (j_)*m)
    #define d_dinvA(i_) (d_dinvA + (i_)*TRI_NB)

*/
    const magma_tally2FloatComplex c_neg_one = MAGMA_tally2_C_NEG_ONE;
    const magma_tally2FloatComplex c_one     = MAGMA_tally2_C_ONE;
    const magma_tally2FloatComplex c_zero    = MAGMA_tally2_C_ZERO;

    magma_tally2_int_t i, jb;
    magma_tally2_int_t nrowA = (side == Magma_tally2Left ? m : n);

    magma_tally2_int_t info = 0;
    if ( side != Magma_tally2Left && side != Magma_tally2Right ) {
        info = -1;
    } else if ( uplo != Magma_tally2Upper && uplo != Magma_tally2Lower ) {
        info = -2;
    } else if ( transA != Magma_tally2NoTrans && transA != Magma_tally2Trans && transA != Magma_tally2ConjTrans ) {
        info = -3;
    } else if ( diag != Magma_tally2Unit && diag != Magma_tally2NonUnit ) {
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
    magma_tally2_int_t size_dinvA;
    if ( side == Magma_tally2Left ) {
        size_dinvA = ((m+TRI_NB-1)/TRI_NB)*TRI_NB*TRI_NB;
    }
    else {
        size_dinvA = ((n+TRI_NB-1)/TRI_NB)*TRI_NB*TRI_NB;
    }
    if(dinvA_length < size_dinvA) info = -19;

    if (info != 0) {
        magma_tally2_xerbla( __func__, -(info) );
        return;
    }

    // quick return if possible.
    if (m == 0 || n == 0)
        return;

    magma_tally2_cdisplace_pointers(dA_displ,       dA_array,    ldda,    0, 0, batchCount, queue); 
    magma_tally2_cdisplace_pointers(dB_displ,       dB_array,    lddb,    0, 0, batchCount, queue); 
    magma_tally2_cdisplace_pointers(dX_displ,       dX_array,    lddx,    0, 0, batchCount, queue); 
    magma_tally2_cdisplace_pointers(dinvA_displ, dinvA_array,  TRI_NB,    0, 0, batchCount, queue); 

    if (side == Magma_tally2Left) {
        // invert diagonal blocks
        if (flag)
            magma_tally2blas_ctrtri_diag_batched( uplo, diag, m, dA_displ, ldda, dinvA_displ, resetozero, batchCount, queue );

        if (transA == Magma_tally2NoTrans) {
            if (uplo == Magma_tally2Lower) {
                // left, lower no-transpose
                // handle first block seperately with alpha
                jb = min(TRI_NB, m);                
                magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, jb, n, jb, alpha, dinvA_displ, TRI_NB, dB_displ, lddb, c_zero, dX_displ, lddx, batchCount, queue );

                if (TRI_NB < m) {
                    magma_tally2_cdisplace_pointers(dA_displ,    dA_array, ldda, TRI_NB, 0, batchCount, queue); 
                    magma_tally2_cdisplace_pointers(dB_displ,    dB_array, lddb, TRI_NB, 0, batchCount, queue);                    
                    magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, m-TRI_NB, n, TRI_NB, c_neg_one, dA_displ, ldda, dX_displ, lddx, alpha, dB_displ, lddb, batchCount, queue );

                    // remaining blocks
                    for( i=TRI_NB; i < m; i += TRI_NB ) {
                        jb = min(m-i, TRI_NB);
                        magma_tally2_cdisplace_pointers(dinvA_displ, dinvA_array, TRI_NB, 0, i, batchCount, queue);
                        magma_tally2_cdisplace_pointers(dB_displ,       dB_array, lddb,    i, 0, batchCount, queue);
                        magma_tally2_cdisplace_pointers(dX_displ,       dX_array, lddx,    i, 0, batchCount, queue);
                        magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, jb, n, jb, c_one, dinvA_displ, TRI_NB, dB_displ, lddb, c_zero, dX_displ, lddx, batchCount, queue );
                        if (i+TRI_NB >= m)
                            break;
                        
                        magma_tally2_cdisplace_pointers(dA_displ,    dA_array,    ldda,  i+TRI_NB, i, batchCount, queue); 
                        magma_tally2_cdisplace_pointers(dB_displ,    dB_array,    lddb,  i+TRI_NB, 0, batchCount, queue);
                        magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, m-i-TRI_NB, n, TRI_NB, c_neg_one, dA_displ, ldda, dX_displ, lddx, c_one, dB_displ, lddb, batchCount, queue );
                    }
                }
            }
            else {
                // left, upper no-transpose
                // handle first block seperately with alpha
                jb = (m % TRI_NB == 0) ? TRI_NB : (m % TRI_NB);
                i = m-jb;
                magma_tally2_cdisplace_pointers(dinvA_displ,    dinvA_array, TRI_NB, 0, i, batchCount, queue);
                magma_tally2_cdisplace_pointers(dB_displ,          dB_array,    lddb, i, 0, batchCount, queue);
                magma_tally2_cdisplace_pointers(dX_displ,          dX_array,    lddx, i, 0, batchCount, queue);
                magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, jb, n, jb, alpha, dinvA_displ, TRI_NB, dB_displ, lddb, c_zero, dX_displ, lddx, batchCount, queue );

                if (i-TRI_NB >= 0) {
                    magma_tally2_cdisplace_pointers(dA_displ,    dA_array,    ldda, 0, i, batchCount, queue); 
                    magma_tally2_cdisplace_pointers(dB_displ,    dB_array,    lddb, 0, 0, batchCount, queue);
                    magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, i, n, jb, c_neg_one, dA_displ, ldda, dX_displ, lddx, alpha, dB_displ, lddb, batchCount, queue );

                    // remaining blocks
                    for( i=m-jb-TRI_NB; i >= 0; i -= TRI_NB ) {
                        magma_tally2_cdisplace_pointers(dinvA_displ, dinvA_array, TRI_NB, 0, i, batchCount, queue);
                        magma_tally2_cdisplace_pointers(dB_displ,       dB_array,    lddb, i, 0, batchCount, queue);
                        magma_tally2_cdisplace_pointers(dX_displ,       dX_array,    lddx, i, 0, batchCount, queue);
                        magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, TRI_NB, n, TRI_NB, c_one, dinvA_displ, TRI_NB, dB_displ, lddb, c_zero, dX_displ, lddx, batchCount, queue );
                        if (i-TRI_NB < 0)
                            break;

                        magma_tally2_cdisplace_pointers(dA_displ,    dA_array,    ldda, 0, i, batchCount, queue); 
                        magma_tally2_cdisplace_pointers(dB_displ,    dB_array,    lddb, 0, 0, batchCount, queue);
                        magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, i, n, TRI_NB, c_neg_one, dA_displ, ldda, dX_displ, lddx, c_one, dB_displ, lddb, batchCount, queue);
                    }
                }
            }
        }
        else {  // transA == Magma_tally2Trans || transA == Magma_tally2ConjTrans
            if (uplo == Magma_tally2Lower) {
                // left, lower transpose
                // handle first block seperately with alpha
                jb = (m % TRI_NB == 0) ? TRI_NB : (m % TRI_NB);
                i = m-jb;
                magma_tally2_cdisplace_pointers(dinvA_displ,    dinvA_array, TRI_NB, 0, i, batchCount, queue);
                magma_tally2_cdisplace_pointers(dB_displ,          dB_array,    lddb, i, 0, batchCount, queue);
                magma_tally2_cdisplace_pointers(dX_displ,          dX_array,    lddx, i, 0, batchCount, queue);
                magma_tally2blas_cgemm_batched( transA, Magma_tally2NoTrans, jb, n, jb, alpha, dinvA_displ, TRI_NB, dB_displ, lddb, c_zero, dX_displ, lddx, batchCount, queue );

                if (i-TRI_NB >= 0) {
                    magma_tally2_cdisplace_pointers(dA_displ,    dA_array,    ldda, i, 0, batchCount, queue); 
                    magma_tally2_cdisplace_pointers(dB_displ,    dB_array,    lddb, 0, 0, batchCount, queue);
                    magma_tally2blas_cgemm_batched( transA, Magma_tally2NoTrans, i, n, jb, c_neg_one, dA_displ, ldda, dX_displ, lddx, alpha, dB_displ, lddb, batchCount, queue );

                    // remaining blocks
                    for( i=m-jb-TRI_NB; i >= 0; i -= TRI_NB ) {
                        magma_tally2_cdisplace_pointers(dinvA_displ, dinvA_array, TRI_NB, 0, i, batchCount, queue);
                        magma_tally2_cdisplace_pointers(dX_displ,       dX_array,    lddx, i, 0, batchCount, queue);
                        magma_tally2_cdisplace_pointers(dB_displ,       dB_array,    lddb, i, 0, batchCount, queue);
                        magma_tally2blas_cgemm_batched( transA, Magma_tally2NoTrans, TRI_NB, n, TRI_NB, c_one, dinvA_displ, TRI_NB, dB_displ, lddb, c_zero, dX_displ, lddx, batchCount, queue );
                        if (i-TRI_NB < 0)
                            break;

                        magma_tally2_cdisplace_pointers(dA_displ,    dA_array,    ldda, i, 0,  batchCount, queue); 
                        magma_tally2_cdisplace_pointers(dB_displ,    dB_array,    lddb, 0, 0, batchCount, queue);
                        magma_tally2blas_cgemm_batched( transA, Magma_tally2NoTrans, i, n, TRI_NB, c_neg_one, dA_displ, ldda, dX_displ, lddx, c_one, dB_displ, lddb, batchCount, queue );
                    }
                }
            }
            else {
                // left, upper transpose
                // handle first block seperately with alpha
                jb = min(TRI_NB, m);
                magma_tally2blas_cgemm_batched( transA, Magma_tally2NoTrans, jb, n, jb, alpha, dinvA_displ, TRI_NB, dB_displ, lddb, c_zero, dX_displ, lddx, batchCount, queue );

                if (TRI_NB < m) {
                    magma_tally2_cdisplace_pointers(dA_displ,    dA_array,    ldda,      0,   TRI_NB, batchCount, queue); 
                    magma_tally2_cdisplace_pointers(dB_displ,    dB_array,    lddb, TRI_NB,        0, batchCount, queue);
                    magma_tally2blas_cgemm_batched( transA, Magma_tally2NoTrans, m-TRI_NB, n, TRI_NB, c_neg_one, dA_displ, ldda, dX_displ, lddx, alpha, dB_displ, lddb, batchCount, queue );

                    // remaining blocks
                    for( i=TRI_NB; i < m; i += TRI_NB ) {
                        jb = min(m-i, TRI_NB);
                        magma_tally2_cdisplace_pointers(dinvA_displ, dinvA_array, TRI_NB, 0, i, batchCount, queue);
                        magma_tally2_cdisplace_pointers(dB_displ,       dB_array,    lddb, i, 0, batchCount, queue);
                        magma_tally2_cdisplace_pointers(dX_displ,       dX_array,    lddx, i, 0, batchCount, queue);
                        magma_tally2blas_cgemm_batched( transA, Magma_tally2NoTrans, jb, n, jb, c_one, dinvA_displ, TRI_NB, dB_displ, lddb, c_zero, dX_displ, lddx, batchCount, queue );
                        if (i+TRI_NB >= m)
                            break;

                        magma_tally2_cdisplace_pointers(dA_displ,    dA_array,    ldda,  i,        i+TRI_NB, batchCount, queue); 
                        magma_tally2_cdisplace_pointers(dB_displ,    dB_array,    lddb,  i+TRI_NB,        0, batchCount, queue);
                        magma_tally2blas_cgemm_batched( transA, Magma_tally2NoTrans, m-i-TRI_NB, n, TRI_NB, c_neg_one, dA_displ, ldda, dX_displ, lddx, c_one, dB_displ, lddb, batchCount, queue );
                    }
                }
            }
        }
    }
    else {  // side == Magma_tally2Right
        // invert diagonal blocks
        if (flag)
            magma_tally2blas_ctrtri_diag_batched( uplo, diag, n, dA_displ, ldda, dinvA_displ, resetozero, batchCount, queue);

        if (transA == Magma_tally2NoTrans) {
            if (uplo == Magma_tally2Lower) {
                // right, lower no-transpose
                // handle first block seperately with alpha
                jb = (n % TRI_NB == 0) ? TRI_NB : (n % TRI_NB);
                i = n-jb;                
                magma_tally2_cdisplace_pointers(dinvA_displ, dinvA_array, TRI_NB, 0, i, batchCount, queue);
                magma_tally2_cdisplace_pointers(dB_displ,       dB_array,    lddb, 0, i, batchCount, queue);
                magma_tally2_cdisplace_pointers(dX_displ,       dX_array,    lddx, 0, i, batchCount, queue);
                magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, m, jb, jb, alpha, dB_displ, lddb, dinvA_displ, TRI_NB, c_zero, dX_displ, lddx, batchCount, queue );

                if (i-TRI_NB >= 0) {
                    magma_tally2_cdisplace_pointers(dA_displ,    dA_array,    ldda, i, 0, batchCount, queue);                        
                    magma_tally2_cdisplace_pointers(dB_displ,    dB_array,    lddb, 0, 0, batchCount, queue);
                    magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, m, i, jb, c_neg_one, dX_displ, lddx, dA_displ, ldda, alpha, dB_displ, lddb, batchCount, queue );

                    // remaining blocks
                    for( i=n-jb-TRI_NB; i >= 0; i -= TRI_NB ) {
                        magma_tally2_cdisplace_pointers(dinvA_displ, dinvA_array, TRI_NB, 0, i, batchCount, queue);
                        magma_tally2_cdisplace_pointers(dB_displ,       dB_array,    lddb, 0, i, batchCount, queue);
                        magma_tally2_cdisplace_pointers(dX_displ,       dX_array,    lddx, 0, i, batchCount, queue);
                        magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, m, TRI_NB, TRI_NB, c_one, dB_displ, lddb, dinvA_displ, TRI_NB, c_zero, dX_displ, lddx, batchCount, queue );
                        if (i-TRI_NB < 0)
                            break;

                        magma_tally2_cdisplace_pointers(dA_displ,    dA_array,    ldda, i, 0, batchCount, queue); 
                        magma_tally2_cdisplace_pointers(dB_displ,    dB_array,    lddb, 0, 0, batchCount, queue);
                        magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, m, i, TRI_NB, c_neg_one, dX_displ, lddx, dA_displ, ldda, c_one, dB_displ, lddb, batchCount, queue );
                    }
                }
            }
            else {
                // right, upper no-transpose
                // handle first block seperately with alpha
                jb = min(TRI_NB, n);
                magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, m, jb, jb, alpha, dB_displ, lddb, dinvA_displ, TRI_NB, c_zero, dX_displ, lddx, batchCount, queue );
                if (TRI_NB < n) {
                    magma_tally2_cdisplace_pointers(dA_displ,    dA_array,    ldda, 0, TRI_NB, batchCount, queue); 
                    magma_tally2_cdisplace_pointers(dB_displ,    dB_array,    lddb, 0, TRI_NB, batchCount, queue);
                    magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, m, n-TRI_NB, TRI_NB, c_neg_one, dX_displ, lddx, dA_displ, ldda, alpha, dB_displ, lddb, batchCount, queue );

                    // remaining blocks
                    for( i=TRI_NB; i < n; i += TRI_NB ) {
                        jb = min(TRI_NB, n-i);
                        magma_tally2_cdisplace_pointers(dinvA_displ, dinvA_array, TRI_NB, 0, i, batchCount, queue);
                        magma_tally2_cdisplace_pointers(dB_displ,       dB_array,    lddb, 0, i, batchCount, queue);
                        magma_tally2_cdisplace_pointers(dX_displ,       dX_array,    lddx, 0, i, batchCount, queue);
                        magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, m, jb, jb, c_one, dB_displ, lddb, dinvA_displ, TRI_NB, c_zero, dX_displ, lddx, batchCount, queue );
                        if (i+TRI_NB >= n)
                            break;

                        magma_tally2_cdisplace_pointers(dA_displ,    dA_array,    ldda, i, i+TRI_NB, batchCount, queue); 
                        magma_tally2_cdisplace_pointers(dB_displ,    dB_array,    lddb, 0, i+TRI_NB, batchCount, queue);
                        magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, m, n-i-TRI_NB, TRI_NB, c_neg_one, dX_displ, lddx, dA_displ, ldda, c_one, dB_displ, lddb, batchCount, queue );
                    }
                }
            }
        }
        else { // transA == Magma_tally2Trans || transA == Magma_tally2ConjTrans
            if (uplo == Magma_tally2Lower) {
                // right, lower transpose
                // handle first block seperately with alpha
                jb = min(TRI_NB, n);
                magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, transA, m, jb, jb, alpha, dB_displ, lddb, dinvA_displ, TRI_NB, c_zero, dX_displ, lddx, batchCount, queue );
                if (TRI_NB < n) {
                    magma_tally2_cdisplace_pointers(dA_displ,    dA_array,    ldda,  TRI_NB,      0, batchCount, queue); 
                    magma_tally2_cdisplace_pointers(dB_displ,    dB_array,    lddb,       0, TRI_NB, batchCount, queue);
                    magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, transA, m, n-TRI_NB, TRI_NB, c_neg_one, dX_displ, lddx, dA_displ, ldda, alpha, dB_displ, lddb, batchCount, queue );

                    // remaining blocks
                    for( i=TRI_NB; i < n; i += TRI_NB ) {
                        jb = min(TRI_NB, n-i);
                        magma_tally2_cdisplace_pointers(dinvA_displ, dinvA_array, TRI_NB, 0, i, batchCount, queue);
                        magma_tally2_cdisplace_pointers(dB_displ,       dB_array,    lddb, 0, i, batchCount, queue);
                        magma_tally2_cdisplace_pointers(dX_displ,       dX_array,    lddx, 0, i, batchCount, queue);
                        magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, transA, m, jb, jb, c_one, dB_displ, lddb, dinvA_displ, TRI_NB, c_zero, dX_displ, lddx, batchCount, queue );
                        if (i+TRI_NB >= n)
                            break;

                        magma_tally2_cdisplace_pointers(dA_displ,    dA_array,    ldda,  TRI_NB+i,        i, batchCount, queue); 
                        magma_tally2_cdisplace_pointers(dB_displ,    dB_array,    lddb,         0, i+TRI_NB, batchCount, queue);
                        magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, transA, m, n-i-TRI_NB, TRI_NB, c_neg_one, dX_displ, lddx, dA_displ, ldda, c_one, dB_displ, lddb, batchCount, queue );
                    }
                }
            }
            else {
                // right, upper transpose
                // handle first block seperately with alpha
                jb = (n % TRI_NB == 0) ? TRI_NB : (n % TRI_NB);
                i = n-jb;
                magma_tally2_cdisplace_pointers(dinvA_displ, dinvA_array, TRI_NB, 0, i, batchCount, queue);
                magma_tally2_cdisplace_pointers(dB_displ,       dB_array,    lddb, 0, i, batchCount, queue);
                magma_tally2_cdisplace_pointers(dX_displ,       dX_array,    lddx, 0, i, batchCount, queue);
                magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, transA, m, jb, jb, alpha, dB_displ, lddb, dinvA_displ, TRI_NB, c_zero, dX_displ, lddx, batchCount, queue );

                if (i-TRI_NB >= 0) {
                    magma_tally2_cdisplace_pointers(dA_displ,    dA_array,    ldda, 0, i, batchCount, queue); 
                    magma_tally2_cdisplace_pointers(dB_displ,    dB_array,    lddb, 0, 0, batchCount, queue);
                    magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, transA, m, i, jb, c_neg_one, dX_displ, lddx, dA_displ, ldda, alpha, dB_displ, lddb, batchCount, queue );

                    // remaining blocks
                    for( i=n-jb-TRI_NB; i >= 0; i -= TRI_NB ) {
                        magma_tally2_cdisplace_pointers(dinvA_displ, dinvA_array, TRI_NB, 0, i, batchCount, queue);
                        magma_tally2_cdisplace_pointers(dB_displ,       dB_array,    lddb, 0, i, batchCount, queue);
                        magma_tally2_cdisplace_pointers(dX_displ,       dX_array,    lddx, 0, i, batchCount, queue);
                        magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, transA, m, TRI_NB, TRI_NB, c_one, dB_displ, lddb, dinvA_displ, TRI_NB, c_zero, dX_displ, lddx, batchCount, queue );
                        if (i-TRI_NB < 0)
                            break;

                        magma_tally2_cdisplace_pointers(dA_displ,    dA_array,    ldda, 0, i, batchCount, queue); 
                        magma_tally2_cdisplace_pointers(dB_displ,    dB_array,    lddb, 0, 0, batchCount, queue);
                        magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, transA, m, i, TRI_NB, c_neg_one, dX_displ, lddx, dA_displ, ldda, c_one, dB_displ, lddb, batchCount, queue );
                    }
                }
            }
        }
    }
}

/**
    @see magma_tally2blas_ctrsm_outofplace_batched
    @ingroup magma_tally2_cblas3
    ********************************************************************/
extern "C"
void magma_tally2blas_ctrsm_work_batched(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t transA, magma_tally2_diag_t diag,
    magma_tally2_int_t flag, magma_tally2_int_t m, magma_tally2_int_t n, 
    magma_tally2FloatComplex alpha, 
    magma_tally2FloatComplex** dA_array,    magma_tally2_int_t ldda,
    magma_tally2FloatComplex** dB_array,    magma_tally2_int_t lddb,
    magma_tally2FloatComplex** dX_array,    magma_tally2_int_t lddx, 
    magma_tally2FloatComplex** dinvA_array, magma_tally2_int_t dinvA_length,
    magma_tally2FloatComplex** dA_displ, magma_tally2FloatComplex** dB_displ, 
    magma_tally2FloatComplex** dX_displ, magma_tally2FloatComplex** dinvA_displ,
    magma_tally2_int_t resetozero, magma_tally2_int_t batchCount, magma_tally2_queue_t queue)
{
    magma_tally2_int_t nrowA = (side == Magma_tally2Left ? m : n);

    magma_tally2_int_t info = 0;
    if ( side != Magma_tally2Left && side != Magma_tally2Right ) {
        info = -1;
    } else if ( uplo != Magma_tally2Upper && uplo != Magma_tally2Lower ) {
        info = -2;
    } else if ( transA != Magma_tally2NoTrans && transA != Magma_tally2Trans && transA != Magma_tally2ConjTrans ) {
        info = -3;
    } else if ( diag != Magma_tally2Unit && diag != Magma_tally2NonUnit ) {
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
        magma_tally2_xerbla( __func__, -(info) );
        return;
    }



    magma_tally2blas_ctrsm_outofplace_batched( 
                    side, uplo, transA, diag, flag,
                    m, n, alpha,
                    dA_array,    ldda,
                    dB_array,    lddb,
                    dX_array,    lddx, 
                    dinvA_array, dinvA_length,
                    dA_displ, dB_displ, 
                    dX_displ, dinvA_displ,
                    resetozero, batchCount, queue );
    // copy X to B
    magma_tally2_cdisplace_pointers(dX_displ,    dX_array, lddx, 0, 0, batchCount, queue);
    magma_tally2_cdisplace_pointers(dB_displ,    dB_array, lddb, 0, 0, batchCount, queue);
    magma_tally2blas_clacpy_batched( Magma_tally2Full, m, n, dX_displ, lddx, dB_displ, lddb, batchCount, queue );
}

/**
    @see magma_tally2blas_ctrsm_work
    @ingroup magma_tally2_cblas3
    ********************************************************************/
extern "C"
void magma_tally2blas_ctrsm_batched(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t transA, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex** dA_array,    magma_tally2_int_t ldda,
    magma_tally2FloatComplex** dB_array,    magma_tally2_int_t lddb,
    magma_tally2_int_t batchCount, magma_tally2_queue_t queue)
{
    magma_tally2_int_t nrowA = (side == Magma_tally2Left ? m : n);

    magma_tally2_int_t info = 0;
    if ( side != Magma_tally2Left && side != Magma_tally2Right ) {
        info = -1;
    } else if ( uplo != Magma_tally2Upper && uplo != Magma_tally2Lower ) {
        info = -2;
    } else if ( transA != Magma_tally2NoTrans && transA != Magma_tally2Trans && transA != Magma_tally2ConjTrans ) {
        info = -3;
    } else if ( diag != Magma_tally2Unit && diag != Magma_tally2NonUnit ) {
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
        magma_tally2_xerbla( __func__, -(info) );
        return;
    }

    magma_tally2FloatComplex **dA_displ     = NULL;
    magma_tally2FloatComplex **dB_displ     = NULL;
    magma_tally2FloatComplex **dX_displ     = NULL;
    magma_tally2FloatComplex **dinvA_displ  = NULL;
    magma_tally2FloatComplex **dX_array     = NULL;
    magma_tally2FloatComplex **dinvA_array  = NULL;

    magma_tally2_malloc((void**)&dA_displ,  batchCount * sizeof(*dA_displ));
    magma_tally2_malloc((void**)&dB_displ,  batchCount * sizeof(*dB_displ));
    magma_tally2_malloc((void**)&dX_displ,  batchCount * sizeof(*dX_displ));
    magma_tally2_malloc((void**)&dinvA_displ,  batchCount * sizeof(*dinvA_displ));
    magma_tally2_malloc((void**)&dinvA_array, batchCount * sizeof(*dinvA_array));
    magma_tally2_malloc((void**)&dX_array, batchCount * sizeof(*dX_array));





    magma_tally2_int_t size_dinvA;
    magma_tally2_int_t lddx = m;
    magma_tally2_int_t size_x = lddx*n;

    if ( side == Magma_tally2Left ) {
        size_dinvA = ((m+TRI_NB-1)/TRI_NB)*TRI_NB*TRI_NB;
    }
    else {
        size_dinvA = ((n+TRI_NB-1)/TRI_NB)*TRI_NB*TRI_NB;
    }
    magma_tally2FloatComplex *dinvA=NULL, *dX=NULL;
    magma_tally2_int_t resetozero = 0;
    magma_tally2_cmalloc( &dinvA, size_dinvA*batchCount );
    magma_tally2_cmalloc( &dX, size_x*batchCount );
    if ( dinvA == NULL || dX == NULL ) {
        info = MAGMA_tally2_ERR_DEVICE_ALLOC;
        magma_tally2_xerbla( __func__, -(info) );
        return;
    }
    magma_tally2blas_claset_q(Magma_tally2Full, size_dinvA, batchCount, MAGMA_tally2_C_ZERO, MAGMA_tally2_C_ZERO, dinvA, size_dinvA, queue);
    magma_tally2blas_claset_q(Magma_tally2Full, lddx, n*batchCount, MAGMA_tally2_C_ZERO, MAGMA_tally2_C_ZERO, dX, lddx, queue);

    cset_pointer(dX_array, dX, lddx, 0, 0, size_x, batchCount, queue);
    cset_pointer(dinvA_array, dinvA, TRI_NB, 0, 0, size_dinvA, batchCount, queue);

    magma_tally2blas_ctrsm_work_batched( 
                    side, uplo, transA, diag, 1, 
                    m, n, alpha,
                    dA_array,    ldda,
                    dB_array,    lddb,
                    dX_array,    lddx, 
                    dinvA_array, size_dinvA,
                    dA_displ, dB_displ, 
                    dX_displ, dinvA_displ,
                    resetozero, batchCount, queue );


    magma_tally2_free( dinvA );
    magma_tally2_free( dX );
    magma_tally2_free(dA_displ);
    magma_tally2_free(dB_displ);
    magma_tally2_free(dX_displ);
    magma_tally2_free(dinvA_displ);
    magma_tally2_free(dinvA_array);
    magma_tally2_free(dX_array);


}

