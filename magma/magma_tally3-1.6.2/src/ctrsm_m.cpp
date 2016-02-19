/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca

       @generated from ztrsm_m.cpp normal z -> c, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_tally3.h"

extern "C" magma_tally3_int_t
magma_tally3_get_ctrsm_m_nb() { return 128; }

/**
    Purpose
    -------
    CTRSM solves one of the matrix equations
       op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
    where alpha is a scalar, X and B are m by n matrices, A is a unit, or
    non-unit, upper or lower triangular matrix and op( A ) is one of

       op( A ) = A      or
       op( A ) = A**T   or
       op( A ) = A**H.

    The matrix X is overwritten on B.

    Arguments
    ---------
    @param[in]
    ngpu    INTEGER
            Number of GPUs to use. ngpu > 0.

    @param[in]
    side    magma_tally3_side_t.
            On entry, SIDE specifies whether op( A ) appears on the left
            or right of X as follows:
         -     = Magma_tally3Left:        op( A )*X = alpha*B.
         -     = Magma_tally3Right:       X*op( A ) = alpha*B.

    @param[in]
    uplo    magma_tally3_uplo_t.
            On entry, UPLO specifies whether the matrix A is an upper or
            lower triangular matrix as follows:
         -     = Magma_tally3Upper:   A is an upper triangular matrix.
         -     = Magma_tally3Lower:   A is a lower triangular matrix.

    @param[in]
    transa  magma_tally3_trans_t.
            On entry, TRANSA specifies the form of op( A ) to be used in
            the matrix multiplication as follows:
         -     = Magma_tally3NoTrans:     op( A ) = A.
         -     = Magma_tally3Trans:       op( A ) = A**T.
         -     = Magma_tally3ConjTrans:   op( A ) = A**H.

    @param[in]
    diag    magma_tally3_diag_t.
            On entry, DIAG specifies whether or not A is unit triangular
            as follows:
         -     = Magma_tally3Unit:      A is assumed to be unit triangular.
         -     = Magma_tally3NonUnit:   A is not assumed to be unit triangular.

    @param[in]
    m       INTEGER.
            On entry, M specifies the number of rows of B. M must be at
            least zero.

    @param[in]
    n       INTEGER.
            On entry, N specifies the number of columns of B. N must be
            at least zero.

    @param[in]
    alpha   COMPLEX.
            On entry, ALPHA specifies the scalar alpha. When alpha is
            zero then A is not referenced and B need not be set before
            entry.

    @param[in]
    A       COMPLEX array of DIMENSION ( LDA, k ), where k is m
            when SIDE = Magma_tally3Left and is n when SIDE = Magma_tally3Right.
            Before entry with UPLO = Magma_tally3Upper, the leading k by k
            upper triangular part of the array A must contain the upper
            triangular matrix and the strictly lower triangular part of
            A is not referenced.
            Before entry with UPLO = Magma_tally3Lower, the leading k by k
            lower triangular part of the array A must contain the lower
            triangular matrix and the strictly upper triangular part of
            A is not referenced.
            Note that when DIAG = Magma_tally3Unit, the diagonal elements of
            A are not referenced either, but are assumed to be unity.

    @param[in]
    lda     INTEGER.
            On entry, LDA specifies the first dimension of A as declared
            in the calling (sub) program.
            When SIDE = Magma_tally3Left  then LDA >= max( 1, m ),
            when SIDE = Magma_tally3Right then LDA >= max( 1, n ).

    @param[in,out]
    B       COMPLEX array of DIMENSION ( LDB, n ).
            Before entry, the leading m by n part of the array B must
            contain the right-hand side matrix B, and on exit is
            overwritten by the solution matrix X.

    @param[in]
    ldb     INTEGER.
            On entry, LDB specifies the first dimension of B as declared
            in the calling (sub) program. LDB must be at least
            max( 1, m ).

    @ingroup magma_tally3_cblas3
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_ctrsm_m(
    magma_tally3_int_t ngpu,
    magma_tally3_side_t side, magma_tally3_uplo_t uplo,
    magma_tally3_trans_t transa, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *B, magma_tally3_int_t ldb)
{
#define A(i, j) (A + (j)*nb*lda + (i)*nb)
#define B(i, j) (B + (j)*nb*ldb + (i)*nb)

#define dB(gpui, i, j) (dw[gpui] + (j)*nb*lddb + (i)*nb)

#define dA(gpui, i, j) (dw[gpui] + dimb*lddb + (i)*nb + (j)*nb*ldda)

    magma_tally3FloatComplex  c_one     = MAGMA_tally3_C_ONE;
    magma_tally3FloatComplex  c_neg_one = MAGMA_tally3_C_NEG_ONE;
    magma_tally3FloatComplex  alpha_;
    magma_tally3FloatComplex* dw[Magma_tally3MaxGPUs];
    magma_tally3_queue_t stream [Magma_tally3MaxGPUs][3];
    magma_tally3_int_t lside;
    magma_tally3_int_t upper;
    magma_tally3_int_t notransp;
    magma_tally3_int_t nrowa;
    magma_tally3_int_t nb = magma_tally3_get_ctrsm_m_nb();
    magma_tally3_int_t igpu = 0;
    magma_tally3_int_t info;
    magma_tally3_int_t k, j, kb, jb;
    magma_tally3_int_t ldda, dima, lddb, dimb;
    
    magma_tally3_device_t orig_dev;
    magma_tally3_getdevice( &orig_dev );
    magma_tally3_queue_t orig_stream;
    magma_tally3blasGetKernelStream( &orig_stream );

    lside = (side == Magma_tally3Left);
    if (lside) {
        nrowa = m;
    } else {
        nrowa = n;
    }
    upper = (uplo == Magma_tally3Upper);
    notransp = (transa == Magma_tally3NoTrans);

    info = 0;
    if (! lside && side != Magma_tally3Right) {
        info = 1;
    } else if (! upper && uplo != Magma_tally3Lower) {
        info = 2;
    } else if (! notransp && transa != Magma_tally3Trans
               && transa != Magma_tally3ConjTrans) {
        info = 3;
    } else if (diag != Magma_tally3Unit && diag != Magma_tally3NonUnit) {
        info = 4;
    } else if (m < 0) {
        info = 5;
    } else if (n < 0) {
        info = 6;
    } else if (lda < max(1,nrowa)) {
        info = 9;
    } else if (ldb < max(1,m)) {
        info = 11;
    }

    // alpha = 0 case not done
    if (MAGMA_tally3_C_REAL(alpha) == 0. && MAGMA_tally3_C_IMAG(alpha) == 0.) {
        info = MAGMA_tally3_ERR_NOT_IMPLEMENTED;
    }

    if (info != 0) {
        magma_tally3_xerbla( __func__, -(info) );
        return info;
    }

    // Quick return if possible.
    if (n == 0) {
        return info;
    }

    magma_tally3_int_t nbl = (n-1)/nb+1; // number of blocks in a row
    magma_tally3_int_t mbl = (m-1)/nb+1; // number of blocks in a column

    if (lside) {
        lddb = m;
        dimb = ((nbl-1)/ngpu+1)*nb;
        if ( notransp ) {
            ldda = m;
            dima = 2 * nb;
        } else {
            ldda = 2 * nb;
            dima = m;
        }
    } else {
        lddb = ((mbl-1)/ngpu+1)*nb;
        dimb = n;
        if ( !notransp ) {
            ldda = n;
            dima = 2 * nb;
        } else {
            ldda = 2 * nb;
            dima = n;
        }
    }

    for (igpu = 0; igpu < ngpu; ++igpu) {
        magma_tally3_setdevice(igpu);
        if (MAGMA_tally3_SUCCESS != magma_tally3_cmalloc( &dw[igpu], (dimb*lddb + dima*ldda) )) {
            info = MAGMA_tally3_ERR_DEVICE_ALLOC;
            return info;
        }
        magma_tally3_queue_create( &stream[igpu][0] );
        magma_tally3_queue_create( &stream[igpu][1] );
        magma_tally3_queue_create( &stream[igpu][2] );
    }

    if (lside) {
        if (notransp) {
            // Form  B := alpha*inv( A )*B
            if (upper) {
                // left upper notranspose
                magma_tally3_int_t nloc[Magma_tally3MaxGPUs];
                for (igpu = 0; igpu < ngpu; ++igpu)
                    nloc[igpu] = 0;

                // copy B to mgpus
                for (k = 0; k < nbl; ++k) {
                    igpu = k%ngpu;
                    magma_tally3_setdevice(igpu);
                    kb = min(nb, n-k*nb);
                    nloc[igpu] += kb;
                    magma_tally3_csetmatrix_async( m, kb,
                                            B(0, k),              ldb,
                                            dB(igpu, 0, k/ngpu), lddb, stream[igpu][(mbl+1)%2] );
                }
                jb = min(nb, m-(mbl-1)*nb);
                for (igpu = 0; igpu < ngpu; ++igpu) {
                    magma_tally3_setdevice(igpu);
                    magma_tally3_csetmatrix_async( m, jb,
                                            A(0, mbl-1),            lda,
                                            dA(igpu, 0, (mbl-1)%2), ldda, stream[igpu][(mbl+1)%2] );
                }
                for (j = mbl-1; j >= 0; --j) {
                    if (j > 0) {
                        jb = nb;
                        for (igpu = 0; igpu < ngpu; ++igpu) {
                            magma_tally3_setdevice(igpu);
                            magma_tally3_csetmatrix_async( j*nb, jb,
                                                    A(0, j-1),            lda,
                                                    dA(igpu, 0, (j+1)%2), ldda, stream[igpu][(j+1)%2] );
                        }
                    }
                    if (j == mbl-1)
                        alpha_=alpha;
                    else
                        alpha_= c_one;

                    jb = min(nb, m-j*nb);

                    for (igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally3_setdevice(igpu);
                        magma_tally3blasSetKernelStream(stream[igpu][j%2]);
                        magma_tally3_ctrsm(side, uplo, transa, diag, jb, nloc[igpu], alpha_, dA(igpu, j, j%2), ldda,
                                    dB(igpu, j, 0), lddb );
                    }

                    if (j > 0) {
                        for (igpu = 0; igpu < ngpu; ++igpu) {
                            magma_tally3_setdevice(igpu);
                            magma_tally3blasSetKernelStream(stream[igpu][j%2]);
                            magma_tally3_cgemm(transa, Magma_tally3NoTrans, j*nb, nloc[igpu], jb, c_neg_one, dA(igpu, 0, j%2), ldda,
                                        dB(igpu, j, 0), lddb, alpha_, dB(igpu, 0, 0), lddb );
                        }
                    }

                    for (igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally3_queue_sync( stream[igpu][j%2] );
                    }

                    for (k = 0; k < nbl; ++k) {
                        igpu = k%ngpu;
                        magma_tally3_setdevice(igpu);
                        kb = min(nb, n-k*nb);
                        magma_tally3_cgetmatrix_async( jb, kb,
                                                dB(igpu, j, k/ngpu), lddb,
                                                B(j, k),              ldb, stream[igpu][2] );
                    }
                }
            }
            else {
                // left lower notranspose
                magma_tally3_int_t nloc[Magma_tally3MaxGPUs];
                for (igpu = 0; igpu < ngpu; ++igpu)
                    nloc[igpu] = 0;

                // copy B to mgpus
                for (k = 0; k < nbl; ++k) {
                    igpu = k%ngpu;
                    magma_tally3_setdevice(igpu);
                    kb = min(nb, n-k*nb);
                    nloc[igpu] += kb;
                    magma_tally3_csetmatrix_async( m, kb,
                                            B(0, k),              ldb,
                                            dB(igpu, 0, k/ngpu), lddb, stream[igpu][0] );
                }
                jb = min(nb, m);
                for (igpu = 0; igpu < ngpu; ++igpu) {
                    magma_tally3_setdevice(igpu);
                    magma_tally3_csetmatrix_async( m, jb,
                                            A(0, 0),        lda,
                                            dA(igpu, 0, 0), ldda, stream[igpu][0] );
                }
                for (j = 0; j < mbl; ++j) {
                    if ((j+1)*nb < m) {
                        jb = min(nb, m-(j+1)*nb);
                        for (igpu = 0; igpu < ngpu; ++igpu) {
                            magma_tally3_setdevice(igpu);
                            magma_tally3_csetmatrix_async( (m-(j+1)*nb), jb,
                                                    A(j+1, j+1),            lda,
                                                    dA(igpu, j+1, (j+1)%2), ldda, stream[igpu][(j+1)%2] );
                        }
                    }
                    jb = min(nb, m-j*nb);

                    if (j == 0)
                        alpha_=alpha;
                    else
                        alpha_= c_one;

                    for (igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally3_setdevice(igpu);
                        magma_tally3blasSetKernelStream(stream[igpu][j%2]);
                        magma_tally3_ctrsm(side, uplo, transa, diag, jb, nloc[igpu], alpha_, dA(igpu, j, j%2), ldda,
                                    dB(igpu, j, 0), lddb );
                    }

                    if ( j < mbl-1 ) {
                        for (igpu = 0; igpu < ngpu; ++igpu) {
                            magma_tally3_setdevice(igpu);
                            magma_tally3blasSetKernelStream(stream[igpu][j%2]);
                            magma_tally3_cgemm(transa, Magma_tally3NoTrans, m-(j+1)*nb, nloc[igpu], nb, c_neg_one, dA(igpu, j+1, j%2), ldda,
                                        dB(igpu, j, 0), lddb, alpha_, dB(igpu, j+1, 0), lddb );
                        }
                    }

                    for (igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally3_queue_sync( stream[igpu][j%2] );
                    }

                    for (k = 0; k < nbl; ++k) {
                        igpu = k%ngpu;
                        magma_tally3_setdevice(igpu);
                        kb = min(nb, n-k*nb);
                        magma_tally3_cgetmatrix_async( jb, kb,
                                                dB(igpu, j, k/ngpu), lddb,
                                                B(j, k),              ldb, stream[igpu][2] );
                    }
                }
            }
        }
        else {
            // Form  B := alpha*inv( A**T )*B
            if (upper) {
                // left upper transpose or conj transpose
                magma_tally3_int_t nloc[Magma_tally3MaxGPUs];
                for (igpu = 0; igpu < ngpu; ++igpu)
                    nloc[igpu] = 0;

                // copy B to mgpus
                for (k = 0; k < nbl; ++k) {
                    igpu = k%ngpu;
                    magma_tally3_setdevice(igpu);
                    kb = min(nb, n-k*nb);
                    nloc[igpu] += kb;
                    magma_tally3_csetmatrix_async( m, kb,
                                            B(0, k),              ldb,
                                            dB(igpu, 0, k/ngpu), lddb, stream[igpu][0] );
                }
                jb = min(nb, m);
                for (igpu = 0; igpu < ngpu; ++igpu) {
                    magma_tally3_setdevice(igpu);
                    magma_tally3_csetmatrix_async( jb, m,
                                            A(0, 0),        lda,
                                            dA(igpu, 0, 0), ldda, stream[igpu][0] );
                }
                for (j = 0; j < mbl; ++j) {
                    if ((j+1)*nb < m) {
                        jb = min(nb, m-(j+1)*nb);
                        for (igpu = 0; igpu < ngpu; ++igpu) {
                            magma_tally3_setdevice(igpu);
                            magma_tally3_csetmatrix_async( jb, m-(j+1)*nb,
                                                    A(j+1, j+1),            lda,
                                                    dA(igpu, (j+1)%2, j+1), ldda, stream[igpu][(j+1)%2] );
                        }
                    }
                    jb = min(nb, m-j*nb);

                    if (j == 0)
                        alpha_=alpha;
                    else
                        alpha_= c_one;

                    for (igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally3_setdevice(igpu);
                        magma_tally3blasSetKernelStream(stream[igpu][j%2]);
                        magma_tally3_ctrsm(side, uplo, transa, diag, jb, nloc[igpu], alpha_, dA(igpu, j%2, j), ldda,
                                    dB(igpu, j, 0), lddb );
                    }

                    if ( j < mbl-1 ) {
                        for (igpu = 0; igpu < ngpu; ++igpu) {
                            magma_tally3_setdevice(igpu);
                            magma_tally3blasSetKernelStream(stream[igpu][j%2]);
                            magma_tally3_cgemm(transa, Magma_tally3NoTrans, m-(j+1)*nb, nloc[igpu], nb, c_neg_one, dA(igpu, j%2, j+1), ldda,
                                        dB(igpu, j, 0), lddb, alpha_, dB(igpu, j+1, 0), lddb );
                        }
                    }

                    for (igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally3_queue_sync( stream[igpu][j%2] );
                    }

                    for (k = 0; k < nbl; ++k) {
                        igpu = k%ngpu;
                        magma_tally3_setdevice(igpu);
                        kb = min(nb, n-k*nb);
                        magma_tally3_cgetmatrix_async( jb, kb,
                                                dB(igpu, j, k/ngpu), lddb,
                                                B(j, k),              ldb, stream[igpu][2] );
                    }
                }
            }
            else {
                // left lower transpose or conj transpose
                magma_tally3_int_t nloc[Magma_tally3MaxGPUs];
                for (igpu = 0; igpu < ngpu; ++igpu)
                    nloc[igpu] = 0;

                // copy B to mgpus
                for (k = 0; k < nbl; ++k) {
                    igpu = k%ngpu;
                    magma_tally3_setdevice(igpu);
                    kb = min(nb, n-k*nb);
                    nloc[igpu] += kb;
                    magma_tally3_csetmatrix_async( m, kb,
                                            B(0, k),              ldb,
                                            dB(igpu, 0, k/ngpu), lddb, stream[igpu][(mbl+1)%2] );
                }
                jb = min(nb, m-(mbl-1)*nb);
                for (igpu = 0; igpu < ngpu; ++igpu) {
                    magma_tally3_setdevice(igpu);
                    magma_tally3_csetmatrix_async( jb, m,
                                            A(mbl-1, 0),            lda,
                                            dA(igpu, (mbl-1)%2, 0), ldda, stream[igpu][(mbl+1)%2] );
                }
                for (j = mbl-1; j >= 0; --j) {
                    if (j > 0) {
                        jb = nb;
                        for (igpu = 0; igpu < ngpu; ++igpu) {
                            magma_tally3_setdevice(igpu);
                            magma_tally3_csetmatrix_async( jb, j*nb,
                                                    A(j-1, 0),            lda,
                                                    dA(igpu, (j+1)%2, 0), ldda, stream[igpu][(j+1)%2] );
                        }
                    }
                    if (j == mbl-1)
                        alpha_=alpha;
                    else
                        alpha_= c_one;

                    jb = min(nb, m-j*nb);

                    for (igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally3_setdevice(igpu);
                        magma_tally3blasSetKernelStream(stream[igpu][j%2]);
                        magma_tally3_ctrsm(side, uplo, transa, diag, jb, nloc[igpu], alpha_, dA(igpu, j%2, j), ldda,
                                    dB(igpu, j, 0), lddb );
                    }

                    if (j > 0) {
                        for (igpu = 0; igpu < ngpu; ++igpu) {
                            magma_tally3_setdevice(igpu);
                            magma_tally3blasSetKernelStream(stream[igpu][j%2]);
                            magma_tally3_cgemm(transa, Magma_tally3NoTrans, j*nb, nloc[igpu], jb, c_neg_one, dA(igpu, j%2, 0), ldda,
                                        dB(igpu, j, 0), lddb, alpha_, dB(igpu, 0, 0), lddb );
                        }
                    }

                    for (igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally3_queue_sync( stream[igpu][j%2] );
                    }

                    for (k = 0; k < nbl; ++k) {
                        igpu = k%ngpu;
                        magma_tally3_setdevice(igpu);
                        kb = min(nb, n-k*nb);
                        magma_tally3_cgetmatrix_async( jb, kb,
                                                dB(igpu, j, k/ngpu), lddb,
                                                B(j, k),              ldb, stream[igpu][2] );
                    }
                }
            }
        }
    }
    else {
        if (notransp) {
            // Form  B := alpha*B*inv( A ).
            if (upper) {
                // right upper notranspose
                magma_tally3_int_t mloc[Magma_tally3MaxGPUs];
                for (igpu = 0; igpu < ngpu; ++igpu)
                    mloc[igpu] = 0;

                // copy B to mgpus
                for (j = 0; j < mbl; ++j) {
                    igpu = j%ngpu;
                    magma_tally3_setdevice(igpu);
                    jb = min(nb, m-j*nb);
                    mloc[igpu] += jb;
                    magma_tally3_csetmatrix_async( jb, n,
                                            B(j, 0),              ldb,
                                            dB(igpu, j/ngpu, 0), lddb, stream[igpu][0] );
                }
                kb = min(nb, n);
                for (igpu = 0; igpu < ngpu; ++igpu) {
                    magma_tally3_setdevice(igpu);
                    magma_tally3_csetmatrix_async( kb, n,
                                            A(0, 0),        lda,
                                            dA(igpu, 0, 0), ldda, stream[igpu][0] );
                }
                for (k = 0; k < nbl; ++k) {
                    if ((k+1)*nb < n) {
                        kb = min(nb, n-(k+1)*nb);
                        for (igpu = 0; igpu < ngpu; ++igpu) {
                            magma_tally3_setdevice(igpu);
                            magma_tally3_csetmatrix_async( kb, n-(k+1)*nb,
                                                    A(k+1, k+1),            lda,
                                                    dA(igpu, (k+1)%2, k+1), ldda, stream[igpu][(k+1)%2] );
                        }
                    }
                    kb = min(nb, n-k*nb);

                    if (k == 0)
                        alpha_=alpha;
                    else
                        alpha_= c_one;

                    for (igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally3_setdevice(igpu);
                        magma_tally3blasSetKernelStream(stream[igpu][k%2]);
                        magma_tally3_ctrsm(side, uplo, transa, diag, mloc[igpu], kb, alpha_, dA(igpu, k%2, k), ldda,
                                    dB(igpu, 0, k), lddb );
                    }

                    if ( k < nbl-1 ) {
                        for (igpu = 0; igpu < ngpu; ++igpu) {
                            magma_tally3_setdevice(igpu);
                            magma_tally3blasSetKernelStream(stream[igpu][k%2]);
                            magma_tally3_cgemm(Magma_tally3NoTrans, transa, mloc[igpu], n-(k+1)*nb, nb, c_neg_one, dB(igpu, 0, k), lddb,
                                        dA(igpu, k%2, k+1), ldda, alpha_, dB(igpu, 0, k+1), lddb );
                        }
                    }

                    for (igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally3_queue_sync( stream[igpu][k%2] );
                    }

                    for (j = 0; j < mbl; ++j) {
                        igpu = j%ngpu;
                        magma_tally3_setdevice(igpu);
                        jb = min(nb, m-j*nb);
                        magma_tally3_cgetmatrix_async( jb, kb,
                                                dB(igpu, j/ngpu, k), lddb,
                                                B(j, k),              ldb, stream[igpu][2] );
                    }
                }
            }
            else {
                // right lower notranspose
                magma_tally3_int_t mloc[Magma_tally3MaxGPUs];
                for (igpu = 0; igpu < ngpu; ++igpu)
                    mloc[igpu] = 0;

                // copy B to mgpus
                for (j = 0; j < mbl; ++j) {
                    igpu = j%ngpu;
                    magma_tally3_setdevice(igpu);
                    jb = min(nb, m-j*nb);
                    mloc[igpu] += jb;
                    magma_tally3_csetmatrix_async( jb, n,
                                            B(j, 0),              ldb,
                                            dB(igpu, j/ngpu, 0), lddb, stream[igpu][(nbl+1)%2] );
                }
                kb = min(nb, n-(nbl-1)*nb);
                for (igpu = 0; igpu < ngpu; ++igpu) {
                    magma_tally3_setdevice(igpu);
                    magma_tally3_csetmatrix_async( kb, n,
                                            A(nbl-1, 0),            lda,
                                            dA(igpu, (nbl-1)%2, 0), ldda, stream[igpu][(nbl+1)%2] );
                }
                for (k = nbl-1; k >= 0; --k) {
                    if (k > 0) {
                        kb = nb;
                        for (igpu = 0; igpu < ngpu; ++igpu) {
                            magma_tally3_setdevice(igpu);
                            magma_tally3_csetmatrix_async( kb, k*nb,
                                                    A(k-1, 0),            lda,
                                                    dA(igpu, (k+1)%2, 0), ldda, stream[igpu][(k+1)%2] );
                        }
                    }
                    if (k == nbl-1)
                        alpha_=alpha;
                    else
                        alpha_= c_one;

                    kb = min(nb, n-k*nb);

                    for (igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally3_setdevice(igpu);
                        magma_tally3blasSetKernelStream(stream[igpu][k%2]);
                        magma_tally3_ctrsm(side, uplo, transa, diag, mloc[igpu], kb, alpha_, dA(igpu, k%2, k), ldda,
                                    dB(igpu, 0, k), lddb );
                    }

                    if (k > 0) {
                        for (igpu = 0; igpu < ngpu; ++igpu) {
                            magma_tally3_setdevice(igpu);
                            magma_tally3blasSetKernelStream(stream[igpu][k%2]);
                            magma_tally3_cgemm(Magma_tally3NoTrans, transa, mloc[igpu], k*nb, kb, c_neg_one, dB(igpu, 0, k), lddb,
                                        dA(igpu, k%2, 0), ldda, alpha_, dB(igpu, 0, 0), lddb );
                        }
                    }

                    for (igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally3_queue_sync( stream[igpu][k%2] );
                    }

                    for (j = 0; j < mbl; ++j) {
                        igpu = j%ngpu;
                        magma_tally3_setdevice(igpu);
                        jb = min(nb, m-j*nb);
                        magma_tally3_cgetmatrix_async( jb, kb,
                                                dB(igpu, j/ngpu, k), lddb,
                                                B(j, k),              ldb, stream[igpu][2] );
                    }
                }
            }
        }
        else {
            // Form  B := alpha*B*inv( A**T ).
            if (upper) {
                // right upper transpose or conj transpose
                magma_tally3_int_t mloc[Magma_tally3MaxGPUs];
                for (igpu = 0; igpu < ngpu; ++igpu)
                    mloc[igpu] = 0;

                // copy B to mgpus
                for (j = 0; j < mbl; ++j) {
                    igpu = j%ngpu;
                    magma_tally3_setdevice(igpu);
                    jb = min(nb, m-j*nb);
                    mloc[igpu] += jb;
                    magma_tally3_csetmatrix_async( jb, n,
                                            B(j, 0),              ldb,
                                            dB(igpu, j/ngpu, 0), lddb, stream[igpu][(nbl+1)%2] );
                }
                kb = min(nb, n-(nbl-1)*nb);
                for (igpu = 0; igpu < ngpu; ++igpu) {
                    magma_tally3_setdevice(igpu);
                    magma_tally3_csetmatrix_async( n, kb,
                                            A(0, nbl-1),            lda,
                                            dA(igpu, 0, (nbl-1)%2), ldda, stream[igpu][(nbl+1)%2] );
                }
                for (k = nbl-1; k >= 0; --k) {
                    if (k > 0) {
                        kb = nb;
                        for (igpu = 0; igpu < ngpu; ++igpu) {
                            magma_tally3_setdevice(igpu);
                            magma_tally3_csetmatrix_async( k*nb, kb,
                                                    A(0, k-1),            lda,
                                                    dA(igpu, 0, (k+1)%2), ldda, stream[igpu][(k+1)%2] );
                        }
                    }
                    if (k == nbl-1)
                        alpha_=alpha;
                    else
                        alpha_= c_one;

                    kb = min(nb, n-k*nb);

                    for (igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally3_setdevice(igpu);
                        magma_tally3blasSetKernelStream(stream[igpu][k%2]);
                        magma_tally3_ctrsm(side, uplo, transa, diag, mloc[igpu], kb, alpha_, dA(igpu, k, k%2), ldda,
                                    dB(igpu, 0, k), lddb );
                    }

                    if (k > 0) {
                        for (igpu = 0; igpu < ngpu; ++igpu) {
                            magma_tally3_setdevice(igpu);
                            magma_tally3blasSetKernelStream(stream[igpu][k%2]);
                            magma_tally3_cgemm(Magma_tally3NoTrans, transa, mloc[igpu], k*nb, kb, c_neg_one, dB(igpu, 0, k), lddb,
                                        dA(igpu, 0, k%2), ldda, alpha_, dB(igpu, 0, 0), lddb );
                        }
                    }

                    for (igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally3_queue_sync( stream[igpu][k%2] );
                    }

                    for (j = 0; j < mbl; ++j) {
                        igpu = j%ngpu;
                        magma_tally3_setdevice(igpu);
                        jb = min(nb, m-j*nb);
                        magma_tally3_cgetmatrix_async( jb, kb,
                                                dB(igpu, j/ngpu, k), lddb,
                                                B(j, k),              ldb, stream[igpu][2] );
                    }
                }
            }
            else {
                // right lower transpose or conj transpose
                magma_tally3_int_t mloc[Magma_tally3MaxGPUs];
                for (igpu = 0; igpu < ngpu; ++igpu)
                    mloc[igpu] = 0;

                // copy B to mgpus
                for (j = 0; j < mbl; ++j) {
                    igpu = j%ngpu;
                    magma_tally3_setdevice(igpu);
                    jb = min(nb, m-j*nb);
                    mloc[igpu] += jb;
                    magma_tally3_csetmatrix_async( jb, n,
                                            B(j, 0),              ldb,
                                            dB(igpu, j/ngpu, 0), lddb, stream[igpu][0] );
                }
                kb = min(nb, n);
                for (igpu = 0; igpu < ngpu; ++igpu) {
                    magma_tally3_setdevice(igpu);
                    magma_tally3_csetmatrix_async( n, kb,
                                            A(0, 0),        lda,
                                            dA(igpu, 0, 0), ldda, stream[igpu][0] );
                }
                for (k = 0; k < nbl; ++k) {
                    if ((k+1)*nb < n) {
                        kb = min(nb, n-(k+1)*nb);
                        for (igpu = 0; igpu < ngpu; ++igpu) {
                            magma_tally3_setdevice(igpu);
                            magma_tally3_csetmatrix_async( (n-(k+1)*nb), kb,
                                                    A(k+1, k+1),            lda,
                                                    dA(igpu, k+1, (k+1)%2), ldda, stream[igpu][(k+1)%2] );
                        }
                    }
                    kb = min(nb, n-k*nb);

                    if (k == 0)
                        alpha_=alpha;
                    else
                        alpha_= c_one;

                    for (igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally3_setdevice(igpu);
                        magma_tally3blasSetKernelStream(stream[igpu][k%2]);
                        magma_tally3_ctrsm(side, uplo, transa, diag, mloc[igpu], kb, alpha_, dA(igpu, k, k%2), ldda,
                                    dB(igpu, 0, k), lddb );
                    }

                    if ( k < nbl-1 ) {
                        for (igpu = 0; igpu < ngpu; ++igpu) {
                            magma_tally3_setdevice(igpu);
                            magma_tally3blasSetKernelStream(stream[igpu][k%2]);
                            magma_tally3_cgemm(Magma_tally3NoTrans, transa, mloc[igpu], n-(k+1)*nb, nb, c_neg_one, dB(igpu, 0, k), lddb,
                                        dA(igpu, k+1, k%2), ldda, alpha_, dB(igpu, 0, k+1), lddb );
                        }
                    }

                    for (igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally3_queue_sync( stream[igpu][k%2] );
                    }

                    for (j = 0; j < mbl; ++j) {
                        igpu = j%ngpu;
                        magma_tally3_setdevice(igpu);
                        jb = min(nb, m-j*nb);
                        magma_tally3_cgetmatrix_async( jb, kb,
                                                dB(igpu, j/ngpu, k), lddb,
                                                B(j, k),              ldb, stream[igpu][2] );
                    }
                }
            }
        }
    }


    for (igpu = 0; igpu < ngpu; ++igpu) {
        magma_tally3_setdevice(igpu);
        magma_tally3_queue_sync( stream[igpu][2] );
        magma_tally3_queue_destroy( stream[igpu][0] );
        magma_tally3_queue_destroy( stream[igpu][1] );
        magma_tally3_queue_destroy( stream[igpu][2] );
        magma_tally3_free( dw[igpu] );
    }
    magma_tally3_setdevice( orig_dev );
    magma_tally3blasSetKernelStream( orig_stream );

    return info;
} /* magma_tally3_ctrsm_m */
