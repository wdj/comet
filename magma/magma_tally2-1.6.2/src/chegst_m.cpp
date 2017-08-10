/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Azzam Haidar

       @generated from zhegst_m.cpp normal z -> c, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_tally2.h"


static void magma_tally2_chegst_m_1_L_col_update(magma_tally2_int_t nk, magma_tally2_int_t nb, magma_tally2FloatComplex* dA_col, magma_tally2_int_t ldda,
                                          magma_tally2FloatComplex* dC1, magma_tally2_int_t lddc1, magma_tally2FloatComplex* dC2, magma_tally2_int_t lddc2);

static void magma_tally2_chegst_m_1_U_row_update(magma_tally2_int_t nk, magma_tally2_int_t nb, magma_tally2FloatComplex* dA_row, magma_tally2_int_t ldda,
                                          magma_tally2FloatComplex* dC1, magma_tally2_int_t lddc1, magma_tally2FloatComplex* dC2, magma_tally2_int_t lddc2);

/**
    Purpose
    -------
    CHEGST_M reduces a complex Hermitian-definite generalized
    eigenproblem to standard form.
    
    If ITYPE = 1, the problem is A*x = lambda*B*x,
    and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
    
    If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
    B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.
    
    B must have been previously factorized as U**H*U or L*L**H by CPOTRF.
    
    Arguments
    ---------
    @param[in]
    ngpu    INTEGER
            Number of GPUs to use. ngpu > 0.

    @param[in]
    itype   INTEGER
            = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);
            = 2 or 3: compute U*A*U**H or L**H*A*L.
    
    @param[in]
    uplo    magma_tally2_uplo_t
      -     = Magma_tally2Upper:  Upper triangle of A is stored and B is factored as U**H*U;
      -     = Magma_tally2Lower:  Lower triangle of A is stored and B is factored as L*L**H.
    
    @param[in]
    n       INTEGER
            The order of the matrices A and B.  N >= 0.
    
    @param[in,out]
    A       COMPLEX array, dimension (LDA,N)
            On entry, the Hermitian matrix A.  If UPLO = Magma_tally2Upper, the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_tally2Lower, the
            leading N-by-N lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
    \n
            On exit, if INFO = 0, the transformed matrix, stored in the
            same format as A.
    
    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).
    
    @param[in]
    B       COMPLEX array, dimension (LDB,N)
            The triangular factor from the Cholesky factorization of B,
            as returned by CPOTRF.
    
    @param[in]
    ldb     INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).
    
    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
    

    @ingroup magma_tally2_cheev_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_chegst_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info)
{
#define A(i, j) (A + (j)*nb*lda + (i)*nb)
#define B(i, j) (B + (j)*nb*ldb + (i)*nb)
#define dA(gpui, i, j)    (dw[gpui] + (j)*nb*ldda + (i)*nb)
#define dB_c(gpui, i, j)  (dw[gpui] + dima*ldda + (i)*nb + (j)*nb*lddbc)
#define dB_r(gpui, i, j)  (dw[gpui] + dima*ldda + (i)*nb + (j)*nb*lddbr)
#define dwork(gpui, i, j) (dw[gpui] + dima*ldda + lddbc*lddbr + (j)*nb*nb + (i)*nb)

    const char* uplo_ = lapack_uplo_const_tally2( uplo );

    float             d_one      = 1.0;
    magma_tally2FloatComplex    c_one      = MAGMA_tally2_C_ONE;
    magma_tally2FloatComplex    c_neg_one  = MAGMA_tally2_C_NEG_ONE;
    magma_tally2FloatComplex    c_half     = MAGMA_tally2_C_HALF;
    magma_tally2FloatComplex    c_neg_half = MAGMA_tally2_C_NEG_HALF;
    magma_tally2FloatComplex* dw[Magma_tally2MaxGPUs];

    magma_tally2_queue_t stream [Magma_tally2MaxGPUs][3];
    magma_tally2_event_t  event  [Magma_tally2MaxGPUs][2];

    int upper = (uplo == Magma_tally2Upper);

    magma_tally2_int_t nb = magma_tally2_get_chegst_nb_m(n);

    /* Test the input parameters. */
    *info = 0;
    if (itype < 1 || itype > 3) {
        *info = -1;
    } else if (! upper && uplo != Magma_tally2Lower) {
        *info = -2;
    } else if (n < 0) {
        *info = -3;
    } else if (lda < max(1,n)) {
        *info = -5;
    } else if (ldb < max(1,n)) {
        *info = -7;
    }
    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return */
    if ( n == 0 )
        return *info;

    magma_tally2_device_t orig_dev;
    magma_tally2_getdevice( &orig_dev );
    magma_tally2_queue_t orig_stream;
    magma_tally2blasGetKernelStream( &orig_stream );

    magma_tally2_int_t nbl = (n-1)/nb+1; // number of blocks

    magma_tally2_int_t ldda = 0;
    magma_tally2_int_t dima = 0;

    if ( (itype == 1 && upper) || (itype != 1 && !upper) ) {
        ldda = ((nbl-1)/ngpu+1)*nb;
        dima = n;
    } else {
        ldda = n;
        dima = ((nbl-1)/ngpu+1)*nb;
    }
    magma_tally2_int_t lddbr = 2 * nb;
    magma_tally2_int_t lddbc = n;

    for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
        magma_tally2_setdevice(igpu);

        if (MAGMA_tally2_SUCCESS != magma_tally2_cmalloc( &dw[igpu], (dima*ldda + lddbc*lddbr + n*nb) )) {
            *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
            return *info;
        }
        magma_tally2_queue_create( &stream[igpu][0] );
        magma_tally2_queue_create( &stream[igpu][1] );
        magma_tally2_queue_create( &stream[igpu][2] );
        magma_tally2_event_create( &event[igpu][0] );
        magma_tally2_event_create( &event[igpu][1] );
    }

    /* Use hybrid blocked code */

    if (itype == 1) {
        if (upper) {
            /* Compute inv(U')*A*inv(U) */

            //copy A to mgpu
            for (magma_tally2_int_t k = 0; k < nbl; ++k) {
                magma_tally2_int_t igpu = k%ngpu;
                magma_tally2_setdevice(igpu);
                magma_tally2_int_t kb = min(nb, n-k*nb);
                magma_tally2_csetmatrix_async(kb, n-k*nb,
                                       A(k, k),              lda,
                                       dA(igpu, k/ngpu, k), ldda, stream[igpu][0] );
            }

            for (magma_tally2_int_t k = 0; k < nbl; ++k) {
                magma_tally2_int_t ind_k  =   k   % 2;
                magma_tally2_int_t ind_k1 = (k+1) % 2;
                magma_tally2_int_t kb= min(n-k*nb,nb);

                // Copy B panel
                for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                    magma_tally2_setdevice(igpu);
                    magma_tally2_queue_sync( stream[igpu][0] ); // sync previous B panel copy

                    // sync dwork copy and update (uses B panel of the next copy)
                    magma_tally2_queue_wait_event( stream[igpu][0], event[igpu][1] );
                    magma_tally2_queue_wait_event( stream[igpu][2], event[igpu][1] );
                    magma_tally2_csetmatrix_async(kb, n-k*nb,
                                           B(k, k),              ldb,
                                           dB_r(igpu, ind_k, k), lddbr, stream[igpu][0] );
                }

                magma_tally2_int_t igpu_p = k%ngpu;

                if (k > 0) {
                    // Update the next panel
                    magma_tally2_setdevice(igpu_p);
                    magma_tally2blasSetKernelStream(stream[igpu_p][2]);

                    magma_tally2_int_t nk = n - k*nb;

                    magma_tally2_cher2k(Magma_tally2Upper, Magma_tally2ConjTrans, kb, nb,
                                 c_neg_one, dwork(igpu_p, 0, k), nb, dB_r(igpu_p, ind_k1, k), lddbr,
                                 d_one, dA(igpu_p, k/ngpu, k), ldda);

                    // copy Akk block on the CPU
                    magma_tally2_cgetmatrix_async(kb, kb,
                                           dA(igpu_p, k/ngpu, k), ldda,
                                           A(k, k),                lda, stream[igpu_p][2] );

                    magma_tally2_event_record( event[igpu_p][0], stream[igpu_p][2]);

                    magma_tally2_cgemm(Magma_tally2ConjTrans, Magma_tally2NoTrans, kb, nk-kb, nb, c_neg_one, dwork(igpu_p, 0, k), nb,
                                dB_r(igpu_p, ind_k1, k+1), lddbr, c_one, dA(igpu_p, k/ngpu, k+1), ldda );

                    magma_tally2_cgemm(Magma_tally2ConjTrans, Magma_tally2NoTrans, kb, nk-kb, nb, c_neg_one, dB_r(igpu_p, ind_k1, k), lddbr,
                                dwork(igpu_p, 0, k+1), nb, c_one, dA(igpu_p, k/ngpu, k+1), ldda );

                    // Update the panels of the other GPUs
                    for (magma_tally2_int_t j=k+1; j < nbl; ++j) {
                        magma_tally2_int_t igpu = j%ngpu;
                        if (igpu != igpu_p) {
                            magma_tally2_setdevice(igpu);
                            magma_tally2blasSetKernelStream(stream[igpu][1]);

                            magma_tally2_chegst_m_1_U_row_update(n-j*nb, nb, dA(igpu, j/ngpu, j), ldda,
                                                          dwork(igpu, 0, j), nb, dB_r(igpu, ind_k1, j), lddbr); // cher2k on j-th row
                        }
                    }
                }
                // compute next panel
                magma_tally2_setdevice(igpu_p);

                if (k+1 < nbl) {
                    magma_tally2_queue_sync( stream[igpu_p][0] ); // sync B panel copy
                    magma_tally2blasSetKernelStream(stream[igpu_p][2]);

                    magma_tally2_ctrsm(Magma_tally2Left, uplo, Magma_tally2ConjTrans, Magma_tally2NonUnit,
                                kb, n-(k+1)*nb,
                                c_one, dB_r(igpu_p, ind_k, k), lddbr,
                                dA(igpu_p, k/ngpu, k+1), ldda);
                }

                magma_tally2_event_sync( event[igpu_p][0] ); // sync Akk copy
                lapackf77_chegst( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);

                if (k+1 < nbl) {
                    magma_tally2_csetmatrix_async(kb, kb,
                                           A(k, k),              lda,
                                           dA(igpu_p, k/ngpu, k), ldda, stream[igpu_p][2] );

                    magma_tally2_chemm(Magma_tally2Left, uplo,
                                kb, n-(k+1)*nb,
                                c_neg_half, dA(igpu_p, k/ngpu, k), ldda,
                                dB_r(igpu_p, ind_k, k+1), lddbr,
                                c_one, dA(igpu_p, k/ngpu, k+1), ldda);

                    magma_tally2_cgetmatrix_async(kb, n-(k+1)*nb,
                                           dA(igpu_p, k/ngpu, k+1), ldda,
                                           A(k, k+1),                lda, stream[igpu_p][2] );
                }

                if (k > 0) {
                    // Update the remaining panels of GPU igpu_p
                    for (magma_tally2_int_t j=k+ngpu; j < nbl; j += ngpu) {
                        magma_tally2_setdevice(igpu_p);
                        magma_tally2blasSetKernelStream(stream[igpu_p][1]);

                        magma_tally2_chegst_m_1_U_row_update(n-j*nb, nb, dA(igpu_p, j/ngpu, j), ldda,
                                                      dwork(igpu_p, 0, j), nb, dB_r(igpu_p, ind_k1, j), lddbr); // cher2k on j-th row
                    }
                }

                if (k+1 < nbl) {
                    // send the partially updated panel of dA to each gpu in dwork block
                    magma_tally2_setdevice(igpu_p);
                    magma_tally2_queue_sync( stream[igpu_p][2] );

                    for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally2_setdevice(igpu);
                        magma_tally2_csetmatrix_async(kb, n-(k+1)*nb,
                                               A(k, k+1),          lda,
                                               dwork(igpu, 0, k+1), nb, stream[igpu][1] );

                        magma_tally2_event_record( event[igpu][1], stream[igpu][1]);
                    }

                    magma_tally2_setdevice(igpu_p);
                    magma_tally2blasSetKernelStream(stream[igpu_p][1]);

                    magma_tally2_chemm(Magma_tally2Left, uplo,
                                kb, n-(k+1)*nb,
                                c_neg_half, dA(igpu_p, k/ngpu, k), ldda,
                                dB_r(igpu_p, ind_k, k+1), lddbr,
                                c_one, dA(igpu_p, k/ngpu, k+1), ldda);
                }
            }

            for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                magma_tally2_queue_sync( stream[igpu][1] );
            }

            if (n > nb) {
                magma_tally2_int_t nloc[Magma_tally2MaxGPUs] = { 0 };

                for (magma_tally2_int_t j = 1; j < nbl; ++j) {
                    nloc[(j-1)%ngpu] += nb;

                    magma_tally2_int_t jb = min(nb, n-j*nb);

                    for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally2_setdevice(igpu);
                        if (nloc[igpu] > 0) {
                            magma_tally2_csetmatrix_async(jb, n-j*nb,
                                                   B(j, j),            ldb,
                                                   dB_r(igpu, j%2, j), lddbr, stream[igpu][j%2] );

                            magma_tally2_queue_wait_event( stream[igpu][j%2], event[igpu][0] );

                            magma_tally2blasSetKernelStream(stream[igpu][j%2]);
                            magma_tally2_ctrsm(Magma_tally2Right, uplo, Magma_tally2NoTrans, Magma_tally2NonUnit, nloc[igpu], jb, c_one, dB_r(igpu, j%2, j), lddbr,
                                        dA(igpu, 0, j), ldda );

                            if ( j < nbl-1 ) {
                                magma_tally2_cgemm(Magma_tally2NoTrans, Magma_tally2NoTrans, nloc[igpu], n-(j+1)*nb, nb, c_neg_one, dA(igpu, 0, j), ldda,
                                            dB_r(igpu, j%2, j+1), lddbr, c_one, dA(igpu, 0, j+1), ldda );
                            }
                            magma_tally2_event_record( event[igpu][0], stream[igpu][j%2]);
                        }
                    }

                    for (magma_tally2_int_t k = 0; k < j; ++k) {
                        magma_tally2_int_t igpu = k%ngpu;
                        magma_tally2_setdevice(igpu);
                        magma_tally2_int_t kb = min(nb, n-k*nb);
                        magma_tally2_cgetmatrix_async(kb, jb,
                                               dA(igpu, k/ngpu, j), ldda,
                                               A(k, j),              lda, stream[igpu][j%2] );
                    }
                }
            }
        } else {
            /* Compute inv(L)*A*inv(L') */

            // Copy A to mgpu
            for (magma_tally2_int_t k = 0; k < nbl; ++k) {
                magma_tally2_int_t igpu = k%ngpu;
                magma_tally2_setdevice(igpu);
                magma_tally2_int_t kb = min(nb, n-k*nb);
                magma_tally2_csetmatrix_async((n-k*nb), kb,
                                       A(k, k),              lda,
                                       dA(igpu, k, k/ngpu), ldda, stream[igpu][0] );
            }

            for (magma_tally2_int_t k = 0; k < nbl; ++k) {
                magma_tally2_int_t ind_k  =   k   % 2;
                magma_tally2_int_t ind_k1 = (k+1) % 2;
                magma_tally2_int_t kb= min(n-k*nb,nb);

                // Copy B panel
                for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                    magma_tally2_setdevice(igpu);
                    magma_tally2_queue_sync( stream[igpu][0] ); // sync previous B panel copy

                    // sync dwork copy and update (uses B panel of the next copy)
                    magma_tally2_queue_wait_event( stream[igpu][0], event[igpu][1] );
                    magma_tally2_queue_wait_event( stream[igpu][2], event[igpu][1] );
                    magma_tally2_csetmatrix_async((n-k*nb), kb,
                                           B(k, k),          ldb,
                                           dB_c(igpu, k, ind_k), lddbc, stream[igpu][0] );
                }

                magma_tally2_int_t igpu_p = k%ngpu;

                if (k > 0) {
                    // Update the next panel
                    magma_tally2_setdevice(igpu_p);
                    magma_tally2blasSetKernelStream(stream[igpu_p][2]);

                    magma_tally2_int_t nk = n-k*nb;

                    magma_tally2_cher2k(Magma_tally2Lower, Magma_tally2NoTrans, kb, nb,
                                 c_neg_one, dwork(igpu_p, k, 0), n, dB_c(igpu_p, k, ind_k1), lddbc,
                                 d_one, dA(igpu_p, k, k/ngpu), ldda);

                    // copy Akk block on the CPU
                    magma_tally2_cgetmatrix_async(kb, kb,
                                           dA(igpu_p, k, k/ngpu), ldda,
                                           A(k, k),                lda, stream[igpu_p][2] );

                    magma_tally2_event_record( event[igpu_p][0], stream[igpu_p][2]);

                    magma_tally2_cgemm(Magma_tally2NoTrans, Magma_tally2ConjTrans, nk-kb, kb, nb, c_neg_one, dwork(igpu_p, k+1, 0), n,
                                dB_c(igpu_p, k, ind_k1), lddbc, c_one, dA(igpu_p, k+1, k/ngpu), ldda );

                    magma_tally2_cgemm(Magma_tally2NoTrans, Magma_tally2ConjTrans, nk-kb, kb, nb, c_neg_one, dB_c(igpu_p, k+1, ind_k1), lddbc,
                                dwork(igpu_p, k, 0), n, c_one, dA(igpu_p, k+1, k/ngpu), ldda );

                    // Update the panels of the other GPUs
                    for (magma_tally2_int_t j=k+1; j < nbl; ++j) {
                        magma_tally2_int_t igpu = j%ngpu;
                        if (igpu != igpu_p) {
                            magma_tally2_setdevice(igpu);
                            magma_tally2blasSetKernelStream(stream[igpu][1]);

                            magma_tally2_chegst_m_1_L_col_update(n-j*nb, nb, dA(igpu, j, j/ngpu), ldda,
                                                          dwork(igpu, j, 0), n, dB_c(igpu, j, ind_k1), lddbc); // cher2k on j-th column
                        }
                    }
                }
                // compute next panel
                magma_tally2_setdevice(igpu_p);

                if (k+1 < nbl) {
                    magma_tally2_queue_sync( stream[igpu_p][0] ); // sync B panel copy
                    magma_tally2blasSetKernelStream(stream[igpu_p][2]);

                    magma_tally2_ctrsm(Magma_tally2Right, uplo, Magma_tally2ConjTrans, Magma_tally2NonUnit,
                                n-(k+1)*nb, kb,
                                c_one, dB_c(igpu_p, k, ind_k), lddbc,
                                dA(igpu_p, k+1, k/ngpu), ldda);
                }

                magma_tally2_event_sync( event[igpu_p][0] ); // sync Akk copy
                lapackf77_chegst( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);

                if (k+1 < nbl) {
                    magma_tally2_csetmatrix_async(kb, kb,
                                           A(k, k),                lda,
                                           dA(igpu_p, k, k/ngpu), ldda, stream[igpu_p][2] );

                    magma_tally2_chemm(Magma_tally2Right, uplo,
                                n-(k+1)*nb, kb,
                                c_neg_half, dA(igpu_p, k, k/ngpu), ldda,
                                dB_c(igpu_p, k+1, ind_k), lddbc,
                                c_one, dA(igpu_p, k+1, k/ngpu), ldda);

                    magma_tally2_cgetmatrix_async(n-(k+1)*nb, kb,
                                           dA(igpu_p, k+1, k/ngpu), ldda,
                                           A(k+1, k),                lda, stream[igpu_p][2] );
                }

                if (k > 0) {
                    // Update the remaining panels of GPU igpu_p
                    for (magma_tally2_int_t j=k+ngpu; j < nbl; j += ngpu) {
                        magma_tally2_setdevice(igpu_p);
                        magma_tally2blasSetKernelStream(stream[igpu_p][1]);

                        magma_tally2_chegst_m_1_L_col_update(n-j*nb, nb, dA(igpu_p, j, j/ngpu), ldda,
                                                      dwork(igpu_p, j, 0), n, dB_c(igpu_p, j, ind_k1), lddbc); // cher2k on j-th column
                    }
                }

                if (k+1 < nbl) {
                    // send the partially updated panel of dA to each gpu in dwork block
                    magma_tally2_setdevice(igpu_p);
                    magma_tally2_queue_sync( stream[igpu_p][2] );

                    for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally2_setdevice(igpu);
                        magma_tally2_csetmatrix_async((n-(k+1)*nb), kb,
                                               A(k+1, k),          lda,
                                               dwork(igpu, k+1, 0), n, stream[igpu][1] );

                        magma_tally2_event_record( event[igpu][1], stream[igpu][1]);
                    }

                    magma_tally2_setdevice(igpu_p);
                    magma_tally2blasSetKernelStream(stream[igpu_p][1]);

                    magma_tally2_chemm(Magma_tally2Right, uplo,
                                n-(k+1)*nb, kb,
                                c_neg_half, dA(igpu_p, k, k/ngpu), ldda,
                                dB_c(igpu_p, k+1, ind_k), lddbc,
                                c_one, dA(igpu_p, k+1, k/ngpu), ldda);
                }
            }

            for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                magma_tally2_queue_sync( stream[igpu][1] );
            }

            if (n > nb) {
                magma_tally2_int_t nloc[Magma_tally2MaxGPUs] = { 0 };

                for (magma_tally2_int_t j = 1; j < nbl; ++j) {
                    nloc[(j-1)%ngpu] += nb;

                    magma_tally2_int_t jb = min(nb, n-j*nb);

                    for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally2_setdevice(igpu);
                        if (nloc[igpu] > 0) {
                            magma_tally2_csetmatrix_async((n-j*nb), jb,
                                                   B(j, j),            ldb,
                                                   dB_c(igpu, j, j%2), lddbc, stream[igpu][j%2] );

                            magma_tally2_queue_wait_event( stream[igpu][j%2], event[igpu][0] );

                            magma_tally2blasSetKernelStream(stream[igpu][j%2]);
                            magma_tally2_ctrsm(Magma_tally2Left, uplo, Magma_tally2NoTrans, Magma_tally2NonUnit, jb, nloc[igpu], c_one, dB_c(igpu, j, j%2), lddbc,
                                        dA(igpu, j, 0), ldda );

                            if ( j < nbl-1 ) {
                                magma_tally2_cgemm(Magma_tally2NoTrans, Magma_tally2NoTrans, n-(j+1)*nb, nloc[igpu], nb, c_neg_one, dB_c(igpu, j+1, j%2), lddbc,
                                            dA(igpu, j, 0), ldda, c_one, dA(igpu, j+1, 0), ldda );
                            }
                            magma_tally2_event_record( event[igpu][0], stream[igpu][j%2]);
                        }
                    }

                    for (magma_tally2_int_t k = 0; k < j; ++k) {
                        magma_tally2_int_t igpu = k%ngpu;
                        magma_tally2_setdevice(igpu);
                        magma_tally2_int_t kb = min(nb, n-k*nb);
                        magma_tally2_cgetmatrix_async(jb, kb,
                                               dA(igpu, j, k/ngpu), ldda,
                                               A(j, k),              lda, stream[igpu][j%2] );
                    }
                }
            }
        }
    } else {
        if (upper) {
            /* Compute U*A*U' */

            if (n > nb) {
                magma_tally2_int_t nloc[Magma_tally2MaxGPUs] = { 0 };
                magma_tally2_int_t iloc[Magma_tally2MaxGPUs] = { 0 };

                for (magma_tally2_int_t j = 0; j < nbl; ++j) {
                    magma_tally2_int_t jb = min(nb, n-j*nb);
                    nloc[j%ngpu] += jb;
                }

                for (magma_tally2_int_t k = 0; k < nbl; ++k) {
                    magma_tally2_int_t kb = min(nb, n-k*nb);

                    for (magma_tally2_int_t j = k; j < nbl; ++j) {
                        magma_tally2_int_t igpu = j%ngpu;
                        magma_tally2_setdevice(igpu);
                        magma_tally2_int_t jb = min(nb, n-j*nb);
                        magma_tally2_csetmatrix_async(kb, jb,
                                               A(k, j),              lda,
                                               dA(igpu, k, j/ngpu), ldda, stream[igpu][k%2] );
                    }

                    magma_tally2_int_t igpu_p = k % ngpu;

                    ++iloc[igpu_p];
                    nloc[igpu_p] -= kb;

                    for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally2_setdevice(igpu);
                        magma_tally2_csetmatrix_async(k*nb + kb, kb,
                                               B(0, k),            ldb,
                                               dB_c(igpu, 0, k%2), lddbc, stream[igpu][k%2] );

                        magma_tally2_queue_wait_event( stream[igpu][k%2], event[igpu][0] );

                        magma_tally2blasSetKernelStream(stream[igpu][k%2]);

                        if (igpu == igpu_p) {
                            magma_tally2_chemm(Magma_tally2Right, uplo,
                                        k*nb, kb,
                                        c_half, dA(igpu, k, k/ngpu), ldda,
                                        dB_c(igpu, 0, k%2), lddbc,
                                        c_one, dA(igpu, 0, k/ngpu), ldda);
                        }

                        magma_tally2_cgemm(Magma_tally2NoTrans, Magma_tally2NoTrans, k*nb, nloc[igpu], kb, c_one, dB_c(igpu, 0, k%2), lddbc,
                                    dA(igpu, k, iloc[igpu]), ldda, c_one, dA(igpu, 0, iloc[igpu]), ldda );

                        magma_tally2_ctrmm(Magma_tally2Left, uplo, Magma_tally2NoTrans, Magma_tally2NonUnit, kb, nloc[igpu], c_one, dB_c(igpu, k, k%2), lddbc,
                                    dA(igpu, k, iloc[igpu]), ldda );

                        magma_tally2_event_record( event[igpu][0], stream[igpu][k%2]);

                        if (igpu == igpu_p) {
                            magma_tally2_cgetmatrix_async(k*nb, kb,
                                                   dA(igpu, 0, k/ngpu), ldda,
                                                   A(0, k),              lda, stream[igpu][k%2] );
                        }
                    }
                }
            }

            for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                magma_tally2_queue_sync( stream[igpu][0] );
                magma_tally2_queue_sync( stream[igpu][1] );
            }

            for (magma_tally2_int_t k = 0; k < nbl; ++k) {
                magma_tally2_int_t ind_k = k % 2;
                magma_tally2_int_t ind_k1 = (k+1) % 2;
                magma_tally2_int_t kb= min(n-k*nb,nb);

                magma_tally2_int_t igpu_p = k%ngpu;

                if (k > 0) {
                    for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally2_setdevice(igpu);

                        magma_tally2_queue_wait_event( stream[igpu][0], event[igpu][ind_k] ); // sync computation that use the B panel of next copy

                        magma_tally2_csetmatrix_async(k*nb+kb, kb,
                                               B(0, k),              ldb,
                                               dB_c(igpu, 0, ind_k), lddbc, stream[igpu][0] );

                        magma_tally2_event_record(event[igpu][ind_k], stream[igpu][0]);

                        magma_tally2_csetmatrix_async(k*nb, kb,
                                               A(0, k),           lda,
                                               dwork(igpu, 0, 0), n, stream[igpu][1] );

                        magma_tally2_queue_wait_event( stream[igpu][1], event[igpu][ind_k] ); // sync B copy
                    }

                    magma_tally2_setdevice(igpu_p);
                    magma_tally2blasSetKernelStream(stream[igpu_p][2]);

                    magma_tally2_queue_wait_event( stream[igpu_p][2], event[igpu_p][ind_k1] ); // sync update of previous step
                    magma_tally2_queue_wait_event( stream[igpu_p][2], event[igpu_p][ind_k] ); // sync B copy

                    magma_tally2_chemm(Magma_tally2Right, uplo,
                                k*nb, kb,
                                c_half, dA(igpu_p, k, k/ngpu), ldda,
                                dB_c(igpu_p, 0, ind_k), lddbc,
                                c_one, dA(igpu_p, 0, k/ngpu), ldda);

                    magma_tally2_ctrmm(Magma_tally2Right, uplo, Magma_tally2ConjTrans, Magma_tally2NonUnit,
                                k*nb, kb,
                                c_one, dB_c(igpu_p, k, ind_k), lddbc,
                                dA(igpu_p, 0, k/ngpu), ldda);

                    magma_tally2_event_record(event[igpu_p][ind_k], stream[igpu_p][2]);
                    magma_tally2_queue_wait_event(stream[igpu_p][1], event[igpu_p][ind_k]);

                    for (magma_tally2_int_t j = 0; j < k; ++j) {
                        magma_tally2_int_t igpu = j%ngpu;
                        magma_tally2_setdevice(igpu);
                        magma_tally2blasSetKernelStream(stream[igpu][1]);

                        magma_tally2_cher2k(uplo, Magma_tally2NoTrans,
                                     nb, kb,
                                     c_one, dwork(igpu, j, 0), n,
                                     dB_c(igpu, j, ind_k), lddbc,
                                     d_one, dA(igpu, j, j/ngpu), ldda);

                        magma_tally2_cgemm(Magma_tally2NoTrans, Magma_tally2ConjTrans, j*nb, nb, kb, c_one, dB_c(igpu, 0, ind_k), lddbc,
                                    dwork(igpu, j, 0), n, c_one, dA(igpu, 0, j/ngpu), ldda );

                        magma_tally2_cgemm(Magma_tally2NoTrans, Magma_tally2ConjTrans, j*nb, nb, kb, c_one, dwork(igpu, 0, 0), n,
                                    dB_c(igpu, j, ind_k), lddbc, c_one, dA(igpu, 0, j/ngpu), ldda );
                    }
                }

                for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                    magma_tally2_setdevice(igpu);
                    magma_tally2_event_record(event[igpu][ind_k], stream[igpu][1]);
                    magma_tally2_queue_sync(stream[igpu][0]); // sync B copy (conflicts with chegst)
                }

                lapackf77_chegst( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);

                magma_tally2_setdevice(igpu_p);

                magma_tally2_csetmatrix_async(kb, kb,
                                       A(k, k),                lda,
                                       dA(igpu_p, k, k/ngpu), ldda, stream[igpu_p][1] );
            }

            for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                magma_tally2_queue_sync( stream[igpu][1] );
            }

            //copy A from mgpus
            for (magma_tally2_int_t j = 0; j < nbl; ++j) {
                magma_tally2_int_t igpu = j%ngpu;
                magma_tally2_setdevice(igpu);
                magma_tally2_int_t jb = min(nb, n-j*nb);
                magma_tally2_cgetmatrix_async(j*nb+jb, jb,
                                       dA(igpu, 0, j/ngpu), ldda,
                                       A(0, j),              lda, stream[igpu][0] );
            }
        } else {
            /* Compute L'*A*L */

            if (n > nb) {
                magma_tally2_int_t nloc[Magma_tally2MaxGPUs] = { 0 };
                magma_tally2_int_t iloc[Magma_tally2MaxGPUs] = { 0 };

                for (magma_tally2_int_t j = 0; j < nbl; ++j) {
                    magma_tally2_int_t jb = min(nb, n-j*nb);
                    nloc[j%ngpu] += jb;
                }

                for (magma_tally2_int_t k = 0; k < nbl; ++k) {
                    magma_tally2_int_t kb = min(nb, n-k*nb);

                    for (magma_tally2_int_t j = k; j < nbl; ++j) {
                        magma_tally2_int_t igpu = j%ngpu;
                        magma_tally2_setdevice(igpu);
                        magma_tally2_int_t jb = min(nb, n-j*nb);

                        magma_tally2_csetmatrix_async(jb, kb,
                                               A(j, k),              lda,
                                               dA(igpu, j/ngpu, k), ldda, stream[igpu][k%2] );
                    }

                    magma_tally2_int_t igpu_p = k % ngpu;

                    ++iloc[igpu_p];
                    nloc[igpu_p] -= kb;

                    for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally2_setdevice(igpu);
                        magma_tally2_csetmatrix_async(kb, k*nb +kb,
                                               B(k, 0),            ldb,
                                               dB_r(igpu, k%2, 0), lddbr, stream[igpu][k%2] );

                        magma_tally2_queue_wait_event( stream[igpu][k%2], event[igpu][0] );

                        magma_tally2blasSetKernelStream(stream[igpu][k%2]);

                        if (igpu == igpu_p) {
                            magma_tally2_chemm(Magma_tally2Left, uplo,
                                        kb, k*nb,
                                        c_half, dA(igpu, k/ngpu, k), ldda,
                                        dB_r(igpu, k%2, 0), lddbr,
                                        c_one, dA(igpu, k/ngpu, 0), ldda);
                        }

                        magma_tally2_cgemm(Magma_tally2NoTrans, Magma_tally2NoTrans, nloc[igpu], k*nb, kb, c_one, dA(igpu, iloc[igpu], k), ldda,
                                    dB_r(igpu, k%2, 0), lddbr, c_one, dA(igpu, iloc[igpu], 0), ldda );

                        magma_tally2_ctrmm(Magma_tally2Right, uplo, Magma_tally2NoTrans, Magma_tally2NonUnit, nloc[igpu], kb, c_one, dB_r(igpu, k%2, k), lddbr,
                                    dA(igpu, iloc[igpu], k), ldda );

                        magma_tally2_event_record( event[igpu][0], stream[igpu][k%2]);

                        if (igpu == igpu_p) {
                            magma_tally2_cgetmatrix_async(kb, k*nb,
                                                   dA(igpu, k/ngpu, 0), ldda,
                                                   A(k, 0),              lda, stream[igpu][k%2] );
                        }
                    }
                }
            }

            for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                magma_tally2_queue_sync( stream[igpu][0] );
                magma_tally2_queue_sync( stream[igpu][1] );
            }

            for (magma_tally2_int_t k = 0; k < nbl; ++k) {
                magma_tally2_int_t ind_k = k % 2;
                magma_tally2_int_t ind_k1 = (k+1) % 2;
                magma_tally2_int_t kb= min(n-k*nb,nb);

                magma_tally2_int_t igpu_p = k%ngpu;

                if (k > 0) {
                    for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                        magma_tally2_setdevice(igpu);

                        magma_tally2_queue_wait_event( stream[igpu][0], event[igpu][ind_k] ); // sync computation that use the B panel of next copy

                        magma_tally2_csetmatrix_async(kb, k*nb+kb,
                                               B(k, 0),              ldb,
                                               dB_r(igpu, ind_k, 0), lddbr, stream[igpu][0] );

                        magma_tally2_event_record(event[igpu][ind_k], stream[igpu][0]);

                        magma_tally2_csetmatrix_async(kb, k*nb,
                                               A(k, 0),           lda,
                                               dwork(igpu, 0, 0), nb, stream[igpu][1] );

                        magma_tally2_queue_wait_event( stream[igpu][1], event[igpu][ind_k] ); // sync B copy
                    }

                    magma_tally2_setdevice(igpu_p);
                    magma_tally2blasSetKernelStream(stream[igpu_p][2]);

                    magma_tally2_queue_wait_event( stream[igpu_p][2], event[igpu_p][ind_k1] ); // sync update of previous step
                    magma_tally2_queue_wait_event( stream[igpu_p][2], event[igpu_p][ind_k] ); // sync B copy

                    magma_tally2_chemm(Magma_tally2Left, uplo,
                                kb, k*nb,
                                c_half, dA(igpu_p, k/ngpu, k), ldda,
                                dB_r(igpu_p, ind_k, 0), lddbr,
                                c_one, dA(igpu_p, k/ngpu, 0), ldda);

                    magma_tally2_ctrmm(Magma_tally2Left, uplo, Magma_tally2ConjTrans, Magma_tally2NonUnit,
                                kb, k*nb,
                                c_one, dB_r(igpu_p, ind_k, k), lddbr,
                                dA(igpu_p, k/ngpu, 0), ldda);

                    magma_tally2_event_record(event[igpu_p][ind_k], stream[igpu_p][2]);
                    magma_tally2_queue_wait_event(stream[igpu_p][1], event[igpu_p][ind_k]);

                    for (magma_tally2_int_t j = 0; j < k; ++j) {
                        magma_tally2_int_t igpu = j%ngpu;
                        magma_tally2_setdevice(igpu);
                        magma_tally2blasSetKernelStream(stream[igpu][1]);

                        magma_tally2_cher2k(uplo, Magma_tally2ConjTrans,
                                     nb, kb,
                                     c_one, dwork(igpu, 0, j), nb,
                                     dB_r(igpu, ind_k, j), lddbr,
                                     d_one, dA(igpu, j/ngpu, j), ldda);

                        magma_tally2_cgemm(Magma_tally2ConjTrans, Magma_tally2NoTrans, nb, j*nb, kb, c_one, dwork(igpu, 0, j), nb,
                                    dB_r(igpu, ind_k, 0), lddbr, c_one, dA(igpu, j/ngpu, 0), ldda );

                        magma_tally2_cgemm(Magma_tally2ConjTrans, Magma_tally2NoTrans, nb, j*nb, kb, c_one, dB_r(igpu, ind_k, j), lddbr,
                                    dwork(igpu, 0, 0), nb, c_one, dA(igpu, j/ngpu, 0), ldda );
                    }
                }

                for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                    magma_tally2_setdevice(igpu);
                    magma_tally2_event_record(event[igpu][ind_k], stream[igpu][1]);
                    magma_tally2_queue_sync(stream[igpu][0]); // sync B copy (conflicts with chegst)
                }

                lapackf77_chegst( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);

                magma_tally2_setdevice(igpu_p);

                magma_tally2_csetmatrix_async(kb, kb,
                                       A(k, k),                lda,
                                       dA(igpu_p, k/ngpu, k), ldda, stream[igpu_p][1] );
            }

            for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
                magma_tally2_queue_sync( stream[igpu][1] );
            }

            //copy A from mgpus
            for (magma_tally2_int_t j = 0; j < nbl; ++j) {
                magma_tally2_int_t igpu = j%ngpu;
                magma_tally2_setdevice(igpu);
                magma_tally2_int_t jb = min(nb, n-j*nb);
                magma_tally2_cgetmatrix_async(jb, j*nb+jb,
                                       dA(igpu, j/ngpu, 0), ldda,
                                       A(j, 0),              lda, stream[igpu][0] );
            }
        }
    }

    for (magma_tally2_int_t igpu = 0; igpu < ngpu; ++igpu) {
        magma_tally2_setdevice(igpu);
        magma_tally2_event_destroy( event[igpu][0] );
        magma_tally2_event_destroy( event[igpu][1] );
        magma_tally2_queue_destroy( stream[igpu][0] );
        magma_tally2_queue_destroy( stream[igpu][1] );
        magma_tally2_queue_destroy( stream[igpu][2] );
        magma_tally2_free( dw[igpu] );
    }
    magma_tally2_setdevice( orig_dev );
    magma_tally2blasSetKernelStream( orig_stream );

    return *info;
} /* magma_tally2_chegst_gpu */


inline static void magma_tally2_chegst_m_1_U_row_update(
    magma_tally2_int_t nk, magma_tally2_int_t nb,
    magma_tally2FloatComplex* dA_row, magma_tally2_int_t ldda,
    magma_tally2FloatComplex* dC1, magma_tally2_int_t lddc1,
    magma_tally2FloatComplex* dC2, magma_tally2_int_t lddc2)
{
    // update 1 rowblock (rowwise cher2k) for itype=1 Upper case
    float             d_one      = 1.0;
    magma_tally2FloatComplex c_one      = MAGMA_tally2_C_ONE;
    magma_tally2FloatComplex c_neg_one  = MAGMA_tally2_C_NEG_ONE;

    magma_tally2_int_t kb = min(nk, nb);

    magma_tally2_cher2k(Magma_tally2Upper, Magma_tally2ConjTrans, kb, nb,
                 c_neg_one, dC1, lddc1, dC2, lddc2,
                 d_one, dA_row, ldda);

    magma_tally2_cgemm(Magma_tally2ConjTrans, Magma_tally2NoTrans, kb, nk-kb, nb, c_neg_one, dC1, lddc1,
                dC2+kb*lddc2, lddc2, c_one, dA_row+kb*ldda, ldda );

    magma_tally2_cgemm(Magma_tally2ConjTrans, Magma_tally2NoTrans, kb, nk-kb, nb, c_neg_one, dC2, lddc2,
                dC1+kb*lddc1, lddc1, c_one, dA_row+kb*ldda, ldda );
}


inline static void magma_tally2_chegst_m_1_L_col_update(
    magma_tally2_int_t nk, magma_tally2_int_t nb,
    magma_tally2FloatComplex* dA_col, magma_tally2_int_t ldda,
    magma_tally2FloatComplex* dC1, magma_tally2_int_t lddc1,
    magma_tally2FloatComplex* dC2, magma_tally2_int_t lddc2)
{
    // update 1 columnblock (columnwise cher2k) for itype=1 Lower case
    float             d_one      = 1.0;
    magma_tally2FloatComplex c_one      = MAGMA_tally2_C_ONE;
    magma_tally2FloatComplex c_neg_one  = MAGMA_tally2_C_NEG_ONE;

    magma_tally2_int_t kb = min(nk, nb);

    magma_tally2_cher2k(Magma_tally2Lower, Magma_tally2NoTrans, kb, nb,
                 c_neg_one, dC1, lddc1, dC2, lddc2,
                 d_one, dA_col, ldda);

    magma_tally2_cgemm(Magma_tally2NoTrans, Magma_tally2ConjTrans, nk-kb, kb, nb, c_neg_one, dC1+kb, lddc1,
                dC2, lddc2, c_one, dA_col+kb, ldda );

    magma_tally2_cgemm(Magma_tally2NoTrans, Magma_tally2ConjTrans, nk-kb, kb, nb, c_neg_one, dC2+kb, lddc2,
                dC1, lddc1, c_one, dA_col+kb, ldda );
}