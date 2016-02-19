/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
  
       @generated from zgeqp3_gpu.cpp normal z -> d, Fri Jan 30 19:00:15 2015

*/
#include "common_magma_tally3.h"

#define PRECISION_d
#define REAL

/**
    Purpose
    -------
    DGEQP3 computes a QR factorization with column pivoting of a
    matrix A:  A*P = Q*R  using Level 3 BLAS.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A. M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in,out]
    dA      DOUBLE_PRECISION array on the GPU, dimension (LDDA,N)
            On entry, the M-by-N matrix A.
            On exit, the upper triangle of the array contains the
            min(M,N)-by-N upper trapezoidal matrix R; the elements below
            the diagonal, together with the array TAU, represent the
            unitary matrix Q as a product of min(M,N) elementary
            reflectors.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A. LDDA >= max(1,M).

    @param[in,out]
    jpvt    INTEGER array, dimension (N)
            On entry, if JPVT(J).ne.0, the J-th column of A is permuted
            to the front of A*P (a leading column); if JPVT(J)=0,
            the J-th column of A is a free column.
            On exit, if JPVT(J)=K, then the J-th column of A*P was the
            the K-th column of A.

    @param[out]
    tau     DOUBLE_PRECISION array, dimension (min(M,N))
            The scalar factors of the elementary reflectors.

    @param[out]
    dwork   (workspace) DOUBLE_PRECISION array on the GPU, dimension (MAX(1,LWORK))
            On exit, if INFO=0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK.
            For [sd]geqp3, LWORK >= (N+1)*NB + 2*N;
            for [cz]geqp3, LWORK >= (N+1)*NB,
            where NB is the optimal blocksize.
    \n
            Note: unlike the CPU interface of this routine, the GPU interface
            does not support a workspace query.

    @param
    rwork   (workspace, for [cz]geqp3 only) DOUBLE PRECISION array, dimension (2*N)

    @param[out]
    info    INTEGER
      -     = 0: successful exit.
      -     < 0: if INFO = -i, the i-th argument had an illegal value.

    Further Details
    ---------------
    The matrix Q is represented as a product of elementary reflectors

        Q = H(1) H(2) . . . H(k), where k = min(m,n).

    Each H(i) has the form

        H(i) = I - tau * v * v'

    where tau is a real scalar, and v is a real vector
    with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in
    A(i+1:m,i), and tau in TAU(i).

    @ingroup magma_tally3_dgeqp3_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_dgeqp3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt, double *tau,
    magma_tally3Double_ptr dwork, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *info )
{
    #define dA(i_, j_) (dA + (i_) + (j_)*ldda)

    magma_tally3_int_t ione = 1;

    //magma_tally3_int_t na;
    magma_tally3_int_t n_j;
    magma_tally3_int_t j, jb, nb, sm, sn, fjb, nfxd, minmn;
    magma_tally3_int_t topbmn, lwkopt;
    
    *info = 0;
    if (m < 0) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (ldda < max(1,m)) {
        *info = -4;
    }
    
    nb = magma_tally3_get_dgeqp3_nb(min(m, n));
    minmn = min(m,n);
    if (*info == 0) {
        if (minmn == 0) {
            lwkopt = 1;
        } else {
            lwkopt = (n + 1)*nb;
            #ifdef REAL
            lwkopt += 2*n;
            #endif
        }
        if (lwork < lwkopt) {
            *info = -8;
        }
    }

    if (*info != 0) {
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }

    if (minmn == 0)
        return *info;

    #ifdef REAL
    double *rwork = dwork + (n + 1)*nb;
    #endif
    magma_tally3Double_ptr df;
    if (MAGMA_tally3_SUCCESS != magma_tally3_dmalloc( &df, (n+1)*nb )) {
        *info = MAGMA_tally3_ERR_DEVICE_ALLOC;
        return *info;
    }
    cudaMemset( df, 0, (n+1)*nb*sizeof(double) );

    nfxd = 0;
    /* Move initial columns up front.
     * Note jpvt uses 1-based indices for historical compatibility. */
    for (j = 0; j < n; ++j) {
        if (jpvt[j] != 0) {
            if (j != nfxd) {
                blasf77_dswap(&m, dA(0, j), &ione, dA(0, nfxd), &ione);  // TODO: ERROR, matrix not on CPU!
                jpvt[j]    = jpvt[nfxd];
                jpvt[nfxd] = j + 1;
            }
            else {
                jpvt[j] = j + 1;
            }
            ++nfxd;
        }
        else {
            jpvt[j] = j + 1;
        }
    }

    /*
        // TODO:
           Factorize fixed columns
           =======================
           Compute the QR factorization of fixed columns and update
           remaining columns.
    if (nfxd > 0) {
        na = min(m,nfxd);
        lapackf77_dgeqrf(&m, &na, dA, &ldda, tau, dwork, &lwork, info);
        if (na < n) {
            n_j = n - na;
            lapackf77_dormqr( Magma_tally3LeftStr, Magma_tally3ConjTransStr, &m, &n_j, &na,
                              dA, &ldda, tau, dA(0, na), &ldda,
                              dwork, &lwork, info );
        }
    }*/
    
    /*  Factorize free columns */
    if (nfxd < minmn) {
        sm = m - nfxd;
        sn = n - nfxd;
        //sminmn = minmn - nfxd;
        
        /* Initialize partial column norms. */
        magma_tally3blas_dnrm2_cols(sm, sn, dA(nfxd,nfxd), ldda, &rwork[nfxd]);
        magma_tally3_dcopymatrix( sn, 1, &rwork[nfxd], sn, &rwork[n+nfxd], sn);
        
        j = nfxd;
        //if (nb < sminmn)
        {
            /* Use blocked code initially. */
            
            /* Compute factorization: while loop. */
            topbmn = minmn; // - nb;
            while(j < topbmn) {
                jb = min(nb, topbmn - j);
                
                /* Factorize JB columns among columns J:N. */
                n_j = n - j;
                
                //magma_tally3_dlaqps_gpu    // this is a cpp-file
                magma_tally3_dlaqps2_gpu   // this is a cuda-file
                    ( m, n_j, j, jb, &fjb,
                      dA(0, j), ldda,
                      &jpvt[j], &tau[j], &rwork[j], &rwork[n + j],
                      dwork,
                      &df[jb],   n_j );
                
                j += fjb;  /* fjb is actual number of columns factored */
            }
        }
        
        /*
        // Use unblocked code to factor the last or only block.
        if (j < minmn) {
            n_j = n - j;
            if (j > nfxd) {
                magma_tally3_dgetmatrix( m-j, n_j,
                                  dA(j,j), ldda,
                                   A(j,j), lda );
            }
            lapackf77_dlaqp2(&m, &n_j, &j, dA(0, j), &ldda, &jpvt[j],
                             &tau[j], &rwork[j], &rwork[n+j], dwork );
        }*/
    }

    magma_tally3_free(df);

    return *info;
} /* magma_tally3_dgeqp3_gpu */
