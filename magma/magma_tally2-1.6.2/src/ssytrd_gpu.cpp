/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Stan Tomov
       @author Mark Gates

       @generated from zhetrd_gpu.cpp normal z -> s, Fri Jan 30 19:00:17 2015

*/
#include "common_magma_tally2.h"

/**
    Purpose
    -------
    SSYTRD_GPU reduces a real symmetric matrix A to real symmetric
    tridiagonal form T by an orthogonal similarity transformation:
    Q**H * A * Q = T.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally2_uplo_t
      -     = Magma_tally2Upper:  Upper triangle of A is stored;
      -     = Magma_tally2Lower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    dA      REAL array on the GPU, dimension (LDDA,N)
            On entry, the symmetric matrix A.  If UPLO = Magma_tally2Upper, the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_tally2Lower, the
            leading N-by-N lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
            On exit, if UPLO = Magma_tally2Upper, the diagonal and first superdiagonal
            of A are overwritten by the corresponding elements of the
            tridiagonal matrix T, and the elements above the first
            superdiagonal, with the array TAU, represent the orthogonal
            matrix Q as a product of elementary reflectors; if UPLO
            = Magma_tally2Lower, the diagonal and first subdiagonal of A are over-
            written by the corresponding elements of the tridiagonal
            matrix T, and the elements below the first subdiagonal, with
            the array TAU, represent the orthogonal matrix Q as a product
            of elementary reflectors. See Further Details.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[out]
    d       REAL array, dimension (N)
            The diagonal elements of the tridiagonal matrix T:
            D(i) = A(i,i).

    @param[out]
    e       REAL array, dimension (N-1)
            The off-diagonal elements of the tridiagonal matrix T:
            E(i) = A(i,i+1) if UPLO = Magma_tally2Upper, E(i) = A(i+1,i) if UPLO = Magma_tally2Lower.

    @param[out]
    tau     REAL array, dimension (N-1)
            The scalar factors of the elementary reflectors (see Further
            Details).

    @param[out]
    A       (workspace) REAL array, dimension (LDA,N)
            On exit the diagonal, the  upper part (UPLO=Magma_tally2Upper)
            or the lower part (UPLO=Magma_tally2Lower) are copies of DA

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[out]
    work    (workspace) REAL array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK.  LWORK >= N*NB, where NB is the
            optimal blocksize given by magma_tally2_get_ssytrd_nb().
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued by XERBLA.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    Further Details
    ---------------
    If UPLO = Magma_tally2Upper, the matrix Q is represented as a product of elementary
    reflectors

       Q = H(n-1) . . . H(2) H(1).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a real scalar, and v is a real vector with
    v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
    A(1:i-1,i+1), and tau in TAU(i).

    If UPLO = Magma_tally2Lower, the matrix Q is represented as a product of elementary
    reflectors

       Q = H(1) H(2) . . . H(n-1).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a real scalar, and v is a real vector with
    v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
    and tau in TAU(i).

    The contents of A on exit are illustrated by the following examples
    with n = 5:

    if UPLO = Magma_tally2Upper:                if UPLO = Magma_tally2Lower:

      (  d   e   v2  v3  v4 )              (  d                  )
      (      d   e   v3  v4 )              (  e   d              )
      (          d   e   v4 )              (  v1  e   d          )
      (              d   e  )              (  v1  v2  e   d      )
      (                  d  )              (  v1  v2  v3  e   d  )

    where d and e denote diagonal and off-diagonal elements of T, and vi
    denotes an element of the vector defining H(i).

    @ingroup magma_tally2_ssyev_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_ssytrd_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *d, float *e, float *tau,
    float *A,  magma_tally2_int_t lda,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info)
{
    #define  A(i_, j_) ( A + (i_) + (j_)*lda )
    #define dA(i_, j_) (dA + (i_) + (j_)*ldda)

    const char* uplo_ = lapack_uplo_const_tally2( uplo );

    magma_tally2_int_t nb = magma_tally2_get_ssytrd_nb( n );

    const float c_zero    = MAGMA_tally2_S_ZERO;
    const float c_neg_one = MAGMA_tally2_S_NEG_ONE;
    const float c_one     = MAGMA_tally2_S_ONE;
    const float             d_one     = MAGMA_tally2_D_ONE;
    
    magma_tally2_int_t kk, nx;
    magma_tally2_int_t i, j, i_n;
    magma_tally2_int_t iinfo;
    magma_tally2_int_t ldw, lddw, lwkopt;
    magma_tally2_int_t lquery;

    *info = 0;
    int upper = (uplo == Magma_tally2Upper);
    lquery = (lwork == -1);
    if (! upper && uplo != Magma_tally2Lower) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (ldda < max(1,n)) {
        *info = -4;
    } else if (lda < max(1,n)) {
        *info = -9;
    } else if (lwork < nb*n && ! lquery) {
        *info = -11;
    }

    /* Determine the block size. */
    ldw = n;
    lddw = roundup( n, 32 );
    lwkopt = n * nb;
    if (*info == 0) {
        work[0] = MAGMA_tally2_S_MAKE( lwkopt, 0 );
    }

    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery)
        return *info;

    /* Quick return if possible */
    if (n == 0) {
        work[0] = c_one;
        return *info;
    }

    //if (n < 2048)
    //    nx = n;
    //else
    //    nx = 512;
    nx = min( 128, n );  // nx <= n is required
    
    magma_tally2Float_ptr dwork;
    if (MAGMA_tally2_SUCCESS != magma_tally2_smalloc( &dwork, lddw*nb )) {
        *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
        return *info;
    }

    // clear out dwork in case it has NANs (used as y in ssymv)
    // rest of dwork (used as work in magma_tally2blas_ssymv) doesn't need to be cleared
    magma_tally2blas_slaset( Magma_tally2Full, n, nb, c_zero, c_zero, dwork, lddw );
    
    if (upper) {
        /* Reduce the upper triangle of A.
           Columns 1:kk are handled by the unblocked method. */
        kk = n - (n - nx + nb - 1) / nb * nb;
        
        for (i = n - nb; i >= kk; i -= nb) {
            /* Reduce columns i:i+nb-1 to tridiagonal form and form the
               matrix W which is needed to update the unreduced part of
               the matrix */
            
            /* Get the current panel */
            magma_tally2_sgetmatrix( i+nb, nb, dA(0, i), ldda, A(0, i), lda );
            
            magma_tally2_slatrd( uplo, i+nb, nb, A(0, 0), lda, e, tau,
                          work, ldw, dA(0, 0), ldda, dwork, lddw );
            
            /* Update the unreduced submatrix A(0:i-2,0:i-2), using an
               update of the form:  A := A - V*W' - W*V' */
            magma_tally2_ssetmatrix( i + nb, nb, work, ldw, dwork, lddw );
            
            magma_tally2_ssyr2k( uplo, Magma_tally2NoTrans, i, nb, c_neg_one,
                          dA(0, i), ldda, dwork, lddw,
                          d_one, dA(0, 0), ldda );
            
            /* Copy superdiagonal elements back into A, and diagonal
               elements into D */
            for (j = i; j < i+nb; ++j) {
                *A(j-1,j) = MAGMA_tally2_S_MAKE( e[j - 1], 0 );
                d[j] = MAGMA_tally2_S_REAL( *A(j, j) );
            }
        }
        
        magma_tally2_sgetmatrix( kk, kk, dA(0, 0), ldda, A(0, 0), lda );
        
        /* Use CPU code to reduce the last or only block */
        lapackf77_ssytrd( uplo_, &kk, A(0, 0), &lda, d, e, tau, work, &lwork, &iinfo );
        
        magma_tally2_ssetmatrix( kk, kk, A(0, 0), lda, dA(0, 0), ldda );
    }
    else {
        /* Reduce the lower triangle of A */
        for (i = 0; i < n-nx; i += nb) {
            /* Reduce columns i:i+nb-1 to tridiagonal form and form the
               matrix W which is needed to update the unreduced part of
               the matrix */
            
            /* Get the current panel */
            magma_tally2_sgetmatrix( n-i, nb, dA(i, i), ldda, A(i, i), lda );
            
            magma_tally2_slatrd( uplo, n-i, nb, A(i, i), lda, &e[i], &tau[i],
                          work, ldw, dA(i, i), ldda, dwork, lddw );
            
            /* Update the unreduced submatrix A(i+ib:n,i+ib:n), using
               an update of the form:  A := A - V*W' - W*V' */
            magma_tally2_ssetmatrix( n-i, nb, work, ldw, dwork, lddw );
            
            // cublas 6.5 crashes here if lddw % 32 != 0, e.g., N=250.
            magma_tally2_ssyr2k( Magma_tally2Lower, Magma_tally2NoTrans, n-i-nb, nb, c_neg_one,
                          dA(i+nb, i), ldda, &dwork[nb], lddw,
                          d_one, dA(i+nb, i+nb), ldda );
            
            /* Copy subdiagonal elements back into A, and diagonal
               elements into D */
            for (j = i; j < i+nb; ++j) {
                *A(j+1,j) = MAGMA_tally2_S_MAKE( e[j], 0 );
                d[j] = MAGMA_tally2_S_REAL( *A(j, j) );
            }
        }
        
        /* Use CPU code to reduce the last or only block */
        magma_tally2_sgetmatrix( n-i, n-i, dA(i, i), ldda, A(i, i), lda );
        
        i_n = n-i;
        lapackf77_ssytrd( uplo_, &i_n, A(i, i), &lda, &d[i], &e[i],
                          &tau[i], work, &lwork, &iinfo );
        
        magma_tally2_ssetmatrix( n-i, n-i, A(i, i), lda, dA(i, i), ldda );
    }
    
    magma_tally2_free( dwork );
    
    work[0] = MAGMA_tally2_S_MAKE( lwkopt, 0 );
    
    return *info;
} /* magma_tally2_ssytrd_gpu */
