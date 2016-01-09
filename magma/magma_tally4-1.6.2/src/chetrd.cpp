/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Stan Tomov
       @author Mark Gates

       @generated from zhetrd.cpp normal z -> c, Fri Jan 30 19:00:17 2015

*/
#include "common_magma_tally4.h"

#define FAST_HEMV

/**
    Purpose
    -------
    CHETRD reduces a complex Hermitian matrix A to real symmetric
    tridiagonal form T by an orthogonal similarity transformation:
    Q**H * A * Q = T.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally4_uplo_t
      -     = Magma_tally4Upper:  Upper triangle of A is stored;
      -     = Magma_tally4Lower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    A       COMPLEX array, dimension (LDA,N)
            On entry, the Hermitian matrix A.  If UPLO = Magma_tally4Upper, the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_tally4Lower, the
            leading N-by-N lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
            On exit, if UPLO = Magma_tally4Upper, the diagonal and first superdiagonal
            of A are overwritten by the corresponding elements of the
            tridiagonal matrix T, and the elements above the first
            superdiagonal, with the array TAU, represent the orthogonal
            matrix Q as a product of elementary reflectors; if UPLO
            = Magma_tally4Lower, the diagonal and first subdiagonal of A are over-
            written by the corresponding elements of the tridiagonal
            matrix T, and the elements below the first subdiagonal, with
            the array TAU, represent the orthogonal matrix Q as a product
            of elementary reflectors. See Further Details.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[out]
    d       COMPLEX array, dimension (N)
            The diagonal elements of the tridiagonal matrix T:
            D(i) = A(i,i).

    @param[out]
    e       COMPLEX array, dimension (N-1)
            The off-diagonal elements of the tridiagonal matrix T:
            E(i) = A(i,i+1) if UPLO = Magma_tally4Upper, E(i) = A(i+1,i) if UPLO = Magma_tally4Lower.

    @param[out]
    tau     COMPLEX array, dimension (N-1)
            The scalar factors of the elementary reflectors (see Further
            Details).

    @param[out]
    work    (workspace) COMPLEX array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK.  LWORK >= N*NB, where NB is the
            optimal blocksize given by magma_tally4_get_chetrd_nb().
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
    If UPLO = Magma_tally4Upper, the matrix Q is represented as a product of elementary
    reflectors

       Q = H(n-1) . . . H(2) H(1).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a complex scalar, and v is a complex vector with
    v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
    A(1:i-1,i+1), and tau in TAU(i).

    If UPLO = Magma_tally4Lower, the matrix Q is represented as a product of elementary
    reflectors

       Q = H(1) H(2) . . . H(n-1).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a complex scalar, and v is a complex vector with
    v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
    and tau in TAU(i).

    The contents of A on exit are illustrated by the following examples
    with n = 5:

    if UPLO = Magma_tally4Upper:                if UPLO = Magma_tally4Lower:

      (  d   e   v2  v3  v4 )              (  d                  )
      (      d   e   v3  v4 )              (  e   d              )
      (          d   e   v4 )              (  v1  e   d          )
      (              d   e  )              (  v1  v2  e   d      )
      (                  d  )              (  v1  v2  v3  e   d  )

    where d and e denote diagonal and off-diagonal elements of T, and vi
    denotes an element of the vector defining H(i).

    @ingroup magma_tally4_cheev_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_chetrd(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    float *d, float *e, magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info)
{
    #define  A(i_, j_) ( A + (i_) + (j_)*lda )
    #define dA(i_, j_) (dA + (i_) + (j_)*ldda)

    const char* uplo_ = lapack_uplo_const_tally4( uplo );

    magma_tally4_int_t ldda = roundup( n, 32 );
    magma_tally4_int_t nb = magma_tally4_get_chetrd_nb( n );

    const magma_tally4FloatComplex c_zero    = MAGMA_tally4_C_ZERO;
    const magma_tally4FloatComplex c_neg_one = MAGMA_tally4_C_NEG_ONE;
    const magma_tally4FloatComplex c_one     = MAGMA_tally4_C_ONE;
    const float             d_one     = MAGMA_tally4_D_ONE;
    
    magma_tally4_int_t kk, nx;
    magma_tally4_int_t i, j, i_n;
    magma_tally4_int_t iinfo;
    magma_tally4_int_t ldw, lddw, lwkopt;
    magma_tally4_int_t lquery;

    *info = 0;
    int upper = (uplo == Magma_tally4Upper);
    lquery = (lwork == -1);
    if (! upper && uplo != Magma_tally4Lower) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (lda < max(1,n)) {
        *info = -4;
    } else if (lwork < nb*n && ! lquery) {
        *info = -9;
    }

    /* Determine the block size. */
    ldw = n;
    lddw = ldda;
    lwkopt = n * nb;
    if (*info == 0) {
        work[0] = MAGMA_tally4_C_MAKE( lwkopt, 0 );
    }

    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery)
        return *info;

    /* Quick return if possible */
    if (n == 0) {
        work[0] = c_one;
        return *info;
    }

    magma_tally4FloatComplex *dA;
    #ifdef FAST_HEMV
    magma_tally4_int_t ldwork2 = ldda*ceildiv(n,64);
    #else
    magma_tally4_int_t ldwork2 = 0;
    #endif
    if (MAGMA_tally4_SUCCESS != magma_tally4_cmalloc( &dA, ldda*n + 2*lddw*nb + ldwork2 )) {
        *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
        return *info;
    }

    magma_tally4FloatComplex *dwork  = dA + ldda*n;
    #ifdef FAST_HEMV
    magma_tally4FloatComplex *dwork2 = dwork + 2*lddw*nb;
    #endif

    //if (n < 2048)
    //    nx = n;
    //else
    //    nx = 512;
    nx = min( 128, n );  // nx <= n is required

    // clear out dwork in case it has NANs (used as y in chemv)
    // rest of dwork (used as work in magma_tally4blas_chemv) doesn't need to be cleared
    magma_tally4blas_claset( Magma_tally4Full, n, nb, c_zero, c_zero, dwork, lddw );

    if (upper) {
        /* Copy the matrix to the GPU */
        magma_tally4_csetmatrix( n, n, A(0, 0), lda, dA(0, 0), ldda );

        /* Reduce the upper triangle of A.
           Columns 1:kk are handled by the unblocked method. */
        kk = n - (n - nx + nb - 1) / nb * nb;

        for (i = n - nb; i >= kk; i -= nb) {
            /* Reduce columns i:i+nb-1 to tridiagonal form and form the
               matrix W which is needed to update the unreduced part of
               the matrix */
            
            /* Get the current panel (no need for the 1st iteration) */
            if (i != n-nb)
                magma_tally4_cgetmatrix( i+nb, nb, dA(0, i), ldda, A(0, i), lda );
            
            #ifdef FAST_HEMV
            magma_tally4_clatrd2( uplo, i+nb, nb, A(0, 0), lda, e, tau,
                           work, ldw, dA(0, 0), ldda, dwork, lddw,
                           dwork2, ldwork2 );
            #else
            magma_tally4_clatrd(  uplo, i+nb, nb, A(0, 0), lda, e, tau,
                           work, ldw, dA(0, 0), ldda, dwork, lddw );
            #endif

            /* Update the unreduced submatrix A(0:i-2,0:i-2), using an
               update of the form:  A := A - V*W' - W*V' */
            magma_tally4_csetmatrix( i + nb, nb, work, ldw, dwork, lddw );

            magma_tally4_cher2k( uplo, Magma_tally4NoTrans, i, nb, c_neg_one,
                          dA(0, i), ldda, dwork, lddw,
                          d_one, dA(0, 0), ldda );
            
            /* Copy superdiagonal elements back into A, and diagonal
               elements into D */
            for (j = i; j < i+nb; ++j) {
                *A(j-1,j) = MAGMA_tally4_C_MAKE( e[j - 1], 0 );
                d[j] = MAGMA_tally4_C_REAL( *A(j, j) );
            }
        }
        
        magma_tally4_cgetmatrix( kk, kk, dA(0, 0), ldda, A(0, 0), lda );
        
        /* Use CPU code to reduce the last or only block */
        lapackf77_chetrd( uplo_, &kk, A(0, 0), &lda, d, e, tau, work, &lwork, &iinfo );
    }
    else {
        /* Copy the matrix to the GPU */
        if (1 <= n-nx)
            magma_tally4_csetmatrix( n, n, A(0,0), lda, dA(0,0), ldda );

        /* Reduce the lower triangle of A */
        for (i = 0; i < n-nx; i += nb) {
            /* Reduce columns i:i+nb-1 to tridiagonal form and form the
               matrix W which is needed to update the unreduced part of
               the matrix */

            /* Get the current panel (no need for the 1st iteration) */
            if (i != 0)
                magma_tally4_cgetmatrix( n-i, nb, dA(i, i), ldda, A(i, i), lda );
            
            #ifdef FAST_HEMV
            magma_tally4_clatrd2( uplo, n-i, nb, A(i, i), lda, &e[i], &tau[i],
                           work, ldw, dA(i, i), ldda, dwork, lddw,
                           dwork2, ldwork2 );
            #else
            magma_tally4_clatrd(  uplo, n-i, nb, A(i, i), lda, &e[i], &tau[i],
                           work, ldw, dA(i, i), ldda, dwork, lddw );
            #endif
            
            /* Update the unreduced submatrix A(i+ib:n,i+ib:n), using
               an update of the form:  A := A - V*W' - W*V' */
            magma_tally4_csetmatrix( n-i, nb, work, ldw, dwork, lddw );

            magma_tally4_cher2k( Magma_tally4Lower, Magma_tally4NoTrans, n-i-nb, nb, c_neg_one,
                          dA(i+nb, i), ldda, &dwork[nb], lddw,
                          d_one, dA(i+nb, i+nb), ldda );
            
            /* Copy subdiagonal elements back into A, and diagonal
               elements into D */
            for (j = i; j < i+nb; ++j) {
                *A(j+1,j) = MAGMA_tally4_C_MAKE( e[j], 0 );
                d[j] = MAGMA_tally4_C_REAL( *A(j, j) );
            }
        }

        /* Use CPU code to reduce the last or only block */
        if (1 <= n-nx)
            magma_tally4_cgetmatrix( n-i, n-i, dA(i, i), ldda, A(i, i), lda );
        
        i_n = n-i;
        lapackf77_chetrd( uplo_, &i_n, A(i, i), &lda, &d[i], &e[i],
                          &tau[i], work, &lwork, &iinfo );
    }
    
    magma_tally4_free( dA );
    
    work[0] = MAGMA_tally4_C_MAKE( lwkopt, 0 );

    return *info;
} /* magma_tally4_chetrd */
