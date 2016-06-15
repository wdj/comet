/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zungqr.cpp normal z -> d, Fri Jan 30 19:00:16 2015

       @author Stan Tomov
       @author Mark Gates
*/
#include "common_magma_tally3.h"

/**
    Purpose
    -------
    DORGQR generates an M-by-N DOUBLE_PRECISION matrix Q with orthonormal columns,
    which is defined as the first N columns of a product of K elementary
    reflectors of order M

          Q  =  H(1) H(2) . . . H(k)

    as returned by DGEQRF.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix Q. M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix Q. M >= N >= 0.

    @param[in]
    k       INTEGER
            The number of elementary reflectors whose product defines the
            matrix Q. N >= K >= 0.

    @param[in,out]
    A       DOUBLE_PRECISION array A, dimension (LDDA,N).
            On entry, the i-th column must contain the vector
            which defines the elementary reflector H(i), for
            i = 1,2,...,k, as returned by DGEQRF_GPU in the
            first k columns of its array argument A.
            On exit, the M-by-N matrix Q.

    @param[in]
    lda     INTEGER
            The first dimension of the array A. LDA >= max(1,M).

    @param[in]
    tau     DOUBLE_PRECISION array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by DGEQRF_GPU.

    @param[in]
    dT      DOUBLE_PRECISION array on the GPU device.
            DT contains the T matrices used in blocking the elementary
            reflectors H(i), e.g., this can be the 6th argument of
            magma_tally3_dgeqrf_gpu.

    @param[in]
    nb      INTEGER
            This is the block size used in DGEQRF_GPU, and correspondingly
            the size of the T matrices, used in the factorization, and
            stored in DT.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument has an illegal value

    @ingroup magma_tally3_dgeqrf_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_dorgqr(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    double *A, magma_tally3_int_t lda,
    double *tau,
    magma_tally3Double_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info)
{
#define  A(i,j) ( A + (i) + (j)*lda )
#define dA(i,j) (dA + (i) + (j)*ldda)
#define dT(j)   (dT + (j)*nb)

    double c_zero = MAGMA_tally3_D_ZERO;
    double c_one  = MAGMA_tally3_D_ONE;

    magma_tally3_int_t  m_kk, n_kk, k_kk, mi;
    magma_tally3_int_t lwork, ldda;
    magma_tally3_int_t i, ib, ki, kk;  //, iinfo;
    magma_tally3_int_t lddwork;
    double *dA, *dV, *dW;
    double *work;

    *info = 0;
    if (m < 0) {
        *info = -1;
    } else if ((n < 0) || (n > m)) {
        *info = -2;
    } else if ((k < 0) || (k > n)) {
        *info = -3;
    } else if (lda < max(1,m)) {
        *info = -5;
    }
    if (*info != 0) {
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }

    if (n <= 0) {
        return *info;
    }

    magma_tally3_queue_t orig_stream;
    magma_tally3blasGetKernelStream( &orig_stream );
    
    // first kk columns are handled by blocked method.
    // ki is start of 2nd-to-last block
    if ((nb > 1) && (nb < k)) {
        ki = (k - nb - 1) / nb * nb;
        kk = min(k, ki + nb);
    } else {
        ki = 0;
        kk = 0;
    }

    // Allocate GPU work space
    // ldda*n     for matrix dA
    // ldda*nb    for dV
    // lddwork*nb for dW larfb workspace
    ldda    = ((m + 31) / 32) * 32;
    lddwork = ((n + 31) / 32) * 32;
    if (MAGMA_tally3_SUCCESS != magma_tally3_dmalloc( &dA, ldda*n + ldda*nb + lddwork*nb )) {
        *info = MAGMA_tally3_ERR_DEVICE_ALLOC;
        return *info;
    }
    dV = dA + ldda*n;
    dW = dA + ldda*n + ldda*nb;

    // Allocate CPU work space
    lwork = (n+m+nb) * nb;
    magma_tally3_dmalloc_cpu( &work, lwork );
    if (work == NULL) {
        magma_tally3_free( dA );
        *info = MAGMA_tally3_ERR_HOST_ALLOC;
        return *info;
    }
    double *V = work + (n+nb)*nb;

    magma_tally3_queue_t stream;
    magma_tally3_queue_create( &stream );

    // Use unblocked code for the last or only block.
    if (kk < n) {
        m_kk = m - kk;
        n_kk = n - kk;
        k_kk = k - kk;
        /*
            // Replacing this with the following 4 routines works but dorgqr is slow for
            // k smaller than the dorgqr's blocking size (new version can be up to 60x faster)
            lapackf77_dorgqr( &m_kk, &n_kk, &k_kk,
                              A(kk, kk), &lda,
                              &tau[kk], work, &lwork, &iinfo );
        */
        lapackf77_dlacpy( Magma_tally3UpperLowerStr, &m_kk, &k_kk, A(kk,kk), &lda, V, &m_kk);
        lapackf77_dlaset( Magma_tally3UpperLowerStr, &m_kk, &n_kk, &c_zero, &c_one, A(kk, kk), &lda );

        lapackf77_dlarft( Magma_tally3ForwardStr, Magma_tally3ColumnwiseStr,
                          &m_kk, &k_kk,
                          V, &m_kk, &tau[kk], work, &k_kk);
        lapackf77_dlarfb( Magma_tally3LeftStr, Magma_tally3NoTransStr, Magma_tally3ForwardStr, Magma_tally3ColumnwiseStr,
                          &m_kk, &n_kk, &k_kk,
                          V, &m_kk, work, &k_kk, A(kk, kk), &lda, work+k_kk*k_kk, &n_kk );
        
        if (kk > 0) {
            magma_tally3_dsetmatrix( m_kk, n_kk,
                              A(kk, kk),  lda,
                              dA(kk, kk), ldda );
        
            // Set A(1:kk,kk+1:n) to zero.
            magma_tally3blas_dlaset( Magma_tally3Full, kk, n - kk, c_zero, c_zero, dA(0, kk), ldda );
        }
    }

    if (kk > 0) {
        // Use blocked code
        // stream: set Aii (V) --> laset --> laset --> larfb --> [next]
        // CPU has no computation
        magma_tally3blasSetKernelStream( stream );
        
        for (i = ki; i >= 0; i -= nb) {
            ib = min(nb, k - i);

            // Send current panel to the GPU
            mi = m - i;
            lapackf77_dlaset( "Upper", &ib, &ib, &c_zero, &c_one, A(i, i), &lda );
            magma_tally3_dsetmatrix_async( mi, ib,
                                    A(i, i), lda,
                                    dV,      ldda, stream );

            // set panel to identity
            magma_tally3blas_dlaset( Magma_tally3Full, i,  ib, c_zero, c_zero, dA(0, i), ldda );
            magma_tally3blas_dlaset( Magma_tally3Full, mi, ib, c_zero, c_one,  dA(i, i), ldda );
            
            if (i < n) {
                // Apply H to A(i:m,i:n) from the left
                magma_tally3_dlarfb_gpu( Magma_tally3Left, Magma_tally3NoTrans, Magma_tally3Forward, Magma_tally3Columnwise,
                                  mi, n-i, ib,
                                  dV,       ldda, dT(i), nb,
                                  dA(i, i), ldda, dW, lddwork );
            }
        }
    
        // copy result back to CPU
        magma_tally3_dgetmatrix( m, n,
                          dA(0, 0), ldda, A(0, 0), lda);
    }

    magma_tally3_queue_destroy( stream );
    magma_tally3_free( dA );
    magma_tally3_free_cpu( work );

    magma_tally3blasSetKernelStream( orig_stream );
    
    return *info;
} /* magma_tally3_dorgqr */