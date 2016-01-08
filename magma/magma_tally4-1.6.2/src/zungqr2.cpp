/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

       @author Stan Tomov
       @author Mark Gates
*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    ZUNGQR generates an M-by-N COMPLEX_16 matrix Q with orthonormal columns,
    which is defined as the first N columns of a product of K elementary
    reflectors of order M

          Q  =  H(1) H(2) . . . H(k)

    as returned by ZGEQRF.

    This version recomputes the T matrices on the CPU and sends them to the GPU.

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
    A       COMPLEX_16 array A, dimension (LDDA,N).
            On entry, the i-th column must contain the vector
            which defines the elementary reflector H(i), for
            i = 1,2,...,k, as returned by ZGEQRF_GPU in the
            first k columns of its array argument A.
            On exit, the M-by-N matrix Q.

    @param[in]
    lda     INTEGER
            The first dimension of the array A. LDA >= max(1,M).

    @param[in]
    tau     COMPLEX_16 array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by ZGEQRF_GPU.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument has an illegal value

    @ingroup magma_tally4_zgeqrf_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_zungqr2(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4_int_t *info)
{
#define  A(i,j) ( A + (i) + (j)*lda )
#define dA(i,j) (dA + (i) + (j)*ldda)

    magma_tally4DoubleComplex c_zero = MAGMA_tally4_Z_ZERO;
    magma_tally4DoubleComplex c_one  = MAGMA_tally4_Z_ONE;

    magma_tally4_int_t nb = magma_tally4_get_zgeqrf_nb(min(m, n));

    magma_tally4_int_t  m_kk, n_kk, k_kk, mi;
    magma_tally4_int_t lwork, ldda;
    magma_tally4_int_t i, ib, ki, kk;  //, iinfo;
    magma_tally4_int_t lddwork;
    magma_tally4DoubleComplex *dA, *dV, *dW, *dT, *T;
    magma_tally4DoubleComplex *work;

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
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }

    if (n <= 0) {
        return *info;
    }

    magma_tally4_queue_t orig_stream;
    magma_tally4blasGetKernelStream( &orig_stream );
    
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
    if (MAGMA_tally4_SUCCESS != magma_tally4_zmalloc( &dA, ldda*n + ldda*nb + lddwork*nb + nb*nb)) {
        *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
        return *info;
    }
    dV = dA + ldda*n;
    dW = dA + ldda*n + ldda*nb;
    dT = dA + ldda*n + ldda*nb + lddwork*nb;

    // Allocate CPU work space
    lwork = (n+m+nb) * nb;
    magma_tally4_zmalloc_cpu( &work, lwork );

    T = work;

    if (work == NULL) {
        magma_tally4_free( dA );
        magma_tally4_free_cpu( work );
        *info = MAGMA_tally4_ERR_HOST_ALLOC;
        return *info;
    }
    magma_tally4DoubleComplex *V = work + (n+nb)*nb;

    magma_tally4_queue_t stream;
    magma_tally4_queue_create( &stream );

    // Use unblocked code for the last or only block.
    if (kk < n) {
        m_kk = m - kk;
        n_kk = n - kk;
        k_kk = k - kk;
        /*
            lapackf77_zungqr( &m_kk, &n_kk, &k_kk,
                              A(kk, kk), &lda,
                              &tau[kk], work, &lwork, &iinfo );
        */
        lapackf77_zlacpy( Magma_tally4UpperLowerStr, &m_kk, &k_kk, A(kk,kk), &lda, V, &m_kk);
        lapackf77_zlaset( Magma_tally4UpperLowerStr, &m_kk, &n_kk, &c_zero, &c_one, A(kk, kk), &lda );

        lapackf77_zlarft( Magma_tally4ForwardStr, Magma_tally4ColumnwiseStr,
                          &m_kk, &k_kk,
                          V, &m_kk, &tau[kk], work, &k_kk);
        lapackf77_zlarfb( Magma_tally4LeftStr, Magma_tally4NoTransStr, Magma_tally4ForwardStr, Magma_tally4ColumnwiseStr,
                          &m_kk, &n_kk, &k_kk,
                          V, &m_kk, work, &k_kk, A(kk, kk), &lda, work+k_kk*k_kk, &n_kk );
        
        if (kk > 0) {
            magma_tally4_zsetmatrix( m_kk, n_kk,
                              A(kk, kk),  lda,
                              dA(kk, kk), ldda );
        
            // Set A(1:kk,kk+1:n) to zero.
            magma_tally4blas_zlaset( Magma_tally4Full, kk, n - kk, c_zero, c_zero, dA(0, kk), ldda );
        }
    }

    if (kk > 0) {
        // Use blocked code
        // stream: set Aii (V) --> laset --> laset --> larfb --> [next]
        // CPU has no computation
        magma_tally4blasSetKernelStream( stream );
        
        for (i = ki; i >= 0; i -= nb) {
            ib = min(nb, k - i);

            // Send current panel to the GPU
            mi = m - i;
            lapackf77_zlaset( "Upper", &ib, &ib, &c_zero, &c_one, A(i, i), &lda );
            magma_tally4_zsetmatrix_async( mi, ib,
                                    A(i, i), lda,
                                    dV,      ldda, stream );
            lapackf77_zlarft( Magma_tally4ForwardStr, Magma_tally4ColumnwiseStr,
                              &mi, &ib,
                              A(i,i), &lda, &tau[i], T, &nb);
            magma_tally4_zsetmatrix_async( ib, ib,
                                    T, nb,
                                    dT, nb, stream );

            // set panel to identity
            magma_tally4blas_zlaset( Magma_tally4Full, i,  ib, c_zero, c_zero, dA(0, i), ldda );
            magma_tally4blas_zlaset( Magma_tally4Full, mi, ib, c_zero, c_one,  dA(i, i), ldda );
            
            magma_tally4_queue_sync( stream );
            if (i < n) {
                // Apply H to A(i:m,i:n) from the left
                magma_tally4_zlarfb_gpu( Magma_tally4Left, Magma_tally4NoTrans, Magma_tally4Forward, Magma_tally4Columnwise,
                                  mi, n-i, ib,
                                  dV,       ldda, dT, nb,
                                  dA(i, i), ldda, dW, lddwork );
            }
        }
    
        // copy result back to CPU
        magma_tally4_zgetmatrix( m, n,
                          dA(0, 0), ldda, A(0, 0), lda);
    }

    magma_tally4_queue_destroy( stream );
    magma_tally4_free( dA );
    magma_tally4_free_cpu( work );

    magma_tally4blasSetKernelStream( orig_stream );
    
    return *info;
} /* magma_tally4_zungqr */
