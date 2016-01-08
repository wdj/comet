/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zungqr_gpu.cpp normal z -> d, Fri Jan 30 19:00:15 2015

       @author Stan Tomov
       @author Mark Gates
*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    DORGQR generates an M-by-N DOUBLE_PRECISION matrix Q with orthonormal columns,
    which is defined as the first N columns of a product of K elementary
    reflectors of order M

          Q  =  H(1) H(2) . . . H(k)

    as returned by DGEQRF_GPU.

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
    dA      DOUBLE_PRECISION array A on the GPU, dimension (LDDA,N).
            On entry, the i-th column must contain the vector
            which defines the elementary reflector H(i), for
            i = 1,2,...,k, as returned by DGEQRF_GPU in the
            first k columns of its array argument A.
            On exit, the M-by-N matrix Q.

    @param[in]
    ldda    INTEGER
            The first dimension of the array A. LDDA >= max(1,M).

    @param[in]
    tau     DOUBLE_PRECISION array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by DGEQRF_GPU.

    @param[in]
    dT      (workspace) DOUBLE_PRECISION work space array on the GPU,
            dimension (2*MIN(M, N) + (N+31)/32*32 )*NB.
            This must be the 6th argument of magma_tally4_dgeqrf_gpu
            [ note that if N here is bigger than N in magma_tally4_dgeqrf_gpu,
              the workspace requirement DT in magma_tally4_dgeqrf_gpu must be
              as specified in this routine ].

    @param[in]
    nb      INTEGER
            This is the block size used in DGEQRF_GPU, and correspondingly
            the size of the T matrices, used in the factorization, and
            stored in DT.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument has an illegal value

    @ingroup magma_tally4_dgeqrf_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_dorgqr_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *tau,
    magma_tally4Double_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info)
{
#define dA(i,j) (dA + (i) + (j)*ldda)
#define dT(j)   (dT + (j)*nb)

    double c_zero = MAGMA_tally4_D_ZERO;
    double c_one  = MAGMA_tally4_D_ONE;
    
    magma_tally4_int_t m_kk, n_kk, k_kk, mi;
    magma_tally4_int_t lwork, lpanel;
    magma_tally4_int_t i, ib, ki, kk, iinfo;
    magma_tally4_int_t lddwork;
    magma_tally4Double_ptr dV, dW;
    double *work, *panel;

    *info = 0;
    if (m < 0) {
        *info = -1;
    } else if ((n < 0) || (n > m)) {
        *info = -2;
    } else if ((k < 0) || (k > n)) {
        *info = -3;
    } else if (ldda < max(1,m)) {
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
        kk = min( k, ki+nb );
    } else {
        ki = 0;
        kk = 0;
    }

    // Allocate CPU work space
    // n*nb for dorgqr workspace
    // (m - kk)*(n - kk) for last block's panel
    lwork = n*nb;
    lpanel = (m - kk)*(n - kk);
    magma_tally4_dmalloc_cpu( &work, lwork + lpanel );
    if ( work == NULL ) {
        *info = MAGMA_tally4_ERR_HOST_ALLOC;
        return *info;
    }
    panel = work + lwork;
    
    // Allocate work space on GPU
    if (MAGMA_tally4_SUCCESS != magma_tally4_dmalloc( &dV, ldda*nb )) {
        magma_tally4_free_cpu( work );
        *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
        return *info;
    }
    
    // dT workspace has:
    // 2*min(m,n)*nb      for T and R^{-1} matrices from geqrf
    // ((n+31)/32*32 )*nb for dW larfb workspace.
    lddwork = min(m,n);
    dW = dT + 2*lddwork*nb;

    magma_tally4_queue_t stream;
    magma_tally4_queue_create( &stream );

    // Use unblocked code for the last or only block.
    if (kk < n) {
        m_kk = m - kk;
        n_kk = n - kk;
        k_kk = k - kk;
        magma_tally4_dgetmatrix( m_kk, k_kk,
                          dA(kk, kk), ldda, panel, m_kk );
        
        lapackf77_dorgqr( &m_kk, &n_kk, &k_kk,
                          panel, &m_kk,
                          &tau[kk], work, &lwork, &iinfo );
        
        magma_tally4_dsetmatrix( m_kk, n_kk,
                          panel, m_kk, dA(kk, kk), ldda );
        
        // Set A(1:kk,kk+1:n) to zero.
        magma_tally4blas_dlaset( Magma_tally4Full, kk, n - kk, c_zero, c_zero, dA(0, kk), ldda );
    }

    if (kk > 0) {
        // Use blocked code
        // stream:  copy Aii to V --> laset --> laset --> larfb --> [next]
        // CPU has no computation
        magma_tally4blasSetKernelStream( stream );
        
        for (i = ki; i >= 0; i -= nb) {
            ib = min( nb, k-i );
            mi = m - i;
            
            // Copy current panel on the GPU from dA to dV
            magma_tally4_dcopymatrix_async( mi, ib,
                                     dA(i,i), ldda,
                                     dV,      ldda, stream );

            // set panel to identity
            magma_tally4blas_dlaset( Magma_tally4Full, i,  ib, c_zero, c_zero, dA(0, i), ldda );
            magma_tally4blas_dlaset( Magma_tally4Full, mi, ib, c_zero, c_one,  dA(i, i), ldda );
            
            if (i < n) {
                // Apply H to A(i:m,i:n) from the left
                magma_tally4_dlarfb_gpu( Magma_tally4Left, Magma_tally4NoTrans, Magma_tally4Forward, Magma_tally4Columnwise,
                                  mi, n-i, ib,
                                  dV,       ldda, dT(i), nb,
                                  dA(i, i), ldda, dW, lddwork );
            }
        }
    }
    magma_tally4_queue_sync( stream );

    magma_tally4_free( dV );
    magma_tally4_free_cpu( work );
    magma_tally4_queue_destroy( stream );
    
    magma_tally4blasSetKernelStream( orig_stream );

    return *info;
} /* magma_tally4_dorgqr_gpu */
