/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zungqr_gpu.cpp normal z -> s, Fri Jan 30 19:00:15 2015

       @author Stan Tomov
       @author Mark Gates
*/
#include "common_magma_tally3.h"

/**
    Purpose
    -------
    SORGQR generates an M-by-N REAL matrix Q with orthonormal columns,
    which is defined as the first N columns of a product of K elementary
    reflectors of order M

          Q  =  H(1) H(2) . . . H(k)

    as returned by SGEQRF_GPU.

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
    dA      REAL array A on the GPU, dimension (LDDA,N).
            On entry, the i-th column must contain the vector
            which defines the elementary reflector H(i), for
            i = 1,2,...,k, as returned by SGEQRF_GPU in the
            first k columns of its array argument A.
            On exit, the M-by-N matrix Q.

    @param[in]
    ldda    INTEGER
            The first dimension of the array A. LDDA >= max(1,M).

    @param[in]
    tau     REAL array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by SGEQRF_GPU.

    @param[in]
    dT      (workspace) REAL work space array on the GPU,
            dimension (2*MIN(M, N) + (N+31)/32*32 )*NB.
            This must be the 6th argument of magma_tally3_sgeqrf_gpu
            [ note that if N here is bigger than N in magma_tally3_sgeqrf_gpu,
              the workspace requirement DT in magma_tally3_sgeqrf_gpu must be
              as specified in this routine ].

    @param[in]
    nb      INTEGER
            This is the block size used in SGEQRF_GPU, and correspondingly
            the size of the T matrices, used in the factorization, and
            stored in DT.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument has an illegal value

    @ingroup magma_tally3_sgeqrf_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_sorgqr_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *tau,
    magma_tally3Float_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info)
{
#define dA(i,j) (dA + (i) + (j)*ldda)
#define dT(j)   (dT + (j)*nb)

    float c_zero = MAGMA_tally3_S_ZERO;
    float c_one  = MAGMA_tally3_S_ONE;
    
    magma_tally3_int_t m_kk, n_kk, k_kk, mi;
    magma_tally3_int_t lwork, lpanel;
    magma_tally3_int_t i, ib, ki, kk, iinfo;
    magma_tally3_int_t lddwork;
    magma_tally3Float_ptr dV, dW;
    float *work, *panel;

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
        kk = min( k, ki+nb );
    } else {
        ki = 0;
        kk = 0;
    }

    // Allocate CPU work space
    // n*nb for sorgqr workspace
    // (m - kk)*(n - kk) for last block's panel
    lwork = n*nb;
    lpanel = (m - kk)*(n - kk);
    magma_tally3_smalloc_cpu( &work, lwork + lpanel );
    if ( work == NULL ) {
        *info = MAGMA_tally3_ERR_HOST_ALLOC;
        return *info;
    }
    panel = work + lwork;
    
    // Allocate work space on GPU
    if (MAGMA_tally3_SUCCESS != magma_tally3_smalloc( &dV, ldda*nb )) {
        magma_tally3_free_cpu( work );
        *info = MAGMA_tally3_ERR_DEVICE_ALLOC;
        return *info;
    }
    
    // dT workspace has:
    // 2*min(m,n)*nb      for T and R^{-1} matrices from geqrf
    // ((n+31)/32*32 )*nb for dW larfb workspace.
    lddwork = min(m,n);
    dW = dT + 2*lddwork*nb;

    magma_tally3_queue_t stream;
    magma_tally3_queue_create( &stream );

    // Use unblocked code for the last or only block.
    if (kk < n) {
        m_kk = m - kk;
        n_kk = n - kk;
        k_kk = k - kk;
        magma_tally3_sgetmatrix( m_kk, k_kk,
                          dA(kk, kk), ldda, panel, m_kk );
        
        lapackf77_sorgqr( &m_kk, &n_kk, &k_kk,
                          panel, &m_kk,
                          &tau[kk], work, &lwork, &iinfo );
        
        magma_tally3_ssetmatrix( m_kk, n_kk,
                          panel, m_kk, dA(kk, kk), ldda );
        
        // Set A(1:kk,kk+1:n) to zero.
        magma_tally3blas_slaset( Magma_tally3Full, kk, n - kk, c_zero, c_zero, dA(0, kk), ldda );
    }

    if (kk > 0) {
        // Use blocked code
        // stream:  copy Aii to V --> laset --> laset --> larfb --> [next]
        // CPU has no computation
        magma_tally3blasSetKernelStream( stream );
        
        for (i = ki; i >= 0; i -= nb) {
            ib = min( nb, k-i );
            mi = m - i;
            
            // Copy current panel on the GPU from dA to dV
            magma_tally3_scopymatrix_async( mi, ib,
                                     dA(i,i), ldda,
                                     dV,      ldda, stream );

            // set panel to identity
            magma_tally3blas_slaset( Magma_tally3Full, i,  ib, c_zero, c_zero, dA(0, i), ldda );
            magma_tally3blas_slaset( Magma_tally3Full, mi, ib, c_zero, c_one,  dA(i, i), ldda );
            
            if (i < n) {
                // Apply H to A(i:m,i:n) from the left
                magma_tally3_slarfb_gpu( Magma_tally3Left, Magma_tally3NoTrans, Magma_tally3Forward, Magma_tally3Columnwise,
                                  mi, n-i, ib,
                                  dV,       ldda, dT(i), nb,
                                  dA(i, i), ldda, dW, lddwork );
            }
        }
    }
    magma_tally3_queue_sync( stream );

    magma_tally3_free( dV );
    magma_tally3_free_cpu( work );
    magma_tally3_queue_destroy( stream );
    
    magma_tally3blasSetKernelStream( orig_stream );

    return *info;
} /* magma_tally3_sorgqr_gpu */
