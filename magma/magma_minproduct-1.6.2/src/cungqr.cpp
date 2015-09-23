/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zungqr.cpp normal z -> c, Fri Jan 30 19:00:15 2015

       @author Stan Tomov
       @author Mark Gates
*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    CUNGQR generates an M-by-N COMPLEX matrix Q with orthonormal columns,
    which is defined as the first N columns of a product of K elementary
    reflectors of order M

          Q  =  H(1) H(2) . . . H(k)

    as returned by CGEQRF.

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
    A       COMPLEX array A, dimension (LDDA,N).
            On entry, the i-th column must contain the vector
            which defines the elementary reflector H(i), for
            i = 1,2,...,k, as returned by CGEQRF_GPU in the
            first k columns of its array argument A.
            On exit, the M-by-N matrix Q.

    @param[in]
    lda     INTEGER
            The first dimension of the array A. LDA >= max(1,M).

    @param[in]
    tau     COMPLEX array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by CGEQRF_GPU.

    @param[in]
    dT      COMPLEX array on the GPU device.
            DT contains the T matrices used in blocking the elementary
            reflectors H(i), e.g., this can be the 6th argument of
            magma_minproduct_cgeqrf_gpu.

    @param[in]
    nb      INTEGER
            This is the block size used in CGEQRF_GPU, and correspondingly
            the size of the T matrices, used in the factorization, and
            stored in DT.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument has an illegal value

    @ingroup magma_minproduct_cgeqrf_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_cungqr(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info)
{
#define  A(i,j) ( A + (i) + (j)*lda )
#define dA(i,j) (dA + (i) + (j)*ldda)
#define dT(j)   (dT + (j)*nb)

    magma_minproductFloatComplex c_zero = MAGMA_minproduct_C_ZERO;
    magma_minproductFloatComplex c_one  = MAGMA_minproduct_C_ONE;

    magma_minproduct_int_t  m_kk, n_kk, k_kk, mi;
    magma_minproduct_int_t lwork, ldda;
    magma_minproduct_int_t i, ib, ki, kk;  //, iinfo;
    magma_minproduct_int_t lddwork;
    magma_minproductFloatComplex *dA, *dV, *dW;
    magma_minproductFloatComplex *work;

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
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }

    if (n <= 0) {
        return *info;
    }

    magma_minproduct_queue_t orig_stream;
    magma_minproductblasGetKernelStream( &orig_stream );
    
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
    if (MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc( &dA, ldda*n + ldda*nb + lddwork*nb )) {
        *info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
        return *info;
    }
    dV = dA + ldda*n;
    dW = dA + ldda*n + ldda*nb;

    // Allocate CPU work space
    lwork = (n+m+nb) * nb;
    magma_minproduct_cmalloc_cpu( &work, lwork );
    if (work == NULL) {
        magma_minproduct_free( dA );
        *info = MAGMA_minproduct_ERR_HOST_ALLOC;
        return *info;
    }
    magma_minproductFloatComplex *V = work + (n+nb)*nb;

    magma_minproduct_queue_t stream;
    magma_minproduct_queue_create( &stream );

    // Use unblocked code for the last or only block.
    if (kk < n) {
        m_kk = m - kk;
        n_kk = n - kk;
        k_kk = k - kk;
        /*
            // Replacing this with the following 4 routines works but cungqr is slow for
            // k smaller than the cungqr's blocking size (new version can be up to 60x faster)
            lapackf77_cungqr( &m_kk, &n_kk, &k_kk,
                              A(kk, kk), &lda,
                              &tau[kk], work, &lwork, &iinfo );
        */
        lapackf77_clacpy( Magma_minproductUpperLowerStr, &m_kk, &k_kk, A(kk,kk), &lda, V, &m_kk);
        lapackf77_claset( Magma_minproductUpperLowerStr, &m_kk, &n_kk, &c_zero, &c_one, A(kk, kk), &lda );

        lapackf77_clarft( Magma_minproductForwardStr, Magma_minproductColumnwiseStr,
                          &m_kk, &k_kk,
                          V, &m_kk, &tau[kk], work, &k_kk);
        lapackf77_clarfb( Magma_minproductLeftStr, Magma_minproductNoTransStr, Magma_minproductForwardStr, Magma_minproductColumnwiseStr,
                          &m_kk, &n_kk, &k_kk,
                          V, &m_kk, work, &k_kk, A(kk, kk), &lda, work+k_kk*k_kk, &n_kk );
        
        if (kk > 0) {
            magma_minproduct_csetmatrix( m_kk, n_kk,
                              A(kk, kk),  lda,
                              dA(kk, kk), ldda );
        
            // Set A(1:kk,kk+1:n) to zero.
            magma_minproductblas_claset( Magma_minproductFull, kk, n - kk, c_zero, c_zero, dA(0, kk), ldda );
        }
    }

    if (kk > 0) {
        // Use blocked code
        // stream: set Aii (V) --> laset --> laset --> larfb --> [next]
        // CPU has no computation
        magma_minproductblasSetKernelStream( stream );
        
        for (i = ki; i >= 0; i -= nb) {
            ib = min(nb, k - i);

            // Send current panel to the GPU
            mi = m - i;
            lapackf77_claset( "Upper", &ib, &ib, &c_zero, &c_one, A(i, i), &lda );
            magma_minproduct_csetmatrix_async( mi, ib,
                                    A(i, i), lda,
                                    dV,      ldda, stream );

            // set panel to identity
            magma_minproductblas_claset( Magma_minproductFull, i,  ib, c_zero, c_zero, dA(0, i), ldda );
            magma_minproductblas_claset( Magma_minproductFull, mi, ib, c_zero, c_one,  dA(i, i), ldda );
            
            if (i < n) {
                // Apply H to A(i:m,i:n) from the left
                magma_minproduct_clarfb_gpu( Magma_minproductLeft, Magma_minproductNoTrans, Magma_minproductForward, Magma_minproductColumnwise,
                                  mi, n-i, ib,
                                  dV,       ldda, dT(i), nb,
                                  dA(i, i), ldda, dW, lddwork );
            }
        }
    
        // copy result back to CPU
        magma_minproduct_cgetmatrix( m, n,
                          dA(0, 0), ldda, A(0, 0), lda);
    }

    magma_minproduct_queue_destroy( stream );
    magma_minproduct_free( dA );
    magma_minproduct_free_cpu( work );

    magma_minproductblasSetKernelStream( orig_stream );
    
    return *info;
} /* magma_minproduct_cungqr */
