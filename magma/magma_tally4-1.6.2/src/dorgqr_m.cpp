/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zungqr_m.cpp normal z -> d, Fri Jan 30 19:00:16 2015

       @author Mark Gates
*/
#include "common_magma_tally4.h"
#include "trace.h"

#define PRECISION_d

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
    T       DOUBLE_PRECISION array, dimension (NB, min(M,N)).
            T contains the T matrices used in blocking the elementary
            reflectors H(i), e.g., this can be the 6th argument of
            magma_tally4_dgeqrf_gpu (except stored on the CPU, not the GPU).

    @param[in]
    nb      INTEGER
            This is the block size used in DGEQRF_GPU, and correspondingly
            the size of the T matrices, used in the factorization, and
            stored in T.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument has an illegal value

    @ingroup magma_tally4_dgeqrf_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_dorgqr_m(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    double *A, magma_tally4_int_t lda,
    double *tau,
    double *T, magma_tally4_int_t nb,
    magma_tally4_int_t *info)
{
#define  A(i,j)   ( A    + (i) + (j)*lda )
#define dA(d,i,j) (dA[d] + (i) + (j)*ldda)
#define dT(d,i,j) (dT[d] + (i) + (j)*nb)

    double c_zero = MAGMA_tally4_D_ZERO;
    double c_one  = MAGMA_tally4_D_ONE;

    magma_tally4_int_t m_kk, n_kk, k_kk, mi;
    magma_tally4_int_t lwork, ldwork;
    magma_tally4_int_t i, ib, ki, kk, iinfo;
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
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }

    if (n <= 0) {
        return *info;
    }
    
    magma_tally4_int_t di, dn;
    magma_tally4_int_t dpanel;

    magma_tally4_int_t ngpu = magma_tally4_num_gpus();
    
    magma_tally4_device_t orig_dev;
    magma_tally4_getdevice( &orig_dev );
    magma_tally4_queue_t orig_stream;
    magma_tally4blasGetKernelStream( &orig_stream );
    
    // Allocate memory on GPUs for A and workspaces
    magma_tally4_int_t ldda    = ((m + 31) / 32) * 32;
    magma_tally4_int_t lddwork = ((n + 31) / 32) * 32;
    magma_tally4_int_t min_lblocks = (n / nb) / ngpu;  // min. blocks per gpu
    magma_tally4_int_t last_dev    = (n / nb) % ngpu;  // device with last block
    
    magma_tally4_int_t  nlocal[ Magma_tally4MaxGPUs ] = { 0 };
    double *dA[ Magma_tally4MaxGPUs ] = { NULL };
    double *dT[ Magma_tally4MaxGPUs ] = { NULL };
    double *dV[ Magma_tally4MaxGPUs ] = { NULL };
    double *dW[ Magma_tally4MaxGPUs ] = { NULL };
    magma_tally4_queue_t stream[ Magma_tally4MaxGPUs ] = { NULL };
    
    for( int d = 0; d < ngpu; ++d ) {
        // example with n = 75, nb = 10, ngpu = 3
        // min_lblocks = 2
        // last_dev    = 1
        // gpu 0: 2  blocks, cols:  0- 9, 30-39, 60-69
        // gpu 1: 1+ blocks, cols: 10-19, 40-49, 70-74 (partial)
        // gpu 2: 1  block,  cols: 20-29, 50-59
        magma_tally4_setdevice( d );
        nlocal[d] = min_lblocks*nb;
        if ( d < last_dev ) {
            nlocal[d] += nb;
        }
        else if ( d == last_dev ) {
            nlocal[d] += (n % nb);
        }
        
        ldwork = nlocal[d]*ldda  // dA
               + nb*m            // dT
               + nb*ldda         // dV
               + nb*lddwork;     // dW
        if ( MAGMA_tally4_SUCCESS != magma_tally4_dmalloc( &dA[d], ldwork )) {
            *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
            goto CLEANUP;
        }
        dT[d] = dA[d] + nlocal[d]*ldda;
        dV[d] = dT[d] + nb*m;
        dW[d] = dV[d] + nb*ldda;
        
        magma_tally4_queue_create( &stream[d] );
    }
    
    trace_init( 1, ngpu, 1, stream );
    
    // first kk columns are handled by blocked method.
    // ki is start of 2nd-to-last block
    if ((nb > 1) && (nb < k)) {
        ki = (k - nb - 1) / nb * nb;
        kk = min(k, ki + nb);
    } else {
        ki = 0;
        kk = 0;
    }

    // Allocate CPU work space
    // n*nb for dorgqr workspace
    lwork = n * nb;
    magma_tally4_dmalloc_cpu( &work, lwork );
    if (work == NULL) {
        *info = MAGMA_tally4_ERR_HOST_ALLOC;
        goto CLEANUP;
    }

    // Use unblocked code for the last or only block.
    if (kk < n) {
        trace_cpu_start( 0, "ungqr", "ungqr last block" );
        m_kk = m - kk;
        n_kk = n - kk;
        k_kk = k - kk;
        dpanel =  (kk / nb) % ngpu;
        di     = ((kk / nb) / ngpu) * nb;
        magma_tally4_setdevice( dpanel );
        
        lapackf77_dorgqr( &m_kk, &n_kk, &k_kk,
                          A(kk, kk), &lda,
                          &tau[kk], work, &lwork, &iinfo );

        magma_tally4_dsetmatrix( m_kk, n_kk,
                          A(kk, kk),  lda,
                          dA(dpanel, kk, di), ldda );
        
        // Set A(1:kk,kk+1:n) to zero.
        magma_tally4blas_dlaset( Magma_tally4Full, kk, n - kk, c_zero, c_zero, dA(dpanel, 0, di), ldda );
        trace_cpu_end( 0 );
    }

    if (kk > 0) {
        // Use blocked code
        // send T to all GPUs
        for( int d = 0; d < ngpu; ++d ) {
            magma_tally4_setdevice( d );
            trace_gpu_start( d, 0, "set", "set T" );
            magma_tally4_dsetmatrix_async( nb, min(m,n), T, nb, dT[d], nb, stream[d] );
            trace_gpu_end( d, 0 );
        }
        
        // stream: set Aii (V) --> laset --> laset --> larfb --> [next]
        // CPU has no computation
        for( i = ki; i >= 0; i -= nb ) {
            ib = min(nb, k - i);
            mi = m - i;
            dpanel =  (i / nb) % ngpu;
            di     = ((i / nb) / ngpu) * nb;

            // Send current panel to the GPUs
            lapackf77_dlaset( "Upper", &ib, &ib, &c_zero, &c_one, A(i, i), &lda );
            for( int d = 0; d < ngpu; ++d ) {
                magma_tally4_setdevice( d );
                trace_gpu_start( d, 0, "set", "set V" );
                magma_tally4_dsetmatrix_async( mi, ib,
                                        A(i, i), lda,
                                        dV[d],   ldda, stream[d] );
                trace_gpu_end( d, 0 );
            }
            
            // set panel to identity
            magma_tally4_setdevice( dpanel );
            magma_tally4blasSetKernelStream( stream[dpanel] );
            trace_gpu_start( dpanel, 0, "laset", "laset" );
            magma_tally4blas_dlaset( Magma_tally4Full, i,  ib, c_zero, c_zero, dA(dpanel, 0, di), ldda );
            magma_tally4blas_dlaset( Magma_tally4Full, mi, ib, c_zero, c_one,  dA(dpanel, i, di), ldda );
            trace_gpu_end( dpanel, 0 );
            
            if (i < n) {
                // Apply H to A(i:m,i:n) from the left
                for( int d = 0; d < ngpu; ++d ) {
                    magma_tally4_setdevice( d );
                    magma_tally4blasSetKernelStream( stream[d] );
                    magma_tally4_indices_1D_bcyclic( nb, ngpu, d, i, n, &di, &dn );
                    trace_gpu_start( d, 0, "larfb", "larfb" );
                    magma_tally4_dlarfb_gpu( Magma_tally4Left, Magma_tally4NoTrans, Magma_tally4Forward, Magma_tally4Columnwise,
                                      mi, dn-di, ib,
                                      dV[d],        ldda, dT(d,0,i), nb,
                                      dA(d, i, di), ldda, dW[d], lddwork );
                    trace_gpu_end( d, 0 );
                }
            }
        }
    }
    
    // copy result back to CPU
    trace_cpu_start( 0, "get", "get A" );
    magma_tally4_dgetmatrix_1D_col_bcyclic( m, n, dA, ldda, A, lda, ngpu, nb );
    trace_cpu_end( 0 );
    
    #ifdef TRACING
    char name[80];
    snprintf( name, sizeof(name), "dorgqr-n%d-ngpu%d.svg", m, ngpu );
    trace_finalize( name, "trace.css" );
    #endif
    
CLEANUP:
    for( int d = 0; d < ngpu; ++d ) {
        magma_tally4_setdevice( d );
        magma_tally4_free( dA[d] );
        magma_tally4_queue_destroy( stream[d] );
    }
    magma_tally4_free_cpu( work );
    magma_tally4_setdevice( orig_dev );
    magma_tally4blasSetKernelStream( orig_stream );
    
    return *info;
} /* magma_tally4_dorgqr */
