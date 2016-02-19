/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Azzam Haidar

       @generated from zgetrs_batched.cpp normal z -> c, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_tally3.h"
#include "batched_kernel_param.h"
/**
    Purpose
    -------
    Solves a system of linear equations
      A * X = B,  A**T * X = B,  or  A**H * X = B
    with a general N-by-N matrix A using the LU factorization computed by CGETRF_GPU.

    Arguments
    ---------
    @param[in]
    trans   magma_tally3_trans_t
            Specifies the form of the system of equations:
      -     = Magma_tally3NoTrans:    A    * X = B  (No transpose)
      -     = Magma_tally3Trans:      A**T * X = B  (Transpose)
      -     = Magma_tally3ConjTrans:  A**H * X = B  (Conjugate transpose)

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in]
    dA      COMPLEX array on the GPU, dimension (LDA,N)
            The factors L and U from the factorization A = P*L*U as computed
            by CGETRF_GPU.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[in]
    ipiv    INTEGER array, dimension (N)
            The pivot indices from CGETRF; for 1 <= i <= N, row i of the
            matrix was interchanged with row IPIV(i).

    @param[in,out]
    dB      COMPLEX array on the GPU, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally3_cgesv_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_cgetrs_batched(
                  magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  magma_tally3FloatComplex **dA_array, magma_tally3_int_t ldda,
                  magma_tally3_int_t **dipiv_array, 
                  magma_tally3FloatComplex **dB_array, magma_tally3_int_t lddb,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue)
{
    magma_tally3_int_t notran = (trans == Magma_tally3NoTrans);
    magma_tally3_int_t info = 0;
    if ( (! notran) &&
         (trans != Magma_tally3Trans) &&
         (trans != Magma_tally3ConjTrans) ) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (nrhs < 0) {
        info = -3;
    } else if (ldda < max(1,n)) {
        info = -5;
    } else if (lddb < max(1,n)) {
        info = -8;
    }
    if (info != 0) {
        magma_tally3_xerbla( __func__, -(info) );
        return info;
    }

    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return info;
    }



    magma_tally3FloatComplex **dW1_displ  = NULL;
    magma_tally3FloatComplex **dW2_displ  = NULL;
    magma_tally3FloatComplex **dW3_displ  = NULL;
    magma_tally3FloatComplex **dW4_displ  = NULL;
    magma_tally3FloatComplex **dinvA_array = NULL;
    magma_tally3FloatComplex **dwork_array = NULL;

    magma_tally3_malloc((void**)&dW1_displ,  batchCount * sizeof(*dW1_displ));
    magma_tally3_malloc((void**)&dW2_displ,  batchCount * sizeof(*dW2_displ));
    magma_tally3_malloc((void**)&dW3_displ,  batchCount * sizeof(*dW3_displ));
    magma_tally3_malloc((void**)&dW4_displ,  batchCount * sizeof(*dW4_displ));
    magma_tally3_malloc((void**)&dinvA_array, batchCount * sizeof(*dinvA_array));
    magma_tally3_malloc((void**)&dwork_array, batchCount * sizeof(*dwork_array));

    magma_tally3_int_t invA_msize = ((n+TRI_NB-1)/TRI_NB)*TRI_NB*TRI_NB;
    magma_tally3_int_t dwork_msize = n*nrhs;
    magma_tally3FloatComplex* dinvA      = NULL;
    magma_tally3FloatComplex* dwork      = NULL;// dinvA and dwork are workspace in ctrsm
    magma_tally3_cmalloc( &dinvA, invA_msize * batchCount);
    magma_tally3_cmalloc( &dwork, dwork_msize * batchCount );
   /* check allocation */
    if ( dW1_displ == NULL || dW2_displ == NULL || dW3_displ   == NULL || dW4_displ   == NULL || 
         dinvA_array == NULL || dwork_array == NULL || dinvA     == NULL || dwork     == NULL ) {
        magma_tally3_free(dW1_displ);
        magma_tally3_free(dW2_displ);
        magma_tally3_free(dW3_displ);
        magma_tally3_free(dW4_displ);
        magma_tally3_free(dinvA_array);
        magma_tally3_free(dwork_array);
        magma_tally3_free( dinvA );
        magma_tally3_free( dwork );
        info = MAGMA_tally3_ERR_DEVICE_ALLOC;
        magma_tally3_xerbla( __func__, -(info) );
        return info;
    }

    magma_tally3blas_claset_q(Magma_tally3Full, invA_msize, batchCount, MAGMA_tally3_C_ZERO, MAGMA_tally3_C_ZERO, dinvA, invA_msize, queue);
    magma_tally3blas_claset_q(Magma_tally3Full, dwork_msize, batchCount, MAGMA_tally3_C_ZERO, MAGMA_tally3_C_ZERO, dwork, dwork_msize, queue);
    cset_pointer(dwork_array, dwork, n, 0, 0, dwork_msize, batchCount, queue);
    cset_pointer(dinvA_array, dinvA, TRI_NB, 0, 0, invA_msize, batchCount, queue);


    magma_tally3_queue_t cstream;
    magma_tally3blasGetKernelStream(&cstream);


    if (notran) {
        magma_tally3_claswp_rowserial_batched(nrhs, dB_array, lddb, 1, n, dipiv_array, batchCount, queue);
        // solve dwork = L^-1 * NRHS
        magma_tally3blas_ctrsm_outofplace_batched(Magma_tally3Left, Magma_tally3Lower, Magma_tally3NoTrans, Magma_tally3Unit, 1,
                n, nrhs,
                MAGMA_tally3_C_ONE,
                dA_array,       ldda, // dA
                dB_array,      lddb, // dB
                dwork_array,        n, // dX //output
                dinvA_array,  invA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue);

        // solve X = U^-1 * dwork
        magma_tally3blas_ctrsm_outofplace_batched(Magma_tally3Left, Magma_tally3Upper, Magma_tally3NoTrans, Magma_tally3NonUnit, 1,
                n, nrhs,
                MAGMA_tally3_C_ONE,
                dA_array,       ldda, // dA
                dwork_array,        n, // dB 
                dB_array,   lddb, // dX //output
                dinvA_array,  invA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue);
    }
    else{
        /* Solve A**T * X = B  or  A**H * X = B. */
        // solve 
        magma_tally3blas_ctrsm_outofplace_batched(Magma_tally3Left, Magma_tally3Upper, trans, Magma_tally3Unit, 1,
                n, nrhs,
                MAGMA_tally3_C_ONE,
                dA_array,       ldda, // dA
                dB_array,      lddb, // dB
                dwork_array,        n, // dX //output
                dinvA_array,  invA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue);

        // solve 
        magma_tally3blas_ctrsm_outofplace_batched(Magma_tally3Left, Magma_tally3Lower, trans, Magma_tally3NonUnit, 1,
                n, nrhs,
                MAGMA_tally3_C_ONE,
                dA_array,       ldda, // dA
                dwork_array,        n, // dB 
                dB_array,   lddb, // dX //output
                dinvA_array,  invA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue);
        magma_tally3_claswp_rowserial_batched(nrhs, dB_array, lddb, 1, n, dipiv_array, batchCount, queue);
    }




    magma_tally3_queue_sync(cstream);

    magma_tally3_free(dW1_displ);
    magma_tally3_free(dW2_displ);
    magma_tally3_free(dW3_displ);
    magma_tally3_free(dW4_displ);
    magma_tally3_free(dinvA_array);
    magma_tally3_free(dwork_array);
    magma_tally3_free( dinvA );
    magma_tally3_free( dwork );

    
    return info;
}
