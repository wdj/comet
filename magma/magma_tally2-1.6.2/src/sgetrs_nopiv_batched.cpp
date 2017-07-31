/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2013
       
       @author Azzam Haidar

       @generated from zgetrs_nopiv_batched.cpp normal z -> s, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_tally2.h"
#include "batched_kernel_param.h"
/**
    Purpose
    -------
    Solves a system of linear equations
      A * X = B,  A**T * X = B,  or  A**H * X = B
    with a general N-by-N matrix A using the LU factorization computed by SGETRF_GPU.

    Arguments
    ---------
    @param[in]
    trans   magma_tally2_trans_t
            Specifies the form of the system of equations:
      -     = Magma_tally2NoTrans:    A    * X = B  (No transpose)
      -     = Magma_tally2Trans:      A**T * X = B  (Transpose)
      -     = Magma_tally2ConjTrans:  A**H * X = B  (Conjugate transpose)

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in]
    dA      REAL array on the GPU, dimension (LDA,N)
            The factors L and U from the factorization A = P*L*U as computed
            by SGETRF_GPU.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[in]
    ipiv    INTEGER array, dimension (N)
            The pivot indices from SGETRF; for 1 <= i <= N, row i of the
            matrix was interchanged with row IPIV(i).

    @param[in,out]
    dB      REAL array on the GPU, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally2_sgesv_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_sgetrs_nopiv_batched(
                  magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
                  float **dA_array, magma_tally2_int_t ldda,
                  float **dB_array, magma_tally2_int_t lddb,
                  magma_tally2_int_t *info_array,
                  magma_tally2_int_t batchCount,  magma_tally2_queue_t queue)
{
    /* Local variables */

    magma_tally2_int_t notran = (trans == Magma_tally2NoTrans);

    magma_tally2_int_t info = 0;
    if ( (! notran) &&
         (trans != Magma_tally2Trans) &&
         (trans != Magma_tally2ConjTrans) ) {
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
        magma_tally2_xerbla( __func__, -(info) );
        return info;
    }


    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return info;
    }

    float **dA_displ   = NULL;
    float **dB_displ  = NULL;
    float **dW1_displ  = NULL;
    float **dW2_displ  = NULL;
    float **dW3_displ  = NULL;
    float **dW4_displ  = NULL;
    float **dinvA_array = NULL;
    float **dwork_array = NULL;

    magma_tally2_malloc((void**)&dA_displ,   batchCount * sizeof(*dA_displ));
    magma_tally2_malloc((void**)&dB_displ,  batchCount * sizeof(*dB_displ));
    magma_tally2_malloc((void**)&dW1_displ,  batchCount * sizeof(*dW1_displ));
    magma_tally2_malloc((void**)&dW2_displ,  batchCount * sizeof(*dW2_displ));
    magma_tally2_malloc((void**)&dW3_displ,  batchCount * sizeof(*dW3_displ));
    magma_tally2_malloc((void**)&dW4_displ,  batchCount * sizeof(*dW4_displ));
    magma_tally2_malloc((void**)&dinvA_array, batchCount * sizeof(*dinvA_array));
    magma_tally2_malloc((void**)&dwork_array, batchCount * sizeof(*dwork_array));




    magma_tally2_int_t invA_msize = ((n+TRI_NB-1)/TRI_NB)*TRI_NB*TRI_NB;
    magma_tally2_int_t dwork_msize = n*nrhs;
    float* dinvA      = NULL;
    float* dwork      = NULL;// dinvA and dwork are workspace in strsm
    magma_tally2_smalloc( &dinvA, invA_msize * batchCount);
    magma_tally2_smalloc( &dwork, dwork_msize * batchCount );
   /* check allocation */
    if ( dW1_displ == NULL || dW2_displ == NULL || dW3_displ   == NULL || dW4_displ   == NULL || 
         dinvA_array == NULL || dwork_array == NULL || dinvA     == NULL || dwork     == NULL ||
         dA_displ == NULL || dB_displ == NULL ) {
        magma_tally2_free(dA_displ);
        magma_tally2_free(dB_displ);
        magma_tally2_free(dW1_displ);
        magma_tally2_free(dW2_displ);
        magma_tally2_free(dW3_displ);
        magma_tally2_free(dW4_displ);
        magma_tally2_free(dinvA_array);
        magma_tally2_free(dwork_array);
        magma_tally2_free( dinvA );
        magma_tally2_free( dwork );
        info = MAGMA_tally2_ERR_DEVICE_ALLOC;
        magma_tally2_xerbla( __func__, -(info) );
        return info;
    }

    magma_tally2blas_slaset_q(Magma_tally2Full, invA_msize, batchCount, MAGMA_tally2_S_ZERO, MAGMA_tally2_S_ZERO, dinvA, invA_msize, queue);
    magma_tally2blas_slaset_q(Magma_tally2Full, dwork_msize, batchCount, MAGMA_tally2_S_ZERO, MAGMA_tally2_S_ZERO, dwork, dwork_msize, queue);
    sset_pointer(dwork_array, dwork, n, 0, 0, dwork_msize, batchCount, queue);
    sset_pointer(dinvA_array, dinvA, TRI_NB, 0, 0, invA_msize, batchCount, queue);

    magma_tally2_sdisplace_pointers(dA_displ, dA_array, ldda, 0, 0, batchCount, queue);
    magma_tally2_sdisplace_pointers(dB_displ, dB_array, lddb, 0, 0, batchCount, queue);

    magma_tally2_queue_t cstream;
    magma_tally2blasGetKernelStream(&cstream);

    //printf(" I am after malloc getri\n");

    if (notran) {
        // solve dwork = L^-1 * NRHS
        magma_tally2blas_strsm_outofplace_batched(Magma_tally2Left, Magma_tally2Lower, Magma_tally2NoTrans, Magma_tally2Unit, 1,
                n, nrhs,
                MAGMA_tally2_S_ONE,
                dA_displ,       ldda, // dA
                dB_displ,      lddb, // dB
                dwork_array,        n, // dX //output
                dinvA_array,  invA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue);

        // solve X = U^-1 * dwork
        magma_tally2blas_strsm_outofplace_batched(Magma_tally2Left, Magma_tally2Upper, Magma_tally2NoTrans, Magma_tally2NonUnit, 1,
                n, nrhs,
                MAGMA_tally2_S_ONE,
                dA_displ,       ldda, // dA
                dwork_array,        n, // dB 
                dB_displ,   lddb, // dX //output
                dinvA_array,  invA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue);
    }
    else{
        /* Solve A**T * X = B  or  A**H * X = B. */
        // solve 
        magma_tally2blas_strsm_outofplace_batched(Magma_tally2Left, Magma_tally2Upper, trans, Magma_tally2Unit, 1,
                n, nrhs,
                MAGMA_tally2_S_ONE,
                dA_displ,       ldda, // dA
                dB_displ,      lddb, // dB
                dwork_array,        n, // dX //output
                dinvA_array,  invA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue);

        // solve 
        magma_tally2blas_strsm_outofplace_batched(Magma_tally2Left, Magma_tally2Lower, trans, Magma_tally2NonUnit, 1,
                n, nrhs,
                MAGMA_tally2_S_ONE,
                dA_displ,       ldda, // dA
                dwork_array,        n, // dB 
                dB_displ,   lddb, // dX //output
                dinvA_array,  invA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue);

    }




    magma_tally2_queue_sync(cstream);

    magma_tally2_free(dA_displ);
    magma_tally2_free(dB_displ);
    magma_tally2_free(dW1_displ);
    magma_tally2_free(dW2_displ);
    magma_tally2_free(dW3_displ);
    magma_tally2_free(dW4_displ);
    magma_tally2_free(dinvA_array);
    magma_tally2_free(dwork_array);
    magma_tally2_free( dinvA );
    magma_tally2_free( dwork );
    
    return info;
}
