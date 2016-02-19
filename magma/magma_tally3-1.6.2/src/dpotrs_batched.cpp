/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Azzam Haidar

       @generated from zpotrs_batched.cpp normal z -> d, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_tally3.h"
#include "batched_kernel_param.h"
/**
    Purpose
    -------
    DPOTRS solves a system of linear equations A*X = B with a symmetric
    positive definite matrix A using the Cholesky factorization
    A = U**H*U or A = L*L**H computed by DPOTRF.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally3_uplo_t
      -     = Magma_tally3Upper:  Upper triangle of A is stored;
      -     = Magma_tally3Lower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in]
    dA      DOUBLE_PRECISION array on the GPU, dimension (LDDA,N)
            The triangular factor U or L from the Cholesky factorization
            A = U**H*U or A = L*L**H, as computed by DPOTRF.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[in,out]
    dB      DOUBLE_PRECISION array on the GPU, dimension (LDDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally3_dposv_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_dpotrs_batched(
                  magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  double **dA_array, magma_tally3_int_t ldda,
                  double **dB_array, magma_tally3_int_t lddb,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue)
{
    double c_one = MAGMA_tally3_D_ONE;
    magma_tally3_int_t info = 0;
    if ( uplo != Magma_tally3Upper && uplo != Magma_tally3Lower )
        info = -1;
    if ( n < 0 )
        info = -2;
    if ( nrhs < 0)
        info = -3;
    if ( ldda < max(1, n) )
        info = -5;
    if ( lddb < max(1, n) )
        info = -7;
    if (info != 0) {
        magma_tally3_xerbla( __func__, -(info) );
        return info;
    }

    /* Quick return if possible */
    if ( (n == 0) || (nrhs == 0) ) {
        return info;
    }


    double **dW1_displ  = NULL;
    double **dW2_displ  = NULL;
    double **dW3_displ  = NULL;
    double **dW4_displ  = NULL;
    double **dinvA_array = NULL;
    double **dwork_array = NULL;



    magma_tally3_malloc((void**)&dW1_displ,  batchCount * sizeof(*dW1_displ));
    magma_tally3_malloc((void**)&dW2_displ,  batchCount * sizeof(*dW2_displ));
    magma_tally3_malloc((void**)&dW3_displ,  batchCount * sizeof(*dW3_displ));
    magma_tally3_malloc((void**)&dW4_displ,  batchCount * sizeof(*dW4_displ));
    magma_tally3_malloc((void**)&dinvA_array, batchCount * sizeof(*dinvA_array));
    magma_tally3_malloc((void**)&dwork_array, batchCount * sizeof(*dwork_array));


    magma_tally3_int_t invA_msize = ((n+TRI_NB-1)/TRI_NB)*TRI_NB*TRI_NB;
    magma_tally3_int_t dwork_msize = n*nrhs;
    double* dinvA      = NULL;
    double* dwork      = NULL;// dinvA and dwork are workspace in dtrsm
    magma_tally3_dmalloc( &dinvA, invA_msize * batchCount);
    magma_tally3_dmalloc( &dwork, dwork_msize * batchCount );
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

    magma_tally3blas_dlaset_q(Magma_tally3Full, invA_msize, batchCount, MAGMA_tally3_D_ZERO, MAGMA_tally3_D_ZERO, dinvA, invA_msize, queue);
    magma_tally3blas_dlaset_q(Magma_tally3Full, dwork_msize, batchCount, MAGMA_tally3_D_ZERO, MAGMA_tally3_D_ZERO, dwork, dwork_msize, queue);
    dset_pointer(dwork_array, dwork, n, 0, 0, dwork_msize, batchCount, queue);
    dset_pointer(dinvA_array, dinvA, TRI_NB, 0, 0, invA_msize, batchCount, queue);

    magma_tally3_queue_t cstream;
    magma_tally3blasGetKernelStream(&cstream);


    if ( uplo == Magma_tally3Upper) {
        // A = U^T U
        // solve U^{T}X = B ==> dworkX = U^-T * B
        magma_tally3blas_dtrsm_outofplace_batched(Magma_tally3Left, Magma_tally3Upper, Magma_tally3ConjTrans, Magma_tally3NonUnit, 1,
                n, nrhs,
                c_one,
                dA_array,       ldda, // dA
                dB_array,      lddb, // dB
                dwork_array,        n, // dX //output
                dinvA_array,  invA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue);

        // solve U X = dwork ==> X = U^-1 * dwork
        magma_tally3blas_dtrsm_outofplace_batched(Magma_tally3Left, Magma_tally3Upper, Magma_tally3NoTrans, Magma_tally3NonUnit, 1,
                n, nrhs,
                c_one,
                dA_array,       ldda, // dA
                dwork_array,        n, // dB 
                dB_array,   lddb, // dX //output
                dinvA_array,  invA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue);
    }
    else{
        // A = L L^T
        // solve LX= B ==> dwork = L^{-1} B
        magma_tally3blas_dtrsm_outofplace_batched(Magma_tally3Left, Magma_tally3Lower, Magma_tally3NoTrans, Magma_tally3NonUnit, 1,
                n, nrhs,
                c_one,
                dA_array,       ldda, // dA
                dB_array,      lddb, // dB
                dwork_array,        n, // dX //output
                dinvA_array,  invA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue);

        // solve L^{T}X= dwork ==> X = L^{-T} dwork
        magma_tally3blas_dtrsm_outofplace_batched(Magma_tally3Left, Magma_tally3Lower, Magma_tally3ConjTrans, Magma_tally3NonUnit, 1,
                n, nrhs,
                c_one,
                dA_array,       ldda, // dA
                dwork_array,        n, // dB 
                dB_array,   lddb, // dX //output
                dinvA_array,  invA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue);
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
