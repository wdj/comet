/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Azzam Haidar

       @generated from zpotrs_batched.cpp normal z -> c, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_tally2.h"
#include "batched_kernel_param.h"
/**
    Purpose
    -------
    CPOTRS solves a system of linear equations A*X = B with a Hermitian
    positive definite matrix A using the Cholesky factorization
    A = U**H*U or A = L*L**H computed by CPOTRF.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally2_uplo_t
      -     = Magma_tally2Upper:  Upper triangle of A is stored;
      -     = Magma_tally2Lower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in]
    dA      COMPLEX array on the GPU, dimension (LDDA,N)
            The triangular factor U or L from the Cholesky factorization
            A = U**H*U or A = L*L**H, as computed by CPOTRF.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[in,out]
    dB      COMPLEX array on the GPU, dimension (LDDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally2_cposv_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_cpotrs_batched(
                  magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
                  magma_tally2FloatComplex **dA_array, magma_tally2_int_t ldda,
                  magma_tally2FloatComplex **dB_array, magma_tally2_int_t lddb,
                  magma_tally2_int_t batchCount, magma_tally2_queue_t queue)
{
    magma_tally2FloatComplex c_one = MAGMA_tally2_C_ONE;
    magma_tally2_int_t info = 0;
    if ( uplo != Magma_tally2Upper && uplo != Magma_tally2Lower )
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
        magma_tally2_xerbla( __func__, -(info) );
        return info;
    }

    /* Quick return if possible */
    if ( (n == 0) || (nrhs == 0) ) {
        return info;
    }


    magma_tally2FloatComplex **dW1_displ  = NULL;
    magma_tally2FloatComplex **dW2_displ  = NULL;
    magma_tally2FloatComplex **dW3_displ  = NULL;
    magma_tally2FloatComplex **dW4_displ  = NULL;
    magma_tally2FloatComplex **dinvA_array = NULL;
    magma_tally2FloatComplex **dwork_array = NULL;



    magma_tally2_malloc((void**)&dW1_displ,  batchCount * sizeof(*dW1_displ));
    magma_tally2_malloc((void**)&dW2_displ,  batchCount * sizeof(*dW2_displ));
    magma_tally2_malloc((void**)&dW3_displ,  batchCount * sizeof(*dW3_displ));
    magma_tally2_malloc((void**)&dW4_displ,  batchCount * sizeof(*dW4_displ));
    magma_tally2_malloc((void**)&dinvA_array, batchCount * sizeof(*dinvA_array));
    magma_tally2_malloc((void**)&dwork_array, batchCount * sizeof(*dwork_array));


    magma_tally2_int_t invA_msize = ((n+TRI_NB-1)/TRI_NB)*TRI_NB*TRI_NB;
    magma_tally2_int_t dwork_msize = n*nrhs;
    magma_tally2FloatComplex* dinvA      = NULL;
    magma_tally2FloatComplex* dwork      = NULL;// dinvA and dwork are workspace in ctrsm
    magma_tally2_cmalloc( &dinvA, invA_msize * batchCount);
    magma_tally2_cmalloc( &dwork, dwork_msize * batchCount );
   /* check allocation */
    if ( dW1_displ == NULL || dW2_displ == NULL || dW3_displ   == NULL || dW4_displ   == NULL || 
         dinvA_array == NULL || dwork_array == NULL || dinvA     == NULL || dwork     == NULL ) {
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

    magma_tally2blas_claset_q(Magma_tally2Full, invA_msize, batchCount, MAGMA_tally2_C_ZERO, MAGMA_tally2_C_ZERO, dinvA, invA_msize, queue);
    magma_tally2blas_claset_q(Magma_tally2Full, dwork_msize, batchCount, MAGMA_tally2_C_ZERO, MAGMA_tally2_C_ZERO, dwork, dwork_msize, queue);
    cset_pointer(dwork_array, dwork, n, 0, 0, dwork_msize, batchCount, queue);
    cset_pointer(dinvA_array, dinvA, TRI_NB, 0, 0, invA_msize, batchCount, queue);

    magma_tally2_queue_t cstream;
    magma_tally2blasGetKernelStream(&cstream);


    if ( uplo == Magma_tally2Upper) {
        // A = U^T U
        // solve U^{T}X = B ==> dworkX = U^-T * B
        magma_tally2blas_ctrsm_outofplace_batched(Magma_tally2Left, Magma_tally2Upper, Magma_tally2ConjTrans, Magma_tally2NonUnit, 1,
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
        magma_tally2blas_ctrsm_outofplace_batched(Magma_tally2Left, Magma_tally2Upper, Magma_tally2NoTrans, Magma_tally2NonUnit, 1,
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
        magma_tally2blas_ctrsm_outofplace_batched(Magma_tally2Left, Magma_tally2Lower, Magma_tally2NoTrans, Magma_tally2NonUnit, 1,
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
        magma_tally2blas_ctrsm_outofplace_batched(Magma_tally2Left, Magma_tally2Lower, Magma_tally2ConjTrans, Magma_tally2NonUnit, 1,
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




    magma_tally2_queue_sync(cstream);

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