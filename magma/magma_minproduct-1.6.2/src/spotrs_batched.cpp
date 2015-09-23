/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Azzam Haidar

       @generated from zpotrs_batched.cpp normal z -> s, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_minproduct.h"
#include "batched_kernel_param.h"
/**
    Purpose
    -------
    SPOTRS solves a system of linear equations A*X = B with a symmetric
    positive definite matrix A using the Cholesky factorization
    A = U**H*U or A = L*L**H computed by SPOTRF.

    Arguments
    ---------
    @param[in]
    uplo    magma_minproduct_uplo_t
      -     = Magma_minproductUpper:  Upper triangle of A is stored;
      -     = Magma_minproductLower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in]
    dA      REAL array on the GPU, dimension (LDDA,N)
            The triangular factor U or L from the Cholesky factorization
            A = U**H*U or A = L*L**H, as computed by SPOTRF.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[in,out]
    dB      REAL array on the GPU, dimension (LDDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_minproduct_sposv_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_spotrs_batched(
                  magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  float **dA_array, magma_minproduct_int_t ldda,
                  float **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue)
{
    float c_one = MAGMA_minproduct_S_ONE;
    magma_minproduct_int_t info = 0;
    if ( uplo != Magma_minproductUpper && uplo != Magma_minproductLower )
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
        magma_minproduct_xerbla( __func__, -(info) );
        return info;
    }

    /* Quick return if possible */
    if ( (n == 0) || (nrhs == 0) ) {
        return info;
    }


    float **dW1_displ  = NULL;
    float **dW2_displ  = NULL;
    float **dW3_displ  = NULL;
    float **dW4_displ  = NULL;
    float **dinvA_array = NULL;
    float **dwork_array = NULL;



    magma_minproduct_malloc((void**)&dW1_displ,  batchCount * sizeof(*dW1_displ));
    magma_minproduct_malloc((void**)&dW2_displ,  batchCount * sizeof(*dW2_displ));
    magma_minproduct_malloc((void**)&dW3_displ,  batchCount * sizeof(*dW3_displ));
    magma_minproduct_malloc((void**)&dW4_displ,  batchCount * sizeof(*dW4_displ));
    magma_minproduct_malloc((void**)&dinvA_array, batchCount * sizeof(*dinvA_array));
    magma_minproduct_malloc((void**)&dwork_array, batchCount * sizeof(*dwork_array));


    magma_minproduct_int_t invA_msize = ((n+TRI_NB-1)/TRI_NB)*TRI_NB*TRI_NB;
    magma_minproduct_int_t dwork_msize = n*nrhs;
    float* dinvA      = NULL;
    float* dwork      = NULL;// dinvA and dwork are workspace in strsm
    magma_minproduct_smalloc( &dinvA, invA_msize * batchCount);
    magma_minproduct_smalloc( &dwork, dwork_msize * batchCount );
   /* check allocation */
    if ( dW1_displ == NULL || dW2_displ == NULL || dW3_displ   == NULL || dW4_displ   == NULL || 
         dinvA_array == NULL || dwork_array == NULL || dinvA     == NULL || dwork     == NULL ) {
        magma_minproduct_free(dW1_displ);
        magma_minproduct_free(dW2_displ);
        magma_minproduct_free(dW3_displ);
        magma_minproduct_free(dW4_displ);
        magma_minproduct_free(dinvA_array);
        magma_minproduct_free(dwork_array);
        magma_minproduct_free( dinvA );
        magma_minproduct_free( dwork );
        info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
        magma_minproduct_xerbla( __func__, -(info) );
        return info;
    }

    magma_minproductblas_slaset_q(Magma_minproductFull, invA_msize, batchCount, MAGMA_minproduct_S_ZERO, MAGMA_minproduct_S_ZERO, dinvA, invA_msize, queue);
    magma_minproductblas_slaset_q(Magma_minproductFull, dwork_msize, batchCount, MAGMA_minproduct_S_ZERO, MAGMA_minproduct_S_ZERO, dwork, dwork_msize, queue);
    sset_pointer(dwork_array, dwork, n, 0, 0, dwork_msize, batchCount, queue);
    sset_pointer(dinvA_array, dinvA, TRI_NB, 0, 0, invA_msize, batchCount, queue);

    magma_minproduct_queue_t cstream;
    magma_minproductblasGetKernelStream(&cstream);


    if ( uplo == Magma_minproductUpper) {
        // A = U^T U
        // solve U^{T}X = B ==> dworkX = U^-T * B
        magma_minproductblas_strsm_outofplace_batched(Magma_minproductLeft, Magma_minproductUpper, Magma_minproductConjTrans, Magma_minproductNonUnit, 1,
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
        magma_minproductblas_strsm_outofplace_batched(Magma_minproductLeft, Magma_minproductUpper, Magma_minproductNoTrans, Magma_minproductNonUnit, 1,
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
        magma_minproductblas_strsm_outofplace_batched(Magma_minproductLeft, Magma_minproductLower, Magma_minproductNoTrans, Magma_minproductNonUnit, 1,
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
        magma_minproductblas_strsm_outofplace_batched(Magma_minproductLeft, Magma_minproductLower, Magma_minproductConjTrans, Magma_minproductNonUnit, 1,
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




    magma_minproduct_queue_sync(cstream);

    magma_minproduct_free(dW1_displ);
    magma_minproduct_free(dW2_displ);
    magma_minproduct_free(dW3_displ);
    magma_minproduct_free(dW4_displ);
    magma_minproduct_free(dinvA_array);
    magma_minproduct_free(dwork_array);
    magma_minproduct_free( dinvA );
    magma_minproduct_free( dwork );

    return info;
}
