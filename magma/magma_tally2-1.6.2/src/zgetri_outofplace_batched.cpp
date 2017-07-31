/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Azzam Haidar
       @author Tingxing Dong

       @precisions normal z -> s d c
*/
#include "common_magma_tally2.h"
#include "batched_kernel_param.h"

/**
    Purpose
    -------
    ZGETRI computes the inverse of a matrix using the LU factorization
    computed by ZGETRF. This method inverts U and then computes inv(A) by
    solving the system inv(A)*L = inv(U) for inv(A).
    
    Note that it is generally both faster and more accurate to use ZGESV,
    or ZGETRF and ZGETRS, to solve the system AX = B, rather than inverting
    the matrix and multiplying to form X = inv(A)*B. Only in special
    instances should an explicit inverse be computed with this routine.

    Arguments
    ---------
    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    dA      COMPLEX_16 array on the GPU, dimension (LDDA,N)
            On entry, the factors L and U from the factorization
            A = P*L*U as computed by ZGETRF_GPU.
            On exit, if INFO = 0, the inverse of the original matrix A.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[in]
    ipiv    INTEGER array, dimension (N)
            The pivot indices from ZGETRF; for 1 <= i <= N, row i of the
            matrix was interchanged with row IPIV(i).

    @param[out]
    dwork   (workspace) COMPLEX_16 array on the GPU, dimension (MAX(1,LWORK))
  
    @param[in]
    lwork   INTEGER
            The dimension of the array DWORK.  LWORK >= N*NB, where NB is
            the optimal blocksize returned by magma_tally2_get_zgetri_nb(n).
    \n
            Unlike LAPACK, this version does not currently support a
            workspace query, because the workspace is on the GPU.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
                  singular and its cannot be computed.

    @ingroup magma_tally2_zgesv_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_zgetri_outofplace_batched( magma_tally2_int_t n, 
                  magma_tally2DoubleComplex **dA_array, magma_tally2_int_t ldda,
                  magma_tally2_int_t **dipiv_array, 
                  magma_tally2DoubleComplex **dinvA_array, magma_tally2_int_t lddia,
                  magma_tally2_int_t *info_array,
                  magma_tally2_int_t batchCount, magma_tally2_queue_t queue)
       
{
    /* Local variables */
  
    magma_tally2_int_t info = 0;
    if (n < 0)
        info = -1;
    else if (ldda < max(1,n))
        info = -3;
    else if (lddia < max(1,n))
        info = -6;

    if (info != 0) {
        magma_tally2_xerbla( __func__, -(info) );
        return info;
    }

    /* Quick return if possible */
    if ( n == 0 )
        return info;




    magma_tally2_int_t ib, j;
    magma_tally2_int_t nb = 256;//256;// BATRF_NB;



    magma_tally2DoubleComplex **dA_displ   = NULL;
    magma_tally2DoubleComplex **dW0_displ  = NULL;
    magma_tally2DoubleComplex **dW1_displ  = NULL;
    magma_tally2DoubleComplex **dW2_displ  = NULL;
    magma_tally2DoubleComplex **dW3_displ  = NULL;
    magma_tally2DoubleComplex **dW4_displ  = NULL;
    magma_tally2DoubleComplex **dinvdiagA_array = NULL;
    magma_tally2DoubleComplex **dwork_array = NULL;
    magma_tally2DoubleComplex **dW5_displ   = NULL;
    magma_tally2_malloc((void**)&dA_displ,   batchCount * sizeof(*dA_displ));
    magma_tally2_malloc((void**)&dW0_displ,  batchCount * sizeof(*dW0_displ));
    magma_tally2_malloc((void**)&dW1_displ,  batchCount * sizeof(*dW1_displ));
    magma_tally2_malloc((void**)&dW2_displ,  batchCount * sizeof(*dW2_displ));
    magma_tally2_malloc((void**)&dW3_displ,  batchCount * sizeof(*dW3_displ));
    magma_tally2_malloc((void**)&dW4_displ,  batchCount * sizeof(*dW4_displ));
    magma_tally2_malloc((void**)&dinvdiagA_array, batchCount * sizeof(*dinvdiagA_array));
    magma_tally2_malloc((void**)&dwork_array, batchCount * sizeof(*dwork_array));
    magma_tally2_malloc((void**)&dW5_displ,  batchCount * sizeof(*dW5_displ));

    magma_tally2DoubleComplex* dinvdiagA;
    magma_tally2DoubleComplex* dwork;// dinvdiagA and dwork are workspace in ztrsm
    //magma_tally2_int_t invdiagA_msize =  BATRI_NB*((nb/BATRI_NB)+(nb % BATRI_NB != 0))* BATRI_NB ;
    magma_tally2_int_t invdiagA_msize = ((n+TRI_NB-1)/TRI_NB)*TRI_NB*TRI_NB;
    magma_tally2_int_t dwork_msize = n*nb;
    magma_tally2_zmalloc( &dinvdiagA, invdiagA_msize * batchCount);
    magma_tally2_zmalloc( &dwork, dwork_msize * batchCount );
    /* check allocation */
    if ( dA_displ  == NULL || dW1_displ == NULL || dW2_displ       == NULL || dW3_displ   == NULL || 
         dW4_displ == NULL || dW5_displ  == NULL || dinvdiagA_array == NULL || dwork_array == NULL || 
         dinvdiagA == NULL || dwork     == NULL ) {
        magma_tally2_free(dA_displ);
        magma_tally2_free(dW1_displ);
        magma_tally2_free(dW2_displ);
        magma_tally2_free(dW3_displ);
        magma_tally2_free(dW4_displ);
        magma_tally2_free(dW5_displ);
        magma_tally2_free(dinvdiagA_array);
        magma_tally2_free(dwork_array);
        magma_tally2_free(dinvdiagA);
        magma_tally2_free( dwork );
        info = MAGMA_tally2_ERR_DEVICE_ALLOC;
        magma_tally2_xerbla( __func__, -(info) );
        return info;
    }

    magma_tally2blas_zlaset_q(Magma_tally2Full, invdiagA_msize, batchCount, MAGMA_tally2_Z_ZERO, MAGMA_tally2_Z_ZERO, dinvdiagA, invdiagA_msize, queue);
    magma_tally2blas_zlaset_q(Magma_tally2Full, dwork_msize, batchCount, MAGMA_tally2_Z_ZERO, MAGMA_tally2_Z_ZERO, dwork, dwork_msize, queue);
    zset_pointer(dwork_array, dwork, n, 0, 0, dwork_msize, batchCount, queue);
    zset_pointer(dinvdiagA_array, dinvdiagA, TRI_NB, 0, 0, invdiagA_msize, batchCount, queue);

    magma_tally2_queue_t cstream;
    magma_tally2blasGetKernelStream(&cstream);

    //printf(" I am after malloc getri\n");


    magma_tally2_zdisplace_pointers(dA_displ, dA_array, ldda, 0, 0, batchCount, queue);
    // set dinvdiagA to identity
    magma_tally2blas_zlaset_batched(Magma_tally2UpperLower, n, n, MAGMA_tally2_Z_ZERO, MAGMA_tally2_Z_ONE, dinvA_array, lddia, batchCount, queue);

    for(j = 0; j < n; j+=nb) {
        ib = min(nb, n-j);
        // dinvdiagA * Piv' = I * U^-1 * L^-1 = U^-1 * L^-1 * I
        // Azzam : optimization can be done:
        //          2- compute invdiagL invdiagU only one time


        //magma_tally2_queue_sync(NULL);
        //printf(" @ step %d calling solve 1 \n",j);
        // solve dwork = L^-1 * I
        magma_tally2blas_zlaset_batched(Magma_tally2UpperLower, j, ib, MAGMA_tally2_Z_ZERO, MAGMA_tally2_Z_ZERO, dwork_array, n, batchCount, queue);
        magma_tally2_zdisplace_pointers(dW5_displ, dwork_array, n, j, 0, batchCount, queue);
        magma_tally2_zdisplace_pointers(dW0_displ, dinvA_array, lddia, j, j, batchCount, queue);
        magma_tally2_zdisplace_pointers(dA_displ, dA_array, ldda, j, j, batchCount, queue);
        
        magma_tally2blas_ztrsm_outofplace_batched(Magma_tally2Left, Magma_tally2Lower, Magma_tally2NoTrans, Magma_tally2Unit, 1,
                n-j, ib,
                MAGMA_tally2_Z_ONE,
                dA_displ,       ldda, // dA
                dW0_displ,   lddia, // dB
                dW5_displ,        n, // dX //output
                dinvdiagA_array,  invdiagA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue);
        
        //magma_tally2_queue_sync(NULL);
        //printf(" @ step %d calling solve 2 \n",j);
        // solve dinvdiagA = U^-1 * dwork
        magma_tally2_zdisplace_pointers(dW5_displ, dwork_array, n, 0, 0, batchCount, queue);
        magma_tally2_zdisplace_pointers(dW0_displ, dinvA_array, lddia, 0, j, batchCount, queue);
        magma_tally2_zdisplace_pointers(dA_displ, dA_array, ldda, 0, 0, batchCount, queue);
        magma_tally2blas_ztrsm_outofplace_batched(Magma_tally2Left, Magma_tally2Upper, Magma_tally2NoTrans, Magma_tally2NonUnit, 1,
                n, ib,
                MAGMA_tally2_Z_ONE,
                dA_displ,       ldda, // dA
                dW5_displ,        n, // dB 
                dW0_displ,   lddia, // dX //output
                dinvdiagA_array,  invdiagA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue);
    }

    // Apply column interchanges
    magma_tally2_zlaswp_columnserial_batched( n, dinvA_array, lddia, max(1,n-1), 1, dipiv_array, batchCount, queue);

    magma_tally2_queue_sync(cstream);

    magma_tally2_free(dA_displ);
    magma_tally2_free(dW1_displ);
    magma_tally2_free(dW2_displ);
    magma_tally2_free(dW3_displ);
    magma_tally2_free(dW4_displ);
    magma_tally2_free(dW5_displ);
    magma_tally2_free(dinvdiagA_array);
    magma_tally2_free(dwork_array);
    magma_tally2_free(dinvdiagA);
    magma_tally2_free( dwork );

    
    return info;
}
