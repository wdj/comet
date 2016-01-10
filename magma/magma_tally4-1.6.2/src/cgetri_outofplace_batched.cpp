/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Azzam Haidar
       @author Tingxing Dong

       @generated from zgetri_outofplace_batched.cpp normal z -> c, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_tally4.h"
#include "batched_kernel_param.h"

/**
    Purpose
    -------
    CGETRI computes the inverse of a matrix using the LU factorization
    computed by CGETRF. This method inverts U and then computes inv(A) by
    solving the system inv(A)*L = inv(U) for inv(A).
    
    Note that it is generally both faster and more accurate to use CGESV,
    or CGETRF and CGETRS, to solve the system AX = B, rather than inverting
    the matrix and multiplying to form X = inv(A)*B. Only in special
    instances should an explicit inverse be computed with this routine.

    Arguments
    ---------
    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    dA      COMPLEX array on the GPU, dimension (LDDA,N)
            On entry, the factors L and U from the factorization
            A = P*L*U as computed by CGETRF_GPU.
            On exit, if INFO = 0, the inverse of the original matrix A.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[in]
    ipiv    INTEGER array, dimension (N)
            The pivot indices from CGETRF; for 1 <= i <= N, row i of the
            matrix was interchanged with row IPIV(i).

    @param[out]
    dwork   (workspace) COMPLEX array on the GPU, dimension (MAX(1,LWORK))
  
    @param[in]
    lwork   INTEGER
            The dimension of the array DWORK.  LWORK >= N*NB, where NB is
            the optimal blocksize returned by magma_tally4_get_cgetri_nb(n).
    \n
            Unlike LAPACK, this version does not currently support a
            workspace query, because the workspace is on the GPU.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
                  singular and its cannot be computed.

    @ingroup magma_tally4_cgesv_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_cgetri_outofplace_batched( magma_tally4_int_t n, 
                  magma_tally4FloatComplex **dA_array, magma_tally4_int_t ldda,
                  magma_tally4_int_t **dipiv_array, 
                  magma_tally4FloatComplex **dinvA_array, magma_tally4_int_t lddia,
                  magma_tally4_int_t *info_array,
                  magma_tally4_int_t batchCount, magma_tally4_queue_t queue)
       
{
    /* Local variables */
  
    magma_tally4_int_t info = 0;
    if (n < 0)
        info = -1;
    else if (ldda < max(1,n))
        info = -3;
    else if (lddia < max(1,n))
        info = -6;

    if (info != 0) {
        magma_tally4_xerbla( __func__, -(info) );
        return info;
    }

    /* Quick return if possible */
    if ( n == 0 )
        return info;




    magma_tally4_int_t ib, j;
    magma_tally4_int_t nb = 256;//256;// BATRF_NB;



    magma_tally4FloatComplex **dA_displ   = NULL;
    magma_tally4FloatComplex **dW0_displ  = NULL;
    magma_tally4FloatComplex **dW1_displ  = NULL;
    magma_tally4FloatComplex **dW2_displ  = NULL;
    magma_tally4FloatComplex **dW3_displ  = NULL;
    magma_tally4FloatComplex **dW4_displ  = NULL;
    magma_tally4FloatComplex **dinvdiagA_array = NULL;
    magma_tally4FloatComplex **dwork_array = NULL;
    magma_tally4FloatComplex **dW5_displ   = NULL;
    magma_tally4_malloc((void**)&dA_displ,   batchCount * sizeof(*dA_displ));
    magma_tally4_malloc((void**)&dW0_displ,  batchCount * sizeof(*dW0_displ));
    magma_tally4_malloc((void**)&dW1_displ,  batchCount * sizeof(*dW1_displ));
    magma_tally4_malloc((void**)&dW2_displ,  batchCount * sizeof(*dW2_displ));
    magma_tally4_malloc((void**)&dW3_displ,  batchCount * sizeof(*dW3_displ));
    magma_tally4_malloc((void**)&dW4_displ,  batchCount * sizeof(*dW4_displ));
    magma_tally4_malloc((void**)&dinvdiagA_array, batchCount * sizeof(*dinvdiagA_array));
    magma_tally4_malloc((void**)&dwork_array, batchCount * sizeof(*dwork_array));
    magma_tally4_malloc((void**)&dW5_displ,  batchCount * sizeof(*dW5_displ));

    magma_tally4FloatComplex* dinvdiagA;
    magma_tally4FloatComplex* dwork;// dinvdiagA and dwork are workspace in ctrsm
    //magma_tally4_int_t invdiagA_msize =  BATRI_NB*((nb/BATRI_NB)+(nb % BATRI_NB != 0))* BATRI_NB ;
    magma_tally4_int_t invdiagA_msize = ((n+TRI_NB-1)/TRI_NB)*TRI_NB*TRI_NB;
    magma_tally4_int_t dwork_msize = n*nb;
    magma_tally4_cmalloc( &dinvdiagA, invdiagA_msize * batchCount);
    magma_tally4_cmalloc( &dwork, dwork_msize * batchCount );
    /* check allocation */
    if ( dA_displ  == NULL || dW1_displ == NULL || dW2_displ       == NULL || dW3_displ   == NULL || 
         dW4_displ == NULL || dW5_displ  == NULL || dinvdiagA_array == NULL || dwork_array == NULL || 
         dinvdiagA == NULL || dwork     == NULL ) {
        magma_tally4_free(dA_displ);
        magma_tally4_free(dW1_displ);
        magma_tally4_free(dW2_displ);
        magma_tally4_free(dW3_displ);
        magma_tally4_free(dW4_displ);
        magma_tally4_free(dW5_displ);
        magma_tally4_free(dinvdiagA_array);
        magma_tally4_free(dwork_array);
        magma_tally4_free(dinvdiagA);
        magma_tally4_free( dwork );
        info = MAGMA_tally4_ERR_DEVICE_ALLOC;
        magma_tally4_xerbla( __func__, -(info) );
        return info;
    }

    magma_tally4blas_claset_q(Magma_tally4Full, invdiagA_msize, batchCount, MAGMA_tally4_C_ZERO, MAGMA_tally4_C_ZERO, dinvdiagA, invdiagA_msize, queue);
    magma_tally4blas_claset_q(Magma_tally4Full, dwork_msize, batchCount, MAGMA_tally4_C_ZERO, MAGMA_tally4_C_ZERO, dwork, dwork_msize, queue);
    cset_pointer(dwork_array, dwork, n, 0, 0, dwork_msize, batchCount, queue);
    cset_pointer(dinvdiagA_array, dinvdiagA, TRI_NB, 0, 0, invdiagA_msize, batchCount, queue);

    magma_tally4_queue_t cstream;
    magma_tally4blasGetKernelStream(&cstream);

    //printf(" I am after malloc getri\n");


    magma_tally4_cdisplace_pointers(dA_displ, dA_array, ldda, 0, 0, batchCount, queue);
    // set dinvdiagA to identity
    magma_tally4blas_claset_batched(Magma_tally4UpperLower, n, n, MAGMA_tally4_C_ZERO, MAGMA_tally4_C_ONE, dinvA_array, lddia, batchCount, queue);

    for(j = 0; j < n; j+=nb) {
        ib = min(nb, n-j);
        // dinvdiagA * Piv' = I * U^-1 * L^-1 = U^-1 * L^-1 * I
        // Azzam : optimization can be done:
        //          2- compute invdiagL invdiagU only one time


        //magma_tally4_queue_sync(NULL);
        //printf(" @ step %d calling solve 1 \n",j);
        // solve dwork = L^-1 * I
        magma_tally4blas_claset_batched(Magma_tally4UpperLower, j, ib, MAGMA_tally4_C_ZERO, MAGMA_tally4_C_ZERO, dwork_array, n, batchCount, queue);
        magma_tally4_cdisplace_pointers(dW5_displ, dwork_array, n, j, 0, batchCount, queue);
        magma_tally4_cdisplace_pointers(dW0_displ, dinvA_array, lddia, j, j, batchCount, queue);
        magma_tally4_cdisplace_pointers(dA_displ, dA_array, ldda, j, j, batchCount, queue);
        
        magma_tally4blas_ctrsm_outofplace_batched(Magma_tally4Left, Magma_tally4Lower, Magma_tally4NoTrans, Magma_tally4Unit, 1,
                n-j, ib,
                MAGMA_tally4_C_ONE,
                dA_displ,       ldda, // dA
                dW0_displ,   lddia, // dB
                dW5_displ,        n, // dX //output
                dinvdiagA_array,  invdiagA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue);
        
        //magma_tally4_queue_sync(NULL);
        //printf(" @ step %d calling solve 2 \n",j);
        // solve dinvdiagA = U^-1 * dwork
        magma_tally4_cdisplace_pointers(dW5_displ, dwork_array, n, 0, 0, batchCount, queue);
        magma_tally4_cdisplace_pointers(dW0_displ, dinvA_array, lddia, 0, j, batchCount, queue);
        magma_tally4_cdisplace_pointers(dA_displ, dA_array, ldda, 0, 0, batchCount, queue);
        magma_tally4blas_ctrsm_outofplace_batched(Magma_tally4Left, Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4NonUnit, 1,
                n, ib,
                MAGMA_tally4_C_ONE,
                dA_displ,       ldda, // dA
                dW5_displ,        n, // dB 
                dW0_displ,   lddia, // dX //output
                dinvdiagA_array,  invdiagA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue);
    }

    // Apply column interchanges
    magma_tally4_claswp_columnserial_batched( n, dinvA_array, lddia, max(1,n-1), 1, dipiv_array, batchCount, queue);

    magma_tally4_queue_sync(cstream);

    magma_tally4_free(dA_displ);
    magma_tally4_free(dW1_displ);
    magma_tally4_free(dW2_displ);
    magma_tally4_free(dW3_displ);
    magma_tally4_free(dW4_displ);
    magma_tally4_free(dW5_displ);
    magma_tally4_free(dinvdiagA_array);
    magma_tally4_free(dwork_array);
    magma_tally4_free(dinvdiagA);
    magma_tally4_free( dwork );

    
    return info;
}