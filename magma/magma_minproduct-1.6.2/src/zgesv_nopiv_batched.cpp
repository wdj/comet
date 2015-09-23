/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2013
       
       @author Azzam Haidar

       @precisions normal z -> s d c
*/
#include "common_magma_minproduct.h"
#include "batched_kernel_param.h"
/**
    Purpose
    -------
    Solves a system of linear equations
       A * X = B
    where A is a general N-by-N matrix and X and B are N-by-NRHS matrices.
    The LU decomposition with partial pivoting and row interchanges is
    used to factor A as
       A = P * L * U,
    where P is a permutation matrix, L is unit lower triangular, and U is
    upper triangular.  The factored form of A is then used to solve the
    system of equations A * X = B.

    Arguments
    ---------
    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in,out]
    dA      COMPLEX_16 array on the GPU, dimension (LDDA,N).
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[out]
    ipiv    INTEGER array, dimension (min(M,N))
            The pivot indices; for 1 <= i <= min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    @param[in,out]
    dB      COMPLEX_16 array on the GPU, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_minproduct_zgesv_driver
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_zgesv_nopiv_batched(
                  magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproductDoubleComplex **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue)
{
    /* Local variables */
    
    magma_minproduct_int_t info;
    info = 0;
    if (n < 0) {
        info = -1;
    } else if (nrhs < 0) {
        info = -2;
    } else if (ldda < max(1,n)) {
        info = -4;
    } else if (lddb < max(1,n)) {
        info = -6;
    }
    if (info != 0) {
        magma_minproduct_xerbla( __func__, -(info) );
        return info;
    }


    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return info;
    }
    info = magma_minproduct_zgetrf_nopiv_batched( n, n, dA_array, ldda, info_array, batchCount, queue);
    if ( (info != MAGMA_minproduct_SUCCESS) ){
        return info;
    }

#ifdef CHECK_INFO
    // check correctness of results throught "dinfo_magma_minproduct" and correctness of argument throught "info"
    magma_minproduct_int_t *cpu_info = (magma_minproduct_int_t*) malloc(batchCount*sizeof(magma_minproduct_int_t));
    magma_minproduct_getvector( batchCount, sizeof(magma_minproduct_int_t), dinfo_array, 1, cpu_info, 1);
    for(int i=0; i<batchCount; i++)
    {
        if(cpu_info[i] != 0 ){
            printf("magma_minproduct_zgetrf_batched matrix %d returned error %d\n",i, (int)cpu_info[i] );
            info = cpu_info[i];
            free (cpu_info);
            return info;
        }
    }
    free (cpu_info);
#endif
           
    info = magma_minproduct_zgetrs_nopiv_batched( Magma_minproductNoTrans, n, nrhs, dA_array, ldda, dB_array, lddb, info_array, batchCount, queue );
            
    return info;
}
