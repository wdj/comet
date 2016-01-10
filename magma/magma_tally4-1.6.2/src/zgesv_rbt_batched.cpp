/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2013
       
       @author Azzam Haidar

       @precisions normal z -> s d c
*/
#include "common_magma_tally4.h"
#include "batched_kernel_param.h"
/**
    Purpose
    -------
    Solves a system of linear equations
      A * X = B,  A**T * X = B,  or  A**H * X = B
    with a general N-by-N matrix A using the LU factorization computed by ZGETRF_GPU.

    Arguments
    ---------
    @param[in]
    trans   magma_tally4_trans_t
            Specifies the form of the system of equations:
      -     = Magma_tally4NoTrans:    A    * X = B  (No transpose)
      -     = Magma_tally4Trans:      A**T * X = B  (Transpose)
      -     = Magma_tally4ConjTrans:  A**H * X = B  (Conjugate transpose)

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in]
    dA      COMPLEX_16 array on the GPU, dimension (LDA,N)
            The factors L and U from the factorization A = P*L*U as computed
            by ZGETRF_GPU.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[in]
    ipiv    INTEGER array, dimension (N)
            The pivot indices from ZGETRF; for 1 <= i <= N, row i of the
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

    @ingroup magma_tally4_zgesv_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_zgesv_rbt_batched(
                  magma_tally4_int_t n, magma_tally4_int_t nrhs,
                  magma_tally4DoubleComplex **dA_array, magma_tally4_int_t ldda,
                  magma_tally4DoubleComplex **dB_array, magma_tally4_int_t lddb,
                  magma_tally4_int_t *info_array,
                  magma_tally4_int_t batchCount, magma_tally4_queue_t queue)
{
    /* Local variables */
    
    magma_tally4_int_t info;
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
        magma_tally4_xerbla( __func__, -(info) );
        return info;
    }


    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return info;
    }

    magma_tally4DoubleComplex *hu, *hv;
    if (MAGMA_tally4_SUCCESS != magma_tally4_zmalloc_cpu( &hu, 2*n )) {
        info = MAGMA_tally4_ERR_HOST_ALLOC;
        return info;
    }

    if (MAGMA_tally4_SUCCESS != magma_tally4_zmalloc_cpu( &hv, 2*n )) {
        info = MAGMA_tally4_ERR_HOST_ALLOC;
        return info;
    }



    info = magma_tally4_zgerbt_batched(Magma_tally4True, n, nrhs, dA_array, n, dB_array, n, hu, hv, &info, batchCount, queue);
    if (info != MAGMA_tally4_SUCCESS)  {
        return info;
    }


    info = magma_tally4_zgetrf_nopiv_batched( n, n, dA_array, ldda, info_array, batchCount, queue);
    if ( (info != MAGMA_tally4_SUCCESS) ){
        return info;
    }

#ifdef CHECK_INFO
    // check correctness of results throught "dinfo_magma_tally4" and correctness of argument throught "info"
    magma_tally4_int_t *cpu_info = (magma_tally4_int_t*) malloc(batchCount*sizeof(magma_tally4_int_t));
    magma_tally4_getvector( batchCount, sizeof(magma_tally4_int_t), dinfo_array, 1, cpu_info, 1);
    for(int i=0; i<batchCount; i++)
    {
        if(cpu_info[i] != 0 ){
            printf("magma_tally4_zgetrf_batched matrix %d returned error %d\n",i, (int)cpu_info[i] );
            info = cpu_info[i];
            free (cpu_info);
            return info;
        }
    }
    free (cpu_info);
#endif

    info = magma_tally4_zgetrs_nopiv_batched( Magma_tally4NoTrans, n, nrhs, dA_array, ldda, dB_array, lddb, info_array, batchCount, queue );


    /* The solution of A.x = b is Vy computed on the GPU */
    magma_tally4DoubleComplex *dv;

    if (MAGMA_tally4_SUCCESS != magma_tally4_zmalloc( &dv, 2*n )) {
        info = MAGMA_tally4_ERR_DEVICE_ALLOC;
        return info;
    }

    magma_tally4_zsetvector(2*n, hv, 1, dv, 1);

    for(int i = 0; i < nrhs; i++)
        magma_tally4blas_zprbt_mv_batched(n, dv, dB_array+(i), batchCount, queue);

 //   magma_tally4_zgetmatrix(n, nrhs, db, nn, B, ldb);


    return info;
}