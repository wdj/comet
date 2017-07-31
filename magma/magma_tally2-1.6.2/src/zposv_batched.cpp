/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Azzam Haidar

       @precisions normal z -> s d c
*/
#include "common_magma_tally2.h"
/**
    Purpose
    -------
    ZPOSV computes the solution to a complex system of linear equations
       A * X = B,
    where A is an N-by-N Hermitian positive definite matrix and X and B
    are N-by-NRHS matrices.
    The Cholesky decomposition is used to factor A as
       A = U**H * U,  if UPLO = Magma_tally2Upper, or
       A = L * L**H,  if UPLO = Magma_tally2Lower,
    where U is an upper triangular matrix and  L is a lower triangular
    matrix.  The factored form of A is then used to solve the system of
    equations A * X = B.

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

    @param[in,out]
    dA      COMPLEX_16 array on the GPU, dimension (LDDA,N)
            On entry, the Hermitian matrix dA.  If UPLO = Magma_tally2Upper, the leading
            N-by-N upper triangular part of dA contains the upper
            triangular part of the matrix dA, and the strictly lower
            triangular part of dA is not referenced.  If UPLO = Magma_tally2Lower, the
            leading N-by-N lower triangular part of dA contains the lower
            triangular part of the matrix dA, and the strictly upper
            triangular part of dA is not referenced.
    \n
            On exit, if INFO = 0, the factor U or L from the Cholesky
            factorization dA = U**H*U or dA = L*L**H.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

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

    @ingroup magma_tally2_zposv_driver
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_zposv_batched(
                  magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
                  magma_tally2DoubleComplex **dA_array, magma_tally2_int_t ldda,
                  magma_tally2DoubleComplex **dB_array, magma_tally2_int_t lddb,
                  magma_tally2_int_t *dinfo_array,
                  magma_tally2_int_t batchCount, magma_tally2_queue_t queue)
{
    /* Local variables */
    
    magma_tally2_int_t info = 0;

    if ( uplo != Magma_tally2Upper && uplo != Magma_tally2Lower )
        info = -1;
    if ( n < 0 )
        info = -2;
    if ( nrhs < 0 )
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

    info = magma_tally2_zpotrf_batched( uplo, n, dA_array, ldda, dinfo_array, batchCount, queue);
    if ( (info != MAGMA_tally2_SUCCESS) ){
        return info;
    }


#ifdef CHECK_INFO
    // check correctness of results throught "dinfo_magma_tally2" and correctness of argument throught "info"
    magma_tally2_int_t *cpu_info = (magma_tally2_int_t*) malloc(batchCount*sizeof(magma_tally2_int_t));
    magma_tally2_getvector( batchCount, sizeof(magma_tally2_int_t), dinfo_array, 1, cpu_info, 1);
    for(int i=0; i<batchCount; i++)
    {
        if(cpu_info[i] != 0 ){
            printf("magma_tally2_zpotrf_batched matrix %d returned error %d\n",i, (int)cpu_info[i] );
            info = cpu_info[i];
            free (cpu_info);
            return info;
        }
    }
    free (cpu_info);
#endif

    info = magma_tally2_zpotrs_batched( uplo, n, nrhs, dA_array, ldda, dB_array, lddb,  batchCount, queue );
    return info;
}
