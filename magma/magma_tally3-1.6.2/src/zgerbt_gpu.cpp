/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
       @author Adrien REMY
*/
#include "common_magma_tally3.h"

#define PRECISION_z
#define COMPLEX



static void
init_butterfly(
        magma_tally3_int_t n,
        magma_tally3DoubleComplex* u, magma_tally3DoubleComplex* v)
{
    magma_tally3_int_t idx;
    double u1, v1;
    for (idx=0; idx<n; idx++){
        u1 = exp((((rand() * 1.0)/RAND_MAX)-0.5)/10);
        v1 = exp((((rand() * 1.0)/RAND_MAX)-0.5)/10);
        u[idx] = MAGMA_tally3_Z_MAKE(u1,u1);

        v[idx] = MAGMA_tally3_Z_MAKE(v1,v1);
    }
}


/**
    Purpose
    -------
    Solves a system of linear equations
       A * X = B
    where A is a general n-by-n matrix and X and B are n-by-nrhs matrices.
    Random Butterfly Tranformation is applied on A and B, then
    the LU decomposition with no pivoting is
    used to factor A as
       A = L * U,
    where L is unit lower triangular, and U is
    upper triangular.  The factored form of A is then used to solve the
    system of equations A * X = B.

    Arguments
    ---------
    
    @param[in]
    gen     magma_tally3_bool_t
     -         = Magma_tally3True:     new matrices are generated for U and V
     -         = Magma_tally3False:    matrices U and V given as parameter are used

    
    @param[in]
    n       INTEGER
            The order of the matrix A.  n >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  nrhs >= 0.

    @param[in,out]
    dA      COMPLEX_16 array, dimension (LDA,n).
            On entry, the M-by-n matrix to be factored.
            On exit, the factors L and U from the factorization
            A = L*U; the unit diagonal elements of L are not stored.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,n).

    @param[in,out]
    dB      COMPLEX_16 array, dimension (LDB,nrhs)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDB >= max(1,n).

    @param[in,out]
    U       COMPLEX_16 array, dimension (2,n)
            Random butterfly matrix, if gen = Magma_tally3True U is generated and returned as output;
            else we use U given as input.
            CPU memory

    @param[in,out]
    V       COMPLEX_16 array, dimension (2,n)
            Random butterfly matrix, if gen = Magma_tally3True V is generated and returned as output;
            else we use U given as input.
            CPU memory

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
                  or another error occured, such as memory allocation failed.

    @ingroup magma_tally3_zgesv_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_zgerbt_gpu(
    magma_tally3_bool_t gen, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3DoubleComplex *U, magma_tally3DoubleComplex *V,
    magma_tally3_int_t *info)
{
    /* Function Body */
    *info = 0;
    if ( ! (gen == Magma_tally3True) &&
            ! (gen == Magma_tally3False) ) {
        *info = -1;
    }
    else if (n < 0) {
        *info = -2;
    } else if (nrhs < 0) {
        *info = -3;
    } else if (ldda < max(1,n)) {
        *info = -5;
    } else if (lddb < max(1,n)) {
        *info = -7;
    }
    if (*info != 0) {
        magma_tally3_xerbla( __func__, -(*info) );

        return *info;
    }

    /* Quick return if possible */
    if (nrhs == 0 || n == 0)
        return *info;

    magma_tally3DoubleComplex *du, *dv;

    /* Allocate memory for the buterfly matrices */
    if (MAGMA_tally3_SUCCESS != magma_tally3_zmalloc( &du, 2*n )) {
        *info = MAGMA_tally3_ERR_DEVICE_ALLOC;
        return *info;
    }
    if (MAGMA_tally3_SUCCESS != magma_tally3_zmalloc( &dv, 2*n )) {
        *info = MAGMA_tally3_ERR_DEVICE_ALLOC;
        return *info;
    }

    /* Initialize Butterfly matrix on the CPU*/
    if(gen == Magma_tally3True)
        init_butterfly(2*n, U, V);

    /* Copy the butterfly to the GPU */
    magma_tally3_zsetvector( 2*n, U, 1, du, 1);
    magma_tally3_zsetvector( 2*n, V, 1, dv, 1);

    /* Perform Partial Random Butterfly Transformation on the GPU*/
    magma_tally3blas_zprbt(n, dA, ldda, du, dv);

    /* Compute U^T.b on the GPU*/
    for(int i= 0; i < nrhs; i++)
        magma_tally3blas_zprbt_mtv(n, du, dB+(i*lddb));

    magma_tally3_free( du );
    magma_tally3_free( dv );

    return *info;
}
