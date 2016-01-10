/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgerbt_gpu.cpp normal z -> d, Fri Jan 30 19:00:14 2015
       @author Adrien REMY
*/
#include "common_magma_tally4.h"

#define PRECISION_d
#define REAL



static void
init_butterfly(
        magma_tally4_int_t n,
        double* u, double* v)
{
    magma_tally4_int_t idx;
    double u1, v1;
    for (idx=0; idx<n; idx++){
        u1 = exp((((rand() * 1.0)/RAND_MAX)-0.5)/10);
        v1 = exp((((rand() * 1.0)/RAND_MAX)-0.5)/10);
        u[idx] = MAGMA_tally4_D_MAKE(u1,u1);

        v[idx] = MAGMA_tally4_D_MAKE(v1,v1);
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
    gen     magma_tally4_bool_t
     -         = Magma_tally4True:     new matrices are generated for U and V
     -         = Magma_tally4False:    matrices U and V given as parameter are used

    
    @param[in]
    n       INTEGER
            The order of the matrix A.  n >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  nrhs >= 0.

    @param[in,out]
    dA      DOUBLE_PRECISION array, dimension (LDA,n).
            On entry, the M-by-n matrix to be factored.
            On exit, the factors L and U from the factorization
            A = L*U; the unit diagonal elements of L are not stored.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,n).

    @param[in,out]
    dB      DOUBLE_PRECISION array, dimension (LDB,nrhs)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDB >= max(1,n).

    @param[in,out]
    U       DOUBLE_PRECISION array, dimension (2,n)
            Random butterfly matrix, if gen = Magma_tally4True U is generated and returned as output;
            else we use U given as input.
            CPU memory

    @param[in,out]
    V       DOUBLE_PRECISION array, dimension (2,n)
            Random butterfly matrix, if gen = Magma_tally4True V is generated and returned as output;
            else we use U given as input.
            CPU memory

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
                  or another error occured, such as memory allocation failed.

    @ingroup magma_tally4_dgesv_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_dgerbt_gpu(
    magma_tally4_bool_t gen, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    double *U, double *V,
    magma_tally4_int_t *info)
{
    /* Function Body */
    *info = 0;
    if ( ! (gen == Magma_tally4True) &&
            ! (gen == Magma_tally4False) ) {
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
        magma_tally4_xerbla( __func__, -(*info) );

        return *info;
    }

    /* Quick return if possible */
    if (nrhs == 0 || n == 0)
        return *info;

    double *du, *dv;

    /* Allocate memory for the buterfly matrices */
    if (MAGMA_tally4_SUCCESS != magma_tally4_dmalloc( &du, 2*n )) {
        *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
        return *info;
    }
    if (MAGMA_tally4_SUCCESS != magma_tally4_dmalloc( &dv, 2*n )) {
        *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
        return *info;
    }

    /* Initialize Butterfly matrix on the CPU*/
    if(gen == Magma_tally4True)
        init_butterfly(2*n, U, V);

    /* Copy the butterfly to the GPU */
    magma_tally4_dsetvector( 2*n, U, 1, du, 1);
    magma_tally4_dsetvector( 2*n, V, 1, dv, 1);

    /* Perform Partial Random Butterfly Transformation on the GPU*/
    magma_tally4blas_dprbt(n, dA, ldda, du, dv);

    /* Compute U^T.b on the GPU*/
    for(int i= 0; i < nrhs; i++)
        magma_tally4blas_dprbt_mtv(n, du, dB+(i*lddb));

    magma_tally4_free( du );
    magma_tally4_free( dv );

    return *info;
}