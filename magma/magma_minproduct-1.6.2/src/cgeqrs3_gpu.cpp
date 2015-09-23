/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgeqrs3_gpu.cpp normal z -> c, Fri Jan 30 19:00:15 2015

*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    Solves the least squares problem
           min || A*X - C ||
    using the QR factorization A = Q*R computed by CGEQRF3_GPU.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A. M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A. M >= N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of columns of the matrix C. NRHS >= 0.

    @param[in]
    dA      COMPLEX array on the GPU, dimension (LDDA,N)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,n, as returned by
            CGEQRF3_GPU in the first n columns of its array argument A.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A, LDDA >= M.

    @param[in]
    tau     COMPLEX array, dimension (N)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by MAGMA_minproduct_CGEQRF_GPU.

    @param[in,out]
    dB      COMPLEX array on the GPU, dimension (LDDB,NRHS)
            On entry, the M-by-NRHS matrix C.
            On exit, the N-by-NRHS solution matrix X.

    @param[in]
    dT      COMPLEX array that is the output (the 6th argument)
            of magma_minproduct_cgeqrf_gpu of size
            2*MIN(M, N)*NB + ((N+31)/32*32 )* MAX(NB, NRHS).
            The array starts with a block of size MIN(M,N)*NB that stores
            the triangular T matrices used in the QR factorization,
            followed by MIN(M,N)*NB block storing the diagonal block
            matrices for the R matrix, followed by work space of size
            ((N+31)/32*32 )* MAX(NB, NRHS).

    @param[in]
    lddb    INTEGER
            The leading dimension of the array dB. LDDB >= M.

    @param[out]
    hwork   (workspace) COMPLEX array, dimension (LWORK)
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK,
            LWORK >= (M - N + NB)*(NRHS + NB) + NRHS*NB,
            where NB is the blocksize given by magma_minproduct_get_cgeqrf_nb( M ).
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the HWORK array, returns
            this value as the first entry of the WORK array.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_minproduct_cgels_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_cgeqrs3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex_ptr dA,    magma_minproduct_int_t ldda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex_ptr dT,
    magma_minproductFloatComplex_ptr dB,    magma_minproduct_int_t lddb,
    magma_minproductFloatComplex *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info)
{
    #define dA(a_1,a_2) (dA + (a_2)*(ldda) + (a_1))
    #define dT(a_1)     (dT + (lddwork+(a_1))*nb)

    magma_minproductFloatComplex c_one     = MAGMA_minproduct_C_ONE;
    magma_minproduct_int_t k, lddwork;

    magma_minproduct_int_t nb     = magma_minproduct_get_cgeqrf_nb(m);
    magma_minproduct_int_t lwkopt = (m - n + nb)*(nrhs + nb) + nrhs*nb;
    int lquery = (lwork == -1);

    hwork[0] = MAGMA_minproduct_C_MAKE( (float)lwkopt, 0. );

    *info = 0;
    if (m < 0)
        *info = -1;
    else if (n < 0 || m < n)
        *info = -2;
    else if (nrhs < 0)
        *info = -3;
    else if (ldda < max(1,m))
        *info = -5;
    else if (lddb < max(1,m))
        *info = -8;
    else if (lwork < lwkopt && ! lquery)
        *info = -10;

    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery)
        return *info;

    k = min(m,n);
    if (k == 0) {
        hwork[0] = c_one;
        return *info;
    }
    lddwork = k;

    /* B := Q' * B */
    magma_minproduct_cunmqr_gpu( Magma_minproductLeft, Magma_minproduct_ConjTrans,
                      m, nrhs, n,
                      dA(0,0), ldda, tau,
                      dB, lddb, hwork, lwork, dT, nb, info );
    if ( *info != 0 ) {
        return *info;
    }

    /* Solve R*X = B(1:n,:)
       1. Move the (k-1)/nb block diagonal submatrices from dT to R
       2. Solve
       3. Restore the data format moving data from R back to dT
    */
    magma_minproductblas_cswapdblk(k-1, nb, dA(0,0), ldda, 1, dT(0), nb, 0);
    if ( nrhs == 1 ) {
        magma_minproduct_ctrsv(Magma_minproductUpper, Magma_minproductNoTrans, Magma_minproductNonUnit,
                    n, dA(0,0), ldda, dB, 1);
    } else {
        magma_minproduct_ctrsm(Magma_minproductLeft, Magma_minproductUpper, Magma_minproductNoTrans, Magma_minproductNonUnit,
                    n, nrhs, c_one, dA(0,0), ldda, dB, lddb);
    }
    magma_minproductblas_cswapdblk(k-1, nb, dT(0), nb, 0, dA(0,0), ldda, 1);

    return *info;
}

#undef dA
#undef dT
