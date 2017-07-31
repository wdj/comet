/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
       
       @author Stan Tomov
       @author Mathieu Faverge
       @author Ichitaro Yamazaki
       @author Mark Gates
*/
#include "common_magma_tally2.h"

// MAX_PIVOTS is maximum number of pivots to apply in each kernel launch
// NTHREADS is number of threads in a block
// 64 and 256 are better on Kepler; 
//#define MAX_PIVOTS 64
//#define NTHREADS   256
#define MAX_PIVOTS 32
#define NTHREADS   64

typedef struct {
    int npivots;
    int ipiv[MAX_PIVOTS];
} zlaswp_params_t;


// Matrix A is stored row-wise in dAT.
// Divide matrix A into block-columns of NTHREADS columns each.
// Each GPU block processes one block-column of A.
// Each thread goes down a column of A,
// swapping rows according to pivots stored in params.
__global__ void zlaswp_kernel(
    int n, magma_tally2DoubleComplex *dAT, int ldda, zlaswp_params_t params )
{
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if ( tid < n ) {
        dAT += tid;
        magma_tally2DoubleComplex *A1  = dAT;
        
        for( int i1 = 0; i1 < params.npivots; ++i1 ) {
            int i2 = params.ipiv[i1];
            magma_tally2DoubleComplex *A2 = dAT + i2*ldda;
            magma_tally2DoubleComplex temp = *A1;
            *A1 = *A2;
            *A2 = temp;
            A1 += ldda;  // A1 = dA + i1*ldx
        }
    }
}


/**
    Purpose:
    =============
    ZLASWP performs a series of row interchanges on the matrix A.
    One row interchange is initiated for each of rows K1 through K2 of A.
    
    ** Unlike LAPACK, here A is stored row-wise (hence dAT). **
    Otherwise, this is identical to LAPACK's interface.
    
    Arguments:
    ==========
    \param[in]
    n        INTEGER
             The number of columns of the matrix A.
    
    \param[in,out]
    dAT      COMPLEX*16 array on GPU, stored row-wise, dimension (LDDA,N)
             On entry, the matrix of column dimension N to which the row
             interchanges will be applied.
             On exit, the permuted matrix.
    
    \param[in]
    ldda     INTEGER
             The leading dimension of the array A. ldda >= n.
    
    \param[in]
    k1       INTEGER
             The first element of IPIV for which a row interchange will
             be done. (Fortran one-based index: 1 <= k1 .)
    
    \param[in]
    k2       INTEGER
             The last element of IPIV for which a row interchange will
             be done. (Fortran one-based index: 1 <= k2 .)
    
    \param[in]
    ipiv     INTEGER array, on CPU, dimension (K2*abs(INCI))
             The vector of pivot indices.  Only the elements in positions
             K1 through K2 of IPIV are accessed.
             IPIV(K) = L implies rows K and L are to be interchanged.
    
    \param[in]
    inci     INTEGER
             The increment between successive values of IPIV.
             Currently, INCI > 0.
             TODO: If INCI is negative, the pivots are applied in reverse order.

    @param[in]
    queue   magma_tally2_queue_t
            Queue to execute in.

    @ingroup magma_tally2_zaux2
    ********************************************************************/
// It is used in zgessm, zgetrf_incpiv.
extern "C" void
magma_tally2blas_zlaswp_q(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dAT, magma_tally2_int_t ldda,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci,
    magma_tally2_queue_t queue )
{
    #define dAT(i_, j_) (dAT + (i_)*ldda + (j_))
    
    magma_tally2_int_t info = 0;
    if ( n < 0 )
        info = -1;
    else if ( n > ldda )
        info = -3;
    else if ( k1 < 1 )
        info = -4;
    else if ( k2 < 1 )
        info = -5;
    else if ( inci <= 0 )
        info = -7;

    if (info != 0) {
        magma_tally2_xerbla( __func__, -(info) );
        return;  //info;
    }
  
    dim3 grid( (n + NTHREADS - 1) / NTHREADS );
    dim3 threads( NTHREADS );
    zlaswp_params_t params;
    
    for( int k = k1-1; k < k2; k += MAX_PIVOTS ) {
        int npivots = min( MAX_PIVOTS, k2-k );
        params.npivots = npivots;
        for( int j = 0; j < npivots; ++j ) {
            params.ipiv[j] = ipiv[(k+j)*inci] - k - 1;
        }
        zlaswp_kernel<<< grid, threads, 0, queue >>>( n, dAT(k,0), ldda, params );
    }
    
    #undef dAT
}


/**
    @see magma_tally2blas_zlaswp_q
    @ingroup magma_tally2_zaux2
    ********************************************************************/
extern "C" void
magma_tally2blas_zlaswp(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dAT, magma_tally2_int_t ldda,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci )
{
    magma_tally2blas_zlaswp_q( n, dAT, ldda, k1, k2, ipiv, inci, magma_tally2_stream );
}






// ------------------------------------------------------------
// Extended version has stride in both directions (ldx, ldy)
// to handle both row-wise and column-wise storage.

// Matrix A is stored row or column-wise in dA.
// Divide matrix A into block-columns of NTHREADS columns each.
// Each GPU block processes one block-column of A.
// Each thread goes down a column of A,
// swapping rows according to pivots stored in params.
__global__ void zlaswpx_kernel(
    int n, magma_tally2DoubleComplex *dA, int ldx, int ldy, zlaswp_params_t params )
{
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if ( tid < n ) {
        dA += tid*ldy;
        magma_tally2DoubleComplex *A1  = dA;
        
        for( int i1 = 0; i1 < params.npivots; ++i1 ) {
            int i2 = params.ipiv[i1];
            magma_tally2DoubleComplex *A2 = dA + i2*ldx;
            magma_tally2DoubleComplex temp = *A1;
            *A1 = *A2;
            *A2 = temp;
            A1 += ldx;  // A1 = dA + i1*ldx
        }
    }
}


/**
    Purpose:
    =============
    ZLASWPX performs a series of row interchanges on the matrix A.
    One row interchange is initiated for each of rows K1 through K2 of A.
    
    ** Unlike LAPACK, here A is stored either row-wise or column-wise,
       depending on ldx and ldy. **
    Otherwise, this is identical to LAPACK's interface.
    
    Arguments:
    ==========
    \param[in]
    n        INTEGER
             The number of columns of the matrix A.
    
    \param[in,out]
    dA       COMPLEX*16 array on GPU, dimension (*,*)
             On entry, the matrix of column dimension N to which the row
             interchanges will be applied.
             On exit, the permuted matrix.
    
    \param[in]
    ldx      INTEGER
             Stride between elements in same column.
    
    \param[in]
    ldy      INTEGER
             Stride between elements in same row.
             For A stored row-wise,    set ldx=ldda and ldy=1.
             For A stored column-wise, set ldx=1    and ldy=ldda.
    
    \param[in]
    k1       INTEGER
             The first element of IPIV for which a row interchange will
             be done. (One based index.)
    
    \param[in]
    k2       INTEGER
             The last element of IPIV for which a row interchange will
             be done. (One based index.)
    
    \param[in]
    ipiv     INTEGER array, on CPU, dimension (K2*abs(INCI))
             The vector of pivot indices.  Only the elements in positions
             K1 through K2 of IPIV are accessed.
             IPIV(K) = L implies rows K and L are to be interchanged.
    
    \param[in]
    inci     INTEGER
             The increment between successive values of IPIV.
             Currently, IPIV > 0.
             TODO: If IPIV is negative, the pivots are applied in reverse order.

    @param[in]
    queue   magma_tally2_queue_t
            Queue to execute in.

    @ingroup magma_tally2_zaux2
    ********************************************************************/
extern "C" void
magma_tally2blas_zlaswpx_q(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldx, magma_tally2_int_t ldy,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci,
    magma_tally2_queue_t queue )
{
    #define dA(i_, j_) (dA + (i_)*ldx + (j_)*ldy)
    
    magma_tally2_int_t info = 0;
    if ( n < 0 )
        info = -1;
    else if ( k1 < 0 )
        info = -4;  
    else if ( k2 < 0 || k2 < k1 )
        info = -5;
    else if ( inci <= 0 )
        info = -7;

    if (info != 0) {
        magma_tally2_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    dim3 grid( (n + NTHREADS - 1) / NTHREADS );
    dim3 threads( NTHREADS );
    zlaswp_params_t params;
    
    for( int k = k1-1; k < k2; k += MAX_PIVOTS ) {
        int npivots = min( MAX_PIVOTS, k2-k );
        params.npivots = npivots;
        for( int j = 0; j < npivots; ++j ) {
            params.ipiv[j] = ipiv[(k+j)*inci] - k - 1;
        }
        zlaswpx_kernel<<< grid, threads, 0, queue >>>( n, dA(k,0), ldx, ldy, params );
    }
    
    #undef dA
}


/**
    @see magma_tally2blas_zlaswpx_q
    @ingroup magma_tally2_zaux2
    ********************************************************************/
extern "C" void
magma_tally2blas_zlaswpx(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldx, magma_tally2_int_t ldy,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    const magma_tally2_int_t *ipiv, magma_tally2_int_t inci )
{
    return magma_tally2blas_zlaswpx_q( n, dA, ldx, ldy, k1, k2, ipiv, inci, magma_tally2_stream );
}






// ------------------------------------------------------------
// This version takes d_ipiv on the GPU. Thus it does not pass pivots
// as an argument using a structure, avoiding all the argument size
// limitations of CUDA and OpenCL. It also needs just one kernel launch
// with all the pivots, instead of multiple kernel launches with small
// batches of pivots. On Fermi, it is faster than magma_tally2blas_zlaswp
// (including copying pivots to the GPU).

__global__ void zlaswp2_kernel(
    int n, magma_tally2DoubleComplex *dAT, int ldda, int npivots,
    const magma_tally2_int_t *d_ipiv, int inci )
{
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if ( tid < n ) {
        dAT += tid;
        magma_tally2DoubleComplex *A1  = dAT;
        
        for( int i1 = 0; i1 < npivots; ++i1 ) {
            int i2 = d_ipiv[i1*inci] - 1;  // Fortran index
            magma_tally2DoubleComplex *A2 = dAT + i2*ldda;
            magma_tally2DoubleComplex temp = *A1;
            *A1 = *A2;
            *A2 = temp;
            A1 += ldda;  // A1 = dA + i1*ldx
        }
    }
}


/**
    Purpose:
    =============
    ZLASWP2 performs a series of row interchanges on the matrix A.
    One row interchange is initiated for each of rows K1 through K2 of A.
    
    ** Unlike LAPACK, here A is stored row-wise (hence dAT). **
    Otherwise, this is identical to LAPACK's interface.
    
    Here, d_ipiv is passed in GPU memory.
    
    Arguments:
    ==========
    \param[in]
    n        INTEGER
             The number of columns of the matrix A.
    
    \param[in,out]
    dAT      COMPLEX*16 array on GPU, stored row-wise, dimension (LDDA,*)
             On entry, the matrix of column dimension N to which the row
             interchanges will be applied.
             On exit, the permuted matrix.
    
    \param[in]
    ldda     INTEGER
             The leading dimension of the array A.
             (I.e., stride between elements in a column.)
    
    \param[in]
    k1       INTEGER
             The first element of IPIV for which a row interchange will
             be done. (One based index.)
    
    \param[in]
    k2       INTEGER
             The last element of IPIV for which a row interchange will
             be done. (One based index.)
    
    \param[in]
    d_ipiv   INTEGER array, on GPU, dimension (K2*abs(INCI))
             The vector of pivot indices.  Only the elements in positions
             K1 through K2 of IPIV are accessed.
             IPIV(K) = L implies rows K and L are to be interchanged.
    
    \param[in]
    inci     INTEGER
             The increment between successive values of IPIV.
             Currently, IPIV > 0.
             TODO: If IPIV is negative, the pivots are applied in reverse order.

    @param[in]
    queue   magma_tally2_queue_t
            Queue to execute in.

    @ingroup magma_tally2_zaux2
    ********************************************************************/
extern "C" void
magma_tally2blas_zlaswp2_q(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dAT, magma_tally2_int_t ldda,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    magma_tally2Int_const_ptr d_ipiv, magma_tally2_int_t inci,
    magma_tally2_queue_t queue )
{
    #define dAT(i_, j_) (dAT + (i_)*ldda + (j_))
    
    magma_tally2_int_t info = 0;
    if ( n < 0 )
        info = -1;
    else if ( k1 < 0 )
        info = -4;  
    else if ( k2 < 0 || k2 < k1 )
        info = -5;
    else if ( inci <= 0 )
        info = -7;

    if (info != 0) {
        magma_tally2_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    magma_tally2_int_t nb = k2-(k1-1);
    
    dim3 grid( (n + NTHREADS - 1) / NTHREADS );
    dim3 threads( NTHREADS );
    zlaswp2_kernel<<< grid, threads, 0, queue >>>
        ( n, dAT(k1-1,0), ldda, nb, d_ipiv, inci );
}


/**
    @see magma_tally2blas_zlaswp2_q
    @ingroup magma_tally2_zaux2
    ********************************************************************/
extern "C" void
magma_tally2blas_zlaswp2(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dAT, magma_tally2_int_t ldda,
    magma_tally2_int_t k1, magma_tally2_int_t k2,
    magma_tally2Int_const_ptr d_ipiv, magma_tally2_int_t inci )
{
    magma_tally2blas_zlaswp2_q( n, dAT, ldda, k1, k2, d_ipiv, inci, magma_tally2_stream );
}
