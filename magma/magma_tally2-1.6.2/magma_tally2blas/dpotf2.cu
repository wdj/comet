/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @generated from zpotf2.cu normal z -> d, Fri Jan 30 19:00:10 2015
*/
#include "common_magma_tally2.h"

#define PRECISION_d

//#if (GPUSHMEM < 200)
#define ddot_max_bs 512  // 512 is max threads for 1.x cards
//#else
//#define ddot_max_bs 1024
//#endif

void dpotf2_dscal(magma_tally2_int_t n, double *x, magma_tally2_int_t incx);
void dpotf2_ddot(magma_tally2_int_t n, double *x, magma_tally2_int_t incx);

#if defined(PRECISION_z) || defined(PRECISION_c)
void dlacgv(magma_tally2_int_t n, double *x, magma_tally2_int_t incx);
#endif

/**
    Purpose
    -------

    dpotf2 computes the Cholesky factorization of a real symmetric
    positive definite matrix A.

    The factorization has the form
        A = U**H * U,  if UPLO = Magma_tally2Upper, or
        A = L  * L**H, if UPLO = Magma_tally2Lower,
    where U is an upper triangular matrix and L is lower triangular.

    This is the unblocked version of the algorithm, calling Level 2 BLAS.

    Arguments
    ---------

    @param[in]
    uplo    magma_tally2_uplo_t
            Specifies whether the upper or lower triangular part of the
            symmetric matrix A is stored.
      -     = Magma_tally2Upper:  Upper triangular
      -     = Magma_tally2Lower:  Lower triangular

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0 and N <= 512.

    @param[in,out]
    dA      DOUBLE_PRECISION array, dimension (LDDA,N)
            On entry, the symmetric matrix A.  If UPLO = Magma_tally2Upper, the leading
            n by n upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_tally2Lower, the
            leading n by n lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
    \n
            On exit, if INFO = 0, the factor U or L from the Cholesky
            factorization A = U**H * U  or A = L * L**H.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0: successful exit
      -     < 0: if INFO = -k, the k-th argument had an illegal value
      -     > 0: if INFO = k, the leading minor of order k is not
                 positive definite, and the factorization could not be
                 completed.

    @ingroup magma_tally2_dposv_aux
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_dpotf2_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info )
{
#define dA(i_, j_)  (dA + (i_) + (j_)*ldda)

    magma_tally2_int_t j;

    *info = 0;
    if ( uplo != Magma_tally2Upper && uplo != Magma_tally2Lower) {
        *info = -1;
    } else if (n < 0 || n > ddot_max_bs) {
        *info = -2;
    } else if (ldda < max(1,n)) {
        *info = -4;
    }

    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }

    // Quick return if possible
    if (n == 0) {
        return *info;
    }

    double alpha = MAGMA_tally2_D_NEG_ONE;
    double beta  = MAGMA_tally2_D_ONE;

    if (uplo == Magma_tally2Upper) {
        for(j = 0; j < n; j++) {
            dpotf2_ddot(j, dA(0,j), 1); // including ddot product and update a(j,j)
            if (j < n) {
                #if defined(PRECISION_z) || defined(PRECISION_c)
                dlacgv(j, dA(0, j), 1);
                #endif
                magma_tally2_dgemv( Magma_tally2Trans, j, n-j-1,
                             alpha, dA(0, j+1), ldda,
                                    dA(0, j),   1,
                             beta,  dA(j, j+1), ldda);

                #if defined(PRECISION_z) || defined(PRECISION_c)
                dlacgv(j, dA(0, j), 1);
                #endif
                dpotf2_dscal(n-j, dA(j,j), ldda);
            }
        }
    }
    else {
        for(j = 0; j < n; j++) {
            dpotf2_ddot(j, dA(j,0), ldda); // including ddot product and update a(j,j)
            if (j < n) {
                #if defined(PRECISION_z) || defined(PRECISION_c)
                dlacgv(j, dA(j, 0), ldda);
                #endif
                magma_tally2_dgemv( Magma_tally2NoTrans, n-j-1, j,
                             alpha, dA(j+1, 0), ldda,
                                    dA(j,0),    ldda,
                             beta,  dA(j+1, j), 1 );

                #if defined(PRECISION_z) || defined(PRECISION_c)
                dlacgv(j, dA(j, 0), ldda);
                #endif
                dpotf2_dscal(n-j, dA(j,j), 1);
            }
        }
    }

    return *info;
}

#define dscal_bs  32
#define ddot_bs  512
#define dlacgv_bs 512

// dynamically allocated shared memory, set to size number of threads when the kernel is launched.
// See CUDA Guide B.2.3
extern __shared__ double shared_data[];

__global__ void kernel_ddot(int n, double *x, int incx, int threadSize)
{
    int tx = threadIdx.x;

    double *sdata = shared_data;

    double res = MAGMA_tally2_D_ZERO;

    if (tx < n) {
       res = x[tx*incx];
    }

    sdata[tx] = MAGMA_tally2_D_REAL(res * MAGMA_tally2_D_CNJG(res));

    __syncthreads();

    for(int s = blockDim.x/2; s > 32; s >>= 1 ) {
        if (tx < s) {
            sdata[tx] += sdata[tx+s];
        }
        __syncthreads();
    }

    if (tx < 32) {
        volatile double* smem = sdata;
        smem[tx] += smem[tx+32];
        smem[tx] += smem[tx+16];
        smem[tx] += smem[tx+8];
        smem[tx] += smem[tx+4];
        smem[tx] += smem[tx+2];
        smem[tx] += smem[tx+1];
    }

    if (tx == 0) {
        double xreal = MAGMA_tally2_D_REAL(x[n*incx]);
        x[n*incx] = MAGMA_tally2_D_MAKE( sqrt(xreal - sdata[0]), 0 );
    }
}

void dpotf2_ddot(magma_tally2_int_t n, double *x, magma_tally2_int_t incx)
{
/*
    Specialized Ddot
    1) performs ddot sum = x[0:n-1]*conj(x[0:n-1])
    2) updates x[n] = sqrt(x[n]-sum);

*/
    if (n > ddot_max_bs) {
        fprintf( stderr, "n = %d > %d is not supported in dpotf2_ddot\n", (int) n, (int) ddot_max_bs);
        return;
    }
    int threadSize;

    if (n <= 1024 && n > 512) {
        threadSize = 1024;
    }
    else if (n <= 512 && n > 256 ) {
        threadSize = 512;
    }
    else if (n <= 256 && n > 128) {
        threadSize = 256;
    }
    else if (n <= 128 && n > 64) {
        threadSize = 128;
    }
    else {
        threadSize = 64;
    }

    kernel_ddot<<< 1, threadSize, threadSize * sizeof(double), magma_tally2_stream>>> (n, x, incx, threadSize);
}

__global__ void kernel_dscal(int n, double *x, int incx)
{
    int id = blockIdx.x * dscal_bs + threadIdx.x;

    __shared__ double factor;

    if (threadIdx.x == 0) {
        factor = MAGMA_tally2_D_MAKE(1.0/MAGMA_tally2_D_REAL(x[0]), 0.0);
    }

    __syncthreads();

    if ( id < n && id >0) {
        x[id*incx] = x[id*incx] * factor;
    }
}


void dpotf2_dscal(magma_tally2_int_t n, double *x, magma_tally2_int_t incx)
{
/*
    Specialized Dscal perform x[1:n-1]/x[0]

*/
    dim3 threads(dscal_bs, 1, 1);
    int num_blocks = (n - 1)/dscal_bs + 1;
    dim3 grid(num_blocks,1);
    kernel_dscal<<< grid, threads, 0, magma_tally2_stream >>> (n, x, incx);
}


#if defined(PRECISION_z) || defined(PRECISION_c)

__global__ void kernel_dlacgv(int n, double *x, int incx)
{
    int id = blockIdx.x * dlacgv_bs + threadIdx.x;

    if ( id < n ) {
        x[id*incx] = MAGMA_tally2_D_CNJG(x[id*incx]);
    }
}


/**
    Purpose
    -------

    DLACGV conjugates a real vector of length N.

    Arguments
    ---------

    @param[in]
    n       INTEGER
            The length of the vector X.  N >= 0.

    @param[in,out]
    x       DOUBLE PRECISION array, dimension
                           (1+(N-1)*abs(INCX))
            On entry, the vector of length N to be conjugated.
            On exit, X is overwritten with conjg(X).

    @param[in]
    incx    INTEGER
            The spacing between successive elements of X.

    @ingroup magma_tally2_dposv_aux
    ********************************************************************/
void dlacgv(magma_tally2_int_t n, double *x, magma_tally2_int_t incx)
{
    dim3 threads(dlacgv_bs, 1, 1);
    int num_blocks = (n - 1)/dlacgv_bs + 1;
    dim3 grid(num_blocks,1);
    kernel_dlacgv<<< grid, threads, 0, magma_tally2_stream >>> (n, x, incx);
}

#endif // defined(PRECISION_z) || defined(PRECISION_c)
