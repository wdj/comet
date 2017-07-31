/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
*/

#ifndef COMMONBLAS_Z_H
#define COMMONBLAS_Z_H

#ifdef __cplusplus
extern "C" {
#endif

/* ======================================================================
 * Internal prototypes
 */

// Tesla GEMM kernels
#define MAGMA_tally2BLAS_ZGEMM( name ) \
void magma_tally2blas_zgemm_##name( \
    magma_tally2DoubleComplex *C, const magma_tally2DoubleComplex *A, const magma_tally2DoubleComplex *B, \
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k, \
    magma_tally2_int_t lda, magma_tally2_int_t ldb, magma_tally2_int_t ldc, \
    magma_tally2DoubleComplex alpha, magma_tally2DoubleComplex beta )

MAGMA_tally2BLAS_ZGEMM( a_0  );
MAGMA_tally2BLAS_ZGEMM( ab_0 );
MAGMA_tally2BLAS_ZGEMM( N_N_64_16_16_16_4_special );
MAGMA_tally2BLAS_ZGEMM( N_N_64_16_16_16_4         );
MAGMA_tally2BLAS_ZGEMM( N_T_64_16_4_16_4          );
MAGMA_tally2BLAS_ZGEMM( T_N_32_32_8_8_8           );
MAGMA_tally2BLAS_ZGEMM( T_T_64_16_16_16_4_special );
MAGMA_tally2BLAS_ZGEMM( T_T_64_16_16_16_4         );
                   
void magma_tally2blas_zgemm_tesla(
    magma_tally2_trans_t transA, magma_tally2_trans_t transB, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex alpha,
    const magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    const magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex *C, magma_tally2_int_t ldc );

void magma_tally2blas_zgemv_tesla(
    magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    const magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    const magma_tally2DoubleComplex *x, magma_tally2_int_t incx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex *y, magma_tally2_int_t incy );


// kernels used in dznrm2, zgeqr2x-v4, laqps2_gpu, laqps3_gpu, zlarfbx, zlarfgx-v2, zlarfx
__global__ void
magma_tally2_zgemv_kernel1(int m, const magma_tally2DoubleComplex * __restrict__ V, int ldv,
                    const magma_tally2DoubleComplex * __restrict__ c,
                    magma_tally2DoubleComplex *dwork);

__global__ void
magma_tally2_zgemv_kernel2(int m, int n, const magma_tally2DoubleComplex * __restrict__ V, int ldv,
                    const magma_tally2DoubleComplex * __restrict__ x, magma_tally2DoubleComplex *c);

__global__ void
magma_tally2_zgemv_kernel3(int m, const magma_tally2DoubleComplex * __restrict__ V, int ldv,
                    magma_tally2DoubleComplex *c, magma_tally2DoubleComplex *dwork,
                    magma_tally2DoubleComplex *tau);

__global__ void
magma_tally2_ztrmv_tkernel(magma_tally2DoubleComplex *T, int ldt, magma_tally2DoubleComplex *v,
                                    magma_tally2DoubleComplex *y);

__global__ void
magma_tally2_ztrmv_kernel2(const magma_tally2DoubleComplex *T, int ldt,
                    magma_tally2DoubleComplex *v, magma_tally2DoubleComplex *y, magma_tally2DoubleComplex *tau);

__global__ void
magma_tally2_dznrm2_adjust_kernel(double *xnorm, magma_tally2DoubleComplex *c);


// kernels used in zhemv
__global__ void
zhemv_kernel_U(
    int n,
    magma_tally2DoubleComplex const * __restrict__ A, int lda,
    magma_tally2DoubleComplex const * __restrict__ x, int incx,
    magma_tally2DoubleComplex       * __restrict__ work);

__global__ void
zhemv_kernel_U_sum(
    int n,
    magma_tally2DoubleComplex alpha,
    int lda,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex       * __restrict__ y, int incy,
    magma_tally2DoubleComplex const * __restrict__ work );

// kernels used in zsymv
__global__ void
zsymv_kernel_U(
    int n,
    magma_tally2DoubleComplex const * __restrict__ A, int lda,
    magma_tally2DoubleComplex const * __restrict__ x, int incx,
    magma_tally2DoubleComplex       * __restrict__ work);

__global__ void
zsymv_kernel_U_sum(
    int n,
    magma_tally2DoubleComplex alpha,
    int lda,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex       * __restrict__ y, int incy,
    magma_tally2DoubleComplex const * __restrict__ work );

// kernels used in zhemv_mgpu
__global__ void
zhemv_kernel_U_mgpu(
    int n,
    magma_tally2DoubleComplex const * __restrict__ A, int lda,
    magma_tally2DoubleComplex const * __restrict__ x, int incx,
    magma_tally2DoubleComplex       * __restrict__ work,
    int my_gpu_id,
    int ngpu,
    int block_offset );

__global__ void
zhemv_kernel_U_mgpu_sum(
    int n,
    magma_tally2DoubleComplex alpha,
    int lda,
    magma_tally2DoubleComplex       * __restrict__ y, int incy,
    magma_tally2DoubleComplex const * __restrict__ work,
    int my_gpu_id,
    int ngpu,
    int block_offset);

#ifdef __cplusplus
}
#endif

#endif /* COMMONBLAS_Z_H */
