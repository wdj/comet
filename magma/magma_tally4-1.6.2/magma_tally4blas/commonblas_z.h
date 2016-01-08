/*
    -- MAGMA_tally4 (version 1.6.1) --
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
#define MAGMA_tally4BLAS_ZGEMM( name ) \
void magma_tally4blas_zgemm_##name( \
    magma_tally4DoubleComplex *C, const magma_tally4DoubleComplex *A, const magma_tally4DoubleComplex *B, \
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k, \
    magma_tally4_int_t lda, magma_tally4_int_t ldb, magma_tally4_int_t ldc, \
    magma_tally4DoubleComplex alpha, magma_tally4DoubleComplex beta )

MAGMA_tally4BLAS_ZGEMM( a_0  );
MAGMA_tally4BLAS_ZGEMM( ab_0 );
MAGMA_tally4BLAS_ZGEMM( N_N_64_16_16_16_4_special );
MAGMA_tally4BLAS_ZGEMM( N_N_64_16_16_16_4         );
MAGMA_tally4BLAS_ZGEMM( N_T_64_16_4_16_4          );
MAGMA_tally4BLAS_ZGEMM( T_N_32_32_8_8_8           );
MAGMA_tally4BLAS_ZGEMM( T_T_64_16_16_16_4_special );
MAGMA_tally4BLAS_ZGEMM( T_T_64_16_16_16_4         );
                   
void magma_tally4blas_zgemm_tesla(
    magma_tally4_trans_t transA, magma_tally4_trans_t transB, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex alpha,
    const magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    const magma_tally4DoubleComplex *B, magma_tally4_int_t ldb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex *C, magma_tally4_int_t ldc );

void magma_tally4blas_zgemv_tesla(
    magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    const magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    const magma_tally4DoubleComplex *x, magma_tally4_int_t incx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex *y, magma_tally4_int_t incy );


// kernels used in dznrm2, zgeqr2x-v4, laqps2_gpu, laqps3_gpu, zlarfbx, zlarfgx-v2, zlarfx
__global__ void
magma_tally4_zgemv_kernel1(int m, const magma_tally4DoubleComplex * __restrict__ V, int ldv,
                    const magma_tally4DoubleComplex * __restrict__ c,
                    magma_tally4DoubleComplex *dwork);

__global__ void
magma_tally4_zgemv_kernel2(int m, int n, const magma_tally4DoubleComplex * __restrict__ V, int ldv,
                    const magma_tally4DoubleComplex * __restrict__ x, magma_tally4DoubleComplex *c);

__global__ void
magma_tally4_zgemv_kernel3(int m, const magma_tally4DoubleComplex * __restrict__ V, int ldv,
                    magma_tally4DoubleComplex *c, magma_tally4DoubleComplex *dwork,
                    magma_tally4DoubleComplex *tau);

__global__ void
magma_tally4_ztrmv_tkernel(magma_tally4DoubleComplex *T, int ldt, magma_tally4DoubleComplex *v,
                                    magma_tally4DoubleComplex *y);

__global__ void
magma_tally4_ztrmv_kernel2(const magma_tally4DoubleComplex *T, int ldt,
                    magma_tally4DoubleComplex *v, magma_tally4DoubleComplex *y, magma_tally4DoubleComplex *tau);

__global__ void
magma_tally4_dznrm2_adjust_kernel(double *xnorm, magma_tally4DoubleComplex *c);


// kernels used in zhemv
__global__ void
zhemv_kernel_U(
    int n,
    magma_tally4DoubleComplex const * __restrict__ A, int lda,
    magma_tally4DoubleComplex const * __restrict__ x, int incx,
    magma_tally4DoubleComplex       * __restrict__ work);

__global__ void
zhemv_kernel_U_sum(
    int n,
    magma_tally4DoubleComplex alpha,
    int lda,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex       * __restrict__ y, int incy,
    magma_tally4DoubleComplex const * __restrict__ work );

// kernels used in zsymv
__global__ void
zsymv_kernel_U(
    int n,
    magma_tally4DoubleComplex const * __restrict__ A, int lda,
    magma_tally4DoubleComplex const * __restrict__ x, int incx,
    magma_tally4DoubleComplex       * __restrict__ work);

__global__ void
zsymv_kernel_U_sum(
    int n,
    magma_tally4DoubleComplex alpha,
    int lda,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex       * __restrict__ y, int incy,
    magma_tally4DoubleComplex const * __restrict__ work );

// kernels used in zhemv_mgpu
__global__ void
zhemv_kernel_U_mgpu(
    int n,
    magma_tally4DoubleComplex const * __restrict__ A, int lda,
    magma_tally4DoubleComplex const * __restrict__ x, int incx,
    magma_tally4DoubleComplex       * __restrict__ work,
    int my_gpu_id,
    int ngpu,
    int block_offset );

__global__ void
zhemv_kernel_U_mgpu_sum(
    int n,
    magma_tally4DoubleComplex alpha,
    int lda,
    magma_tally4DoubleComplex       * __restrict__ y, int incy,
    magma_tally4DoubleComplex const * __restrict__ work,
    int my_gpu_id,
    int ngpu,
    int block_offset);

#ifdef __cplusplus
}
#endif

#endif /* COMMONBLAS_Z_H */
