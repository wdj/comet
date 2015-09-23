/*
    -- MAGMA_minproduct (version 1.6.1) --
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
#define MAGMA_minproductBLAS_ZGEMM( name ) \
void magma_minproductblas_zgemm_##name( \
    magma_minproductDoubleComplex *C, const magma_minproductDoubleComplex *A, const magma_minproductDoubleComplex *B, \
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k, \
    magma_minproduct_int_t lda, magma_minproduct_int_t ldb, magma_minproduct_int_t ldc, \
    magma_minproductDoubleComplex alpha, magma_minproductDoubleComplex beta )

MAGMA_minproductBLAS_ZGEMM( a_0  );
MAGMA_minproductBLAS_ZGEMM( ab_0 );
MAGMA_minproductBLAS_ZGEMM( N_N_64_16_16_16_4_special );
MAGMA_minproductBLAS_ZGEMM( N_N_64_16_16_16_4         );
MAGMA_minproductBLAS_ZGEMM( N_T_64_16_4_16_4          );
MAGMA_minproductBLAS_ZGEMM( T_N_32_32_8_8_8           );
MAGMA_minproductBLAS_ZGEMM( T_T_64_16_16_16_4_special );
MAGMA_minproductBLAS_ZGEMM( T_T_64_16_16_16_4         );
                   
void magma_minproductblas_zgemm_tesla(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex alpha,
    const magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    const magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex *C, magma_minproduct_int_t ldc );

void magma_minproductblas_zgemv_tesla(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    const magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    const magma_minproductDoubleComplex *x, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex *y, magma_minproduct_int_t incy );


// kernels used in dznrm2, zgeqr2x-v4, laqps2_gpu, laqps3_gpu, zlarfbx, zlarfgx-v2, zlarfx
__global__ void
magma_minproduct_zgemv_kernel1(int m, const magma_minproductDoubleComplex * __restrict__ V, int ldv,
                    const magma_minproductDoubleComplex * __restrict__ c,
                    magma_minproductDoubleComplex *dwork);

__global__ void
magma_minproduct_zgemv_kernel2(int m, int n, const magma_minproductDoubleComplex * __restrict__ V, int ldv,
                    const magma_minproductDoubleComplex * __restrict__ x, magma_minproductDoubleComplex *c);

__global__ void
magma_minproduct_zgemv_kernel3(int m, const magma_minproductDoubleComplex * __restrict__ V, int ldv,
                    magma_minproductDoubleComplex *c, magma_minproductDoubleComplex *dwork,
                    magma_minproductDoubleComplex *tau);

__global__ void
magma_minproduct_ztrmv_tkernel(magma_minproductDoubleComplex *T, int ldt, magma_minproductDoubleComplex *v,
                                    magma_minproductDoubleComplex *y);

__global__ void
magma_minproduct_ztrmv_kernel2(const magma_minproductDoubleComplex *T, int ldt,
                    magma_minproductDoubleComplex *v, magma_minproductDoubleComplex *y, magma_minproductDoubleComplex *tau);

__global__ void
magma_minproduct_dznrm2_adjust_kernel(double *xnorm, magma_minproductDoubleComplex *c);


// kernels used in zhemv
__global__ void
zhemv_kernel_U(
    int n,
    magma_minproductDoubleComplex const * __restrict__ A, int lda,
    magma_minproductDoubleComplex const * __restrict__ x, int incx,
    magma_minproductDoubleComplex       * __restrict__ work);

__global__ void
zhemv_kernel_U_sum(
    int n,
    magma_minproductDoubleComplex alpha,
    int lda,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex       * __restrict__ y, int incy,
    magma_minproductDoubleComplex const * __restrict__ work );

// kernels used in zsymv
__global__ void
zsymv_kernel_U(
    int n,
    magma_minproductDoubleComplex const * __restrict__ A, int lda,
    magma_minproductDoubleComplex const * __restrict__ x, int incx,
    magma_minproductDoubleComplex       * __restrict__ work);

__global__ void
zsymv_kernel_U_sum(
    int n,
    magma_minproductDoubleComplex alpha,
    int lda,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex       * __restrict__ y, int incy,
    magma_minproductDoubleComplex const * __restrict__ work );

// kernels used in zhemv_mgpu
__global__ void
zhemv_kernel_U_mgpu(
    int n,
    magma_minproductDoubleComplex const * __restrict__ A, int lda,
    magma_minproductDoubleComplex const * __restrict__ x, int incx,
    magma_minproductDoubleComplex       * __restrict__ work,
    int my_gpu_id,
    int ngpu,
    int block_offset );

__global__ void
zhemv_kernel_U_mgpu_sum(
    int n,
    magma_minproductDoubleComplex alpha,
    int lda,
    magma_minproductDoubleComplex       * __restrict__ y, int incy,
    magma_minproductDoubleComplex const * __restrict__ work,
    int my_gpu_id,
    int ngpu,
    int block_offset);

#ifdef __cplusplus
}
#endif

#endif /* COMMONBLAS_Z_H */
