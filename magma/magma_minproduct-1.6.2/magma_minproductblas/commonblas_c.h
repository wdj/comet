/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from commonblas_z.h normal z -> c, Fri Jan 30 19:00:10 2015
*/

#ifndef COMMONBLAS_C_H
#define COMMONBLAS_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* ======================================================================
 * Internal prototypes
 */

// Tesla GEMM kernels
#define MAGMA_minproductBLAS_CGEMM( name ) \
void magma_minproductblas_cgemm_##name( \
    magma_minproductFloatComplex *C, const magma_minproductFloatComplex *A, const magma_minproductFloatComplex *B, \
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k, \
    magma_minproduct_int_t lda, magma_minproduct_int_t ldb, magma_minproduct_int_t ldc, \
    magma_minproductFloatComplex alpha, magma_minproductFloatComplex beta )

MAGMA_minproductBLAS_CGEMM( a_0  );
MAGMA_minproductBLAS_CGEMM( ab_0 );
MAGMA_minproductBLAS_CGEMM( N_N_64_16_16_16_4_special );
MAGMA_minproductBLAS_CGEMM( N_N_64_16_16_16_4         );
MAGMA_minproductBLAS_CGEMM( N_T_64_16_4_16_4          );
MAGMA_minproductBLAS_CGEMM( T_N_32_32_8_8_8           );
MAGMA_minproductBLAS_CGEMM( T_T_64_16_16_16_4_special );
MAGMA_minproductBLAS_CGEMM( T_T_64_16_16_16_4         );
                   
void magma_minproductblas_cgemm_tesla(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    const magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    const magma_minproductFloatComplex *B, magma_minproduct_int_t ldb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex *C, magma_minproduct_int_t ldc );

void magma_minproductblas_cgemv_tesla(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    const magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    const magma_minproductFloatComplex *x, magma_minproduct_int_t incx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex *y, magma_minproduct_int_t incy );


// kernels used in scnrm2, cgeqr2x-v4, laqps2_gpu, laqps3_gpu, clarfbx, clarfgx-v2, clarfx
__global__ void
magma_minproduct_cgemv_kernel1(int m, const magma_minproductFloatComplex * __restrict__ V, int ldv,
                    const magma_minproductFloatComplex * __restrict__ c,
                    magma_minproductFloatComplex *dwork);

__global__ void
magma_minproduct_cgemv_kernel2(int m, int n, const magma_minproductFloatComplex * __restrict__ V, int ldv,
                    const magma_minproductFloatComplex * __restrict__ x, magma_minproductFloatComplex *c);

__global__ void
magma_minproduct_cgemv_kernel3(int m, const magma_minproductFloatComplex * __restrict__ V, int ldv,
                    magma_minproductFloatComplex *c, magma_minproductFloatComplex *dwork,
                    magma_minproductFloatComplex *tau);

__global__ void
magma_minproduct_ctrmv_tkernel(magma_minproductFloatComplex *T, int ldt, magma_minproductFloatComplex *v,
                                    magma_minproductFloatComplex *y);

__global__ void
magma_minproduct_ctrmv_kernel2(const magma_minproductFloatComplex *T, int ldt,
                    magma_minproductFloatComplex *v, magma_minproductFloatComplex *y, magma_minproductFloatComplex *tau);

__global__ void
magma_minproduct_scnrm2_adjust_kernel(float *xnorm, magma_minproductFloatComplex *c);


// kernels used in chemv
__global__ void
chemv_kernel_U(
    int n,
    magma_minproductFloatComplex const * __restrict__ A, int lda,
    magma_minproductFloatComplex const * __restrict__ x, int incx,
    magma_minproductFloatComplex       * __restrict__ work);

__global__ void
chemv_kernel_U_sum(
    int n,
    magma_minproductFloatComplex alpha,
    int lda,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex       * __restrict__ y, int incy,
    magma_minproductFloatComplex const * __restrict__ work );

// kernels used in csymv
__global__ void
csymv_kernel_U(
    int n,
    magma_minproductFloatComplex const * __restrict__ A, int lda,
    magma_minproductFloatComplex const * __restrict__ x, int incx,
    magma_minproductFloatComplex       * __restrict__ work);

__global__ void
csymv_kernel_U_sum(
    int n,
    magma_minproductFloatComplex alpha,
    int lda,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex       * __restrict__ y, int incy,
    magma_minproductFloatComplex const * __restrict__ work );

// kernels used in chemv_mgpu
__global__ void
chemv_kernel_U_mgpu(
    int n,
    magma_minproductFloatComplex const * __restrict__ A, int lda,
    magma_minproductFloatComplex const * __restrict__ x, int incx,
    magma_minproductFloatComplex       * __restrict__ work,
    int my_gpu_id,
    int ngpu,
    int block_offset );

__global__ void
chemv_kernel_U_mgpu_sum(
    int n,
    magma_minproductFloatComplex alpha,
    int lda,
    magma_minproductFloatComplex       * __restrict__ y, int incy,
    magma_minproductFloatComplex const * __restrict__ work,
    int my_gpu_id,
    int ngpu,
    int block_offset);

#ifdef __cplusplus
}
#endif

#endif /* COMMONBLAS_C_H */
