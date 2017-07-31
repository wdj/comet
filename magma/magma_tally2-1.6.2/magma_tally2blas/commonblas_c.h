/*
    -- MAGMA_tally2 (version 1.6.1) --
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
#define MAGMA_tally2BLAS_CGEMM( name ) \
void magma_tally2blas_cgemm_##name( \
    magma_tally2FloatComplex *C, const magma_tally2FloatComplex *A, const magma_tally2FloatComplex *B, \
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k, \
    magma_tally2_int_t lda, magma_tally2_int_t ldb, magma_tally2_int_t ldc, \
    magma_tally2FloatComplex alpha, magma_tally2FloatComplex beta )

MAGMA_tally2BLAS_CGEMM( a_0  );
MAGMA_tally2BLAS_CGEMM( ab_0 );
MAGMA_tally2BLAS_CGEMM( N_N_64_16_16_16_4_special );
MAGMA_tally2BLAS_CGEMM( N_N_64_16_16_16_4         );
MAGMA_tally2BLAS_CGEMM( N_T_64_16_4_16_4          );
MAGMA_tally2BLAS_CGEMM( T_N_32_32_8_8_8           );
MAGMA_tally2BLAS_CGEMM( T_T_64_16_16_16_4_special );
MAGMA_tally2BLAS_CGEMM( T_T_64_16_16_16_4         );
                   
void magma_tally2blas_cgemm_tesla(
    magma_tally2_trans_t transA, magma_tally2_trans_t transB, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex alpha,
    const magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    const magma_tally2FloatComplex *B, magma_tally2_int_t ldb,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex *C, magma_tally2_int_t ldc );

void magma_tally2blas_cgemv_tesla(
    magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex alpha,
    const magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    const magma_tally2FloatComplex *x, magma_tally2_int_t incx,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex *y, magma_tally2_int_t incy );


// kernels used in scnrm2, cgeqr2x-v4, laqps2_gpu, laqps3_gpu, clarfbx, clarfgx-v2, clarfx
__global__ void
magma_tally2_cgemv_kernel1(int m, const magma_tally2FloatComplex * __restrict__ V, int ldv,
                    const magma_tally2FloatComplex * __restrict__ c,
                    magma_tally2FloatComplex *dwork);

__global__ void
magma_tally2_cgemv_kernel2(int m, int n, const magma_tally2FloatComplex * __restrict__ V, int ldv,
                    const magma_tally2FloatComplex * __restrict__ x, magma_tally2FloatComplex *c);

__global__ void
magma_tally2_cgemv_kernel3(int m, const magma_tally2FloatComplex * __restrict__ V, int ldv,
                    magma_tally2FloatComplex *c, magma_tally2FloatComplex *dwork,
                    magma_tally2FloatComplex *tau);

__global__ void
magma_tally2_ctrmv_tkernel(magma_tally2FloatComplex *T, int ldt, magma_tally2FloatComplex *v,
                                    magma_tally2FloatComplex *y);

__global__ void
magma_tally2_ctrmv_kernel2(const magma_tally2FloatComplex *T, int ldt,
                    magma_tally2FloatComplex *v, magma_tally2FloatComplex *y, magma_tally2FloatComplex *tau);

__global__ void
magma_tally2_scnrm2_adjust_kernel(float *xnorm, magma_tally2FloatComplex *c);


// kernels used in chemv
__global__ void
chemv_kernel_U(
    int n,
    magma_tally2FloatComplex const * __restrict__ A, int lda,
    magma_tally2FloatComplex const * __restrict__ x, int incx,
    magma_tally2FloatComplex       * __restrict__ work);

__global__ void
chemv_kernel_U_sum(
    int n,
    magma_tally2FloatComplex alpha,
    int lda,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex       * __restrict__ y, int incy,
    magma_tally2FloatComplex const * __restrict__ work );

// kernels used in csymv
__global__ void
csymv_kernel_U(
    int n,
    magma_tally2FloatComplex const * __restrict__ A, int lda,
    magma_tally2FloatComplex const * __restrict__ x, int incx,
    magma_tally2FloatComplex       * __restrict__ work);

__global__ void
csymv_kernel_U_sum(
    int n,
    magma_tally2FloatComplex alpha,
    int lda,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex       * __restrict__ y, int incy,
    magma_tally2FloatComplex const * __restrict__ work );

// kernels used in chemv_mgpu
__global__ void
chemv_kernel_U_mgpu(
    int n,
    magma_tally2FloatComplex const * __restrict__ A, int lda,
    magma_tally2FloatComplex const * __restrict__ x, int incx,
    magma_tally2FloatComplex       * __restrict__ work,
    int my_gpu_id,
    int ngpu,
    int block_offset );

__global__ void
chemv_kernel_U_mgpu_sum(
    int n,
    magma_tally2FloatComplex alpha,
    int lda,
    magma_tally2FloatComplex       * __restrict__ y, int incy,
    magma_tally2FloatComplex const * __restrict__ work,
    int my_gpu_id,
    int ngpu,
    int block_offset);

#ifdef __cplusplus
}
#endif

#endif /* COMMONBLAS_C_H */
