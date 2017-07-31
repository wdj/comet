/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from commonblas_z.h normal z -> s, Fri Jan 30 19:00:10 2015
*/

#ifndef COMMONBLAS_S_H
#define COMMONBLAS_S_H

#ifdef __cplusplus
extern "C" {
#endif

/* ======================================================================
 * Internal prototypes
 */

// Tesla GEMM kernels
#define MAGMA_tally2BLAS_SGEMM( name ) \
void magma_tally2blas_sgemm_##name( \
    float *C, const float *A, const float *B, \
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k, \
    magma_tally2_int_t lda, magma_tally2_int_t ldb, magma_tally2_int_t ldc, \
    float alpha, float beta )

MAGMA_tally2BLAS_SGEMM( a_0  );
MAGMA_tally2BLAS_SGEMM( ab_0 );
MAGMA_tally2BLAS_SGEMM( N_N_64_16_16_16_4_special );
MAGMA_tally2BLAS_SGEMM( N_N_64_16_16_16_4         );
MAGMA_tally2BLAS_SGEMM( N_T_64_16_4_16_4          );
MAGMA_tally2BLAS_SGEMM( T_N_32_32_8_8_8           );
MAGMA_tally2BLAS_SGEMM( T_T_64_16_16_16_4_special );
MAGMA_tally2BLAS_SGEMM( T_T_64_16_16_16_4         );
                   
void magma_tally2blas_sgemm_tesla(
    magma_tally2_trans_t transA, magma_tally2_trans_t transB, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    float alpha,
    const float *A, magma_tally2_int_t lda,
    const float *B, magma_tally2_int_t ldb,
    float beta,
    float *C, magma_tally2_int_t ldc );

void magma_tally2blas_sgemv_tesla(
    magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n,
    float alpha,
    const float *A, magma_tally2_int_t lda,
    const float *x, magma_tally2_int_t incx,
    float beta,
    float *y, magma_tally2_int_t incy );


// kernels used in snrm2, sgeqr2x-v4, laqps2_gpu, laqps3_gpu, slarfbx, slarfgx-v2, slarfx
__global__ void
magma_tally2_sgemv_kernel1(int m, const float * __restrict__ V, int ldv,
                    const float * __restrict__ c,
                    float *dwork);

__global__ void
magma_tally2_sgemv_kernel2(int m, int n, const float * __restrict__ V, int ldv,
                    const float * __restrict__ x, float *c);

__global__ void
magma_tally2_sgemv_kernel3(int m, const float * __restrict__ V, int ldv,
                    float *c, float *dwork,
                    float *tau);

__global__ void
magma_tally2_strmv_tkernel(float *T, int ldt, float *v,
                                    float *y);

__global__ void
magma_tally2_strmv_kernel2(const float *T, int ldt,
                    float *v, float *y, float *tau);

__global__ void
magma_tally2_snrm2_adjust_kernel(float *xnorm, float *c);


// kernels used in ssymv
__global__ void
ssymv_kernel_U(
    int n,
    float const * __restrict__ A, int lda,
    float const * __restrict__ x, int incx,
    float       * __restrict__ work);

__global__ void
ssymv_kernel_U_sum(
    int n,
    float alpha,
    int lda,
    float beta,
    float       * __restrict__ y, int incy,
    float const * __restrict__ work );

// kernels used in ssymv
__global__ void
ssymv_kernel_U(
    int n,
    float const * __restrict__ A, int lda,
    float const * __restrict__ x, int incx,
    float       * __restrict__ work);

__global__ void
ssymv_kernel_U_sum(
    int n,
    float alpha,
    int lda,
    float beta,
    float       * __restrict__ y, int incy,
    float const * __restrict__ work );

// kernels used in ssymv_mgpu
__global__ void
ssymv_kernel_U_mgpu(
    int n,
    float const * __restrict__ A, int lda,
    float const * __restrict__ x, int incx,
    float       * __restrict__ work,
    int my_gpu_id,
    int ngpu,
    int block_offset );

__global__ void
ssymv_kernel_U_mgpu_sum(
    int n,
    float alpha,
    int lda,
    float       * __restrict__ y, int incy,
    float const * __restrict__ work,
    int my_gpu_id,
    int ngpu,
    int block_offset);

#ifdef __cplusplus
}
#endif

#endif /* COMMONBLAS_S_H */
