/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from commonblas_z.h normal z -> d, Fri Jan 30 19:00:10 2015
*/

#ifndef COMMONBLAS_D_H
#define COMMONBLAS_D_H

#ifdef __cplusplus
extern "C" {
#endif

/* ======================================================================
 * Internal prototypes
 */

// Tesla GEMM kernels
#define MAGMA_tally3BLAS_DGEMM( name ) \
void magma_tally3blas_dgemm_##name( \
    double *C, const double *A, const double *B, \
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k, \
    magma_tally3_int_t lda, magma_tally3_int_t ldb, magma_tally3_int_t ldc, \
    double alpha, double beta )

MAGMA_tally3BLAS_DGEMM( a_0  );
MAGMA_tally3BLAS_DGEMM( ab_0 );
MAGMA_tally3BLAS_DGEMM( N_N_64_16_16_16_4_special );
MAGMA_tally3BLAS_DGEMM( N_N_64_16_16_16_4         );
MAGMA_tally3BLAS_DGEMM( N_T_64_16_4_16_4          );
MAGMA_tally3BLAS_DGEMM( T_N_32_32_8_8_8           );
MAGMA_tally3BLAS_DGEMM( T_T_64_16_16_16_4_special );
MAGMA_tally3BLAS_DGEMM( T_T_64_16_16_16_4         );
                   
void magma_tally3blas_dgemm_tesla(
    magma_tally3_trans_t transA, magma_tally3_trans_t transB, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    double alpha,
    const double *A, magma_tally3_int_t lda,
    const double *B, magma_tally3_int_t ldb,
    double beta,
    double *C, magma_tally3_int_t ldc );

void magma_tally3blas_dgemv_tesla(
    magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n,
    double alpha,
    const double *A, magma_tally3_int_t lda,
    const double *x, magma_tally3_int_t incx,
    double beta,
    double *y, magma_tally3_int_t incy );


// kernels used in dnrm2, dgeqr2x-v4, laqps2_gpu, laqps3_gpu, dlarfbx, dlarfgx-v2, dlarfx
__global__ void
magma_tally3_dgemv_kernel1(int m, const double * __restrict__ V, int ldv,
                    const double * __restrict__ c,
                    double *dwork);

__global__ void
magma_tally3_dgemv_kernel2(int m, int n, const double * __restrict__ V, int ldv,
                    const double * __restrict__ x, double *c);

__global__ void
magma_tally3_dgemv_kernel3(int m, const double * __restrict__ V, int ldv,
                    double *c, double *dwork,
                    double *tau);

__global__ void
magma_tally3_dtrmv_tkernel(double *T, int ldt, double *v,
                                    double *y);

__global__ void
magma_tally3_dtrmv_kernel2(const double *T, int ldt,
                    double *v, double *y, double *tau);

__global__ void
magma_tally3_dnrm2_adjust_kernel(double *xnorm, double *c);


// kernels used in dsymv
__global__ void
dsymv_kernel_U(
    int n,
    double const * __restrict__ A, int lda,
    double const * __restrict__ x, int incx,
    double       * __restrict__ work);

__global__ void
dsymv_kernel_U_sum(
    int n,
    double alpha,
    int lda,
    double beta,
    double       * __restrict__ y, int incy,
    double const * __restrict__ work );

// kernels used in dsymv
__global__ void
dsymv_kernel_U(
    int n,
    double const * __restrict__ A, int lda,
    double const * __restrict__ x, int incx,
    double       * __restrict__ work);

__global__ void
dsymv_kernel_U_sum(
    int n,
    double alpha,
    int lda,
    double beta,
    double       * __restrict__ y, int incy,
    double const * __restrict__ work );

// kernels used in dsymv_mgpu
__global__ void
dsymv_kernel_U_mgpu(
    int n,
    double const * __restrict__ A, int lda,
    double const * __restrict__ x, int incx,
    double       * __restrict__ work,
    int my_gpu_id,
    int ngpu,
    int block_offset );

__global__ void
dsymv_kernel_U_mgpu_sum(
    int n,
    double alpha,
    int lda,
    double       * __restrict__ y, int incy,
    double const * __restrict__ work,
    int my_gpu_id,
    int ngpu,
    int block_offset);

#ifdef __cplusplus
}
#endif

#endif /* COMMONBLAS_D_H */
