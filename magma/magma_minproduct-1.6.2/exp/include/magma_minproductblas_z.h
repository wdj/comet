/*
 *   -- MAGMA_minproduct (version 1.6.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      @date January 2015
 *
 * @precisions normal z -> s d c
 */

#ifndef _MAGMA_minproductBLAS_Z_H_
#define _MAGMA_minproductBLAS_Z_H_

#define PRECISION_z
#ifdef __cplusplus
extern "C" {
#endif

  /*
   * Interface to clean
   */
double cpu_gpu_zdiff(             int M, int N, 
                  cuDoubleComplex * a, int lda, 
                  cuDoubleComplex *da, int ldda);
void zzero_32x32_block(           cuDoubleComplex *, magma_minproduct_int_t);
void zzero_nbxnb_block(           magma_minproduct_int_t, cuDoubleComplex *, magma_minproduct_int_t);
void magma_minproductblas_zinplace_transpose(cuDoubleComplex *, magma_minproduct_int_t, magma_minproduct_int_t);
void magma_minproductblas_zpermute_long(     cuDoubleComplex *, magma_minproduct_int_t, 
                  magma_minproduct_int_t *, magma_minproduct_int_t, magma_minproduct_int_t);
void magma_minproductblas_zpermute_long2(    cuDoubleComplex *, magma_minproduct_int_t, 
                  magma_minproduct_int_t *, magma_minproduct_int_t, magma_minproduct_int_t);
void magma_minproductblas_ztranspose(        cuDoubleComplex *, magma_minproduct_int_t, 
                  cuDoubleComplex *, magma_minproduct_int_t, 
                  magma_minproduct_int_t, magma_minproduct_int_t);
void magma_minproductblas_ztranspose2(       cuDoubleComplex *, magma_minproduct_int_t, 
                  cuDoubleComplex *, magma_minproduct_int_t, 
                  magma_minproduct_int_t, magma_minproduct_int_t);
void magma_minproductblas_zhtodt(            cuDoubleComplex  *ha, int lda,
                          cuDoubleComplex *dat, int ldda,
                          cuDoubleComplex  *dB, int lddb,
                                  int m, int n , int nb);
void magma_minproductblas_zdtoht(            cuDoubleComplex *dat, int ldda,
                  cuDoubleComplex  *ha, int lda,
                       cuDoubleComplex  *dB, int lddb,
                  int m, int n , int nb);
  
  /*
   * LAPACK auxiliary functions
   */
void   magma_minproductblas_zlacpy( char uplo, 
             magma_minproduct_int_t m, magma_minproduct_int_t n, 
             cuDoubleComplex *A, magma_minproduct_int_t lda, 
             cuDoubleComplex *B, magma_minproduct_int_t ldb);
double magma_minproductblas_zlange( char norm, 
             magma_minproduct_int_t m, magma_minproduct_int_t n, 
             cuDoubleComplex *A, magma_minproduct_int_t lda, double *WORK);
double magma_minproductblas_zlanhe( char norm, char uplo, 
             magma_minproduct_int_t n,
             cuDoubleComplex *A, magma_minproduct_int_t lda, double *WORK);
double magma_minproductblas_zlansy( char norm, char uplo,
             magma_minproduct_int_t n, 
             cuDoubleComplex *A, magma_minproduct_int_t lda, double *WORK);
void   magma_minproductblas_zlascl( char type, int kl, int ku,
             double cfrom, double cto,
             int m, int n,
             cuDoubleComplex *A, int lda, int *info );
void   magma_minproductblas_zlaset( magma_minproduct_int_t m, magma_minproduct_int_t n,
             cuDoubleComplex *A, magma_minproduct_int_t lda);
void   magma_minproductblas_zlaswp( magma_minproduct_int_t N, 
             cuDoubleComplex *dAT, magma_minproduct_int_t lda, 
             magma_minproduct_int_t i1,  magma_minproduct_int_t i2, 
             magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci );
void   magma_minproductblas_zlaswpx(magma_minproduct_int_t N, 
             cuDoubleComplex *dAT, magma_minproduct_int_t ldx, magma_minproduct_int_t ldy, 
             magma_minproduct_int_t i1, magma_minproduct_int_t i2,
             magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci );

  /*
   * Level 1 BLAS
   */
void   magma_minproductblas_zswap(   magma_minproduct_int_t N, 
              cuDoubleComplex *dA1, magma_minproduct_int_t lda1, 
              cuDoubleComplex *dA2, magma_minproduct_int_t lda2 );
void   magma_minproductblas_zswapblk(char storev, 
              magma_minproduct_int_t N, 
              cuDoubleComplex *dA1, magma_minproduct_int_t lda1, 
              cuDoubleComplex *dA2, magma_minproduct_int_t lda2,
              magma_minproduct_int_t i1, magma_minproduct_int_t i2, 
              magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci, 
              magma_minproduct_int_t offset);
void magma_minproductblas_zswapdblk(magma_minproduct_int_t n, magma_minproduct_int_t nb,
             cuDoubleComplex *dA1, magma_minproduct_int_t ldda1, magma_minproduct_int_t inca1,
             cuDoubleComplex *dA2, magma_minproduct_int_t ldda2, magma_minproduct_int_t inca2 );

  /*
   * Level 2 BLAS
   */
void magma_minproductblas_zgemv(char t, magma_minproduct_int_t M, magma_minproduct_int_t N, 
             cuDoubleComplex alpha,
             cuDoubleComplex *A, magma_minproduct_int_t lda, 
             cuDoubleComplex * X, magma_minproduct_int_t incX, 
             cuDoubleComplex beta, 
             cuDoubleComplex *Y, magma_minproduct_int_t incY);
#if defined(PRECISION_z) || defined(PRECISION_c)
magma_minproduct_int_t magma_minproductblas_zhemv(char u, magma_minproduct_int_t N, 
                            cuDoubleComplex alpha, 
                            cuDoubleComplex *A, magma_minproduct_int_t lda, 
                            cuDoubleComplex *X, magma_minproduct_int_t incX, 
                            cuDoubleComplex beta, 
                            cuDoubleComplex *Y, magma_minproduct_int_t incY);
#endif
magma_minproduct_int_t magma_minproductblas_zsymv(char u, magma_minproduct_int_t N, 
                            cuDoubleComplex alpha, 
                            cuDoubleComplex *A, magma_minproduct_int_t lda, 
                            cuDoubleComplex *X, magma_minproduct_int_t incX, 
                            cuDoubleComplex beta, 
                            cuDoubleComplex *Y, magma_minproduct_int_t incY);

  /*
   * Level 3 BLAS
   */
void magma_minproductblas_zgemm(char tA, char tB,
             magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k, 
             cuDoubleComplex alpha,
             const cuDoubleComplex *A, magma_minproduct_int_t lda, 
             const cuDoubleComplex *B, magma_minproduct_int_t ldb, 
             cuDoubleComplex beta,
             cuDoubleComplex *C, magma_minproduct_int_t ldc);
void magma_minproductblas_zgemm_fermi80(char tA, char tB, 
                 magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
                 cuDoubleComplex alpha, 
                 const cuDoubleComplex *A, magma_minproduct_int_t lda, 
                 const cuDoubleComplex *B, magma_minproduct_int_t ldb,
                 cuDoubleComplex beta, 
                 cuDoubleComplex *C, magma_minproduct_int_t ldc);
void magma_minproductblas_zgemm_fermi64(char tA, char tB, 
                 magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
                 cuDoubleComplex alpha, 
                 const cuDoubleComplex *A, magma_minproduct_int_t lda, 
                 const cuDoubleComplex *B, magma_minproduct_int_t ldb, 
                 cuDoubleComplex beta, 
                 cuDoubleComplex *C, magma_minproduct_int_t ldc);
void magma_minproductblas_zhemm(char s, char u,          
             magma_minproduct_int_t m, magma_minproduct_int_t n,
             cuDoubleComplex alpha, 
             const cuDoubleComplex *A, magma_minproduct_int_t lda,
             const cuDoubleComplex *B, magma_minproduct_int_t ldb,
             cuDoubleComplex beta, 
             cuDoubleComplex *C, magma_minproduct_int_t ldc);
void magma_minproductblas_zsymm(char s, char u,
             magma_minproduct_int_t m, magma_minproduct_int_t n,
             cuDoubleComplex alpha, 
             const cuDoubleComplex *A, magma_minproduct_int_t lda, 
             const cuDoubleComplex *B, magma_minproduct_int_t ldb,
             cuDoubleComplex beta,
             cuDoubleComplex *C, magma_minproduct_int_t ldc);
void magma_minproductblas_zsyrk(char u, char t,
             magma_minproduct_int_t n, magma_minproduct_int_t k, 
             cuDoubleComplex alpha, 
             const cuDoubleComplex *A, magma_minproduct_int_t lda,
             cuDoubleComplex beta,
             cuDoubleComplex *C, magma_minproduct_int_t ldc);
void magma_minproductblas_zherk(char u, char t,
             magma_minproduct_int_t n, magma_minproduct_int_t k, 
             double  alpha, 
             const cuDoubleComplex *A, magma_minproduct_int_t lda,
             double  beta, 
             cuDoubleComplex *C, magma_minproduct_int_t ldc);
void magma_minproductblas_zsyr2k(char u, char t,
              magma_minproduct_int_t n, magma_minproduct_int_t k,
              cuDoubleComplex alpha, 
              const cuDoubleComplex *A, magma_minproduct_int_t lda,
              const cuDoubleComplex *B, magma_minproduct_int_t ldb, 
              cuDoubleComplex beta, 
              cuDoubleComplex *C, magma_minproduct_int_t ldc);
void magma_minproductblas_zher2k(char u, char t,
              magma_minproduct_int_t n, magma_minproduct_int_t k, 
              cuDoubleComplex alpha, 
              const cuDoubleComplex *A, magma_minproduct_int_t lda, 
              const cuDoubleComplex *B, magma_minproduct_int_t ldb,
              double  beta,
              cuDoubleComplex *C, magma_minproduct_int_t ldc);
void magma_minproductblas_ztrmm(char s, char u, char t,  char d, 
             magma_minproduct_int_t m, magma_minproduct_int_t n,
             cuDoubleComplex alpha,
             const cuDoubleComplex *A, magma_minproduct_int_t lda,
             cuDoubleComplex *B, magma_minproduct_int_t ldb);
void magma_minproductblas_ztrsm(char s, char u, char t, char d,
             magma_minproduct_int_t m, magma_minproduct_int_t n,
             cuDoubleComplex alpha,
             /*const*/ cuDoubleComplex *A, magma_minproduct_int_t lda,
             cuDoubleComplex *B, magma_minproduct_int_t ldb);


  /*
   * Workspace interface (alphabetical order)
   */
magma_minproduct_int_t magma_minproductblasw_zsymv(char u, magma_minproduct_int_t N, 
                 cuDoubleComplex alpha, 
                 cuDoubleComplex *A, magma_minproduct_int_t lda, 
                 cuDoubleComplex *X, magma_minproduct_int_t incX, 
                 cuDoubleComplex beta, 
                 cuDoubleComplex *Y, magma_minproduct_int_t incY,
                 cuDoubleComplex *dWork);

#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif
