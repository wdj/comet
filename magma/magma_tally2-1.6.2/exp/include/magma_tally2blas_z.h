/*
 *   -- MAGMA_tally2 (version 1.6.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      @date January 2015
 *
 * @precisions normal z -> s d c
 */

#ifndef _MAGMA_tally2BLAS_Z_H_
#define _MAGMA_tally2BLAS_Z_H_

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
void zzero_32x32_block(           cuDoubleComplex *, magma_tally2_int_t);
void zzero_nbxnb_block(           magma_tally2_int_t, cuDoubleComplex *, magma_tally2_int_t);
void magma_tally2blas_zinplace_transpose(cuDoubleComplex *, magma_tally2_int_t, magma_tally2_int_t);
void magma_tally2blas_zpermute_long(     cuDoubleComplex *, magma_tally2_int_t, 
                  magma_tally2_int_t *, magma_tally2_int_t, magma_tally2_int_t);
void magma_tally2blas_zpermute_long2(    cuDoubleComplex *, magma_tally2_int_t, 
                  magma_tally2_int_t *, magma_tally2_int_t, magma_tally2_int_t);
void magma_tally2blas_ztranspose(        cuDoubleComplex *, magma_tally2_int_t, 
                  cuDoubleComplex *, magma_tally2_int_t, 
                  magma_tally2_int_t, magma_tally2_int_t);
void magma_tally2blas_ztranspose2(       cuDoubleComplex *, magma_tally2_int_t, 
                  cuDoubleComplex *, magma_tally2_int_t, 
                  magma_tally2_int_t, magma_tally2_int_t);
void magma_tally2blas_zhtodt(            cuDoubleComplex  *ha, int lda,
                          cuDoubleComplex *dat, int ldda,
                          cuDoubleComplex  *dB, int lddb,
                                  int m, int n , int nb);
void magma_tally2blas_zdtoht(            cuDoubleComplex *dat, int ldda,
                  cuDoubleComplex  *ha, int lda,
                       cuDoubleComplex  *dB, int lddb,
                  int m, int n , int nb);
  
  /*
   * LAPACK auxiliary functions
   */
void   magma_tally2blas_zlacpy( char uplo, 
             magma_tally2_int_t m, magma_tally2_int_t n, 
             cuDoubleComplex *A, magma_tally2_int_t lda, 
             cuDoubleComplex *B, magma_tally2_int_t ldb);
double magma_tally2blas_zlange( char norm, 
             magma_tally2_int_t m, magma_tally2_int_t n, 
             cuDoubleComplex *A, magma_tally2_int_t lda, double *WORK);
double magma_tally2blas_zlanhe( char norm, char uplo, 
             magma_tally2_int_t n,
             cuDoubleComplex *A, magma_tally2_int_t lda, double *WORK);
double magma_tally2blas_zlansy( char norm, char uplo,
             magma_tally2_int_t n, 
             cuDoubleComplex *A, magma_tally2_int_t lda, double *WORK);
void   magma_tally2blas_zlascl( char type, int kl, int ku,
             double cfrom, double cto,
             int m, int n,
             cuDoubleComplex *A, int lda, int *info );
void   magma_tally2blas_zlaset( magma_tally2_int_t m, magma_tally2_int_t n,
             cuDoubleComplex *A, magma_tally2_int_t lda);
void   magma_tally2blas_zlaswp( magma_tally2_int_t N, 
             cuDoubleComplex *dAT, magma_tally2_int_t lda, 
             magma_tally2_int_t i1,  magma_tally2_int_t i2, 
             magma_tally2_int_t *ipiv, magma_tally2_int_t inci );
void   magma_tally2blas_zlaswpx(magma_tally2_int_t N, 
             cuDoubleComplex *dAT, magma_tally2_int_t ldx, magma_tally2_int_t ldy, 
             magma_tally2_int_t i1, magma_tally2_int_t i2,
             magma_tally2_int_t *ipiv, magma_tally2_int_t inci );

  /*
   * Level 1 BLAS
   */
void   magma_tally2blas_zswap(   magma_tally2_int_t N, 
              cuDoubleComplex *dA1, magma_tally2_int_t lda1, 
              cuDoubleComplex *dA2, magma_tally2_int_t lda2 );
void   magma_tally2blas_zswapblk(char storev, 
              magma_tally2_int_t N, 
              cuDoubleComplex *dA1, magma_tally2_int_t lda1, 
              cuDoubleComplex *dA2, magma_tally2_int_t lda2,
              magma_tally2_int_t i1, magma_tally2_int_t i2, 
              magma_tally2_int_t *ipiv, magma_tally2_int_t inci, 
              magma_tally2_int_t offset);
void magma_tally2blas_zswapdblk(magma_tally2_int_t n, magma_tally2_int_t nb,
             cuDoubleComplex *dA1, magma_tally2_int_t ldda1, magma_tally2_int_t inca1,
             cuDoubleComplex *dA2, magma_tally2_int_t ldda2, magma_tally2_int_t inca2 );

  /*
   * Level 2 BLAS
   */
void magma_tally2blas_zgemv(char t, magma_tally2_int_t M, magma_tally2_int_t N, 
             cuDoubleComplex alpha,
             cuDoubleComplex *A, magma_tally2_int_t lda, 
             cuDoubleComplex * X, magma_tally2_int_t incX, 
             cuDoubleComplex beta, 
             cuDoubleComplex *Y, magma_tally2_int_t incY);
#if defined(PRECISION_z) || defined(PRECISION_c)
magma_tally2_int_t magma_tally2blas_zhemv(char u, magma_tally2_int_t N, 
                            cuDoubleComplex alpha, 
                            cuDoubleComplex *A, magma_tally2_int_t lda, 
                            cuDoubleComplex *X, magma_tally2_int_t incX, 
                            cuDoubleComplex beta, 
                            cuDoubleComplex *Y, magma_tally2_int_t incY);
#endif
magma_tally2_int_t magma_tally2blas_zsymv(char u, magma_tally2_int_t N, 
                            cuDoubleComplex alpha, 
                            cuDoubleComplex *A, magma_tally2_int_t lda, 
                            cuDoubleComplex *X, magma_tally2_int_t incX, 
                            cuDoubleComplex beta, 
                            cuDoubleComplex *Y, magma_tally2_int_t incY);

  /*
   * Level 3 BLAS
   */
void magma_tally2blas_zgemm(char tA, char tB,
             magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k, 
             cuDoubleComplex alpha,
             const cuDoubleComplex *A, magma_tally2_int_t lda, 
             const cuDoubleComplex *B, magma_tally2_int_t ldb, 
             cuDoubleComplex beta,
             cuDoubleComplex *C, magma_tally2_int_t ldc);
void magma_tally2blas_zgemm_fermi80(char tA, char tB, 
                 magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
                 cuDoubleComplex alpha, 
                 const cuDoubleComplex *A, magma_tally2_int_t lda, 
                 const cuDoubleComplex *B, magma_tally2_int_t ldb,
                 cuDoubleComplex beta, 
                 cuDoubleComplex *C, magma_tally2_int_t ldc);
void magma_tally2blas_zgemm_fermi64(char tA, char tB, 
                 magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
                 cuDoubleComplex alpha, 
                 const cuDoubleComplex *A, magma_tally2_int_t lda, 
                 const cuDoubleComplex *B, magma_tally2_int_t ldb, 
                 cuDoubleComplex beta, 
                 cuDoubleComplex *C, magma_tally2_int_t ldc);
void magma_tally2blas_zhemm(char s, char u,          
             magma_tally2_int_t m, magma_tally2_int_t n,
             cuDoubleComplex alpha, 
             const cuDoubleComplex *A, magma_tally2_int_t lda,
             const cuDoubleComplex *B, magma_tally2_int_t ldb,
             cuDoubleComplex beta, 
             cuDoubleComplex *C, magma_tally2_int_t ldc);
void magma_tally2blas_zsymm(char s, char u,
             magma_tally2_int_t m, magma_tally2_int_t n,
             cuDoubleComplex alpha, 
             const cuDoubleComplex *A, magma_tally2_int_t lda, 
             const cuDoubleComplex *B, magma_tally2_int_t ldb,
             cuDoubleComplex beta,
             cuDoubleComplex *C, magma_tally2_int_t ldc);
void magma_tally2blas_zsyrk(char u, char t,
             magma_tally2_int_t n, magma_tally2_int_t k, 
             cuDoubleComplex alpha, 
             const cuDoubleComplex *A, magma_tally2_int_t lda,
             cuDoubleComplex beta,
             cuDoubleComplex *C, magma_tally2_int_t ldc);
void magma_tally2blas_zherk(char u, char t,
             magma_tally2_int_t n, magma_tally2_int_t k, 
             double  alpha, 
             const cuDoubleComplex *A, magma_tally2_int_t lda,
             double  beta, 
             cuDoubleComplex *C, magma_tally2_int_t ldc);
void magma_tally2blas_zsyr2k(char u, char t,
              magma_tally2_int_t n, magma_tally2_int_t k,
              cuDoubleComplex alpha, 
              const cuDoubleComplex *A, magma_tally2_int_t lda,
              const cuDoubleComplex *B, magma_tally2_int_t ldb, 
              cuDoubleComplex beta, 
              cuDoubleComplex *C, magma_tally2_int_t ldc);
void magma_tally2blas_zher2k(char u, char t,
              magma_tally2_int_t n, magma_tally2_int_t k, 
              cuDoubleComplex alpha, 
              const cuDoubleComplex *A, magma_tally2_int_t lda, 
              const cuDoubleComplex *B, magma_tally2_int_t ldb,
              double  beta,
              cuDoubleComplex *C, magma_tally2_int_t ldc);
void magma_tally2blas_ztrmm(char s, char u, char t,  char d, 
             magma_tally2_int_t m, magma_tally2_int_t n,
             cuDoubleComplex alpha,
             const cuDoubleComplex *A, magma_tally2_int_t lda,
             cuDoubleComplex *B, magma_tally2_int_t ldb);
void magma_tally2blas_ztrsm(char s, char u, char t, char d,
             magma_tally2_int_t m, magma_tally2_int_t n,
             cuDoubleComplex alpha,
             /*const*/ cuDoubleComplex *A, magma_tally2_int_t lda,
             cuDoubleComplex *B, magma_tally2_int_t ldb);


  /*
   * Workspace interface (alphabetical order)
   */
magma_tally2_int_t magma_tally2blasw_zsymv(char u, magma_tally2_int_t N, 
                 cuDoubleComplex alpha, 
                 cuDoubleComplex *A, magma_tally2_int_t lda, 
                 cuDoubleComplex *X, magma_tally2_int_t incX, 
                 cuDoubleComplex beta, 
                 cuDoubleComplex *Y, magma_tally2_int_t incY,
                 cuDoubleComplex *dWork);

#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif
