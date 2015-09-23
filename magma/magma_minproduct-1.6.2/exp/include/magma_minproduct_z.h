/*
 *   -- MAGMA_minproduct (version 1.6.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      @date January 2015
 *
 * @precisions normal z -> s d c
 */

#ifndef _MAGMA_minproduct_Z_H_
#define _MAGMA_minproduct_Z_H_
#define PRECISION_z

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct function definitions / Data on CPU
*/
magma_minproduct_int_t magma_minproduct_zgebrd( magma_minproduct_int_t m, magma_minproduct_int_t n, cuDoubleComplex *A, 
              magma_minproduct_int_t lda, double *d, double *e,
              cuDoubleComplex *tauq,  cuDoubleComplex *taup, 
              cuDoubleComplex *work, magma_minproduct_int_t lwork, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgehrd2(magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
              cuDoubleComplex *A, magma_minproduct_int_t lda, cuDoubleComplex *tau, 
              cuDoubleComplex *work, magma_minproduct_int_t *lwork, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgehrd( magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
              cuDoubleComplex *A, magma_minproduct_int_t lda, cuDoubleComplex *tau,
              cuDoubleComplex *work, magma_minproduct_int_t lwork,
              cuDoubleComplex *d_T, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgelqf( magma_minproduct_int_t m, magma_minproduct_int_t n, 
                          cuDoubleComplex *A,    magma_minproduct_int_t lda,   cuDoubleComplex *tau, 
                          cuDoubleComplex *work, magma_minproduct_int_t lwork, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgeqlf( magma_minproduct_int_t m, magma_minproduct_int_t n, 
                          cuDoubleComplex *A,    magma_minproduct_int_t lda,   cuDoubleComplex *tau, 
                          cuDoubleComplex *work, magma_minproduct_int_t lwork, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgeqrf( magma_minproduct_int_t m, magma_minproduct_int_t n, cuDoubleComplex *A, 
              magma_minproduct_int_t lda, cuDoubleComplex *tau, cuDoubleComplex *work, 
              magma_minproduct_int_t lwork, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgeqrf_ooc( magma_minproduct_int_t m, magma_minproduct_int_t n, cuDoubleComplex *A,
              magma_minproduct_int_t lda, cuDoubleComplex *tau, cuDoubleComplex *work,
                  magma_minproduct_int_t lwork, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgesv ( magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
              cuDoubleComplex *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv, 
              cuDoubleComplex *B, magma_minproduct_int_t ldb, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgetrf( magma_minproduct_int_t m, magma_minproduct_int_t n, cuDoubleComplex *A, 
              magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv, 
              magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgetrf_mc( magma_minproduct_context *cntxt, magma_minproduct_int_t *m, magma_minproduct_int_t *n, cuDoubleComplex *A,
                          magma_minproduct_int_t *lda, magma_minproduct_int_t *ipiv, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgeqrf_mc(magma_minproduct_context *cntxt, magma_minproduct_int_t *m, magma_minproduct_int_t *n,
                          cuDoubleComplex *A, magma_minproduct_int_t *lda,
                          cuDoubleComplex *tau, cuDoubleComplex *work,
                          magma_minproduct_int_t *lwork, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgetrf2(magma_minproduct_int_t m, magma_minproduct_int_t n, cuDoubleComplex *a, 
                          magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zlatrd( char uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb, cuDoubleComplex *a, 
                          magma_minproduct_int_t lda, double *e, cuDoubleComplex *tau, 
              cuDoubleComplex *w, magma_minproduct_int_t ldw,
                          cuDoubleComplex *da, magma_minproduct_int_t ldda, 
              cuDoubleComplex *dw, magma_minproduct_int_t lddw);
magma_minproduct_int_t magma_minproduct_zlatrd2(char uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
              cuDoubleComplex *a,  magma_minproduct_int_t lda,
              double *e, cuDoubleComplex *tau,
              cuDoubleComplex *w,  magma_minproduct_int_t ldw,
              cuDoubleComplex *da, magma_minproduct_int_t ldda,
              cuDoubleComplex *dw, magma_minproduct_int_t lddw,
              cuDoubleComplex *dwork, magma_minproduct_int_t ldwork);
magma_minproduct_int_t magma_minproduct_zlahr2( magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb, 
              cuDoubleComplex *da, cuDoubleComplex *dv, cuDoubleComplex *a, 
              magma_minproduct_int_t lda, cuDoubleComplex *tau, cuDoubleComplex *t, 
              magma_minproduct_int_t ldt, cuDoubleComplex *y, magma_minproduct_int_t ldy);
magma_minproduct_int_t magma_minproduct_zlahru( magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb, 
              cuDoubleComplex *a, magma_minproduct_int_t lda, 
              cuDoubleComplex *da, cuDoubleComplex *y, 
              cuDoubleComplex *v, cuDoubleComplex *t, 
              cuDoubleComplex *dwork);
magma_minproduct_int_t magma_minproduct_zposv ( char uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
              cuDoubleComplex *A, magma_minproduct_int_t lda, 
              cuDoubleComplex *B, magma_minproduct_int_t ldb, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zpotrf( char uplo, magma_minproduct_int_t n, cuDoubleComplex *A, 
              magma_minproduct_int_t lda, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zpotrf_mc( magma_minproduct_context *cntxt, char *uplo, magma_minproduct_int_t *n, cuDoubleComplex *A,
                          magma_minproduct_int_t *lda, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zpotri( char uplo, magma_minproduct_int_t n, cuDoubleComplex *A,
                  magma_minproduct_int_t lda, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zlauum( char uplo, magma_minproduct_int_t n, cuDoubleComplex *A,
                  magma_minproduct_int_t lda, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_ztrtri( char uplo, char diag, magma_minproduct_int_t n, cuDoubleComplex *A, 
                  magma_minproduct_int_t lda, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zhetrd( char uplo, magma_minproduct_int_t n, cuDoubleComplex *A, 
              magma_minproduct_int_t lda, double *d, double *e, 
              cuDoubleComplex *tau, cuDoubleComplex *work, magma_minproduct_int_t lwork, 
              magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgeqrf2(magma_minproduct_context *cntxt, magma_minproduct_int_t m, magma_minproduct_int_t n,
                          cuDoubleComplex *a, magma_minproduct_int_t lda, cuDoubleComplex *tau,
                          cuDoubleComplex *work, magma_minproduct_int_t lwork,
                          magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgeqrf3(magma_minproduct_context *cntxt, magma_minproduct_int_t m, magma_minproduct_int_t n,
                          cuDoubleComplex *a, magma_minproduct_int_t lda, cuDoubleComplex *tau,
                          cuDoubleComplex *work, magma_minproduct_int_t lwork,
                          magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zungqr( magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
              cuDoubleComplex *a, magma_minproduct_int_t lda,
              cuDoubleComplex *tau, cuDoubleComplex *dwork,
              magma_minproduct_int_t nb, magma_minproduct_int_t *info );
magma_minproduct_int_t magma_minproduct_zunmql( const char side, const char trans,
              magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
              cuDoubleComplex *a, magma_minproduct_int_t lda,
              cuDoubleComplex *tau,
              cuDoubleComplex *c, magma_minproduct_int_t ldc,
              cuDoubleComplex *work, magma_minproduct_int_t lwork,
              magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zunmqr( char side, char trans, 
                          magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k, 
                          cuDoubleComplex *a, magma_minproduct_int_t lda, cuDoubleComplex *tau, 
                          cuDoubleComplex *c, magma_minproduct_int_t ldc, 
                          cuDoubleComplex *work, magma_minproduct_int_t lwork, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zunmtr( char side, char uplo, char trans,
              magma_minproduct_int_t m, magma_minproduct_int_t n,
              cuDoubleComplex *a,    magma_minproduct_int_t lda,
              cuDoubleComplex *tau,
              cuDoubleComplex *c,    magma_minproduct_int_t ldc,
              cuDoubleComplex *work, magma_minproduct_int_t lwork,
              magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zunghr( magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
              cuDoubleComplex *a, magma_minproduct_int_t lda,
              cuDoubleComplex *tau,
              cuDoubleComplex *dT, magma_minproduct_int_t nb,
              magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zheev( char jobz, char uplo, magma_minproduct_int_t n,
             cuDoubleComplex *a, magma_minproduct_int_t lda, double *w,
             cuDoubleComplex *work, magma_minproduct_int_t lwork,
                 double *rwork, magma_minproduct_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
magma_minproduct_int_t  magma_minproduct_zgeev( char jobvl, char jobvr, magma_minproduct_int_t n,
              cuDoubleComplex *a, magma_minproduct_int_t lda,
              cuDoubleComplex *w,
              cuDoubleComplex *vl, magma_minproduct_int_t ldvl,
              cuDoubleComplex *vr, magma_minproduct_int_t ldvr,
              cuDoubleComplex *work, magma_minproduct_int_t lwork,
              double *rwork, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgesvd( char jobu, char jobvt, magma_minproduct_int_t m, magma_minproduct_int_t n,
              cuDoubleComplex *a,    magma_minproduct_int_t lda, double *s, 
                          cuDoubleComplex *u,    magma_minproduct_int_t ldu, 
              cuDoubleComplex *vt,   magma_minproduct_int_t ldvt,
              cuDoubleComplex *work, magma_minproduct_int_t lwork,
              double *rwork, magma_minproduct_int_t *info );
magma_minproduct_int_t magma_minproduct_zheevd( char jobz, char uplo, magma_minproduct_int_t n,
              cuDoubleComplex *a, magma_minproduct_int_t lda, double *w,
              cuDoubleComplex *work, magma_minproduct_int_t lwork,
                  double *rwork, magma_minproduct_int_t lrwork,
                  magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zhegvd( magma_minproduct_int_t itype, char jobz, char uplo, magma_minproduct_int_t n,
              cuDoubleComplex *a, magma_minproduct_int_t lda,
              cuDoubleComplex *b, magma_minproduct_int_t ldb,
              double *w, cuDoubleComplex *work, magma_minproduct_int_t lwork,
              double *rwork, magma_minproduct_int_t lrwork, magma_minproduct_int_t *iwork,
              magma_minproduct_int_t liwork, magma_minproduct_int_t *info);
#else
magma_minproduct_int_t  magma_minproduct_zgeev( char jobvl, char jobvr, magma_minproduct_int_t n,
              cuDoubleComplex *a,    magma_minproduct_int_t lda,
              cuDoubleComplex *wr, cuDoubleComplex *wi,
              cuDoubleComplex *vl,   magma_minproduct_int_t ldvl,
              cuDoubleComplex *vr,   magma_minproduct_int_t ldvr,
              cuDoubleComplex *work, magma_minproduct_int_t lwork,
              magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgesvd( char jobu, char jobvt, magma_minproduct_int_t m, magma_minproduct_int_t n,
              cuDoubleComplex *a,    magma_minproduct_int_t lda, double *s, 
                          cuDoubleComplex *u,    magma_minproduct_int_t ldu, 
              cuDoubleComplex *vt,   magma_minproduct_int_t ldvt,
              cuDoubleComplex *work, magma_minproduct_int_t lwork,
              magma_minproduct_int_t *info );
magma_minproduct_int_t magma_minproduct_zheevd( char jobz, char uplo, magma_minproduct_int_t n,
              cuDoubleComplex *a, magma_minproduct_int_t lda, double *w,
              cuDoubleComplex *work, magma_minproduct_int_t lwork,
              magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zhegvd( magma_minproduct_int_t itype, char jobz, char uplo, magma_minproduct_int_t n,
              cuDoubleComplex *a, magma_minproduct_int_t lda,
              cuDoubleComplex *b, magma_minproduct_int_t ldb,
              double *w, cuDoubleComplex *work, magma_minproduct_int_t lwork,
              magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork, magma_minproduct_int_t *info);
#endif

magma_minproduct_int_t magma_minproduct_zhegst( magma_minproduct_int_t itype, char uplo, magma_minproduct_int_t n,
              cuDoubleComplex *a, magma_minproduct_int_t lda,
              cuDoubleComplex *b, magma_minproduct_int_t ldb, magma_minproduct_int_t *info);

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA_minproduct function definitions / Data on GPU
*/
magma_minproduct_int_t magma_minproduct_zgels_gpu(  char trans, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  cuDoubleComplex *dA,    magma_minproduct_int_t ldda, 
                  cuDoubleComplex *dB,    magma_minproduct_int_t lddb, 
                  cuDoubleComplex *hwork, magma_minproduct_int_t lwork, 
                  magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgels3_gpu( char trans, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  cuDoubleComplex *dA,    magma_minproduct_int_t ldda,
                  cuDoubleComplex *dB,    magma_minproduct_int_t lddb,
                  cuDoubleComplex *hwork, magma_minproduct_int_t lwork,
                  magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgelqf_gpu( magma_minproduct_int_t m, magma_minproduct_int_t n,
                  cuDoubleComplex *dA,    magma_minproduct_int_t ldda,   cuDoubleComplex *tau,
                  cuDoubleComplex *work, magma_minproduct_int_t lwork, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgeqrf_gpu( magma_minproduct_int_t m, magma_minproduct_int_t n, 
                  cuDoubleComplex *dA,  magma_minproduct_int_t ldda, 
                  cuDoubleComplex *tau, cuDoubleComplex *dT, 
                  magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgeqrf2_gpu(magma_minproduct_int_t m, magma_minproduct_int_t n, 
                  cuDoubleComplex *dA,  magma_minproduct_int_t ldda, 
                  cuDoubleComplex *tau, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgeqrf3_gpu(magma_minproduct_int_t m, magma_minproduct_int_t n, 
                  cuDoubleComplex *dA,  magma_minproduct_int_t ldda, 
                  cuDoubleComplex *tau, cuDoubleComplex *dT, 
                  magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgeqrs_gpu( magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
                  cuDoubleComplex *dA,     magma_minproduct_int_t ldda, 
                  cuDoubleComplex *tau,   cuDoubleComplex *dT,
                  cuDoubleComplex *dB,    magma_minproduct_int_t lddb,
                  cuDoubleComplex *hwork, magma_minproduct_int_t lhwork, 
                  magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgeqrs3_gpu( magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
                  cuDoubleComplex *dA,     magma_minproduct_int_t ldda, 
                  cuDoubleComplex *tau,   cuDoubleComplex *dT,
                  cuDoubleComplex *dB,    magma_minproduct_int_t lddb,
                  cuDoubleComplex *hwork, magma_minproduct_int_t lhwork, 
                  magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgessm_gpu( char storev, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t ib, 
                              magma_minproduct_int_t *ipiv, 
                              cuDoubleComplex *dL1, magma_minproduct_int_t lddl1, 
                              cuDoubleComplex *dL,  magma_minproduct_int_t lddl, 
                              cuDoubleComplex *dA,  magma_minproduct_int_t ldda, 
                              magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgesv_gpu(  magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
                  cuDoubleComplex *dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv, 
                  cuDoubleComplex *dB, magma_minproduct_int_t lddb, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgetrf_incpiv_gpu( char storev, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t ib,
                              cuDoubleComplex *hA, magma_minproduct_int_t ldha, cuDoubleComplex *dA, magma_minproduct_int_t ldda,
                              cuDoubleComplex *hL, magma_minproduct_int_t ldhl, cuDoubleComplex *dL, magma_minproduct_int_t lddl,
                              magma_minproduct_int_t *ipiv, 
                              cuDoubleComplex *dwork, magma_minproduct_int_t lddwork,
                              magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgetrf_gpu( magma_minproduct_int_t m, magma_minproduct_int_t n, 
                  cuDoubleComplex *dA, magma_minproduct_int_t ldda, 
                  magma_minproduct_int_t *ipiv, magma_minproduct_int_t *info);
magma_minproduct_int_t 
magma_minproduct_zgetrf_nopiv_gpu      ( magma_minproduct_int_t m, magma_minproduct_int_t n,
                  cuDoubleComplex *dA, magma_minproduct_int_t ldda,
                  magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zgetrs_gpu( char trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
                  cuDoubleComplex *dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv, 
                  cuDoubleComplex *dB, magma_minproduct_int_t lddb, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zlabrd_gpu( magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb, 
                              cuDoubleComplex *a, magma_minproduct_int_t lda, cuDoubleComplex *da, magma_minproduct_int_t ldda,
                              double *d, double *e, cuDoubleComplex *tauq, cuDoubleComplex *taup,  
                              cuDoubleComplex *x, magma_minproduct_int_t ldx, cuDoubleComplex *dx, magma_minproduct_int_t lddx, 
                              cuDoubleComplex *y, magma_minproduct_int_t ldy, cuDoubleComplex *dy, magma_minproduct_int_t lddy);
magma_minproduct_int_t magma_minproduct_zlarfb_gpu( char side, char trans, char direct, char storev, 
                  magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
                  cuDoubleComplex *dv, magma_minproduct_int_t ldv, cuDoubleComplex *dt,    magma_minproduct_int_t ldt, 
                  cuDoubleComplex *dc, magma_minproduct_int_t ldc, cuDoubleComplex *dowrk, magma_minproduct_int_t ldwork );
magma_minproduct_int_t magma_minproduct_zposv_gpu(  char uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
                  cuDoubleComplex *dA, magma_minproduct_int_t ldda, 
                  cuDoubleComplex *dB, magma_minproduct_int_t lddb, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zpotrf_gpu( char uplo,  magma_minproduct_int_t n, 
                  cuDoubleComplex *dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zpotri_gpu( char uplo,  magma_minproduct_int_t n,
                      cuDoubleComplex *dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zlauum_gpu( char uplo,  magma_minproduct_int_t n,
                      cuDoubleComplex *dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_ztrtri_gpu( char uplo,  char diag, magma_minproduct_int_t n,
                      cuDoubleComplex *dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zhetrd_gpu( char uplo, magma_minproduct_int_t n,
                  cuDoubleComplex *da, magma_minproduct_int_t ldda,
                  double *d, double *e, cuDoubleComplex *tau,
                  cuDoubleComplex *wa,  magma_minproduct_int_t ldwa,
                  cuDoubleComplex *work, magma_minproduct_int_t lwork,
                  magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zhetrd2_gpu(char uplo, magma_minproduct_int_t n,
                  cuDoubleComplex *da, magma_minproduct_int_t ldda,
                  double *d, double *e, cuDoubleComplex *tau,
                  cuDoubleComplex *wa,  magma_minproduct_int_t ldwa,
                  cuDoubleComplex *work, magma_minproduct_int_t lwork,
                  cuDoubleComplex *dwork, magma_minproduct_int_t ldwork,
                  magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zpotrs_gpu( char uplo,  magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
                  cuDoubleComplex *dA, magma_minproduct_int_t ldda, 
                  cuDoubleComplex *dB, magma_minproduct_int_t lddb, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zssssm_gpu( char storev, magma_minproduct_int_t m1, magma_minproduct_int_t n1, 
                              magma_minproduct_int_t m2, magma_minproduct_int_t n2, magma_minproduct_int_t k, magma_minproduct_int_t ib, 
                              cuDoubleComplex *dA1, magma_minproduct_int_t ldda1, 
                              cuDoubleComplex *dA2, magma_minproduct_int_t ldda2, 
                              cuDoubleComplex *dL1, magma_minproduct_int_t lddl1, 
                              cuDoubleComplex *dL2, magma_minproduct_int_t lddl2,
                              magma_minproduct_int_t *IPIV, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_ztstrf_gpu( char storev, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t ib, magma_minproduct_int_t nb,
                              cuDoubleComplex *hU, magma_minproduct_int_t ldhu, cuDoubleComplex *dU, magma_minproduct_int_t lddu, 
                              cuDoubleComplex *hA, magma_minproduct_int_t ldha, cuDoubleComplex *dA, magma_minproduct_int_t ldda, 
                              cuDoubleComplex *hL, magma_minproduct_int_t ldhl, cuDoubleComplex *dL, magma_minproduct_int_t lddl,
                              magma_minproduct_int_t *ipiv, 
                              cuDoubleComplex *hwork, magma_minproduct_int_t ldhwork, 
                  cuDoubleComplex *dwork, magma_minproduct_int_t lddwork,
                              magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zungqr_gpu( magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k, 
                              cuDoubleComplex *da, magma_minproduct_int_t ldda, 
                              cuDoubleComplex *tau, cuDoubleComplex *dwork, 
                              magma_minproduct_int_t nb, magma_minproduct_int_t *info );
magma_minproduct_int_t magma_minproduct_zunmql2_gpu(const char side, const char trans,
                  magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
                  cuDoubleComplex *da, magma_minproduct_int_t ldda,
                  cuDoubleComplex *tau,
                  cuDoubleComplex *dc, magma_minproduct_int_t lddc,
                  cuDoubleComplex *wa, magma_minproduct_int_t ldwa,
                  magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zunmqr_gpu( char side, char trans, 
                              magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
                              cuDoubleComplex *a,    magma_minproduct_int_t lda, cuDoubleComplex *tau, 
                              cuDoubleComplex *c,    magma_minproduct_int_t ldc,
                              cuDoubleComplex *work, magma_minproduct_int_t lwork, 
                              cuDoubleComplex *td,   magma_minproduct_int_t nb, magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zunmqr2_gpu(const char side, const char trans,
                  magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
                  cuDoubleComplex *da,   magma_minproduct_int_t ldda,
                  cuDoubleComplex *tau,
                  cuDoubleComplex *dc,    magma_minproduct_int_t lddc,
                  cuDoubleComplex *wa,    magma_minproduct_int_t ldwa,
                  magma_minproduct_int_t *info);
magma_minproduct_int_t magma_minproduct_zunmtr_gpu( char side, char uplo, char trans,
                  magma_minproduct_int_t m, magma_minproduct_int_t n,
                  cuDoubleComplex *da,    magma_minproduct_int_t ldda,
                  cuDoubleComplex *tau,
                  cuDoubleComplex *dc,    magma_minproduct_int_t lddc,
                  cuDoubleComplex *wa,    magma_minproduct_int_t ldwa,
                  magma_minproduct_int_t *info);

#if defined(PRECISION_z) || defined(PRECISION_c)
magma_minproduct_int_t magma_minproduct_zheevd_gpu( char jobz, char uplo,
                  magma_minproduct_int_t n,
                  cuDoubleComplex *da, magma_minproduct_int_t ldda,
                  double *w,
                  cuDoubleComplex *wa,  magma_minproduct_int_t ldwa,
                  cuDoubleComplex *work, magma_minproduct_int_t lwork,
                  double *rwork, magma_minproduct_int_t lrwork,
                  magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
                  magma_minproduct_int_t *info);
#else
magma_minproduct_int_t magma_minproduct_zheevd_gpu( char jobz, char uplo,
                  magma_minproduct_int_t n,
                  cuDoubleComplex *da, magma_minproduct_int_t ldda,
                  cuDoubleComplex *w,
                  cuDoubleComplex *wa,  magma_minproduct_int_t ldwa,
                  cuDoubleComplex *work, magma_minproduct_int_t lwork,
                  magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
                  magma_minproduct_int_t *info);
#endif

magma_minproduct_int_t magma_minproduct_zhegst_gpu(magma_minproduct_int_t itype, char uplo, magma_minproduct_int_t n,
                 cuDoubleComplex *da, magma_minproduct_int_t ldda,
                 cuDoubleComplex *db, magma_minproduct_int_t lddb, magma_minproduct_int_t *info);


#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* _MAGMA_minproduct_Z_H_ */

