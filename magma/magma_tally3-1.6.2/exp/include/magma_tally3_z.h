/*
 *   -- MAGMA_tally3 (version 1.6.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      @date January 2015
 *
 * @precisions normal z -> s d c
 */

#ifndef _MAGMA_tally3_Z_H_
#define _MAGMA_tally3_Z_H_
#define PRECISION_z

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 function definitions / Data on CPU
*/
magma_tally3_int_t magma_tally3_zgebrd( magma_tally3_int_t m, magma_tally3_int_t n, cuDoubleComplex *A, 
              magma_tally3_int_t lda, double *d, double *e,
              cuDoubleComplex *tauq,  cuDoubleComplex *taup, 
              cuDoubleComplex *work, magma_tally3_int_t lwork, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgehrd2(magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
              cuDoubleComplex *A, magma_tally3_int_t lda, cuDoubleComplex *tau, 
              cuDoubleComplex *work, magma_tally3_int_t *lwork, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgehrd( magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
              cuDoubleComplex *A, magma_tally3_int_t lda, cuDoubleComplex *tau,
              cuDoubleComplex *work, magma_tally3_int_t lwork,
              cuDoubleComplex *d_T, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgelqf( magma_tally3_int_t m, magma_tally3_int_t n, 
                          cuDoubleComplex *A,    magma_tally3_int_t lda,   cuDoubleComplex *tau, 
                          cuDoubleComplex *work, magma_tally3_int_t lwork, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgeqlf( magma_tally3_int_t m, magma_tally3_int_t n, 
                          cuDoubleComplex *A,    magma_tally3_int_t lda,   cuDoubleComplex *tau, 
                          cuDoubleComplex *work, magma_tally3_int_t lwork, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgeqrf( magma_tally3_int_t m, magma_tally3_int_t n, cuDoubleComplex *A, 
              magma_tally3_int_t lda, cuDoubleComplex *tau, cuDoubleComplex *work, 
              magma_tally3_int_t lwork, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgeqrf_ooc( magma_tally3_int_t m, magma_tally3_int_t n, cuDoubleComplex *A,
              magma_tally3_int_t lda, cuDoubleComplex *tau, cuDoubleComplex *work,
                  magma_tally3_int_t lwork, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgesv ( magma_tally3_int_t n, magma_tally3_int_t nrhs, 
              cuDoubleComplex *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv, 
              cuDoubleComplex *B, magma_tally3_int_t ldb, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgetrf( magma_tally3_int_t m, magma_tally3_int_t n, cuDoubleComplex *A, 
              magma_tally3_int_t lda, magma_tally3_int_t *ipiv, 
              magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgetrf_mc( magma_tally3_context *cntxt, magma_tally3_int_t *m, magma_tally3_int_t *n, cuDoubleComplex *A,
                          magma_tally3_int_t *lda, magma_tally3_int_t *ipiv, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgeqrf_mc(magma_tally3_context *cntxt, magma_tally3_int_t *m, magma_tally3_int_t *n,
                          cuDoubleComplex *A, magma_tally3_int_t *lda,
                          cuDoubleComplex *tau, cuDoubleComplex *work,
                          magma_tally3_int_t *lwork, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgetrf2(magma_tally3_int_t m, magma_tally3_int_t n, cuDoubleComplex *a, 
                          magma_tally3_int_t lda, magma_tally3_int_t *ipiv, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zlatrd( char uplo, magma_tally3_int_t n, magma_tally3_int_t nb, cuDoubleComplex *a, 
                          magma_tally3_int_t lda, double *e, cuDoubleComplex *tau, 
              cuDoubleComplex *w, magma_tally3_int_t ldw,
                          cuDoubleComplex *da, magma_tally3_int_t ldda, 
              cuDoubleComplex *dw, magma_tally3_int_t lddw);
magma_tally3_int_t magma_tally3_zlatrd2(char uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
              cuDoubleComplex *a,  magma_tally3_int_t lda,
              double *e, cuDoubleComplex *tau,
              cuDoubleComplex *w,  magma_tally3_int_t ldw,
              cuDoubleComplex *da, magma_tally3_int_t ldda,
              cuDoubleComplex *dw, magma_tally3_int_t lddw,
              cuDoubleComplex *dwork, magma_tally3_int_t ldwork);
magma_tally3_int_t magma_tally3_zlahr2( magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb, 
              cuDoubleComplex *da, cuDoubleComplex *dv, cuDoubleComplex *a, 
              magma_tally3_int_t lda, cuDoubleComplex *tau, cuDoubleComplex *t, 
              magma_tally3_int_t ldt, cuDoubleComplex *y, magma_tally3_int_t ldy);
magma_tally3_int_t magma_tally3_zlahru( magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb, 
              cuDoubleComplex *a, magma_tally3_int_t lda, 
              cuDoubleComplex *da, cuDoubleComplex *y, 
              cuDoubleComplex *v, cuDoubleComplex *t, 
              cuDoubleComplex *dwork);
magma_tally3_int_t magma_tally3_zposv ( char uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs, 
              cuDoubleComplex *A, magma_tally3_int_t lda, 
              cuDoubleComplex *B, magma_tally3_int_t ldb, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zpotrf( char uplo, magma_tally3_int_t n, cuDoubleComplex *A, 
              magma_tally3_int_t lda, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zpotrf_mc( magma_tally3_context *cntxt, char *uplo, magma_tally3_int_t *n, cuDoubleComplex *A,
                          magma_tally3_int_t *lda, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zpotri( char uplo, magma_tally3_int_t n, cuDoubleComplex *A,
                  magma_tally3_int_t lda, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zlauum( char uplo, magma_tally3_int_t n, cuDoubleComplex *A,
                  magma_tally3_int_t lda, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_ztrtri( char uplo, char diag, magma_tally3_int_t n, cuDoubleComplex *A, 
                  magma_tally3_int_t lda, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zhetrd( char uplo, magma_tally3_int_t n, cuDoubleComplex *A, 
              magma_tally3_int_t lda, double *d, double *e, 
              cuDoubleComplex *tau, cuDoubleComplex *work, magma_tally3_int_t lwork, 
              magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgeqrf2(magma_tally3_context *cntxt, magma_tally3_int_t m, magma_tally3_int_t n,
                          cuDoubleComplex *a, magma_tally3_int_t lda, cuDoubleComplex *tau,
                          cuDoubleComplex *work, magma_tally3_int_t lwork,
                          magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgeqrf3(magma_tally3_context *cntxt, magma_tally3_int_t m, magma_tally3_int_t n,
                          cuDoubleComplex *a, magma_tally3_int_t lda, cuDoubleComplex *tau,
                          cuDoubleComplex *work, magma_tally3_int_t lwork,
                          magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zungqr( magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
              cuDoubleComplex *a, magma_tally3_int_t lda,
              cuDoubleComplex *tau, cuDoubleComplex *dwork,
              magma_tally3_int_t nb, magma_tally3_int_t *info );
magma_tally3_int_t magma_tally3_zunmql( const char side, const char trans,
              magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
              cuDoubleComplex *a, magma_tally3_int_t lda,
              cuDoubleComplex *tau,
              cuDoubleComplex *c, magma_tally3_int_t ldc,
              cuDoubleComplex *work, magma_tally3_int_t lwork,
              magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zunmqr( char side, char trans, 
                          magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k, 
                          cuDoubleComplex *a, magma_tally3_int_t lda, cuDoubleComplex *tau, 
                          cuDoubleComplex *c, magma_tally3_int_t ldc, 
                          cuDoubleComplex *work, magma_tally3_int_t lwork, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zunmtr( char side, char uplo, char trans,
              magma_tally3_int_t m, magma_tally3_int_t n,
              cuDoubleComplex *a,    magma_tally3_int_t lda,
              cuDoubleComplex *tau,
              cuDoubleComplex *c,    magma_tally3_int_t ldc,
              cuDoubleComplex *work, magma_tally3_int_t lwork,
              magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zunghr( magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
              cuDoubleComplex *a, magma_tally3_int_t lda,
              cuDoubleComplex *tau,
              cuDoubleComplex *dT, magma_tally3_int_t nb,
              magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zheev( char jobz, char uplo, magma_tally3_int_t n,
             cuDoubleComplex *a, magma_tally3_int_t lda, double *w,
             cuDoubleComplex *work, magma_tally3_int_t lwork,
                 double *rwork, magma_tally3_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
magma_tally3_int_t  magma_tally3_zgeev( char jobvl, char jobvr, magma_tally3_int_t n,
              cuDoubleComplex *a, magma_tally3_int_t lda,
              cuDoubleComplex *w,
              cuDoubleComplex *vl, magma_tally3_int_t ldvl,
              cuDoubleComplex *vr, magma_tally3_int_t ldvr,
              cuDoubleComplex *work, magma_tally3_int_t lwork,
              double *rwork, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgesvd( char jobu, char jobvt, magma_tally3_int_t m, magma_tally3_int_t n,
              cuDoubleComplex *a,    magma_tally3_int_t lda, double *s, 
                          cuDoubleComplex *u,    magma_tally3_int_t ldu, 
              cuDoubleComplex *vt,   magma_tally3_int_t ldvt,
              cuDoubleComplex *work, magma_tally3_int_t lwork,
              double *rwork, magma_tally3_int_t *info );
magma_tally3_int_t magma_tally3_zheevd( char jobz, char uplo, magma_tally3_int_t n,
              cuDoubleComplex *a, magma_tally3_int_t lda, double *w,
              cuDoubleComplex *work, magma_tally3_int_t lwork,
                  double *rwork, magma_tally3_int_t lrwork,
                  magma_tally3_int_t *iwork, magma_tally3_int_t liwork, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zhegvd( magma_tally3_int_t itype, char jobz, char uplo, magma_tally3_int_t n,
              cuDoubleComplex *a, magma_tally3_int_t lda,
              cuDoubleComplex *b, magma_tally3_int_t ldb,
              double *w, cuDoubleComplex *work, magma_tally3_int_t lwork,
              double *rwork, magma_tally3_int_t lrwork, magma_tally3_int_t *iwork,
              magma_tally3_int_t liwork, magma_tally3_int_t *info);
#else
magma_tally3_int_t  magma_tally3_zgeev( char jobvl, char jobvr, magma_tally3_int_t n,
              cuDoubleComplex *a,    magma_tally3_int_t lda,
              cuDoubleComplex *wr, cuDoubleComplex *wi,
              cuDoubleComplex *vl,   magma_tally3_int_t ldvl,
              cuDoubleComplex *vr,   magma_tally3_int_t ldvr,
              cuDoubleComplex *work, magma_tally3_int_t lwork,
              magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgesvd( char jobu, char jobvt, magma_tally3_int_t m, magma_tally3_int_t n,
              cuDoubleComplex *a,    magma_tally3_int_t lda, double *s, 
                          cuDoubleComplex *u,    magma_tally3_int_t ldu, 
              cuDoubleComplex *vt,   magma_tally3_int_t ldvt,
              cuDoubleComplex *work, magma_tally3_int_t lwork,
              magma_tally3_int_t *info );
magma_tally3_int_t magma_tally3_zheevd( char jobz, char uplo, magma_tally3_int_t n,
              cuDoubleComplex *a, magma_tally3_int_t lda, double *w,
              cuDoubleComplex *work, magma_tally3_int_t lwork,
              magma_tally3_int_t *iwork, magma_tally3_int_t liwork, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zhegvd( magma_tally3_int_t itype, char jobz, char uplo, magma_tally3_int_t n,
              cuDoubleComplex *a, magma_tally3_int_t lda,
              cuDoubleComplex *b, magma_tally3_int_t ldb,
              double *w, cuDoubleComplex *work, magma_tally3_int_t lwork,
              magma_tally3_int_t *iwork, magma_tally3_int_t liwork, magma_tally3_int_t *info);
#endif

magma_tally3_int_t magma_tally3_zhegst( magma_tally3_int_t itype, char uplo, magma_tally3_int_t n,
              cuDoubleComplex *a, magma_tally3_int_t lda,
              cuDoubleComplex *b, magma_tally3_int_t ldb, magma_tally3_int_t *info);

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA_tally3 function definitions / Data on GPU
*/
magma_tally3_int_t magma_tally3_zgels_gpu(  char trans, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  cuDoubleComplex *dA,    magma_tally3_int_t ldda, 
                  cuDoubleComplex *dB,    magma_tally3_int_t lddb, 
                  cuDoubleComplex *hwork, magma_tally3_int_t lwork, 
                  magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgels3_gpu( char trans, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  cuDoubleComplex *dA,    magma_tally3_int_t ldda,
                  cuDoubleComplex *dB,    magma_tally3_int_t lddb,
                  cuDoubleComplex *hwork, magma_tally3_int_t lwork,
                  magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgelqf_gpu( magma_tally3_int_t m, magma_tally3_int_t n,
                  cuDoubleComplex *dA,    magma_tally3_int_t ldda,   cuDoubleComplex *tau,
                  cuDoubleComplex *work, magma_tally3_int_t lwork, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgeqrf_gpu( magma_tally3_int_t m, magma_tally3_int_t n, 
                  cuDoubleComplex *dA,  magma_tally3_int_t ldda, 
                  cuDoubleComplex *tau, cuDoubleComplex *dT, 
                  magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgeqrf2_gpu(magma_tally3_int_t m, magma_tally3_int_t n, 
                  cuDoubleComplex *dA,  magma_tally3_int_t ldda, 
                  cuDoubleComplex *tau, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgeqrf3_gpu(magma_tally3_int_t m, magma_tally3_int_t n, 
                  cuDoubleComplex *dA,  magma_tally3_int_t ldda, 
                  cuDoubleComplex *tau, cuDoubleComplex *dT, 
                  magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgeqrs_gpu( magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs, 
                  cuDoubleComplex *dA,     magma_tally3_int_t ldda, 
                  cuDoubleComplex *tau,   cuDoubleComplex *dT,
                  cuDoubleComplex *dB,    magma_tally3_int_t lddb,
                  cuDoubleComplex *hwork, magma_tally3_int_t lhwork, 
                  magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgeqrs3_gpu( magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs, 
                  cuDoubleComplex *dA,     magma_tally3_int_t ldda, 
                  cuDoubleComplex *tau,   cuDoubleComplex *dT,
                  cuDoubleComplex *dB,    magma_tally3_int_t lddb,
                  cuDoubleComplex *hwork, magma_tally3_int_t lhwork, 
                  magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgessm_gpu( char storev, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3_int_t ib, 
                              magma_tally3_int_t *ipiv, 
                              cuDoubleComplex *dL1, magma_tally3_int_t lddl1, 
                              cuDoubleComplex *dL,  magma_tally3_int_t lddl, 
                              cuDoubleComplex *dA,  magma_tally3_int_t ldda, 
                              magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgesv_gpu(  magma_tally3_int_t n, magma_tally3_int_t nrhs, 
                  cuDoubleComplex *dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv, 
                  cuDoubleComplex *dB, magma_tally3_int_t lddb, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgetrf_incpiv_gpu( char storev, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t ib,
                              cuDoubleComplex *hA, magma_tally3_int_t ldha, cuDoubleComplex *dA, magma_tally3_int_t ldda,
                              cuDoubleComplex *hL, magma_tally3_int_t ldhl, cuDoubleComplex *dL, magma_tally3_int_t lddl,
                              magma_tally3_int_t *ipiv, 
                              cuDoubleComplex *dwork, magma_tally3_int_t lddwork,
                              magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgetrf_gpu( magma_tally3_int_t m, magma_tally3_int_t n, 
                  cuDoubleComplex *dA, magma_tally3_int_t ldda, 
                  magma_tally3_int_t *ipiv, magma_tally3_int_t *info);
magma_tally3_int_t 
magma_tally3_zgetrf_nopiv_gpu      ( magma_tally3_int_t m, magma_tally3_int_t n,
                  cuDoubleComplex *dA, magma_tally3_int_t ldda,
                  magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zgetrs_gpu( char trans, magma_tally3_int_t n, magma_tally3_int_t nrhs, 
                  cuDoubleComplex *dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv, 
                  cuDoubleComplex *dB, magma_tally3_int_t lddb, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zlabrd_gpu( magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb, 
                              cuDoubleComplex *a, magma_tally3_int_t lda, cuDoubleComplex *da, magma_tally3_int_t ldda,
                              double *d, double *e, cuDoubleComplex *tauq, cuDoubleComplex *taup,  
                              cuDoubleComplex *x, magma_tally3_int_t ldx, cuDoubleComplex *dx, magma_tally3_int_t lddx, 
                              cuDoubleComplex *y, magma_tally3_int_t ldy, cuDoubleComplex *dy, magma_tally3_int_t lddy);
magma_tally3_int_t magma_tally3_zlarfb_gpu( char side, char trans, char direct, char storev, 
                  magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
                  cuDoubleComplex *dv, magma_tally3_int_t ldv, cuDoubleComplex *dt,    magma_tally3_int_t ldt, 
                  cuDoubleComplex *dc, magma_tally3_int_t ldc, cuDoubleComplex *dowrk, magma_tally3_int_t ldwork );
magma_tally3_int_t magma_tally3_zposv_gpu(  char uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs, 
                  cuDoubleComplex *dA, magma_tally3_int_t ldda, 
                  cuDoubleComplex *dB, magma_tally3_int_t lddb, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zpotrf_gpu( char uplo,  magma_tally3_int_t n, 
                  cuDoubleComplex *dA, magma_tally3_int_t ldda, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zpotri_gpu( char uplo,  magma_tally3_int_t n,
                      cuDoubleComplex *dA, magma_tally3_int_t ldda, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zlauum_gpu( char uplo,  magma_tally3_int_t n,
                      cuDoubleComplex *dA, magma_tally3_int_t ldda, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_ztrtri_gpu( char uplo,  char diag, magma_tally3_int_t n,
                      cuDoubleComplex *dA, magma_tally3_int_t ldda, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zhetrd_gpu( char uplo, magma_tally3_int_t n,
                  cuDoubleComplex *da, magma_tally3_int_t ldda,
                  double *d, double *e, cuDoubleComplex *tau,
                  cuDoubleComplex *wa,  magma_tally3_int_t ldwa,
                  cuDoubleComplex *work, magma_tally3_int_t lwork,
                  magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zhetrd2_gpu(char uplo, magma_tally3_int_t n,
                  cuDoubleComplex *da, magma_tally3_int_t ldda,
                  double *d, double *e, cuDoubleComplex *tau,
                  cuDoubleComplex *wa,  magma_tally3_int_t ldwa,
                  cuDoubleComplex *work, magma_tally3_int_t lwork,
                  cuDoubleComplex *dwork, magma_tally3_int_t ldwork,
                  magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zpotrs_gpu( char uplo,  magma_tally3_int_t n, magma_tally3_int_t nrhs, 
                  cuDoubleComplex *dA, magma_tally3_int_t ldda, 
                  cuDoubleComplex *dB, magma_tally3_int_t lddb, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zssssm_gpu( char storev, magma_tally3_int_t m1, magma_tally3_int_t n1, 
                              magma_tally3_int_t m2, magma_tally3_int_t n2, magma_tally3_int_t k, magma_tally3_int_t ib, 
                              cuDoubleComplex *dA1, magma_tally3_int_t ldda1, 
                              cuDoubleComplex *dA2, magma_tally3_int_t ldda2, 
                              cuDoubleComplex *dL1, magma_tally3_int_t lddl1, 
                              cuDoubleComplex *dL2, magma_tally3_int_t lddl2,
                              magma_tally3_int_t *IPIV, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_ztstrf_gpu( char storev, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t ib, magma_tally3_int_t nb,
                              cuDoubleComplex *hU, magma_tally3_int_t ldhu, cuDoubleComplex *dU, magma_tally3_int_t lddu, 
                              cuDoubleComplex *hA, magma_tally3_int_t ldha, cuDoubleComplex *dA, magma_tally3_int_t ldda, 
                              cuDoubleComplex *hL, magma_tally3_int_t ldhl, cuDoubleComplex *dL, magma_tally3_int_t lddl,
                              magma_tally3_int_t *ipiv, 
                              cuDoubleComplex *hwork, magma_tally3_int_t ldhwork, 
                  cuDoubleComplex *dwork, magma_tally3_int_t lddwork,
                              magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zungqr_gpu( magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k, 
                              cuDoubleComplex *da, magma_tally3_int_t ldda, 
                              cuDoubleComplex *tau, cuDoubleComplex *dwork, 
                              magma_tally3_int_t nb, magma_tally3_int_t *info );
magma_tally3_int_t magma_tally3_zunmql2_gpu(const char side, const char trans,
                  magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
                  cuDoubleComplex *da, magma_tally3_int_t ldda,
                  cuDoubleComplex *tau,
                  cuDoubleComplex *dc, magma_tally3_int_t lddc,
                  cuDoubleComplex *wa, magma_tally3_int_t ldwa,
                  magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zunmqr_gpu( char side, char trans, 
                              magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
                              cuDoubleComplex *a,    magma_tally3_int_t lda, cuDoubleComplex *tau, 
                              cuDoubleComplex *c,    magma_tally3_int_t ldc,
                              cuDoubleComplex *work, magma_tally3_int_t lwork, 
                              cuDoubleComplex *td,   magma_tally3_int_t nb, magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zunmqr2_gpu(const char side, const char trans,
                  magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
                  cuDoubleComplex *da,   magma_tally3_int_t ldda,
                  cuDoubleComplex *tau,
                  cuDoubleComplex *dc,    magma_tally3_int_t lddc,
                  cuDoubleComplex *wa,    magma_tally3_int_t ldwa,
                  magma_tally3_int_t *info);
magma_tally3_int_t magma_tally3_zunmtr_gpu( char side, char uplo, char trans,
                  magma_tally3_int_t m, magma_tally3_int_t n,
                  cuDoubleComplex *da,    magma_tally3_int_t ldda,
                  cuDoubleComplex *tau,
                  cuDoubleComplex *dc,    magma_tally3_int_t lddc,
                  cuDoubleComplex *wa,    magma_tally3_int_t ldwa,
                  magma_tally3_int_t *info);

#if defined(PRECISION_z) || defined(PRECISION_c)
magma_tally3_int_t magma_tally3_zheevd_gpu( char jobz, char uplo,
                  magma_tally3_int_t n,
                  cuDoubleComplex *da, magma_tally3_int_t ldda,
                  double *w,
                  cuDoubleComplex *wa,  magma_tally3_int_t ldwa,
                  cuDoubleComplex *work, magma_tally3_int_t lwork,
                  double *rwork, magma_tally3_int_t lrwork,
                  magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
                  magma_tally3_int_t *info);
#else
magma_tally3_int_t magma_tally3_zheevd_gpu( char jobz, char uplo,
                  magma_tally3_int_t n,
                  cuDoubleComplex *da, magma_tally3_int_t ldda,
                  cuDoubleComplex *w,
                  cuDoubleComplex *wa,  magma_tally3_int_t ldwa,
                  cuDoubleComplex *work, magma_tally3_int_t lwork,
                  magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
                  magma_tally3_int_t *info);
#endif

magma_tally3_int_t magma_tally3_zhegst_gpu(magma_tally3_int_t itype, char uplo, magma_tally3_int_t n,
                 cuDoubleComplex *da, magma_tally3_int_t ldda,
                 cuDoubleComplex *db, magma_tally3_int_t lddb, magma_tally3_int_t *info);


#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* _MAGMA_tally3_Z_H_ */

