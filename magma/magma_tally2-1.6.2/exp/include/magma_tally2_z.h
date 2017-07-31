/*
 *   -- MAGMA_tally2 (version 1.6.1) --
 *      Univ. of Tennessee, Knoxville
 *      Univ. of California, Berkeley
 *      Univ. of Colorado, Denver
 *      @date January 2015
 *
 * @precisions normal z -> s d c
 */

#ifndef _MAGMA_tally2_Z_H_
#define _MAGMA_tally2_Z_H_
#define PRECISION_z

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 function definitions / Data on CPU
*/
magma_tally2_int_t magma_tally2_zgebrd( magma_tally2_int_t m, magma_tally2_int_t n, cuDoubleComplex *A, 
              magma_tally2_int_t lda, double *d, double *e,
              cuDoubleComplex *tauq,  cuDoubleComplex *taup, 
              cuDoubleComplex *work, magma_tally2_int_t lwork, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgehrd2(magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
              cuDoubleComplex *A, magma_tally2_int_t lda, cuDoubleComplex *tau, 
              cuDoubleComplex *work, magma_tally2_int_t *lwork, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgehrd( magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
              cuDoubleComplex *A, magma_tally2_int_t lda, cuDoubleComplex *tau,
              cuDoubleComplex *work, magma_tally2_int_t lwork,
              cuDoubleComplex *d_T, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgelqf( magma_tally2_int_t m, magma_tally2_int_t n, 
                          cuDoubleComplex *A,    magma_tally2_int_t lda,   cuDoubleComplex *tau, 
                          cuDoubleComplex *work, magma_tally2_int_t lwork, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgeqlf( magma_tally2_int_t m, magma_tally2_int_t n, 
                          cuDoubleComplex *A,    magma_tally2_int_t lda,   cuDoubleComplex *tau, 
                          cuDoubleComplex *work, magma_tally2_int_t lwork, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgeqrf( magma_tally2_int_t m, magma_tally2_int_t n, cuDoubleComplex *A, 
              magma_tally2_int_t lda, cuDoubleComplex *tau, cuDoubleComplex *work, 
              magma_tally2_int_t lwork, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgeqrf_ooc( magma_tally2_int_t m, magma_tally2_int_t n, cuDoubleComplex *A,
              magma_tally2_int_t lda, cuDoubleComplex *tau, cuDoubleComplex *work,
                  magma_tally2_int_t lwork, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgesv ( magma_tally2_int_t n, magma_tally2_int_t nrhs, 
              cuDoubleComplex *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv, 
              cuDoubleComplex *B, magma_tally2_int_t ldb, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgetrf( magma_tally2_int_t m, magma_tally2_int_t n, cuDoubleComplex *A, 
              magma_tally2_int_t lda, magma_tally2_int_t *ipiv, 
              magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgetrf_mc( magma_tally2_context *cntxt, magma_tally2_int_t *m, magma_tally2_int_t *n, cuDoubleComplex *A,
                          magma_tally2_int_t *lda, magma_tally2_int_t *ipiv, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgeqrf_mc(magma_tally2_context *cntxt, magma_tally2_int_t *m, magma_tally2_int_t *n,
                          cuDoubleComplex *A, magma_tally2_int_t *lda,
                          cuDoubleComplex *tau, cuDoubleComplex *work,
                          magma_tally2_int_t *lwork, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgetrf2(magma_tally2_int_t m, magma_tally2_int_t n, cuDoubleComplex *a, 
                          magma_tally2_int_t lda, magma_tally2_int_t *ipiv, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zlatrd( char uplo, magma_tally2_int_t n, magma_tally2_int_t nb, cuDoubleComplex *a, 
                          magma_tally2_int_t lda, double *e, cuDoubleComplex *tau, 
              cuDoubleComplex *w, magma_tally2_int_t ldw,
                          cuDoubleComplex *da, magma_tally2_int_t ldda, 
              cuDoubleComplex *dw, magma_tally2_int_t lddw);
magma_tally2_int_t magma_tally2_zlatrd2(char uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
              cuDoubleComplex *a,  magma_tally2_int_t lda,
              double *e, cuDoubleComplex *tau,
              cuDoubleComplex *w,  magma_tally2_int_t ldw,
              cuDoubleComplex *da, magma_tally2_int_t ldda,
              cuDoubleComplex *dw, magma_tally2_int_t lddw,
              cuDoubleComplex *dwork, magma_tally2_int_t ldwork);
magma_tally2_int_t magma_tally2_zlahr2( magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb, 
              cuDoubleComplex *da, cuDoubleComplex *dv, cuDoubleComplex *a, 
              magma_tally2_int_t lda, cuDoubleComplex *tau, cuDoubleComplex *t, 
              magma_tally2_int_t ldt, cuDoubleComplex *y, magma_tally2_int_t ldy);
magma_tally2_int_t magma_tally2_zlahru( magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb, 
              cuDoubleComplex *a, magma_tally2_int_t lda, 
              cuDoubleComplex *da, cuDoubleComplex *y, 
              cuDoubleComplex *v, cuDoubleComplex *t, 
              cuDoubleComplex *dwork);
magma_tally2_int_t magma_tally2_zposv ( char uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs, 
              cuDoubleComplex *A, magma_tally2_int_t lda, 
              cuDoubleComplex *B, magma_tally2_int_t ldb, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zpotrf( char uplo, magma_tally2_int_t n, cuDoubleComplex *A, 
              magma_tally2_int_t lda, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zpotrf_mc( magma_tally2_context *cntxt, char *uplo, magma_tally2_int_t *n, cuDoubleComplex *A,
                          magma_tally2_int_t *lda, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zpotri( char uplo, magma_tally2_int_t n, cuDoubleComplex *A,
                  magma_tally2_int_t lda, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zlauum( char uplo, magma_tally2_int_t n, cuDoubleComplex *A,
                  magma_tally2_int_t lda, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_ztrtri( char uplo, char diag, magma_tally2_int_t n, cuDoubleComplex *A, 
                  magma_tally2_int_t lda, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zhetrd( char uplo, magma_tally2_int_t n, cuDoubleComplex *A, 
              magma_tally2_int_t lda, double *d, double *e, 
              cuDoubleComplex *tau, cuDoubleComplex *work, magma_tally2_int_t lwork, 
              magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgeqrf2(magma_tally2_context *cntxt, magma_tally2_int_t m, magma_tally2_int_t n,
                          cuDoubleComplex *a, magma_tally2_int_t lda, cuDoubleComplex *tau,
                          cuDoubleComplex *work, magma_tally2_int_t lwork,
                          magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgeqrf3(magma_tally2_context *cntxt, magma_tally2_int_t m, magma_tally2_int_t n,
                          cuDoubleComplex *a, magma_tally2_int_t lda, cuDoubleComplex *tau,
                          cuDoubleComplex *work, magma_tally2_int_t lwork,
                          magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zungqr( magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
              cuDoubleComplex *a, magma_tally2_int_t lda,
              cuDoubleComplex *tau, cuDoubleComplex *dwork,
              magma_tally2_int_t nb, magma_tally2_int_t *info );
magma_tally2_int_t magma_tally2_zunmql( const char side, const char trans,
              magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
              cuDoubleComplex *a, magma_tally2_int_t lda,
              cuDoubleComplex *tau,
              cuDoubleComplex *c, magma_tally2_int_t ldc,
              cuDoubleComplex *work, magma_tally2_int_t lwork,
              magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zunmqr( char side, char trans, 
                          magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k, 
                          cuDoubleComplex *a, magma_tally2_int_t lda, cuDoubleComplex *tau, 
                          cuDoubleComplex *c, magma_tally2_int_t ldc, 
                          cuDoubleComplex *work, magma_tally2_int_t lwork, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zunmtr( char side, char uplo, char trans,
              magma_tally2_int_t m, magma_tally2_int_t n,
              cuDoubleComplex *a,    magma_tally2_int_t lda,
              cuDoubleComplex *tau,
              cuDoubleComplex *c,    magma_tally2_int_t ldc,
              cuDoubleComplex *work, magma_tally2_int_t lwork,
              magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zunghr( magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
              cuDoubleComplex *a, magma_tally2_int_t lda,
              cuDoubleComplex *tau,
              cuDoubleComplex *dT, magma_tally2_int_t nb,
              magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zheev( char jobz, char uplo, magma_tally2_int_t n,
             cuDoubleComplex *a, magma_tally2_int_t lda, double *w,
             cuDoubleComplex *work, magma_tally2_int_t lwork,
                 double *rwork, magma_tally2_int_t *info);
#if defined(PRECISION_z) || defined(PRECISION_c)
magma_tally2_int_t  magma_tally2_zgeev( char jobvl, char jobvr, magma_tally2_int_t n,
              cuDoubleComplex *a, magma_tally2_int_t lda,
              cuDoubleComplex *w,
              cuDoubleComplex *vl, magma_tally2_int_t ldvl,
              cuDoubleComplex *vr, magma_tally2_int_t ldvr,
              cuDoubleComplex *work, magma_tally2_int_t lwork,
              double *rwork, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgesvd( char jobu, char jobvt, magma_tally2_int_t m, magma_tally2_int_t n,
              cuDoubleComplex *a,    magma_tally2_int_t lda, double *s, 
                          cuDoubleComplex *u,    magma_tally2_int_t ldu, 
              cuDoubleComplex *vt,   magma_tally2_int_t ldvt,
              cuDoubleComplex *work, magma_tally2_int_t lwork,
              double *rwork, magma_tally2_int_t *info );
magma_tally2_int_t magma_tally2_zheevd( char jobz, char uplo, magma_tally2_int_t n,
              cuDoubleComplex *a, magma_tally2_int_t lda, double *w,
              cuDoubleComplex *work, magma_tally2_int_t lwork,
                  double *rwork, magma_tally2_int_t lrwork,
                  magma_tally2_int_t *iwork, magma_tally2_int_t liwork, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zhegvd( magma_tally2_int_t itype, char jobz, char uplo, magma_tally2_int_t n,
              cuDoubleComplex *a, magma_tally2_int_t lda,
              cuDoubleComplex *b, magma_tally2_int_t ldb,
              double *w, cuDoubleComplex *work, magma_tally2_int_t lwork,
              double *rwork, magma_tally2_int_t lrwork, magma_tally2_int_t *iwork,
              magma_tally2_int_t liwork, magma_tally2_int_t *info);
#else
magma_tally2_int_t  magma_tally2_zgeev( char jobvl, char jobvr, magma_tally2_int_t n,
              cuDoubleComplex *a,    magma_tally2_int_t lda,
              cuDoubleComplex *wr, cuDoubleComplex *wi,
              cuDoubleComplex *vl,   magma_tally2_int_t ldvl,
              cuDoubleComplex *vr,   magma_tally2_int_t ldvr,
              cuDoubleComplex *work, magma_tally2_int_t lwork,
              magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgesvd( char jobu, char jobvt, magma_tally2_int_t m, magma_tally2_int_t n,
              cuDoubleComplex *a,    magma_tally2_int_t lda, double *s, 
                          cuDoubleComplex *u,    magma_tally2_int_t ldu, 
              cuDoubleComplex *vt,   magma_tally2_int_t ldvt,
              cuDoubleComplex *work, magma_tally2_int_t lwork,
              magma_tally2_int_t *info );
magma_tally2_int_t magma_tally2_zheevd( char jobz, char uplo, magma_tally2_int_t n,
              cuDoubleComplex *a, magma_tally2_int_t lda, double *w,
              cuDoubleComplex *work, magma_tally2_int_t lwork,
              magma_tally2_int_t *iwork, magma_tally2_int_t liwork, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zhegvd( magma_tally2_int_t itype, char jobz, char uplo, magma_tally2_int_t n,
              cuDoubleComplex *a, magma_tally2_int_t lda,
              cuDoubleComplex *b, magma_tally2_int_t ldb,
              double *w, cuDoubleComplex *work, magma_tally2_int_t lwork,
              magma_tally2_int_t *iwork, magma_tally2_int_t liwork, magma_tally2_int_t *info);
#endif

magma_tally2_int_t magma_tally2_zhegst( magma_tally2_int_t itype, char uplo, magma_tally2_int_t n,
              cuDoubleComplex *a, magma_tally2_int_t lda,
              cuDoubleComplex *b, magma_tally2_int_t ldb, magma_tally2_int_t *info);

/* //////////////////////////////////////////////////////////////////////////// 
 -- MAGMA_tally2 function definitions / Data on GPU
*/
magma_tally2_int_t magma_tally2_zgels_gpu(  char trans, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
                  cuDoubleComplex *dA,    magma_tally2_int_t ldda, 
                  cuDoubleComplex *dB,    magma_tally2_int_t lddb, 
                  cuDoubleComplex *hwork, magma_tally2_int_t lwork, 
                  magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgels3_gpu( char trans, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
                  cuDoubleComplex *dA,    magma_tally2_int_t ldda,
                  cuDoubleComplex *dB,    magma_tally2_int_t lddb,
                  cuDoubleComplex *hwork, magma_tally2_int_t lwork,
                  magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgelqf_gpu( magma_tally2_int_t m, magma_tally2_int_t n,
                  cuDoubleComplex *dA,    magma_tally2_int_t ldda,   cuDoubleComplex *tau,
                  cuDoubleComplex *work, magma_tally2_int_t lwork, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgeqrf_gpu( magma_tally2_int_t m, magma_tally2_int_t n, 
                  cuDoubleComplex *dA,  magma_tally2_int_t ldda, 
                  cuDoubleComplex *tau, cuDoubleComplex *dT, 
                  magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgeqrf2_gpu(magma_tally2_int_t m, magma_tally2_int_t n, 
                  cuDoubleComplex *dA,  magma_tally2_int_t ldda, 
                  cuDoubleComplex *tau, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgeqrf3_gpu(magma_tally2_int_t m, magma_tally2_int_t n, 
                  cuDoubleComplex *dA,  magma_tally2_int_t ldda, 
                  cuDoubleComplex *tau, cuDoubleComplex *dT, 
                  magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgeqrs_gpu( magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs, 
                  cuDoubleComplex *dA,     magma_tally2_int_t ldda, 
                  cuDoubleComplex *tau,   cuDoubleComplex *dT,
                  cuDoubleComplex *dB,    magma_tally2_int_t lddb,
                  cuDoubleComplex *hwork, magma_tally2_int_t lhwork, 
                  magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgeqrs3_gpu( magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs, 
                  cuDoubleComplex *dA,     magma_tally2_int_t ldda, 
                  cuDoubleComplex *tau,   cuDoubleComplex *dT,
                  cuDoubleComplex *dB,    magma_tally2_int_t lddb,
                  cuDoubleComplex *hwork, magma_tally2_int_t lhwork, 
                  magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgessm_gpu( char storev, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2_int_t ib, 
                              magma_tally2_int_t *ipiv, 
                              cuDoubleComplex *dL1, magma_tally2_int_t lddl1, 
                              cuDoubleComplex *dL,  magma_tally2_int_t lddl, 
                              cuDoubleComplex *dA,  magma_tally2_int_t ldda, 
                              magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgesv_gpu(  magma_tally2_int_t n, magma_tally2_int_t nrhs, 
                  cuDoubleComplex *dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv, 
                  cuDoubleComplex *dB, magma_tally2_int_t lddb, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgetrf_incpiv_gpu( char storev, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t ib,
                              cuDoubleComplex *hA, magma_tally2_int_t ldha, cuDoubleComplex *dA, magma_tally2_int_t ldda,
                              cuDoubleComplex *hL, magma_tally2_int_t ldhl, cuDoubleComplex *dL, magma_tally2_int_t lddl,
                              magma_tally2_int_t *ipiv, 
                              cuDoubleComplex *dwork, magma_tally2_int_t lddwork,
                              magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgetrf_gpu( magma_tally2_int_t m, magma_tally2_int_t n, 
                  cuDoubleComplex *dA, magma_tally2_int_t ldda, 
                  magma_tally2_int_t *ipiv, magma_tally2_int_t *info);
magma_tally2_int_t 
magma_tally2_zgetrf_nopiv_gpu      ( magma_tally2_int_t m, magma_tally2_int_t n,
                  cuDoubleComplex *dA, magma_tally2_int_t ldda,
                  magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zgetrs_gpu( char trans, magma_tally2_int_t n, magma_tally2_int_t nrhs, 
                  cuDoubleComplex *dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv, 
                  cuDoubleComplex *dB, magma_tally2_int_t lddb, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zlabrd_gpu( magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb, 
                              cuDoubleComplex *a, magma_tally2_int_t lda, cuDoubleComplex *da, magma_tally2_int_t ldda,
                              double *d, double *e, cuDoubleComplex *tauq, cuDoubleComplex *taup,  
                              cuDoubleComplex *x, magma_tally2_int_t ldx, cuDoubleComplex *dx, magma_tally2_int_t lddx, 
                              cuDoubleComplex *y, magma_tally2_int_t ldy, cuDoubleComplex *dy, magma_tally2_int_t lddy);
magma_tally2_int_t magma_tally2_zlarfb_gpu( char side, char trans, char direct, char storev, 
                  magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
                  cuDoubleComplex *dv, magma_tally2_int_t ldv, cuDoubleComplex *dt,    magma_tally2_int_t ldt, 
                  cuDoubleComplex *dc, magma_tally2_int_t ldc, cuDoubleComplex *dowrk, magma_tally2_int_t ldwork );
magma_tally2_int_t magma_tally2_zposv_gpu(  char uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs, 
                  cuDoubleComplex *dA, magma_tally2_int_t ldda, 
                  cuDoubleComplex *dB, magma_tally2_int_t lddb, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zpotrf_gpu( char uplo,  magma_tally2_int_t n, 
                  cuDoubleComplex *dA, magma_tally2_int_t ldda, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zpotri_gpu( char uplo,  magma_tally2_int_t n,
                      cuDoubleComplex *dA, magma_tally2_int_t ldda, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zlauum_gpu( char uplo,  magma_tally2_int_t n,
                      cuDoubleComplex *dA, magma_tally2_int_t ldda, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_ztrtri_gpu( char uplo,  char diag, magma_tally2_int_t n,
                      cuDoubleComplex *dA, magma_tally2_int_t ldda, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zhetrd_gpu( char uplo, magma_tally2_int_t n,
                  cuDoubleComplex *da, magma_tally2_int_t ldda,
                  double *d, double *e, cuDoubleComplex *tau,
                  cuDoubleComplex *wa,  magma_tally2_int_t ldwa,
                  cuDoubleComplex *work, magma_tally2_int_t lwork,
                  magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zhetrd2_gpu(char uplo, magma_tally2_int_t n,
                  cuDoubleComplex *da, magma_tally2_int_t ldda,
                  double *d, double *e, cuDoubleComplex *tau,
                  cuDoubleComplex *wa,  magma_tally2_int_t ldwa,
                  cuDoubleComplex *work, magma_tally2_int_t lwork,
                  cuDoubleComplex *dwork, magma_tally2_int_t ldwork,
                  magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zpotrs_gpu( char uplo,  magma_tally2_int_t n, magma_tally2_int_t nrhs, 
                  cuDoubleComplex *dA, magma_tally2_int_t ldda, 
                  cuDoubleComplex *dB, magma_tally2_int_t lddb, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zssssm_gpu( char storev, magma_tally2_int_t m1, magma_tally2_int_t n1, 
                              magma_tally2_int_t m2, magma_tally2_int_t n2, magma_tally2_int_t k, magma_tally2_int_t ib, 
                              cuDoubleComplex *dA1, magma_tally2_int_t ldda1, 
                              cuDoubleComplex *dA2, magma_tally2_int_t ldda2, 
                              cuDoubleComplex *dL1, magma_tally2_int_t lddl1, 
                              cuDoubleComplex *dL2, magma_tally2_int_t lddl2,
                              magma_tally2_int_t *IPIV, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_ztstrf_gpu( char storev, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t ib, magma_tally2_int_t nb,
                              cuDoubleComplex *hU, magma_tally2_int_t ldhu, cuDoubleComplex *dU, magma_tally2_int_t lddu, 
                              cuDoubleComplex *hA, magma_tally2_int_t ldha, cuDoubleComplex *dA, magma_tally2_int_t ldda, 
                              cuDoubleComplex *hL, magma_tally2_int_t ldhl, cuDoubleComplex *dL, magma_tally2_int_t lddl,
                              magma_tally2_int_t *ipiv, 
                              cuDoubleComplex *hwork, magma_tally2_int_t ldhwork, 
                  cuDoubleComplex *dwork, magma_tally2_int_t lddwork,
                              magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zungqr_gpu( magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k, 
                              cuDoubleComplex *da, magma_tally2_int_t ldda, 
                              cuDoubleComplex *tau, cuDoubleComplex *dwork, 
                              magma_tally2_int_t nb, magma_tally2_int_t *info );
magma_tally2_int_t magma_tally2_zunmql2_gpu(const char side, const char trans,
                  magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
                  cuDoubleComplex *da, magma_tally2_int_t ldda,
                  cuDoubleComplex *tau,
                  cuDoubleComplex *dc, magma_tally2_int_t lddc,
                  cuDoubleComplex *wa, magma_tally2_int_t ldwa,
                  magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zunmqr_gpu( char side, char trans, 
                              magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
                              cuDoubleComplex *a,    magma_tally2_int_t lda, cuDoubleComplex *tau, 
                              cuDoubleComplex *c,    magma_tally2_int_t ldc,
                              cuDoubleComplex *work, magma_tally2_int_t lwork, 
                              cuDoubleComplex *td,   magma_tally2_int_t nb, magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zunmqr2_gpu(const char side, const char trans,
                  magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
                  cuDoubleComplex *da,   magma_tally2_int_t ldda,
                  cuDoubleComplex *tau,
                  cuDoubleComplex *dc,    magma_tally2_int_t lddc,
                  cuDoubleComplex *wa,    magma_tally2_int_t ldwa,
                  magma_tally2_int_t *info);
magma_tally2_int_t magma_tally2_zunmtr_gpu( char side, char uplo, char trans,
                  magma_tally2_int_t m, magma_tally2_int_t n,
                  cuDoubleComplex *da,    magma_tally2_int_t ldda,
                  cuDoubleComplex *tau,
                  cuDoubleComplex *dc,    magma_tally2_int_t lddc,
                  cuDoubleComplex *wa,    magma_tally2_int_t ldwa,
                  magma_tally2_int_t *info);

#if defined(PRECISION_z) || defined(PRECISION_c)
magma_tally2_int_t magma_tally2_zheevd_gpu( char jobz, char uplo,
                  magma_tally2_int_t n,
                  cuDoubleComplex *da, magma_tally2_int_t ldda,
                  double *w,
                  cuDoubleComplex *wa,  magma_tally2_int_t ldwa,
                  cuDoubleComplex *work, magma_tally2_int_t lwork,
                  double *rwork, magma_tally2_int_t lrwork,
                  magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
                  magma_tally2_int_t *info);
#else
magma_tally2_int_t magma_tally2_zheevd_gpu( char jobz, char uplo,
                  magma_tally2_int_t n,
                  cuDoubleComplex *da, magma_tally2_int_t ldda,
                  cuDoubleComplex *w,
                  cuDoubleComplex *wa,  magma_tally2_int_t ldwa,
                  cuDoubleComplex *work, magma_tally2_int_t lwork,
                  magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
                  magma_tally2_int_t *info);
#endif

magma_tally2_int_t magma_tally2_zhegst_gpu(magma_tally2_int_t itype, char uplo, magma_tally2_int_t n,
                 cuDoubleComplex *da, magma_tally2_int_t ldda,
                 cuDoubleComplex *db, magma_tally2_int_t lddb, magma_tally2_int_t *info);


#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* _MAGMA_tally2_Z_H_ */

