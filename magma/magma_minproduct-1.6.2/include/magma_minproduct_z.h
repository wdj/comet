/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
*/

#ifndef MAGMA_minproduct_Z_H
#define MAGMA_minproduct_Z_H

#include "magma_minproduct_types.h"
#include "magma_minproduct_zgehrd_m.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct Auxiliary functions to get the NB used
*/
#ifdef REAL
magma_minproduct_int_t magma_minproduct_get_dlaex3_m_nb();       // defined in dlaex3_m.cpp
#endif

magma_minproduct_int_t magma_minproduct_get_zpotrf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_zgetrf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_zgetri_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_zgeqp3_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_zgeqrf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_zgeqlf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_zgehrd_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_zhetrd_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_zhetrf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_zhetrf_nopiv_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_zgelqf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_zgebrd_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_zhegst_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_zgesvd_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_zhegst_nb_m( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_zbulge_nb( magma_minproduct_int_t m, magma_minproduct_int_t nbthreads );
magma_minproduct_int_t magma_minproduct_get_zbulge_nb_mgpu( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_zbulge_get_Vblksiz( magma_minproduct_int_t m, magma_minproduct_int_t nb, magma_minproduct_int_t nbthreads );
magma_minproduct_int_t magma_minproduct_get_zbulge_gcperf();


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct function definitions / Data on CPU (alphabetical order)
*/

#ifdef REAL
// only applicable to real [sd] precisions
void
magma_minproduct_dmove_eig(
    magma_minproduct_range_t range, magma_minproduct_int_t n, double *w,
    magma_minproduct_int_t *il, magma_minproduct_int_t *iu, double vl, double vu, magma_minproduct_int_t *m);

// defined in dlaex3.cpp
void
magma_minproduct_zvrange(
    magma_minproduct_int_t k, double *d, magma_minproduct_int_t *il, magma_minproduct_int_t *iu, double vl, double vu);

void
magma_minproduct_zirange(
    magma_minproduct_int_t k, magma_minproduct_int_t *indxq, magma_minproduct_int_t *iil, magma_minproduct_int_t *iiu, magma_minproduct_int_t il, magma_minproduct_int_t iu);
#endif

magma_minproduct_int_t
magma_minproduct_zgebrd(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    double *d, double *e,
    magma_minproductDoubleComplex *tauq, magma_minproductDoubleComplex *taup,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeev(
    magma_minproduct_vec_t jobvl, magma_minproduct_vec_t jobvr, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    #ifdef COMPLEX
    magma_minproductDoubleComplex *w,
    #else
    double *wr, double *wi,
    #endif
    magma_minproductDoubleComplex *VL, magma_minproduct_int_t ldvl,
    magma_minproductDoubleComplex *VR, magma_minproduct_int_t ldvr,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgehrd(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproductDoubleComplex_ptr dT,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgehrd2(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgelqf(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A,    magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqlf(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A,    magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqp3(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *jpvt, magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqrf(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqrf_ooc(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqrf4(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A,    magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgesdd(
    magma_minproduct_vec_t jobz, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    double *s,
    magma_minproductDoubleComplex *U, magma_minproduct_int_t ldu,
    magma_minproductDoubleComplex *VT, magma_minproduct_int_t ldvt,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_minproduct_int_t *iwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgesv(
    magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgesv_rbt(
    magma_minproduct_bool_t ref, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda, 
    magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgesvd(
    magma_minproduct_vec_t jobu, magma_minproduct_vec_t jobvt, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A,    magma_minproduct_int_t lda, double *s,
    magma_minproductDoubleComplex *U,    magma_minproduct_int_t ldu,
    magma_minproductDoubleComplex *VT,   magma_minproduct_int_t ldvt,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgetf2_nopiv(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgetrf(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgetrf_nopiv(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgetrf_piv(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t NB,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgetrf2(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zheevd(
    magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    double *w,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zheevdx(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zheevdx_2stage(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_minproduct_int_t
magma_minproduct_zheevr(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    double vl, double vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, double abstol, magma_minproduct_int_t *m,
    double *w,
    magma_minproductDoubleComplex *Z, magma_minproduct_int_t ldz,
    magma_minproduct_int_t *isuppz,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    double *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

// no real [sd] precisions available
magma_minproduct_int_t
magma_minproduct_zheevx(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    double vl, double vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, double abstol, magma_minproduct_int_t *m,
    double *w,
    magma_minproductDoubleComplex *Z, magma_minproduct_int_t ldz,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    double *rwork, magma_minproduct_int_t *iwork,
    magma_minproduct_int_t *ifail,
    magma_minproduct_int_t *info);
#endif  // COMPLEX

magma_minproduct_int_t
magma_minproduct_zhegst(
    magma_minproduct_int_t itype, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhegvd(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
    double *w, magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhegvdx(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n, magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhegvdx_2stage(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhegvr(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    double abstol, magma_minproduct_int_t *m, double *w,
    magma_minproductDoubleComplex *Z, magma_minproduct_int_t ldz,
    magma_minproduct_int_t *isuppz, magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    double *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhegvx(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n, magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    double abstol, magma_minproduct_int_t *m, double *w,
    magma_minproductDoubleComplex *Z, magma_minproduct_int_t ldz,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork, double *rwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t *ifail,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhesv(magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
            magma_minproductDoubleComplex *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
            magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
            magma_minproduct_int_t *info );

magma_minproduct_int_t
magma_minproduct_zhetrd(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    double *d, double *e, magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhetrf(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhetrf_nopiv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhetrd_hb2st(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t Vblksiz,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    double *d, double *e,
    magma_minproductDoubleComplex *V, magma_minproduct_int_t ldv,
    magma_minproductDoubleComplex *TAU, magma_minproduct_int_t compT,
    magma_minproductDoubleComplex *T, magma_minproduct_int_t ldt);

magma_minproduct_int_t
magma_minproduct_zhetrd_he2hb(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproductDoubleComplex_ptr dT,
    magma_minproduct_int_t *info);

#ifdef REAL
// only applicable to real [sd] precisions
magma_minproduct_int_t
magma_minproduct_dlaex0(
    magma_minproduct_int_t n, double *d, double *e,
    double *Q, magma_minproduct_int_t ldq,
    double *work, magma_minproduct_int_t *iwork,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_range_t range, double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dlaex1(
    magma_minproduct_int_t n, double *d,
    double *Q, magma_minproduct_int_t ldq,
    magma_minproduct_int_t *indxq, double rho, magma_minproduct_int_t cutpnt,
    double *work, magma_minproduct_int_t *iwork,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_range_t range, double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dlaex3(
    magma_minproduct_int_t k, magma_minproduct_int_t n, magma_minproduct_int_t n1, double *d,
    double *Q, magma_minproduct_int_t ldq,
    double rho,
    double *dlamda, double *Q2, magma_minproduct_int_t *indx,
    magma_minproduct_int_t *ctot, double *w, double *s, magma_minproduct_int_t *indxq,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_range_t range, double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *info);
#endif  // REAL

magma_minproduct_int_t
magma_minproduct_zlahef_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    magma_minproductDoubleComplex    *hA, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductDoubleComplex_ptr dW, magma_minproduct_int_t lddw,
    magma_minproduct_queue_t queues[], magma_minproduct_event_t event[],
    magma_minproduct_int_t *info);

magma_minproduct_int_t
zhetrf_nopiv_cpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t ib,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhetrs_nopiv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhesv_nopiv_gpu(
    magma_minproduct_uplo_t uplo,  magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb, 
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zlahr2(
    magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductDoubleComplex *A,  magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *T,  magma_minproduct_int_t ldt,
    magma_minproductDoubleComplex *Y,  magma_minproduct_int_t ldy);

magma_minproduct_int_t
magma_minproduct_zlahru(
    magma_minproduct_int_t n, magma_minproduct_int_t ihi, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex     *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dY, magma_minproduct_int_t lddy,
    magma_minproductDoubleComplex_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductDoubleComplex_ptr dT,
    magma_minproductDoubleComplex_ptr dwork);

#ifdef REAL
magma_minproduct_int_t
magma_minproduct_dlaln2(
    magma_minproduct_int_t trans, magma_minproduct_int_t na, magma_minproduct_int_t nw,
    double smin, double ca, const double *A, magma_minproduct_int_t lda,
    double d1, double d2,   const double *B, magma_minproduct_int_t ldb,
    double wr, double wi, double *X, magma_minproduct_int_t ldx,
    double *scale, double *xnorm,
    magma_minproduct_int_t *info);
#endif

magma_minproduct_int_t
magma_minproduct_zlaqps(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t offset,
    magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    magma_minproductDoubleComplex *A,  magma_minproduct_int_t lda,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt, magma_minproductDoubleComplex *tau, double *vn1, double *vn2,
    magma_minproductDoubleComplex *auxv,
    magma_minproductDoubleComplex *F,  magma_minproduct_int_t ldf,
    magma_minproductDoubleComplex_ptr dF, magma_minproduct_int_t lddf );

#ifdef REAL
magma_minproduct_int_t
magma_minproduct_zlaqtrsd(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n,
    const double *T, magma_minproduct_int_t ldt,
    double *x,       magma_minproduct_int_t ldx,
    const double *cnorm,
    magma_minproduct_int_t *info);
#endif

magma_minproduct_int_t
magma_minproduct_zlatrd(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    double *e, magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *W, magma_minproduct_int_t ldw,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dW, magma_minproduct_int_t lddw);

magma_minproduct_int_t
magma_minproduct_zlatrd2(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex *A,  magma_minproduct_int_t lda,
    double *e, magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *W,  magma_minproduct_int_t ldw,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dW, magma_minproduct_int_t lddw,
    magma_minproductDoubleComplex_ptr dwork, magma_minproduct_int_t ldwork);

#ifdef COMPLEX
magma_minproduct_int_t
magma_minproduct_zlatrsd(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_diag_t diag, magma_minproduct_bool_t normin,
    magma_minproduct_int_t n, const magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex lambda,
    magma_minproductDoubleComplex *x,
    double *scale, double *cnorm,
    magma_minproduct_int_t *info);
#endif

magma_minproduct_int_t
magma_minproduct_zlauum(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zposv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zpotrf(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zpotri(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zstedx(
    magma_minproduct_range_t range, magma_minproduct_int_t n, double vl, double vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, double *d, double *e,
    magma_minproductDoubleComplex *Z, magma_minproduct_int_t ldz,
    double *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ztrevc3(
    magma_minproduct_side_t side, magma_minproduct_vec_t howmany,
    magma_minproduct_int_t *select, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *T,  magma_minproduct_int_t ldt,
    magma_minproductDoubleComplex *VL, magma_minproduct_int_t ldvl,
    magma_minproductDoubleComplex *VR, magma_minproduct_int_t ldvr,
    magma_minproduct_int_t mm, magma_minproduct_int_t *mout,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ztrevc3_mt(
    magma_minproduct_side_t side, magma_minproduct_vec_t howmany,
    magma_minproduct_int_t *select, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *T,  magma_minproduct_int_t ldt,
    magma_minproductDoubleComplex *VL, magma_minproduct_int_t ldvl,
    magma_minproductDoubleComplex *VR, magma_minproduct_int_t ldvr,
    magma_minproduct_int_t mm, magma_minproduct_int_t *mout,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ztrtri(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zunghr(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zungqr(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zungqr2(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zunmbr(
    magma_minproduct_vect_t vect, magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *C, magma_minproduct_int_t ldc,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zunmlq(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *C, magma_minproduct_int_t ldc,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

// not yet implemented
//magma_minproduct_int_t magma_minproduct_zunmrq( magma_minproduct_side_t side, magma_minproduct_trans_t trans,
//                          magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
//                          magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
//                          magma_minproductDoubleComplex *tau,
//                          magma_minproductDoubleComplex *C, magma_minproduct_int_t ldc,
//                          magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
//                          magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zunmql(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *C, magma_minproduct_int_t ldc,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zunmqr(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *C, magma_minproduct_int_t ldc,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zunmtr(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A,    magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *C,    magma_minproduct_int_t ldc,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct function definitions / Data on CPU / Multi-GPU (alphabetical order)
*/
magma_minproduct_int_t
magma_minproduct_zgeev_m(
    magma_minproduct_vec_t jobvl, magma_minproduct_vec_t jobvr, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    #ifdef COMPLEX
    magma_minproductDoubleComplex *w,
    #else
    double *wr, double *wi,
    #endif
    magma_minproductDoubleComplex *VL, magma_minproduct_int_t ldvl,
    magma_minproductDoubleComplex *VR, magma_minproduct_int_t ldvr,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgehrd_m(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproductDoubleComplex *T,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgetrf_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zheevd_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    double *w,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zheevdx_2stage_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zheevdx_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhegst_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t itype, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhegvd_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
    double *w,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhegvdx_2stage_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhegvdx_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

#ifdef REAL
magma_minproduct_int_t
magma_minproduct_dlaex0_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t n, double *d, double *e,
    double *Q, magma_minproduct_int_t ldq,
    double *work, magma_minproduct_int_t *iwork,
    magma_minproduct_range_t range, double vl, double vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dlaex1_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t n, double *d,
    double *Q, magma_minproduct_int_t ldq,
    magma_minproduct_int_t *indxq, double rho, magma_minproduct_int_t cutpnt,
    double *work, magma_minproduct_int_t *iwork,
    magma_minproductDouble_ptr dwork[],
    magma_minproduct_queue_t queues[Magma_minproductMaxGPUs][2],
    magma_minproduct_range_t range, double vl, double vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dlaex3_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t k, magma_minproduct_int_t n, magma_minproduct_int_t n1, double *d,
    double *Q, magma_minproduct_int_t ldq, double rho,
    double *dlamda, double *Q2, magma_minproduct_int_t *indx,
    magma_minproduct_int_t *ctot, double *w, double *s, magma_minproduct_int_t *indxq,
    magma_minproductDouble_ptr dwork[],
    magma_minproduct_queue_t queues[Magma_minproductMaxGPUs][2],
    magma_minproduct_range_t range, double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *info);
#endif

magma_minproduct_int_t
magma_minproduct_zlahr2_m(
    magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *T, magma_minproduct_int_t ldt,
    magma_minproductDoubleComplex *Y, magma_minproduct_int_t ldy,
    struct zgehrd_data *data );

magma_minproduct_int_t
magma_minproduct_zlahru_m(
    magma_minproduct_int_t n, magma_minproduct_int_t ihi, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    struct zgehrd_data *data );

magma_minproduct_int_t
magma_minproduct_zpotrf_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zstedx_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_range_t range, magma_minproduct_int_t n, double vl, double vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, double *d, double *e,
    magma_minproductDoubleComplex *Z, magma_minproduct_int_t ldz,
    double *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ztrsm_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transa, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb);

magma_minproduct_int_t
magma_minproduct_zunghr_m(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *T, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zungqr_m(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *T, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zunmqr_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex *A,    magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *C,    magma_minproduct_int_t ldc,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zunmtr_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A,    magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *C,    magma_minproduct_int_t ldc,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct function definitions / Data on GPU (alphabetical order)
*/
magma_minproduct_int_t
magma_minproduct_zgegqr_gpu(
    magma_minproduct_int_t ikind, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dwork, magma_minproductDoubleComplex *work,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgelqf_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgels_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgels3_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqp3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt, magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex_ptr dwork, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqr2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dtau,
    magma_minproductDouble_ptr        dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqr2x_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dtau,
    magma_minproductDoubleComplex_ptr dT, magma_minproductDoubleComplex_ptr ddA,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqr2x2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dtau,
    magma_minproductDoubleComplex_ptr dT, magma_minproductDoubleComplex_ptr ddA,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqr2x3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dtau,
    magma_minproductDoubleComplex_ptr dT,
    magma_minproductDoubleComplex_ptr ddA,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqr2x4_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dtau,
    magma_minproductDoubleComplex_ptr dT, magma_minproductDoubleComplex_ptr ddA,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_queue_t queue,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqrf_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex_ptr dT,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqrf2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex *tau,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqrf2_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dlA[], magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex *tau,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqrf3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex_ptr dT,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqrs_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex_ptr dT,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgeqrs3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex_ptr dT,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgerbt_gpu(
    magma_minproduct_bool_t gen, magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda, 
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb, 
    magma_minproductDoubleComplex *U, magma_minproductDoubleComplex *V,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgessm_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t ib,
    magma_minproduct_int_t *ipiv,
    magma_minproductDoubleComplex_ptr dL1, magma_minproduct_int_t lddl1,
    magma_minproductDoubleComplex_ptr dL,  magma_minproduct_int_t lddl,
    magma_minproductDoubleComplex_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgesv_gpu(
    magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgesv_nopiv_gpu( 
    magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb, 
                 magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgetf2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgetrf_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgetrf_incpiv_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t ib,
    magma_minproductDoubleComplex    *hA, magma_minproduct_int_t ldha,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex    *hL, magma_minproduct_int_t ldhl,
    magma_minproductDoubleComplex_ptr dL, magma_minproduct_int_t lddl,
    magma_minproduct_int_t *ipiv,
    magma_minproductDoubleComplex_ptr dwork, magma_minproduct_int_t lddwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgetrf_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr d_lA[], magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgetrf_nopiv_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgetrf2_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t offset,
    magma_minproductDoubleComplex_ptr d_lAT[], magma_minproduct_int_t lddat, magma_minproduct_int_t *ipiv,
    magma_minproductDoubleComplex_ptr d_lAP[],
    magma_minproductDoubleComplex *W, magma_minproduct_int_t ldw,
    magma_minproduct_queue_t queues[][2],
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgetri_gpu(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductDoubleComplex_ptr dwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgetrs_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zgetrs_nopiv_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zheevd_gpu(
    magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    double *w,
    magma_minproductDoubleComplex *wA,  magma_minproduct_int_t ldwa,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zheevdx_gpu(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    double vl, double vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    magma_minproductDoubleComplex *wA,  magma_minproduct_int_t ldwa,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_minproduct_int_t
magma_minproduct_zheevr_gpu(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    double vl, double vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, double abstol, magma_minproduct_int_t *m,
    double *w,
    magma_minproductDoubleComplex_ptr dZ, magma_minproduct_int_t lddz,
    magma_minproduct_int_t *isuppz,
    magma_minproductDoubleComplex *wA, magma_minproduct_int_t ldwa,
    magma_minproductDoubleComplex *wZ, magma_minproduct_int_t ldwz,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    double *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zheevx_gpu(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    double abstol, magma_minproduct_int_t *m,
    double *w,
    magma_minproductDoubleComplex_ptr dZ, magma_minproduct_int_t lddz,
    magma_minproductDoubleComplex *wA, magma_minproduct_int_t ldwa,
    magma_minproductDoubleComplex *wZ, magma_minproduct_int_t ldwz,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    double *rwork, magma_minproduct_int_t *iwork,
    magma_minproduct_int_t *ifail,
    magma_minproduct_int_t *info);
#endif  // COMPLEX

magma_minproduct_int_t
magma_minproduct_zhegst_gpu(
    magma_minproduct_int_t itype, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhetrd_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    double *d, double *e, magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *wA,  magma_minproduct_int_t ldwa,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhetrd_he2hb_mgpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproductDoubleComplex_ptr dAmgpu[], magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dTmgpu[], magma_minproduct_int_t lddt,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t distblk,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhetrd_he2hb_mgpu_spec(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproductDoubleComplex_ptr dAmgpu[], magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dTmgpu[], magma_minproduct_int_t lddt,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t distblk,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhetrd_mgpu(
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nqueue,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    double *d, double *e, magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhetrd2_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    double *d, double *e, magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *wA,  magma_minproduct_int_t ldwa,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t lwork,
    magma_minproductDoubleComplex_ptr dwork, magma_minproduct_int_t ldwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zhetrf_nopiv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zlabrd_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex     *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    double *d, double *e, magma_minproductDoubleComplex *tauq, magma_minproductDoubleComplex *taup,
    magma_minproductDoubleComplex     *X, magma_minproduct_int_t ldx,
    magma_minproductDoubleComplex_ptr dX, magma_minproduct_int_t lddx,
    magma_minproductDoubleComplex     *Y, magma_minproduct_int_t ldy,
    magma_minproductDoubleComplex_ptr dY, magma_minproduct_int_t lddy );

magma_minproduct_int_t
magma_minproduct_zlaqps_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t offset,
    magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    magma_minproductDoubleComplex_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt, magma_minproductDoubleComplex *tau,
    double *vn1, double *vn2,
    magma_minproductDoubleComplex_ptr dauxv,
    magma_minproductDoubleComplex_ptr dF, magma_minproduct_int_t lddf);

magma_minproduct_int_t
magma_minproduct_zlaqps2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t offset,
    magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    magma_minproductDoubleComplex_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt,
    magma_minproductDoubleComplex_ptr dtau,
    magma_minproductDouble_ptr dvn1, magma_minproductDouble_ptr dvn2,
    magma_minproductDoubleComplex_ptr dauxv,
    magma_minproductDoubleComplex_ptr dF, magma_minproduct_int_t lddf);

magma_minproduct_int_t
magma_minproduct_zlaqps3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t offset,
    magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    magma_minproductDoubleComplex_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt,
    magma_minproductDoubleComplex_ptr dtau,
    magma_minproductDouble_ptr dvn1, magma_minproductDouble_ptr dvn2,
    magma_minproductDoubleComplex_ptr dauxv,
    magma_minproductDoubleComplex_ptr dF, magma_minproduct_int_t lddf);

magma_minproduct_int_t
magma_minproduct_zlarf_gpu(
    magma_minproduct_int_t m,  magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dv, magma_minproductDoubleComplex_const_ptr dtau,
    magma_minproductDoubleComplex_ptr dC,  magma_minproduct_int_t lddc);

magma_minproduct_int_t
magma_minproduct_zlarfb_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_direct_t direct, magma_minproduct_storev_t storev,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex_const_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductDoubleComplex_const_ptr dT, magma_minproduct_int_t lddt,
    magma_minproductDoubleComplex_ptr dC,       magma_minproduct_int_t lddc,
    magma_minproductDoubleComplex_ptr dwork,    magma_minproduct_int_t ldwork );

magma_minproduct_int_t
magma_minproduct_zlarfb_gpu_gemm(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_direct_t direct, magma_minproduct_storev_t storev,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex_const_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductDoubleComplex_const_ptr dT, magma_minproduct_int_t lddt,
    magma_minproductDoubleComplex_ptr dC,       magma_minproduct_int_t lddc,
    magma_minproductDoubleComplex_ptr dwork,    magma_minproduct_int_t ldwork,
    magma_minproductDoubleComplex_ptr dworkvt,  magma_minproduct_int_t ldworkvt);

magma_minproduct_int_t
magma_minproduct_zlarfb2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex_const_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductDoubleComplex_const_ptr dT, magma_minproduct_int_t lddt,
    magma_minproductDoubleComplex_ptr dC,       magma_minproduct_int_t lddc,
    magma_minproductDoubleComplex_ptr dwork,    magma_minproduct_int_t ldwork );

magma_minproduct_int_t
magma_minproduct_zlatrd_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t nb0,
    magma_minproductDoubleComplex *A,  magma_minproduct_int_t lda,
    double *e, magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex    *W,       magma_minproduct_int_t ldw,
    magma_minproductDoubleComplex_ptr dA[],    magma_minproduct_int_t ldda, magma_minproduct_int_t offset,
    magma_minproductDoubleComplex_ptr dW[],    magma_minproduct_int_t lddw,
    magma_minproductDoubleComplex    *hwork,   magma_minproduct_int_t lhwork,
    magma_minproductDoubleComplex_ptr dwork[], magma_minproduct_int_t ldwork,
    magma_minproduct_queue_t queues[] );

magma_minproduct_int_t
magma_minproduct_zlauum_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zposv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zpotf2_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zpotrf_gpu(
    magma_minproduct_uplo_t uplo,  magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zpotrf_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr d_lA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zpotrf_mgpu_right(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr d_lA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zpotrf3_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t off_i, magma_minproduct_int_t off_j, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex_ptr d_lA[],  magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr d_lP[],  magma_minproduct_int_t lddp,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda, magma_minproduct_int_t h,
    magma_minproduct_queue_t queues[][3], magma_minproduct_event_t events[][5],
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zpotri_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zpotrs_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zssssm_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m1, magma_minproduct_int_t n1,
    magma_minproduct_int_t m2, magma_minproduct_int_t n2, magma_minproduct_int_t k, magma_minproduct_int_t ib,
    magma_minproductDoubleComplex_ptr dA1, magma_minproduct_int_t ldda1,
    magma_minproductDoubleComplex_ptr dA2, magma_minproduct_int_t ldda2,
    magma_minproductDoubleComplex_ptr dL1, magma_minproduct_int_t lddl1,
    magma_minproductDoubleComplex_ptr dL2, magma_minproduct_int_t lddl2,
    magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ztrtri_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ztsqrt_gpu(
    magma_minproduct_int_t *m, magma_minproduct_int_t *n,
    magma_minproductDoubleComplex *A1, magma_minproductDoubleComplex *A2, magma_minproduct_int_t *lda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex *work, magma_minproduct_int_t *lwork,
    magma_minproductDoubleComplex_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ztstrf_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t ib, magma_minproduct_int_t nb,
    magma_minproductDoubleComplex    *hU, magma_minproduct_int_t ldhu,
    magma_minproductDoubleComplex_ptr dU, magma_minproduct_int_t lddu,
    magma_minproductDoubleComplex    *hA, magma_minproduct_int_t ldha,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex    *hL, magma_minproduct_int_t ldhl,
    magma_minproductDoubleComplex_ptr dL, magma_minproduct_int_t lddl,
    magma_minproduct_int_t *ipiv,
    magma_minproductDoubleComplex *hwork, magma_minproduct_int_t ldhwork,
    magma_minproductDoubleComplex_ptr dwork, magma_minproduct_int_t lddwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zungqr_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zunmql2_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex_ptr dC, magma_minproduct_int_t lddc,
    magma_minproductDoubleComplex *wA, magma_minproduct_int_t ldwa,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zunmqr_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex_ptr dC, magma_minproduct_int_t lddc,
    magma_minproductDoubleComplex *hwork, magma_minproduct_int_t lwork,
    magma_minproductDoubleComplex_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zunmqr2_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex_ptr dC, magma_minproduct_int_t lddc,
    magma_minproductDoubleComplex    *wA, magma_minproduct_int_t ldwa,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_zunmtr_gpu(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex *tau,
    magma_minproductDoubleComplex_ptr dC, magma_minproduct_int_t lddc,
    magma_minproductDoubleComplex    *wA, magma_minproduct_int_t ldwa,
    magma_minproduct_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct utility function definitions
*/

extern const magma_minproductDoubleComplex MAGMA_minproduct_Z_NAN;
extern const magma_minproductDoubleComplex MAGMA_minproduct_Z_INF;

magma_minproduct_int_t
magma_minproduct_znan_inf(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    const magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *cnt_nan,
    magma_minproduct_int_t *cnt_inf );

magma_minproduct_int_t
magma_minproduct_znan_inf_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *cnt_nan,
    magma_minproduct_int_t *cnt_inf );

void magma_minproduct_zprint(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const magma_minproductDoubleComplex *A, magma_minproduct_int_t lda );

void magma_minproduct_zprint_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dA, magma_minproduct_int_t ldda );

void zpanel_to_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t ib,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *work );

void zq_to_panel(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t ib,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *work );

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif /* MAGMA_minproduct_Z_H */
