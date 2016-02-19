/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
*/

#ifndef MAGMA_tally3_Z_H
#define MAGMA_tally3_Z_H

#include "magma_tally3_types.h"
#include "magma_tally3_zgehrd_m.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 Auxiliary functions to get the NB used
*/
#ifdef REAL
magma_tally3_int_t magma_tally3_get_dlaex3_m_nb();       // defined in dlaex3_m.cpp
#endif

magma_tally3_int_t magma_tally3_get_zpotrf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_zgetrf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_zgetri_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_zgeqp3_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_zgeqrf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_zgeqlf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_zgehrd_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_zhetrd_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_zhetrf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_zhetrf_nopiv_tally3_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_zgelqf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_zgebrd_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_zhegst_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_zgesvd_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_zhegst_nb_m( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_zbulge_nb( magma_tally3_int_t m, magma_tally3_int_t nbthreads );
magma_tally3_int_t magma_tally3_get_zbulge_nb_mgpu( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_zbulge_get_Vblksiz( magma_tally3_int_t m, magma_tally3_int_t nb, magma_tally3_int_t nbthreads );
magma_tally3_int_t magma_tally3_get_zbulge_gcperf();


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 function definitions / Data on CPU (alphabetical order)
*/

#ifdef REAL
// only applicable to real [sd] precisions
void
magma_tally3_dmove_eig(
    magma_tally3_range_t range, magma_tally3_int_t n, double *w,
    magma_tally3_int_t *il, magma_tally3_int_t *iu, double vl, double vu, magma_tally3_int_t *m);

// defined in dlaex3.cpp
void
magma_tally3_zvrange(
    magma_tally3_int_t k, double *d, magma_tally3_int_t *il, magma_tally3_int_t *iu, double vl, double vu);

void
magma_tally3_zirange(
    magma_tally3_int_t k, magma_tally3_int_t *indxq, magma_tally3_int_t *iil, magma_tally3_int_t *iiu, magma_tally3_int_t il, magma_tally3_int_t iu);
#endif

magma_tally3_int_t
magma_tally3_zgebrd(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    double *d, double *e,
    magma_tally3DoubleComplex *tauq, magma_tally3DoubleComplex *taup,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeev(
    magma_tally3_vec_t jobvl, magma_tally3_vec_t jobvr, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    #ifdef COMPLEX
    magma_tally3DoubleComplex *w,
    #else
    double *wr, double *wi,
    #endif
    magma_tally3DoubleComplex *VL, magma_tally3_int_t ldvl,
    magma_tally3DoubleComplex *VR, magma_tally3_int_t ldvr,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgehrd(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3DoubleComplex_ptr dT,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgehrd2(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgelqf(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A,    magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqlf(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A,    magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqp3(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *jpvt, magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqrf(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqrf_ooc(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqrf4(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A,    magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgesdd(
    magma_tally3_vec_t jobz, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    double *s,
    magma_tally3DoubleComplex *U, magma_tally3_int_t ldu,
    magma_tally3DoubleComplex *VT, magma_tally3_int_t ldvt,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *iwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgesv(
    magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3DoubleComplex *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgesv_rbt(
    magma_tally3_bool_t ref, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda, 
    magma_tally3DoubleComplex *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgesvd(
    magma_tally3_vec_t jobu, magma_tally3_vec_t jobvt, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A,    magma_tally3_int_t lda, double *s,
    magma_tally3DoubleComplex *U,    magma_tally3_int_t ldu,
    magma_tally3DoubleComplex *VT,   magma_tally3_int_t ldvt,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgetf2_nopiv(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgetrf(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgetrf_nopiv(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgetrf_piv(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t NB,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgetrf2(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zheevd(
    magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    double *w,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zheevdx(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zheevdx_2stage(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally3_int_t
magma_tally3_zheevr(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    double vl, double vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, double abstol, magma_tally3_int_t *m,
    double *w,
    magma_tally3DoubleComplex *Z, magma_tally3_int_t ldz,
    magma_tally3_int_t *isuppz,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    double *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

// no real [sd] precisions available
magma_tally3_int_t
magma_tally3_zheevx(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    double vl, double vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, double abstol, magma_tally3_int_t *m,
    double *w,
    magma_tally3DoubleComplex *Z, magma_tally3_int_t ldz,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    double *rwork, magma_tally3_int_t *iwork,
    magma_tally3_int_t *ifail,
    magma_tally3_int_t *info);
#endif  // COMPLEX

magma_tally3_int_t
magma_tally3_zhegst(
    magma_tally3_int_t itype, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhegvd(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *B, magma_tally3_int_t ldb,
    double *w, magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhegvdx(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n, magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *B, magma_tally3_int_t ldb,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhegvdx_2stage(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *B, magma_tally3_int_t ldb,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhegvr(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *B, magma_tally3_int_t ldb,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    double abstol, magma_tally3_int_t *m, double *w,
    magma_tally3DoubleComplex *Z, magma_tally3_int_t ldz,
    magma_tally3_int_t *isuppz, magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    double *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhegvx(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n, magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *B, magma_tally3_int_t ldb,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    double abstol, magma_tally3_int_t *m, double *w,
    magma_tally3DoubleComplex *Z, magma_tally3_int_t ldz,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork, double *rwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t *ifail,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhesv(magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
            magma_tally3DoubleComplex *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
            magma_tally3DoubleComplex *B, magma_tally3_int_t ldb,
            magma_tally3_int_t *info );

magma_tally3_int_t
magma_tally3_zhetrd(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    double *d, double *e, magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhetrf(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhetrf_nopiv_tally3(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhetrd_hb2st(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t Vblksiz,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    double *d, double *e,
    magma_tally3DoubleComplex *V, magma_tally3_int_t ldv,
    magma_tally3DoubleComplex *TAU, magma_tally3_int_t compT,
    magma_tally3DoubleComplex *T, magma_tally3_int_t ldt);

magma_tally3_int_t
magma_tally3_zhetrd_he2hb(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3DoubleComplex_ptr dT,
    magma_tally3_int_t *info);

#ifdef REAL
// only applicable to real [sd] precisions
magma_tally3_int_t
magma_tally3_dlaex0(
    magma_tally3_int_t n, double *d, double *e,
    double *Q, magma_tally3_int_t ldq,
    double *work, magma_tally3_int_t *iwork,
    magma_tally3Double_ptr dwork,
    magma_tally3_range_t range, double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dlaex1(
    magma_tally3_int_t n, double *d,
    double *Q, magma_tally3_int_t ldq,
    magma_tally3_int_t *indxq, double rho, magma_tally3_int_t cutpnt,
    double *work, magma_tally3_int_t *iwork,
    magma_tally3Double_ptr dwork,
    magma_tally3_range_t range, double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dlaex3(
    magma_tally3_int_t k, magma_tally3_int_t n, magma_tally3_int_t n1, double *d,
    double *Q, magma_tally3_int_t ldq,
    double rho,
    double *dlamda, double *Q2, magma_tally3_int_t *indx,
    magma_tally3_int_t *ctot, double *w, double *s, magma_tally3_int_t *indxq,
    magma_tally3Double_ptr dwork,
    magma_tally3_range_t range, double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *info);
#endif  // REAL

magma_tally3_int_t
magma_tally3_zlahef_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3DoubleComplex    *hA, magma_tally3_int_t lda,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3DoubleComplex_ptr dW, magma_tally3_int_t lddw,
    magma_tally3_queue_t queues[], magma_tally3_event_t event[],
    magma_tally3_int_t *info);

magma_tally3_int_t
zhetrf_nopiv_tally3_cpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t ib,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhetrs_nopiv_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhesv_nopiv_gpu(
    magma_tally3_uplo_t uplo,  magma_tally3_int_t n, magma_tally3_int_t nrhs, 
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb, 
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zlahr2(
    magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3_int_t nb,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dV, magma_tally3_int_t lddv,
    magma_tally3DoubleComplex *A,  magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *T,  magma_tally3_int_t ldt,
    magma_tally3DoubleComplex *Y,  magma_tally3_int_t ldy);

magma_tally3_int_t
magma_tally3_zlahru(
    magma_tally3_int_t n, magma_tally3_int_t ihi, magma_tally3_int_t k, magma_tally3_int_t nb,
    magma_tally3DoubleComplex     *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dY, magma_tally3_int_t lddy,
    magma_tally3DoubleComplex_ptr dV, magma_tally3_int_t lddv,
    magma_tally3DoubleComplex_ptr dT,
    magma_tally3DoubleComplex_ptr dwork);

#ifdef REAL
magma_tally3_int_t
magma_tally3_dlaln2(
    magma_tally3_int_t trans, magma_tally3_int_t na, magma_tally3_int_t nw,
    double smin, double ca, const double *A, magma_tally3_int_t lda,
    double d1, double d2,   const double *B, magma_tally3_int_t ldb,
    double wr, double wi, double *X, magma_tally3_int_t ldx,
    double *scale, double *xnorm,
    magma_tally3_int_t *info);
#endif

magma_tally3_int_t
magma_tally3_zlaqps(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3DoubleComplex *A,  magma_tally3_int_t lda,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt, magma_tally3DoubleComplex *tau, double *vn1, double *vn2,
    magma_tally3DoubleComplex *auxv,
    magma_tally3DoubleComplex *F,  magma_tally3_int_t ldf,
    magma_tally3DoubleComplex_ptr dF, magma_tally3_int_t lddf );

#ifdef REAL
magma_tally3_int_t
magma_tally3_zlaqtrsd(
    magma_tally3_trans_t trans, magma_tally3_int_t n,
    const double *T, magma_tally3_int_t ldt,
    double *x,       magma_tally3_int_t ldx,
    const double *cnorm,
    magma_tally3_int_t *info);
#endif

magma_tally3_int_t
magma_tally3_zlatrd(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    double *e, magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *W, magma_tally3_int_t ldw,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dW, magma_tally3_int_t lddw);

magma_tally3_int_t
magma_tally3_zlatrd2(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3DoubleComplex *A,  magma_tally3_int_t lda,
    double *e, magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *W,  magma_tally3_int_t ldw,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dW, magma_tally3_int_t lddw,
    magma_tally3DoubleComplex_ptr dwork, magma_tally3_int_t ldwork);

#ifdef COMPLEX
magma_tally3_int_t
magma_tally3_zlatrsd(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_diag_t diag, magma_tally3_bool_t normin,
    magma_tally3_int_t n, const magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex lambda,
    magma_tally3DoubleComplex *x,
    double *scale, double *cnorm,
    magma_tally3_int_t *info);
#endif

magma_tally3_int_t
magma_tally3_zlauum(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zposv(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zpotrf(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zpotri(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zstedx(
    magma_tally3_range_t range, magma_tally3_int_t n, double vl, double vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, double *d, double *e,
    magma_tally3DoubleComplex *Z, magma_tally3_int_t ldz,
    double *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3Double_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ztrevc3(
    magma_tally3_side_t side, magma_tally3_vec_t howmany,
    magma_tally3_int_t *select, magma_tally3_int_t n,
    magma_tally3DoubleComplex *T,  magma_tally3_int_t ldt,
    magma_tally3DoubleComplex *VL, magma_tally3_int_t ldvl,
    magma_tally3DoubleComplex *VR, magma_tally3_int_t ldvr,
    magma_tally3_int_t mm, magma_tally3_int_t *mout,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ztrevc3_mt(
    magma_tally3_side_t side, magma_tally3_vec_t howmany,
    magma_tally3_int_t *select, magma_tally3_int_t n,
    magma_tally3DoubleComplex *T,  magma_tally3_int_t ldt,
    magma_tally3DoubleComplex *VL, magma_tally3_int_t ldvl,
    magma_tally3DoubleComplex *VR, magma_tally3_int_t ldvr,
    magma_tally3_int_t mm, magma_tally3_int_t *mout,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ztrtri(
    magma_tally3_uplo_t uplo, magma_tally3_diag_t diag, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zunghr(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zungqr(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zungqr2(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zunmbr(
    magma_tally3_vect_t vect, magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *C, magma_tally3_int_t ldc,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zunmlq(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *C, magma_tally3_int_t ldc,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

// not yet implemented
//magma_tally3_int_t magma_tally3_zunmrq( magma_tally3_side_t side, magma_tally3_trans_t trans,
//                          magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
//                          magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
//                          magma_tally3DoubleComplex *tau,
//                          magma_tally3DoubleComplex *C, magma_tally3_int_t ldc,
//                          magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
//                          magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zunmql(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *C, magma_tally3_int_t ldc,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zunmqr(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *C, magma_tally3_int_t ldc,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zunmtr(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A,    magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *C,    magma_tally3_int_t ldc,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 function definitions / Data on CPU / Multi-GPU (alphabetical order)
*/
magma_tally3_int_t
magma_tally3_zgeev_m(
    magma_tally3_vec_t jobvl, magma_tally3_vec_t jobvr, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    #ifdef COMPLEX
    magma_tally3DoubleComplex *w,
    #else
    double *wr, double *wi,
    #endif
    magma_tally3DoubleComplex *VL, magma_tally3_int_t ldvl,
    magma_tally3DoubleComplex *VR, magma_tally3_int_t ldvr,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgehrd_m(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3DoubleComplex *T,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgetrf_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zheevd_m(
    magma_tally3_int_t ngpu,
    magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    double *w,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zheevdx_2stage_m(
    magma_tally3_int_t ngpu,
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zheevdx_m(
    magma_tally3_int_t ngpu,
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhegst_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhegvd_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *B, magma_tally3_int_t ldb,
    double *w,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhegvdx_2stage_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *B, magma_tally3_int_t ldb,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhegvdx_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *B, magma_tally3_int_t ldb,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

#ifdef REAL
magma_tally3_int_t
magma_tally3_dlaex0_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t n, double *d, double *e,
    double *Q, magma_tally3_int_t ldq,
    double *work, magma_tally3_int_t *iwork,
    magma_tally3_range_t range, double vl, double vu,
    magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dlaex1_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t n, double *d,
    double *Q, magma_tally3_int_t ldq,
    magma_tally3_int_t *indxq, double rho, magma_tally3_int_t cutpnt,
    double *work, magma_tally3_int_t *iwork,
    magma_tally3Double_ptr dwork[],
    magma_tally3_queue_t queues[Magma_tally3MaxGPUs][2],
    magma_tally3_range_t range, double vl, double vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dlaex3_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t k, magma_tally3_int_t n, magma_tally3_int_t n1, double *d,
    double *Q, magma_tally3_int_t ldq, double rho,
    double *dlamda, double *Q2, magma_tally3_int_t *indx,
    magma_tally3_int_t *ctot, double *w, double *s, magma_tally3_int_t *indxq,
    magma_tally3Double_ptr dwork[],
    magma_tally3_queue_t queues[Magma_tally3MaxGPUs][2],
    magma_tally3_range_t range, double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *info);
#endif

magma_tally3_int_t
magma_tally3_zlahr2_m(
    magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3_int_t nb,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *T, magma_tally3_int_t ldt,
    magma_tally3DoubleComplex *Y, magma_tally3_int_t ldy,
    struct zgehrd_data_tally3 *data );

magma_tally3_int_t
magma_tally3_zlahru_m(
    magma_tally3_int_t n, magma_tally3_int_t ihi, magma_tally3_int_t k, magma_tally3_int_t nb,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    struct zgehrd_data_tally3 *data );

magma_tally3_int_t
magma_tally3_zpotrf_m(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zstedx_m(
    magma_tally3_int_t ngpu,
    magma_tally3_range_t range, magma_tally3_int_t n, double vl, double vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, double *d, double *e,
    magma_tally3DoubleComplex *Z, magma_tally3_int_t ldz,
    double *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ztrsm_m(
    magma_tally3_int_t ngpu,
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transa, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3DoubleComplex alpha,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *B, magma_tally3_int_t ldb);

magma_tally3_int_t
magma_tally3_zunghr_m(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *T, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zungqr_m(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *T, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zunmqr_m(
    magma_tally3_int_t ngpu,
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex *A,    magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *C,    magma_tally3_int_t ldc,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zunmtr_m(
    magma_tally3_int_t ngpu,
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A,    magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *C,    magma_tally3_int_t ldc,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 function definitions / Data on GPU (alphabetical order)
*/
magma_tally3_int_t
magma_tally3_zgegqr_gpu(
    magma_tally3_int_t ikind, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dwork, magma_tally3DoubleComplex *work,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgelqf_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgels_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3DoubleComplex *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgels3_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3DoubleComplex *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqp3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt, magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex_ptr dwork, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqr2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dtau,
    magma_tally3Double_ptr        dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqr2x_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dtau,
    magma_tally3DoubleComplex_ptr dT, magma_tally3DoubleComplex_ptr ddA,
    magma_tally3Double_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqr2x2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dtau,
    magma_tally3DoubleComplex_ptr dT, magma_tally3DoubleComplex_ptr ddA,
    magma_tally3Double_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqr2x3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dtau,
    magma_tally3DoubleComplex_ptr dT,
    magma_tally3DoubleComplex_ptr ddA,
    magma_tally3Double_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqr2x4_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dtau,
    magma_tally3DoubleComplex_ptr dT, magma_tally3DoubleComplex_ptr ddA,
    magma_tally3Double_ptr dwork,
    magma_tally3_queue_t queue,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqrf_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex_ptr dT,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqrf2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex *tau,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqrf2_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dlA[], magma_tally3_int_t ldda,
    magma_tally3DoubleComplex *tau,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqrf3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex_ptr dT,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqrs_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex_ptr dT,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3DoubleComplex *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgeqrs3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex_ptr dT,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3DoubleComplex *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgerbt_gpu(
    magma_tally3_bool_t gen, magma_tally3_int_t n, magma_tally3_int_t nrhs, 
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda, 
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb, 
    magma_tally3DoubleComplex *U, magma_tally3DoubleComplex *V,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgessm_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3_int_t ib,
    magma_tally3_int_t *ipiv,
    magma_tally3DoubleComplex_ptr dL1, magma_tally3_int_t lddl1,
    magma_tally3DoubleComplex_ptr dL,  magma_tally3_int_t lddl,
    magma_tally3DoubleComplex_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgesv_gpu(
    magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgesv_nopiv_gpu( 
    magma_tally3_int_t n, magma_tally3_int_t nrhs, 
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb, 
                 magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgetf2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgetrf_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgetrf_incpiv_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t ib,
    magma_tally3DoubleComplex    *hA, magma_tally3_int_t ldha,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex    *hL, magma_tally3_int_t ldhl,
    magma_tally3DoubleComplex_ptr dL, magma_tally3_int_t lddl,
    magma_tally3_int_t *ipiv,
    magma_tally3DoubleComplex_ptr dwork, magma_tally3_int_t lddwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgetrf_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr d_lA[], magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgetrf_nopiv_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgetrf2_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t offset,
    magma_tally3DoubleComplex_ptr d_lAT[], magma_tally3_int_t lddat, magma_tally3_int_t *ipiv,
    magma_tally3DoubleComplex_ptr d_lAP[],
    magma_tally3DoubleComplex *W, magma_tally3_int_t ldw,
    magma_tally3_queue_t queues[][2],
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgetri_gpu(
    magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3DoubleComplex_ptr dwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgetrs_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zgetrs_nopiv_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zheevd_gpu(
    magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    double *w,
    magma_tally3DoubleComplex *wA,  magma_tally3_int_t ldwa,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zheevdx_gpu(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    double vl, double vu,
    magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    magma_tally3DoubleComplex *wA,  magma_tally3_int_t ldwa,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally3_int_t
magma_tally3_zheevr_gpu(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    double vl, double vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, double abstol, magma_tally3_int_t *m,
    double *w,
    magma_tally3DoubleComplex_ptr dZ, magma_tally3_int_t lddz,
    magma_tally3_int_t *isuppz,
    magma_tally3DoubleComplex *wA, magma_tally3_int_t ldwa,
    magma_tally3DoubleComplex *wZ, magma_tally3_int_t ldwz,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    double *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zheevx_gpu(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    double abstol, magma_tally3_int_t *m,
    double *w,
    magma_tally3DoubleComplex_ptr dZ, magma_tally3_int_t lddz,
    magma_tally3DoubleComplex *wA, magma_tally3_int_t ldwa,
    magma_tally3DoubleComplex *wZ, magma_tally3_int_t ldwz,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    double *rwork, magma_tally3_int_t *iwork,
    magma_tally3_int_t *ifail,
    magma_tally3_int_t *info);
#endif  // COMPLEX

magma_tally3_int_t
magma_tally3_zhegst_gpu(
    magma_tally3_int_t itype, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhetrd_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    double *d, double *e, magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *wA,  magma_tally3_int_t ldwa,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhetrd_he2hb_mgpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3DoubleComplex_ptr dAmgpu[], magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dTmgpu[], magma_tally3_int_t lddt,
    magma_tally3_int_t ngpu, magma_tally3_int_t distblk,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhetrd_he2hb_mgpu_spec(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3DoubleComplex_ptr dAmgpu[], magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dTmgpu[], magma_tally3_int_t lddt,
    magma_tally3_int_t ngpu, magma_tally3_int_t distblk,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhetrd_mgpu(
    magma_tally3_int_t ngpu, magma_tally3_int_t nqueue,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    double *d, double *e, magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhetrd2_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    double *d, double *e, magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *wA,  magma_tally3_int_t ldwa,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3DoubleComplex_ptr dwork, magma_tally3_int_t ldwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zhetrf_nopiv_tally3_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zlabrd_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3DoubleComplex     *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    double *d, double *e, magma_tally3DoubleComplex *tauq, magma_tally3DoubleComplex *taup,
    magma_tally3DoubleComplex     *X, magma_tally3_int_t ldx,
    magma_tally3DoubleComplex_ptr dX, magma_tally3_int_t lddx,
    magma_tally3DoubleComplex     *Y, magma_tally3_int_t ldy,
    magma_tally3DoubleComplex_ptr dY, magma_tally3_int_t lddy );

magma_tally3_int_t
magma_tally3_zlaqps_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3DoubleComplex_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt, magma_tally3DoubleComplex *tau,
    double *vn1, double *vn2,
    magma_tally3DoubleComplex_ptr dauxv,
    magma_tally3DoubleComplex_ptr dF, magma_tally3_int_t lddf);

magma_tally3_int_t
magma_tally3_zlaqps2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3DoubleComplex_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt,
    magma_tally3DoubleComplex_ptr dtau,
    magma_tally3Double_ptr dvn1, magma_tally3Double_ptr dvn2,
    magma_tally3DoubleComplex_ptr dauxv,
    magma_tally3DoubleComplex_ptr dF, magma_tally3_int_t lddf);

magma_tally3_int_t
magma_tally3_zlaqps3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3DoubleComplex_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt,
    magma_tally3DoubleComplex_ptr dtau,
    magma_tally3Double_ptr dvn1, magma_tally3Double_ptr dvn2,
    magma_tally3DoubleComplex_ptr dauxv,
    magma_tally3DoubleComplex_ptr dF, magma_tally3_int_t lddf);

magma_tally3_int_t
magma_tally3_zlarf_gpu(
    magma_tally3_int_t m,  magma_tally3_int_t n,
    magma_tally3DoubleComplex_const_ptr dv, magma_tally3DoubleComplex_const_ptr dtau,
    magma_tally3DoubleComplex_ptr dC,  magma_tally3_int_t lddc);

magma_tally3_int_t
magma_tally3_zlarfb_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans, magma_tally3_direct_t direct, magma_tally3_storev_t storev,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex_const_ptr dV, magma_tally3_int_t lddv,
    magma_tally3DoubleComplex_const_ptr dT, magma_tally3_int_t lddt,
    magma_tally3DoubleComplex_ptr dC,       magma_tally3_int_t lddc,
    magma_tally3DoubleComplex_ptr dwork,    magma_tally3_int_t ldwork );

magma_tally3_int_t
magma_tally3_zlarfb_gpu_gemm(
    magma_tally3_side_t side, magma_tally3_trans_t trans, magma_tally3_direct_t direct, magma_tally3_storev_t storev,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex_const_ptr dV, magma_tally3_int_t lddv,
    magma_tally3DoubleComplex_const_ptr dT, magma_tally3_int_t lddt,
    magma_tally3DoubleComplex_ptr dC,       magma_tally3_int_t lddc,
    magma_tally3DoubleComplex_ptr dwork,    magma_tally3_int_t ldwork,
    magma_tally3DoubleComplex_ptr dworkvt,  magma_tally3_int_t ldworkvt);

magma_tally3_int_t
magma_tally3_zlarfb2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex_const_ptr dV, magma_tally3_int_t lddv,
    magma_tally3DoubleComplex_const_ptr dT, magma_tally3_int_t lddt,
    magma_tally3DoubleComplex_ptr dC,       magma_tally3_int_t lddc,
    magma_tally3DoubleComplex_ptr dwork,    magma_tally3_int_t ldwork );

magma_tally3_int_t
magma_tally3_zlatrd_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t nb0,
    magma_tally3DoubleComplex *A,  magma_tally3_int_t lda,
    double *e, magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex    *W,       magma_tally3_int_t ldw,
    magma_tally3DoubleComplex_ptr dA[],    magma_tally3_int_t ldda, magma_tally3_int_t offset,
    magma_tally3DoubleComplex_ptr dW[],    magma_tally3_int_t lddw,
    magma_tally3DoubleComplex    *hwork,   magma_tally3_int_t lhwork,
    magma_tally3DoubleComplex_ptr dwork[], magma_tally3_int_t ldwork,
    magma_tally3_queue_t queues[] );

magma_tally3_int_t
magma_tally3_zlauum_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zposv_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zpotf2_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zpotrf_gpu(
    magma_tally3_uplo_t uplo,  magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zpotrf_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr d_lA[], magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zpotrf_mgpu_right(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr d_lA[], magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zpotrf3_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t off_i, magma_tally3_int_t off_j, magma_tally3_int_t nb,
    magma_tally3DoubleComplex_ptr d_lA[],  magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr d_lP[],  magma_tally3_int_t lddp,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda, magma_tally3_int_t h,
    magma_tally3_queue_t queues[][3], magma_tally3_event_t events[][5],
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zpotri_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zpotrs_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zssssm_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m1, magma_tally3_int_t n1,
    magma_tally3_int_t m2, magma_tally3_int_t n2, magma_tally3_int_t k, magma_tally3_int_t ib,
    magma_tally3DoubleComplex_ptr dA1, magma_tally3_int_t ldda1,
    magma_tally3DoubleComplex_ptr dA2, magma_tally3_int_t ldda2,
    magma_tally3DoubleComplex_ptr dL1, magma_tally3_int_t lddl1,
    magma_tally3DoubleComplex_ptr dL2, magma_tally3_int_t lddl2,
    magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ztrtri_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_diag_t diag, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ztsqrt_gpu(
    magma_tally3_int_t *m, magma_tally3_int_t *n,
    magma_tally3DoubleComplex *A1, magma_tally3DoubleComplex *A2, magma_tally3_int_t *lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *work, magma_tally3_int_t *lwork,
    magma_tally3DoubleComplex_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ztstrf_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t ib, magma_tally3_int_t nb,
    magma_tally3DoubleComplex    *hU, magma_tally3_int_t ldhu,
    magma_tally3DoubleComplex_ptr dU, magma_tally3_int_t lddu,
    magma_tally3DoubleComplex    *hA, magma_tally3_int_t ldha,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex    *hL, magma_tally3_int_t ldhl,
    magma_tally3DoubleComplex_ptr dL, magma_tally3_int_t lddl,
    magma_tally3_int_t *ipiv,
    magma_tally3DoubleComplex *hwork, magma_tally3_int_t ldhwork,
    magma_tally3DoubleComplex_ptr dwork, magma_tally3_int_t lddwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zungqr_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zunmql2_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex_ptr dC, magma_tally3_int_t lddc,
    magma_tally3DoubleComplex *wA, magma_tally3_int_t ldwa,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zunmqr_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex_ptr dC, magma_tally3_int_t lddc,
    magma_tally3DoubleComplex *hwork, magma_tally3_int_t lwork,
    magma_tally3DoubleComplex_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zunmqr2_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex_ptr dC, magma_tally3_int_t lddc,
    magma_tally3DoubleComplex    *wA, magma_tally3_int_t ldwa,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_zunmtr_gpu(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex_ptr dC, magma_tally3_int_t lddc,
    magma_tally3DoubleComplex    *wA, magma_tally3_int_t ldwa,
    magma_tally3_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 utility function definitions
*/

extern const magma_tally3DoubleComplex MAGMA_tally3_Z_NAN;
extern const magma_tally3DoubleComplex MAGMA_tally3_Z_INF;

magma_tally3_int_t
magma_tally3_znan_inf(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *cnt_nan,
    magma_tally3_int_t *cnt_inf );

magma_tally3_int_t
magma_tally3_znan_inf_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *cnt_nan,
    magma_tally3_int_t *cnt_inf );

void magma_tally3_zprint(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3DoubleComplex *A, magma_tally3_int_t lda );

void magma_tally3_zprint_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_const_ptr dA, magma_tally3_int_t ldda );

void zpanel_to_q_tally3(
    magma_tally3_uplo_t uplo, magma_tally3_int_t ib,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *work );

void zq_to_panel_tally3(
    magma_tally3_uplo_t uplo, magma_tally3_int_t ib,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *work );

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif /* MAGMA_tally3_Z_H */
