/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
*/

#ifndef MAGMA_tally2_Z_H
#define MAGMA_tally2_Z_H

#include "magma_tally2_types.h"
#include "magma_tally2_zgehrd_m.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 Auxiliary functions to get the NB used
*/
#ifdef REAL
magma_tally2_int_t magma_tally2_get_dlaex3_m_nb();       // defined in dlaex3_m.cpp
#endif

magma_tally2_int_t magma_tally2_get_zpotrf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_zgetrf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_zgetri_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_zgeqp3_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_zgeqrf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_zgeqlf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_zgehrd_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_zhetrd_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_zhetrf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_zhetrf_nopiv_tally2_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_zgelqf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_zgebrd_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_zhegst_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_zgesvd_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_zhegst_nb_m( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_zbulge_nb( magma_tally2_int_t m, magma_tally2_int_t nbthreads );
magma_tally2_int_t magma_tally2_get_zbulge_nb_mgpu( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_zbulge_get_Vblksiz( magma_tally2_int_t m, magma_tally2_int_t nb, magma_tally2_int_t nbthreads );
magma_tally2_int_t magma_tally2_get_zbulge_gcperf();


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 function definitions / Data on CPU (alphabetical order)
*/

#ifdef REAL
// only applicable to real [sd] precisions
void
magma_tally2_dmove_eig(
    magma_tally2_range_t range, magma_tally2_int_t n, double *w,
    magma_tally2_int_t *il, magma_tally2_int_t *iu, double vl, double vu, magma_tally2_int_t *m);

// defined in dlaex3.cpp
void
magma_tally2_zvrange(
    magma_tally2_int_t k, double *d, magma_tally2_int_t *il, magma_tally2_int_t *iu, double vl, double vu);

void
magma_tally2_zirange(
    magma_tally2_int_t k, magma_tally2_int_t *indxq, magma_tally2_int_t *iil, magma_tally2_int_t *iiu, magma_tally2_int_t il, magma_tally2_int_t iu);
#endif

magma_tally2_int_t
magma_tally2_zgebrd(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    double *d, double *e,
    magma_tally2DoubleComplex *tauq, magma_tally2DoubleComplex *taup,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeev(
    magma_tally2_vec_t jobvl, magma_tally2_vec_t jobvr, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    #ifdef COMPLEX
    magma_tally2DoubleComplex *w,
    #else
    double *wr, double *wi,
    #endif
    magma_tally2DoubleComplex *VL, magma_tally2_int_t ldvl,
    magma_tally2DoubleComplex *VR, magma_tally2_int_t ldvr,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgehrd(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2DoubleComplex_ptr dT,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgehrd2(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgelqf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A,    magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqlf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A,    magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqp3(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *jpvt, magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqrf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqrf_ooc(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqrf4(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A,    magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgesdd(
    magma_tally2_vec_t jobz, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    double *s,
    magma_tally2DoubleComplex *U, magma_tally2_int_t ldu,
    magma_tally2DoubleComplex *VT, magma_tally2_int_t ldvt,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally2_int_t *iwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgesv(
    magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgesv_rbt(
    magma_tally2_bool_t ref, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda, 
    magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgesvd(
    magma_tally2_vec_t jobu, magma_tally2_vec_t jobvt, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A,    magma_tally2_int_t lda, double *s,
    magma_tally2DoubleComplex *U,    magma_tally2_int_t ldu,
    magma_tally2DoubleComplex *VT,   magma_tally2_int_t ldvt,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgetf2_nopiv(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgetrf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgetrf_nopiv(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgetrf_piv(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t NB,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgetrf2(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zheevd(
    magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    double *w,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zheevdx(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zheevdx_2stage(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally2_int_t
magma_tally2_zheevr(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    double vl, double vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, double abstol, magma_tally2_int_t *m,
    double *w,
    magma_tally2DoubleComplex *Z, magma_tally2_int_t ldz,
    magma_tally2_int_t *isuppz,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    double *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

// no real [sd] precisions available
magma_tally2_int_t
magma_tally2_zheevx(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    double vl, double vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, double abstol, magma_tally2_int_t *m,
    double *w,
    magma_tally2DoubleComplex *Z, magma_tally2_int_t ldz,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    double *rwork, magma_tally2_int_t *iwork,
    magma_tally2_int_t *ifail,
    magma_tally2_int_t *info);
#endif  // COMPLEX

magma_tally2_int_t
magma_tally2_zhegst(
    magma_tally2_int_t itype, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhegvd(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
    double *w, magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhegvdx(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n, magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhegvdx_2stage(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhegvr(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    double abstol, magma_tally2_int_t *m, double *w,
    magma_tally2DoubleComplex *Z, magma_tally2_int_t ldz,
    magma_tally2_int_t *isuppz, magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    double *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhegvx(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n, magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    double abstol, magma_tally2_int_t *m, double *w,
    magma_tally2DoubleComplex *Z, magma_tally2_int_t ldz,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork, double *rwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t *ifail,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhesv(magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
            magma_tally2DoubleComplex *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
            magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
            magma_tally2_int_t *info );

magma_tally2_int_t
magma_tally2_zhetrd(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    double *d, double *e, magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhetrf(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhetrf_nopiv_tally2(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhetrd_hb2st(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t Vblksiz,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    double *d, double *e,
    magma_tally2DoubleComplex *V, magma_tally2_int_t ldv,
    magma_tally2DoubleComplex *TAU, magma_tally2_int_t compT,
    magma_tally2DoubleComplex *T, magma_tally2_int_t ldt);

magma_tally2_int_t
magma_tally2_zhetrd_he2hb(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2DoubleComplex_ptr dT,
    magma_tally2_int_t *info);

#ifdef REAL
// only applicable to real [sd] precisions
magma_tally2_int_t
magma_tally2_dlaex0(
    magma_tally2_int_t n, double *d, double *e,
    double *Q, magma_tally2_int_t ldq,
    double *work, magma_tally2_int_t *iwork,
    magma_tally2Double_ptr dwork,
    magma_tally2_range_t range, double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dlaex1(
    magma_tally2_int_t n, double *d,
    double *Q, magma_tally2_int_t ldq,
    magma_tally2_int_t *indxq, double rho, magma_tally2_int_t cutpnt,
    double *work, magma_tally2_int_t *iwork,
    magma_tally2Double_ptr dwork,
    magma_tally2_range_t range, double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dlaex3(
    magma_tally2_int_t k, magma_tally2_int_t n, magma_tally2_int_t n1, double *d,
    double *Q, magma_tally2_int_t ldq,
    double rho,
    double *dlamda, double *Q2, magma_tally2_int_t *indx,
    magma_tally2_int_t *ctot, double *w, double *s, magma_tally2_int_t *indxq,
    magma_tally2Double_ptr dwork,
    magma_tally2_range_t range, double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *info);
#endif  // REAL

magma_tally2_int_t
magma_tally2_zlahef_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t *kb,
    magma_tally2DoubleComplex    *hA, magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2DoubleComplex_ptr dW, magma_tally2_int_t lddw,
    magma_tally2_queue_t queues[], magma_tally2_event_t event[],
    magma_tally2_int_t *info);

magma_tally2_int_t
zhetrf_nopiv_tally2_cpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t ib,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhetrs_nopiv_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhesv_nopiv_gpu(
    magma_tally2_uplo_t uplo,  magma_tally2_int_t n, magma_tally2_int_t nrhs, 
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb, 
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zlahr2(
    magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2_int_t nb,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dV, magma_tally2_int_t lddv,
    magma_tally2DoubleComplex *A,  magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *T,  magma_tally2_int_t ldt,
    magma_tally2DoubleComplex *Y,  magma_tally2_int_t ldy);

magma_tally2_int_t
magma_tally2_zlahru(
    magma_tally2_int_t n, magma_tally2_int_t ihi, magma_tally2_int_t k, magma_tally2_int_t nb,
    magma_tally2DoubleComplex     *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dY, magma_tally2_int_t lddy,
    magma_tally2DoubleComplex_ptr dV, magma_tally2_int_t lddv,
    magma_tally2DoubleComplex_ptr dT,
    magma_tally2DoubleComplex_ptr dwork);

#ifdef REAL
magma_tally2_int_t
magma_tally2_dlaln2(
    magma_tally2_int_t trans, magma_tally2_int_t na, magma_tally2_int_t nw,
    double smin, double ca, const double *A, magma_tally2_int_t lda,
    double d1, double d2,   const double *B, magma_tally2_int_t ldb,
    double wr, double wi, double *X, magma_tally2_int_t ldx,
    double *scale, double *xnorm,
    magma_tally2_int_t *info);
#endif

magma_tally2_int_t
magma_tally2_zlaqps(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t offset,
    magma_tally2_int_t nb, magma_tally2_int_t *kb,
    magma_tally2DoubleComplex *A,  magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt, magma_tally2DoubleComplex *tau, double *vn1, double *vn2,
    magma_tally2DoubleComplex *auxv,
    magma_tally2DoubleComplex *F,  magma_tally2_int_t ldf,
    magma_tally2DoubleComplex_ptr dF, magma_tally2_int_t lddf );

#ifdef REAL
magma_tally2_int_t
magma_tally2_zlaqtrsd(
    magma_tally2_trans_t trans, magma_tally2_int_t n,
    const double *T, magma_tally2_int_t ldt,
    double *x,       magma_tally2_int_t ldx,
    const double *cnorm,
    magma_tally2_int_t *info);
#endif

magma_tally2_int_t
magma_tally2_zlatrd(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    double *e, magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *W, magma_tally2_int_t ldw,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dW, magma_tally2_int_t lddw);

magma_tally2_int_t
magma_tally2_zlatrd2(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2DoubleComplex *A,  magma_tally2_int_t lda,
    double *e, magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *W,  magma_tally2_int_t ldw,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dW, magma_tally2_int_t lddw,
    magma_tally2DoubleComplex_ptr dwork, magma_tally2_int_t ldwork);

#ifdef COMPLEX
magma_tally2_int_t
magma_tally2_zlatrsd(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_diag_t diag, magma_tally2_bool_t normin,
    magma_tally2_int_t n, const magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex lambda,
    magma_tally2DoubleComplex *x,
    double *scale, double *cnorm,
    magma_tally2_int_t *info);
#endif

magma_tally2_int_t
magma_tally2_zlauum(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zposv(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zpotrf(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zpotri(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zstedx(
    magma_tally2_range_t range, magma_tally2_int_t n, double vl, double vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, double *d, double *e,
    magma_tally2DoubleComplex *Z, magma_tally2_int_t ldz,
    double *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2Double_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ztrevc3(
    magma_tally2_side_t side, magma_tally2_vec_t howmany,
    magma_tally2_int_t *select, magma_tally2_int_t n,
    magma_tally2DoubleComplex *T,  magma_tally2_int_t ldt,
    magma_tally2DoubleComplex *VL, magma_tally2_int_t ldvl,
    magma_tally2DoubleComplex *VR, magma_tally2_int_t ldvr,
    magma_tally2_int_t mm, magma_tally2_int_t *mout,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ztrevc3_mt(
    magma_tally2_side_t side, magma_tally2_vec_t howmany,
    magma_tally2_int_t *select, magma_tally2_int_t n,
    magma_tally2DoubleComplex *T,  magma_tally2_int_t ldt,
    magma_tally2DoubleComplex *VL, magma_tally2_int_t ldvl,
    magma_tally2DoubleComplex *VR, magma_tally2_int_t ldvr,
    magma_tally2_int_t mm, magma_tally2_int_t *mout,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ztrtri(
    magma_tally2_uplo_t uplo, magma_tally2_diag_t diag, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zunghr(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zungqr(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zungqr2(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zunmbr(
    magma_tally2_vect_t vect, magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *C, magma_tally2_int_t ldc,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zunmlq(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *C, magma_tally2_int_t ldc,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

// not yet implemented
//magma_tally2_int_t magma_tally2_zunmrq( magma_tally2_side_t side, magma_tally2_trans_t trans,
//                          magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
//                          magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
//                          magma_tally2DoubleComplex *tau,
//                          magma_tally2DoubleComplex *C, magma_tally2_int_t ldc,
//                          magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
//                          magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zunmql(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *C, magma_tally2_int_t ldc,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zunmqr(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *C, magma_tally2_int_t ldc,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zunmtr(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A,    magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *C,    magma_tally2_int_t ldc,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 function definitions / Data on CPU / Multi-GPU (alphabetical order)
*/
magma_tally2_int_t
magma_tally2_zgeev_m(
    magma_tally2_vec_t jobvl, magma_tally2_vec_t jobvr, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    #ifdef COMPLEX
    magma_tally2DoubleComplex *w,
    #else
    double *wr, double *wi,
    #endif
    magma_tally2DoubleComplex *VL, magma_tally2_int_t ldvl,
    magma_tally2DoubleComplex *VR, magma_tally2_int_t ldvr,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgehrd_m(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2DoubleComplex *T,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgetrf_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zheevd_m(
    magma_tally2_int_t ngpu,
    magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    double *w,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zheevdx_2stage_m(
    magma_tally2_int_t ngpu,
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zheevdx_m(
    magma_tally2_int_t ngpu,
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhegst_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhegvd_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
    double *w,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhegvdx_2stage_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhegvdx_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

#ifdef REAL
magma_tally2_int_t
magma_tally2_dlaex0_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t n, double *d, double *e,
    double *Q, magma_tally2_int_t ldq,
    double *work, magma_tally2_int_t *iwork,
    magma_tally2_range_t range, double vl, double vu,
    magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dlaex1_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t n, double *d,
    double *Q, magma_tally2_int_t ldq,
    magma_tally2_int_t *indxq, double rho, magma_tally2_int_t cutpnt,
    double *work, magma_tally2_int_t *iwork,
    magma_tally2Double_ptr dwork[],
    magma_tally2_queue_t queues[Magma_tally2MaxGPUs][2],
    magma_tally2_range_t range, double vl, double vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dlaex3_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t k, magma_tally2_int_t n, magma_tally2_int_t n1, double *d,
    double *Q, magma_tally2_int_t ldq, double rho,
    double *dlamda, double *Q2, magma_tally2_int_t *indx,
    magma_tally2_int_t *ctot, double *w, double *s, magma_tally2_int_t *indxq,
    magma_tally2Double_ptr dwork[],
    magma_tally2_queue_t queues[Magma_tally2MaxGPUs][2],
    magma_tally2_range_t range, double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *info);
#endif

magma_tally2_int_t
magma_tally2_zlahr2_m(
    magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2_int_t nb,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *T, magma_tally2_int_t ldt,
    magma_tally2DoubleComplex *Y, magma_tally2_int_t ldy,
    struct zgehrd_data_tally2 *data );

magma_tally2_int_t
magma_tally2_zlahru_m(
    magma_tally2_int_t n, magma_tally2_int_t ihi, magma_tally2_int_t k, magma_tally2_int_t nb,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    struct zgehrd_data_tally2 *data );

magma_tally2_int_t
magma_tally2_zpotrf_m(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zstedx_m(
    magma_tally2_int_t ngpu,
    magma_tally2_range_t range, magma_tally2_int_t n, double vl, double vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, double *d, double *e,
    magma_tally2DoubleComplex *Z, magma_tally2_int_t ldz,
    double *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ztrsm_m(
    magma_tally2_int_t ngpu,
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t transa, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *B, magma_tally2_int_t ldb);

magma_tally2_int_t
magma_tally2_zunghr_m(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *T, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zungqr_m(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *T, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zunmqr_m(
    magma_tally2_int_t ngpu,
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex *A,    magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *C,    magma_tally2_int_t ldc,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zunmtr_m(
    magma_tally2_int_t ngpu,
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A,    magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *C,    magma_tally2_int_t ldc,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 function definitions / Data on GPU (alphabetical order)
*/
magma_tally2_int_t
magma_tally2_zgegqr_gpu(
    magma_tally2_int_t ikind, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dwork, magma_tally2DoubleComplex *work,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgelqf_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgels_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgels3_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqp3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt, magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex_ptr dwork, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqr2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dtau,
    magma_tally2Double_ptr        dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqr2x_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dtau,
    magma_tally2DoubleComplex_ptr dT, magma_tally2DoubleComplex_ptr ddA,
    magma_tally2Double_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqr2x2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dtau,
    magma_tally2DoubleComplex_ptr dT, magma_tally2DoubleComplex_ptr ddA,
    magma_tally2Double_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqr2x3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dtau,
    magma_tally2DoubleComplex_ptr dT,
    magma_tally2DoubleComplex_ptr ddA,
    magma_tally2Double_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqr2x4_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dtau,
    magma_tally2DoubleComplex_ptr dT, magma_tally2DoubleComplex_ptr ddA,
    magma_tally2Double_ptr dwork,
    magma_tally2_queue_t queue,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqrf_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex_ptr dT,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqrf2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex *tau,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqrf2_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dlA[], magma_tally2_int_t ldda,
    magma_tally2DoubleComplex *tau,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqrf3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex_ptr dT,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqrs_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex_ptr dT,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgeqrs3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex_ptr dT,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgerbt_gpu(
    magma_tally2_bool_t gen, magma_tally2_int_t n, magma_tally2_int_t nrhs, 
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda, 
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb, 
    magma_tally2DoubleComplex *U, magma_tally2DoubleComplex *V,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgessm_gpu(
    magma_tally2_order_t order, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2_int_t ib,
    magma_tally2_int_t *ipiv,
    magma_tally2DoubleComplex_ptr dL1, magma_tally2_int_t lddl1,
    magma_tally2DoubleComplex_ptr dL,  magma_tally2_int_t lddl,
    magma_tally2DoubleComplex_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgesv_gpu(
    magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgesv_nopiv_gpu( 
    magma_tally2_int_t n, magma_tally2_int_t nrhs, 
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb, 
                 magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgetf2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgetrf_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgetrf_incpiv_gpu(
    magma_tally2_order_t order, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t ib,
    magma_tally2DoubleComplex    *hA, magma_tally2_int_t ldha,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex    *hL, magma_tally2_int_t ldhl,
    magma_tally2DoubleComplex_ptr dL, magma_tally2_int_t lddl,
    magma_tally2_int_t *ipiv,
    magma_tally2DoubleComplex_ptr dwork, magma_tally2_int_t lddwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgetrf_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr d_lA[], magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgetrf_nopiv_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgetrf2_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t offset,
    magma_tally2DoubleComplex_ptr d_lAT[], magma_tally2_int_t lddat, magma_tally2_int_t *ipiv,
    magma_tally2DoubleComplex_ptr d_lAP[],
    magma_tally2DoubleComplex *W, magma_tally2_int_t ldw,
    magma_tally2_queue_t queues[][2],
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgetri_gpu(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2DoubleComplex_ptr dwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgetrs_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zgetrs_nopiv_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zheevd_gpu(
    magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    double *w,
    magma_tally2DoubleComplex *wA,  magma_tally2_int_t ldwa,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zheevdx_gpu(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    double vl, double vu,
    magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    magma_tally2DoubleComplex *wA,  magma_tally2_int_t ldwa,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally2_int_t
magma_tally2_zheevr_gpu(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    double vl, double vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, double abstol, magma_tally2_int_t *m,
    double *w,
    magma_tally2DoubleComplex_ptr dZ, magma_tally2_int_t lddz,
    magma_tally2_int_t *isuppz,
    magma_tally2DoubleComplex *wA, magma_tally2_int_t ldwa,
    magma_tally2DoubleComplex *wZ, magma_tally2_int_t ldwz,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    double *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zheevx_gpu(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    double abstol, magma_tally2_int_t *m,
    double *w,
    magma_tally2DoubleComplex_ptr dZ, magma_tally2_int_t lddz,
    magma_tally2DoubleComplex *wA, magma_tally2_int_t ldwa,
    magma_tally2DoubleComplex *wZ, magma_tally2_int_t ldwz,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    double *rwork, magma_tally2_int_t *iwork,
    magma_tally2_int_t *ifail,
    magma_tally2_int_t *info);
#endif  // COMPLEX

magma_tally2_int_t
magma_tally2_zhegst_gpu(
    magma_tally2_int_t itype, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhetrd_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    double *d, double *e, magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *wA,  magma_tally2_int_t ldwa,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhetrd_he2hb_mgpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2DoubleComplex_ptr dAmgpu[], magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dTmgpu[], magma_tally2_int_t lddt,
    magma_tally2_int_t ngpu, magma_tally2_int_t distblk,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhetrd_he2hb_mgpu_spec(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2DoubleComplex_ptr dAmgpu[], magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dTmgpu[], magma_tally2_int_t lddt,
    magma_tally2_int_t ngpu, magma_tally2_int_t distblk,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhetrd_mgpu(
    magma_tally2_int_t ngpu, magma_tally2_int_t nqueue,
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    double *d, double *e, magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhetrd2_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    double *d, double *e, magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *wA,  magma_tally2_int_t ldwa,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2DoubleComplex_ptr dwork, magma_tally2_int_t ldwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zhetrf_nopiv_tally2_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zlabrd_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2DoubleComplex     *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    double *d, double *e, magma_tally2DoubleComplex *tauq, magma_tally2DoubleComplex *taup,
    magma_tally2DoubleComplex     *X, magma_tally2_int_t ldx,
    magma_tally2DoubleComplex_ptr dX, magma_tally2_int_t lddx,
    magma_tally2DoubleComplex     *Y, magma_tally2_int_t ldy,
    magma_tally2DoubleComplex_ptr dY, magma_tally2_int_t lddy );

magma_tally2_int_t
magma_tally2_zlaqps_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t offset,
    magma_tally2_int_t nb, magma_tally2_int_t *kb,
    magma_tally2DoubleComplex_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt, magma_tally2DoubleComplex *tau,
    double *vn1, double *vn2,
    magma_tally2DoubleComplex_ptr dauxv,
    magma_tally2DoubleComplex_ptr dF, magma_tally2_int_t lddf);

magma_tally2_int_t
magma_tally2_zlaqps2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t offset,
    magma_tally2_int_t nb, magma_tally2_int_t *kb,
    magma_tally2DoubleComplex_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt,
    magma_tally2DoubleComplex_ptr dtau,
    magma_tally2Double_ptr dvn1, magma_tally2Double_ptr dvn2,
    magma_tally2DoubleComplex_ptr dauxv,
    magma_tally2DoubleComplex_ptr dF, magma_tally2_int_t lddf);

magma_tally2_int_t
magma_tally2_zlaqps3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t offset,
    magma_tally2_int_t nb, magma_tally2_int_t *kb,
    magma_tally2DoubleComplex_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt,
    magma_tally2DoubleComplex_ptr dtau,
    magma_tally2Double_ptr dvn1, magma_tally2Double_ptr dvn2,
    magma_tally2DoubleComplex_ptr dauxv,
    magma_tally2DoubleComplex_ptr dF, magma_tally2_int_t lddf);

magma_tally2_int_t
magma_tally2_zlarf_gpu(
    magma_tally2_int_t m,  magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dv, magma_tally2DoubleComplex_const_ptr dtau,
    magma_tally2DoubleComplex_ptr dC,  magma_tally2_int_t lddc);

magma_tally2_int_t
magma_tally2_zlarfb_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans, magma_tally2_direct_t direct, magma_tally2_storev_t storev,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex_const_ptr dV, magma_tally2_int_t lddv,
    magma_tally2DoubleComplex_const_ptr dT, magma_tally2_int_t lddt,
    magma_tally2DoubleComplex_ptr dC,       magma_tally2_int_t lddc,
    magma_tally2DoubleComplex_ptr dwork,    magma_tally2_int_t ldwork );

magma_tally2_int_t
magma_tally2_zlarfb_gpu_gemm(
    magma_tally2_side_t side, magma_tally2_trans_t trans, magma_tally2_direct_t direct, magma_tally2_storev_t storev,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex_const_ptr dV, magma_tally2_int_t lddv,
    magma_tally2DoubleComplex_const_ptr dT, magma_tally2_int_t lddt,
    magma_tally2DoubleComplex_ptr dC,       magma_tally2_int_t lddc,
    magma_tally2DoubleComplex_ptr dwork,    magma_tally2_int_t ldwork,
    magma_tally2DoubleComplex_ptr dworkvt,  magma_tally2_int_t ldworkvt);

magma_tally2_int_t
magma_tally2_zlarfb2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex_const_ptr dV, magma_tally2_int_t lddv,
    magma_tally2DoubleComplex_const_ptr dT, magma_tally2_int_t lddt,
    magma_tally2DoubleComplex_ptr dC,       magma_tally2_int_t lddc,
    magma_tally2DoubleComplex_ptr dwork,    magma_tally2_int_t ldwork );

magma_tally2_int_t
magma_tally2_zlatrd_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t nb0,
    magma_tally2DoubleComplex *A,  magma_tally2_int_t lda,
    double *e, magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex    *W,       magma_tally2_int_t ldw,
    magma_tally2DoubleComplex_ptr dA[],    magma_tally2_int_t ldda, magma_tally2_int_t offset,
    magma_tally2DoubleComplex_ptr dW[],    magma_tally2_int_t lddw,
    magma_tally2DoubleComplex    *hwork,   magma_tally2_int_t lhwork,
    magma_tally2DoubleComplex_ptr dwork[], magma_tally2_int_t ldwork,
    magma_tally2_queue_t queues[] );

magma_tally2_int_t
magma_tally2_zlauum_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zposv_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zpotf2_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zpotrf_gpu(
    magma_tally2_uplo_t uplo,  magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zpotrf_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr d_lA[], magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zpotrf_mgpu_right(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr d_lA[], magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zpotrf3_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t off_i, magma_tally2_int_t off_j, magma_tally2_int_t nb,
    magma_tally2DoubleComplex_ptr d_lA[],  magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr d_lP[],  magma_tally2_int_t lddp,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda, magma_tally2_int_t h,
    magma_tally2_queue_t queues[][3], magma_tally2_event_t events[][5],
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zpotri_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zpotrs_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zssssm_gpu(
    magma_tally2_order_t order, magma_tally2_int_t m1, magma_tally2_int_t n1,
    magma_tally2_int_t m2, magma_tally2_int_t n2, magma_tally2_int_t k, magma_tally2_int_t ib,
    magma_tally2DoubleComplex_ptr dA1, magma_tally2_int_t ldda1,
    magma_tally2DoubleComplex_ptr dA2, magma_tally2_int_t ldda2,
    magma_tally2DoubleComplex_ptr dL1, magma_tally2_int_t lddl1,
    magma_tally2DoubleComplex_ptr dL2, magma_tally2_int_t lddl2,
    magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ztrtri_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_diag_t diag, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ztsqrt_gpu(
    magma_tally2_int_t *m, magma_tally2_int_t *n,
    magma_tally2DoubleComplex *A1, magma_tally2DoubleComplex *A2, magma_tally2_int_t *lda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t *lwork,
    magma_tally2DoubleComplex_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ztstrf_gpu(
    magma_tally2_order_t order, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t ib, magma_tally2_int_t nb,
    magma_tally2DoubleComplex    *hU, magma_tally2_int_t ldhu,
    magma_tally2DoubleComplex_ptr dU, magma_tally2_int_t lddu,
    magma_tally2DoubleComplex    *hA, magma_tally2_int_t ldha,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex    *hL, magma_tally2_int_t ldhl,
    magma_tally2DoubleComplex_ptr dL, magma_tally2_int_t lddl,
    magma_tally2_int_t *ipiv,
    magma_tally2DoubleComplex *hwork, magma_tally2_int_t ldhwork,
    magma_tally2DoubleComplex_ptr dwork, magma_tally2_int_t lddwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zungqr_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zunmql2_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex_ptr dC, magma_tally2_int_t lddc,
    magma_tally2DoubleComplex *wA, magma_tally2_int_t ldwa,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zunmqr_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex_ptr dC, magma_tally2_int_t lddc,
    magma_tally2DoubleComplex *hwork, magma_tally2_int_t lwork,
    magma_tally2DoubleComplex_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zunmqr2_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex_ptr dC, magma_tally2_int_t lddc,
    magma_tally2DoubleComplex    *wA, magma_tally2_int_t ldwa,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_zunmtr_gpu(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex_ptr dC, magma_tally2_int_t lddc,
    magma_tally2DoubleComplex    *wA, magma_tally2_int_t ldwa,
    magma_tally2_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 utility function definitions
*/

extern const magma_tally2DoubleComplex MAGMA_tally2_Z_NAN;
extern const magma_tally2DoubleComplex MAGMA_tally2_Z_INF;

magma_tally2_int_t
magma_tally2_znan_inf(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    const magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *cnt_nan,
    magma_tally2_int_t *cnt_inf );

magma_tally2_int_t
magma_tally2_znan_inf_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *cnt_nan,
    magma_tally2_int_t *cnt_inf );

void magma_tally2_zprint(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const magma_tally2DoubleComplex *A, magma_tally2_int_t lda );

void magma_tally2_zprint_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr dA, magma_tally2_int_t ldda );

void zpanel_to_q_tally2(
    magma_tally2_uplo_t uplo, magma_tally2_int_t ib,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *work );

void zq_to_panel_tally2(
    magma_tally2_uplo_t uplo, magma_tally2_int_t ib,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *work );

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif /* MAGMA_tally2_Z_H */
