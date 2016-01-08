/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
*/

#ifndef MAGMA_tally4_Z_H
#define MAGMA_tally4_Z_H

#include "magma_tally4_types.h"
#include "magma_tally4_zgehrd_m.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 Auxiliary functions to get the NB used
*/
#ifdef REAL
magma_tally4_int_t magma_tally4_get_dlaex3_m_nb();       // defined in dlaex3_m.cpp
#endif

magma_tally4_int_t magma_tally4_get_zpotrf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_zgetrf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_zgetri_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_zgeqp3_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_zgeqrf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_zgeqlf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_zgehrd_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_zhetrd_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_zhetrf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_zhetrf_nopiv_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_zgelqf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_zgebrd_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_zhegst_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_zgesvd_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_zhegst_nb_m( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_zbulge_nb( magma_tally4_int_t m, magma_tally4_int_t nbthreads );
magma_tally4_int_t magma_tally4_get_zbulge_nb_mgpu( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_zbulge_get_Vblksiz( magma_tally4_int_t m, magma_tally4_int_t nb, magma_tally4_int_t nbthreads );
magma_tally4_int_t magma_tally4_get_zbulge_gcperf();


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 function definitions / Data on CPU (alphabetical order)
*/

#ifdef REAL
// only applicable to real [sd] precisions
void
magma_tally4_dmove_eig(
    magma_tally4_range_t range, magma_tally4_int_t n, double *w,
    magma_tally4_int_t *il, magma_tally4_int_t *iu, double vl, double vu, magma_tally4_int_t *m);

// defined in dlaex3.cpp
void
magma_tally4_zvrange(
    magma_tally4_int_t k, double *d, magma_tally4_int_t *il, magma_tally4_int_t *iu, double vl, double vu);

void
magma_tally4_zirange(
    magma_tally4_int_t k, magma_tally4_int_t *indxq, magma_tally4_int_t *iil, magma_tally4_int_t *iiu, magma_tally4_int_t il, magma_tally4_int_t iu);
#endif

magma_tally4_int_t
magma_tally4_zgebrd(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    double *d, double *e,
    magma_tally4DoubleComplex *tauq, magma_tally4DoubleComplex *taup,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeev(
    magma_tally4_vec_t jobvl, magma_tally4_vec_t jobvr, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    #ifdef COMPLEX
    magma_tally4DoubleComplex *w,
    #else
    double *wr, double *wi,
    #endif
    magma_tally4DoubleComplex *VL, magma_tally4_int_t ldvl,
    magma_tally4DoubleComplex *VR, magma_tally4_int_t ldvr,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgehrd(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4DoubleComplex_ptr dT,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgehrd2(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgelqf(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A,    magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqlf(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A,    magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqp3(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *jpvt, magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqrf(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqrf_ooc(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqrf4(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A,    magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgesdd(
    magma_tally4_vec_t jobz, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    double *s,
    magma_tally4DoubleComplex *U, magma_tally4_int_t ldu,
    magma_tally4DoubleComplex *VT, magma_tally4_int_t ldvt,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally4_int_t *iwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgesv(
    magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4DoubleComplex *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgesv_rbt(
    magma_tally4_bool_t ref, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda, 
    magma_tally4DoubleComplex *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgesvd(
    magma_tally4_vec_t jobu, magma_tally4_vec_t jobvt, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A,    magma_tally4_int_t lda, double *s,
    magma_tally4DoubleComplex *U,    magma_tally4_int_t ldu,
    magma_tally4DoubleComplex *VT,   magma_tally4_int_t ldvt,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgetf2_nopiv(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgetrf(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgetrf_nopiv(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgetrf_piv(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t NB,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgetrf2(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zheevd(
    magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    double *w,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zheevdx(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zheevdx_2stage(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally4_int_t
magma_tally4_zheevr(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    double vl, double vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, double abstol, magma_tally4_int_t *m,
    double *w,
    magma_tally4DoubleComplex *Z, magma_tally4_int_t ldz,
    magma_tally4_int_t *isuppz,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    double *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

// no real [sd] precisions available
magma_tally4_int_t
magma_tally4_zheevx(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    double vl, double vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, double abstol, magma_tally4_int_t *m,
    double *w,
    magma_tally4DoubleComplex *Z, magma_tally4_int_t ldz,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    double *rwork, magma_tally4_int_t *iwork,
    magma_tally4_int_t *ifail,
    magma_tally4_int_t *info);
#endif  // COMPLEX

magma_tally4_int_t
magma_tally4_zhegst(
    magma_tally4_int_t itype, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhegvd(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *B, magma_tally4_int_t ldb,
    double *w, magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhegvdx(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n, magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *B, magma_tally4_int_t ldb,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhegvdx_2stage(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *B, magma_tally4_int_t ldb,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhegvr(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *B, magma_tally4_int_t ldb,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    double abstol, magma_tally4_int_t *m, double *w,
    magma_tally4DoubleComplex *Z, magma_tally4_int_t ldz,
    magma_tally4_int_t *isuppz, magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    double *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhegvx(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n, magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *B, magma_tally4_int_t ldb,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    double abstol, magma_tally4_int_t *m, double *w,
    magma_tally4DoubleComplex *Z, magma_tally4_int_t ldz,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork, double *rwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t *ifail,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhesv(magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
            magma_tally4DoubleComplex *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
            magma_tally4DoubleComplex *B, magma_tally4_int_t ldb,
            magma_tally4_int_t *info );

magma_tally4_int_t
magma_tally4_zhetrd(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    double *d, double *e, magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhetrf(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhetrf_nopiv(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhetrd_hb2st(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t Vblksiz,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    double *d, double *e,
    magma_tally4DoubleComplex *V, magma_tally4_int_t ldv,
    magma_tally4DoubleComplex *TAU, magma_tally4_int_t compT,
    magma_tally4DoubleComplex *T, magma_tally4_int_t ldt);

magma_tally4_int_t
magma_tally4_zhetrd_he2hb(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4DoubleComplex_ptr dT,
    magma_tally4_int_t *info);

#ifdef REAL
// only applicable to real [sd] precisions
magma_tally4_int_t
magma_tally4_dlaex0(
    magma_tally4_int_t n, double *d, double *e,
    double *Q, magma_tally4_int_t ldq,
    double *work, magma_tally4_int_t *iwork,
    magma_tally4Double_ptr dwork,
    magma_tally4_range_t range, double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dlaex1(
    magma_tally4_int_t n, double *d,
    double *Q, magma_tally4_int_t ldq,
    magma_tally4_int_t *indxq, double rho, magma_tally4_int_t cutpnt,
    double *work, magma_tally4_int_t *iwork,
    magma_tally4Double_ptr dwork,
    magma_tally4_range_t range, double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dlaex3(
    magma_tally4_int_t k, magma_tally4_int_t n, magma_tally4_int_t n1, double *d,
    double *Q, magma_tally4_int_t ldq,
    double rho,
    double *dlamda, double *Q2, magma_tally4_int_t *indx,
    magma_tally4_int_t *ctot, double *w, double *s, magma_tally4_int_t *indxq,
    magma_tally4Double_ptr dwork,
    magma_tally4_range_t range, double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *info);
#endif  // REAL

magma_tally4_int_t
magma_tally4_zlahef_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t *kb,
    magma_tally4DoubleComplex    *hA, magma_tally4_int_t lda,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4DoubleComplex_ptr dW, magma_tally4_int_t lddw,
    magma_tally4_queue_t queues[], magma_tally4_event_t event[],
    magma_tally4_int_t *info);

magma_tally4_int_t
zhetrf_nopiv_cpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t ib,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhetrs_nopiv_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhesv_nopiv_gpu(
    magma_tally4_uplo_t uplo,  magma_tally4_int_t n, magma_tally4_int_t nrhs, 
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb, 
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zlahr2(
    magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4_int_t nb,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dV, magma_tally4_int_t lddv,
    magma_tally4DoubleComplex *A,  magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *T,  magma_tally4_int_t ldt,
    magma_tally4DoubleComplex *Y,  magma_tally4_int_t ldy);

magma_tally4_int_t
magma_tally4_zlahru(
    magma_tally4_int_t n, magma_tally4_int_t ihi, magma_tally4_int_t k, magma_tally4_int_t nb,
    magma_tally4DoubleComplex     *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dY, magma_tally4_int_t lddy,
    magma_tally4DoubleComplex_ptr dV, magma_tally4_int_t lddv,
    magma_tally4DoubleComplex_ptr dT,
    magma_tally4DoubleComplex_ptr dwork);

#ifdef REAL
magma_tally4_int_t
magma_tally4_dlaln2(
    magma_tally4_int_t trans, magma_tally4_int_t na, magma_tally4_int_t nw,
    double smin, double ca, const double *A, magma_tally4_int_t lda,
    double d1, double d2,   const double *B, magma_tally4_int_t ldb,
    double wr, double wi, double *X, magma_tally4_int_t ldx,
    double *scale, double *xnorm,
    magma_tally4_int_t *info);
#endif

magma_tally4_int_t
magma_tally4_zlaqps(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t offset,
    magma_tally4_int_t nb, magma_tally4_int_t *kb,
    magma_tally4DoubleComplex *A,  magma_tally4_int_t lda,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt, magma_tally4DoubleComplex *tau, double *vn1, double *vn2,
    magma_tally4DoubleComplex *auxv,
    magma_tally4DoubleComplex *F,  magma_tally4_int_t ldf,
    magma_tally4DoubleComplex_ptr dF, magma_tally4_int_t lddf );

#ifdef REAL
magma_tally4_int_t
magma_tally4_zlaqtrsd(
    magma_tally4_trans_t trans, magma_tally4_int_t n,
    const double *T, magma_tally4_int_t ldt,
    double *x,       magma_tally4_int_t ldx,
    const double *cnorm,
    magma_tally4_int_t *info);
#endif

magma_tally4_int_t
magma_tally4_zlatrd(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    double *e, magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *W, magma_tally4_int_t ldw,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dW, magma_tally4_int_t lddw);

magma_tally4_int_t
magma_tally4_zlatrd2(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4DoubleComplex *A,  magma_tally4_int_t lda,
    double *e, magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *W,  magma_tally4_int_t ldw,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dW, magma_tally4_int_t lddw,
    magma_tally4DoubleComplex_ptr dwork, magma_tally4_int_t ldwork);

#ifdef COMPLEX
magma_tally4_int_t
magma_tally4_zlatrsd(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_diag_t diag, magma_tally4_bool_t normin,
    magma_tally4_int_t n, const magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex lambda,
    magma_tally4DoubleComplex *x,
    double *scale, double *cnorm,
    magma_tally4_int_t *info);
#endif

magma_tally4_int_t
magma_tally4_zlauum(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zposv(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zpotrf(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zpotri(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zstedx(
    magma_tally4_range_t range, magma_tally4_int_t n, double vl, double vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, double *d, double *e,
    magma_tally4DoubleComplex *Z, magma_tally4_int_t ldz,
    double *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4Double_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ztrevc3(
    magma_tally4_side_t side, magma_tally4_vec_t howmany,
    magma_tally4_int_t *select, magma_tally4_int_t n,
    magma_tally4DoubleComplex *T,  magma_tally4_int_t ldt,
    magma_tally4DoubleComplex *VL, magma_tally4_int_t ldvl,
    magma_tally4DoubleComplex *VR, magma_tally4_int_t ldvr,
    magma_tally4_int_t mm, magma_tally4_int_t *mout,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ztrevc3_mt(
    magma_tally4_side_t side, magma_tally4_vec_t howmany,
    magma_tally4_int_t *select, magma_tally4_int_t n,
    magma_tally4DoubleComplex *T,  magma_tally4_int_t ldt,
    magma_tally4DoubleComplex *VL, magma_tally4_int_t ldvl,
    magma_tally4DoubleComplex *VR, magma_tally4_int_t ldvr,
    magma_tally4_int_t mm, magma_tally4_int_t *mout,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ztrtri(
    magma_tally4_uplo_t uplo, magma_tally4_diag_t diag, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zunghr(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zungqr(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zungqr2(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zunmbr(
    magma_tally4_vect_t vect, magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *C, magma_tally4_int_t ldc,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zunmlq(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *C, magma_tally4_int_t ldc,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

// not yet implemented
//magma_tally4_int_t magma_tally4_zunmrq( magma_tally4_side_t side, magma_tally4_trans_t trans,
//                          magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
//                          magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
//                          magma_tally4DoubleComplex *tau,
//                          magma_tally4DoubleComplex *C, magma_tally4_int_t ldc,
//                          magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
//                          magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zunmql(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *C, magma_tally4_int_t ldc,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zunmqr(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *C, magma_tally4_int_t ldc,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zunmtr(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A,    magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *C,    magma_tally4_int_t ldc,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 function definitions / Data on CPU / Multi-GPU (alphabetical order)
*/
magma_tally4_int_t
magma_tally4_zgeev_m(
    magma_tally4_vec_t jobvl, magma_tally4_vec_t jobvr, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    #ifdef COMPLEX
    magma_tally4DoubleComplex *w,
    #else
    double *wr, double *wi,
    #endif
    magma_tally4DoubleComplex *VL, magma_tally4_int_t ldvl,
    magma_tally4DoubleComplex *VR, magma_tally4_int_t ldvr,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgehrd_m(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4DoubleComplex *T,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgetrf_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zheevd_m(
    magma_tally4_int_t ngpu,
    magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    double *w,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zheevdx_2stage_m(
    magma_tally4_int_t ngpu,
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zheevdx_m(
    magma_tally4_int_t ngpu,
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhegst_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t itype, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhegvd_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *B, magma_tally4_int_t ldb,
    double *w,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhegvdx_2stage_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *B, magma_tally4_int_t ldb,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhegvdx_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *B, magma_tally4_int_t ldb,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

#ifdef REAL
magma_tally4_int_t
magma_tally4_dlaex0_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t n, double *d, double *e,
    double *Q, magma_tally4_int_t ldq,
    double *work, magma_tally4_int_t *iwork,
    magma_tally4_range_t range, double vl, double vu,
    magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dlaex1_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t n, double *d,
    double *Q, magma_tally4_int_t ldq,
    magma_tally4_int_t *indxq, double rho, magma_tally4_int_t cutpnt,
    double *work, magma_tally4_int_t *iwork,
    magma_tally4Double_ptr dwork[],
    magma_tally4_queue_t queues[Magma_tally4MaxGPUs][2],
    magma_tally4_range_t range, double vl, double vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dlaex3_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t k, magma_tally4_int_t n, magma_tally4_int_t n1, double *d,
    double *Q, magma_tally4_int_t ldq, double rho,
    double *dlamda, double *Q2, magma_tally4_int_t *indx,
    magma_tally4_int_t *ctot, double *w, double *s, magma_tally4_int_t *indxq,
    magma_tally4Double_ptr dwork[],
    magma_tally4_queue_t queues[Magma_tally4MaxGPUs][2],
    magma_tally4_range_t range, double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *info);
#endif

magma_tally4_int_t
magma_tally4_zlahr2_m(
    magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4_int_t nb,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *T, magma_tally4_int_t ldt,
    magma_tally4DoubleComplex *Y, magma_tally4_int_t ldy,
    struct zgehrd_data *data );

magma_tally4_int_t
magma_tally4_zlahru_m(
    magma_tally4_int_t n, magma_tally4_int_t ihi, magma_tally4_int_t k, magma_tally4_int_t nb,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    struct zgehrd_data *data );

magma_tally4_int_t
magma_tally4_zpotrf_m(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zstedx_m(
    magma_tally4_int_t ngpu,
    magma_tally4_range_t range, magma_tally4_int_t n, double vl, double vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, double *d, double *e,
    magma_tally4DoubleComplex *Z, magma_tally4_int_t ldz,
    double *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ztrsm_m(
    magma_tally4_int_t ngpu,
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t transa, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *B, magma_tally4_int_t ldb);

magma_tally4_int_t
magma_tally4_zunghr_m(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *T, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zungqr_m(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *T, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zunmqr_m(
    magma_tally4_int_t ngpu,
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex *A,    magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *C,    magma_tally4_int_t ldc,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zunmtr_m(
    magma_tally4_int_t ngpu,
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A,    magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *C,    magma_tally4_int_t ldc,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 function definitions / Data on GPU (alphabetical order)
*/
magma_tally4_int_t
magma_tally4_zgegqr_gpu(
    magma_tally4_int_t ikind, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dwork, magma_tally4DoubleComplex *work,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgelqf_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgels_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgels3_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqp3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt, magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex_ptr dwork, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqr2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dtau,
    magma_tally4Double_ptr        dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqr2x_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dtau,
    magma_tally4DoubleComplex_ptr dT, magma_tally4DoubleComplex_ptr ddA,
    magma_tally4Double_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqr2x2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dtau,
    magma_tally4DoubleComplex_ptr dT, magma_tally4DoubleComplex_ptr ddA,
    magma_tally4Double_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqr2x3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dtau,
    magma_tally4DoubleComplex_ptr dT,
    magma_tally4DoubleComplex_ptr ddA,
    magma_tally4Double_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqr2x4_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dtau,
    magma_tally4DoubleComplex_ptr dT, magma_tally4DoubleComplex_ptr ddA,
    magma_tally4Double_ptr dwork,
    magma_tally4_queue_t queue,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqrf_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex_ptr dT,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqrf2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex *tau,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqrf2_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dlA[], magma_tally4_int_t ldda,
    magma_tally4DoubleComplex *tau,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqrf3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex_ptr dT,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqrs_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex_ptr dT,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgeqrs3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex_ptr dT,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgerbt_gpu(
    magma_tally4_bool_t gen, magma_tally4_int_t n, magma_tally4_int_t nrhs, 
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda, 
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb, 
    magma_tally4DoubleComplex *U, magma_tally4DoubleComplex *V,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgessm_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4_int_t ib,
    magma_tally4_int_t *ipiv,
    magma_tally4DoubleComplex_ptr dL1, magma_tally4_int_t lddl1,
    magma_tally4DoubleComplex_ptr dL,  magma_tally4_int_t lddl,
    magma_tally4DoubleComplex_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgesv_gpu(
    magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgesv_nopiv_gpu( 
    magma_tally4_int_t n, magma_tally4_int_t nrhs, 
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb, 
                 magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgetf2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgetrf_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgetrf_incpiv_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t ib,
    magma_tally4DoubleComplex    *hA, magma_tally4_int_t ldha,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex    *hL, magma_tally4_int_t ldhl,
    magma_tally4DoubleComplex_ptr dL, magma_tally4_int_t lddl,
    magma_tally4_int_t *ipiv,
    magma_tally4DoubleComplex_ptr dwork, magma_tally4_int_t lddwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgetrf_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr d_lA[], magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgetrf_nopiv_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgetrf2_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t offset,
    magma_tally4DoubleComplex_ptr d_lAT[], magma_tally4_int_t lddat, magma_tally4_int_t *ipiv,
    magma_tally4DoubleComplex_ptr d_lAP[],
    magma_tally4DoubleComplex *W, magma_tally4_int_t ldw,
    magma_tally4_queue_t queues[][2],
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgetri_gpu(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4DoubleComplex_ptr dwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgetrs_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zgetrs_nopiv_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zheevd_gpu(
    magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    double *w,
    magma_tally4DoubleComplex *wA,  magma_tally4_int_t ldwa,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zheevdx_gpu(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    double vl, double vu,
    magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    magma_tally4DoubleComplex *wA,  magma_tally4_int_t ldwa,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally4_int_t
magma_tally4_zheevr_gpu(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    double vl, double vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, double abstol, magma_tally4_int_t *m,
    double *w,
    magma_tally4DoubleComplex_ptr dZ, magma_tally4_int_t lddz,
    magma_tally4_int_t *isuppz,
    magma_tally4DoubleComplex *wA, magma_tally4_int_t ldwa,
    magma_tally4DoubleComplex *wZ, magma_tally4_int_t ldwz,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    double *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zheevx_gpu(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    double abstol, magma_tally4_int_t *m,
    double *w,
    magma_tally4DoubleComplex_ptr dZ, magma_tally4_int_t lddz,
    magma_tally4DoubleComplex *wA, magma_tally4_int_t ldwa,
    magma_tally4DoubleComplex *wZ, magma_tally4_int_t ldwz,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    double *rwork, magma_tally4_int_t *iwork,
    magma_tally4_int_t *ifail,
    magma_tally4_int_t *info);
#endif  // COMPLEX

magma_tally4_int_t
magma_tally4_zhegst_gpu(
    magma_tally4_int_t itype, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhetrd_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    double *d, double *e, magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *wA,  magma_tally4_int_t ldwa,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhetrd_he2hb_mgpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4DoubleComplex_ptr dAmgpu[], magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dTmgpu[], magma_tally4_int_t lddt,
    magma_tally4_int_t ngpu, magma_tally4_int_t distblk,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhetrd_he2hb_mgpu_spec(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4DoubleComplex_ptr dAmgpu[], magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dTmgpu[], magma_tally4_int_t lddt,
    magma_tally4_int_t ngpu, magma_tally4_int_t distblk,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhetrd_mgpu(
    magma_tally4_int_t ngpu, magma_tally4_int_t nqueue,
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    double *d, double *e, magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhetrd2_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    double *d, double *e, magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *wA,  magma_tally4_int_t ldwa,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4DoubleComplex_ptr dwork, magma_tally4_int_t ldwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zhetrf_nopiv_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zlabrd_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4DoubleComplex     *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    double *d, double *e, magma_tally4DoubleComplex *tauq, magma_tally4DoubleComplex *taup,
    magma_tally4DoubleComplex     *X, magma_tally4_int_t ldx,
    magma_tally4DoubleComplex_ptr dX, magma_tally4_int_t lddx,
    magma_tally4DoubleComplex     *Y, magma_tally4_int_t ldy,
    magma_tally4DoubleComplex_ptr dY, magma_tally4_int_t lddy );

magma_tally4_int_t
magma_tally4_zlaqps_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t offset,
    magma_tally4_int_t nb, magma_tally4_int_t *kb,
    magma_tally4DoubleComplex_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt, magma_tally4DoubleComplex *tau,
    double *vn1, double *vn2,
    magma_tally4DoubleComplex_ptr dauxv,
    magma_tally4DoubleComplex_ptr dF, magma_tally4_int_t lddf);

magma_tally4_int_t
magma_tally4_zlaqps2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t offset,
    magma_tally4_int_t nb, magma_tally4_int_t *kb,
    magma_tally4DoubleComplex_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt,
    magma_tally4DoubleComplex_ptr dtau,
    magma_tally4Double_ptr dvn1, magma_tally4Double_ptr dvn2,
    magma_tally4DoubleComplex_ptr dauxv,
    magma_tally4DoubleComplex_ptr dF, magma_tally4_int_t lddf);

magma_tally4_int_t
magma_tally4_zlaqps3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t offset,
    magma_tally4_int_t nb, magma_tally4_int_t *kb,
    magma_tally4DoubleComplex_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt,
    magma_tally4DoubleComplex_ptr dtau,
    magma_tally4Double_ptr dvn1, magma_tally4Double_ptr dvn2,
    magma_tally4DoubleComplex_ptr dauxv,
    magma_tally4DoubleComplex_ptr dF, magma_tally4_int_t lddf);

magma_tally4_int_t
magma_tally4_zlarf_gpu(
    magma_tally4_int_t m,  magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dv, magma_tally4DoubleComplex_const_ptr dtau,
    magma_tally4DoubleComplex_ptr dC,  magma_tally4_int_t lddc);

magma_tally4_int_t
magma_tally4_zlarfb_gpu(
    magma_tally4_side_t side, magma_tally4_trans_t trans, magma_tally4_direct_t direct, magma_tally4_storev_t storev,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex_const_ptr dV, magma_tally4_int_t lddv,
    magma_tally4DoubleComplex_const_ptr dT, magma_tally4_int_t lddt,
    magma_tally4DoubleComplex_ptr dC,       magma_tally4_int_t lddc,
    magma_tally4DoubleComplex_ptr dwork,    magma_tally4_int_t ldwork );

magma_tally4_int_t
magma_tally4_zlarfb_gpu_gemm(
    magma_tally4_side_t side, magma_tally4_trans_t trans, magma_tally4_direct_t direct, magma_tally4_storev_t storev,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex_const_ptr dV, magma_tally4_int_t lddv,
    magma_tally4DoubleComplex_const_ptr dT, magma_tally4_int_t lddt,
    magma_tally4DoubleComplex_ptr dC,       magma_tally4_int_t lddc,
    magma_tally4DoubleComplex_ptr dwork,    magma_tally4_int_t ldwork,
    magma_tally4DoubleComplex_ptr dworkvt,  magma_tally4_int_t ldworkvt);

magma_tally4_int_t
magma_tally4_zlarfb2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex_const_ptr dV, magma_tally4_int_t lddv,
    magma_tally4DoubleComplex_const_ptr dT, magma_tally4_int_t lddt,
    magma_tally4DoubleComplex_ptr dC,       magma_tally4_int_t lddc,
    magma_tally4DoubleComplex_ptr dwork,    magma_tally4_int_t ldwork );

magma_tally4_int_t
magma_tally4_zlatrd_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t nb0,
    magma_tally4DoubleComplex *A,  magma_tally4_int_t lda,
    double *e, magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex    *W,       magma_tally4_int_t ldw,
    magma_tally4DoubleComplex_ptr dA[],    magma_tally4_int_t ldda, magma_tally4_int_t offset,
    magma_tally4DoubleComplex_ptr dW[],    magma_tally4_int_t lddw,
    magma_tally4DoubleComplex    *hwork,   magma_tally4_int_t lhwork,
    magma_tally4DoubleComplex_ptr dwork[], magma_tally4_int_t ldwork,
    magma_tally4_queue_t queues[] );

magma_tally4_int_t
magma_tally4_zlauum_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zposv_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zpotf2_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zpotrf_gpu(
    magma_tally4_uplo_t uplo,  magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zpotrf_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr d_lA[], magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zpotrf_mgpu_right(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr d_lA[], magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zpotrf3_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t off_i, magma_tally4_int_t off_j, magma_tally4_int_t nb,
    magma_tally4DoubleComplex_ptr d_lA[],  magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr d_lP[],  magma_tally4_int_t lddp,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda, magma_tally4_int_t h,
    magma_tally4_queue_t queues[][3], magma_tally4_event_t events[][5],
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zpotri_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zpotrs_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zssssm_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m1, magma_tally4_int_t n1,
    magma_tally4_int_t m2, magma_tally4_int_t n2, magma_tally4_int_t k, magma_tally4_int_t ib,
    magma_tally4DoubleComplex_ptr dA1, magma_tally4_int_t ldda1,
    magma_tally4DoubleComplex_ptr dA2, magma_tally4_int_t ldda2,
    magma_tally4DoubleComplex_ptr dL1, magma_tally4_int_t lddl1,
    magma_tally4DoubleComplex_ptr dL2, magma_tally4_int_t lddl2,
    magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ztrtri_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_diag_t diag, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ztsqrt_gpu(
    magma_tally4_int_t *m, magma_tally4_int_t *n,
    magma_tally4DoubleComplex *A1, magma_tally4DoubleComplex *A2, magma_tally4_int_t *lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t *lwork,
    magma_tally4DoubleComplex_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ztstrf_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t ib, magma_tally4_int_t nb,
    magma_tally4DoubleComplex    *hU, magma_tally4_int_t ldhu,
    magma_tally4DoubleComplex_ptr dU, magma_tally4_int_t lddu,
    magma_tally4DoubleComplex    *hA, magma_tally4_int_t ldha,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex    *hL, magma_tally4_int_t ldhl,
    magma_tally4DoubleComplex_ptr dL, magma_tally4_int_t lddl,
    magma_tally4_int_t *ipiv,
    magma_tally4DoubleComplex *hwork, magma_tally4_int_t ldhwork,
    magma_tally4DoubleComplex_ptr dwork, magma_tally4_int_t lddwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zungqr_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zunmql2_gpu(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex_ptr dC, magma_tally4_int_t lddc,
    magma_tally4DoubleComplex *wA, magma_tally4_int_t ldwa,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zunmqr_gpu(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex_ptr dC, magma_tally4_int_t lddc,
    magma_tally4DoubleComplex *hwork, magma_tally4_int_t lwork,
    magma_tally4DoubleComplex_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zunmqr2_gpu(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex_ptr dC, magma_tally4_int_t lddc,
    magma_tally4DoubleComplex    *wA, magma_tally4_int_t ldwa,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_zunmtr_gpu(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex_ptr dC, magma_tally4_int_t lddc,
    magma_tally4DoubleComplex    *wA, magma_tally4_int_t ldwa,
    magma_tally4_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 utility function definitions
*/

extern const magma_tally4DoubleComplex MAGMA_tally4_Z_NAN;
extern const magma_tally4DoubleComplex MAGMA_tally4_Z_INF;

magma_tally4_int_t
magma_tally4_znan_inf(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    const magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *cnt_nan,
    magma_tally4_int_t *cnt_inf );

magma_tally4_int_t
magma_tally4_znan_inf_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *cnt_nan,
    magma_tally4_int_t *cnt_inf );

void magma_tally4_zprint(
    magma_tally4_int_t m, magma_tally4_int_t n,
    const magma_tally4DoubleComplex *A, magma_tally4_int_t lda );

void magma_tally4_zprint_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr dA, magma_tally4_int_t ldda );

void zpanel_to_q(
    magma_tally4_uplo_t uplo, magma_tally4_int_t ib,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *work );

void zq_to_panel(
    magma_tally4_uplo_t uplo, magma_tally4_int_t ib,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *work );

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif /* MAGMA_tally4_Z_H */
