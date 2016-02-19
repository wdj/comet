/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally3_z.h normal z -> c, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally3_C_H
#define MAGMA_tally3_C_H

#include "magma_tally3_types.h"
#include "magma_tally3_cgehrd_m.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 Auxiliary functions to get the NB used
*/
#ifdef REAL
magma_tally3_int_t magma_tally3_get_slaex3_m_nb();       // defined in slaex3_m.cpp
#endif

magma_tally3_int_t magma_tally3_get_cpotrf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_cgetrf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_cgetri_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_cgeqp3_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_cgeqrf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_cgeqlf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_cgehrd_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_chetrd_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_chetrf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_chetrf_nopiv_tally3_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_cgelqf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_cgebrd_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_chegst_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_cgesvd_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_chegst_nb_m( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_cbulge_nb( magma_tally3_int_t m, magma_tally3_int_t nbthreads );
magma_tally3_int_t magma_tally3_get_cbulge_nb_mgpu( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_cbulge_get_Vblksiz( magma_tally3_int_t m, magma_tally3_int_t nb, magma_tally3_int_t nbthreads );
magma_tally3_int_t magma_tally3_get_cbulge_gcperf();


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 function definitions / Data on CPU (alphabetical order)
*/

#ifdef REAL
// only applicable to real [sd] precisions
void
magma_tally3_smove_eig(
    magma_tally3_range_t range, magma_tally3_int_t n, float *w,
    magma_tally3_int_t *il, magma_tally3_int_t *iu, float vl, float vu, magma_tally3_int_t *m);

// defined in slaex3.cpp
void
magma_tally3_cvrange(
    magma_tally3_int_t k, float *d, magma_tally3_int_t *il, magma_tally3_int_t *iu, float vl, float vu);

void
magma_tally3_cirange(
    magma_tally3_int_t k, magma_tally3_int_t *indxq, magma_tally3_int_t *iil, magma_tally3_int_t *iiu, magma_tally3_int_t il, magma_tally3_int_t iu);
#endif

magma_tally3_int_t
magma_tally3_cgebrd(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    float *d, float *e,
    magma_tally3FloatComplex *tauq, magma_tally3FloatComplex *taup,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeev(
    magma_tally3_vec_t jobvl, magma_tally3_vec_t jobvr, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    #ifdef COMPLEX
    magma_tally3FloatComplex *w,
    #else
    float *wr, float *wi,
    #endif
    magma_tally3FloatComplex *VL, magma_tally3_int_t ldvl,
    magma_tally3FloatComplex *VR, magma_tally3_int_t ldvr,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgehrd(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3FloatComplex_ptr dT,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgehrd2(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgelqf(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex *A,    magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqlf(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex *A,    magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqp3(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *jpvt, magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqrf(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqrf_ooc(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqrf4(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex *A,    magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgesdd(
    magma_tally3_vec_t jobz, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    float *s,
    magma_tally3FloatComplex *U, magma_tally3_int_t ldu,
    magma_tally3FloatComplex *VT, magma_tally3_int_t ldvt,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally3_int_t *iwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgesv(
    magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3FloatComplex *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgesv_rbt(
    magma_tally3_bool_t ref, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda, 
    magma_tally3FloatComplex *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgesvd(
    magma_tally3_vec_t jobu, magma_tally3_vec_t jobvt, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex *A,    magma_tally3_int_t lda, float *s,
    magma_tally3FloatComplex *U,    magma_tally3_int_t ldu,
    magma_tally3FloatComplex *VT,   magma_tally3_int_t ldvt,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgetf2_nopiv(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgetrf(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgetrf_nopiv(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgetrf_piv(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t NB,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgetrf2(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cheevd(
    magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    float *w,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cheevdx(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cheevdx_2stage(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally3_int_t
magma_tally3_cheevr(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    float vl, float vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, float abstol, magma_tally3_int_t *m,
    float *w,
    magma_tally3FloatComplex *Z, magma_tally3_int_t ldz,
    magma_tally3_int_t *isuppz,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    float *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

// no real [sd] precisions available
magma_tally3_int_t
magma_tally3_cheevx(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    float vl, float vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, float abstol, magma_tally3_int_t *m,
    float *w,
    magma_tally3FloatComplex *Z, magma_tally3_int_t ldz,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    float *rwork, magma_tally3_int_t *iwork,
    magma_tally3_int_t *ifail,
    magma_tally3_int_t *info);
#endif  // COMPLEX

magma_tally3_int_t
magma_tally3_chegst(
    magma_tally3_int_t itype, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chegvd(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *B, magma_tally3_int_t ldb,
    float *w, magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chegvdx(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n, magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *B, magma_tally3_int_t ldb,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chegvdx_2stage(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *B, magma_tally3_int_t ldb,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chegvr(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *B, magma_tally3_int_t ldb,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    float abstol, magma_tally3_int_t *m, float *w,
    magma_tally3FloatComplex *Z, magma_tally3_int_t ldz,
    magma_tally3_int_t *isuppz, magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    float *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chegvx(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n, magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *B, magma_tally3_int_t ldb,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    float abstol, magma_tally3_int_t *m, float *w,
    magma_tally3FloatComplex *Z, magma_tally3_int_t ldz,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork, float *rwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t *ifail,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chesv(magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
            magma_tally3FloatComplex *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
            magma_tally3FloatComplex *B, magma_tally3_int_t ldb,
            magma_tally3_int_t *info );

magma_tally3_int_t
magma_tally3_chetrd(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    float *d, float *e, magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chetrf(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chetrf_nopiv_tally3(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chetrd_hb2st(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t Vblksiz,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    float *d, float *e,
    magma_tally3FloatComplex *V, magma_tally3_int_t ldv,
    magma_tally3FloatComplex *TAU, magma_tally3_int_t compT,
    magma_tally3FloatComplex *T, magma_tally3_int_t ldt);

magma_tally3_int_t
magma_tally3_chetrd_he2hb(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3FloatComplex_ptr dT,
    magma_tally3_int_t *info);

#ifdef REAL
// only applicable to real [sd] precisions
magma_tally3_int_t
magma_tally3_slaex0(
    magma_tally3_int_t n, float *d, float *e,
    float *Q, magma_tally3_int_t ldq,
    float *work, magma_tally3_int_t *iwork,
    magma_tally3Float_ptr dwork,
    magma_tally3_range_t range, float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_slaex1(
    magma_tally3_int_t n, float *d,
    float *Q, magma_tally3_int_t ldq,
    magma_tally3_int_t *indxq, float rho, magma_tally3_int_t cutpnt,
    float *work, magma_tally3_int_t *iwork,
    magma_tally3Float_ptr dwork,
    magma_tally3_range_t range, float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_slaex3(
    magma_tally3_int_t k, magma_tally3_int_t n, magma_tally3_int_t n1, float *d,
    float *Q, magma_tally3_int_t ldq,
    float rho,
    float *dlamda, float *Q2, magma_tally3_int_t *indx,
    magma_tally3_int_t *ctot, float *w, float *s, magma_tally3_int_t *indxq,
    magma_tally3Float_ptr dwork,
    magma_tally3_range_t range, float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *info);
#endif  // REAL

magma_tally3_int_t
magma_tally3_clahef_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3FloatComplex    *hA, magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3FloatComplex_ptr dW, magma_tally3_int_t lddw,
    magma_tally3_queue_t queues[], magma_tally3_event_t event[],
    magma_tally3_int_t *info);

magma_tally3_int_t
chetrf_nopiv_tally3_cpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t ib,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chetrs_nopiv_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chesv_nopiv_gpu(
    magma_tally3_uplo_t uplo,  magma_tally3_int_t n, magma_tally3_int_t nrhs, 
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dB, magma_tally3_int_t lddb, 
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_clahr2(
    magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3_int_t nb,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dV, magma_tally3_int_t lddv,
    magma_tally3FloatComplex *A,  magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *T,  magma_tally3_int_t ldt,
    magma_tally3FloatComplex *Y,  magma_tally3_int_t ldy);

magma_tally3_int_t
magma_tally3_clahru(
    magma_tally3_int_t n, magma_tally3_int_t ihi, magma_tally3_int_t k, magma_tally3_int_t nb,
    magma_tally3FloatComplex     *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dY, magma_tally3_int_t lddy,
    magma_tally3FloatComplex_ptr dV, magma_tally3_int_t lddv,
    magma_tally3FloatComplex_ptr dT,
    magma_tally3FloatComplex_ptr dwork);

#ifdef REAL
magma_tally3_int_t
magma_tally3_slaln2(
    magma_tally3_int_t trans, magma_tally3_int_t na, magma_tally3_int_t nw,
    float smin, float ca, const float *A, magma_tally3_int_t lda,
    float d1, float d2,   const float *B, magma_tally3_int_t ldb,
    float wr, float wi, float *X, magma_tally3_int_t ldx,
    float *scale, float *xnorm,
    magma_tally3_int_t *info);
#endif

magma_tally3_int_t
magma_tally3_claqps(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3FloatComplex *A,  magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt, magma_tally3FloatComplex *tau, float *vn1, float *vn2,
    magma_tally3FloatComplex *auxv,
    magma_tally3FloatComplex *F,  magma_tally3_int_t ldf,
    magma_tally3FloatComplex_ptr dF, magma_tally3_int_t lddf );

#ifdef REAL
magma_tally3_int_t
magma_tally3_claqtrsd(
    magma_tally3_trans_t trans, magma_tally3_int_t n,
    const float *T, magma_tally3_int_t ldt,
    float *x,       magma_tally3_int_t ldx,
    const float *cnorm,
    magma_tally3_int_t *info);
#endif

magma_tally3_int_t
magma_tally3_clatrd(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    float *e, magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *W, magma_tally3_int_t ldw,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dW, magma_tally3_int_t lddw);

magma_tally3_int_t
magma_tally3_clatrd2(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3FloatComplex *A,  magma_tally3_int_t lda,
    float *e, magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *W,  magma_tally3_int_t ldw,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dW, magma_tally3_int_t lddw,
    magma_tally3FloatComplex_ptr dwork, magma_tally3_int_t ldwork);

#ifdef COMPLEX
magma_tally3_int_t
magma_tally3_clatrsd(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_diag_t diag, magma_tally3_bool_t normin,
    magma_tally3_int_t n, const magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex lambda,
    magma_tally3FloatComplex *x,
    float *scale, float *cnorm,
    magma_tally3_int_t *info);
#endif

magma_tally3_int_t
magma_tally3_clauum(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cposv(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cpotrf(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cpotri(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cstedx(
    magma_tally3_range_t range, magma_tally3_int_t n, float vl, float vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, float *d, float *e,
    magma_tally3FloatComplex *Z, magma_tally3_int_t ldz,
    float *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3Float_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ctrevc3(
    magma_tally3_side_t side, magma_tally3_vec_t howmany,
    magma_tally3_int_t *select, magma_tally3_int_t n,
    magma_tally3FloatComplex *T,  magma_tally3_int_t ldt,
    magma_tally3FloatComplex *VL, magma_tally3_int_t ldvl,
    magma_tally3FloatComplex *VR, magma_tally3_int_t ldvr,
    magma_tally3_int_t mm, magma_tally3_int_t *mout,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ctrevc3_mt(
    magma_tally3_side_t side, magma_tally3_vec_t howmany,
    magma_tally3_int_t *select, magma_tally3_int_t n,
    magma_tally3FloatComplex *T,  magma_tally3_int_t ldt,
    magma_tally3FloatComplex *VL, magma_tally3_int_t ldvl,
    magma_tally3FloatComplex *VR, magma_tally3_int_t ldvr,
    magma_tally3_int_t mm, magma_tally3_int_t *mout,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ctrtri(
    magma_tally3_uplo_t uplo, magma_tally3_diag_t diag, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cunghr(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cungqr(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cungqr2(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cunmbr(
    magma_tally3_vect_t vect, magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *C, magma_tally3_int_t ldc,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cunmlq(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *C, magma_tally3_int_t ldc,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

// not yet implemented
//magma_tally3_int_t magma_tally3_cunmrq( magma_tally3_side_t side, magma_tally3_trans_t trans,
//                          magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
//                          magma_tally3FloatComplex *A, magma_tally3_int_t lda,
//                          magma_tally3FloatComplex *tau,
//                          magma_tally3FloatComplex *C, magma_tally3_int_t ldc,
//                          magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
//                          magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cunmql(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *C, magma_tally3_int_t ldc,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cunmqr(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *C, magma_tally3_int_t ldc,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cunmtr(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex *A,    magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *C,    magma_tally3_int_t ldc,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 function definitions / Data on CPU / Multi-GPU (alphabetical order)
*/
magma_tally3_int_t
magma_tally3_cgeev_m(
    magma_tally3_vec_t jobvl, magma_tally3_vec_t jobvr, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    #ifdef COMPLEX
    magma_tally3FloatComplex *w,
    #else
    float *wr, float *wi,
    #endif
    magma_tally3FloatComplex *VL, magma_tally3_int_t ldvl,
    magma_tally3FloatComplex *VR, magma_tally3_int_t ldvr,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgehrd_m(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3FloatComplex *T,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgetrf_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cheevd_m(
    magma_tally3_int_t ngpu,
    magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    float *w,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cheevdx_2stage_m(
    magma_tally3_int_t ngpu,
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cheevdx_m(
    magma_tally3_int_t ngpu,
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chegst_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chegvd_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *B, magma_tally3_int_t ldb,
    float *w,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chegvdx_2stage_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *B, magma_tally3_int_t ldb,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chegvdx_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *B, magma_tally3_int_t ldb,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

#ifdef REAL
magma_tally3_int_t
magma_tally3_slaex0_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t n, float *d, float *e,
    float *Q, magma_tally3_int_t ldq,
    float *work, magma_tally3_int_t *iwork,
    magma_tally3_range_t range, float vl, float vu,
    magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_slaex1_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t n, float *d,
    float *Q, magma_tally3_int_t ldq,
    magma_tally3_int_t *indxq, float rho, magma_tally3_int_t cutpnt,
    float *work, magma_tally3_int_t *iwork,
    magma_tally3Float_ptr dwork[],
    magma_tally3_queue_t queues[Magma_tally3MaxGPUs][2],
    magma_tally3_range_t range, float vl, float vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_slaex3_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t k, magma_tally3_int_t n, magma_tally3_int_t n1, float *d,
    float *Q, magma_tally3_int_t ldq, float rho,
    float *dlamda, float *Q2, magma_tally3_int_t *indx,
    magma_tally3_int_t *ctot, float *w, float *s, magma_tally3_int_t *indxq,
    magma_tally3Float_ptr dwork[],
    magma_tally3_queue_t queues[Magma_tally3MaxGPUs][2],
    magma_tally3_range_t range, float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *info);
#endif

magma_tally3_int_t
magma_tally3_clahr2_m(
    magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3_int_t nb,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *T, magma_tally3_int_t ldt,
    magma_tally3FloatComplex *Y, magma_tally3_int_t ldy,
    struct cgehrd_data_tally3 *data );

magma_tally3_int_t
magma_tally3_clahru_m(
    magma_tally3_int_t n, magma_tally3_int_t ihi, magma_tally3_int_t k, magma_tally3_int_t nb,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    struct cgehrd_data_tally3 *data );

magma_tally3_int_t
magma_tally3_cpotrf_m(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cstedx_m(
    magma_tally3_int_t ngpu,
    magma_tally3_range_t range, magma_tally3_int_t n, float vl, float vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, float *d, float *e,
    magma_tally3FloatComplex *Z, magma_tally3_int_t ldz,
    float *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ctrsm_m(
    magma_tally3_int_t ngpu,
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transa, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *B, magma_tally3_int_t ldb);

magma_tally3_int_t
magma_tally3_cunghr_m(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *T, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cungqr_m(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *T, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cunmqr_m(
    magma_tally3_int_t ngpu,
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex *A,    magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *C,    magma_tally3_int_t ldc,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cunmtr_m(
    magma_tally3_int_t ngpu,
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex *A,    magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *C,    magma_tally3_int_t ldc,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 function definitions / Data on GPU (alphabetical order)
*/
magma_tally3_int_t
magma_tally3_cgegqr_gpu(
    magma_tally3_int_t ikind, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dwork, magma_tally3FloatComplex *work,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgelqf_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgels_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3FloatComplex *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgels3_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3FloatComplex *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqp3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt, magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex_ptr dwork, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqr2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dtau,
    magma_tally3Float_ptr        dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqr2x_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dtau,
    magma_tally3FloatComplex_ptr dT, magma_tally3FloatComplex_ptr ddA,
    magma_tally3Float_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqr2x2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dtau,
    magma_tally3FloatComplex_ptr dT, magma_tally3FloatComplex_ptr ddA,
    magma_tally3Float_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqr2x3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dtau,
    magma_tally3FloatComplex_ptr dT,
    magma_tally3FloatComplex_ptr ddA,
    magma_tally3Float_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqr2x4_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dtau,
    magma_tally3FloatComplex_ptr dT, magma_tally3FloatComplex_ptr ddA,
    magma_tally3Float_ptr dwork,
    magma_tally3_queue_t queue,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqrf_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex_ptr dT,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqrf2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex *tau,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqrf2_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dlA[], magma_tally3_int_t ldda,
    magma_tally3FloatComplex *tau,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqrf3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex_ptr dT,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqrs_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex_ptr dT,
    magma_tally3FloatComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3FloatComplex *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgeqrs3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex_ptr dT,
    magma_tally3FloatComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3FloatComplex *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgerbt_gpu(
    magma_tally3_bool_t gen, magma_tally3_int_t n, magma_tally3_int_t nrhs, 
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda, 
    magma_tally3FloatComplex_ptr dB, magma_tally3_int_t lddb, 
    magma_tally3FloatComplex *U, magma_tally3FloatComplex *V,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgessm_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3_int_t ib,
    magma_tally3_int_t *ipiv,
    magma_tally3FloatComplex_ptr dL1, magma_tally3_int_t lddl1,
    magma_tally3FloatComplex_ptr dL,  magma_tally3_int_t lddl,
    magma_tally3FloatComplex_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgesv_gpu(
    magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3FloatComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgesv_nopiv_gpu( 
    magma_tally3_int_t n, magma_tally3_int_t nrhs, 
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dB, magma_tally3_int_t lddb, 
                 magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgetf2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgetrf_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgetrf_incpiv_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t ib,
    magma_tally3FloatComplex    *hA, magma_tally3_int_t ldha,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex    *hL, magma_tally3_int_t ldhl,
    magma_tally3FloatComplex_ptr dL, magma_tally3_int_t lddl,
    magma_tally3_int_t *ipiv,
    magma_tally3FloatComplex_ptr dwork, magma_tally3_int_t lddwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgetrf_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr d_lA[], magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgetrf_nopiv_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgetrf2_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t offset,
    magma_tally3FloatComplex_ptr d_lAT[], magma_tally3_int_t lddat, magma_tally3_int_t *ipiv,
    magma_tally3FloatComplex_ptr d_lAP[],
    magma_tally3FloatComplex *W, magma_tally3_int_t ldw,
    magma_tally3_queue_t queues[][2],
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgetri_gpu(
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3FloatComplex_ptr dwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgetrs_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3FloatComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cgetrs_nopiv_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cheevd_gpu(
    magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    float *w,
    magma_tally3FloatComplex *wA,  magma_tally3_int_t ldwa,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cheevdx_gpu(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    float vl, float vu,
    magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    magma_tally3FloatComplex *wA,  magma_tally3_int_t ldwa,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally3_int_t
magma_tally3_cheevr_gpu(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    float vl, float vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, float abstol, magma_tally3_int_t *m,
    float *w,
    magma_tally3FloatComplex_ptr dZ, magma_tally3_int_t lddz,
    magma_tally3_int_t *isuppz,
    magma_tally3FloatComplex *wA, magma_tally3_int_t ldwa,
    magma_tally3FloatComplex *wZ, magma_tally3_int_t ldwz,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    float *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cheevx_gpu(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    float abstol, magma_tally3_int_t *m,
    float *w,
    magma_tally3FloatComplex_ptr dZ, magma_tally3_int_t lddz,
    magma_tally3FloatComplex *wA, magma_tally3_int_t ldwa,
    magma_tally3FloatComplex *wZ, magma_tally3_int_t ldwz,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    float *rwork, magma_tally3_int_t *iwork,
    magma_tally3_int_t *ifail,
    magma_tally3_int_t *info);
#endif  // COMPLEX

magma_tally3_int_t
magma_tally3_chegst_gpu(
    magma_tally3_int_t itype, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chetrd_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    float *d, float *e, magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *wA,  magma_tally3_int_t ldwa,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chetrd_he2hb_mgpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3FloatComplex_ptr dAmgpu[], magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dTmgpu[], magma_tally3_int_t lddt,
    magma_tally3_int_t ngpu, magma_tally3_int_t distblk,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chetrd_he2hb_mgpu_spec(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3FloatComplex_ptr dAmgpu[], magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dTmgpu[], magma_tally3_int_t lddt,
    magma_tally3_int_t ngpu, magma_tally3_int_t distblk,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chetrd_mgpu(
    magma_tally3_int_t ngpu, magma_tally3_int_t nqueue,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    float *d, float *e, magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chetrd2_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    float *d, float *e, magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *wA,  magma_tally3_int_t ldwa,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3FloatComplex_ptr dwork, magma_tally3_int_t ldwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_chetrf_nopiv_tally3_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_clabrd_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3FloatComplex     *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    float *d, float *e, magma_tally3FloatComplex *tauq, magma_tally3FloatComplex *taup,
    magma_tally3FloatComplex     *X, magma_tally3_int_t ldx,
    magma_tally3FloatComplex_ptr dX, magma_tally3_int_t lddx,
    magma_tally3FloatComplex     *Y, magma_tally3_int_t ldy,
    magma_tally3FloatComplex_ptr dY, magma_tally3_int_t lddy );

magma_tally3_int_t
magma_tally3_claqps_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3FloatComplex_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt, magma_tally3FloatComplex *tau,
    float *vn1, float *vn2,
    magma_tally3FloatComplex_ptr dauxv,
    magma_tally3FloatComplex_ptr dF, magma_tally3_int_t lddf);

magma_tally3_int_t
magma_tally3_claqps2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3FloatComplex_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt,
    magma_tally3FloatComplex_ptr dtau,
    magma_tally3Float_ptr dvn1, magma_tally3Float_ptr dvn2,
    magma_tally3FloatComplex_ptr dauxv,
    magma_tally3FloatComplex_ptr dF, magma_tally3_int_t lddf);

magma_tally3_int_t
magma_tally3_claqps3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3FloatComplex_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt,
    magma_tally3FloatComplex_ptr dtau,
    magma_tally3Float_ptr dvn1, magma_tally3Float_ptr dvn2,
    magma_tally3FloatComplex_ptr dauxv,
    magma_tally3FloatComplex_ptr dF, magma_tally3_int_t lddf);

magma_tally3_int_t
magma_tally3_clarf_gpu(
    magma_tally3_int_t m,  magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dv, magma_tally3FloatComplex_const_ptr dtau,
    magma_tally3FloatComplex_ptr dC,  magma_tally3_int_t lddc);

magma_tally3_int_t
magma_tally3_clarfb_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans, magma_tally3_direct_t direct, magma_tally3_storev_t storev,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex_const_ptr dV, magma_tally3_int_t lddv,
    magma_tally3FloatComplex_const_ptr dT, magma_tally3_int_t lddt,
    magma_tally3FloatComplex_ptr dC,       magma_tally3_int_t lddc,
    magma_tally3FloatComplex_ptr dwork,    magma_tally3_int_t ldwork );

magma_tally3_int_t
magma_tally3_clarfb_gpu_gemm(
    magma_tally3_side_t side, magma_tally3_trans_t trans, magma_tally3_direct_t direct, magma_tally3_storev_t storev,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex_const_ptr dV, magma_tally3_int_t lddv,
    magma_tally3FloatComplex_const_ptr dT, magma_tally3_int_t lddt,
    magma_tally3FloatComplex_ptr dC,       magma_tally3_int_t lddc,
    magma_tally3FloatComplex_ptr dwork,    magma_tally3_int_t ldwork,
    magma_tally3FloatComplex_ptr dworkvt,  magma_tally3_int_t ldworkvt);

magma_tally3_int_t
magma_tally3_clarfb2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex_const_ptr dV, magma_tally3_int_t lddv,
    magma_tally3FloatComplex_const_ptr dT, magma_tally3_int_t lddt,
    magma_tally3FloatComplex_ptr dC,       magma_tally3_int_t lddc,
    magma_tally3FloatComplex_ptr dwork,    magma_tally3_int_t ldwork );

magma_tally3_int_t
magma_tally3_clatrd_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t nb0,
    magma_tally3FloatComplex *A,  magma_tally3_int_t lda,
    float *e, magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex    *W,       magma_tally3_int_t ldw,
    magma_tally3FloatComplex_ptr dA[],    magma_tally3_int_t ldda, magma_tally3_int_t offset,
    magma_tally3FloatComplex_ptr dW[],    magma_tally3_int_t lddw,
    magma_tally3FloatComplex    *hwork,   magma_tally3_int_t lhwork,
    magma_tally3FloatComplex_ptr dwork[], magma_tally3_int_t ldwork,
    magma_tally3_queue_t queues[] );

magma_tally3_int_t
magma_tally3_clauum_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cposv_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cpotf2_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cpotrf_gpu(
    magma_tally3_uplo_t uplo,  magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cpotrf_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr d_lA[], magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cpotrf_mgpu_right(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr d_lA[], magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cpotrf3_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t off_i, magma_tally3_int_t off_j, magma_tally3_int_t nb,
    magma_tally3FloatComplex_ptr d_lA[],  magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr d_lP[],  magma_tally3_int_t lddp,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda, magma_tally3_int_t h,
    magma_tally3_queue_t queues[][3], magma_tally3_event_t events[][5],
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cpotri_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cpotrs_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cssssm_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m1, magma_tally3_int_t n1,
    magma_tally3_int_t m2, magma_tally3_int_t n2, magma_tally3_int_t k, magma_tally3_int_t ib,
    magma_tally3FloatComplex_ptr dA1, magma_tally3_int_t ldda1,
    magma_tally3FloatComplex_ptr dA2, magma_tally3_int_t ldda2,
    magma_tally3FloatComplex_ptr dL1, magma_tally3_int_t lddl1,
    magma_tally3FloatComplex_ptr dL2, magma_tally3_int_t lddl2,
    magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ctrtri_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_diag_t diag, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ctsqrt_gpu(
    magma_tally3_int_t *m, magma_tally3_int_t *n,
    magma_tally3FloatComplex *A1, magma_tally3FloatComplex *A2, magma_tally3_int_t *lda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t *lwork,
    magma_tally3FloatComplex_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ctstrf_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t ib, magma_tally3_int_t nb,
    magma_tally3FloatComplex    *hU, magma_tally3_int_t ldhu,
    magma_tally3FloatComplex_ptr dU, magma_tally3_int_t lddu,
    magma_tally3FloatComplex    *hA, magma_tally3_int_t ldha,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex    *hL, magma_tally3_int_t ldhl,
    magma_tally3FloatComplex_ptr dL, magma_tally3_int_t lddl,
    magma_tally3_int_t *ipiv,
    magma_tally3FloatComplex *hwork, magma_tally3_int_t ldhwork,
    magma_tally3FloatComplex_ptr dwork, magma_tally3_int_t lddwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cungqr_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cunmql2_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex_ptr dC, magma_tally3_int_t lddc,
    magma_tally3FloatComplex *wA, magma_tally3_int_t ldwa,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cunmqr_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex_ptr dC, magma_tally3_int_t lddc,
    magma_tally3FloatComplex *hwork, magma_tally3_int_t lwork,
    magma_tally3FloatComplex_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cunmqr2_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex_ptr dC, magma_tally3_int_t lddc,
    magma_tally3FloatComplex    *wA, magma_tally3_int_t ldwa,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_cunmtr_gpu(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex_ptr dC, magma_tally3_int_t lddc,
    magma_tally3FloatComplex    *wA, magma_tally3_int_t ldwa,
    magma_tally3_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 utility function definitions
*/

extern const magma_tally3FloatComplex MAGMA_tally3_C_NAN;
extern const magma_tally3FloatComplex MAGMA_tally3_C_INF;

magma_tally3_int_t
magma_tally3_cnan_inf(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3_int_t *cnt_nan,
    magma_tally3_int_t *cnt_inf );

magma_tally3_int_t
magma_tally3_cnan_inf_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *cnt_nan,
    magma_tally3_int_t *cnt_inf );

void magma_tally3_cprint(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const magma_tally3FloatComplex *A, magma_tally3_int_t lda );

void magma_tally3_cprint_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr dA, magma_tally3_int_t ldda );

void cpanel_to_q_tally3(
    magma_tally3_uplo_t uplo, magma_tally3_int_t ib,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *work );

void cq_to_panel_tally3(
    magma_tally3_uplo_t uplo, magma_tally3_int_t ib,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex *work );

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif /* MAGMA_tally3_C_H */
