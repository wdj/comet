/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally4_z.h normal z -> c, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally4_C_H
#define MAGMA_tally4_C_H

#include "magma_tally4_types.h"
#include "magma_tally4_cgehrd_m.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 Auxiliary functions to get the NB used
*/
#ifdef REAL
magma_tally4_int_t magma_tally4_get_slaex3_m_nb();       // defined in slaex3_m.cpp
#endif

magma_tally4_int_t magma_tally4_get_cpotrf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_cgetrf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_cgetri_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_cgeqp3_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_cgeqrf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_cgeqlf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_cgehrd_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_chetrd_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_chetrf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_chetrf_nopiv_tally4_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_cgelqf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_cgebrd_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_chegst_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_cgesvd_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_chegst_nb_m( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_cbulge_nb( magma_tally4_int_t m, magma_tally4_int_t nbthreads );
magma_tally4_int_t magma_tally4_get_cbulge_nb_mgpu( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_cbulge_get_Vblksiz( magma_tally4_int_t m, magma_tally4_int_t nb, magma_tally4_int_t nbthreads );
magma_tally4_int_t magma_tally4_get_cbulge_gcperf();


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 function definitions / Data on CPU (alphabetical order)
*/

#ifdef REAL
// only applicable to real [sd] precisions
void
magma_tally4_smove_eig(
    magma_tally4_range_t range, magma_tally4_int_t n, float *w,
    magma_tally4_int_t *il, magma_tally4_int_t *iu, float vl, float vu, magma_tally4_int_t *m);

// defined in slaex3.cpp
void
magma_tally4_cvrange(
    magma_tally4_int_t k, float *d, magma_tally4_int_t *il, magma_tally4_int_t *iu, float vl, float vu);

void
magma_tally4_cirange(
    magma_tally4_int_t k, magma_tally4_int_t *indxq, magma_tally4_int_t *iil, magma_tally4_int_t *iiu, magma_tally4_int_t il, magma_tally4_int_t iu);
#endif

magma_tally4_int_t
magma_tally4_cgebrd(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    float *d, float *e,
    magma_tally4FloatComplex *tauq, magma_tally4FloatComplex *taup,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeev(
    magma_tally4_vec_t jobvl, magma_tally4_vec_t jobvr, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    #ifdef COMPLEX
    magma_tally4FloatComplex *w,
    #else
    float *wr, float *wi,
    #endif
    magma_tally4FloatComplex *VL, magma_tally4_int_t ldvl,
    magma_tally4FloatComplex *VR, magma_tally4_int_t ldvr,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgehrd(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4FloatComplex_ptr dT,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgehrd2(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgelqf(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A,    magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqlf(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A,    magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqp3(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *jpvt, magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqrf(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqrf_ooc(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqrf4(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A,    magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgesdd(
    magma_tally4_vec_t jobz, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    float *s,
    magma_tally4FloatComplex *U, magma_tally4_int_t ldu,
    magma_tally4FloatComplex *VT, magma_tally4_int_t ldvt,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally4_int_t *iwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgesv(
    magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4FloatComplex *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgesv_rbt(
    magma_tally4_bool_t ref, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda, 
    magma_tally4FloatComplex *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgesvd(
    magma_tally4_vec_t jobu, magma_tally4_vec_t jobvt, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A,    magma_tally4_int_t lda, float *s,
    magma_tally4FloatComplex *U,    magma_tally4_int_t ldu,
    magma_tally4FloatComplex *VT,   magma_tally4_int_t ldvt,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgetf2_nopiv(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgetrf(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgetrf_nopiv(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgetrf_piv(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t NB,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgetrf2(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cheevd(
    magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    float *w,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cheevdx(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cheevdx_2stage(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally4_int_t
magma_tally4_cheevr(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    float vl, float vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, float abstol, magma_tally4_int_t *m,
    float *w,
    magma_tally4FloatComplex *Z, magma_tally4_int_t ldz,
    magma_tally4_int_t *isuppz,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    float *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

// no real [sd] precisions available
magma_tally4_int_t
magma_tally4_cheevx(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    float vl, float vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, float abstol, magma_tally4_int_t *m,
    float *w,
    magma_tally4FloatComplex *Z, magma_tally4_int_t ldz,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    float *rwork, magma_tally4_int_t *iwork,
    magma_tally4_int_t *ifail,
    magma_tally4_int_t *info);
#endif  // COMPLEX

magma_tally4_int_t
magma_tally4_chegst(
    magma_tally4_int_t itype, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chegvd(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *B, magma_tally4_int_t ldb,
    float *w, magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chegvdx(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n, magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *B, magma_tally4_int_t ldb,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chegvdx_2stage(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *B, magma_tally4_int_t ldb,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chegvr(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *B, magma_tally4_int_t ldb,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    float abstol, magma_tally4_int_t *m, float *w,
    magma_tally4FloatComplex *Z, magma_tally4_int_t ldz,
    magma_tally4_int_t *isuppz, magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    float *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chegvx(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n, magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *B, magma_tally4_int_t ldb,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    float abstol, magma_tally4_int_t *m, float *w,
    magma_tally4FloatComplex *Z, magma_tally4_int_t ldz,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork, float *rwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t *ifail,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chesv(magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
            magma_tally4FloatComplex *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
            magma_tally4FloatComplex *B, magma_tally4_int_t ldb,
            magma_tally4_int_t *info );

magma_tally4_int_t
magma_tally4_chetrd(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    float *d, float *e, magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chetrf(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chetrf_nopiv_tally4(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chetrd_hb2st(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t Vblksiz,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    float *d, float *e,
    magma_tally4FloatComplex *V, magma_tally4_int_t ldv,
    magma_tally4FloatComplex *TAU, magma_tally4_int_t compT,
    magma_tally4FloatComplex *T, magma_tally4_int_t ldt);

magma_tally4_int_t
magma_tally4_chetrd_he2hb(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4FloatComplex_ptr dT,
    magma_tally4_int_t *info);

#ifdef REAL
// only applicable to real [sd] precisions
magma_tally4_int_t
magma_tally4_slaex0(
    magma_tally4_int_t n, float *d, float *e,
    float *Q, magma_tally4_int_t ldq,
    float *work, magma_tally4_int_t *iwork,
    magma_tally4Float_ptr dwork,
    magma_tally4_range_t range, float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_slaex1(
    magma_tally4_int_t n, float *d,
    float *Q, magma_tally4_int_t ldq,
    magma_tally4_int_t *indxq, float rho, magma_tally4_int_t cutpnt,
    float *work, magma_tally4_int_t *iwork,
    magma_tally4Float_ptr dwork,
    magma_tally4_range_t range, float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_slaex3(
    magma_tally4_int_t k, magma_tally4_int_t n, magma_tally4_int_t n1, float *d,
    float *Q, magma_tally4_int_t ldq,
    float rho,
    float *dlamda, float *Q2, magma_tally4_int_t *indx,
    magma_tally4_int_t *ctot, float *w, float *s, magma_tally4_int_t *indxq,
    magma_tally4Float_ptr dwork,
    magma_tally4_range_t range, float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *info);
#endif  // REAL

magma_tally4_int_t
magma_tally4_clahef_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t *kb,
    magma_tally4FloatComplex    *hA, magma_tally4_int_t lda,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4FloatComplex_ptr dW, magma_tally4_int_t lddw,
    magma_tally4_queue_t queues[], magma_tally4_event_t event[],
    magma_tally4_int_t *info);

magma_tally4_int_t
chetrf_nopiv_tally4_cpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t ib,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chetrs_nopiv_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chesv_nopiv_gpu(
    magma_tally4_uplo_t uplo,  magma_tally4_int_t n, magma_tally4_int_t nrhs, 
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb, 
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_clahr2(
    magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4_int_t nb,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dV, magma_tally4_int_t lddv,
    magma_tally4FloatComplex *A,  magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *T,  magma_tally4_int_t ldt,
    magma_tally4FloatComplex *Y,  magma_tally4_int_t ldy);

magma_tally4_int_t
magma_tally4_clahru(
    magma_tally4_int_t n, magma_tally4_int_t ihi, magma_tally4_int_t k, magma_tally4_int_t nb,
    magma_tally4FloatComplex     *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dY, magma_tally4_int_t lddy,
    magma_tally4FloatComplex_ptr dV, magma_tally4_int_t lddv,
    magma_tally4FloatComplex_ptr dT,
    magma_tally4FloatComplex_ptr dwork);

#ifdef REAL
magma_tally4_int_t
magma_tally4_slaln2(
    magma_tally4_int_t trans, magma_tally4_int_t na, magma_tally4_int_t nw,
    float smin, float ca, const float *A, magma_tally4_int_t lda,
    float d1, float d2,   const float *B, magma_tally4_int_t ldb,
    float wr, float wi, float *X, magma_tally4_int_t ldx,
    float *scale, float *xnorm,
    magma_tally4_int_t *info);
#endif

magma_tally4_int_t
magma_tally4_claqps(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t offset,
    magma_tally4_int_t nb, magma_tally4_int_t *kb,
    magma_tally4FloatComplex *A,  magma_tally4_int_t lda,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt, magma_tally4FloatComplex *tau, float *vn1, float *vn2,
    magma_tally4FloatComplex *auxv,
    magma_tally4FloatComplex *F,  magma_tally4_int_t ldf,
    magma_tally4FloatComplex_ptr dF, magma_tally4_int_t lddf );

#ifdef REAL
magma_tally4_int_t
magma_tally4_claqtrsd(
    magma_tally4_trans_t trans, magma_tally4_int_t n,
    const float *T, magma_tally4_int_t ldt,
    float *x,       magma_tally4_int_t ldx,
    const float *cnorm,
    magma_tally4_int_t *info);
#endif

magma_tally4_int_t
magma_tally4_clatrd(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    float *e, magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *W, magma_tally4_int_t ldw,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dW, magma_tally4_int_t lddw);

magma_tally4_int_t
magma_tally4_clatrd2(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4FloatComplex *A,  magma_tally4_int_t lda,
    float *e, magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *W,  magma_tally4_int_t ldw,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dW, magma_tally4_int_t lddw,
    magma_tally4FloatComplex_ptr dwork, magma_tally4_int_t ldwork);

#ifdef COMPLEX
magma_tally4_int_t
magma_tally4_clatrsd(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_diag_t diag, magma_tally4_bool_t normin,
    magma_tally4_int_t n, const magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex lambda,
    magma_tally4FloatComplex *x,
    float *scale, float *cnorm,
    magma_tally4_int_t *info);
#endif

magma_tally4_int_t
magma_tally4_clauum(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cposv(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cpotrf(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cpotri(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cstedx(
    magma_tally4_range_t range, magma_tally4_int_t n, float vl, float vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, float *d, float *e,
    magma_tally4FloatComplex *Z, magma_tally4_int_t ldz,
    float *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4Float_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ctrevc3(
    magma_tally4_side_t side, magma_tally4_vec_t howmany,
    magma_tally4_int_t *select, magma_tally4_int_t n,
    magma_tally4FloatComplex *T,  magma_tally4_int_t ldt,
    magma_tally4FloatComplex *VL, magma_tally4_int_t ldvl,
    magma_tally4FloatComplex *VR, magma_tally4_int_t ldvr,
    magma_tally4_int_t mm, magma_tally4_int_t *mout,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ctrevc3_mt(
    magma_tally4_side_t side, magma_tally4_vec_t howmany,
    magma_tally4_int_t *select, magma_tally4_int_t n,
    magma_tally4FloatComplex *T,  magma_tally4_int_t ldt,
    magma_tally4FloatComplex *VL, magma_tally4_int_t ldvl,
    magma_tally4FloatComplex *VR, magma_tally4_int_t ldvr,
    magma_tally4_int_t mm, magma_tally4_int_t *mout,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ctrtri(
    magma_tally4_uplo_t uplo, magma_tally4_diag_t diag, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cunghr(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cungqr(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cungqr2(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cunmbr(
    magma_tally4_vect_t vect, magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *C, magma_tally4_int_t ldc,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cunmlq(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *C, magma_tally4_int_t ldc,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

// not yet implemented
//magma_tally4_int_t magma_tally4_cunmrq( magma_tally4_side_t side, magma_tally4_trans_t trans,
//                          magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
//                          magma_tally4FloatComplex *A, magma_tally4_int_t lda,
//                          magma_tally4FloatComplex *tau,
//                          magma_tally4FloatComplex *C, magma_tally4_int_t ldc,
//                          magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
//                          magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cunmql(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *C, magma_tally4_int_t ldc,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cunmqr(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *C, magma_tally4_int_t ldc,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cunmtr(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A,    magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *C,    magma_tally4_int_t ldc,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 function definitions / Data on CPU / Multi-GPU (alphabetical order)
*/
magma_tally4_int_t
magma_tally4_cgeev_m(
    magma_tally4_vec_t jobvl, magma_tally4_vec_t jobvr, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    #ifdef COMPLEX
    magma_tally4FloatComplex *w,
    #else
    float *wr, float *wi,
    #endif
    magma_tally4FloatComplex *VL, magma_tally4_int_t ldvl,
    magma_tally4FloatComplex *VR, magma_tally4_int_t ldvr,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgehrd_m(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4FloatComplex *T,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgetrf_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cheevd_m(
    magma_tally4_int_t ngpu,
    magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    float *w,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cheevdx_2stage_m(
    magma_tally4_int_t ngpu,
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cheevdx_m(
    magma_tally4_int_t ngpu,
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chegst_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t itype, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chegvd_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *B, magma_tally4_int_t ldb,
    float *w,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chegvdx_2stage_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *B, magma_tally4_int_t ldb,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chegvdx_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *B, magma_tally4_int_t ldb,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

#ifdef REAL
magma_tally4_int_t
magma_tally4_slaex0_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t n, float *d, float *e,
    float *Q, magma_tally4_int_t ldq,
    float *work, magma_tally4_int_t *iwork,
    magma_tally4_range_t range, float vl, float vu,
    magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_slaex1_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t n, float *d,
    float *Q, magma_tally4_int_t ldq,
    magma_tally4_int_t *indxq, float rho, magma_tally4_int_t cutpnt,
    float *work, magma_tally4_int_t *iwork,
    magma_tally4Float_ptr dwork[],
    magma_tally4_queue_t queues[Magma_tally4MaxGPUs][2],
    magma_tally4_range_t range, float vl, float vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_slaex3_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t k, magma_tally4_int_t n, magma_tally4_int_t n1, float *d,
    float *Q, magma_tally4_int_t ldq, float rho,
    float *dlamda, float *Q2, magma_tally4_int_t *indx,
    magma_tally4_int_t *ctot, float *w, float *s, magma_tally4_int_t *indxq,
    magma_tally4Float_ptr dwork[],
    magma_tally4_queue_t queues[Magma_tally4MaxGPUs][2],
    magma_tally4_range_t range, float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *info);
#endif

magma_tally4_int_t
magma_tally4_clahr2_m(
    magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4_int_t nb,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *T, magma_tally4_int_t ldt,
    magma_tally4FloatComplex *Y, magma_tally4_int_t ldy,
    struct cgehrd_data_tally4 *data );

magma_tally4_int_t
magma_tally4_clahru_m(
    magma_tally4_int_t n, magma_tally4_int_t ihi, magma_tally4_int_t k, magma_tally4_int_t nb,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    struct cgehrd_data_tally4 *data );

magma_tally4_int_t
magma_tally4_cpotrf_m(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cstedx_m(
    magma_tally4_int_t ngpu,
    magma_tally4_range_t range, magma_tally4_int_t n, float vl, float vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, float *d, float *e,
    magma_tally4FloatComplex *Z, magma_tally4_int_t ldz,
    float *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ctrsm_m(
    magma_tally4_int_t ngpu,
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t transa, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *B, magma_tally4_int_t ldb);

magma_tally4_int_t
magma_tally4_cunghr_m(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *T, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cungqr_m(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *T, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cunmqr_m(
    magma_tally4_int_t ngpu,
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex *A,    magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *C,    magma_tally4_int_t ldc,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cunmtr_m(
    magma_tally4_int_t ngpu,
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex *A,    magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *C,    magma_tally4_int_t ldc,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 function definitions / Data on GPU (alphabetical order)
*/
magma_tally4_int_t
magma_tally4_cgegqr_gpu(
    magma_tally4_int_t ikind, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dwork, magma_tally4FloatComplex *work,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgelqf_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgels_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4FloatComplex *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgels3_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4FloatComplex *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqp3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt, magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex_ptr dwork, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqr2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dtau,
    magma_tally4Float_ptr        dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqr2x_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dtau,
    magma_tally4FloatComplex_ptr dT, magma_tally4FloatComplex_ptr ddA,
    magma_tally4Float_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqr2x2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dtau,
    magma_tally4FloatComplex_ptr dT, magma_tally4FloatComplex_ptr ddA,
    magma_tally4Float_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqr2x3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dtau,
    magma_tally4FloatComplex_ptr dT,
    magma_tally4FloatComplex_ptr ddA,
    magma_tally4Float_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqr2x4_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dtau,
    magma_tally4FloatComplex_ptr dT, magma_tally4FloatComplex_ptr ddA,
    magma_tally4Float_ptr dwork,
    magma_tally4_queue_t queue,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqrf_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex_ptr dT,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqrf2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex *tau,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqrf2_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dlA[], magma_tally4_int_t ldda,
    magma_tally4FloatComplex *tau,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqrf3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex_ptr dT,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqrs_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex_ptr dT,
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4FloatComplex *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgeqrs3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex_ptr dT,
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4FloatComplex *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgerbt_gpu(
    magma_tally4_bool_t gen, magma_tally4_int_t n, magma_tally4_int_t nrhs, 
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda, 
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb, 
    magma_tally4FloatComplex *U, magma_tally4FloatComplex *V,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgessm_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4_int_t ib,
    magma_tally4_int_t *ipiv,
    magma_tally4FloatComplex_ptr dL1, magma_tally4_int_t lddl1,
    magma_tally4FloatComplex_ptr dL,  magma_tally4_int_t lddl,
    magma_tally4FloatComplex_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgesv_gpu(
    magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgesv_nopiv_gpu( 
    magma_tally4_int_t n, magma_tally4_int_t nrhs, 
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb, 
                 magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgetf2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgetrf_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgetrf_incpiv_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t ib,
    magma_tally4FloatComplex    *hA, magma_tally4_int_t ldha,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex    *hL, magma_tally4_int_t ldhl,
    magma_tally4FloatComplex_ptr dL, magma_tally4_int_t lddl,
    magma_tally4_int_t *ipiv,
    magma_tally4FloatComplex_ptr dwork, magma_tally4_int_t lddwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgetrf_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr d_lA[], magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgetrf_nopiv_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgetrf2_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t offset,
    magma_tally4FloatComplex_ptr d_lAT[], magma_tally4_int_t lddat, magma_tally4_int_t *ipiv,
    magma_tally4FloatComplex_ptr d_lAP[],
    magma_tally4FloatComplex *W, magma_tally4_int_t ldw,
    magma_tally4_queue_t queues[][2],
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgetri_gpu(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4FloatComplex_ptr dwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgetrs_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cgetrs_nopiv_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cheevd_gpu(
    magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    float *w,
    magma_tally4FloatComplex *wA,  magma_tally4_int_t ldwa,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cheevdx_gpu(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    float vl, float vu,
    magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    magma_tally4FloatComplex *wA,  magma_tally4_int_t ldwa,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally4_int_t
magma_tally4_cheevr_gpu(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    float vl, float vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, float abstol, magma_tally4_int_t *m,
    float *w,
    magma_tally4FloatComplex_ptr dZ, magma_tally4_int_t lddz,
    magma_tally4_int_t *isuppz,
    magma_tally4FloatComplex *wA, magma_tally4_int_t ldwa,
    magma_tally4FloatComplex *wZ, magma_tally4_int_t ldwz,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    float *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cheevx_gpu(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    float abstol, magma_tally4_int_t *m,
    float *w,
    magma_tally4FloatComplex_ptr dZ, magma_tally4_int_t lddz,
    magma_tally4FloatComplex *wA, magma_tally4_int_t ldwa,
    magma_tally4FloatComplex *wZ, magma_tally4_int_t ldwz,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    float *rwork, magma_tally4_int_t *iwork,
    magma_tally4_int_t *ifail,
    magma_tally4_int_t *info);
#endif  // COMPLEX

magma_tally4_int_t
magma_tally4_chegst_gpu(
    magma_tally4_int_t itype, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chetrd_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    float *d, float *e, magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *wA,  magma_tally4_int_t ldwa,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chetrd_he2hb_mgpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4FloatComplex_ptr dAmgpu[], magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dTmgpu[], magma_tally4_int_t lddt,
    magma_tally4_int_t ngpu, magma_tally4_int_t distblk,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chetrd_he2hb_mgpu_spec(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4FloatComplex_ptr dAmgpu[], magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dTmgpu[], magma_tally4_int_t lddt,
    magma_tally4_int_t ngpu, magma_tally4_int_t distblk,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chetrd_mgpu(
    magma_tally4_int_t ngpu, magma_tally4_int_t nqueue,
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    float *d, float *e, magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chetrd2_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    float *d, float *e, magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *wA,  magma_tally4_int_t ldwa,
    magma_tally4FloatComplex *work, magma_tally4_int_t lwork,
    magma_tally4FloatComplex_ptr dwork, magma_tally4_int_t ldwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_chetrf_nopiv_tally4_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_clabrd_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4FloatComplex     *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    float *d, float *e, magma_tally4FloatComplex *tauq, magma_tally4FloatComplex *taup,
    magma_tally4FloatComplex     *X, magma_tally4_int_t ldx,
    magma_tally4FloatComplex_ptr dX, magma_tally4_int_t lddx,
    magma_tally4FloatComplex     *Y, magma_tally4_int_t ldy,
    magma_tally4FloatComplex_ptr dY, magma_tally4_int_t lddy );

magma_tally4_int_t
magma_tally4_claqps_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t offset,
    magma_tally4_int_t nb, magma_tally4_int_t *kb,
    magma_tally4FloatComplex_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt, magma_tally4FloatComplex *tau,
    float *vn1, float *vn2,
    magma_tally4FloatComplex_ptr dauxv,
    magma_tally4FloatComplex_ptr dF, magma_tally4_int_t lddf);

magma_tally4_int_t
magma_tally4_claqps2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t offset,
    magma_tally4_int_t nb, magma_tally4_int_t *kb,
    magma_tally4FloatComplex_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt,
    magma_tally4FloatComplex_ptr dtau,
    magma_tally4Float_ptr dvn1, magma_tally4Float_ptr dvn2,
    magma_tally4FloatComplex_ptr dauxv,
    magma_tally4FloatComplex_ptr dF, magma_tally4_int_t lddf);

magma_tally4_int_t
magma_tally4_claqps3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t offset,
    magma_tally4_int_t nb, magma_tally4_int_t *kb,
    magma_tally4FloatComplex_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt,
    magma_tally4FloatComplex_ptr dtau,
    magma_tally4Float_ptr dvn1, magma_tally4Float_ptr dvn2,
    magma_tally4FloatComplex_ptr dauxv,
    magma_tally4FloatComplex_ptr dF, magma_tally4_int_t lddf);

magma_tally4_int_t
magma_tally4_clarf_gpu(
    magma_tally4_int_t m,  magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dv, magma_tally4FloatComplex_const_ptr dtau,
    magma_tally4FloatComplex_ptr dC,  magma_tally4_int_t lddc);

magma_tally4_int_t
magma_tally4_clarfb_gpu(
    magma_tally4_side_t side, magma_tally4_trans_t trans, magma_tally4_direct_t direct, magma_tally4_storev_t storev,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex_const_ptr dV, magma_tally4_int_t lddv,
    magma_tally4FloatComplex_const_ptr dT, magma_tally4_int_t lddt,
    magma_tally4FloatComplex_ptr dC,       magma_tally4_int_t lddc,
    magma_tally4FloatComplex_ptr dwork,    magma_tally4_int_t ldwork );

magma_tally4_int_t
magma_tally4_clarfb_gpu_gemm(
    magma_tally4_side_t side, magma_tally4_trans_t trans, magma_tally4_direct_t direct, magma_tally4_storev_t storev,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex_const_ptr dV, magma_tally4_int_t lddv,
    magma_tally4FloatComplex_const_ptr dT, magma_tally4_int_t lddt,
    magma_tally4FloatComplex_ptr dC,       magma_tally4_int_t lddc,
    magma_tally4FloatComplex_ptr dwork,    magma_tally4_int_t ldwork,
    magma_tally4FloatComplex_ptr dworkvt,  magma_tally4_int_t ldworkvt);

magma_tally4_int_t
magma_tally4_clarfb2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex_const_ptr dV, magma_tally4_int_t lddv,
    magma_tally4FloatComplex_const_ptr dT, magma_tally4_int_t lddt,
    magma_tally4FloatComplex_ptr dC,       magma_tally4_int_t lddc,
    magma_tally4FloatComplex_ptr dwork,    magma_tally4_int_t ldwork );

magma_tally4_int_t
magma_tally4_clatrd_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t nb0,
    magma_tally4FloatComplex *A,  magma_tally4_int_t lda,
    float *e, magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex    *W,       magma_tally4_int_t ldw,
    magma_tally4FloatComplex_ptr dA[],    magma_tally4_int_t ldda, magma_tally4_int_t offset,
    magma_tally4FloatComplex_ptr dW[],    magma_tally4_int_t lddw,
    magma_tally4FloatComplex    *hwork,   magma_tally4_int_t lhwork,
    magma_tally4FloatComplex_ptr dwork[], magma_tally4_int_t ldwork,
    magma_tally4_queue_t queues[] );

magma_tally4_int_t
magma_tally4_clauum_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cposv_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cpotf2_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cpotrf_gpu(
    magma_tally4_uplo_t uplo,  magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cpotrf_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr d_lA[], magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cpotrf_mgpu_right(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr d_lA[], magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cpotrf3_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t off_i, magma_tally4_int_t off_j, magma_tally4_int_t nb,
    magma_tally4FloatComplex_ptr d_lA[],  magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr d_lP[],  magma_tally4_int_t lddp,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda, magma_tally4_int_t h,
    magma_tally4_queue_t queues[][3], magma_tally4_event_t events[][5],
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cpotri_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cpotrs_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cssssm_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m1, magma_tally4_int_t n1,
    magma_tally4_int_t m2, magma_tally4_int_t n2, magma_tally4_int_t k, magma_tally4_int_t ib,
    magma_tally4FloatComplex_ptr dA1, magma_tally4_int_t ldda1,
    magma_tally4FloatComplex_ptr dA2, magma_tally4_int_t ldda2,
    magma_tally4FloatComplex_ptr dL1, magma_tally4_int_t lddl1,
    magma_tally4FloatComplex_ptr dL2, magma_tally4_int_t lddl2,
    magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ctrtri_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_diag_t diag, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ctsqrt_gpu(
    magma_tally4_int_t *m, magma_tally4_int_t *n,
    magma_tally4FloatComplex *A1, magma_tally4FloatComplex *A2, magma_tally4_int_t *lda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex *work, magma_tally4_int_t *lwork,
    magma_tally4FloatComplex_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ctstrf_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t ib, magma_tally4_int_t nb,
    magma_tally4FloatComplex    *hU, magma_tally4_int_t ldhu,
    magma_tally4FloatComplex_ptr dU, magma_tally4_int_t lddu,
    magma_tally4FloatComplex    *hA, magma_tally4_int_t ldha,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex    *hL, magma_tally4_int_t ldhl,
    magma_tally4FloatComplex_ptr dL, magma_tally4_int_t lddl,
    magma_tally4_int_t *ipiv,
    magma_tally4FloatComplex *hwork, magma_tally4_int_t ldhwork,
    magma_tally4FloatComplex_ptr dwork, magma_tally4_int_t lddwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cungqr_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cunmql2_gpu(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex_ptr dC, magma_tally4_int_t lddc,
    magma_tally4FloatComplex *wA, magma_tally4_int_t ldwa,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cunmqr_gpu(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex_ptr dC, magma_tally4_int_t lddc,
    magma_tally4FloatComplex *hwork, magma_tally4_int_t lwork,
    magma_tally4FloatComplex_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cunmqr2_gpu(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex_ptr dC, magma_tally4_int_t lddc,
    magma_tally4FloatComplex    *wA, magma_tally4_int_t ldwa,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_cunmtr_gpu(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex_ptr dC, magma_tally4_int_t lddc,
    magma_tally4FloatComplex    *wA, magma_tally4_int_t ldwa,
    magma_tally4_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 utility function definitions
*/

extern const magma_tally4FloatComplex MAGMA_tally4_C_NAN;
extern const magma_tally4FloatComplex MAGMA_tally4_C_INF;

magma_tally4_int_t
magma_tally4_cnan_inf(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    const magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *cnt_nan,
    magma_tally4_int_t *cnt_inf );

magma_tally4_int_t
magma_tally4_cnan_inf_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *cnt_nan,
    magma_tally4_int_t *cnt_inf );

void magma_tally4_cprint(
    magma_tally4_int_t m, magma_tally4_int_t n,
    const magma_tally4FloatComplex *A, magma_tally4_int_t lda );

void magma_tally4_cprint_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda );

void cpanel_to_q_tally4(
    magma_tally4_uplo_t uplo, magma_tally4_int_t ib,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *work );

void cq_to_panel_tally4(
    magma_tally4_uplo_t uplo, magma_tally4_int_t ib,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4FloatComplex *work );

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif /* MAGMA_tally4_C_H */
