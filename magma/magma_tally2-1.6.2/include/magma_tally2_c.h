/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally2_z.h normal z -> c, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally2_C_H
#define MAGMA_tally2_C_H

#include "magma_tally2_types.h"
#include "magma_tally2_cgehrd_m.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 Auxiliary functions to get the NB used
*/
#ifdef REAL
magma_tally2_int_t magma_tally2_get_slaex3_m_nb();       // defined in slaex3_m.cpp
#endif

magma_tally2_int_t magma_tally2_get_cpotrf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_cgetrf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_cgetri_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_cgeqp3_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_cgeqrf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_cgeqlf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_cgehrd_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_chetrd_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_chetrf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_chetrf_nopiv_tally2_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_cgelqf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_cgebrd_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_chegst_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_cgesvd_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_chegst_nb_m( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_cbulge_nb( magma_tally2_int_t m, magma_tally2_int_t nbthreads );
magma_tally2_int_t magma_tally2_get_cbulge_nb_mgpu( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_cbulge_get_Vblksiz( magma_tally2_int_t m, magma_tally2_int_t nb, magma_tally2_int_t nbthreads );
magma_tally2_int_t magma_tally2_get_cbulge_gcperf();


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 function definitions / Data on CPU (alphabetical order)
*/

#ifdef REAL
// only applicable to real [sd] precisions
void
magma_tally2_smove_eig(
    magma_tally2_range_t range, magma_tally2_int_t n, float *w,
    magma_tally2_int_t *il, magma_tally2_int_t *iu, float vl, float vu, magma_tally2_int_t *m);

// defined in slaex3.cpp
void
magma_tally2_cvrange(
    magma_tally2_int_t k, float *d, magma_tally2_int_t *il, magma_tally2_int_t *iu, float vl, float vu);

void
magma_tally2_cirange(
    magma_tally2_int_t k, magma_tally2_int_t *indxq, magma_tally2_int_t *iil, magma_tally2_int_t *iiu, magma_tally2_int_t il, magma_tally2_int_t iu);
#endif

magma_tally2_int_t
magma_tally2_cgebrd(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    float *d, float *e,
    magma_tally2FloatComplex *tauq, magma_tally2FloatComplex *taup,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeev(
    magma_tally2_vec_t jobvl, magma_tally2_vec_t jobvr, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    #ifdef COMPLEX
    magma_tally2FloatComplex *w,
    #else
    float *wr, float *wi,
    #endif
    magma_tally2FloatComplex *VL, magma_tally2_int_t ldvl,
    magma_tally2FloatComplex *VR, magma_tally2_int_t ldvr,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgehrd(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2FloatComplex_ptr dT,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgehrd2(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgelqf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex *A,    magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqlf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex *A,    magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqp3(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *jpvt, magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqrf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqrf_ooc(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqrf4(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex *A,    magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgesdd(
    magma_tally2_vec_t jobz, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    float *s,
    magma_tally2FloatComplex *U, magma_tally2_int_t ldu,
    magma_tally2FloatComplex *VT, magma_tally2_int_t ldvt,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally2_int_t *iwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgesv(
    magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2FloatComplex *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgesv_rbt(
    magma_tally2_bool_t ref, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda, 
    magma_tally2FloatComplex *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgesvd(
    magma_tally2_vec_t jobu, magma_tally2_vec_t jobvt, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex *A,    magma_tally2_int_t lda, float *s,
    magma_tally2FloatComplex *U,    magma_tally2_int_t ldu,
    magma_tally2FloatComplex *VT,   magma_tally2_int_t ldvt,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgetf2_nopiv(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgetrf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgetrf_nopiv(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgetrf_piv(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t NB,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgetrf2(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cheevd(
    magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    float *w,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cheevdx(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cheevdx_2stage(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally2_int_t
magma_tally2_cheevr(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    float vl, float vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, float abstol, magma_tally2_int_t *m,
    float *w,
    magma_tally2FloatComplex *Z, magma_tally2_int_t ldz,
    magma_tally2_int_t *isuppz,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    float *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

// no real [sd] precisions available
magma_tally2_int_t
magma_tally2_cheevx(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    float vl, float vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, float abstol, magma_tally2_int_t *m,
    float *w,
    magma_tally2FloatComplex *Z, magma_tally2_int_t ldz,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    float *rwork, magma_tally2_int_t *iwork,
    magma_tally2_int_t *ifail,
    magma_tally2_int_t *info);
#endif  // COMPLEX

magma_tally2_int_t
magma_tally2_chegst(
    magma_tally2_int_t itype, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chegvd(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *B, magma_tally2_int_t ldb,
    float *w, magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chegvdx(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n, magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *B, magma_tally2_int_t ldb,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chegvdx_2stage(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *B, magma_tally2_int_t ldb,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chegvr(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *B, magma_tally2_int_t ldb,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    float abstol, magma_tally2_int_t *m, float *w,
    magma_tally2FloatComplex *Z, magma_tally2_int_t ldz,
    magma_tally2_int_t *isuppz, magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    float *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chegvx(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n, magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *B, magma_tally2_int_t ldb,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    float abstol, magma_tally2_int_t *m, float *w,
    magma_tally2FloatComplex *Z, magma_tally2_int_t ldz,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork, float *rwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t *ifail,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chesv(magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
            magma_tally2FloatComplex *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
            magma_tally2FloatComplex *B, magma_tally2_int_t ldb,
            magma_tally2_int_t *info );

magma_tally2_int_t
magma_tally2_chetrd(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    float *d, float *e, magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chetrf(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chetrf_nopiv_tally2(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chetrd_hb2st(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t Vblksiz,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    float *d, float *e,
    magma_tally2FloatComplex *V, magma_tally2_int_t ldv,
    magma_tally2FloatComplex *TAU, magma_tally2_int_t compT,
    magma_tally2FloatComplex *T, magma_tally2_int_t ldt);

magma_tally2_int_t
magma_tally2_chetrd_he2hb(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2FloatComplex_ptr dT,
    magma_tally2_int_t *info);

#ifdef REAL
// only applicable to real [sd] precisions
magma_tally2_int_t
magma_tally2_slaex0(
    magma_tally2_int_t n, float *d, float *e,
    float *Q, magma_tally2_int_t ldq,
    float *work, magma_tally2_int_t *iwork,
    magma_tally2Float_ptr dwork,
    magma_tally2_range_t range, float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_slaex1(
    magma_tally2_int_t n, float *d,
    float *Q, magma_tally2_int_t ldq,
    magma_tally2_int_t *indxq, float rho, magma_tally2_int_t cutpnt,
    float *work, magma_tally2_int_t *iwork,
    magma_tally2Float_ptr dwork,
    magma_tally2_range_t range, float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_slaex3(
    magma_tally2_int_t k, magma_tally2_int_t n, magma_tally2_int_t n1, float *d,
    float *Q, magma_tally2_int_t ldq,
    float rho,
    float *dlamda, float *Q2, magma_tally2_int_t *indx,
    magma_tally2_int_t *ctot, float *w, float *s, magma_tally2_int_t *indxq,
    magma_tally2Float_ptr dwork,
    magma_tally2_range_t range, float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *info);
#endif  // REAL

magma_tally2_int_t
magma_tally2_clahef_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t *kb,
    magma_tally2FloatComplex    *hA, magma_tally2_int_t lda,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2FloatComplex_ptr dW, magma_tally2_int_t lddw,
    magma_tally2_queue_t queues[], magma_tally2_event_t event[],
    magma_tally2_int_t *info);

magma_tally2_int_t
chetrf_nopiv_tally2_cpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t ib,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chetrs_nopiv_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chesv_nopiv_gpu(
    magma_tally2_uplo_t uplo,  magma_tally2_int_t n, magma_tally2_int_t nrhs, 
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb, 
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_clahr2(
    magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2_int_t nb,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dV, magma_tally2_int_t lddv,
    magma_tally2FloatComplex *A,  magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *T,  magma_tally2_int_t ldt,
    magma_tally2FloatComplex *Y,  magma_tally2_int_t ldy);

magma_tally2_int_t
magma_tally2_clahru(
    magma_tally2_int_t n, magma_tally2_int_t ihi, magma_tally2_int_t k, magma_tally2_int_t nb,
    magma_tally2FloatComplex     *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dY, magma_tally2_int_t lddy,
    magma_tally2FloatComplex_ptr dV, magma_tally2_int_t lddv,
    magma_tally2FloatComplex_ptr dT,
    magma_tally2FloatComplex_ptr dwork);

#ifdef REAL
magma_tally2_int_t
magma_tally2_slaln2(
    magma_tally2_int_t trans, magma_tally2_int_t na, magma_tally2_int_t nw,
    float smin, float ca, const float *A, magma_tally2_int_t lda,
    float d1, float d2,   const float *B, magma_tally2_int_t ldb,
    float wr, float wi, float *X, magma_tally2_int_t ldx,
    float *scale, float *xnorm,
    magma_tally2_int_t *info);
#endif

magma_tally2_int_t
magma_tally2_claqps(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t offset,
    magma_tally2_int_t nb, magma_tally2_int_t *kb,
    magma_tally2FloatComplex *A,  magma_tally2_int_t lda,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt, magma_tally2FloatComplex *tau, float *vn1, float *vn2,
    magma_tally2FloatComplex *auxv,
    magma_tally2FloatComplex *F,  magma_tally2_int_t ldf,
    magma_tally2FloatComplex_ptr dF, magma_tally2_int_t lddf );

#ifdef REAL
magma_tally2_int_t
magma_tally2_claqtrsd(
    magma_tally2_trans_t trans, magma_tally2_int_t n,
    const float *T, magma_tally2_int_t ldt,
    float *x,       magma_tally2_int_t ldx,
    const float *cnorm,
    magma_tally2_int_t *info);
#endif

magma_tally2_int_t
magma_tally2_clatrd(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    float *e, magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *W, magma_tally2_int_t ldw,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dW, magma_tally2_int_t lddw);

magma_tally2_int_t
magma_tally2_clatrd2(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2FloatComplex *A,  magma_tally2_int_t lda,
    float *e, magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *W,  magma_tally2_int_t ldw,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dW, magma_tally2_int_t lddw,
    magma_tally2FloatComplex_ptr dwork, magma_tally2_int_t ldwork);

#ifdef COMPLEX
magma_tally2_int_t
magma_tally2_clatrsd(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_diag_t diag, magma_tally2_bool_t normin,
    magma_tally2_int_t n, const magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex lambda,
    magma_tally2FloatComplex *x,
    float *scale, float *cnorm,
    magma_tally2_int_t *info);
#endif

magma_tally2_int_t
magma_tally2_clauum(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cposv(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cpotrf(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cpotri(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cstedx(
    magma_tally2_range_t range, magma_tally2_int_t n, float vl, float vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, float *d, float *e,
    magma_tally2FloatComplex *Z, magma_tally2_int_t ldz,
    float *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2Float_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ctrevc3(
    magma_tally2_side_t side, magma_tally2_vec_t howmany,
    magma_tally2_int_t *select, magma_tally2_int_t n,
    magma_tally2FloatComplex *T,  magma_tally2_int_t ldt,
    magma_tally2FloatComplex *VL, magma_tally2_int_t ldvl,
    magma_tally2FloatComplex *VR, magma_tally2_int_t ldvr,
    magma_tally2_int_t mm, magma_tally2_int_t *mout,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ctrevc3_mt(
    magma_tally2_side_t side, magma_tally2_vec_t howmany,
    magma_tally2_int_t *select, magma_tally2_int_t n,
    magma_tally2FloatComplex *T,  magma_tally2_int_t ldt,
    magma_tally2FloatComplex *VL, magma_tally2_int_t ldvl,
    magma_tally2FloatComplex *VR, magma_tally2_int_t ldvr,
    magma_tally2_int_t mm, magma_tally2_int_t *mout,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ctrtri(
    magma_tally2_uplo_t uplo, magma_tally2_diag_t diag, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cunghr(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cungqr(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cungqr2(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cunmbr(
    magma_tally2_vect_t vect, magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *C, magma_tally2_int_t ldc,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cunmlq(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *C, magma_tally2_int_t ldc,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

// not yet implemented
//magma_tally2_int_t magma_tally2_cunmrq( magma_tally2_side_t side, magma_tally2_trans_t trans,
//                          magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
//                          magma_tally2FloatComplex *A, magma_tally2_int_t lda,
//                          magma_tally2FloatComplex *tau,
//                          magma_tally2FloatComplex *C, magma_tally2_int_t ldc,
//                          magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
//                          magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cunmql(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *C, magma_tally2_int_t ldc,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cunmqr(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *C, magma_tally2_int_t ldc,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cunmtr(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex *A,    magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *C,    magma_tally2_int_t ldc,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 function definitions / Data on CPU / Multi-GPU (alphabetical order)
*/
magma_tally2_int_t
magma_tally2_cgeev_m(
    magma_tally2_vec_t jobvl, magma_tally2_vec_t jobvr, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    #ifdef COMPLEX
    magma_tally2FloatComplex *w,
    #else
    float *wr, float *wi,
    #endif
    magma_tally2FloatComplex *VL, magma_tally2_int_t ldvl,
    magma_tally2FloatComplex *VR, magma_tally2_int_t ldvr,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgehrd_m(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2FloatComplex *T,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgetrf_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cheevd_m(
    magma_tally2_int_t ngpu,
    magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    float *w,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cheevdx_2stage_m(
    magma_tally2_int_t ngpu,
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cheevdx_m(
    magma_tally2_int_t ngpu,
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chegst_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chegvd_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *B, magma_tally2_int_t ldb,
    float *w,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chegvdx_2stage_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *B, magma_tally2_int_t ldb,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chegvdx_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *B, magma_tally2_int_t ldb,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

#ifdef REAL
magma_tally2_int_t
magma_tally2_slaex0_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t n, float *d, float *e,
    float *Q, magma_tally2_int_t ldq,
    float *work, magma_tally2_int_t *iwork,
    magma_tally2_range_t range, float vl, float vu,
    magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_slaex1_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t n, float *d,
    float *Q, magma_tally2_int_t ldq,
    magma_tally2_int_t *indxq, float rho, magma_tally2_int_t cutpnt,
    float *work, magma_tally2_int_t *iwork,
    magma_tally2Float_ptr dwork[],
    magma_tally2_queue_t queues[Magma_tally2MaxGPUs][2],
    magma_tally2_range_t range, float vl, float vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_slaex3_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t k, magma_tally2_int_t n, magma_tally2_int_t n1, float *d,
    float *Q, magma_tally2_int_t ldq, float rho,
    float *dlamda, float *Q2, magma_tally2_int_t *indx,
    magma_tally2_int_t *ctot, float *w, float *s, magma_tally2_int_t *indxq,
    magma_tally2Float_ptr dwork[],
    magma_tally2_queue_t queues[Magma_tally2MaxGPUs][2],
    magma_tally2_range_t range, float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *info);
#endif

magma_tally2_int_t
magma_tally2_clahr2_m(
    magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2_int_t nb,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *T, magma_tally2_int_t ldt,
    magma_tally2FloatComplex *Y, magma_tally2_int_t ldy,
    struct cgehrd_data_tally2 *data );

magma_tally2_int_t
magma_tally2_clahru_m(
    magma_tally2_int_t n, magma_tally2_int_t ihi, magma_tally2_int_t k, magma_tally2_int_t nb,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    struct cgehrd_data_tally2 *data );

magma_tally2_int_t
magma_tally2_cpotrf_m(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cstedx_m(
    magma_tally2_int_t ngpu,
    magma_tally2_range_t range, magma_tally2_int_t n, float vl, float vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, float *d, float *e,
    magma_tally2FloatComplex *Z, magma_tally2_int_t ldz,
    float *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ctrsm_m(
    magma_tally2_int_t ngpu,
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t transa, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *B, magma_tally2_int_t ldb);

magma_tally2_int_t
magma_tally2_cunghr_m(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *T, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cungqr_m(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *T, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cunmqr_m(
    magma_tally2_int_t ngpu,
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex *A,    magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *C,    magma_tally2_int_t ldc,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cunmtr_m(
    magma_tally2_int_t ngpu,
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex *A,    magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *C,    magma_tally2_int_t ldc,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 function definitions / Data on GPU (alphabetical order)
*/
magma_tally2_int_t
magma_tally2_cgegqr_gpu(
    magma_tally2_int_t ikind, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dwork, magma_tally2FloatComplex *work,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgelqf_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgels_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2FloatComplex *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgels3_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2FloatComplex *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqp3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt, magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex_ptr dwork, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqr2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dtau,
    magma_tally2Float_ptr        dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqr2x_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dtau,
    magma_tally2FloatComplex_ptr dT, magma_tally2FloatComplex_ptr ddA,
    magma_tally2Float_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqr2x2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dtau,
    magma_tally2FloatComplex_ptr dT, magma_tally2FloatComplex_ptr ddA,
    magma_tally2Float_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqr2x3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dtau,
    magma_tally2FloatComplex_ptr dT,
    magma_tally2FloatComplex_ptr ddA,
    magma_tally2Float_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqr2x4_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dtau,
    magma_tally2FloatComplex_ptr dT, magma_tally2FloatComplex_ptr ddA,
    magma_tally2Float_ptr dwork,
    magma_tally2_queue_t queue,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqrf_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex_ptr dT,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqrf2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex *tau,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqrf2_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dlA[], magma_tally2_int_t ldda,
    magma_tally2FloatComplex *tau,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqrf3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex_ptr dT,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqrs_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex_ptr dT,
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2FloatComplex *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgeqrs3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex_ptr dT,
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2FloatComplex *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgerbt_gpu(
    magma_tally2_bool_t gen, magma_tally2_int_t n, magma_tally2_int_t nrhs, 
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda, 
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb, 
    magma_tally2FloatComplex *U, magma_tally2FloatComplex *V,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgessm_gpu(
    magma_tally2_order_t order, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2_int_t ib,
    magma_tally2_int_t *ipiv,
    magma_tally2FloatComplex_ptr dL1, magma_tally2_int_t lddl1,
    magma_tally2FloatComplex_ptr dL,  magma_tally2_int_t lddl,
    magma_tally2FloatComplex_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgesv_gpu(
    magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgesv_nopiv_gpu( 
    magma_tally2_int_t n, magma_tally2_int_t nrhs, 
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb, 
                 magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgetf2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgetrf_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgetrf_incpiv_gpu(
    magma_tally2_order_t order, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t ib,
    magma_tally2FloatComplex    *hA, magma_tally2_int_t ldha,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex    *hL, magma_tally2_int_t ldhl,
    magma_tally2FloatComplex_ptr dL, magma_tally2_int_t lddl,
    magma_tally2_int_t *ipiv,
    magma_tally2FloatComplex_ptr dwork, magma_tally2_int_t lddwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgetrf_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr d_lA[], magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgetrf_nopiv_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgetrf2_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t offset,
    magma_tally2FloatComplex_ptr d_lAT[], magma_tally2_int_t lddat, magma_tally2_int_t *ipiv,
    magma_tally2FloatComplex_ptr d_lAP[],
    magma_tally2FloatComplex *W, magma_tally2_int_t ldw,
    magma_tally2_queue_t queues[][2],
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgetri_gpu(
    magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2FloatComplex_ptr dwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgetrs_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cgetrs_nopiv_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cheevd_gpu(
    magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    float *w,
    magma_tally2FloatComplex *wA,  magma_tally2_int_t ldwa,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cheevdx_gpu(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    float vl, float vu,
    magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    magma_tally2FloatComplex *wA,  magma_tally2_int_t ldwa,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally2_int_t
magma_tally2_cheevr_gpu(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    float vl, float vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, float abstol, magma_tally2_int_t *m,
    float *w,
    magma_tally2FloatComplex_ptr dZ, magma_tally2_int_t lddz,
    magma_tally2_int_t *isuppz,
    magma_tally2FloatComplex *wA, magma_tally2_int_t ldwa,
    magma_tally2FloatComplex *wZ, magma_tally2_int_t ldwz,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    float *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cheevx_gpu(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    float abstol, magma_tally2_int_t *m,
    float *w,
    magma_tally2FloatComplex_ptr dZ, magma_tally2_int_t lddz,
    magma_tally2FloatComplex *wA, magma_tally2_int_t ldwa,
    magma_tally2FloatComplex *wZ, magma_tally2_int_t ldwz,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    float *rwork, magma_tally2_int_t *iwork,
    magma_tally2_int_t *ifail,
    magma_tally2_int_t *info);
#endif  // COMPLEX

magma_tally2_int_t
magma_tally2_chegst_gpu(
    magma_tally2_int_t itype, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chetrd_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    float *d, float *e, magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *wA,  magma_tally2_int_t ldwa,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chetrd_he2hb_mgpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2FloatComplex_ptr dAmgpu[], magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dTmgpu[], magma_tally2_int_t lddt,
    magma_tally2_int_t ngpu, magma_tally2_int_t distblk,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chetrd_he2hb_mgpu_spec(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2FloatComplex_ptr dAmgpu[], magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dTmgpu[], magma_tally2_int_t lddt,
    magma_tally2_int_t ngpu, magma_tally2_int_t distblk,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chetrd_mgpu(
    magma_tally2_int_t ngpu, magma_tally2_int_t nqueue,
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    float *d, float *e, magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chetrd2_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    float *d, float *e, magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *wA,  magma_tally2_int_t ldwa,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2FloatComplex_ptr dwork, magma_tally2_int_t ldwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_chetrf_nopiv_tally2_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_clabrd_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2FloatComplex     *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    float *d, float *e, magma_tally2FloatComplex *tauq, magma_tally2FloatComplex *taup,
    magma_tally2FloatComplex     *X, magma_tally2_int_t ldx,
    magma_tally2FloatComplex_ptr dX, magma_tally2_int_t lddx,
    magma_tally2FloatComplex     *Y, magma_tally2_int_t ldy,
    magma_tally2FloatComplex_ptr dY, magma_tally2_int_t lddy );

magma_tally2_int_t
magma_tally2_claqps_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t offset,
    magma_tally2_int_t nb, magma_tally2_int_t *kb,
    magma_tally2FloatComplex_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt, magma_tally2FloatComplex *tau,
    float *vn1, float *vn2,
    magma_tally2FloatComplex_ptr dauxv,
    magma_tally2FloatComplex_ptr dF, magma_tally2_int_t lddf);

magma_tally2_int_t
magma_tally2_claqps2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t offset,
    magma_tally2_int_t nb, magma_tally2_int_t *kb,
    magma_tally2FloatComplex_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt,
    magma_tally2FloatComplex_ptr dtau,
    magma_tally2Float_ptr dvn1, magma_tally2Float_ptr dvn2,
    magma_tally2FloatComplex_ptr dauxv,
    magma_tally2FloatComplex_ptr dF, magma_tally2_int_t lddf);

magma_tally2_int_t
magma_tally2_claqps3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t offset,
    magma_tally2_int_t nb, magma_tally2_int_t *kb,
    magma_tally2FloatComplex_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt,
    magma_tally2FloatComplex_ptr dtau,
    magma_tally2Float_ptr dvn1, magma_tally2Float_ptr dvn2,
    magma_tally2FloatComplex_ptr dauxv,
    magma_tally2FloatComplex_ptr dF, magma_tally2_int_t lddf);

magma_tally2_int_t
magma_tally2_clarf_gpu(
    magma_tally2_int_t m,  magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr dv, magma_tally2FloatComplex_const_ptr dtau,
    magma_tally2FloatComplex_ptr dC,  magma_tally2_int_t lddc);

magma_tally2_int_t
magma_tally2_clarfb_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans, magma_tally2_direct_t direct, magma_tally2_storev_t storev,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex_const_ptr dV, magma_tally2_int_t lddv,
    magma_tally2FloatComplex_const_ptr dT, magma_tally2_int_t lddt,
    magma_tally2FloatComplex_ptr dC,       magma_tally2_int_t lddc,
    magma_tally2FloatComplex_ptr dwork,    magma_tally2_int_t ldwork );

magma_tally2_int_t
magma_tally2_clarfb_gpu_gemm(
    magma_tally2_side_t side, magma_tally2_trans_t trans, magma_tally2_direct_t direct, magma_tally2_storev_t storev,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex_const_ptr dV, magma_tally2_int_t lddv,
    magma_tally2FloatComplex_const_ptr dT, magma_tally2_int_t lddt,
    magma_tally2FloatComplex_ptr dC,       magma_tally2_int_t lddc,
    magma_tally2FloatComplex_ptr dwork,    magma_tally2_int_t ldwork,
    magma_tally2FloatComplex_ptr dworkvt,  magma_tally2_int_t ldworkvt);

magma_tally2_int_t
magma_tally2_clarfb2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex_const_ptr dV, magma_tally2_int_t lddv,
    magma_tally2FloatComplex_const_ptr dT, magma_tally2_int_t lddt,
    magma_tally2FloatComplex_ptr dC,       magma_tally2_int_t lddc,
    magma_tally2FloatComplex_ptr dwork,    magma_tally2_int_t ldwork );

magma_tally2_int_t
magma_tally2_clatrd_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t nb0,
    magma_tally2FloatComplex *A,  magma_tally2_int_t lda,
    float *e, magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex    *W,       magma_tally2_int_t ldw,
    magma_tally2FloatComplex_ptr dA[],    magma_tally2_int_t ldda, magma_tally2_int_t offset,
    magma_tally2FloatComplex_ptr dW[],    magma_tally2_int_t lddw,
    magma_tally2FloatComplex    *hwork,   magma_tally2_int_t lhwork,
    magma_tally2FloatComplex_ptr dwork[], magma_tally2_int_t ldwork,
    magma_tally2_queue_t queues[] );

magma_tally2_int_t
magma_tally2_clauum_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cposv_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cpotf2_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cpotrf_gpu(
    magma_tally2_uplo_t uplo,  magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cpotrf_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr d_lA[], magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cpotrf_mgpu_right(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr d_lA[], magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cpotrf3_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t off_i, magma_tally2_int_t off_j, magma_tally2_int_t nb,
    magma_tally2FloatComplex_ptr d_lA[],  magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr d_lP[],  magma_tally2_int_t lddp,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda, magma_tally2_int_t h,
    magma_tally2_queue_t queues[][3], magma_tally2_event_t events[][5],
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cpotri_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cpotrs_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cssssm_gpu(
    magma_tally2_order_t order, magma_tally2_int_t m1, magma_tally2_int_t n1,
    magma_tally2_int_t m2, magma_tally2_int_t n2, magma_tally2_int_t k, magma_tally2_int_t ib,
    magma_tally2FloatComplex_ptr dA1, magma_tally2_int_t ldda1,
    magma_tally2FloatComplex_ptr dA2, magma_tally2_int_t ldda2,
    magma_tally2FloatComplex_ptr dL1, magma_tally2_int_t lddl1,
    magma_tally2FloatComplex_ptr dL2, magma_tally2_int_t lddl2,
    magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ctrtri_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_diag_t diag, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ctsqrt_gpu(
    magma_tally2_int_t *m, magma_tally2_int_t *n,
    magma_tally2FloatComplex *A1, magma_tally2FloatComplex *A2, magma_tally2_int_t *lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t *lwork,
    magma_tally2FloatComplex_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ctstrf_gpu(
    magma_tally2_order_t order, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t ib, magma_tally2_int_t nb,
    magma_tally2FloatComplex    *hU, magma_tally2_int_t ldhu,
    magma_tally2FloatComplex_ptr dU, magma_tally2_int_t lddu,
    magma_tally2FloatComplex    *hA, magma_tally2_int_t ldha,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex    *hL, magma_tally2_int_t ldhl,
    magma_tally2FloatComplex_ptr dL, magma_tally2_int_t lddl,
    magma_tally2_int_t *ipiv,
    magma_tally2FloatComplex *hwork, magma_tally2_int_t ldhwork,
    magma_tally2FloatComplex_ptr dwork, magma_tally2_int_t lddwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cungqr_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cunmql2_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex_ptr dC, magma_tally2_int_t lddc,
    magma_tally2FloatComplex *wA, magma_tally2_int_t ldwa,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cunmqr_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex_ptr dC, magma_tally2_int_t lddc,
    magma_tally2FloatComplex *hwork, magma_tally2_int_t lwork,
    magma_tally2FloatComplex_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cunmqr2_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex_ptr dC, magma_tally2_int_t lddc,
    magma_tally2FloatComplex    *wA, magma_tally2_int_t ldwa,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_cunmtr_gpu(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex_ptr dC, magma_tally2_int_t lddc,
    magma_tally2FloatComplex    *wA, magma_tally2_int_t ldwa,
    magma_tally2_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 utility function definitions
*/

extern const magma_tally2FloatComplex MAGMA_tally2_C_NAN;
extern const magma_tally2FloatComplex MAGMA_tally2_C_INF;

magma_tally2_int_t
magma_tally2_cnan_inf(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    const magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2_int_t *cnt_nan,
    magma_tally2_int_t *cnt_inf );

magma_tally2_int_t
magma_tally2_cnan_inf_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *cnt_nan,
    magma_tally2_int_t *cnt_inf );

void magma_tally2_cprint(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const magma_tally2FloatComplex *A, magma_tally2_int_t lda );

void magma_tally2_cprint_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr dA, magma_tally2_int_t ldda );

void cpanel_to_q_tally2(
    magma_tally2_uplo_t uplo, magma_tally2_int_t ib,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *work );

void cq_to_panel_tally2(
    magma_tally2_uplo_t uplo, magma_tally2_int_t ib,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *work );

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif /* MAGMA_tally2_C_H */
