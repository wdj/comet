/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproduct_z.h normal z -> c, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_minproduct_C_H
#define MAGMA_minproduct_C_H

#include "magma_minproduct_types.h"
#include "magma_minproduct_cgehrd_m.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct Auxiliary functions to get the NB used
*/
#ifdef REAL
magma_minproduct_int_t magma_minproduct_get_slaex3_m_nb();       // defined in slaex3_m.cpp
#endif

magma_minproduct_int_t magma_minproduct_get_cpotrf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_cgetrf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_cgetri_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_cgeqp3_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_cgeqrf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_cgeqlf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_cgehrd_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_chetrd_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_chetrf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_chetrf_nopiv_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_cgelqf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_cgebrd_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_chegst_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_cgesvd_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_chegst_nb_m( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_cbulge_nb( magma_minproduct_int_t m, magma_minproduct_int_t nbthreads );
magma_minproduct_int_t magma_minproduct_get_cbulge_nb_mgpu( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_cbulge_get_Vblksiz( magma_minproduct_int_t m, magma_minproduct_int_t nb, magma_minproduct_int_t nbthreads );
magma_minproduct_int_t magma_minproduct_get_cbulge_gcperf();


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct function definitions / Data on CPU (alphabetical order)
*/

#ifdef REAL
// only applicable to real [sd] precisions
void
magma_minproduct_smove_eig(
    magma_minproduct_range_t range, magma_minproduct_int_t n, float *w,
    magma_minproduct_int_t *il, magma_minproduct_int_t *iu, float vl, float vu, magma_minproduct_int_t *m);

// defined in slaex3.cpp
void
magma_minproduct_cvrange(
    magma_minproduct_int_t k, float *d, magma_minproduct_int_t *il, magma_minproduct_int_t *iu, float vl, float vu);

void
magma_minproduct_cirange(
    magma_minproduct_int_t k, magma_minproduct_int_t *indxq, magma_minproduct_int_t *iil, magma_minproduct_int_t *iiu, magma_minproduct_int_t il, magma_minproduct_int_t iu);
#endif

magma_minproduct_int_t
magma_minproduct_cgebrd(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    float *d, float *e,
    magma_minproductFloatComplex *tauq, magma_minproductFloatComplex *taup,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeev(
    magma_minproduct_vec_t jobvl, magma_minproduct_vec_t jobvr, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    #ifdef COMPLEX
    magma_minproductFloatComplex *w,
    #else
    float *wr, float *wi,
    #endif
    magma_minproductFloatComplex *VL, magma_minproduct_int_t ldvl,
    magma_minproductFloatComplex *VR, magma_minproduct_int_t ldvr,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgehrd(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproductFloatComplex_ptr dT,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgehrd2(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgelqf(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A,    magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqlf(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A,    magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqp3(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *jpvt, magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqrf(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqrf_ooc(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqrf4(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A,    magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgesdd(
    magma_minproduct_vec_t jobz, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    float *s,
    magma_minproductFloatComplex *U, magma_minproduct_int_t ldu,
    magma_minproductFloatComplex *VT, magma_minproduct_int_t ldvt,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_minproduct_int_t *iwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgesv(
    magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproductFloatComplex *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgesv_rbt(
    magma_minproduct_bool_t ref, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda, 
    magma_minproductFloatComplex *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgesvd(
    magma_minproduct_vec_t jobu, magma_minproduct_vec_t jobvt, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A,    magma_minproduct_int_t lda, float *s,
    magma_minproductFloatComplex *U,    magma_minproduct_int_t ldu,
    magma_minproductFloatComplex *VT,   magma_minproduct_int_t ldvt,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgetf2_nopiv(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgetrf(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgetrf_nopiv(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgetrf_piv(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t NB,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgetrf2(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cheevd(
    magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    float *w,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cheevdx(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cheevdx_2stage(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_minproduct_int_t
magma_minproduct_cheevr(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    float vl, float vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, float abstol, magma_minproduct_int_t *m,
    float *w,
    magma_minproductFloatComplex *Z, magma_minproduct_int_t ldz,
    magma_minproduct_int_t *isuppz,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    float *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

// no real [sd] precisions available
magma_minproduct_int_t
magma_minproduct_cheevx(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    float vl, float vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, float abstol, magma_minproduct_int_t *m,
    float *w,
    magma_minproductFloatComplex *Z, magma_minproduct_int_t ldz,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    float *rwork, magma_minproduct_int_t *iwork,
    magma_minproduct_int_t *ifail,
    magma_minproduct_int_t *info);
#endif  // COMPLEX

magma_minproduct_int_t
magma_minproduct_chegst(
    magma_minproduct_int_t itype, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chegvd(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *B, magma_minproduct_int_t ldb,
    float *w, magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chegvdx(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n, magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *B, magma_minproduct_int_t ldb,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chegvdx_2stage(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *B, magma_minproduct_int_t ldb,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chegvr(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *B, magma_minproduct_int_t ldb,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    float abstol, magma_minproduct_int_t *m, float *w,
    magma_minproductFloatComplex *Z, magma_minproduct_int_t ldz,
    magma_minproduct_int_t *isuppz, magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    float *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chegvx(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n, magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *B, magma_minproduct_int_t ldb,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    float abstol, magma_minproduct_int_t *m, float *w,
    magma_minproductFloatComplex *Z, magma_minproduct_int_t ldz,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork, float *rwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t *ifail,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chesv(magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
            magma_minproductFloatComplex *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
            magma_minproductFloatComplex *B, magma_minproduct_int_t ldb,
            magma_minproduct_int_t *info );

magma_minproduct_int_t
magma_minproduct_chetrd(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    float *d, float *e, magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chetrf(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chetrf_nopiv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chetrd_hb2st(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t Vblksiz,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    float *d, float *e,
    magma_minproductFloatComplex *V, magma_minproduct_int_t ldv,
    magma_minproductFloatComplex *TAU, magma_minproduct_int_t compT,
    magma_minproductFloatComplex *T, magma_minproduct_int_t ldt);

magma_minproduct_int_t
magma_minproduct_chetrd_he2hb(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproductFloatComplex_ptr dT,
    magma_minproduct_int_t *info);

#ifdef REAL
// only applicable to real [sd] precisions
magma_minproduct_int_t
magma_minproduct_slaex0(
    magma_minproduct_int_t n, float *d, float *e,
    float *Q, magma_minproduct_int_t ldq,
    float *work, magma_minproduct_int_t *iwork,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_range_t range, float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_slaex1(
    magma_minproduct_int_t n, float *d,
    float *Q, magma_minproduct_int_t ldq,
    magma_minproduct_int_t *indxq, float rho, magma_minproduct_int_t cutpnt,
    float *work, magma_minproduct_int_t *iwork,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_range_t range, float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_slaex3(
    magma_minproduct_int_t k, magma_minproduct_int_t n, magma_minproduct_int_t n1, float *d,
    float *Q, magma_minproduct_int_t ldq,
    float rho,
    float *dlamda, float *Q2, magma_minproduct_int_t *indx,
    magma_minproduct_int_t *ctot, float *w, float *s, magma_minproduct_int_t *indxq,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_range_t range, float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *info);
#endif  // REAL

magma_minproduct_int_t
magma_minproduct_clahef_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    magma_minproductFloatComplex    *hA, magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductFloatComplex_ptr dW, magma_minproduct_int_t lddw,
    magma_minproduct_queue_t queues[], magma_minproduct_event_t event[],
    magma_minproduct_int_t *info);

magma_minproduct_int_t
chetrf_nopiv_cpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t ib,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chetrs_nopiv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chesv_nopiv_gpu(
    magma_minproduct_uplo_t uplo,  magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb, 
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_clahr2(
    magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductFloatComplex *A,  magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *T,  magma_minproduct_int_t ldt,
    magma_minproductFloatComplex *Y,  magma_minproduct_int_t ldy);

magma_minproduct_int_t
magma_minproduct_clahru(
    magma_minproduct_int_t n, magma_minproduct_int_t ihi, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    magma_minproductFloatComplex     *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dY, magma_minproduct_int_t lddy,
    magma_minproductFloatComplex_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductFloatComplex_ptr dT,
    magma_minproductFloatComplex_ptr dwork);

#ifdef REAL
magma_minproduct_int_t
magma_minproduct_slaln2(
    magma_minproduct_int_t trans, magma_minproduct_int_t na, magma_minproduct_int_t nw,
    float smin, float ca, const float *A, magma_minproduct_int_t lda,
    float d1, float d2,   const float *B, magma_minproduct_int_t ldb,
    float wr, float wi, float *X, magma_minproduct_int_t ldx,
    float *scale, float *xnorm,
    magma_minproduct_int_t *info);
#endif

magma_minproduct_int_t
magma_minproduct_claqps(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t offset,
    magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    magma_minproductFloatComplex *A,  magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt, magma_minproductFloatComplex *tau, float *vn1, float *vn2,
    magma_minproductFloatComplex *auxv,
    magma_minproductFloatComplex *F,  magma_minproduct_int_t ldf,
    magma_minproductFloatComplex_ptr dF, magma_minproduct_int_t lddf );

#ifdef REAL
magma_minproduct_int_t
magma_minproduct_claqtrsd(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n,
    const float *T, magma_minproduct_int_t ldt,
    float *x,       magma_minproduct_int_t ldx,
    const float *cnorm,
    magma_minproduct_int_t *info);
#endif

magma_minproduct_int_t
magma_minproduct_clatrd(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    float *e, magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *W, magma_minproduct_int_t ldw,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dW, magma_minproduct_int_t lddw);

magma_minproduct_int_t
magma_minproduct_clatrd2(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductFloatComplex *A,  magma_minproduct_int_t lda,
    float *e, magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *W,  magma_minproduct_int_t ldw,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dW, magma_minproduct_int_t lddw,
    magma_minproductFloatComplex_ptr dwork, magma_minproduct_int_t ldwork);

#ifdef COMPLEX
magma_minproduct_int_t
magma_minproduct_clatrsd(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_diag_t diag, magma_minproduct_bool_t normin,
    magma_minproduct_int_t n, const magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex lambda,
    magma_minproductFloatComplex *x,
    float *scale, float *cnorm,
    magma_minproduct_int_t *info);
#endif

magma_minproduct_int_t
magma_minproduct_clauum(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cposv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cpotrf(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cpotri(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cstedx(
    magma_minproduct_range_t range, magma_minproduct_int_t n, float vl, float vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, float *d, float *e,
    magma_minproductFloatComplex *Z, magma_minproduct_int_t ldz,
    float *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ctrevc3(
    magma_minproduct_side_t side, magma_minproduct_vec_t howmany,
    magma_minproduct_int_t *select, magma_minproduct_int_t n,
    magma_minproductFloatComplex *T,  magma_minproduct_int_t ldt,
    magma_minproductFloatComplex *VL, magma_minproduct_int_t ldvl,
    magma_minproductFloatComplex *VR, magma_minproduct_int_t ldvr,
    magma_minproduct_int_t mm, magma_minproduct_int_t *mout,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ctrevc3_mt(
    magma_minproduct_side_t side, magma_minproduct_vec_t howmany,
    magma_minproduct_int_t *select, magma_minproduct_int_t n,
    magma_minproductFloatComplex *T,  magma_minproduct_int_t ldt,
    magma_minproductFloatComplex *VL, magma_minproduct_int_t ldvl,
    magma_minproductFloatComplex *VR, magma_minproduct_int_t ldvr,
    magma_minproduct_int_t mm, magma_minproduct_int_t *mout,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ctrtri(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cunghr(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cungqr(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cungqr2(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cunmbr(
    magma_minproduct_vect_t vect, magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *C, magma_minproduct_int_t ldc,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cunmlq(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *C, magma_minproduct_int_t ldc,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

// not yet implemented
//magma_minproduct_int_t magma_minproduct_cunmrq( magma_minproduct_side_t side, magma_minproduct_trans_t trans,
//                          magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
//                          magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
//                          magma_minproductFloatComplex *tau,
//                          magma_minproductFloatComplex *C, magma_minproduct_int_t ldc,
//                          magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
//                          magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cunmql(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *C, magma_minproduct_int_t ldc,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cunmqr(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *C, magma_minproduct_int_t ldc,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cunmtr(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A,    magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *C,    magma_minproduct_int_t ldc,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct function definitions / Data on CPU / Multi-GPU (alphabetical order)
*/
magma_minproduct_int_t
magma_minproduct_cgeev_m(
    magma_minproduct_vec_t jobvl, magma_minproduct_vec_t jobvr, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    #ifdef COMPLEX
    magma_minproductFloatComplex *w,
    #else
    float *wr, float *wi,
    #endif
    magma_minproductFloatComplex *VL, magma_minproduct_int_t ldvl,
    magma_minproductFloatComplex *VR, magma_minproduct_int_t ldvr,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgehrd_m(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproductFloatComplex *T,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgetrf_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cheevd_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    float *w,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cheevdx_2stage_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cheevdx_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chegst_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t itype, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chegvd_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *B, magma_minproduct_int_t ldb,
    float *w,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chegvdx_2stage_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *B, magma_minproduct_int_t ldb,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chegvdx_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *B, magma_minproduct_int_t ldb,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

#ifdef REAL
magma_minproduct_int_t
magma_minproduct_slaex0_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t n, float *d, float *e,
    float *Q, magma_minproduct_int_t ldq,
    float *work, magma_minproduct_int_t *iwork,
    magma_minproduct_range_t range, float vl, float vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_slaex1_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t n, float *d,
    float *Q, magma_minproduct_int_t ldq,
    magma_minproduct_int_t *indxq, float rho, magma_minproduct_int_t cutpnt,
    float *work, magma_minproduct_int_t *iwork,
    magma_minproductFloat_ptr dwork[],
    magma_minproduct_queue_t queues[Magma_minproductMaxGPUs][2],
    magma_minproduct_range_t range, float vl, float vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_slaex3_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t k, magma_minproduct_int_t n, magma_minproduct_int_t n1, float *d,
    float *Q, magma_minproduct_int_t ldq, float rho,
    float *dlamda, float *Q2, magma_minproduct_int_t *indx,
    magma_minproduct_int_t *ctot, float *w, float *s, magma_minproduct_int_t *indxq,
    magma_minproductFloat_ptr dwork[],
    magma_minproduct_queue_t queues[Magma_minproductMaxGPUs][2],
    magma_minproduct_range_t range, float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *info);
#endif

magma_minproduct_int_t
magma_minproduct_clahr2_m(
    magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *T, magma_minproduct_int_t ldt,
    magma_minproductFloatComplex *Y, magma_minproduct_int_t ldy,
    struct cgehrd_data *data );

magma_minproduct_int_t
magma_minproduct_clahru_m(
    magma_minproduct_int_t n, magma_minproduct_int_t ihi, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    struct cgehrd_data *data );

magma_minproduct_int_t
magma_minproduct_cpotrf_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cstedx_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_range_t range, magma_minproduct_int_t n, float vl, float vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, float *d, float *e,
    magma_minproductFloatComplex *Z, magma_minproduct_int_t ldz,
    float *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ctrsm_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transa, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *B, magma_minproduct_int_t ldb);

magma_minproduct_int_t
magma_minproduct_cunghr_m(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *T, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cungqr_m(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *T, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cunmqr_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex *A,    magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *C,    magma_minproduct_int_t ldc,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cunmtr_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A,    magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *C,    magma_minproduct_int_t ldc,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct function definitions / Data on GPU (alphabetical order)
*/
magma_minproduct_int_t
magma_minproduct_cgegqr_gpu(
    magma_minproduct_int_t ikind, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dwork, magma_minproductFloatComplex *work,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgelqf_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgels_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgels3_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqp3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt, magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex_ptr dwork, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqr2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dtau,
    magma_minproductFloat_ptr        dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqr2x_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dtau,
    magma_minproductFloatComplex_ptr dT, magma_minproductFloatComplex_ptr ddA,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqr2x2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dtau,
    magma_minproductFloatComplex_ptr dT, magma_minproductFloatComplex_ptr ddA,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqr2x3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dtau,
    magma_minproductFloatComplex_ptr dT,
    magma_minproductFloatComplex_ptr ddA,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqr2x4_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dtau,
    magma_minproductFloatComplex_ptr dT, magma_minproductFloatComplex_ptr ddA,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_queue_t queue,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqrf_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex_ptr dT,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqrf2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex *tau,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqrf2_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dlA[], magma_minproduct_int_t ldda,
    magma_minproductFloatComplex *tau,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqrf3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex_ptr dT,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqrs_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex_ptr dT,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgeqrs3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex_ptr dT,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgerbt_gpu(
    magma_minproduct_bool_t gen, magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda, 
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb, 
    magma_minproductFloatComplex *U, magma_minproductFloatComplex *V,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgessm_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t ib,
    magma_minproduct_int_t *ipiv,
    magma_minproductFloatComplex_ptr dL1, magma_minproduct_int_t lddl1,
    magma_minproductFloatComplex_ptr dL,  magma_minproduct_int_t lddl,
    magma_minproductFloatComplex_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgesv_gpu(
    magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgesv_nopiv_gpu( 
    magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb, 
                 magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgetf2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgetrf_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgetrf_incpiv_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t ib,
    magma_minproductFloatComplex    *hA, magma_minproduct_int_t ldha,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex    *hL, magma_minproduct_int_t ldhl,
    magma_minproductFloatComplex_ptr dL, magma_minproduct_int_t lddl,
    magma_minproduct_int_t *ipiv,
    magma_minproductFloatComplex_ptr dwork, magma_minproduct_int_t lddwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgetrf_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr d_lA[], magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgetrf_nopiv_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgetrf2_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t offset,
    magma_minproductFloatComplex_ptr d_lAT[], magma_minproduct_int_t lddat, magma_minproduct_int_t *ipiv,
    magma_minproductFloatComplex_ptr d_lAP[],
    magma_minproductFloatComplex *W, magma_minproduct_int_t ldw,
    magma_minproduct_queue_t queues[][2],
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgetri_gpu(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductFloatComplex_ptr dwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgetrs_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cgetrs_nopiv_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cheevd_gpu(
    magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    float *w,
    magma_minproductFloatComplex *wA,  magma_minproduct_int_t ldwa,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cheevdx_gpu(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    float vl, float vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    magma_minproductFloatComplex *wA,  magma_minproduct_int_t ldwa,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_minproduct_int_t
magma_minproduct_cheevr_gpu(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    float vl, float vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, float abstol, magma_minproduct_int_t *m,
    float *w,
    magma_minproductFloatComplex_ptr dZ, magma_minproduct_int_t lddz,
    magma_minproduct_int_t *isuppz,
    magma_minproductFloatComplex *wA, magma_minproduct_int_t ldwa,
    magma_minproductFloatComplex *wZ, magma_minproduct_int_t ldwz,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    float *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cheevx_gpu(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    float abstol, magma_minproduct_int_t *m,
    float *w,
    magma_minproductFloatComplex_ptr dZ, magma_minproduct_int_t lddz,
    magma_minproductFloatComplex *wA, magma_minproduct_int_t ldwa,
    magma_minproductFloatComplex *wZ, magma_minproduct_int_t ldwz,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    float *rwork, magma_minproduct_int_t *iwork,
    magma_minproduct_int_t *ifail,
    magma_minproduct_int_t *info);
#endif  // COMPLEX

magma_minproduct_int_t
magma_minproduct_chegst_gpu(
    magma_minproduct_int_t itype, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chetrd_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    float *d, float *e, magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *wA,  magma_minproduct_int_t ldwa,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chetrd_he2hb_mgpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproductFloatComplex_ptr dAmgpu[], magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dTmgpu[], magma_minproduct_int_t lddt,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t distblk,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chetrd_he2hb_mgpu_spec(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproductFloatComplex_ptr dAmgpu[], magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dTmgpu[], magma_minproduct_int_t lddt,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t distblk,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chetrd_mgpu(
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nqueue,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    float *d, float *e, magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chetrd2_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    float *d, float *e, magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *wA,  magma_minproduct_int_t ldwa,
    magma_minproductFloatComplex *work, magma_minproduct_int_t lwork,
    magma_minproductFloatComplex_ptr dwork, magma_minproduct_int_t ldwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_chetrf_nopiv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_clabrd_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductFloatComplex     *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    float *d, float *e, magma_minproductFloatComplex *tauq, magma_minproductFloatComplex *taup,
    magma_minproductFloatComplex     *X, magma_minproduct_int_t ldx,
    magma_minproductFloatComplex_ptr dX, magma_minproduct_int_t lddx,
    magma_minproductFloatComplex     *Y, magma_minproduct_int_t ldy,
    magma_minproductFloatComplex_ptr dY, magma_minproduct_int_t lddy );

magma_minproduct_int_t
magma_minproduct_claqps_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t offset,
    magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    magma_minproductFloatComplex_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt, magma_minproductFloatComplex *tau,
    float *vn1, float *vn2,
    magma_minproductFloatComplex_ptr dauxv,
    magma_minproductFloatComplex_ptr dF, magma_minproduct_int_t lddf);

magma_minproduct_int_t
magma_minproduct_claqps2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t offset,
    magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    magma_minproductFloatComplex_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt,
    magma_minproductFloatComplex_ptr dtau,
    magma_minproductFloat_ptr dvn1, magma_minproductFloat_ptr dvn2,
    magma_minproductFloatComplex_ptr dauxv,
    magma_minproductFloatComplex_ptr dF, magma_minproduct_int_t lddf);

magma_minproduct_int_t
magma_minproduct_claqps3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t offset,
    magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    magma_minproductFloatComplex_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt,
    magma_minproductFloatComplex_ptr dtau,
    magma_minproductFloat_ptr dvn1, magma_minproductFloat_ptr dvn2,
    magma_minproductFloatComplex_ptr dauxv,
    magma_minproductFloatComplex_ptr dF, magma_minproduct_int_t lddf);

magma_minproduct_int_t
magma_minproduct_clarf_gpu(
    magma_minproduct_int_t m,  magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dv, magma_minproductFloatComplex_const_ptr dtau,
    magma_minproductFloatComplex_ptr dC,  magma_minproduct_int_t lddc);

magma_minproduct_int_t
magma_minproduct_clarfb_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_direct_t direct, magma_minproduct_storev_t storev,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex_const_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductFloatComplex_const_ptr dT, magma_minproduct_int_t lddt,
    magma_minproductFloatComplex_ptr dC,       magma_minproduct_int_t lddc,
    magma_minproductFloatComplex_ptr dwork,    magma_minproduct_int_t ldwork );

magma_minproduct_int_t
magma_minproduct_clarfb_gpu_gemm(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_direct_t direct, magma_minproduct_storev_t storev,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex_const_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductFloatComplex_const_ptr dT, magma_minproduct_int_t lddt,
    magma_minproductFloatComplex_ptr dC,       magma_minproduct_int_t lddc,
    magma_minproductFloatComplex_ptr dwork,    magma_minproduct_int_t ldwork,
    magma_minproductFloatComplex_ptr dworkvt,  magma_minproduct_int_t ldworkvt);

magma_minproduct_int_t
magma_minproduct_clarfb2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex_const_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductFloatComplex_const_ptr dT, magma_minproduct_int_t lddt,
    magma_minproductFloatComplex_ptr dC,       magma_minproduct_int_t lddc,
    magma_minproductFloatComplex_ptr dwork,    magma_minproduct_int_t ldwork );

magma_minproduct_int_t
magma_minproduct_clatrd_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t nb0,
    magma_minproductFloatComplex *A,  magma_minproduct_int_t lda,
    float *e, magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex    *W,       magma_minproduct_int_t ldw,
    magma_minproductFloatComplex_ptr dA[],    magma_minproduct_int_t ldda, magma_minproduct_int_t offset,
    magma_minproductFloatComplex_ptr dW[],    magma_minproduct_int_t lddw,
    magma_minproductFloatComplex    *hwork,   magma_minproduct_int_t lhwork,
    magma_minproductFloatComplex_ptr dwork[], magma_minproduct_int_t ldwork,
    magma_minproduct_queue_t queues[] );

magma_minproduct_int_t
magma_minproduct_clauum_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cposv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cpotf2_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cpotrf_gpu(
    magma_minproduct_uplo_t uplo,  magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cpotrf_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr d_lA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cpotrf_mgpu_right(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr d_lA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cpotrf3_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t off_i, magma_minproduct_int_t off_j, magma_minproduct_int_t nb,
    magma_minproductFloatComplex_ptr d_lA[],  magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr d_lP[],  magma_minproduct_int_t lddp,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda, magma_minproduct_int_t h,
    magma_minproduct_queue_t queues[][3], magma_minproduct_event_t events[][5],
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cpotri_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cpotrs_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cssssm_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m1, magma_minproduct_int_t n1,
    magma_minproduct_int_t m2, magma_minproduct_int_t n2, magma_minproduct_int_t k, magma_minproduct_int_t ib,
    magma_minproductFloatComplex_ptr dA1, magma_minproduct_int_t ldda1,
    magma_minproductFloatComplex_ptr dA2, magma_minproduct_int_t ldda2,
    magma_minproductFloatComplex_ptr dL1, magma_minproduct_int_t lddl1,
    magma_minproductFloatComplex_ptr dL2, magma_minproduct_int_t lddl2,
    magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ctrtri_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ctsqrt_gpu(
    magma_minproduct_int_t *m, magma_minproduct_int_t *n,
    magma_minproductFloatComplex *A1, magma_minproductFloatComplex *A2, magma_minproduct_int_t *lda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex *work, magma_minproduct_int_t *lwork,
    magma_minproductFloatComplex_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ctstrf_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t ib, magma_minproduct_int_t nb,
    magma_minproductFloatComplex    *hU, magma_minproduct_int_t ldhu,
    magma_minproductFloatComplex_ptr dU, magma_minproduct_int_t lddu,
    magma_minproductFloatComplex    *hA, magma_minproduct_int_t ldha,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex    *hL, magma_minproduct_int_t ldhl,
    magma_minproductFloatComplex_ptr dL, magma_minproduct_int_t lddl,
    magma_minproduct_int_t *ipiv,
    magma_minproductFloatComplex *hwork, magma_minproduct_int_t ldhwork,
    magma_minproductFloatComplex_ptr dwork, magma_minproduct_int_t lddwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cungqr_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cunmql2_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex_ptr dC, magma_minproduct_int_t lddc,
    magma_minproductFloatComplex *wA, magma_minproduct_int_t ldwa,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cunmqr_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex_ptr dC, magma_minproduct_int_t lddc,
    magma_minproductFloatComplex *hwork, magma_minproduct_int_t lwork,
    magma_minproductFloatComplex_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cunmqr2_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex_ptr dC, magma_minproduct_int_t lddc,
    magma_minproductFloatComplex    *wA, magma_minproduct_int_t ldwa,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_cunmtr_gpu(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex *tau,
    magma_minproductFloatComplex_ptr dC, magma_minproduct_int_t lddc,
    magma_minproductFloatComplex    *wA, magma_minproduct_int_t ldwa,
    magma_minproduct_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct utility function definitions
*/

extern const magma_minproductFloatComplex MAGMA_minproduct_C_NAN;
extern const magma_minproductFloatComplex MAGMA_minproduct_C_INF;

magma_minproduct_int_t
magma_minproduct_cnan_inf(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    const magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *cnt_nan,
    magma_minproduct_int_t *cnt_inf );

magma_minproduct_int_t
magma_minproduct_cnan_inf_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *cnt_nan,
    magma_minproduct_int_t *cnt_inf );

void magma_minproduct_cprint(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const magma_minproductFloatComplex *A, magma_minproduct_int_t lda );

void magma_minproduct_cprint_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr dA, magma_minproduct_int_t ldda );

void cpanel_to_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t ib,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *work );

void cq_to_panel(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t ib,
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex *work );

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif /* MAGMA_minproduct_C_H */
