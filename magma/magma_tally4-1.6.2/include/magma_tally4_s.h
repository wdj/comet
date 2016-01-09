/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally4_z.h normal z -> s, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally4_S_H
#define MAGMA_tally4_S_H

#include "magma_tally4_types.h"
#include "magma_tally4_sgehrd_m.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 Auxiliary functions to get the NB used
*/
#ifdef REAL
magma_tally4_int_t magma_tally4_get_slaex3_m_nb();       // defined in slaex3_m.cpp
#endif

magma_tally4_int_t magma_tally4_get_spotrf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_sgetrf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_sgetri_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_sgeqp3_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_sgeqrf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_sgeqlf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_sgehrd_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_ssytrd_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_ssytrf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_ssytrf_nopiv_tally4_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_sgelqf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_sgebrd_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_ssygst_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_sgesvd_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_ssygst_nb_m( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_sbulge_nb( magma_tally4_int_t m, magma_tally4_int_t nbthreads );
magma_tally4_int_t magma_tally4_get_sbulge_nb_mgpu( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_sbulge_get_Vblksiz( magma_tally4_int_t m, magma_tally4_int_t nb, magma_tally4_int_t nbthreads );
magma_tally4_int_t magma_tally4_get_sbulge_gcperf();


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
magma_tally4_svrange(
    magma_tally4_int_t k, float *d, magma_tally4_int_t *il, magma_tally4_int_t *iu, float vl, float vu);

void
magma_tally4_sirange(
    magma_tally4_int_t k, magma_tally4_int_t *indxq, magma_tally4_int_t *iil, magma_tally4_int_t *iiu, magma_tally4_int_t il, magma_tally4_int_t iu);
#endif

magma_tally4_int_t
magma_tally4_sgebrd(
    magma_tally4_int_t m, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float *d, float *e,
    float *tauq, float *taup,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeev(
    magma_tally4_vec_t jobvl, magma_tally4_vec_t jobvr, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    #ifdef COMPLEX
    float *w,
    #else
    float *wr, float *wi,
    #endif
    float *VL, magma_tally4_int_t ldvl,
    float *VR, magma_tally4_int_t ldvr,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgehrd(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    float *A, magma_tally4_int_t lda,
    float *tau,
    float *work, magma_tally4_int_t lwork,
    magma_tally4Float_ptr dT,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgehrd2(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    float *A, magma_tally4_int_t lda,
    float *tau,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgelqf(
    magma_tally4_int_t m, magma_tally4_int_t n,
    float *A,    magma_tally4_int_t lda,
    float *tau,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqlf(
    magma_tally4_int_t m, magma_tally4_int_t n,
    float *A,    magma_tally4_int_t lda,
    float *tau,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqp3(
    magma_tally4_int_t m, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    magma_tally4_int_t *jpvt, float *tau,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqrf(
    magma_tally4_int_t m, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float *tau,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqrf_ooc(
    magma_tally4_int_t m, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float *tau,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqrf4(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n,
    float *A,    magma_tally4_int_t lda,
    float *tau,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgesdd(
    magma_tally4_vec_t jobz, magma_tally4_int_t m, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float *s,
    float *U, magma_tally4_int_t ldu,
    float *VT, magma_tally4_int_t ldvt,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally4_int_t *iwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgesv(
    magma_tally4_int_t n, magma_tally4_int_t nrhs,
    float *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    float *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgesv_rbt(
    magma_tally4_bool_t ref, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    float *A, magma_tally4_int_t lda, 
    float *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgesvd(
    magma_tally4_vec_t jobu, magma_tally4_vec_t jobvt, magma_tally4_int_t m, magma_tally4_int_t n,
    float *A,    magma_tally4_int_t lda, float *s,
    float *U,    magma_tally4_int_t ldu,
    float *VT,   magma_tally4_int_t ldvt,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgetf2_nopiv(
    magma_tally4_int_t m, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgetrf(
    magma_tally4_int_t m, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgetrf_nopiv(
    magma_tally4_int_t m, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgetrf_piv(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t NB,
    float *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgetrf2(
    magma_tally4_int_t m, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssyevd(
    magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float *w,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssyevdx(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssyevdx_2stage(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally4_int_t
magma_tally4_ssyevr(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float vl, float vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, float abstol, magma_tally4_int_t *m,
    float *w,
    float *Z, magma_tally4_int_t ldz,
    magma_tally4_int_t *isuppz,
    float *work, magma_tally4_int_t lwork,
    float *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

// no real [sd] precisions available
magma_tally4_int_t
magma_tally4_ssyevx(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float vl, float vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, float abstol, magma_tally4_int_t *m,
    float *w,
    float *Z, magma_tally4_int_t ldz,
    float *work, magma_tally4_int_t lwork,
    float *rwork, magma_tally4_int_t *iwork,
    magma_tally4_int_t *ifail,
    magma_tally4_int_t *info);
#endif  // COMPLEX

magma_tally4_int_t
magma_tally4_ssygst(
    magma_tally4_int_t itype, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssygvd(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float *B, magma_tally4_int_t ldb,
    float *w, float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssygvdx(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n, float *A, magma_tally4_int_t lda,
    float *B, magma_tally4_int_t ldb,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssygvdx_2stage(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float *B, magma_tally4_int_t ldb,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssygvr(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float *B, magma_tally4_int_t ldb,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    float abstol, magma_tally4_int_t *m, float *w,
    float *Z, magma_tally4_int_t ldz,
    magma_tally4_int_t *isuppz, float *work, magma_tally4_int_t lwork,
    float *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssygvx(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n, float *A, magma_tally4_int_t lda,
    float *B, magma_tally4_int_t ldb,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    float abstol, magma_tally4_int_t *m, float *w,
    float *Z, magma_tally4_int_t ldz,
    float *work, magma_tally4_int_t lwork, float *rwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t *ifail,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssysv(magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
            float *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
            float *B, magma_tally4_int_t ldb,
            magma_tally4_int_t *info );

magma_tally4_int_t
magma_tally4_ssytrd(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float *d, float *e, float *tau,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssytrf(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssytrf_nopiv_tally4(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssytrd_sb2st(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t Vblksiz,
    float *A, magma_tally4_int_t lda,
    float *d, float *e,
    float *V, magma_tally4_int_t ldv,
    float *TAU, magma_tally4_int_t compT,
    float *T, magma_tally4_int_t ldt);

magma_tally4_int_t
magma_tally4_ssytrd_sy2sb(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    float *A, magma_tally4_int_t lda,
    float *tau,
    float *work, magma_tally4_int_t lwork,
    magma_tally4Float_ptr dT,
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
magma_tally4_slasyf_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t *kb,
    float    *hA, magma_tally4_int_t lda,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4Float_ptr dW, magma_tally4_int_t lddw,
    magma_tally4_queue_t queues[], magma_tally4_event_t event[],
    magma_tally4_int_t *info);

magma_tally4_int_t
ssytrf_nopiv_tally4_cpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t ib,
    float *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssytrs_nopiv_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssysv_nopiv_gpu(
    magma_tally4_uplo_t uplo,  magma_tally4_int_t n, magma_tally4_int_t nrhs, 
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb, 
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_slahr2(
    magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4_int_t nb,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dV, magma_tally4_int_t lddv,
    float *A,  magma_tally4_int_t lda,
    float *tau,
    float *T,  magma_tally4_int_t ldt,
    float *Y,  magma_tally4_int_t ldy);

magma_tally4_int_t
magma_tally4_slahru(
    magma_tally4_int_t n, magma_tally4_int_t ihi, magma_tally4_int_t k, magma_tally4_int_t nb,
    float     *A, magma_tally4_int_t lda,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dY, magma_tally4_int_t lddy,
    magma_tally4Float_ptr dV, magma_tally4_int_t lddv,
    magma_tally4Float_ptr dT,
    magma_tally4Float_ptr dwork);

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
magma_tally4_slaqps(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t offset,
    magma_tally4_int_t nb, magma_tally4_int_t *kb,
    float *A,  magma_tally4_int_t lda,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt, float *tau, float *vn1, float *vn2,
    float *auxv,
    float *F,  magma_tally4_int_t ldf,
    magma_tally4Float_ptr dF, magma_tally4_int_t lddf );

#ifdef REAL
magma_tally4_int_t
magma_tally4_slaqtrsd(
    magma_tally4_trans_t trans, magma_tally4_int_t n,
    const float *T, magma_tally4_int_t ldt,
    float *x,       magma_tally4_int_t ldx,
    const float *cnorm,
    magma_tally4_int_t *info);
#endif

magma_tally4_int_t
magma_tally4_slatrd(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    float *A, magma_tally4_int_t lda,
    float *e, float *tau,
    float *W, magma_tally4_int_t ldw,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dW, magma_tally4_int_t lddw);

magma_tally4_int_t
magma_tally4_slatrd2(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    float *A,  magma_tally4_int_t lda,
    float *e, float *tau,
    float *W,  magma_tally4_int_t ldw,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dW, magma_tally4_int_t lddw,
    magma_tally4Float_ptr dwork, magma_tally4_int_t ldwork);

#ifdef COMPLEX
magma_tally4_int_t
magma_tally4_slatrsd(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_diag_t diag, magma_tally4_bool_t normin,
    magma_tally4_int_t n, const float *A, magma_tally4_int_t lda,
    float lambda,
    float *x,
    float *scale, float *cnorm,
    magma_tally4_int_t *info);
#endif

magma_tally4_int_t
magma_tally4_slauum(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sposv(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    float *A, magma_tally4_int_t lda,
    float *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_spotrf(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_spotri(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sstedx(
    magma_tally4_range_t range, magma_tally4_int_t n, float vl, float vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, float *d, float *e,
    float *Z, magma_tally4_int_t ldz,
    float *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4Float_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_strevc3(
    magma_tally4_side_t side, magma_tally4_vec_t howmany,
    magma_tally4_int_t *select, magma_tally4_int_t n,
    float *T,  magma_tally4_int_t ldt,
    float *VL, magma_tally4_int_t ldvl,
    float *VR, magma_tally4_int_t ldvr,
    magma_tally4_int_t mm, magma_tally4_int_t *mout,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_strevc3_mt(
    magma_tally4_side_t side, magma_tally4_vec_t howmany,
    magma_tally4_int_t *select, magma_tally4_int_t n,
    float *T,  magma_tally4_int_t ldt,
    float *VL, magma_tally4_int_t ldvl,
    float *VR, magma_tally4_int_t ldvr,
    magma_tally4_int_t mm, magma_tally4_int_t *mout,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_strtri(
    magma_tally4_uplo_t uplo, magma_tally4_diag_t diag, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sorghr(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    float *A, magma_tally4_int_t lda,
    float *tau,
    magma_tally4Float_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sorgqr(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    float *A, magma_tally4_int_t lda,
    float *tau,
    magma_tally4Float_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sorgqr2(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    float *A, magma_tally4_int_t lda,
    float *tau,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sormbr(
    magma_tally4_vect_t vect, magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    float *A, magma_tally4_int_t lda,
    float *tau,
    float *C, magma_tally4_int_t ldc,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sormlq(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    float *A, magma_tally4_int_t lda,
    float *tau,
    float *C, magma_tally4_int_t ldc,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

// not yet implemented
//magma_tally4_int_t magma_tally4_sunmrq( magma_tally4_side_t side, magma_tally4_trans_t trans,
//                          magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
//                          float *A, magma_tally4_int_t lda,
//                          float *tau,
//                          float *C, magma_tally4_int_t ldc,
//                          float *work, magma_tally4_int_t lwork,
//                          magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sormql(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    float *A, magma_tally4_int_t lda,
    float *tau,
    float *C, magma_tally4_int_t ldc,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sormqr(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    float *A, magma_tally4_int_t lda,
    float *tau,
    float *C, magma_tally4_int_t ldc,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sormtr(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n,
    float *A,    magma_tally4_int_t lda,
    float *tau,
    float *C,    magma_tally4_int_t ldc,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 function definitions / Data on CPU / Multi-GPU (alphabetical order)
*/
magma_tally4_int_t
magma_tally4_sgeev_m(
    magma_tally4_vec_t jobvl, magma_tally4_vec_t jobvr, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    #ifdef COMPLEX
    float *w,
    #else
    float *wr, float *wi,
    #endif
    float *VL, magma_tally4_int_t ldvl,
    float *VR, magma_tally4_int_t ldvr,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgehrd_m(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    float *A, magma_tally4_int_t lda,
    float *tau,
    float *work, magma_tally4_int_t lwork,
    float *T,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgetrf_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssyevd_m(
    magma_tally4_int_t ngpu,
    magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float *w,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssyevdx_2stage_m(
    magma_tally4_int_t ngpu,
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssyevdx_m(
    magma_tally4_int_t ngpu,
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssygst_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t itype, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssygvd_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float *B, magma_tally4_int_t ldb,
    float *w,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssygvdx_2stage_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float *B, magma_tally4_int_t ldb,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssygvdx_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float *B, magma_tally4_int_t ldb,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    float *work, magma_tally4_int_t lwork,
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
magma_tally4_slahr2_m(
    magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4_int_t nb,
    float *A, magma_tally4_int_t lda,
    float *tau,
    float *T, magma_tally4_int_t ldt,
    float *Y, magma_tally4_int_t ldy,
    struct sgehrd_data_tally4 *data );

magma_tally4_int_t
magma_tally4_slahru_m(
    magma_tally4_int_t n, magma_tally4_int_t ihi, magma_tally4_int_t k, magma_tally4_int_t nb,
    float *A, magma_tally4_int_t lda,
    struct sgehrd_data_tally4 *data );

magma_tally4_int_t
magma_tally4_spotrf_m(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sstedx_m(
    magma_tally4_int_t ngpu,
    magma_tally4_range_t range, magma_tally4_int_t n, float vl, float vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, float *d, float *e,
    float *Z, magma_tally4_int_t ldz,
    float *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_strsm_m(
    magma_tally4_int_t ngpu,
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t transa, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n, float alpha,
    float *A, magma_tally4_int_t lda,
    float *B, magma_tally4_int_t ldb);

magma_tally4_int_t
magma_tally4_sorghr_m(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    float *A, magma_tally4_int_t lda,
    float *tau,
    float *T, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sorgqr_m(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    float *A, magma_tally4_int_t lda,
    float *tau,
    float *T, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sormqr_m(
    magma_tally4_int_t ngpu,
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    float *A,    magma_tally4_int_t lda,
    float *tau,
    float *C,    magma_tally4_int_t ldc,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sormtr_m(
    magma_tally4_int_t ngpu,
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n,
    float *A,    magma_tally4_int_t lda,
    float *tau,
    float *C,    magma_tally4_int_t ldc,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 function definitions / Data on GPU (alphabetical order)
*/
magma_tally4_int_t
magma_tally4_sgegqr_gpu(
    magma_tally4_int_t ikind, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dwork, float *work,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgelqf_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *tau,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgels_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb,
    float *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgels3_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb,
    float *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqp3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt, float *tau,
    magma_tally4Float_ptr dwork, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqr2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dtau,
    magma_tally4Float_ptr        dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqr2x_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dtau,
    magma_tally4Float_ptr dT, magma_tally4Float_ptr ddA,
    magma_tally4Float_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqr2x2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dtau,
    magma_tally4Float_ptr dT, magma_tally4Float_ptr ddA,
    magma_tally4Float_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqr2x3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dtau,
    magma_tally4Float_ptr dT,
    magma_tally4Float_ptr ddA,
    magma_tally4Float_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqr2x4_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dtau,
    magma_tally4Float_ptr dT, magma_tally4Float_ptr ddA,
    magma_tally4Float_ptr dwork,
    magma_tally4_queue_t queue,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqrf_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *tau,
    magma_tally4Float_ptr dT,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqrf2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *tau,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqrf2_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dlA[], magma_tally4_int_t ldda,
    float *tau,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqrf3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *tau,
    magma_tally4Float_ptr dT,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqrs_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *tau,
    magma_tally4Float_ptr dT,
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb,
    float *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgeqrs3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *tau,
    magma_tally4Float_ptr dT,
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb,
    float *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgerbt_gpu(
    magma_tally4_bool_t gen, magma_tally4_int_t n, magma_tally4_int_t nrhs, 
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda, 
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb, 
    float *U, float *V,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgessm_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4_int_t ib,
    magma_tally4_int_t *ipiv,
    magma_tally4Float_ptr dL1, magma_tally4_int_t lddl1,
    magma_tally4Float_ptr dL,  magma_tally4_int_t lddl,
    magma_tally4Float_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgesv_gpu(
    magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgesv_nopiv_gpu( 
    magma_tally4_int_t n, magma_tally4_int_t nrhs, 
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb, 
                 magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgetf2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgetrf_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgetrf_incpiv_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t ib,
    float    *hA, magma_tally4_int_t ldha,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float    *hL, magma_tally4_int_t ldhl,
    magma_tally4Float_ptr dL, magma_tally4_int_t lddl,
    magma_tally4_int_t *ipiv,
    magma_tally4Float_ptr dwork, magma_tally4_int_t lddwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgetrf_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr d_lA[], magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgetrf_nopiv_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgetrf2_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t offset,
    magma_tally4Float_ptr d_lAT[], magma_tally4_int_t lddat, magma_tally4_int_t *ipiv,
    magma_tally4Float_ptr d_lAP[],
    float *W, magma_tally4_int_t ldw,
    magma_tally4_queue_t queues[][2],
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgetri_gpu(
    magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4Float_ptr dwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgetrs_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sgetrs_nopiv_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssyevd_gpu(
    magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *w,
    float *wA,  magma_tally4_int_t ldwa,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssyevdx_gpu(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float vl, float vu,
    magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, float *w,
    float *wA,  magma_tally4_int_t ldwa,
    float *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally4_int_t
magma_tally4_ssyevr_gpu(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float vl, float vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, float abstol, magma_tally4_int_t *m,
    float *w,
    magma_tally4Float_ptr dZ, magma_tally4_int_t lddz,
    magma_tally4_int_t *isuppz,
    float *wA, magma_tally4_int_t ldwa,
    float *wZ, magma_tally4_int_t ldwz,
    float *work, magma_tally4_int_t lwork,
    float *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssyevx_gpu(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float vl, float vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    float abstol, magma_tally4_int_t *m,
    float *w,
    magma_tally4Float_ptr dZ, magma_tally4_int_t lddz,
    float *wA, magma_tally4_int_t ldwa,
    float *wZ, magma_tally4_int_t ldwz,
    float *work, magma_tally4_int_t lwork,
    float *rwork, magma_tally4_int_t *iwork,
    magma_tally4_int_t *ifail,
    magma_tally4_int_t *info);
#endif  // COMPLEX

magma_tally4_int_t
magma_tally4_ssygst_gpu(
    magma_tally4_int_t itype, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssytrd_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *d, float *e, float *tau,
    float *wA,  magma_tally4_int_t ldwa,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssytrd_sy2sb_mgpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    float *A, magma_tally4_int_t lda,
    float *tau,
    float *work, magma_tally4_int_t lwork,
    magma_tally4Float_ptr dAmgpu[], magma_tally4_int_t ldda,
    magma_tally4Float_ptr dTmgpu[], magma_tally4_int_t lddt,
    magma_tally4_int_t ngpu, magma_tally4_int_t distblk,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssytrd_sy2sb_mgpu_spec(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    float *A, magma_tally4_int_t lda,
    float *tau,
    float *work, magma_tally4_int_t lwork,
    magma_tally4Float_ptr dAmgpu[], magma_tally4_int_t ldda,
    magma_tally4Float_ptr dTmgpu[], magma_tally4_int_t lddt,
    magma_tally4_int_t ngpu, magma_tally4_int_t distblk,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssytrd_mgpu(
    magma_tally4_int_t ngpu, magma_tally4_int_t nqueue,
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    float *A, magma_tally4_int_t lda,
    float *d, float *e, float *tau,
    float *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssytrd2_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *d, float *e, float *tau,
    float *wA,  magma_tally4_int_t ldwa,
    float *work, magma_tally4_int_t lwork,
    magma_tally4Float_ptr dwork, magma_tally4_int_t ldwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ssytrf_nopiv_tally4_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_slabrd_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nb,
    float     *A, magma_tally4_int_t lda,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *d, float *e, float *tauq, float *taup,
    float     *X, magma_tally4_int_t ldx,
    magma_tally4Float_ptr dX, magma_tally4_int_t lddx,
    float     *Y, magma_tally4_int_t ldy,
    magma_tally4Float_ptr dY, magma_tally4_int_t lddy );

magma_tally4_int_t
magma_tally4_slaqps_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t offset,
    magma_tally4_int_t nb, magma_tally4_int_t *kb,
    magma_tally4Float_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt, float *tau,
    float *vn1, float *vn2,
    magma_tally4Float_ptr dauxv,
    magma_tally4Float_ptr dF, magma_tally4_int_t lddf);

magma_tally4_int_t
magma_tally4_slaqps2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t offset,
    magma_tally4_int_t nb, magma_tally4_int_t *kb,
    magma_tally4Float_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt,
    magma_tally4Float_ptr dtau,
    magma_tally4Float_ptr dvn1, magma_tally4Float_ptr dvn2,
    magma_tally4Float_ptr dauxv,
    magma_tally4Float_ptr dF, magma_tally4_int_t lddf);

magma_tally4_int_t
magma_tally4_slaqps3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t offset,
    magma_tally4_int_t nb, magma_tally4_int_t *kb,
    magma_tally4Float_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt,
    magma_tally4Float_ptr dtau,
    magma_tally4Float_ptr dvn1, magma_tally4Float_ptr dvn2,
    magma_tally4Float_ptr dauxv,
    magma_tally4Float_ptr dF, magma_tally4_int_t lddf);

magma_tally4_int_t
magma_tally4_slarf_gpu(
    magma_tally4_int_t m,  magma_tally4_int_t n,
    magma_tally4Float_const_ptr dv, magma_tally4Float_const_ptr dtau,
    magma_tally4Float_ptr dC,  magma_tally4_int_t lddc);

magma_tally4_int_t
magma_tally4_slarfb_gpu(
    magma_tally4_side_t side, magma_tally4_trans_t trans, magma_tally4_direct_t direct, magma_tally4_storev_t storev,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Float_const_ptr dV, magma_tally4_int_t lddv,
    magma_tally4Float_const_ptr dT, magma_tally4_int_t lddt,
    magma_tally4Float_ptr dC,       magma_tally4_int_t lddc,
    magma_tally4Float_ptr dwork,    magma_tally4_int_t ldwork );

magma_tally4_int_t
magma_tally4_slarfb_gpu_gemm(
    magma_tally4_side_t side, magma_tally4_trans_t trans, magma_tally4_direct_t direct, magma_tally4_storev_t storev,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Float_const_ptr dV, magma_tally4_int_t lddv,
    magma_tally4Float_const_ptr dT, magma_tally4_int_t lddt,
    magma_tally4Float_ptr dC,       magma_tally4_int_t lddc,
    magma_tally4Float_ptr dwork,    magma_tally4_int_t ldwork,
    magma_tally4Float_ptr dworkvt,  magma_tally4_int_t ldworkvt);

magma_tally4_int_t
magma_tally4_slarfb2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Float_const_ptr dV, magma_tally4_int_t lddv,
    magma_tally4Float_const_ptr dT, magma_tally4_int_t lddt,
    magma_tally4Float_ptr dC,       magma_tally4_int_t lddc,
    magma_tally4Float_ptr dwork,    magma_tally4_int_t ldwork );

magma_tally4_int_t
magma_tally4_slatrd_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t nb0,
    float *A,  magma_tally4_int_t lda,
    float *e, float *tau,
    float    *W,       magma_tally4_int_t ldw,
    magma_tally4Float_ptr dA[],    magma_tally4_int_t ldda, magma_tally4_int_t offset,
    magma_tally4Float_ptr dW[],    magma_tally4_int_t lddw,
    float    *hwork,   magma_tally4_int_t lhwork,
    magma_tally4Float_ptr dwork[], magma_tally4_int_t ldwork,
    magma_tally4_queue_t queues[] );

magma_tally4_int_t
magma_tally4_slauum_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sposv_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_spotf2_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_spotrf_gpu(
    magma_tally4_uplo_t uplo,  magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_spotrf_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Float_ptr d_lA[], magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_spotrf_mgpu_right(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Float_ptr d_lA[], magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_spotrf3_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t off_i, magma_tally4_int_t off_j, magma_tally4_int_t nb,
    magma_tally4Float_ptr d_lA[],  magma_tally4_int_t ldda,
    magma_tally4Float_ptr d_lP[],  magma_tally4_int_t lddp,
    float *A, magma_tally4_int_t lda, magma_tally4_int_t h,
    magma_tally4_queue_t queues[][3], magma_tally4_event_t events[][5],
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_spotri_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_spotrs_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sssssm_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m1, magma_tally4_int_t n1,
    magma_tally4_int_t m2, magma_tally4_int_t n2, magma_tally4_int_t k, magma_tally4_int_t ib,
    magma_tally4Float_ptr dA1, magma_tally4_int_t ldda1,
    magma_tally4Float_ptr dA2, magma_tally4_int_t ldda2,
    magma_tally4Float_ptr dL1, magma_tally4_int_t lddl1,
    magma_tally4Float_ptr dL2, magma_tally4_int_t lddl2,
    magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_strtri_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_diag_t diag, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_stsqrt_gpu(
    magma_tally4_int_t *m, magma_tally4_int_t *n,
    float *A1, float *A2, magma_tally4_int_t *lda,
    float *tau,
    float *work, magma_tally4_int_t *lwork,
    magma_tally4Float_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_ststrf_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t ib, magma_tally4_int_t nb,
    float    *hU, magma_tally4_int_t ldhu,
    magma_tally4Float_ptr dU, magma_tally4_int_t lddu,
    float    *hA, magma_tally4_int_t ldha,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float    *hL, magma_tally4_int_t ldhl,
    magma_tally4Float_ptr dL, magma_tally4_int_t lddl,
    magma_tally4_int_t *ipiv,
    float *hwork, magma_tally4_int_t ldhwork,
    magma_tally4Float_ptr dwork, magma_tally4_int_t lddwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sorgqr_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *tau,
    magma_tally4Float_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sormql2_gpu(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *tau,
    magma_tally4Float_ptr dC, magma_tally4_int_t lddc,
    float *wA, magma_tally4_int_t ldwa,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sormqr_gpu(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *tau,
    magma_tally4Float_ptr dC, magma_tally4_int_t lddc,
    float *hwork, magma_tally4_int_t lwork,
    magma_tally4Float_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sormqr2_gpu(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *tau,
    magma_tally4Float_ptr dC, magma_tally4_int_t lddc,
    float    *wA, magma_tally4_int_t ldwa,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_sormtr_gpu(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *tau,
    magma_tally4Float_ptr dC, magma_tally4_int_t lddc,
    float    *wA, magma_tally4_int_t ldwa,
    magma_tally4_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 utility function definitions
*/

extern const float MAGMA_tally4_S_NAN;
extern const float MAGMA_tally4_S_INF;

magma_tally4_int_t
magma_tally4_snan_inf(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    const float *A, magma_tally4_int_t lda,
    magma_tally4_int_t *cnt_nan,
    magma_tally4_int_t *cnt_inf );

magma_tally4_int_t
magma_tally4_snan_inf_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *cnt_nan,
    magma_tally4_int_t *cnt_inf );

void magma_tally4_sprint(
    magma_tally4_int_t m, magma_tally4_int_t n,
    const float *A, magma_tally4_int_t lda );

void magma_tally4_sprint_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr dA, magma_tally4_int_t ldda );

void spanel_to_q_tally4(
    magma_tally4_uplo_t uplo, magma_tally4_int_t ib,
    float *A, magma_tally4_int_t lda,
    float *work );

void sq_to_panel_tally4(
    magma_tally4_uplo_t uplo, magma_tally4_int_t ib,
    float *A, magma_tally4_int_t lda,
    float *work );

#ifdef __cplusplus
}
#endif

#undef REAL

#endif /* MAGMA_tally4_S_H */
