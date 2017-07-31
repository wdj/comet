/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally2_z.h normal z -> s, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally2_S_H
#define MAGMA_tally2_S_H

#include "magma_tally2_types.h"
#include "magma_tally2_sgehrd_m.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 Auxiliary functions to get the NB used
*/
#ifdef REAL
magma_tally2_int_t magma_tally2_get_slaex3_m_nb();       // defined in slaex3_m.cpp
#endif

magma_tally2_int_t magma_tally2_get_spotrf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_sgetrf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_sgetri_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_sgeqp3_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_sgeqrf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_sgeqlf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_sgehrd_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_ssytrd_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_ssytrf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_ssytrf_nopiv_tally2_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_sgelqf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_sgebrd_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_ssygst_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_sgesvd_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_ssygst_nb_m( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_sbulge_nb( magma_tally2_int_t m, magma_tally2_int_t nbthreads );
magma_tally2_int_t magma_tally2_get_sbulge_nb_mgpu( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_sbulge_get_Vblksiz( magma_tally2_int_t m, magma_tally2_int_t nb, magma_tally2_int_t nbthreads );
magma_tally2_int_t magma_tally2_get_sbulge_gcperf();


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
magma_tally2_svrange(
    magma_tally2_int_t k, float *d, magma_tally2_int_t *il, magma_tally2_int_t *iu, float vl, float vu);

void
magma_tally2_sirange(
    magma_tally2_int_t k, magma_tally2_int_t *indxq, magma_tally2_int_t *iil, magma_tally2_int_t *iiu, magma_tally2_int_t il, magma_tally2_int_t iu);
#endif

magma_tally2_int_t
magma_tally2_sgebrd(
    magma_tally2_int_t m, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float *d, float *e,
    float *tauq, float *taup,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeev(
    magma_tally2_vec_t jobvl, magma_tally2_vec_t jobvr, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    #ifdef COMPLEX
    float *w,
    #else
    float *wr, float *wi,
    #endif
    float *VL, magma_tally2_int_t ldvl,
    float *VR, magma_tally2_int_t ldvr,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgehrd(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    float *A, magma_tally2_int_t lda,
    float *tau,
    float *work, magma_tally2_int_t lwork,
    magma_tally2Float_ptr dT,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgehrd2(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    float *A, magma_tally2_int_t lda,
    float *tau,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgelqf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    float *A,    magma_tally2_int_t lda,
    float *tau,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqlf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    float *A,    magma_tally2_int_t lda,
    float *tau,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqp3(
    magma_tally2_int_t m, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    magma_tally2_int_t *jpvt, float *tau,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqrf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float *tau,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqrf_ooc(
    magma_tally2_int_t m, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float *tau,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqrf4(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    float *A,    magma_tally2_int_t lda,
    float *tau,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgesdd(
    magma_tally2_vec_t jobz, magma_tally2_int_t m, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float *s,
    float *U, magma_tally2_int_t ldu,
    float *VT, magma_tally2_int_t ldvt,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally2_int_t *iwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgesv(
    magma_tally2_int_t n, magma_tally2_int_t nrhs,
    float *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    float *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgesv_rbt(
    magma_tally2_bool_t ref, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    float *A, magma_tally2_int_t lda, 
    float *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgesvd(
    magma_tally2_vec_t jobu, magma_tally2_vec_t jobvt, magma_tally2_int_t m, magma_tally2_int_t n,
    float *A,    magma_tally2_int_t lda, float *s,
    float *U,    magma_tally2_int_t ldu,
    float *VT,   magma_tally2_int_t ldvt,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgetf2_nopiv(
    magma_tally2_int_t m, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgetrf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgetrf_nopiv(
    magma_tally2_int_t m, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgetrf_piv(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t NB,
    float *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgetrf2(
    magma_tally2_int_t m, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssyevd(
    magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float *w,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssyevdx(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssyevdx_2stage(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally2_int_t
magma_tally2_ssyevr(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float vl, float vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, float abstol, magma_tally2_int_t *m,
    float *w,
    float *Z, magma_tally2_int_t ldz,
    magma_tally2_int_t *isuppz,
    float *work, magma_tally2_int_t lwork,
    float *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

// no real [sd] precisions available
magma_tally2_int_t
magma_tally2_ssyevx(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float vl, float vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, float abstol, magma_tally2_int_t *m,
    float *w,
    float *Z, magma_tally2_int_t ldz,
    float *work, magma_tally2_int_t lwork,
    float *rwork, magma_tally2_int_t *iwork,
    magma_tally2_int_t *ifail,
    magma_tally2_int_t *info);
#endif  // COMPLEX

magma_tally2_int_t
magma_tally2_ssygst(
    magma_tally2_int_t itype, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssygvd(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float *B, magma_tally2_int_t ldb,
    float *w, float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssygvdx(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n, float *A, magma_tally2_int_t lda,
    float *B, magma_tally2_int_t ldb,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssygvdx_2stage(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float *B, magma_tally2_int_t ldb,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssygvr(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float *B, magma_tally2_int_t ldb,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    float abstol, magma_tally2_int_t *m, float *w,
    float *Z, magma_tally2_int_t ldz,
    magma_tally2_int_t *isuppz, float *work, magma_tally2_int_t lwork,
    float *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssygvx(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n, float *A, magma_tally2_int_t lda,
    float *B, magma_tally2_int_t ldb,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    float abstol, magma_tally2_int_t *m, float *w,
    float *Z, magma_tally2_int_t ldz,
    float *work, magma_tally2_int_t lwork, float *rwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t *ifail,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssysv(magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
            float *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
            float *B, magma_tally2_int_t ldb,
            magma_tally2_int_t *info );

magma_tally2_int_t
magma_tally2_ssytrd(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float *d, float *e, float *tau,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssytrf(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssytrf_nopiv_tally2(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssytrd_sb2st(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t Vblksiz,
    float *A, magma_tally2_int_t lda,
    float *d, float *e,
    float *V, magma_tally2_int_t ldv,
    float *TAU, magma_tally2_int_t compT,
    float *T, magma_tally2_int_t ldt);

magma_tally2_int_t
magma_tally2_ssytrd_sy2sb(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    float *A, magma_tally2_int_t lda,
    float *tau,
    float *work, magma_tally2_int_t lwork,
    magma_tally2Float_ptr dT,
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
magma_tally2_slasyf_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t *kb,
    float    *hA, magma_tally2_int_t lda,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2Float_ptr dW, magma_tally2_int_t lddw,
    magma_tally2_queue_t queues[], magma_tally2_event_t event[],
    magma_tally2_int_t *info);

magma_tally2_int_t
ssytrf_nopiv_tally2_cpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t ib,
    float *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssytrs_nopiv_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssysv_nopiv_gpu(
    magma_tally2_uplo_t uplo,  magma_tally2_int_t n, magma_tally2_int_t nrhs, 
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dB, magma_tally2_int_t lddb, 
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_slahr2(
    magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2_int_t nb,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dV, magma_tally2_int_t lddv,
    float *A,  magma_tally2_int_t lda,
    float *tau,
    float *T,  magma_tally2_int_t ldt,
    float *Y,  magma_tally2_int_t ldy);

magma_tally2_int_t
magma_tally2_slahru(
    magma_tally2_int_t n, magma_tally2_int_t ihi, magma_tally2_int_t k, magma_tally2_int_t nb,
    float     *A, magma_tally2_int_t lda,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dY, magma_tally2_int_t lddy,
    magma_tally2Float_ptr dV, magma_tally2_int_t lddv,
    magma_tally2Float_ptr dT,
    magma_tally2Float_ptr dwork);

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
magma_tally2_slaqps(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t offset,
    magma_tally2_int_t nb, magma_tally2_int_t *kb,
    float *A,  magma_tally2_int_t lda,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt, float *tau, float *vn1, float *vn2,
    float *auxv,
    float *F,  magma_tally2_int_t ldf,
    magma_tally2Float_ptr dF, magma_tally2_int_t lddf );

#ifdef REAL
magma_tally2_int_t
magma_tally2_slaqtrsd(
    magma_tally2_trans_t trans, magma_tally2_int_t n,
    const float *T, magma_tally2_int_t ldt,
    float *x,       magma_tally2_int_t ldx,
    const float *cnorm,
    magma_tally2_int_t *info);
#endif

magma_tally2_int_t
magma_tally2_slatrd(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    float *A, magma_tally2_int_t lda,
    float *e, float *tau,
    float *W, magma_tally2_int_t ldw,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dW, magma_tally2_int_t lddw);

magma_tally2_int_t
magma_tally2_slatrd2(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    float *A,  magma_tally2_int_t lda,
    float *e, float *tau,
    float *W,  magma_tally2_int_t ldw,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dW, magma_tally2_int_t lddw,
    magma_tally2Float_ptr dwork, magma_tally2_int_t ldwork);

#ifdef COMPLEX
magma_tally2_int_t
magma_tally2_slatrsd(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_diag_t diag, magma_tally2_bool_t normin,
    magma_tally2_int_t n, const float *A, magma_tally2_int_t lda,
    float lambda,
    float *x,
    float *scale, float *cnorm,
    magma_tally2_int_t *info);
#endif

magma_tally2_int_t
magma_tally2_slauum(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sposv(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    float *A, magma_tally2_int_t lda,
    float *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_spotrf(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_spotri(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sstedx(
    magma_tally2_range_t range, magma_tally2_int_t n, float vl, float vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, float *d, float *e,
    float *Z, magma_tally2_int_t ldz,
    float *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2Float_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_strevc3(
    magma_tally2_side_t side, magma_tally2_vec_t howmany,
    magma_tally2_int_t *select, magma_tally2_int_t n,
    float *T,  magma_tally2_int_t ldt,
    float *VL, magma_tally2_int_t ldvl,
    float *VR, magma_tally2_int_t ldvr,
    magma_tally2_int_t mm, magma_tally2_int_t *mout,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_strevc3_mt(
    magma_tally2_side_t side, magma_tally2_vec_t howmany,
    magma_tally2_int_t *select, magma_tally2_int_t n,
    float *T,  magma_tally2_int_t ldt,
    float *VL, magma_tally2_int_t ldvl,
    float *VR, magma_tally2_int_t ldvr,
    magma_tally2_int_t mm, magma_tally2_int_t *mout,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_strtri(
    magma_tally2_uplo_t uplo, magma_tally2_diag_t diag, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sorghr(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    float *A, magma_tally2_int_t lda,
    float *tau,
    magma_tally2Float_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sorgqr(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    float *A, magma_tally2_int_t lda,
    float *tau,
    magma_tally2Float_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sorgqr2(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    float *A, magma_tally2_int_t lda,
    float *tau,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sormbr(
    magma_tally2_vect_t vect, magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    float *A, magma_tally2_int_t lda,
    float *tau,
    float *C, magma_tally2_int_t ldc,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sormlq(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    float *A, magma_tally2_int_t lda,
    float *tau,
    float *C, magma_tally2_int_t ldc,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

// not yet implemented
//magma_tally2_int_t magma_tally2_sunmrq( magma_tally2_side_t side, magma_tally2_trans_t trans,
//                          magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
//                          float *A, magma_tally2_int_t lda,
//                          float *tau,
//                          float *C, magma_tally2_int_t ldc,
//                          float *work, magma_tally2_int_t lwork,
//                          magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sormql(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    float *A, magma_tally2_int_t lda,
    float *tau,
    float *C, magma_tally2_int_t ldc,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sormqr(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    float *A, magma_tally2_int_t lda,
    float *tau,
    float *C, magma_tally2_int_t ldc,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sormtr(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n,
    float *A,    magma_tally2_int_t lda,
    float *tau,
    float *C,    magma_tally2_int_t ldc,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 function definitions / Data on CPU / Multi-GPU (alphabetical order)
*/
magma_tally2_int_t
magma_tally2_sgeev_m(
    magma_tally2_vec_t jobvl, magma_tally2_vec_t jobvr, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    #ifdef COMPLEX
    float *w,
    #else
    float *wr, float *wi,
    #endif
    float *VL, magma_tally2_int_t ldvl,
    float *VR, magma_tally2_int_t ldvr,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgehrd_m(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    float *A, magma_tally2_int_t lda,
    float *tau,
    float *work, magma_tally2_int_t lwork,
    float *T,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgetrf_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssyevd_m(
    magma_tally2_int_t ngpu,
    magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float *w,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssyevdx_2stage_m(
    magma_tally2_int_t ngpu,
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssyevdx_m(
    magma_tally2_int_t ngpu,
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssygst_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssygvd_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float *B, magma_tally2_int_t ldb,
    float *w,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssygvdx_2stage_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float *B, magma_tally2_int_t ldb,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssygvdx_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float *B, magma_tally2_int_t ldb,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    float *work, magma_tally2_int_t lwork,
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
magma_tally2_slahr2_m(
    magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2_int_t nb,
    float *A, magma_tally2_int_t lda,
    float *tau,
    float *T, magma_tally2_int_t ldt,
    float *Y, magma_tally2_int_t ldy,
    struct sgehrd_data_tally2 *data );

magma_tally2_int_t
magma_tally2_slahru_m(
    magma_tally2_int_t n, magma_tally2_int_t ihi, magma_tally2_int_t k, magma_tally2_int_t nb,
    float *A, magma_tally2_int_t lda,
    struct sgehrd_data_tally2 *data );

magma_tally2_int_t
magma_tally2_spotrf_m(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sstedx_m(
    magma_tally2_int_t ngpu,
    magma_tally2_range_t range, magma_tally2_int_t n, float vl, float vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, float *d, float *e,
    float *Z, magma_tally2_int_t ldz,
    float *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_strsm_m(
    magma_tally2_int_t ngpu,
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t transa, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n, float alpha,
    float *A, magma_tally2_int_t lda,
    float *B, magma_tally2_int_t ldb);

magma_tally2_int_t
magma_tally2_sorghr_m(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    float *A, magma_tally2_int_t lda,
    float *tau,
    float *T, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sorgqr_m(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    float *A, magma_tally2_int_t lda,
    float *tau,
    float *T, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sormqr_m(
    magma_tally2_int_t ngpu,
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    float *A,    magma_tally2_int_t lda,
    float *tau,
    float *C,    magma_tally2_int_t ldc,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sormtr_m(
    magma_tally2_int_t ngpu,
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n,
    float *A,    magma_tally2_int_t lda,
    float *tau,
    float *C,    magma_tally2_int_t ldc,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 function definitions / Data on GPU (alphabetical order)
*/
magma_tally2_int_t
magma_tally2_sgegqr_gpu(
    magma_tally2_int_t ikind, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dwork, float *work,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgelqf_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *tau,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgels_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dB, magma_tally2_int_t lddb,
    float *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgels3_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dB, magma_tally2_int_t lddb,
    float *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqp3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt, float *tau,
    magma_tally2Float_ptr dwork, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqr2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dtau,
    magma_tally2Float_ptr        dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqr2x_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dtau,
    magma_tally2Float_ptr dT, magma_tally2Float_ptr ddA,
    magma_tally2Float_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqr2x2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dtau,
    magma_tally2Float_ptr dT, magma_tally2Float_ptr ddA,
    magma_tally2Float_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqr2x3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dtau,
    magma_tally2Float_ptr dT,
    magma_tally2Float_ptr ddA,
    magma_tally2Float_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqr2x4_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dtau,
    magma_tally2Float_ptr dT, magma_tally2Float_ptr ddA,
    magma_tally2Float_ptr dwork,
    magma_tally2_queue_t queue,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqrf_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *tau,
    magma_tally2Float_ptr dT,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqrf2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *tau,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqrf2_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dlA[], magma_tally2_int_t ldda,
    float *tau,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqrf3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *tau,
    magma_tally2Float_ptr dT,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqrs_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *tau,
    magma_tally2Float_ptr dT,
    magma_tally2Float_ptr dB, magma_tally2_int_t lddb,
    float *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgeqrs3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *tau,
    magma_tally2Float_ptr dT,
    magma_tally2Float_ptr dB, magma_tally2_int_t lddb,
    float *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgerbt_gpu(
    magma_tally2_bool_t gen, magma_tally2_int_t n, magma_tally2_int_t nrhs, 
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda, 
    magma_tally2Float_ptr dB, magma_tally2_int_t lddb, 
    float *U, float *V,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgessm_gpu(
    magma_tally2_order_t order, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2_int_t ib,
    magma_tally2_int_t *ipiv,
    magma_tally2Float_ptr dL1, magma_tally2_int_t lddl1,
    magma_tally2Float_ptr dL,  magma_tally2_int_t lddl,
    magma_tally2Float_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgesv_gpu(
    magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2Float_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgesv_nopiv_gpu( 
    magma_tally2_int_t n, magma_tally2_int_t nrhs, 
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dB, magma_tally2_int_t lddb, 
                 magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgetf2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgetrf_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgetrf_incpiv_gpu(
    magma_tally2_order_t order, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t ib,
    float    *hA, magma_tally2_int_t ldha,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float    *hL, magma_tally2_int_t ldhl,
    magma_tally2Float_ptr dL, magma_tally2_int_t lddl,
    magma_tally2_int_t *ipiv,
    magma_tally2Float_ptr dwork, magma_tally2_int_t lddwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgetrf_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr d_lA[], magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgetrf_nopiv_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgetrf2_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t offset,
    magma_tally2Float_ptr d_lAT[], magma_tally2_int_t lddat, magma_tally2_int_t *ipiv,
    magma_tally2Float_ptr d_lAP[],
    float *W, magma_tally2_int_t ldw,
    magma_tally2_queue_t queues[][2],
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgetri_gpu(
    magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2Float_ptr dwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgetrs_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2Float_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sgetrs_nopiv_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssyevd_gpu(
    magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *w,
    float *wA,  magma_tally2_int_t ldwa,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssyevdx_gpu(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float vl, float vu,
    magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, float *w,
    float *wA,  magma_tally2_int_t ldwa,
    float *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally2_int_t
magma_tally2_ssyevr_gpu(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float vl, float vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, float abstol, magma_tally2_int_t *m,
    float *w,
    magma_tally2Float_ptr dZ, magma_tally2_int_t lddz,
    magma_tally2_int_t *isuppz,
    float *wA, magma_tally2_int_t ldwa,
    float *wZ, magma_tally2_int_t ldwz,
    float *work, magma_tally2_int_t lwork,
    float *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssyevx_gpu(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float vl, float vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    float abstol, magma_tally2_int_t *m,
    float *w,
    magma_tally2Float_ptr dZ, magma_tally2_int_t lddz,
    float *wA, magma_tally2_int_t ldwa,
    float *wZ, magma_tally2_int_t ldwz,
    float *work, magma_tally2_int_t lwork,
    float *rwork, magma_tally2_int_t *iwork,
    magma_tally2_int_t *ifail,
    magma_tally2_int_t *info);
#endif  // COMPLEX

magma_tally2_int_t
magma_tally2_ssygst_gpu(
    magma_tally2_int_t itype, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssytrd_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *d, float *e, float *tau,
    float *wA,  magma_tally2_int_t ldwa,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssytrd_sy2sb_mgpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    float *A, magma_tally2_int_t lda,
    float *tau,
    float *work, magma_tally2_int_t lwork,
    magma_tally2Float_ptr dAmgpu[], magma_tally2_int_t ldda,
    magma_tally2Float_ptr dTmgpu[], magma_tally2_int_t lddt,
    magma_tally2_int_t ngpu, magma_tally2_int_t distblk,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssytrd_sy2sb_mgpu_spec(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    float *A, magma_tally2_int_t lda,
    float *tau,
    float *work, magma_tally2_int_t lwork,
    magma_tally2Float_ptr dAmgpu[], magma_tally2_int_t ldda,
    magma_tally2Float_ptr dTmgpu[], magma_tally2_int_t lddt,
    magma_tally2_int_t ngpu, magma_tally2_int_t distblk,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssytrd_mgpu(
    magma_tally2_int_t ngpu, magma_tally2_int_t nqueue,
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    float *A, magma_tally2_int_t lda,
    float *d, float *e, float *tau,
    float *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssytrd2_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *d, float *e, float *tau,
    float *wA,  magma_tally2_int_t ldwa,
    float *work, magma_tally2_int_t lwork,
    magma_tally2Float_ptr dwork, magma_tally2_int_t ldwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ssytrf_nopiv_tally2_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_slabrd_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb,
    float     *A, magma_tally2_int_t lda,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *d, float *e, float *tauq, float *taup,
    float     *X, magma_tally2_int_t ldx,
    magma_tally2Float_ptr dX, magma_tally2_int_t lddx,
    float     *Y, magma_tally2_int_t ldy,
    magma_tally2Float_ptr dY, magma_tally2_int_t lddy );

magma_tally2_int_t
magma_tally2_slaqps_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t offset,
    magma_tally2_int_t nb, magma_tally2_int_t *kb,
    magma_tally2Float_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt, float *tau,
    float *vn1, float *vn2,
    magma_tally2Float_ptr dauxv,
    magma_tally2Float_ptr dF, magma_tally2_int_t lddf);

magma_tally2_int_t
magma_tally2_slaqps2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t offset,
    magma_tally2_int_t nb, magma_tally2_int_t *kb,
    magma_tally2Float_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt,
    magma_tally2Float_ptr dtau,
    magma_tally2Float_ptr dvn1, magma_tally2Float_ptr dvn2,
    magma_tally2Float_ptr dauxv,
    magma_tally2Float_ptr dF, magma_tally2_int_t lddf);

magma_tally2_int_t
magma_tally2_slaqps3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t offset,
    magma_tally2_int_t nb, magma_tally2_int_t *kb,
    magma_tally2Float_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt,
    magma_tally2Float_ptr dtau,
    magma_tally2Float_ptr dvn1, magma_tally2Float_ptr dvn2,
    magma_tally2Float_ptr dauxv,
    magma_tally2Float_ptr dF, magma_tally2_int_t lddf);

magma_tally2_int_t
magma_tally2_slarf_gpu(
    magma_tally2_int_t m,  magma_tally2_int_t n,
    magma_tally2Float_const_ptr dv, magma_tally2Float_const_ptr dtau,
    magma_tally2Float_ptr dC,  magma_tally2_int_t lddc);

magma_tally2_int_t
magma_tally2_slarfb_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans, magma_tally2_direct_t direct, magma_tally2_storev_t storev,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2Float_const_ptr dV, magma_tally2_int_t lddv,
    magma_tally2Float_const_ptr dT, magma_tally2_int_t lddt,
    magma_tally2Float_ptr dC,       magma_tally2_int_t lddc,
    magma_tally2Float_ptr dwork,    magma_tally2_int_t ldwork );

magma_tally2_int_t
magma_tally2_slarfb_gpu_gemm(
    magma_tally2_side_t side, magma_tally2_trans_t trans, magma_tally2_direct_t direct, magma_tally2_storev_t storev,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2Float_const_ptr dV, magma_tally2_int_t lddv,
    magma_tally2Float_const_ptr dT, magma_tally2_int_t lddt,
    magma_tally2Float_ptr dC,       magma_tally2_int_t lddc,
    magma_tally2Float_ptr dwork,    magma_tally2_int_t ldwork,
    magma_tally2Float_ptr dworkvt,  magma_tally2_int_t ldworkvt);

magma_tally2_int_t
magma_tally2_slarfb2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2Float_const_ptr dV, magma_tally2_int_t lddv,
    magma_tally2Float_const_ptr dT, magma_tally2_int_t lddt,
    magma_tally2Float_ptr dC,       magma_tally2_int_t lddc,
    magma_tally2Float_ptr dwork,    magma_tally2_int_t ldwork );

magma_tally2_int_t
magma_tally2_slatrd_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t nb0,
    float *A,  magma_tally2_int_t lda,
    float *e, float *tau,
    float    *W,       magma_tally2_int_t ldw,
    magma_tally2Float_ptr dA[],    magma_tally2_int_t ldda, magma_tally2_int_t offset,
    magma_tally2Float_ptr dW[],    magma_tally2_int_t lddw,
    float    *hwork,   magma_tally2_int_t lhwork,
    magma_tally2Float_ptr dwork[], magma_tally2_int_t ldwork,
    magma_tally2_queue_t queues[] );

magma_tally2_int_t
magma_tally2_slauum_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sposv_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_spotf2_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_spotrf_gpu(
    magma_tally2_uplo_t uplo,  magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_spotrf_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Float_ptr d_lA[], magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_spotrf_mgpu_right(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Float_ptr d_lA[], magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_spotrf3_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t off_i, magma_tally2_int_t off_j, magma_tally2_int_t nb,
    magma_tally2Float_ptr d_lA[],  magma_tally2_int_t ldda,
    magma_tally2Float_ptr d_lP[],  magma_tally2_int_t lddp,
    float *A, magma_tally2_int_t lda, magma_tally2_int_t h,
    magma_tally2_queue_t queues[][3], magma_tally2_event_t events[][5],
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_spotri_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_spotrs_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sssssm_gpu(
    magma_tally2_order_t order, magma_tally2_int_t m1, magma_tally2_int_t n1,
    magma_tally2_int_t m2, magma_tally2_int_t n2, magma_tally2_int_t k, magma_tally2_int_t ib,
    magma_tally2Float_ptr dA1, magma_tally2_int_t ldda1,
    magma_tally2Float_ptr dA2, magma_tally2_int_t ldda2,
    magma_tally2Float_ptr dL1, magma_tally2_int_t lddl1,
    magma_tally2Float_ptr dL2, magma_tally2_int_t lddl2,
    magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_strtri_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_diag_t diag, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_stsqrt_gpu(
    magma_tally2_int_t *m, magma_tally2_int_t *n,
    float *A1, float *A2, magma_tally2_int_t *lda,
    float *tau,
    float *work, magma_tally2_int_t *lwork,
    magma_tally2Float_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_ststrf_gpu(
    magma_tally2_order_t order, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t ib, magma_tally2_int_t nb,
    float    *hU, magma_tally2_int_t ldhu,
    magma_tally2Float_ptr dU, magma_tally2_int_t lddu,
    float    *hA, magma_tally2_int_t ldha,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float    *hL, magma_tally2_int_t ldhl,
    magma_tally2Float_ptr dL, magma_tally2_int_t lddl,
    magma_tally2_int_t *ipiv,
    float *hwork, magma_tally2_int_t ldhwork,
    magma_tally2Float_ptr dwork, magma_tally2_int_t lddwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sorgqr_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *tau,
    magma_tally2Float_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sormql2_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *tau,
    magma_tally2Float_ptr dC, magma_tally2_int_t lddc,
    float *wA, magma_tally2_int_t ldwa,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sormqr_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *tau,
    magma_tally2Float_ptr dC, magma_tally2_int_t lddc,
    float *hwork, magma_tally2_int_t lwork,
    magma_tally2Float_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sormqr2_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *tau,
    magma_tally2Float_ptr dC, magma_tally2_int_t lddc,
    float    *wA, magma_tally2_int_t ldwa,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_sormtr_gpu(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *tau,
    magma_tally2Float_ptr dC, magma_tally2_int_t lddc,
    float    *wA, magma_tally2_int_t ldwa,
    magma_tally2_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 utility function definitions
*/

extern const float MAGMA_tally2_S_NAN;
extern const float MAGMA_tally2_S_INF;

magma_tally2_int_t
magma_tally2_snan_inf(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    const float *A, magma_tally2_int_t lda,
    magma_tally2_int_t *cnt_nan,
    magma_tally2_int_t *cnt_inf );

magma_tally2_int_t
magma_tally2_snan_inf_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *cnt_nan,
    magma_tally2_int_t *cnt_inf );

void magma_tally2_sprint(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const float *A, magma_tally2_int_t lda );

void magma_tally2_sprint_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_const_ptr dA, magma_tally2_int_t ldda );

void spanel_to_q_tally2(
    magma_tally2_uplo_t uplo, magma_tally2_int_t ib,
    float *A, magma_tally2_int_t lda,
    float *work );

void sq_to_panel_tally2(
    magma_tally2_uplo_t uplo, magma_tally2_int_t ib,
    float *A, magma_tally2_int_t lda,
    float *work );

#ifdef __cplusplus
}
#endif

#undef REAL

#endif /* MAGMA_tally2_S_H */
