/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally3_z.h normal z -> s, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally3_S_H
#define MAGMA_tally3_S_H

#include "magma_tally3_types.h"
#include "magma_tally3_sgehrd_m.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 Auxiliary functions to get the NB used
*/
#ifdef REAL
magma_tally3_int_t magma_tally3_get_slaex3_m_nb();       // defined in slaex3_m.cpp
#endif

magma_tally3_int_t magma_tally3_get_spotrf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_sgetrf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_sgetri_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_sgeqp3_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_sgeqrf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_sgeqlf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_sgehrd_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_ssytrd_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_ssytrf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_ssytrf_nopiv_tally3_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_sgelqf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_sgebrd_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_ssygst_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_sgesvd_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_ssygst_nb_m( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_sbulge_nb( magma_tally3_int_t m, magma_tally3_int_t nbthreads );
magma_tally3_int_t magma_tally3_get_sbulge_nb_mgpu( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_sbulge_get_Vblksiz( magma_tally3_int_t m, magma_tally3_int_t nb, magma_tally3_int_t nbthreads );
magma_tally3_int_t magma_tally3_get_sbulge_gcperf();


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
magma_tally3_svrange(
    magma_tally3_int_t k, float *d, magma_tally3_int_t *il, magma_tally3_int_t *iu, float vl, float vu);

void
magma_tally3_sirange(
    magma_tally3_int_t k, magma_tally3_int_t *indxq, magma_tally3_int_t *iil, magma_tally3_int_t *iiu, magma_tally3_int_t il, magma_tally3_int_t iu);
#endif

magma_tally3_int_t
magma_tally3_sgebrd(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float *d, float *e,
    float *tauq, float *taup,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeev(
    magma_tally3_vec_t jobvl, magma_tally3_vec_t jobvr, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    #ifdef COMPLEX
    float *w,
    #else
    float *wr, float *wi,
    #endif
    float *VL, magma_tally3_int_t ldvl,
    float *VR, magma_tally3_int_t ldvr,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgehrd(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    float *A, magma_tally3_int_t lda,
    float *tau,
    float *work, magma_tally3_int_t lwork,
    magma_tally3Float_ptr dT,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgehrd2(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    float *A, magma_tally3_int_t lda,
    float *tau,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgelqf(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float *A,    magma_tally3_int_t lda,
    float *tau,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqlf(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float *A,    magma_tally3_int_t lda,
    float *tau,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqp3(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    magma_tally3_int_t *jpvt, float *tau,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqrf(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float *tau,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqrf_ooc(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float *tau,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqrf4(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float *A,    magma_tally3_int_t lda,
    float *tau,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgesdd(
    magma_tally3_vec_t jobz, magma_tally3_int_t m, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float *s,
    float *U, magma_tally3_int_t ldu,
    float *VT, magma_tally3_int_t ldvt,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally3_int_t *iwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgesv(
    magma_tally3_int_t n, magma_tally3_int_t nrhs,
    float *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    float *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgesv_rbt(
    magma_tally3_bool_t ref, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    float *A, magma_tally3_int_t lda, 
    float *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgesvd(
    magma_tally3_vec_t jobu, magma_tally3_vec_t jobvt, magma_tally3_int_t m, magma_tally3_int_t n,
    float *A,    magma_tally3_int_t lda, float *s,
    float *U,    magma_tally3_int_t ldu,
    float *VT,   magma_tally3_int_t ldvt,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgetf2_nopiv(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgetrf(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgetrf_nopiv(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgetrf_piv(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t NB,
    float *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgetrf2(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssyevd(
    magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float *w,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssyevdx(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssyevdx_2stage(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally3_int_t
magma_tally3_ssyevr(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float vl, float vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, float abstol, magma_tally3_int_t *m,
    float *w,
    float *Z, magma_tally3_int_t ldz,
    magma_tally3_int_t *isuppz,
    float *work, magma_tally3_int_t lwork,
    float *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

// no real [sd] precisions available
magma_tally3_int_t
magma_tally3_ssyevx(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float vl, float vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, float abstol, magma_tally3_int_t *m,
    float *w,
    float *Z, magma_tally3_int_t ldz,
    float *work, magma_tally3_int_t lwork,
    float *rwork, magma_tally3_int_t *iwork,
    magma_tally3_int_t *ifail,
    magma_tally3_int_t *info);
#endif  // COMPLEX

magma_tally3_int_t
magma_tally3_ssygst(
    magma_tally3_int_t itype, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssygvd(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float *B, magma_tally3_int_t ldb,
    float *w, float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssygvdx(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n, float *A, magma_tally3_int_t lda,
    float *B, magma_tally3_int_t ldb,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssygvdx_2stage(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float *B, magma_tally3_int_t ldb,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssygvr(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float *B, magma_tally3_int_t ldb,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    float abstol, magma_tally3_int_t *m, float *w,
    float *Z, magma_tally3_int_t ldz,
    magma_tally3_int_t *isuppz, float *work, magma_tally3_int_t lwork,
    float *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssygvx(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n, float *A, magma_tally3_int_t lda,
    float *B, magma_tally3_int_t ldb,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    float abstol, magma_tally3_int_t *m, float *w,
    float *Z, magma_tally3_int_t ldz,
    float *work, magma_tally3_int_t lwork, float *rwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t *ifail,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssysv(magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
            float *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
            float *B, magma_tally3_int_t ldb,
            magma_tally3_int_t *info );

magma_tally3_int_t
magma_tally3_ssytrd(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float *d, float *e, float *tau,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssytrf(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssytrf_nopiv_tally3(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssytrd_sb2st(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t Vblksiz,
    float *A, magma_tally3_int_t lda,
    float *d, float *e,
    float *V, magma_tally3_int_t ldv,
    float *TAU, magma_tally3_int_t compT,
    float *T, magma_tally3_int_t ldt);

magma_tally3_int_t
magma_tally3_ssytrd_sy2sb(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    float *A, magma_tally3_int_t lda,
    float *tau,
    float *work, magma_tally3_int_t lwork,
    magma_tally3Float_ptr dT,
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
magma_tally3_slasyf_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t *kb,
    float    *hA, magma_tally3_int_t lda,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3Float_ptr dW, magma_tally3_int_t lddw,
    magma_tally3_queue_t queues[], magma_tally3_event_t event[],
    magma_tally3_int_t *info);

magma_tally3_int_t
ssytrf_nopiv_tally3_cpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t ib,
    float *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssytrs_nopiv_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssysv_nopiv_gpu(
    magma_tally3_uplo_t uplo,  magma_tally3_int_t n, magma_tally3_int_t nrhs, 
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb, 
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_slahr2(
    magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3_int_t nb,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dV, magma_tally3_int_t lddv,
    float *A,  magma_tally3_int_t lda,
    float *tau,
    float *T,  magma_tally3_int_t ldt,
    float *Y,  magma_tally3_int_t ldy);

magma_tally3_int_t
magma_tally3_slahru(
    magma_tally3_int_t n, magma_tally3_int_t ihi, magma_tally3_int_t k, magma_tally3_int_t nb,
    float     *A, magma_tally3_int_t lda,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dY, magma_tally3_int_t lddy,
    magma_tally3Float_ptr dV, magma_tally3_int_t lddv,
    magma_tally3Float_ptr dT,
    magma_tally3Float_ptr dwork);

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
magma_tally3_slaqps(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    float *A,  magma_tally3_int_t lda,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt, float *tau, float *vn1, float *vn2,
    float *auxv,
    float *F,  magma_tally3_int_t ldf,
    magma_tally3Float_ptr dF, magma_tally3_int_t lddf );

#ifdef REAL
magma_tally3_int_t
magma_tally3_slaqtrsd(
    magma_tally3_trans_t trans, magma_tally3_int_t n,
    const float *T, magma_tally3_int_t ldt,
    float *x,       magma_tally3_int_t ldx,
    const float *cnorm,
    magma_tally3_int_t *info);
#endif

magma_tally3_int_t
magma_tally3_slatrd(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    float *A, magma_tally3_int_t lda,
    float *e, float *tau,
    float *W, magma_tally3_int_t ldw,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dW, magma_tally3_int_t lddw);

magma_tally3_int_t
magma_tally3_slatrd2(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    float *A,  magma_tally3_int_t lda,
    float *e, float *tau,
    float *W,  magma_tally3_int_t ldw,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dW, magma_tally3_int_t lddw,
    magma_tally3Float_ptr dwork, magma_tally3_int_t ldwork);

#ifdef COMPLEX
magma_tally3_int_t
magma_tally3_slatrsd(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_diag_t diag, magma_tally3_bool_t normin,
    magma_tally3_int_t n, const float *A, magma_tally3_int_t lda,
    float lambda,
    float *x,
    float *scale, float *cnorm,
    magma_tally3_int_t *info);
#endif

magma_tally3_int_t
magma_tally3_slauum(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sposv(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    float *A, magma_tally3_int_t lda,
    float *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_spotrf(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_spotri(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sstedx(
    magma_tally3_range_t range, magma_tally3_int_t n, float vl, float vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, float *d, float *e,
    float *Z, magma_tally3_int_t ldz,
    float *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3Float_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_strevc3(
    magma_tally3_side_t side, magma_tally3_vec_t howmany,
    magma_tally3_int_t *select, magma_tally3_int_t n,
    float *T,  magma_tally3_int_t ldt,
    float *VL, magma_tally3_int_t ldvl,
    float *VR, magma_tally3_int_t ldvr,
    magma_tally3_int_t mm, magma_tally3_int_t *mout,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_strevc3_mt(
    magma_tally3_side_t side, magma_tally3_vec_t howmany,
    magma_tally3_int_t *select, magma_tally3_int_t n,
    float *T,  magma_tally3_int_t ldt,
    float *VL, magma_tally3_int_t ldvl,
    float *VR, magma_tally3_int_t ldvr,
    magma_tally3_int_t mm, magma_tally3_int_t *mout,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_strtri(
    magma_tally3_uplo_t uplo, magma_tally3_diag_t diag, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sorghr(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    float *A, magma_tally3_int_t lda,
    float *tau,
    magma_tally3Float_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sorgqr(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    float *A, magma_tally3_int_t lda,
    float *tau,
    magma_tally3Float_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sorgqr2(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    float *A, magma_tally3_int_t lda,
    float *tau,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sormbr(
    magma_tally3_vect_t vect, magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    float *A, magma_tally3_int_t lda,
    float *tau,
    float *C, magma_tally3_int_t ldc,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sormlq(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    float *A, magma_tally3_int_t lda,
    float *tau,
    float *C, magma_tally3_int_t ldc,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

// not yet implemented
//magma_tally3_int_t magma_tally3_sunmrq( magma_tally3_side_t side, magma_tally3_trans_t trans,
//                          magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
//                          float *A, magma_tally3_int_t lda,
//                          float *tau,
//                          float *C, magma_tally3_int_t ldc,
//                          float *work, magma_tally3_int_t lwork,
//                          magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sormql(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    float *A, magma_tally3_int_t lda,
    float *tau,
    float *C, magma_tally3_int_t ldc,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sormqr(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    float *A, magma_tally3_int_t lda,
    float *tau,
    float *C, magma_tally3_int_t ldc,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sormtr(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float *A,    magma_tally3_int_t lda,
    float *tau,
    float *C,    magma_tally3_int_t ldc,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 function definitions / Data on CPU / Multi-GPU (alphabetical order)
*/
magma_tally3_int_t
magma_tally3_sgeev_m(
    magma_tally3_vec_t jobvl, magma_tally3_vec_t jobvr, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    #ifdef COMPLEX
    float *w,
    #else
    float *wr, float *wi,
    #endif
    float *VL, magma_tally3_int_t ldvl,
    float *VR, magma_tally3_int_t ldvr,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgehrd_m(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    float *A, magma_tally3_int_t lda,
    float *tau,
    float *work, magma_tally3_int_t lwork,
    float *T,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgetrf_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssyevd_m(
    magma_tally3_int_t ngpu,
    magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float *w,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssyevdx_2stage_m(
    magma_tally3_int_t ngpu,
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssyevdx_m(
    magma_tally3_int_t ngpu,
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssygst_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssygvd_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float *B, magma_tally3_int_t ldb,
    float *w,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssygvdx_2stage_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float *B, magma_tally3_int_t ldb,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssygvdx_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float *B, magma_tally3_int_t ldb,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    float *work, magma_tally3_int_t lwork,
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
magma_tally3_slahr2_m(
    magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3_int_t nb,
    float *A, magma_tally3_int_t lda,
    float *tau,
    float *T, magma_tally3_int_t ldt,
    float *Y, magma_tally3_int_t ldy,
    struct sgehrd_data_tally3 *data );

magma_tally3_int_t
magma_tally3_slahru_m(
    magma_tally3_int_t n, magma_tally3_int_t ihi, magma_tally3_int_t k, magma_tally3_int_t nb,
    float *A, magma_tally3_int_t lda,
    struct sgehrd_data_tally3 *data );

magma_tally3_int_t
magma_tally3_spotrf_m(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sstedx_m(
    magma_tally3_int_t ngpu,
    magma_tally3_range_t range, magma_tally3_int_t n, float vl, float vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, float *d, float *e,
    float *Z, magma_tally3_int_t ldz,
    float *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_strsm_m(
    magma_tally3_int_t ngpu,
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transa, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n, float alpha,
    float *A, magma_tally3_int_t lda,
    float *B, magma_tally3_int_t ldb);

magma_tally3_int_t
magma_tally3_sorghr_m(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    float *A, magma_tally3_int_t lda,
    float *tau,
    float *T, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sorgqr_m(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    float *A, magma_tally3_int_t lda,
    float *tau,
    float *T, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sormqr_m(
    magma_tally3_int_t ngpu,
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    float *A,    magma_tally3_int_t lda,
    float *tau,
    float *C,    magma_tally3_int_t ldc,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sormtr_m(
    magma_tally3_int_t ngpu,
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float *A,    magma_tally3_int_t lda,
    float *tau,
    float *C,    magma_tally3_int_t ldc,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 function definitions / Data on GPU (alphabetical order)
*/
magma_tally3_int_t
magma_tally3_sgegqr_gpu(
    magma_tally3_int_t ikind, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dwork, float *work,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgelqf_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *tau,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgels_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb,
    float *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgels3_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb,
    float *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqp3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt, float *tau,
    magma_tally3Float_ptr dwork, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqr2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dtau,
    magma_tally3Float_ptr        dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqr2x_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dtau,
    magma_tally3Float_ptr dT, magma_tally3Float_ptr ddA,
    magma_tally3Float_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqr2x2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dtau,
    magma_tally3Float_ptr dT, magma_tally3Float_ptr ddA,
    magma_tally3Float_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqr2x3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dtau,
    magma_tally3Float_ptr dT,
    magma_tally3Float_ptr ddA,
    magma_tally3Float_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqr2x4_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dtau,
    magma_tally3Float_ptr dT, magma_tally3Float_ptr ddA,
    magma_tally3Float_ptr dwork,
    magma_tally3_queue_t queue,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqrf_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *tau,
    magma_tally3Float_ptr dT,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqrf2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *tau,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqrf2_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dlA[], magma_tally3_int_t ldda,
    float *tau,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqrf3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *tau,
    magma_tally3Float_ptr dT,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqrs_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *tau,
    magma_tally3Float_ptr dT,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb,
    float *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgeqrs3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *tau,
    magma_tally3Float_ptr dT,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb,
    float *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgerbt_gpu(
    magma_tally3_bool_t gen, magma_tally3_int_t n, magma_tally3_int_t nrhs, 
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda, 
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb, 
    float *U, float *V,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgessm_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3_int_t ib,
    magma_tally3_int_t *ipiv,
    magma_tally3Float_ptr dL1, magma_tally3_int_t lddl1,
    magma_tally3Float_ptr dL,  magma_tally3_int_t lddl,
    magma_tally3Float_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgesv_gpu(
    magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgesv_nopiv_gpu( 
    magma_tally3_int_t n, magma_tally3_int_t nrhs, 
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb, 
                 magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgetf2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgetrf_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgetrf_incpiv_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t ib,
    float    *hA, magma_tally3_int_t ldha,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float    *hL, magma_tally3_int_t ldhl,
    magma_tally3Float_ptr dL, magma_tally3_int_t lddl,
    magma_tally3_int_t *ipiv,
    magma_tally3Float_ptr dwork, magma_tally3_int_t lddwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgetrf_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr d_lA[], magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgetrf_nopiv_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgetrf2_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t offset,
    magma_tally3Float_ptr d_lAT[], magma_tally3_int_t lddat, magma_tally3_int_t *ipiv,
    magma_tally3Float_ptr d_lAP[],
    float *W, magma_tally3_int_t ldw,
    magma_tally3_queue_t queues[][2],
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgetri_gpu(
    magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3Float_ptr dwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgetrs_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sgetrs_nopiv_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssyevd_gpu(
    magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *w,
    float *wA,  magma_tally3_int_t ldwa,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssyevdx_gpu(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float vl, float vu,
    magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, float *w,
    float *wA,  magma_tally3_int_t ldwa,
    float *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally3_int_t
magma_tally3_ssyevr_gpu(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float vl, float vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, float abstol, magma_tally3_int_t *m,
    float *w,
    magma_tally3Float_ptr dZ, magma_tally3_int_t lddz,
    magma_tally3_int_t *isuppz,
    float *wA, magma_tally3_int_t ldwa,
    float *wZ, magma_tally3_int_t ldwz,
    float *work, magma_tally3_int_t lwork,
    float *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssyevx_gpu(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float vl, float vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    float abstol, magma_tally3_int_t *m,
    float *w,
    magma_tally3Float_ptr dZ, magma_tally3_int_t lddz,
    float *wA, magma_tally3_int_t ldwa,
    float *wZ, magma_tally3_int_t ldwz,
    float *work, magma_tally3_int_t lwork,
    float *rwork, magma_tally3_int_t *iwork,
    magma_tally3_int_t *ifail,
    magma_tally3_int_t *info);
#endif  // COMPLEX

magma_tally3_int_t
magma_tally3_ssygst_gpu(
    magma_tally3_int_t itype, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssytrd_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *d, float *e, float *tau,
    float *wA,  magma_tally3_int_t ldwa,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssytrd_sy2sb_mgpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    float *A, magma_tally3_int_t lda,
    float *tau,
    float *work, magma_tally3_int_t lwork,
    magma_tally3Float_ptr dAmgpu[], magma_tally3_int_t ldda,
    magma_tally3Float_ptr dTmgpu[], magma_tally3_int_t lddt,
    magma_tally3_int_t ngpu, magma_tally3_int_t distblk,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssytrd_sy2sb_mgpu_spec(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    float *A, magma_tally3_int_t lda,
    float *tau,
    float *work, magma_tally3_int_t lwork,
    magma_tally3Float_ptr dAmgpu[], magma_tally3_int_t ldda,
    magma_tally3Float_ptr dTmgpu[], magma_tally3_int_t lddt,
    magma_tally3_int_t ngpu, magma_tally3_int_t distblk,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssytrd_mgpu(
    magma_tally3_int_t ngpu, magma_tally3_int_t nqueue,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float *A, magma_tally3_int_t lda,
    float *d, float *e, float *tau,
    float *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssytrd2_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *d, float *e, float *tau,
    float *wA,  magma_tally3_int_t ldwa,
    float *work, magma_tally3_int_t lwork,
    magma_tally3Float_ptr dwork, magma_tally3_int_t ldwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ssytrf_nopiv_tally3_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_slabrd_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb,
    float     *A, magma_tally3_int_t lda,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *d, float *e, float *tauq, float *taup,
    float     *X, magma_tally3_int_t ldx,
    magma_tally3Float_ptr dX, magma_tally3_int_t lddx,
    float     *Y, magma_tally3_int_t ldy,
    magma_tally3Float_ptr dY, magma_tally3_int_t lddy );

magma_tally3_int_t
magma_tally3_slaqps_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3Float_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt, float *tau,
    float *vn1, float *vn2,
    magma_tally3Float_ptr dauxv,
    magma_tally3Float_ptr dF, magma_tally3_int_t lddf);

magma_tally3_int_t
magma_tally3_slaqps2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3Float_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt,
    magma_tally3Float_ptr dtau,
    magma_tally3Float_ptr dvn1, magma_tally3Float_ptr dvn2,
    magma_tally3Float_ptr dauxv,
    magma_tally3Float_ptr dF, magma_tally3_int_t lddf);

magma_tally3_int_t
magma_tally3_slaqps3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3Float_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt,
    magma_tally3Float_ptr dtau,
    magma_tally3Float_ptr dvn1, magma_tally3Float_ptr dvn2,
    magma_tally3Float_ptr dauxv,
    magma_tally3Float_ptr dF, magma_tally3_int_t lddf);

magma_tally3_int_t
magma_tally3_slarf_gpu(
    magma_tally3_int_t m,  magma_tally3_int_t n,
    magma_tally3Float_const_ptr dv, magma_tally3Float_const_ptr dtau,
    magma_tally3Float_ptr dC,  magma_tally3_int_t lddc);

magma_tally3_int_t
magma_tally3_slarfb_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans, magma_tally3_direct_t direct, magma_tally3_storev_t storev,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Float_const_ptr dV, magma_tally3_int_t lddv,
    magma_tally3Float_const_ptr dT, magma_tally3_int_t lddt,
    magma_tally3Float_ptr dC,       magma_tally3_int_t lddc,
    magma_tally3Float_ptr dwork,    magma_tally3_int_t ldwork );

magma_tally3_int_t
magma_tally3_slarfb_gpu_gemm(
    magma_tally3_side_t side, magma_tally3_trans_t trans, magma_tally3_direct_t direct, magma_tally3_storev_t storev,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Float_const_ptr dV, magma_tally3_int_t lddv,
    magma_tally3Float_const_ptr dT, magma_tally3_int_t lddt,
    magma_tally3Float_ptr dC,       magma_tally3_int_t lddc,
    magma_tally3Float_ptr dwork,    magma_tally3_int_t ldwork,
    magma_tally3Float_ptr dworkvt,  magma_tally3_int_t ldworkvt);

magma_tally3_int_t
magma_tally3_slarfb2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Float_const_ptr dV, magma_tally3_int_t lddv,
    magma_tally3Float_const_ptr dT, magma_tally3_int_t lddt,
    magma_tally3Float_ptr dC,       magma_tally3_int_t lddc,
    magma_tally3Float_ptr dwork,    magma_tally3_int_t ldwork );

magma_tally3_int_t
magma_tally3_slatrd_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t nb0,
    float *A,  magma_tally3_int_t lda,
    float *e, float *tau,
    float    *W,       magma_tally3_int_t ldw,
    magma_tally3Float_ptr dA[],    magma_tally3_int_t ldda, magma_tally3_int_t offset,
    magma_tally3Float_ptr dW[],    magma_tally3_int_t lddw,
    float    *hwork,   magma_tally3_int_t lhwork,
    magma_tally3Float_ptr dwork[], magma_tally3_int_t ldwork,
    magma_tally3_queue_t queues[] );

magma_tally3_int_t
magma_tally3_slauum_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sposv_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_spotf2_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_spotrf_gpu(
    magma_tally3_uplo_t uplo,  magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_spotrf_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Float_ptr d_lA[], magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_spotrf_mgpu_right(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Float_ptr d_lA[], magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_spotrf3_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t off_i, magma_tally3_int_t off_j, magma_tally3_int_t nb,
    magma_tally3Float_ptr d_lA[],  magma_tally3_int_t ldda,
    magma_tally3Float_ptr d_lP[],  magma_tally3_int_t lddp,
    float *A, magma_tally3_int_t lda, magma_tally3_int_t h,
    magma_tally3_queue_t queues[][3], magma_tally3_event_t events[][5],
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_spotri_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_spotrs_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sssssm_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m1, magma_tally3_int_t n1,
    magma_tally3_int_t m2, magma_tally3_int_t n2, magma_tally3_int_t k, magma_tally3_int_t ib,
    magma_tally3Float_ptr dA1, magma_tally3_int_t ldda1,
    magma_tally3Float_ptr dA2, magma_tally3_int_t ldda2,
    magma_tally3Float_ptr dL1, magma_tally3_int_t lddl1,
    magma_tally3Float_ptr dL2, magma_tally3_int_t lddl2,
    magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_strtri_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_diag_t diag, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_stsqrt_gpu(
    magma_tally3_int_t *m, magma_tally3_int_t *n,
    float *A1, float *A2, magma_tally3_int_t *lda,
    float *tau,
    float *work, magma_tally3_int_t *lwork,
    magma_tally3Float_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_ststrf_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t ib, magma_tally3_int_t nb,
    float    *hU, magma_tally3_int_t ldhu,
    magma_tally3Float_ptr dU, magma_tally3_int_t lddu,
    float    *hA, magma_tally3_int_t ldha,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float    *hL, magma_tally3_int_t ldhl,
    magma_tally3Float_ptr dL, magma_tally3_int_t lddl,
    magma_tally3_int_t *ipiv,
    float *hwork, magma_tally3_int_t ldhwork,
    magma_tally3Float_ptr dwork, magma_tally3_int_t lddwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sorgqr_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *tau,
    magma_tally3Float_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sormql2_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *tau,
    magma_tally3Float_ptr dC, magma_tally3_int_t lddc,
    float *wA, magma_tally3_int_t ldwa,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sormqr_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *tau,
    magma_tally3Float_ptr dC, magma_tally3_int_t lddc,
    float *hwork, magma_tally3_int_t lwork,
    magma_tally3Float_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sormqr2_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *tau,
    magma_tally3Float_ptr dC, magma_tally3_int_t lddc,
    float    *wA, magma_tally3_int_t ldwa,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_sormtr_gpu(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *tau,
    magma_tally3Float_ptr dC, magma_tally3_int_t lddc,
    float    *wA, magma_tally3_int_t ldwa,
    magma_tally3_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 utility function definitions
*/

extern const float MAGMA_tally3_S_NAN;
extern const float MAGMA_tally3_S_INF;

magma_tally3_int_t
magma_tally3_snan_inf(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    const float *A, magma_tally3_int_t lda,
    magma_tally3_int_t *cnt_nan,
    magma_tally3_int_t *cnt_inf );

magma_tally3_int_t
magma_tally3_snan_inf_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *cnt_nan,
    magma_tally3_int_t *cnt_inf );

void magma_tally3_sprint(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const float *A, magma_tally3_int_t lda );

void magma_tally3_sprint_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr dA, magma_tally3_int_t ldda );

void spanel_to_q_tally3(
    magma_tally3_uplo_t uplo, magma_tally3_int_t ib,
    float *A, magma_tally3_int_t lda,
    float *work );

void sq_to_panel_tally3(
    magma_tally3_uplo_t uplo, magma_tally3_int_t ib,
    float *A, magma_tally3_int_t lda,
    float *work );

#ifdef __cplusplus
}
#endif

#undef REAL

#endif /* MAGMA_tally3_S_H */
