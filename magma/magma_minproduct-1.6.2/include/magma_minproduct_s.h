/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproduct_z.h normal z -> s, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_minproduct_S_H
#define MAGMA_minproduct_S_H

#include "magma_minproduct_types.h"
#include "magma_minproduct_sgehrd_m.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct Auxiliary functions to get the NB used
*/
#ifdef REAL
magma_minproduct_int_t magma_minproduct_get_slaex3_m_nb();       // defined in slaex3_m.cpp
#endif

magma_minproduct_int_t magma_minproduct_get_spotrf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_sgetrf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_sgetri_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_sgeqp3_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_sgeqrf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_sgeqlf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_sgehrd_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_ssytrd_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_ssytrf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_ssytrf_nopiv_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_sgelqf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_sgebrd_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_ssygst_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_sgesvd_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_ssygst_nb_m( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_sbulge_nb( magma_minproduct_int_t m, magma_minproduct_int_t nbthreads );
magma_minproduct_int_t magma_minproduct_get_sbulge_nb_mgpu( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_sbulge_get_Vblksiz( magma_minproduct_int_t m, magma_minproduct_int_t nb, magma_minproduct_int_t nbthreads );
magma_minproduct_int_t magma_minproduct_get_sbulge_gcperf();


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
magma_minproduct_svrange(
    magma_minproduct_int_t k, float *d, magma_minproduct_int_t *il, magma_minproduct_int_t *iu, float vl, float vu);

void
magma_minproduct_sirange(
    magma_minproduct_int_t k, magma_minproduct_int_t *indxq, magma_minproduct_int_t *iil, magma_minproduct_int_t *iiu, magma_minproduct_int_t il, magma_minproduct_int_t iu);
#endif

magma_minproduct_int_t
magma_minproduct_sgebrd(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float *d, float *e,
    float *tauq, float *taup,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeev(
    magma_minproduct_vec_t jobvl, magma_minproduct_vec_t jobvr, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    #ifdef COMPLEX
    float *w,
    #else
    float *wr, float *wi,
    #endif
    float *VL, magma_minproduct_int_t ldvl,
    float *VR, magma_minproduct_int_t ldvr,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgehrd(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    float *work, magma_minproduct_int_t lwork,
    magma_minproductFloat_ptr dT,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgehrd2(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgelqf(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float *A,    magma_minproduct_int_t lda,
    float *tau,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqlf(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float *A,    magma_minproduct_int_t lda,
    float *tau,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqp3(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *jpvt, float *tau,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqrf(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqrf_ooc(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqrf4(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float *A,    magma_minproduct_int_t lda,
    float *tau,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgesdd(
    magma_minproduct_vec_t jobz, magma_minproduct_int_t m, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float *s,
    float *U, magma_minproduct_int_t ldu,
    float *VT, magma_minproduct_int_t ldvt,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_minproduct_int_t *iwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgesv(
    magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    float *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    float *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgesv_rbt(
    magma_minproduct_bool_t ref, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    float *A, magma_minproduct_int_t lda, 
    float *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgesvd(
    magma_minproduct_vec_t jobu, magma_minproduct_vec_t jobvt, magma_minproduct_int_t m, magma_minproduct_int_t n,
    float *A,    magma_minproduct_int_t lda, float *s,
    float *U,    magma_minproduct_int_t ldu,
    float *VT,   magma_minproduct_int_t ldvt,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgetf2_nopiv(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgetrf(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgetrf_nopiv(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgetrf_piv(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t NB,
    float *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgetrf2(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssyevd(
    magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float *w,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssyevdx(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssyevdx_2stage(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_minproduct_int_t
magma_minproduct_ssyevr(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float vl, float vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, float abstol, magma_minproduct_int_t *m,
    float *w,
    float *Z, magma_minproduct_int_t ldz,
    magma_minproduct_int_t *isuppz,
    float *work, magma_minproduct_int_t lwork,
    float *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

// no real [sd] precisions available
magma_minproduct_int_t
magma_minproduct_ssyevx(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float vl, float vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, float abstol, magma_minproduct_int_t *m,
    float *w,
    float *Z, magma_minproduct_int_t ldz,
    float *work, magma_minproduct_int_t lwork,
    float *rwork, magma_minproduct_int_t *iwork,
    magma_minproduct_int_t *ifail,
    magma_minproduct_int_t *info);
#endif  // COMPLEX

magma_minproduct_int_t
magma_minproduct_ssygst(
    magma_minproduct_int_t itype, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssygvd(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float *B, magma_minproduct_int_t ldb,
    float *w, float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssygvdx(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n, float *A, magma_minproduct_int_t lda,
    float *B, magma_minproduct_int_t ldb,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssygvdx_2stage(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float *B, magma_minproduct_int_t ldb,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssygvr(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float *B, magma_minproduct_int_t ldb,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    float abstol, magma_minproduct_int_t *m, float *w,
    float *Z, magma_minproduct_int_t ldz,
    magma_minproduct_int_t *isuppz, float *work, magma_minproduct_int_t lwork,
    float *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssygvx(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n, float *A, magma_minproduct_int_t lda,
    float *B, magma_minproduct_int_t ldb,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    float abstol, magma_minproduct_int_t *m, float *w,
    float *Z, magma_minproduct_int_t ldz,
    float *work, magma_minproduct_int_t lwork, float *rwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t *ifail,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssysv(magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
            float *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
            float *B, magma_minproduct_int_t ldb,
            magma_minproduct_int_t *info );

magma_minproduct_int_t
magma_minproduct_ssytrd(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float *d, float *e, float *tau,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssytrf(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssytrf_nopiv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssytrd_sb2st(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t Vblksiz,
    float *A, magma_minproduct_int_t lda,
    float *d, float *e,
    float *V, magma_minproduct_int_t ldv,
    float *TAU, magma_minproduct_int_t compT,
    float *T, magma_minproduct_int_t ldt);

magma_minproduct_int_t
magma_minproduct_ssytrd_sy2sb(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    float *work, magma_minproduct_int_t lwork,
    magma_minproductFloat_ptr dT,
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
magma_minproduct_slasyf_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    float    *hA, magma_minproduct_int_t lda,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductFloat_ptr dW, magma_minproduct_int_t lddw,
    magma_minproduct_queue_t queues[], magma_minproduct_event_t event[],
    magma_minproduct_int_t *info);

magma_minproduct_int_t
ssytrf_nopiv_cpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t ib,
    float *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssytrs_nopiv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssysv_nopiv_gpu(
    magma_minproduct_uplo_t uplo,  magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb, 
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_slahr2(
    magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dV, magma_minproduct_int_t lddv,
    float *A,  magma_minproduct_int_t lda,
    float *tau,
    float *T,  magma_minproduct_int_t ldt,
    float *Y,  magma_minproduct_int_t ldy);

magma_minproduct_int_t
magma_minproduct_slahru(
    magma_minproduct_int_t n, magma_minproduct_int_t ihi, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    float     *A, magma_minproduct_int_t lda,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dY, magma_minproduct_int_t lddy,
    magma_minproductFloat_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductFloat_ptr dT,
    magma_minproductFloat_ptr dwork);

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
magma_minproduct_slaqps(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t offset,
    magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    float *A,  magma_minproduct_int_t lda,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt, float *tau, float *vn1, float *vn2,
    float *auxv,
    float *F,  magma_minproduct_int_t ldf,
    magma_minproductFloat_ptr dF, magma_minproduct_int_t lddf );

#ifdef REAL
magma_minproduct_int_t
magma_minproduct_slaqtrsd(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n,
    const float *T, magma_minproduct_int_t ldt,
    float *x,       magma_minproduct_int_t ldx,
    const float *cnorm,
    magma_minproduct_int_t *info);
#endif

magma_minproduct_int_t
magma_minproduct_slatrd(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    float *A, magma_minproduct_int_t lda,
    float *e, float *tau,
    float *W, magma_minproduct_int_t ldw,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dW, magma_minproduct_int_t lddw);

magma_minproduct_int_t
magma_minproduct_slatrd2(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    float *A,  magma_minproduct_int_t lda,
    float *e, float *tau,
    float *W,  magma_minproduct_int_t ldw,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dW, magma_minproduct_int_t lddw,
    magma_minproductFloat_ptr dwork, magma_minproduct_int_t ldwork);

#ifdef COMPLEX
magma_minproduct_int_t
magma_minproduct_slatrsd(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_diag_t diag, magma_minproduct_bool_t normin,
    magma_minproduct_int_t n, const float *A, magma_minproduct_int_t lda,
    float lambda,
    float *x,
    float *scale, float *cnorm,
    magma_minproduct_int_t *info);
#endif

magma_minproduct_int_t
magma_minproduct_slauum(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sposv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    float *A, magma_minproduct_int_t lda,
    float *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_spotrf(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_spotri(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sstedx(
    magma_minproduct_range_t range, magma_minproduct_int_t n, float vl, float vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, float *d, float *e,
    float *Z, magma_minproduct_int_t ldz,
    float *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_strevc3(
    magma_minproduct_side_t side, magma_minproduct_vec_t howmany,
    magma_minproduct_int_t *select, magma_minproduct_int_t n,
    float *T,  magma_minproduct_int_t ldt,
    float *VL, magma_minproduct_int_t ldvl,
    float *VR, magma_minproduct_int_t ldvr,
    magma_minproduct_int_t mm, magma_minproduct_int_t *mout,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_strevc3_mt(
    magma_minproduct_side_t side, magma_minproduct_vec_t howmany,
    magma_minproduct_int_t *select, magma_minproduct_int_t n,
    float *T,  magma_minproduct_int_t ldt,
    float *VL, magma_minproduct_int_t ldvl,
    float *VR, magma_minproduct_int_t ldvr,
    magma_minproduct_int_t mm, magma_minproduct_int_t *mout,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_strtri(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sorghr(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    magma_minproductFloat_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sorgqr(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    magma_minproductFloat_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sorgqr2(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sormbr(
    magma_minproduct_vect_t vect, magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    float *C, magma_minproduct_int_t ldc,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sormlq(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    float *C, magma_minproduct_int_t ldc,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

// not yet implemented
//magma_minproduct_int_t magma_minproduct_sunmrq( magma_minproduct_side_t side, magma_minproduct_trans_t trans,
//                          magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
//                          float *A, magma_minproduct_int_t lda,
//                          float *tau,
//                          float *C, magma_minproduct_int_t ldc,
//                          float *work, magma_minproduct_int_t lwork,
//                          magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sormql(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    float *C, magma_minproduct_int_t ldc,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sormqr(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    float *C, magma_minproduct_int_t ldc,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sormtr(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float *A,    magma_minproduct_int_t lda,
    float *tau,
    float *C,    magma_minproduct_int_t ldc,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct function definitions / Data on CPU / Multi-GPU (alphabetical order)
*/
magma_minproduct_int_t
magma_minproduct_sgeev_m(
    magma_minproduct_vec_t jobvl, magma_minproduct_vec_t jobvr, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    #ifdef COMPLEX
    float *w,
    #else
    float *wr, float *wi,
    #endif
    float *VL, magma_minproduct_int_t ldvl,
    float *VR, magma_minproduct_int_t ldvr,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgehrd_m(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    float *work, magma_minproduct_int_t lwork,
    float *T,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgetrf_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssyevd_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float *w,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssyevdx_2stage_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssyevdx_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssygst_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t itype, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssygvd_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float *B, magma_minproduct_int_t ldb,
    float *w,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssygvdx_2stage_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float *B, magma_minproduct_int_t ldb,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssygvdx_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float *B, magma_minproduct_int_t ldb,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    float *work, magma_minproduct_int_t lwork,
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
magma_minproduct_slahr2_m(
    magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    float *T, magma_minproduct_int_t ldt,
    float *Y, magma_minproduct_int_t ldy,
    struct sgehrd_data *data );

magma_minproduct_int_t
magma_minproduct_slahru_m(
    magma_minproduct_int_t n, magma_minproduct_int_t ihi, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    float *A, magma_minproduct_int_t lda,
    struct sgehrd_data *data );

magma_minproduct_int_t
magma_minproduct_spotrf_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sstedx_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_range_t range, magma_minproduct_int_t n, float vl, float vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, float *d, float *e,
    float *Z, magma_minproduct_int_t ldz,
    float *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_strsm_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transa, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n, float alpha,
    float *A, magma_minproduct_int_t lda,
    float *B, magma_minproduct_int_t ldb);

magma_minproduct_int_t
magma_minproduct_sorghr_m(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    float *T, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sorgqr_m(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    float *T, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sormqr_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float *A,    magma_minproduct_int_t lda,
    float *tau,
    float *C,    magma_minproduct_int_t ldc,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sormtr_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float *A,    magma_minproduct_int_t lda,
    float *tau,
    float *C,    magma_minproduct_int_t ldc,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct function definitions / Data on GPU (alphabetical order)
*/
magma_minproduct_int_t
magma_minproduct_sgegqr_gpu(
    magma_minproduct_int_t ikind, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dwork, float *work,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgelqf_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *tau,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgels_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb,
    float *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgels3_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb,
    float *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqp3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt, float *tau,
    magma_minproductFloat_ptr dwork, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqr2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dtau,
    magma_minproductFloat_ptr        dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqr2x_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dtau,
    magma_minproductFloat_ptr dT, magma_minproductFloat_ptr ddA,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqr2x2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dtau,
    magma_minproductFloat_ptr dT, magma_minproductFloat_ptr ddA,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqr2x3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dtau,
    magma_minproductFloat_ptr dT,
    magma_minproductFloat_ptr ddA,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqr2x4_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dtau,
    magma_minproductFloat_ptr dT, magma_minproductFloat_ptr ddA,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_queue_t queue,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqrf_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *tau,
    magma_minproductFloat_ptr dT,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqrf2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *tau,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqrf2_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dlA[], magma_minproduct_int_t ldda,
    float *tau,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqrf3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *tau,
    magma_minproductFloat_ptr dT,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqrs_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *tau,
    magma_minproductFloat_ptr dT,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb,
    float *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgeqrs3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *tau,
    magma_minproductFloat_ptr dT,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb,
    float *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgerbt_gpu(
    magma_minproduct_bool_t gen, magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda, 
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb, 
    float *U, float *V,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgessm_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t ib,
    magma_minproduct_int_t *ipiv,
    magma_minproductFloat_ptr dL1, magma_minproduct_int_t lddl1,
    magma_minproductFloat_ptr dL,  magma_minproduct_int_t lddl,
    magma_minproductFloat_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgesv_gpu(
    magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgesv_nopiv_gpu( 
    magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb, 
                 magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgetf2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgetrf_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgetrf_incpiv_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t ib,
    float    *hA, magma_minproduct_int_t ldha,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float    *hL, magma_minproduct_int_t ldhl,
    magma_minproductFloat_ptr dL, magma_minproduct_int_t lddl,
    magma_minproduct_int_t *ipiv,
    magma_minproductFloat_ptr dwork, magma_minproduct_int_t lddwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgetrf_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr d_lA[], magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgetrf_nopiv_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgetrf2_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t offset,
    magma_minproductFloat_ptr d_lAT[], magma_minproduct_int_t lddat, magma_minproduct_int_t *ipiv,
    magma_minproductFloat_ptr d_lAP[],
    float *W, magma_minproduct_int_t ldw,
    magma_minproduct_queue_t queues[][2],
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgetri_gpu(
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductFloat_ptr dwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgetrs_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sgetrs_nopiv_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssyevd_gpu(
    magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *w,
    float *wA,  magma_minproduct_int_t ldwa,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssyevdx_gpu(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float vl, float vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, float *w,
    float *wA,  magma_minproduct_int_t ldwa,
    float *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    float *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_minproduct_int_t
magma_minproduct_ssyevr_gpu(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float vl, float vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, float abstol, magma_minproduct_int_t *m,
    float *w,
    magma_minproductFloat_ptr dZ, magma_minproduct_int_t lddz,
    magma_minproduct_int_t *isuppz,
    float *wA, magma_minproduct_int_t ldwa,
    float *wZ, magma_minproduct_int_t ldwz,
    float *work, magma_minproduct_int_t lwork,
    float *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssyevx_gpu(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float vl, float vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    float abstol, magma_minproduct_int_t *m,
    float *w,
    magma_minproductFloat_ptr dZ, magma_minproduct_int_t lddz,
    float *wA, magma_minproduct_int_t ldwa,
    float *wZ, magma_minproduct_int_t ldwz,
    float *work, magma_minproduct_int_t lwork,
    float *rwork, magma_minproduct_int_t *iwork,
    magma_minproduct_int_t *ifail,
    magma_minproduct_int_t *info);
#endif  // COMPLEX

magma_minproduct_int_t
magma_minproduct_ssygst_gpu(
    magma_minproduct_int_t itype, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssytrd_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *d, float *e, float *tau,
    float *wA,  magma_minproduct_int_t ldwa,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssytrd_sy2sb_mgpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    float *work, magma_minproduct_int_t lwork,
    magma_minproductFloat_ptr dAmgpu[], magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dTmgpu[], magma_minproduct_int_t lddt,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t distblk,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssytrd_sy2sb_mgpu_spec(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    float *work, magma_minproduct_int_t lwork,
    magma_minproductFloat_ptr dAmgpu[], magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dTmgpu[], magma_minproduct_int_t lddt,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t distblk,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssytrd_mgpu(
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nqueue,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float *A, magma_minproduct_int_t lda,
    float *d, float *e, float *tau,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssytrd2_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *d, float *e, float *tau,
    float *wA,  magma_minproduct_int_t ldwa,
    float *work, magma_minproduct_int_t lwork,
    magma_minproductFloat_ptr dwork, magma_minproduct_int_t ldwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ssytrf_nopiv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_slabrd_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    float     *A, magma_minproduct_int_t lda,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *d, float *e, float *tauq, float *taup,
    float     *X, magma_minproduct_int_t ldx,
    magma_minproductFloat_ptr dX, magma_minproduct_int_t lddx,
    float     *Y, magma_minproduct_int_t ldy,
    magma_minproductFloat_ptr dY, magma_minproduct_int_t lddy );

magma_minproduct_int_t
magma_minproduct_slaqps_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t offset,
    magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    magma_minproductFloat_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt, float *tau,
    float *vn1, float *vn2,
    magma_minproductFloat_ptr dauxv,
    magma_minproductFloat_ptr dF, magma_minproduct_int_t lddf);

magma_minproduct_int_t
magma_minproduct_slaqps2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t offset,
    magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    magma_minproductFloat_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt,
    magma_minproductFloat_ptr dtau,
    magma_minproductFloat_ptr dvn1, magma_minproductFloat_ptr dvn2,
    magma_minproductFloat_ptr dauxv,
    magma_minproductFloat_ptr dF, magma_minproduct_int_t lddf);

magma_minproduct_int_t
magma_minproduct_slaqps3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t offset,
    magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    magma_minproductFloat_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt,
    magma_minproductFloat_ptr dtau,
    magma_minproductFloat_ptr dvn1, magma_minproductFloat_ptr dvn2,
    magma_minproductFloat_ptr dauxv,
    magma_minproductFloat_ptr dF, magma_minproduct_int_t lddf);

magma_minproduct_int_t
magma_minproduct_slarf_gpu(
    magma_minproduct_int_t m,  magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dv, magma_minproductFloat_const_ptr dtau,
    magma_minproductFloat_ptr dC,  magma_minproduct_int_t lddc);

magma_minproduct_int_t
magma_minproduct_slarfb_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_direct_t direct, magma_minproduct_storev_t storev,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloat_const_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductFloat_const_ptr dT, magma_minproduct_int_t lddt,
    magma_minproductFloat_ptr dC,       magma_minproduct_int_t lddc,
    magma_minproductFloat_ptr dwork,    magma_minproduct_int_t ldwork );

magma_minproduct_int_t
magma_minproduct_slarfb_gpu_gemm(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_direct_t direct, magma_minproduct_storev_t storev,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloat_const_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductFloat_const_ptr dT, magma_minproduct_int_t lddt,
    magma_minproductFloat_ptr dC,       magma_minproduct_int_t lddc,
    magma_minproductFloat_ptr dwork,    magma_minproduct_int_t ldwork,
    magma_minproductFloat_ptr dworkvt,  magma_minproduct_int_t ldworkvt);

magma_minproduct_int_t
magma_minproduct_slarfb2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloat_const_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductFloat_const_ptr dT, magma_minproduct_int_t lddt,
    magma_minproductFloat_ptr dC,       magma_minproduct_int_t lddc,
    magma_minproductFloat_ptr dwork,    magma_minproduct_int_t ldwork );

magma_minproduct_int_t
magma_minproduct_slatrd_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t nb0,
    float *A,  magma_minproduct_int_t lda,
    float *e, float *tau,
    float    *W,       magma_minproduct_int_t ldw,
    magma_minproductFloat_ptr dA[],    magma_minproduct_int_t ldda, magma_minproduct_int_t offset,
    magma_minproductFloat_ptr dW[],    magma_minproduct_int_t lddw,
    float    *hwork,   magma_minproduct_int_t lhwork,
    magma_minproductFloat_ptr dwork[], magma_minproduct_int_t ldwork,
    magma_minproduct_queue_t queues[] );

magma_minproduct_int_t
magma_minproduct_slauum_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sposv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_spotf2_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_spotrf_gpu(
    magma_minproduct_uplo_t uplo,  magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_spotrf_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloat_ptr d_lA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_spotrf_mgpu_right(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloat_ptr d_lA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_spotrf3_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t off_i, magma_minproduct_int_t off_j, magma_minproduct_int_t nb,
    magma_minproductFloat_ptr d_lA[],  magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr d_lP[],  magma_minproduct_int_t lddp,
    float *A, magma_minproduct_int_t lda, magma_minproduct_int_t h,
    magma_minproduct_queue_t queues[][3], magma_minproduct_event_t events[][5],
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_spotri_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_spotrs_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sssssm_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m1, magma_minproduct_int_t n1,
    magma_minproduct_int_t m2, magma_minproduct_int_t n2, magma_minproduct_int_t k, magma_minproduct_int_t ib,
    magma_minproductFloat_ptr dA1, magma_minproduct_int_t ldda1,
    magma_minproductFloat_ptr dA2, magma_minproduct_int_t ldda2,
    magma_minproductFloat_ptr dL1, magma_minproduct_int_t lddl1,
    magma_minproductFloat_ptr dL2, magma_minproduct_int_t lddl2,
    magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_strtri_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_stsqrt_gpu(
    magma_minproduct_int_t *m, magma_minproduct_int_t *n,
    float *A1, float *A2, magma_minproduct_int_t *lda,
    float *tau,
    float *work, magma_minproduct_int_t *lwork,
    magma_minproductFloat_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_ststrf_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t ib, magma_minproduct_int_t nb,
    float    *hU, magma_minproduct_int_t ldhu,
    magma_minproductFloat_ptr dU, magma_minproduct_int_t lddu,
    float    *hA, magma_minproduct_int_t ldha,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float    *hL, magma_minproduct_int_t ldhl,
    magma_minproductFloat_ptr dL, magma_minproduct_int_t lddl,
    magma_minproduct_int_t *ipiv,
    float *hwork, magma_minproduct_int_t ldhwork,
    magma_minproductFloat_ptr dwork, magma_minproduct_int_t lddwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sorgqr_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *tau,
    magma_minproductFloat_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sormql2_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *tau,
    magma_minproductFloat_ptr dC, magma_minproduct_int_t lddc,
    float *wA, magma_minproduct_int_t ldwa,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sormqr_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *tau,
    magma_minproductFloat_ptr dC, magma_minproduct_int_t lddc,
    float *hwork, magma_minproduct_int_t lwork,
    magma_minproductFloat_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sormqr2_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *tau,
    magma_minproductFloat_ptr dC, magma_minproduct_int_t lddc,
    float    *wA, magma_minproduct_int_t ldwa,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_sormtr_gpu(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *tau,
    magma_minproductFloat_ptr dC, magma_minproduct_int_t lddc,
    float    *wA, magma_minproduct_int_t ldwa,
    magma_minproduct_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct utility function definitions
*/

extern const float MAGMA_minproduct_S_NAN;
extern const float MAGMA_minproduct_S_INF;

magma_minproduct_int_t
magma_minproduct_snan_inf(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    const float *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *cnt_nan,
    magma_minproduct_int_t *cnt_inf );

magma_minproduct_int_t
magma_minproduct_snan_inf_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *cnt_nan,
    magma_minproduct_int_t *cnt_inf );

void magma_minproduct_sprint(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const float *A, magma_minproduct_int_t lda );

void magma_minproduct_sprint_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda );

void spanel_to_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t ib,
    float *A, magma_minproduct_int_t lda,
    float *work );

void sq_to_panel(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t ib,
    float *A, magma_minproduct_int_t lda,
    float *work );

#ifdef __cplusplus
}
#endif

#undef REAL

#endif /* MAGMA_minproduct_S_H */
