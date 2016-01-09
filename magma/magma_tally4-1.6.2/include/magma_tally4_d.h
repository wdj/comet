/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally4_z.h normal z -> d, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally4_D_H
#define MAGMA_tally4_D_H

#include "magma_tally4_types.h"
#include "magma_tally4_dgehrd_m.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 Auxiliary functions to get the NB used
*/
#ifdef REAL
magma_tally4_int_t magma_tally4_get_dlaex3_m_nb();       // defined in dlaex3_m.cpp
#endif

magma_tally4_int_t magma_tally4_get_dpotrf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_dgetrf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_dgetri_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_dgeqp3_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_dgeqrf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_dgeqlf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_dgehrd_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_dsytrd_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_dsytrf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_dsytrf_nopiv_tally4_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_dgelqf_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_dgebrd_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_dsygst_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_dgesvd_nb( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_dsygst_nb_m( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_get_dbulge_nb( magma_tally4_int_t m, magma_tally4_int_t nbthreads );
magma_tally4_int_t magma_tally4_get_dbulge_nb_mgpu( magma_tally4_int_t m );
magma_tally4_int_t magma_tally4_dbulge_get_Vblksiz( magma_tally4_int_t m, magma_tally4_int_t nb, magma_tally4_int_t nbthreads );
magma_tally4_int_t magma_tally4_get_dbulge_gcperf();


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
magma_tally4_dvrange(
    magma_tally4_int_t k, double *d, magma_tally4_int_t *il, magma_tally4_int_t *iu, double vl, double vu);

void
magma_tally4_dirange(
    magma_tally4_int_t k, magma_tally4_int_t *indxq, magma_tally4_int_t *iil, magma_tally4_int_t *iiu, magma_tally4_int_t il, magma_tally4_int_t iu);
#endif

magma_tally4_int_t
magma_tally4_dgebrd(
    magma_tally4_int_t m, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *d, double *e,
    double *tauq, double *taup,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeev(
    magma_tally4_vec_t jobvl, magma_tally4_vec_t jobvr, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    #ifdef COMPLEX
    double *w,
    #else
    double *wr, double *wi,
    #endif
    double *VL, magma_tally4_int_t ldvl,
    double *VR, magma_tally4_int_t ldvr,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgehrd(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    double *A, magma_tally4_int_t lda,
    double *tau,
    double *work, magma_tally4_int_t lwork,
    magma_tally4Double_ptr dT,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgehrd2(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    double *A, magma_tally4_int_t lda,
    double *tau,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgelqf(
    magma_tally4_int_t m, magma_tally4_int_t n,
    double *A,    magma_tally4_int_t lda,
    double *tau,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqlf(
    magma_tally4_int_t m, magma_tally4_int_t n,
    double *A,    magma_tally4_int_t lda,
    double *tau,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqp3(
    magma_tally4_int_t m, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    magma_tally4_int_t *jpvt, double *tau,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqrf(
    magma_tally4_int_t m, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *tau,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqrf_ooc(
    magma_tally4_int_t m, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *tau,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqrf4(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n,
    double *A,    magma_tally4_int_t lda,
    double *tau,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgesdd(
    magma_tally4_vec_t jobz, magma_tally4_int_t m, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *s,
    double *U, magma_tally4_int_t ldu,
    double *VT, magma_tally4_int_t ldvt,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally4_int_t *iwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgesv(
    magma_tally4_int_t n, magma_tally4_int_t nrhs,
    double *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    double *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgesv_rbt(
    magma_tally4_bool_t ref, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    double *A, magma_tally4_int_t lda, 
    double *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgesvd(
    magma_tally4_vec_t jobu, magma_tally4_vec_t jobvt, magma_tally4_int_t m, magma_tally4_int_t n,
    double *A,    magma_tally4_int_t lda, double *s,
    double *U,    magma_tally4_int_t ldu,
    double *VT,   magma_tally4_int_t ldvt,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgetf2_nopiv(
    magma_tally4_int_t m, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgetrf(
    magma_tally4_int_t m, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgetrf_nopiv(
    magma_tally4_int_t m, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgetrf_piv(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t NB,
    double *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgetrf2(
    magma_tally4_int_t m, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsyevd(
    magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *w,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsyevdx(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsyevdx_2stage(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally4_int_t
magma_tally4_dsyevr(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double vl, double vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, double abstol, magma_tally4_int_t *m,
    double *w,
    double *Z, magma_tally4_int_t ldz,
    magma_tally4_int_t *isuppz,
    double *work, magma_tally4_int_t lwork,
    double *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

// no real [sd] precisions available
magma_tally4_int_t
magma_tally4_dsyevx(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double vl, double vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, double abstol, magma_tally4_int_t *m,
    double *w,
    double *Z, magma_tally4_int_t ldz,
    double *work, magma_tally4_int_t lwork,
    double *rwork, magma_tally4_int_t *iwork,
    magma_tally4_int_t *ifail,
    magma_tally4_int_t *info);
#endif  // COMPLEX

magma_tally4_int_t
magma_tally4_dsygst(
    magma_tally4_int_t itype, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsygvd(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *B, magma_tally4_int_t ldb,
    double *w, double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsygvdx(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n, double *A, magma_tally4_int_t lda,
    double *B, magma_tally4_int_t ldb,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsygvdx_2stage(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *B, magma_tally4_int_t ldb,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsygvr(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *B, magma_tally4_int_t ldb,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    double abstol, magma_tally4_int_t *m, double *w,
    double *Z, magma_tally4_int_t ldz,
    magma_tally4_int_t *isuppz, double *work, magma_tally4_int_t lwork,
    double *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsygvx(
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n, double *A, magma_tally4_int_t lda,
    double *B, magma_tally4_int_t ldb,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    double abstol, magma_tally4_int_t *m, double *w,
    double *Z, magma_tally4_int_t ldz,
    double *work, magma_tally4_int_t lwork, double *rwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t *ifail,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsysv(magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
            double *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
            double *B, magma_tally4_int_t ldb,
            magma_tally4_int_t *info );

magma_tally4_int_t
magma_tally4_dsytrd(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *d, double *e, double *tau,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsytrf(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsytrf_nopiv_tally4(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsytrd_sb2st(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t Vblksiz,
    double *A, magma_tally4_int_t lda,
    double *d, double *e,
    double *V, magma_tally4_int_t ldv,
    double *TAU, magma_tally4_int_t compT,
    double *T, magma_tally4_int_t ldt);

magma_tally4_int_t
magma_tally4_dsytrd_sy2sb(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    double *A, magma_tally4_int_t lda,
    double *tau,
    double *work, magma_tally4_int_t lwork,
    magma_tally4Double_ptr dT,
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
magma_tally4_dlasyf_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t *kb,
    double    *hA, magma_tally4_int_t lda,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4Double_ptr dW, magma_tally4_int_t lddw,
    magma_tally4_queue_t queues[], magma_tally4_event_t event[],
    magma_tally4_int_t *info);

magma_tally4_int_t
dsytrf_nopiv_tally4_cpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t ib,
    double *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsytrs_nopiv_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsysv_nopiv_gpu(
    magma_tally4_uplo_t uplo,  magma_tally4_int_t n, magma_tally4_int_t nrhs, 
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb, 
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dlahr2(
    magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4_int_t nb,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dV, magma_tally4_int_t lddv,
    double *A,  magma_tally4_int_t lda,
    double *tau,
    double *T,  magma_tally4_int_t ldt,
    double *Y,  magma_tally4_int_t ldy);

magma_tally4_int_t
magma_tally4_dlahru(
    magma_tally4_int_t n, magma_tally4_int_t ihi, magma_tally4_int_t k, magma_tally4_int_t nb,
    double     *A, magma_tally4_int_t lda,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dY, magma_tally4_int_t lddy,
    magma_tally4Double_ptr dV, magma_tally4_int_t lddv,
    magma_tally4Double_ptr dT,
    magma_tally4Double_ptr dwork);

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
magma_tally4_dlaqps(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t offset,
    magma_tally4_int_t nb, magma_tally4_int_t *kb,
    double *A,  magma_tally4_int_t lda,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt, double *tau, double *vn1, double *vn2,
    double *auxv,
    double *F,  magma_tally4_int_t ldf,
    magma_tally4Double_ptr dF, magma_tally4_int_t lddf );

#ifdef REAL
magma_tally4_int_t
magma_tally4_dlaqtrsd(
    magma_tally4_trans_t trans, magma_tally4_int_t n,
    const double *T, magma_tally4_int_t ldt,
    double *x,       magma_tally4_int_t ldx,
    const double *cnorm,
    magma_tally4_int_t *info);
#endif

magma_tally4_int_t
magma_tally4_dlatrd(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    double *A, magma_tally4_int_t lda,
    double *e, double *tau,
    double *W, magma_tally4_int_t ldw,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dW, magma_tally4_int_t lddw);

magma_tally4_int_t
magma_tally4_dlatrd2(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    double *A,  magma_tally4_int_t lda,
    double *e, double *tau,
    double *W,  magma_tally4_int_t ldw,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dW, magma_tally4_int_t lddw,
    magma_tally4Double_ptr dwork, magma_tally4_int_t ldwork);

#ifdef COMPLEX
magma_tally4_int_t
magma_tally4_dlatrsd(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_diag_t diag, magma_tally4_bool_t normin,
    magma_tally4_int_t n, const double *A, magma_tally4_int_t lda,
    double lambda,
    double *x,
    double *scale, double *cnorm,
    magma_tally4_int_t *info);
#endif

magma_tally4_int_t
magma_tally4_dlauum(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dposv(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    double *A, magma_tally4_int_t lda,
    double *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dpotrf(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dpotri(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dstedx(
    magma_tally4_range_t range, magma_tally4_int_t n, double vl, double vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, double *d, double *e,
    double *Z, magma_tally4_int_t ldz,
    double *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4Double_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dtrevc3(
    magma_tally4_side_t side, magma_tally4_vec_t howmany,
    magma_tally4_int_t *select, magma_tally4_int_t n,
    double *T,  magma_tally4_int_t ldt,
    double *VL, magma_tally4_int_t ldvl,
    double *VR, magma_tally4_int_t ldvr,
    magma_tally4_int_t mm, magma_tally4_int_t *mout,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dtrevc3_mt(
    magma_tally4_side_t side, magma_tally4_vec_t howmany,
    magma_tally4_int_t *select, magma_tally4_int_t n,
    double *T,  magma_tally4_int_t ldt,
    double *VL, magma_tally4_int_t ldvl,
    double *VR, magma_tally4_int_t ldvr,
    magma_tally4_int_t mm, magma_tally4_int_t *mout,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dtrtri(
    magma_tally4_uplo_t uplo, magma_tally4_diag_t diag, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dorghr(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    double *A, magma_tally4_int_t lda,
    double *tau,
    magma_tally4Double_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dorgqr(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    double *A, magma_tally4_int_t lda,
    double *tau,
    magma_tally4Double_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dorgqr2(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    double *A, magma_tally4_int_t lda,
    double *tau,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dormbr(
    magma_tally4_vect_t vect, magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    double *A, magma_tally4_int_t lda,
    double *tau,
    double *C, magma_tally4_int_t ldc,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dormlq(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    double *A, magma_tally4_int_t lda,
    double *tau,
    double *C, magma_tally4_int_t ldc,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

// not yet implemented
//magma_tally4_int_t magma_tally4_dunmrq( magma_tally4_side_t side, magma_tally4_trans_t trans,
//                          magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
//                          double *A, magma_tally4_int_t lda,
//                          double *tau,
//                          double *C, magma_tally4_int_t ldc,
//                          double *work, magma_tally4_int_t lwork,
//                          magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dormql(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    double *A, magma_tally4_int_t lda,
    double *tau,
    double *C, magma_tally4_int_t ldc,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dormqr(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    double *A, magma_tally4_int_t lda,
    double *tau,
    double *C, magma_tally4_int_t ldc,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dormtr(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n,
    double *A,    magma_tally4_int_t lda,
    double *tau,
    double *C,    magma_tally4_int_t ldc,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 function definitions / Data on CPU / Multi-GPU (alphabetical order)
*/
magma_tally4_int_t
magma_tally4_dgeev_m(
    magma_tally4_vec_t jobvl, magma_tally4_vec_t jobvr, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    #ifdef COMPLEX
    double *w,
    #else
    double *wr, double *wi,
    #endif
    double *VL, magma_tally4_int_t ldvl,
    double *VR, magma_tally4_int_t ldvr,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgehrd_m(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    double *A, magma_tally4_int_t lda,
    double *tau,
    double *work, magma_tally4_int_t lwork,
    double *T,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgetrf_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsyevd_m(
    magma_tally4_int_t ngpu,
    magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *w,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsyevdx_2stage_m(
    magma_tally4_int_t ngpu,
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsyevdx_m(
    magma_tally4_int_t ngpu,
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsygst_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t itype, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsygvd_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *B, magma_tally4_int_t ldb,
    double *w,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsygvdx_2stage_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *B, magma_tally4_int_t ldb,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsygvdx_m(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t itype, magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *B, magma_tally4_int_t ldb,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    double *work, magma_tally4_int_t lwork,
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
magma_tally4_dlahr2_m(
    magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4_int_t nb,
    double *A, magma_tally4_int_t lda,
    double *tau,
    double *T, magma_tally4_int_t ldt,
    double *Y, magma_tally4_int_t ldy,
    struct dgehrd_data_tally4 *data );

magma_tally4_int_t
magma_tally4_dlahru_m(
    magma_tally4_int_t n, magma_tally4_int_t ihi, magma_tally4_int_t k, magma_tally4_int_t nb,
    double *A, magma_tally4_int_t lda,
    struct dgehrd_data_tally4 *data );

magma_tally4_int_t
magma_tally4_dpotrf_m(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dstedx_m(
    magma_tally4_int_t ngpu,
    magma_tally4_range_t range, magma_tally4_int_t n, double vl, double vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, double *d, double *e,
    double *Z, magma_tally4_int_t ldz,
    double *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dtrsm_m(
    magma_tally4_int_t ngpu,
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t transa, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n, double alpha,
    double *A, magma_tally4_int_t lda,
    double *B, magma_tally4_int_t ldb);

magma_tally4_int_t
magma_tally4_dorghr_m(
    magma_tally4_int_t n, magma_tally4_int_t ilo, magma_tally4_int_t ihi,
    double *A, magma_tally4_int_t lda,
    double *tau,
    double *T, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dorgqr_m(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    double *A, magma_tally4_int_t lda,
    double *tau,
    double *T, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dormqr_m(
    magma_tally4_int_t ngpu,
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    double *A,    magma_tally4_int_t lda,
    double *tau,
    double *C,    magma_tally4_int_t ldc,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dormtr_m(
    magma_tally4_int_t ngpu,
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n,
    double *A,    magma_tally4_int_t lda,
    double *tau,
    double *C,    magma_tally4_int_t ldc,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 function definitions / Data on GPU (alphabetical order)
*/
magma_tally4_int_t
magma_tally4_dgegqr_gpu(
    magma_tally4_int_t ikind, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dwork, double *work,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgelqf_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *tau,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgels_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    double *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgels3_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    double *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqp3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt, double *tau,
    magma_tally4Double_ptr dwork, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqr2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dtau,
    magma_tally4Double_ptr        dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqr2x_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dtau,
    magma_tally4Double_ptr dT, magma_tally4Double_ptr ddA,
    magma_tally4Double_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqr2x2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dtau,
    magma_tally4Double_ptr dT, magma_tally4Double_ptr ddA,
    magma_tally4Double_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqr2x3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dtau,
    magma_tally4Double_ptr dT,
    magma_tally4Double_ptr ddA,
    magma_tally4Double_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqr2x4_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dtau,
    magma_tally4Double_ptr dT, magma_tally4Double_ptr ddA,
    magma_tally4Double_ptr dwork,
    magma_tally4_queue_t queue,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqrf_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *tau,
    magma_tally4Double_ptr dT,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqrf2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *tau,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqrf2_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dlA[], magma_tally4_int_t ldda,
    double *tau,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqrf3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *tau,
    magma_tally4Double_ptr dT,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqrs_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *tau,
    magma_tally4Double_ptr dT,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    double *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgeqrs3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *tau,
    magma_tally4Double_ptr dT,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    double *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgerbt_gpu(
    magma_tally4_bool_t gen, magma_tally4_int_t n, magma_tally4_int_t nrhs, 
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda, 
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb, 
    double *U, double *V,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgessm_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4_int_t ib,
    magma_tally4_int_t *ipiv,
    magma_tally4Double_ptr dL1, magma_tally4_int_t lddl1,
    magma_tally4Double_ptr dL,  magma_tally4_int_t lddl,
    magma_tally4Double_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgesv_gpu(
    magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgesv_nopiv_gpu( 
    magma_tally4_int_t n, magma_tally4_int_t nrhs, 
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb, 
                 magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgetf2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgetrf_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgetrf_incpiv_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t ib,
    double    *hA, magma_tally4_int_t ldha,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double    *hL, magma_tally4_int_t ldhl,
    magma_tally4Double_ptr dL, magma_tally4_int_t lddl,
    magma_tally4_int_t *ipiv,
    magma_tally4Double_ptr dwork, magma_tally4_int_t lddwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgetrf_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr d_lA[], magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgetrf_nopiv_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgetrf2_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t offset,
    magma_tally4Double_ptr d_lAT[], magma_tally4_int_t lddat, magma_tally4_int_t *ipiv,
    magma_tally4Double_ptr d_lAP[],
    double *W, magma_tally4_int_t ldw,
    magma_tally4_queue_t queues[][2],
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgetri_gpu(
    magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4Double_ptr dwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgetrs_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dgetrs_nopiv_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsyevd_gpu(
    magma_tally4_vec_t jobz, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *w,
    double *wA,  magma_tally4_int_t ldwa,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsyevdx_gpu(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double vl, double vu,
    magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *m, double *w,
    double *wA,  magma_tally4_int_t ldwa,
    double *work, magma_tally4_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally4_int_t lrwork,
    #endif
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally4_int_t
magma_tally4_dsyevr_gpu(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double vl, double vu,
    magma_tally4_int_t il, magma_tally4_int_t iu, double abstol, magma_tally4_int_t *m,
    double *w,
    magma_tally4Double_ptr dZ, magma_tally4_int_t lddz,
    magma_tally4_int_t *isuppz,
    double *wA, magma_tally4_int_t ldwa,
    double *wZ, magma_tally4_int_t ldwz,
    double *work, magma_tally4_int_t lwork,
    double *rwork, magma_tally4_int_t lrwork,
    magma_tally4_int_t *iwork, magma_tally4_int_t liwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsyevx_gpu(
    magma_tally4_vec_t jobz, magma_tally4_range_t range, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double vl, double vu, magma_tally4_int_t il, magma_tally4_int_t iu,
    double abstol, magma_tally4_int_t *m,
    double *w,
    magma_tally4Double_ptr dZ, magma_tally4_int_t lddz,
    double *wA, magma_tally4_int_t ldwa,
    double *wZ, magma_tally4_int_t ldwz,
    double *work, magma_tally4_int_t lwork,
    double *rwork, magma_tally4_int_t *iwork,
    magma_tally4_int_t *ifail,
    magma_tally4_int_t *info);
#endif  // COMPLEX

magma_tally4_int_t
magma_tally4_dsygst_gpu(
    magma_tally4_int_t itype, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsytrd_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *d, double *e, double *tau,
    double *wA,  magma_tally4_int_t ldwa,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsytrd_sy2sb_mgpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    double *A, magma_tally4_int_t lda,
    double *tau,
    double *work, magma_tally4_int_t lwork,
    magma_tally4Double_ptr dAmgpu[], magma_tally4_int_t ldda,
    magma_tally4Double_ptr dTmgpu[], magma_tally4_int_t lddt,
    magma_tally4_int_t ngpu, magma_tally4_int_t distblk,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsytrd_sy2sb_mgpu_spec(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    double *A, magma_tally4_int_t lda,
    double *tau,
    double *work, magma_tally4_int_t lwork,
    magma_tally4Double_ptr dAmgpu[], magma_tally4_int_t ldda,
    magma_tally4Double_ptr dTmgpu[], magma_tally4_int_t lddt,
    magma_tally4_int_t ngpu, magma_tally4_int_t distblk,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsytrd_mgpu(
    magma_tally4_int_t ngpu, magma_tally4_int_t nqueue,
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *d, double *e, double *tau,
    double *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsytrd2_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *d, double *e, double *tau,
    double *wA,  magma_tally4_int_t ldwa,
    double *work, magma_tally4_int_t lwork,
    magma_tally4Double_ptr dwork, magma_tally4_int_t ldwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dsytrf_nopiv_tally4_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dlabrd_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nb,
    double     *A, magma_tally4_int_t lda,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *d, double *e, double *tauq, double *taup,
    double     *X, magma_tally4_int_t ldx,
    magma_tally4Double_ptr dX, magma_tally4_int_t lddx,
    double     *Y, magma_tally4_int_t ldy,
    magma_tally4Double_ptr dY, magma_tally4_int_t lddy );

magma_tally4_int_t
magma_tally4_dlaqps_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t offset,
    magma_tally4_int_t nb, magma_tally4_int_t *kb,
    magma_tally4Double_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt, double *tau,
    double *vn1, double *vn2,
    magma_tally4Double_ptr dauxv,
    magma_tally4Double_ptr dF, magma_tally4_int_t lddf);

magma_tally4_int_t
magma_tally4_dlaqps2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t offset,
    magma_tally4_int_t nb, magma_tally4_int_t *kb,
    magma_tally4Double_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt,
    magma_tally4Double_ptr dtau,
    magma_tally4Double_ptr dvn1, magma_tally4Double_ptr dvn2,
    magma_tally4Double_ptr dauxv,
    magma_tally4Double_ptr dF, magma_tally4_int_t lddf);

magma_tally4_int_t
magma_tally4_dlaqps3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t offset,
    magma_tally4_int_t nb, magma_tally4_int_t *kb,
    magma_tally4Double_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *jpvt,
    magma_tally4Double_ptr dtau,
    magma_tally4Double_ptr dvn1, magma_tally4Double_ptr dvn2,
    magma_tally4Double_ptr dauxv,
    magma_tally4Double_ptr dF, magma_tally4_int_t lddf);

magma_tally4_int_t
magma_tally4_dlarf_gpu(
    magma_tally4_int_t m,  magma_tally4_int_t n,
    magma_tally4Double_const_ptr dv, magma_tally4Double_const_ptr dtau,
    magma_tally4Double_ptr dC,  magma_tally4_int_t lddc);

magma_tally4_int_t
magma_tally4_dlarfb_gpu(
    magma_tally4_side_t side, magma_tally4_trans_t trans, magma_tally4_direct_t direct, magma_tally4_storev_t storev,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Double_const_ptr dV, magma_tally4_int_t lddv,
    magma_tally4Double_const_ptr dT, magma_tally4_int_t lddt,
    magma_tally4Double_ptr dC,       magma_tally4_int_t lddc,
    magma_tally4Double_ptr dwork,    magma_tally4_int_t ldwork );

magma_tally4_int_t
magma_tally4_dlarfb_gpu_gemm(
    magma_tally4_side_t side, magma_tally4_trans_t trans, magma_tally4_direct_t direct, magma_tally4_storev_t storev,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Double_const_ptr dV, magma_tally4_int_t lddv,
    magma_tally4Double_const_ptr dT, magma_tally4_int_t lddt,
    magma_tally4Double_ptr dC,       magma_tally4_int_t lddc,
    magma_tally4Double_ptr dwork,    magma_tally4_int_t ldwork,
    magma_tally4Double_ptr dworkvt,  magma_tally4_int_t ldworkvt);

magma_tally4_int_t
magma_tally4_dlarfb2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Double_const_ptr dV, magma_tally4_int_t lddv,
    magma_tally4Double_const_ptr dT, magma_tally4_int_t lddt,
    magma_tally4Double_ptr dC,       magma_tally4_int_t lddc,
    magma_tally4Double_ptr dwork,    magma_tally4_int_t ldwork );

magma_tally4_int_t
magma_tally4_dlatrd_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t nb0,
    double *A,  magma_tally4_int_t lda,
    double *e, double *tau,
    double    *W,       magma_tally4_int_t ldw,
    magma_tally4Double_ptr dA[],    magma_tally4_int_t ldda, magma_tally4_int_t offset,
    magma_tally4Double_ptr dW[],    magma_tally4_int_t lddw,
    double    *hwork,   magma_tally4_int_t lhwork,
    magma_tally4Double_ptr dwork[], magma_tally4_int_t ldwork,
    magma_tally4_queue_t queues[] );

magma_tally4_int_t
magma_tally4_dlauum_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dposv_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dpotf2_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dpotrf_gpu(
    magma_tally4_uplo_t uplo,  magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dpotrf_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Double_ptr d_lA[], magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dpotrf_mgpu_right(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Double_ptr d_lA[], magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dpotrf3_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t off_i, magma_tally4_int_t off_j, magma_tally4_int_t nb,
    magma_tally4Double_ptr d_lA[],  magma_tally4_int_t ldda,
    magma_tally4Double_ptr d_lP[],  magma_tally4_int_t lddp,
    double *A, magma_tally4_int_t lda, magma_tally4_int_t h,
    magma_tally4_queue_t queues[][3], magma_tally4_event_t events[][5],
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dpotri_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dpotrs_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dssssm_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m1, magma_tally4_int_t n1,
    magma_tally4_int_t m2, magma_tally4_int_t n2, magma_tally4_int_t k, magma_tally4_int_t ib,
    magma_tally4Double_ptr dA1, magma_tally4_int_t ldda1,
    magma_tally4Double_ptr dA2, magma_tally4_int_t ldda2,
    magma_tally4Double_ptr dL1, magma_tally4_int_t lddl1,
    magma_tally4Double_ptr dL2, magma_tally4_int_t lddl2,
    magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dtrtri_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_diag_t diag, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dtsqrt_gpu(
    magma_tally4_int_t *m, magma_tally4_int_t *n,
    double *A1, double *A2, magma_tally4_int_t *lda,
    double *tau,
    double *work, magma_tally4_int_t *lwork,
    magma_tally4Double_ptr dwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dtstrf_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t ib, magma_tally4_int_t nb,
    double    *hU, magma_tally4_int_t ldhu,
    magma_tally4Double_ptr dU, magma_tally4_int_t lddu,
    double    *hA, magma_tally4_int_t ldha,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double    *hL, magma_tally4_int_t ldhl,
    magma_tally4Double_ptr dL, magma_tally4_int_t lddl,
    magma_tally4_int_t *ipiv,
    double *hwork, magma_tally4_int_t ldhwork,
    magma_tally4Double_ptr dwork, magma_tally4_int_t lddwork,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dorgqr_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *tau,
    magma_tally4Double_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dormql2_gpu(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *tau,
    magma_tally4Double_ptr dC, magma_tally4_int_t lddc,
    double *wA, magma_tally4_int_t ldwa,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dormqr_gpu(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *tau,
    magma_tally4Double_ptr dC, magma_tally4_int_t lddc,
    double *hwork, magma_tally4_int_t lwork,
    magma_tally4Double_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dormqr2_gpu(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *tau,
    magma_tally4Double_ptr dC, magma_tally4_int_t lddc,
    double    *wA, magma_tally4_int_t ldwa,
    magma_tally4_int_t *info);

magma_tally4_int_t
magma_tally4_dormtr_gpu(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *tau,
    magma_tally4Double_ptr dC, magma_tally4_int_t lddc,
    double    *wA, magma_tally4_int_t ldwa,
    magma_tally4_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally4 utility function definitions
*/

extern const double MAGMA_tally4_D_NAN;
extern const double MAGMA_tally4_D_INF;

magma_tally4_int_t
magma_tally4_dnan_inf(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    const double *A, magma_tally4_int_t lda,
    magma_tally4_int_t *cnt_nan,
    magma_tally4_int_t *cnt_inf );

magma_tally4_int_t
magma_tally4_dnan_inf_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *cnt_nan,
    magma_tally4_int_t *cnt_inf );

void magma_tally4_dprint(
    magma_tally4_int_t m, magma_tally4_int_t n,
    const double *A, magma_tally4_int_t lda );

void magma_tally4_dprint_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_const_ptr dA, magma_tally4_int_t ldda );

void dpanel_to_q_tally4(
    magma_tally4_uplo_t uplo, magma_tally4_int_t ib,
    double *A, magma_tally4_int_t lda,
    double *work );

void dq_to_panel_tally4(
    magma_tally4_uplo_t uplo, magma_tally4_int_t ib,
    double *A, magma_tally4_int_t lda,
    double *work );

#ifdef __cplusplus
}
#endif

#undef REAL

#endif /* MAGMA_tally4_D_H */
