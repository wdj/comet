/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally2_z.h normal z -> d, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally2_D_H
#define MAGMA_tally2_D_H

#include "magma_tally2_types.h"
#include "magma_tally2_dgehrd_m.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 Auxiliary functions to get the NB used
*/
#ifdef REAL
magma_tally2_int_t magma_tally2_get_dlaex3_m_nb();       // defined in dlaex3_m.cpp
#endif

magma_tally2_int_t magma_tally2_get_dpotrf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_dgetrf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_dgetri_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_dgeqp3_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_dgeqrf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_dgeqlf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_dgehrd_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_dsytrd_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_dsytrf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_dsytrf_nopiv_tally2_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_dgelqf_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_dgebrd_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_dsygst_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_dgesvd_nb( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_dsygst_nb_m( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_get_dbulge_nb( magma_tally2_int_t m, magma_tally2_int_t nbthreads );
magma_tally2_int_t magma_tally2_get_dbulge_nb_mgpu( magma_tally2_int_t m );
magma_tally2_int_t magma_tally2_dbulge_get_Vblksiz( magma_tally2_int_t m, magma_tally2_int_t nb, magma_tally2_int_t nbthreads );
magma_tally2_int_t magma_tally2_get_dbulge_gcperf();


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 function definitions / Data on CPU (alphabetical order)
*/

#ifdef REAL
// only applicable to real [sd] precisions
void
magma_tally2_dmove_eig(
    magma_tally2_range_t range, magma_tally2_int_t n, double *w,
    magma_tally2_int_t *il, magma_tally2_int_t *iu, double vl, double vu, magma_tally2_int_t *m);

// defined in dlaex3.cpp
void
magma_tally2_dvrange(
    magma_tally2_int_t k, double *d, magma_tally2_int_t *il, magma_tally2_int_t *iu, double vl, double vu);

void
magma_tally2_dirange(
    magma_tally2_int_t k, magma_tally2_int_t *indxq, magma_tally2_int_t *iil, magma_tally2_int_t *iiu, magma_tally2_int_t il, magma_tally2_int_t iu);
#endif

magma_tally2_int_t
magma_tally2_dgebrd(
    magma_tally2_int_t m, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double *d, double *e,
    double *tauq, double *taup,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeev(
    magma_tally2_vec_t jobvl, magma_tally2_vec_t jobvr, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    #ifdef COMPLEX
    double *w,
    #else
    double *wr, double *wi,
    #endif
    double *VL, magma_tally2_int_t ldvl,
    double *VR, magma_tally2_int_t ldvr,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgehrd(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    double *A, magma_tally2_int_t lda,
    double *tau,
    double *work, magma_tally2_int_t lwork,
    magma_tally2Double_ptr dT,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgehrd2(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    double *A, magma_tally2_int_t lda,
    double *tau,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgelqf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    double *A,    magma_tally2_int_t lda,
    double *tau,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqlf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    double *A,    magma_tally2_int_t lda,
    double *tau,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqp3(
    magma_tally2_int_t m, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    magma_tally2_int_t *jpvt, double *tau,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqrf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double *tau,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqrf_ooc(
    magma_tally2_int_t m, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double *tau,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqrf4(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double *A,    magma_tally2_int_t lda,
    double *tau,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgesdd(
    magma_tally2_vec_t jobz, magma_tally2_int_t m, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double *s,
    double *U, magma_tally2_int_t ldu,
    double *VT, magma_tally2_int_t ldvt,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally2_int_t *iwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgesv(
    magma_tally2_int_t n, magma_tally2_int_t nrhs,
    double *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    double *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgesv_rbt(
    magma_tally2_bool_t ref, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    double *A, magma_tally2_int_t lda, 
    double *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgesvd(
    magma_tally2_vec_t jobu, magma_tally2_vec_t jobvt, magma_tally2_int_t m, magma_tally2_int_t n,
    double *A,    magma_tally2_int_t lda, double *s,
    double *U,    magma_tally2_int_t ldu,
    double *VT,   magma_tally2_int_t ldvt,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgetf2_nopiv(
    magma_tally2_int_t m, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgetrf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgetrf_nopiv(
    magma_tally2_int_t m, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgetrf_piv(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t NB,
    double *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgetrf2(
    magma_tally2_int_t m, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsyevd(
    magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double *w,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsyevdx(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsyevdx_2stage(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally2_int_t
magma_tally2_dsyevr(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double vl, double vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, double abstol, magma_tally2_int_t *m,
    double *w,
    double *Z, magma_tally2_int_t ldz,
    magma_tally2_int_t *isuppz,
    double *work, magma_tally2_int_t lwork,
    double *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

// no real [sd] precisions available
magma_tally2_int_t
magma_tally2_dsyevx(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double vl, double vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, double abstol, magma_tally2_int_t *m,
    double *w,
    double *Z, magma_tally2_int_t ldz,
    double *work, magma_tally2_int_t lwork,
    double *rwork, magma_tally2_int_t *iwork,
    magma_tally2_int_t *ifail,
    magma_tally2_int_t *info);
#endif  // COMPLEX

magma_tally2_int_t
magma_tally2_dsygst(
    magma_tally2_int_t itype, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsygvd(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double *B, magma_tally2_int_t ldb,
    double *w, double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsygvdx(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n, double *A, magma_tally2_int_t lda,
    double *B, magma_tally2_int_t ldb,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsygvdx_2stage(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double *B, magma_tally2_int_t ldb,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsygvr(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double *B, magma_tally2_int_t ldb,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    double abstol, magma_tally2_int_t *m, double *w,
    double *Z, magma_tally2_int_t ldz,
    magma_tally2_int_t *isuppz, double *work, magma_tally2_int_t lwork,
    double *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsygvx(
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n, double *A, magma_tally2_int_t lda,
    double *B, magma_tally2_int_t ldb,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    double abstol, magma_tally2_int_t *m, double *w,
    double *Z, magma_tally2_int_t ldz,
    double *work, magma_tally2_int_t lwork, double *rwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t *ifail,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsysv(magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
            double *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
            double *B, magma_tally2_int_t ldb,
            magma_tally2_int_t *info );

magma_tally2_int_t
magma_tally2_dsytrd(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double *d, double *e, double *tau,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsytrf(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsytrf_nopiv_tally2(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsytrd_sb2st(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t Vblksiz,
    double *A, magma_tally2_int_t lda,
    double *d, double *e,
    double *V, magma_tally2_int_t ldv,
    double *TAU, magma_tally2_int_t compT,
    double *T, magma_tally2_int_t ldt);

magma_tally2_int_t
magma_tally2_dsytrd_sy2sb(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    double *A, magma_tally2_int_t lda,
    double *tau,
    double *work, magma_tally2_int_t lwork,
    magma_tally2Double_ptr dT,
    magma_tally2_int_t *info);

#ifdef REAL
// only applicable to real [sd] precisions
magma_tally2_int_t
magma_tally2_dlaex0(
    magma_tally2_int_t n, double *d, double *e,
    double *Q, magma_tally2_int_t ldq,
    double *work, magma_tally2_int_t *iwork,
    magma_tally2Double_ptr dwork,
    magma_tally2_range_t range, double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dlaex1(
    magma_tally2_int_t n, double *d,
    double *Q, magma_tally2_int_t ldq,
    magma_tally2_int_t *indxq, double rho, magma_tally2_int_t cutpnt,
    double *work, magma_tally2_int_t *iwork,
    magma_tally2Double_ptr dwork,
    magma_tally2_range_t range, double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dlaex3(
    magma_tally2_int_t k, magma_tally2_int_t n, magma_tally2_int_t n1, double *d,
    double *Q, magma_tally2_int_t ldq,
    double rho,
    double *dlamda, double *Q2, magma_tally2_int_t *indx,
    magma_tally2_int_t *ctot, double *w, double *s, magma_tally2_int_t *indxq,
    magma_tally2Double_ptr dwork,
    magma_tally2_range_t range, double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *info);
#endif  // REAL

magma_tally2_int_t
magma_tally2_dlasyf_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t *kb,
    double    *hA, magma_tally2_int_t lda,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2Double_ptr dW, magma_tally2_int_t lddw,
    magma_tally2_queue_t queues[], magma_tally2_event_t event[],
    magma_tally2_int_t *info);

magma_tally2_int_t
dsytrf_nopiv_tally2_cpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t ib,
    double *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsytrs_nopiv_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsysv_nopiv_gpu(
    magma_tally2_uplo_t uplo,  magma_tally2_int_t n, magma_tally2_int_t nrhs, 
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb, 
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dlahr2(
    magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2_int_t nb,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dV, magma_tally2_int_t lddv,
    double *A,  magma_tally2_int_t lda,
    double *tau,
    double *T,  magma_tally2_int_t ldt,
    double *Y,  magma_tally2_int_t ldy);

magma_tally2_int_t
magma_tally2_dlahru(
    magma_tally2_int_t n, magma_tally2_int_t ihi, magma_tally2_int_t k, magma_tally2_int_t nb,
    double     *A, magma_tally2_int_t lda,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dY, magma_tally2_int_t lddy,
    magma_tally2Double_ptr dV, magma_tally2_int_t lddv,
    magma_tally2Double_ptr dT,
    magma_tally2Double_ptr dwork);

#ifdef REAL
magma_tally2_int_t
magma_tally2_dlaln2(
    magma_tally2_int_t trans, magma_tally2_int_t na, magma_tally2_int_t nw,
    double smin, double ca, const double *A, magma_tally2_int_t lda,
    double d1, double d2,   const double *B, magma_tally2_int_t ldb,
    double wr, double wi, double *X, magma_tally2_int_t ldx,
    double *scale, double *xnorm,
    magma_tally2_int_t *info);
#endif

magma_tally2_int_t
magma_tally2_dlaqps(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t offset,
    magma_tally2_int_t nb, magma_tally2_int_t *kb,
    double *A,  magma_tally2_int_t lda,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt, double *tau, double *vn1, double *vn2,
    double *auxv,
    double *F,  magma_tally2_int_t ldf,
    magma_tally2Double_ptr dF, magma_tally2_int_t lddf );

#ifdef REAL
magma_tally2_int_t
magma_tally2_dlaqtrsd(
    magma_tally2_trans_t trans, magma_tally2_int_t n,
    const double *T, magma_tally2_int_t ldt,
    double *x,       magma_tally2_int_t ldx,
    const double *cnorm,
    magma_tally2_int_t *info);
#endif

magma_tally2_int_t
magma_tally2_dlatrd(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    double *A, magma_tally2_int_t lda,
    double *e, double *tau,
    double *W, magma_tally2_int_t ldw,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dW, magma_tally2_int_t lddw);

magma_tally2_int_t
magma_tally2_dlatrd2(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    double *A,  magma_tally2_int_t lda,
    double *e, double *tau,
    double *W,  magma_tally2_int_t ldw,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dW, magma_tally2_int_t lddw,
    magma_tally2Double_ptr dwork, magma_tally2_int_t ldwork);

#ifdef COMPLEX
magma_tally2_int_t
magma_tally2_dlatrsd(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_diag_t diag, magma_tally2_bool_t normin,
    magma_tally2_int_t n, const double *A, magma_tally2_int_t lda,
    double lambda,
    double *x,
    double *scale, double *cnorm,
    magma_tally2_int_t *info);
#endif

magma_tally2_int_t
magma_tally2_dlauum(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dposv(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    double *A, magma_tally2_int_t lda,
    double *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dpotrf(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dpotri(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dstedx(
    magma_tally2_range_t range, magma_tally2_int_t n, double vl, double vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, double *d, double *e,
    double *Z, magma_tally2_int_t ldz,
    double *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2Double_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dtrevc3(
    magma_tally2_side_t side, magma_tally2_vec_t howmany,
    magma_tally2_int_t *select, magma_tally2_int_t n,
    double *T,  magma_tally2_int_t ldt,
    double *VL, magma_tally2_int_t ldvl,
    double *VR, magma_tally2_int_t ldvr,
    magma_tally2_int_t mm, magma_tally2_int_t *mout,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dtrevc3_mt(
    magma_tally2_side_t side, magma_tally2_vec_t howmany,
    magma_tally2_int_t *select, magma_tally2_int_t n,
    double *T,  magma_tally2_int_t ldt,
    double *VL, magma_tally2_int_t ldvl,
    double *VR, magma_tally2_int_t ldvr,
    magma_tally2_int_t mm, magma_tally2_int_t *mout,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dtrtri(
    magma_tally2_uplo_t uplo, magma_tally2_diag_t diag, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dorghr(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    double *A, magma_tally2_int_t lda,
    double *tau,
    magma_tally2Double_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dorgqr(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    double *A, magma_tally2_int_t lda,
    double *tau,
    magma_tally2Double_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dorgqr2(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    double *A, magma_tally2_int_t lda,
    double *tau,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dormbr(
    magma_tally2_vect_t vect, magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    double *A, magma_tally2_int_t lda,
    double *tau,
    double *C, magma_tally2_int_t ldc,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dormlq(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    double *A, magma_tally2_int_t lda,
    double *tau,
    double *C, magma_tally2_int_t ldc,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

// not yet implemented
//magma_tally2_int_t magma_tally2_dunmrq( magma_tally2_side_t side, magma_tally2_trans_t trans,
//                          magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
//                          double *A, magma_tally2_int_t lda,
//                          double *tau,
//                          double *C, magma_tally2_int_t ldc,
//                          double *work, magma_tally2_int_t lwork,
//                          magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dormql(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    double *A, magma_tally2_int_t lda,
    double *tau,
    double *C, magma_tally2_int_t ldc,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dormqr(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    double *A, magma_tally2_int_t lda,
    double *tau,
    double *C, magma_tally2_int_t ldc,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dormtr(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double *A,    magma_tally2_int_t lda,
    double *tau,
    double *C,    magma_tally2_int_t ldc,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 function definitions / Data on CPU / Multi-GPU (alphabetical order)
*/
magma_tally2_int_t
magma_tally2_dgeev_m(
    magma_tally2_vec_t jobvl, magma_tally2_vec_t jobvr, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    #ifdef COMPLEX
    double *w,
    #else
    double *wr, double *wi,
    #endif
    double *VL, magma_tally2_int_t ldvl,
    double *VR, magma_tally2_int_t ldvr,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgehrd_m(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    double *A, magma_tally2_int_t lda,
    double *tau,
    double *work, magma_tally2_int_t lwork,
    double *T,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgetrf_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsyevd_m(
    magma_tally2_int_t ngpu,
    magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double *w,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsyevdx_2stage_m(
    magma_tally2_int_t ngpu,
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsyevdx_m(
    magma_tally2_int_t ngpu,
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsygst_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsygvd_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double *B, magma_tally2_int_t ldb,
    double *w,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsygvdx_2stage_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double *B, magma_tally2_int_t ldb,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsygvdx_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t itype, magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double *B, magma_tally2_int_t ldb,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

#ifdef REAL
magma_tally2_int_t
magma_tally2_dlaex0_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t n, double *d, double *e,
    double *Q, magma_tally2_int_t ldq,
    double *work, magma_tally2_int_t *iwork,
    magma_tally2_range_t range, double vl, double vu,
    magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dlaex1_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t n, double *d,
    double *Q, magma_tally2_int_t ldq,
    magma_tally2_int_t *indxq, double rho, magma_tally2_int_t cutpnt,
    double *work, magma_tally2_int_t *iwork,
    magma_tally2Double_ptr dwork[],
    magma_tally2_queue_t queues[Magma_tally2MaxGPUs][2],
    magma_tally2_range_t range, double vl, double vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dlaex3_m(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t k, magma_tally2_int_t n, magma_tally2_int_t n1, double *d,
    double *Q, magma_tally2_int_t ldq, double rho,
    double *dlamda, double *Q2, magma_tally2_int_t *indx,
    magma_tally2_int_t *ctot, double *w, double *s, magma_tally2_int_t *indxq,
    magma_tally2Double_ptr dwork[],
    magma_tally2_queue_t queues[Magma_tally2MaxGPUs][2],
    magma_tally2_range_t range, double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *info);
#endif

magma_tally2_int_t
magma_tally2_dlahr2_m(
    magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2_int_t nb,
    double *A, magma_tally2_int_t lda,
    double *tau,
    double *T, magma_tally2_int_t ldt,
    double *Y, magma_tally2_int_t ldy,
    struct dgehrd_data_tally2 *data );

magma_tally2_int_t
magma_tally2_dlahru_m(
    magma_tally2_int_t n, magma_tally2_int_t ihi, magma_tally2_int_t k, magma_tally2_int_t nb,
    double *A, magma_tally2_int_t lda,
    struct dgehrd_data_tally2 *data );

magma_tally2_int_t
magma_tally2_dpotrf_m(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dstedx_m(
    magma_tally2_int_t ngpu,
    magma_tally2_range_t range, magma_tally2_int_t n, double vl, double vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, double *d, double *e,
    double *Z, magma_tally2_int_t ldz,
    double *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dtrsm_m(
    magma_tally2_int_t ngpu,
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t transa, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n, double alpha,
    double *A, magma_tally2_int_t lda,
    double *B, magma_tally2_int_t ldb);

magma_tally2_int_t
magma_tally2_dorghr_m(
    magma_tally2_int_t n, magma_tally2_int_t ilo, magma_tally2_int_t ihi,
    double *A, magma_tally2_int_t lda,
    double *tau,
    double *T, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dorgqr_m(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    double *A, magma_tally2_int_t lda,
    double *tau,
    double *T, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dormqr_m(
    magma_tally2_int_t ngpu,
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    double *A,    magma_tally2_int_t lda,
    double *tau,
    double *C,    magma_tally2_int_t ldc,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dormtr_m(
    magma_tally2_int_t ngpu,
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double *A,    magma_tally2_int_t lda,
    double *tau,
    double *C,    magma_tally2_int_t ldc,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 function definitions / Data on GPU (alphabetical order)
*/
magma_tally2_int_t
magma_tally2_dgegqr_gpu(
    magma_tally2_int_t ikind, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dwork, double *work,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgelqf_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double *tau,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgels_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    double *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgels3_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    double *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqp3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt, double *tau,
    magma_tally2Double_ptr dwork, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqr2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dtau,
    magma_tally2Double_ptr        dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqr2x_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dtau,
    magma_tally2Double_ptr dT, magma_tally2Double_ptr ddA,
    magma_tally2Double_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqr2x2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dtau,
    magma_tally2Double_ptr dT, magma_tally2Double_ptr ddA,
    magma_tally2Double_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqr2x3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dtau,
    magma_tally2Double_ptr dT,
    magma_tally2Double_ptr ddA,
    magma_tally2Double_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqr2x4_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dtau,
    magma_tally2Double_ptr dT, magma_tally2Double_ptr ddA,
    magma_tally2Double_ptr dwork,
    magma_tally2_queue_t queue,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqrf_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double *tau,
    magma_tally2Double_ptr dT,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqrf2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double *tau,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqrf2_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dlA[], magma_tally2_int_t ldda,
    double *tau,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqrf3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double *tau,
    magma_tally2Double_ptr dT,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqrs_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double *tau,
    magma_tally2Double_ptr dT,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    double *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgeqrs3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double *tau,
    magma_tally2Double_ptr dT,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    double *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgerbt_gpu(
    magma_tally2_bool_t gen, magma_tally2_int_t n, magma_tally2_int_t nrhs, 
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda, 
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb, 
    double *U, double *V,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgessm_gpu(
    magma_tally2_order_t order, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2_int_t ib,
    magma_tally2_int_t *ipiv,
    magma_tally2Double_ptr dL1, magma_tally2_int_t lddl1,
    magma_tally2Double_ptr dL,  magma_tally2_int_t lddl,
    magma_tally2Double_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgesv_gpu(
    magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgesv_nopiv_gpu( 
    magma_tally2_int_t n, magma_tally2_int_t nrhs, 
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb, 
                 magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgetf2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgetrf_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgetrf_incpiv_gpu(
    magma_tally2_order_t order, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t ib,
    double    *hA, magma_tally2_int_t ldha,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double    *hL, magma_tally2_int_t ldhl,
    magma_tally2Double_ptr dL, magma_tally2_int_t lddl,
    magma_tally2_int_t *ipiv,
    magma_tally2Double_ptr dwork, magma_tally2_int_t lddwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgetrf_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr d_lA[], magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgetrf_nopiv_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgetrf2_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t offset,
    magma_tally2Double_ptr d_lAT[], magma_tally2_int_t lddat, magma_tally2_int_t *ipiv,
    magma_tally2Double_ptr d_lAP[],
    double *W, magma_tally2_int_t ldw,
    magma_tally2_queue_t queues[][2],
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgetri_gpu(
    magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2Double_ptr dwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgetrs_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dgetrs_nopiv_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsyevd_gpu(
    magma_tally2_vec_t jobz, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double *w,
    double *wA,  magma_tally2_int_t ldwa,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsyevdx_gpu(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double vl, double vu,
    magma_tally2_int_t il, magma_tally2_int_t iu,
    magma_tally2_int_t *m, double *w,
    double *wA,  magma_tally2_int_t ldwa,
    double *work, magma_tally2_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally2_int_t lrwork,
    #endif
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally2_int_t
magma_tally2_dsyevr_gpu(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double vl, double vu,
    magma_tally2_int_t il, magma_tally2_int_t iu, double abstol, magma_tally2_int_t *m,
    double *w,
    magma_tally2Double_ptr dZ, magma_tally2_int_t lddz,
    magma_tally2_int_t *isuppz,
    double *wA, magma_tally2_int_t ldwa,
    double *wZ, magma_tally2_int_t ldwz,
    double *work, magma_tally2_int_t lwork,
    double *rwork, magma_tally2_int_t lrwork,
    magma_tally2_int_t *iwork, magma_tally2_int_t liwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsyevx_gpu(
    magma_tally2_vec_t jobz, magma_tally2_range_t range, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double vl, double vu, magma_tally2_int_t il, magma_tally2_int_t iu,
    double abstol, magma_tally2_int_t *m,
    double *w,
    magma_tally2Double_ptr dZ, magma_tally2_int_t lddz,
    double *wA, magma_tally2_int_t ldwa,
    double *wZ, magma_tally2_int_t ldwz,
    double *work, magma_tally2_int_t lwork,
    double *rwork, magma_tally2_int_t *iwork,
    magma_tally2_int_t *ifail,
    magma_tally2_int_t *info);
#endif  // COMPLEX

magma_tally2_int_t
magma_tally2_dsygst_gpu(
    magma_tally2_int_t itype, magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsytrd_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double *d, double *e, double *tau,
    double *wA,  magma_tally2_int_t ldwa,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsytrd_sy2sb_mgpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    double *A, magma_tally2_int_t lda,
    double *tau,
    double *work, magma_tally2_int_t lwork,
    magma_tally2Double_ptr dAmgpu[], magma_tally2_int_t ldda,
    magma_tally2Double_ptr dTmgpu[], magma_tally2_int_t lddt,
    magma_tally2_int_t ngpu, magma_tally2_int_t distblk,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsytrd_sy2sb_mgpu_spec(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    double *A, magma_tally2_int_t lda,
    double *tau,
    double *work, magma_tally2_int_t lwork,
    magma_tally2Double_ptr dAmgpu[], magma_tally2_int_t ldda,
    magma_tally2Double_ptr dTmgpu[], magma_tally2_int_t lddt,
    magma_tally2_int_t ngpu, magma_tally2_int_t distblk,
    magma_tally2_queue_t queues[][20], magma_tally2_int_t nqueue,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsytrd_mgpu(
    magma_tally2_int_t ngpu, magma_tally2_int_t nqueue,
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    double *A, magma_tally2_int_t lda,
    double *d, double *e, double *tau,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsytrd2_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double *d, double *e, double *tau,
    double *wA,  magma_tally2_int_t ldwa,
    double *work, magma_tally2_int_t lwork,
    magma_tally2Double_ptr dwork, magma_tally2_int_t ldwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dsytrf_nopiv_tally2_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dlabrd_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb,
    double     *A, magma_tally2_int_t lda,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double *d, double *e, double *tauq, double *taup,
    double     *X, magma_tally2_int_t ldx,
    magma_tally2Double_ptr dX, magma_tally2_int_t lddx,
    double     *Y, magma_tally2_int_t ldy,
    magma_tally2Double_ptr dY, magma_tally2_int_t lddy );

magma_tally2_int_t
magma_tally2_dlaqps_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t offset,
    magma_tally2_int_t nb, magma_tally2_int_t *kb,
    magma_tally2Double_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt, double *tau,
    double *vn1, double *vn2,
    magma_tally2Double_ptr dauxv,
    magma_tally2Double_ptr dF, magma_tally2_int_t lddf);

magma_tally2_int_t
magma_tally2_dlaqps2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t offset,
    magma_tally2_int_t nb, magma_tally2_int_t *kb,
    magma_tally2Double_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt,
    magma_tally2Double_ptr dtau,
    magma_tally2Double_ptr dvn1, magma_tally2Double_ptr dvn2,
    magma_tally2Double_ptr dauxv,
    magma_tally2Double_ptr dF, magma_tally2_int_t lddf);

magma_tally2_int_t
magma_tally2_dlaqps3_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t offset,
    magma_tally2_int_t nb, magma_tally2_int_t *kb,
    magma_tally2Double_ptr dA,  magma_tally2_int_t ldda,
    magma_tally2_int_t *jpvt,
    magma_tally2Double_ptr dtau,
    magma_tally2Double_ptr dvn1, magma_tally2Double_ptr dvn2,
    magma_tally2Double_ptr dauxv,
    magma_tally2Double_ptr dF, magma_tally2_int_t lddf);

magma_tally2_int_t
magma_tally2_dlarf_gpu(
    magma_tally2_int_t m,  magma_tally2_int_t n,
    magma_tally2Double_const_ptr dv, magma_tally2Double_const_ptr dtau,
    magma_tally2Double_ptr dC,  magma_tally2_int_t lddc);

magma_tally2_int_t
magma_tally2_dlarfb_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans, magma_tally2_direct_t direct, magma_tally2_storev_t storev,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2Double_const_ptr dV, magma_tally2_int_t lddv,
    magma_tally2Double_const_ptr dT, magma_tally2_int_t lddt,
    magma_tally2Double_ptr dC,       magma_tally2_int_t lddc,
    magma_tally2Double_ptr dwork,    magma_tally2_int_t ldwork );

magma_tally2_int_t
magma_tally2_dlarfb_gpu_gemm(
    magma_tally2_side_t side, magma_tally2_trans_t trans, magma_tally2_direct_t direct, magma_tally2_storev_t storev,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2Double_const_ptr dV, magma_tally2_int_t lddv,
    magma_tally2Double_const_ptr dT, magma_tally2_int_t lddt,
    magma_tally2Double_ptr dC,       magma_tally2_int_t lddc,
    magma_tally2Double_ptr dwork,    magma_tally2_int_t ldwork,
    magma_tally2Double_ptr dworkvt,  magma_tally2_int_t ldworkvt);

magma_tally2_int_t
magma_tally2_dlarfb2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2Double_const_ptr dV, magma_tally2_int_t lddv,
    magma_tally2Double_const_ptr dT, magma_tally2_int_t lddt,
    magma_tally2Double_ptr dC,       magma_tally2_int_t lddc,
    magma_tally2Double_ptr dwork,    magma_tally2_int_t ldwork );

magma_tally2_int_t
magma_tally2_dlatrd_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n, magma_tally2_int_t nb, magma_tally2_int_t nb0,
    double *A,  magma_tally2_int_t lda,
    double *e, double *tau,
    double    *W,       magma_tally2_int_t ldw,
    magma_tally2Double_ptr dA[],    magma_tally2_int_t ldda, magma_tally2_int_t offset,
    magma_tally2Double_ptr dW[],    magma_tally2_int_t lddw,
    double    *hwork,   magma_tally2_int_t lhwork,
    magma_tally2Double_ptr dwork[], magma_tally2_int_t ldwork,
    magma_tally2_queue_t queues[] );

magma_tally2_int_t
magma_tally2_dlauum_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dposv_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dpotf2_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dpotrf_gpu(
    magma_tally2_uplo_t uplo,  magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dpotrf_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Double_ptr d_lA[], magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dpotrf_mgpu_right(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Double_ptr d_lA[], magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dpotrf3_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t off_i, magma_tally2_int_t off_j, magma_tally2_int_t nb,
    magma_tally2Double_ptr d_lA[],  magma_tally2_int_t ldda,
    magma_tally2Double_ptr d_lP[],  magma_tally2_int_t lddp,
    double *A, magma_tally2_int_t lda, magma_tally2_int_t h,
    magma_tally2_queue_t queues[][3], magma_tally2_event_t events[][5],
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dpotri_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dpotrs_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dssssm_gpu(
    magma_tally2_order_t order, magma_tally2_int_t m1, magma_tally2_int_t n1,
    magma_tally2_int_t m2, magma_tally2_int_t n2, magma_tally2_int_t k, magma_tally2_int_t ib,
    magma_tally2Double_ptr dA1, magma_tally2_int_t ldda1,
    magma_tally2Double_ptr dA2, magma_tally2_int_t ldda2,
    magma_tally2Double_ptr dL1, magma_tally2_int_t lddl1,
    magma_tally2Double_ptr dL2, magma_tally2_int_t lddl2,
    magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dtrtri_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_diag_t diag, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dtsqrt_gpu(
    magma_tally2_int_t *m, magma_tally2_int_t *n,
    double *A1, double *A2, magma_tally2_int_t *lda,
    double *tau,
    double *work, magma_tally2_int_t *lwork,
    magma_tally2Double_ptr dwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dtstrf_gpu(
    magma_tally2_order_t order, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t ib, magma_tally2_int_t nb,
    double    *hU, magma_tally2_int_t ldhu,
    magma_tally2Double_ptr dU, magma_tally2_int_t lddu,
    double    *hA, magma_tally2_int_t ldha,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double    *hL, magma_tally2_int_t ldhl,
    magma_tally2Double_ptr dL, magma_tally2_int_t lddl,
    magma_tally2_int_t *ipiv,
    double *hwork, magma_tally2_int_t ldhwork,
    magma_tally2Double_ptr dwork, magma_tally2_int_t lddwork,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dorgqr_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double *tau,
    magma_tally2Double_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dormql2_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double *tau,
    magma_tally2Double_ptr dC, magma_tally2_int_t lddc,
    double *wA, magma_tally2_int_t ldwa,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dormqr_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double *tau,
    magma_tally2Double_ptr dC, magma_tally2_int_t lddc,
    double *hwork, magma_tally2_int_t lwork,
    magma_tally2Double_ptr dT, magma_tally2_int_t nb,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dormqr2_gpu(
    magma_tally2_side_t side, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double *tau,
    magma_tally2Double_ptr dC, magma_tally2_int_t lddc,
    double    *wA, magma_tally2_int_t ldwa,
    magma_tally2_int_t *info);

magma_tally2_int_t
magma_tally2_dormtr_gpu(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double *tau,
    magma_tally2Double_ptr dC, magma_tally2_int_t lddc,
    double    *wA, magma_tally2_int_t ldwa,
    magma_tally2_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally2 utility function definitions
*/

extern const double MAGMA_tally2_D_NAN;
extern const double MAGMA_tally2_D_INF;

magma_tally2_int_t
magma_tally2_dnan_inf(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    const double *A, magma_tally2_int_t lda,
    magma_tally2_int_t *cnt_nan,
    magma_tally2_int_t *cnt_inf );

magma_tally2_int_t
magma_tally2_dnan_inf_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *cnt_nan,
    magma_tally2_int_t *cnt_inf );

void magma_tally2_dprint(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const double *A, magma_tally2_int_t lda );

void magma_tally2_dprint_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_const_ptr dA, magma_tally2_int_t ldda );

void dpanel_to_q_tally2(
    magma_tally2_uplo_t uplo, magma_tally2_int_t ib,
    double *A, magma_tally2_int_t lda,
    double *work );

void dq_to_panel_tally2(
    magma_tally2_uplo_t uplo, magma_tally2_int_t ib,
    double *A, magma_tally2_int_t lda,
    double *work );

#ifdef __cplusplus
}
#endif

#undef REAL

#endif /* MAGMA_tally2_D_H */
