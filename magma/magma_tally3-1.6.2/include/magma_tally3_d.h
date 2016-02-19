/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally3_z.h normal z -> d, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally3_D_H
#define MAGMA_tally3_D_H

#include "magma_tally3_types.h"
#include "magma_tally3_dgehrd_m.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 Auxiliary functions to get the NB used
*/
#ifdef REAL
magma_tally3_int_t magma_tally3_get_dlaex3_m_nb();       // defined in dlaex3_m.cpp
#endif

magma_tally3_int_t magma_tally3_get_dpotrf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_dgetrf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_dgetri_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_dgeqp3_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_dgeqrf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_dgeqlf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_dgehrd_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_dsytrd_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_dsytrf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_dsytrf_nopiv_tally3_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_dgelqf_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_dgebrd_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_dsygst_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_dgesvd_nb( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_dsygst_nb_m( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_get_dbulge_nb( magma_tally3_int_t m, magma_tally3_int_t nbthreads );
magma_tally3_int_t magma_tally3_get_dbulge_nb_mgpu( magma_tally3_int_t m );
magma_tally3_int_t magma_tally3_dbulge_get_Vblksiz( magma_tally3_int_t m, magma_tally3_int_t nb, magma_tally3_int_t nbthreads );
magma_tally3_int_t magma_tally3_get_dbulge_gcperf();


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 function definitions / Data on CPU (alphabetical order)
*/

#ifdef REAL
// only applicable to real [sd] precisions
void
magma_tally3_dmove_eig(
    magma_tally3_range_t range, magma_tally3_int_t n, double *w,
    magma_tally3_int_t *il, magma_tally3_int_t *iu, double vl, double vu, magma_tally3_int_t *m);

// defined in dlaex3.cpp
void
magma_tally3_dvrange(
    magma_tally3_int_t k, double *d, magma_tally3_int_t *il, magma_tally3_int_t *iu, double vl, double vu);

void
magma_tally3_dirange(
    magma_tally3_int_t k, magma_tally3_int_t *indxq, magma_tally3_int_t *iil, magma_tally3_int_t *iiu, magma_tally3_int_t il, magma_tally3_int_t iu);
#endif

magma_tally3_int_t
magma_tally3_dgebrd(
    magma_tally3_int_t m, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *d, double *e,
    double *tauq, double *taup,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeev(
    magma_tally3_vec_t jobvl, magma_tally3_vec_t jobvr, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    #ifdef COMPLEX
    double *w,
    #else
    double *wr, double *wi,
    #endif
    double *VL, magma_tally3_int_t ldvl,
    double *VR, magma_tally3_int_t ldvr,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgehrd(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    double *A, magma_tally3_int_t lda,
    double *tau,
    double *work, magma_tally3_int_t lwork,
    magma_tally3Double_ptr dT,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgehrd2(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    double *A, magma_tally3_int_t lda,
    double *tau,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgelqf(
    magma_tally3_int_t m, magma_tally3_int_t n,
    double *A,    magma_tally3_int_t lda,
    double *tau,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqlf(
    magma_tally3_int_t m, magma_tally3_int_t n,
    double *A,    magma_tally3_int_t lda,
    double *tau,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqp3(
    magma_tally3_int_t m, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    magma_tally3_int_t *jpvt, double *tau,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqrf(
    magma_tally3_int_t m, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *tau,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqrf_ooc(
    magma_tally3_int_t m, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *tau,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqrf4(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n,
    double *A,    magma_tally3_int_t lda,
    double *tau,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgesdd(
    magma_tally3_vec_t jobz, magma_tally3_int_t m, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *s,
    double *U, magma_tally3_int_t ldu,
    double *VT, magma_tally3_int_t ldvt,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *iwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgesv(
    magma_tally3_int_t n, magma_tally3_int_t nrhs,
    double *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    double *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgesv_rbt(
    magma_tally3_bool_t ref, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    double *A, magma_tally3_int_t lda, 
    double *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgesvd(
    magma_tally3_vec_t jobu, magma_tally3_vec_t jobvt, magma_tally3_int_t m, magma_tally3_int_t n,
    double *A,    magma_tally3_int_t lda, double *s,
    double *U,    magma_tally3_int_t ldu,
    double *VT,   magma_tally3_int_t ldvt,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgetf2_nopiv(
    magma_tally3_int_t m, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgetrf(
    magma_tally3_int_t m, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgetrf_nopiv(
    magma_tally3_int_t m, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgetrf_piv(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t NB,
    double *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgetrf2(
    magma_tally3_int_t m, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsyevd(
    magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *w,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsyevdx(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsyevdx_2stage(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally3_int_t
magma_tally3_dsyevr(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double vl, double vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, double abstol, magma_tally3_int_t *m,
    double *w,
    double *Z, magma_tally3_int_t ldz,
    magma_tally3_int_t *isuppz,
    double *work, magma_tally3_int_t lwork,
    double *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

// no real [sd] precisions available
magma_tally3_int_t
magma_tally3_dsyevx(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double vl, double vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, double abstol, magma_tally3_int_t *m,
    double *w,
    double *Z, magma_tally3_int_t ldz,
    double *work, magma_tally3_int_t lwork,
    double *rwork, magma_tally3_int_t *iwork,
    magma_tally3_int_t *ifail,
    magma_tally3_int_t *info);
#endif  // COMPLEX

magma_tally3_int_t
magma_tally3_dsygst(
    magma_tally3_int_t itype, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsygvd(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *B, magma_tally3_int_t ldb,
    double *w, double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsygvdx(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n, double *A, magma_tally3_int_t lda,
    double *B, magma_tally3_int_t ldb,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsygvdx_2stage(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *B, magma_tally3_int_t ldb,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsygvr(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *B, magma_tally3_int_t ldb,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    double abstol, magma_tally3_int_t *m, double *w,
    double *Z, magma_tally3_int_t ldz,
    magma_tally3_int_t *isuppz, double *work, magma_tally3_int_t lwork,
    double *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsygvx(
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n, double *A, magma_tally3_int_t lda,
    double *B, magma_tally3_int_t ldb,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    double abstol, magma_tally3_int_t *m, double *w,
    double *Z, magma_tally3_int_t ldz,
    double *work, magma_tally3_int_t lwork, double *rwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t *ifail,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsysv(magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
            double *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
            double *B, magma_tally3_int_t ldb,
            magma_tally3_int_t *info );

magma_tally3_int_t
magma_tally3_dsytrd(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *d, double *e, double *tau,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsytrf(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsytrf_nopiv_tally3(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsytrd_sb2st(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t Vblksiz,
    double *A, magma_tally3_int_t lda,
    double *d, double *e,
    double *V, magma_tally3_int_t ldv,
    double *TAU, magma_tally3_int_t compT,
    double *T, magma_tally3_int_t ldt);

magma_tally3_int_t
magma_tally3_dsytrd_sy2sb(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    double *A, magma_tally3_int_t lda,
    double *tau,
    double *work, magma_tally3_int_t lwork,
    magma_tally3Double_ptr dT,
    magma_tally3_int_t *info);

#ifdef REAL
// only applicable to real [sd] precisions
magma_tally3_int_t
magma_tally3_dlaex0(
    magma_tally3_int_t n, double *d, double *e,
    double *Q, magma_tally3_int_t ldq,
    double *work, magma_tally3_int_t *iwork,
    magma_tally3Double_ptr dwork,
    magma_tally3_range_t range, double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dlaex1(
    magma_tally3_int_t n, double *d,
    double *Q, magma_tally3_int_t ldq,
    magma_tally3_int_t *indxq, double rho, magma_tally3_int_t cutpnt,
    double *work, magma_tally3_int_t *iwork,
    magma_tally3Double_ptr dwork,
    magma_tally3_range_t range, double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dlaex3(
    magma_tally3_int_t k, magma_tally3_int_t n, magma_tally3_int_t n1, double *d,
    double *Q, magma_tally3_int_t ldq,
    double rho,
    double *dlamda, double *Q2, magma_tally3_int_t *indx,
    magma_tally3_int_t *ctot, double *w, double *s, magma_tally3_int_t *indxq,
    magma_tally3Double_ptr dwork,
    magma_tally3_range_t range, double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *info);
#endif  // REAL

magma_tally3_int_t
magma_tally3_dlasyf_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t *kb,
    double    *hA, magma_tally3_int_t lda,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3Double_ptr dW, magma_tally3_int_t lddw,
    magma_tally3_queue_t queues[], magma_tally3_event_t event[],
    magma_tally3_int_t *info);

magma_tally3_int_t
dsytrf_nopiv_tally3_cpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t ib,
    double *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsytrs_nopiv_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsysv_nopiv_gpu(
    magma_tally3_uplo_t uplo,  magma_tally3_int_t n, magma_tally3_int_t nrhs, 
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb, 
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dlahr2(
    magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3_int_t nb,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dV, magma_tally3_int_t lddv,
    double *A,  magma_tally3_int_t lda,
    double *tau,
    double *T,  magma_tally3_int_t ldt,
    double *Y,  magma_tally3_int_t ldy);

magma_tally3_int_t
magma_tally3_dlahru(
    magma_tally3_int_t n, magma_tally3_int_t ihi, magma_tally3_int_t k, magma_tally3_int_t nb,
    double     *A, magma_tally3_int_t lda,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dY, magma_tally3_int_t lddy,
    magma_tally3Double_ptr dV, magma_tally3_int_t lddv,
    magma_tally3Double_ptr dT,
    magma_tally3Double_ptr dwork);

#ifdef REAL
magma_tally3_int_t
magma_tally3_dlaln2(
    magma_tally3_int_t trans, magma_tally3_int_t na, magma_tally3_int_t nw,
    double smin, double ca, const double *A, magma_tally3_int_t lda,
    double d1, double d2,   const double *B, magma_tally3_int_t ldb,
    double wr, double wi, double *X, magma_tally3_int_t ldx,
    double *scale, double *xnorm,
    magma_tally3_int_t *info);
#endif

magma_tally3_int_t
magma_tally3_dlaqps(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    double *A,  magma_tally3_int_t lda,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt, double *tau, double *vn1, double *vn2,
    double *auxv,
    double *F,  magma_tally3_int_t ldf,
    magma_tally3Double_ptr dF, magma_tally3_int_t lddf );

#ifdef REAL
magma_tally3_int_t
magma_tally3_dlaqtrsd(
    magma_tally3_trans_t trans, magma_tally3_int_t n,
    const double *T, magma_tally3_int_t ldt,
    double *x,       magma_tally3_int_t ldx,
    const double *cnorm,
    magma_tally3_int_t *info);
#endif

magma_tally3_int_t
magma_tally3_dlatrd(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    double *A, magma_tally3_int_t lda,
    double *e, double *tau,
    double *W, magma_tally3_int_t ldw,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dW, magma_tally3_int_t lddw);

magma_tally3_int_t
magma_tally3_dlatrd2(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    double *A,  magma_tally3_int_t lda,
    double *e, double *tau,
    double *W,  magma_tally3_int_t ldw,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dW, magma_tally3_int_t lddw,
    magma_tally3Double_ptr dwork, magma_tally3_int_t ldwork);

#ifdef COMPLEX
magma_tally3_int_t
magma_tally3_dlatrsd(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_diag_t diag, magma_tally3_bool_t normin,
    magma_tally3_int_t n, const double *A, magma_tally3_int_t lda,
    double lambda,
    double *x,
    double *scale, double *cnorm,
    magma_tally3_int_t *info);
#endif

magma_tally3_int_t
magma_tally3_dlauum(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dposv(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    double *A, magma_tally3_int_t lda,
    double *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dpotrf(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dpotri(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dstedx(
    magma_tally3_range_t range, magma_tally3_int_t n, double vl, double vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, double *d, double *e,
    double *Z, magma_tally3_int_t ldz,
    double *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3Double_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dtrevc3(
    magma_tally3_side_t side, magma_tally3_vec_t howmany,
    magma_tally3_int_t *select, magma_tally3_int_t n,
    double *T,  magma_tally3_int_t ldt,
    double *VL, magma_tally3_int_t ldvl,
    double *VR, magma_tally3_int_t ldvr,
    magma_tally3_int_t mm, magma_tally3_int_t *mout,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dtrevc3_mt(
    magma_tally3_side_t side, magma_tally3_vec_t howmany,
    magma_tally3_int_t *select, magma_tally3_int_t n,
    double *T,  magma_tally3_int_t ldt,
    double *VL, magma_tally3_int_t ldvl,
    double *VR, magma_tally3_int_t ldvr,
    magma_tally3_int_t mm, magma_tally3_int_t *mout,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dtrtri(
    magma_tally3_uplo_t uplo, magma_tally3_diag_t diag, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dorghr(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    double *A, magma_tally3_int_t lda,
    double *tau,
    magma_tally3Double_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dorgqr(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    double *A, magma_tally3_int_t lda,
    double *tau,
    magma_tally3Double_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dorgqr2(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    double *A, magma_tally3_int_t lda,
    double *tau,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dormbr(
    magma_tally3_vect_t vect, magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    double *A, magma_tally3_int_t lda,
    double *tau,
    double *C, magma_tally3_int_t ldc,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dormlq(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    double *A, magma_tally3_int_t lda,
    double *tau,
    double *C, magma_tally3_int_t ldc,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

// not yet implemented
//magma_tally3_int_t magma_tally3_dunmrq( magma_tally3_side_t side, magma_tally3_trans_t trans,
//                          magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
//                          double *A, magma_tally3_int_t lda,
//                          double *tau,
//                          double *C, magma_tally3_int_t ldc,
//                          double *work, magma_tally3_int_t lwork,
//                          magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dormql(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    double *A, magma_tally3_int_t lda,
    double *tau,
    double *C, magma_tally3_int_t ldc,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dormqr(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    double *A, magma_tally3_int_t lda,
    double *tau,
    double *C, magma_tally3_int_t ldc,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dormtr(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n,
    double *A,    magma_tally3_int_t lda,
    double *tau,
    double *C,    magma_tally3_int_t ldc,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 function definitions / Data on CPU / Multi-GPU (alphabetical order)
*/
magma_tally3_int_t
magma_tally3_dgeev_m(
    magma_tally3_vec_t jobvl, magma_tally3_vec_t jobvr, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    #ifdef COMPLEX
    double *w,
    #else
    double *wr, double *wi,
    #endif
    double *VL, magma_tally3_int_t ldvl,
    double *VR, magma_tally3_int_t ldvr,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgehrd_m(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    double *A, magma_tally3_int_t lda,
    double *tau,
    double *work, magma_tally3_int_t lwork,
    double *T,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgetrf_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsyevd_m(
    magma_tally3_int_t ngpu,
    magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *w,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsyevdx_2stage_m(
    magma_tally3_int_t ngpu,
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsyevdx_m(
    magma_tally3_int_t ngpu,
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsygst_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *B, magma_tally3_int_t ldb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsygvd_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *B, magma_tally3_int_t ldb,
    double *w,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsygvdx_2stage_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *B, magma_tally3_int_t ldb,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsygvdx_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t itype, magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *B, magma_tally3_int_t ldb,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

#ifdef REAL
magma_tally3_int_t
magma_tally3_dlaex0_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t n, double *d, double *e,
    double *Q, magma_tally3_int_t ldq,
    double *work, magma_tally3_int_t *iwork,
    magma_tally3_range_t range, double vl, double vu,
    magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dlaex1_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t n, double *d,
    double *Q, magma_tally3_int_t ldq,
    magma_tally3_int_t *indxq, double rho, magma_tally3_int_t cutpnt,
    double *work, magma_tally3_int_t *iwork,
    magma_tally3Double_ptr dwork[],
    magma_tally3_queue_t queues[Magma_tally3MaxGPUs][2],
    magma_tally3_range_t range, double vl, double vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dlaex3_m(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t k, magma_tally3_int_t n, magma_tally3_int_t n1, double *d,
    double *Q, magma_tally3_int_t ldq, double rho,
    double *dlamda, double *Q2, magma_tally3_int_t *indx,
    magma_tally3_int_t *ctot, double *w, double *s, magma_tally3_int_t *indxq,
    magma_tally3Double_ptr dwork[],
    magma_tally3_queue_t queues[Magma_tally3MaxGPUs][2],
    magma_tally3_range_t range, double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *info);
#endif

magma_tally3_int_t
magma_tally3_dlahr2_m(
    magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3_int_t nb,
    double *A, magma_tally3_int_t lda,
    double *tau,
    double *T, magma_tally3_int_t ldt,
    double *Y, magma_tally3_int_t ldy,
    struct dgehrd_data_tally3 *data );

magma_tally3_int_t
magma_tally3_dlahru_m(
    magma_tally3_int_t n, magma_tally3_int_t ihi, magma_tally3_int_t k, magma_tally3_int_t nb,
    double *A, magma_tally3_int_t lda,
    struct dgehrd_data_tally3 *data );

magma_tally3_int_t
magma_tally3_dpotrf_m(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dstedx_m(
    magma_tally3_int_t ngpu,
    magma_tally3_range_t range, magma_tally3_int_t n, double vl, double vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, double *d, double *e,
    double *Z, magma_tally3_int_t ldz,
    double *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dtrsm_m(
    magma_tally3_int_t ngpu,
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transa, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n, double alpha,
    double *A, magma_tally3_int_t lda,
    double *B, magma_tally3_int_t ldb);

magma_tally3_int_t
magma_tally3_dorghr_m(
    magma_tally3_int_t n, magma_tally3_int_t ilo, magma_tally3_int_t ihi,
    double *A, magma_tally3_int_t lda,
    double *tau,
    double *T, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dorgqr_m(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    double *A, magma_tally3_int_t lda,
    double *tau,
    double *T, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dormqr_m(
    magma_tally3_int_t ngpu,
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    double *A,    magma_tally3_int_t lda,
    double *tau,
    double *C,    magma_tally3_int_t ldc,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dormtr_m(
    magma_tally3_int_t ngpu,
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n,
    double *A,    magma_tally3_int_t lda,
    double *tau,
    double *C,    magma_tally3_int_t ldc,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 function definitions / Data on GPU (alphabetical order)
*/
magma_tally3_int_t
magma_tally3_dgegqr_gpu(
    magma_tally3_int_t ikind, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dwork, double *work,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgelqf_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *tau,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgels_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb,
    double *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgels3_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb,
    double *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqp3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt, double *tau,
    magma_tally3Double_ptr dwork, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqr2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dtau,
    magma_tally3Double_ptr        dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqr2x_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dtau,
    magma_tally3Double_ptr dT, magma_tally3Double_ptr ddA,
    magma_tally3Double_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqr2x2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dtau,
    magma_tally3Double_ptr dT, magma_tally3Double_ptr ddA,
    magma_tally3Double_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqr2x3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dtau,
    magma_tally3Double_ptr dT,
    magma_tally3Double_ptr ddA,
    magma_tally3Double_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqr2x4_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dtau,
    magma_tally3Double_ptr dT, magma_tally3Double_ptr ddA,
    magma_tally3Double_ptr dwork,
    magma_tally3_queue_t queue,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqrf_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *tau,
    magma_tally3Double_ptr dT,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqrf2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *tau,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqrf2_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dlA[], magma_tally3_int_t ldda,
    double *tau,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqrf3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *tau,
    magma_tally3Double_ptr dT,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqrs_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *tau,
    magma_tally3Double_ptr dT,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb,
    double *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgeqrs3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *tau,
    magma_tally3Double_ptr dT,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb,
    double *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgerbt_gpu(
    magma_tally3_bool_t gen, magma_tally3_int_t n, magma_tally3_int_t nrhs, 
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda, 
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb, 
    double *U, double *V,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgessm_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3_int_t ib,
    magma_tally3_int_t *ipiv,
    magma_tally3Double_ptr dL1, magma_tally3_int_t lddl1,
    magma_tally3Double_ptr dL,  magma_tally3_int_t lddl,
    magma_tally3Double_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgesv_gpu(
    magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgesv_nopiv_gpu( 
    magma_tally3_int_t n, magma_tally3_int_t nrhs, 
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb, 
                 magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgetf2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgetrf_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgetrf_incpiv_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t ib,
    double    *hA, magma_tally3_int_t ldha,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double    *hL, magma_tally3_int_t ldhl,
    magma_tally3Double_ptr dL, magma_tally3_int_t lddl,
    magma_tally3_int_t *ipiv,
    magma_tally3Double_ptr dwork, magma_tally3_int_t lddwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgetrf_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr d_lA[], magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgetrf_nopiv_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgetrf2_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t offset,
    magma_tally3Double_ptr d_lAT[], magma_tally3_int_t lddat, magma_tally3_int_t *ipiv,
    magma_tally3Double_ptr d_lAP[],
    double *W, magma_tally3_int_t ldw,
    magma_tally3_queue_t queues[][2],
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgetri_gpu(
    magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3Double_ptr dwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgetrs_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t *ipiv,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dgetrs_nopiv_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsyevd_gpu(
    magma_tally3_vec_t jobz, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *w,
    double *wA,  magma_tally3_int_t ldwa,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsyevdx_gpu(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo,
    magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double vl, double vu,
    magma_tally3_int_t il, magma_tally3_int_t iu,
    magma_tally3_int_t *m, double *w,
    double *wA,  magma_tally3_int_t ldwa,
    double *work, magma_tally3_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_tally3_int_t lrwork,
    #endif
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_tally3_int_t
magma_tally3_dsyevr_gpu(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double vl, double vu,
    magma_tally3_int_t il, magma_tally3_int_t iu, double abstol, magma_tally3_int_t *m,
    double *w,
    magma_tally3Double_ptr dZ, magma_tally3_int_t lddz,
    magma_tally3_int_t *isuppz,
    double *wA, magma_tally3_int_t ldwa,
    double *wZ, magma_tally3_int_t ldwz,
    double *work, magma_tally3_int_t lwork,
    double *rwork, magma_tally3_int_t lrwork,
    magma_tally3_int_t *iwork, magma_tally3_int_t liwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsyevx_gpu(
    magma_tally3_vec_t jobz, magma_tally3_range_t range, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double vl, double vu, magma_tally3_int_t il, magma_tally3_int_t iu,
    double abstol, magma_tally3_int_t *m,
    double *w,
    magma_tally3Double_ptr dZ, magma_tally3_int_t lddz,
    double *wA, magma_tally3_int_t ldwa,
    double *wZ, magma_tally3_int_t ldwz,
    double *work, magma_tally3_int_t lwork,
    double *rwork, magma_tally3_int_t *iwork,
    magma_tally3_int_t *ifail,
    magma_tally3_int_t *info);
#endif  // COMPLEX

magma_tally3_int_t
magma_tally3_dsygst_gpu(
    magma_tally3_int_t itype, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsytrd_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *d, double *e, double *tau,
    double *wA,  magma_tally3_int_t ldwa,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsytrd_sy2sb_mgpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    double *A, magma_tally3_int_t lda,
    double *tau,
    double *work, magma_tally3_int_t lwork,
    magma_tally3Double_ptr dAmgpu[], magma_tally3_int_t ldda,
    magma_tally3Double_ptr dTmgpu[], magma_tally3_int_t lddt,
    magma_tally3_int_t ngpu, magma_tally3_int_t distblk,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsytrd_sy2sb_mgpu_spec(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    double *A, magma_tally3_int_t lda,
    double *tau,
    double *work, magma_tally3_int_t lwork,
    magma_tally3Double_ptr dAmgpu[], magma_tally3_int_t ldda,
    magma_tally3Double_ptr dTmgpu[], magma_tally3_int_t lddt,
    magma_tally3_int_t ngpu, magma_tally3_int_t distblk,
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsytrd_mgpu(
    magma_tally3_int_t ngpu, magma_tally3_int_t nqueue,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    double *A, magma_tally3_int_t lda,
    double *d, double *e, double *tau,
    double *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsytrd2_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *d, double *e, double *tau,
    double *wA,  magma_tally3_int_t ldwa,
    double *work, magma_tally3_int_t lwork,
    magma_tally3Double_ptr dwork, magma_tally3_int_t ldwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dsytrf_nopiv_tally3_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dlabrd_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb,
    double     *A, magma_tally3_int_t lda,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *d, double *e, double *tauq, double *taup,
    double     *X, magma_tally3_int_t ldx,
    magma_tally3Double_ptr dX, magma_tally3_int_t lddx,
    double     *Y, magma_tally3_int_t ldy,
    magma_tally3Double_ptr dY, magma_tally3_int_t lddy );

magma_tally3_int_t
magma_tally3_dlaqps_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3Double_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt, double *tau,
    double *vn1, double *vn2,
    magma_tally3Double_ptr dauxv,
    magma_tally3Double_ptr dF, magma_tally3_int_t lddf);

magma_tally3_int_t
magma_tally3_dlaqps2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3Double_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt,
    magma_tally3Double_ptr dtau,
    magma_tally3Double_ptr dvn1, magma_tally3Double_ptr dvn2,
    magma_tally3Double_ptr dauxv,
    magma_tally3Double_ptr dF, magma_tally3_int_t lddf);

magma_tally3_int_t
magma_tally3_dlaqps3_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3Double_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt,
    magma_tally3Double_ptr dtau,
    magma_tally3Double_ptr dvn1, magma_tally3Double_ptr dvn2,
    magma_tally3Double_ptr dauxv,
    magma_tally3Double_ptr dF, magma_tally3_int_t lddf);

magma_tally3_int_t
magma_tally3_dlarf_gpu(
    magma_tally3_int_t m,  magma_tally3_int_t n,
    magma_tally3Double_const_ptr dv, magma_tally3Double_const_ptr dtau,
    magma_tally3Double_ptr dC,  magma_tally3_int_t lddc);

magma_tally3_int_t
magma_tally3_dlarfb_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans, magma_tally3_direct_t direct, magma_tally3_storev_t storev,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Double_const_ptr dV, magma_tally3_int_t lddv,
    magma_tally3Double_const_ptr dT, magma_tally3_int_t lddt,
    magma_tally3Double_ptr dC,       magma_tally3_int_t lddc,
    magma_tally3Double_ptr dwork,    magma_tally3_int_t ldwork );

magma_tally3_int_t
magma_tally3_dlarfb_gpu_gemm(
    magma_tally3_side_t side, magma_tally3_trans_t trans, magma_tally3_direct_t direct, magma_tally3_storev_t storev,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Double_const_ptr dV, magma_tally3_int_t lddv,
    magma_tally3Double_const_ptr dT, magma_tally3_int_t lddt,
    magma_tally3Double_ptr dC,       magma_tally3_int_t lddc,
    magma_tally3Double_ptr dwork,    magma_tally3_int_t ldwork,
    magma_tally3Double_ptr dworkvt,  magma_tally3_int_t ldworkvt);

magma_tally3_int_t
magma_tally3_dlarfb2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Double_const_ptr dV, magma_tally3_int_t lddv,
    magma_tally3Double_const_ptr dT, magma_tally3_int_t lddt,
    magma_tally3Double_ptr dC,       magma_tally3_int_t lddc,
    magma_tally3Double_ptr dwork,    magma_tally3_int_t ldwork );

magma_tally3_int_t
magma_tally3_dlatrd_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t nb0,
    double *A,  magma_tally3_int_t lda,
    double *e, double *tau,
    double    *W,       magma_tally3_int_t ldw,
    magma_tally3Double_ptr dA[],    magma_tally3_int_t ldda, magma_tally3_int_t offset,
    magma_tally3Double_ptr dW[],    magma_tally3_int_t lddw,
    double    *hwork,   magma_tally3_int_t lhwork,
    magma_tally3Double_ptr dwork[], magma_tally3_int_t ldwork,
    magma_tally3_queue_t queues[] );

magma_tally3_int_t
magma_tally3_dlauum_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dposv_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dpotf2_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dpotrf_gpu(
    magma_tally3_uplo_t uplo,  magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dpotrf_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Double_ptr d_lA[], magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dpotrf_mgpu_right(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Double_ptr d_lA[], magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dpotrf3_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t off_i, magma_tally3_int_t off_j, magma_tally3_int_t nb,
    magma_tally3Double_ptr d_lA[],  magma_tally3_int_t ldda,
    magma_tally3Double_ptr d_lP[],  magma_tally3_int_t lddp,
    double *A, magma_tally3_int_t lda, magma_tally3_int_t h,
    magma_tally3_queue_t queues[][3], magma_tally3_event_t events[][5],
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dpotri_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dpotrs_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dssssm_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m1, magma_tally3_int_t n1,
    magma_tally3_int_t m2, magma_tally3_int_t n2, magma_tally3_int_t k, magma_tally3_int_t ib,
    magma_tally3Double_ptr dA1, magma_tally3_int_t ldda1,
    magma_tally3Double_ptr dA2, magma_tally3_int_t ldda2,
    magma_tally3Double_ptr dL1, magma_tally3_int_t lddl1,
    magma_tally3Double_ptr dL2, magma_tally3_int_t lddl2,
    magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dtrtri_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_diag_t diag, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dtsqrt_gpu(
    magma_tally3_int_t *m, magma_tally3_int_t *n,
    double *A1, double *A2, magma_tally3_int_t *lda,
    double *tau,
    double *work, magma_tally3_int_t *lwork,
    magma_tally3Double_ptr dwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dtstrf_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t ib, magma_tally3_int_t nb,
    double    *hU, magma_tally3_int_t ldhu,
    magma_tally3Double_ptr dU, magma_tally3_int_t lddu,
    double    *hA, magma_tally3_int_t ldha,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double    *hL, magma_tally3_int_t ldhl,
    magma_tally3Double_ptr dL, magma_tally3_int_t lddl,
    magma_tally3_int_t *ipiv,
    double *hwork, magma_tally3_int_t ldhwork,
    magma_tally3Double_ptr dwork, magma_tally3_int_t lddwork,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dorgqr_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *tau,
    magma_tally3Double_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dormql2_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *tau,
    magma_tally3Double_ptr dC, magma_tally3_int_t lddc,
    double *wA, magma_tally3_int_t ldwa,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dormqr_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *tau,
    magma_tally3Double_ptr dC, magma_tally3_int_t lddc,
    double *hwork, magma_tally3_int_t lwork,
    magma_tally3Double_ptr dT, magma_tally3_int_t nb,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dormqr2_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *tau,
    magma_tally3Double_ptr dC, magma_tally3_int_t lddc,
    double    *wA, magma_tally3_int_t ldwa,
    magma_tally3_int_t *info);

magma_tally3_int_t
magma_tally3_dormtr_gpu(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *tau,
    magma_tally3Double_ptr dC, magma_tally3_int_t lddc,
    double    *wA, magma_tally3_int_t ldwa,
    magma_tally3_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_tally3 utility function definitions
*/

extern const double MAGMA_tally3_D_NAN;
extern const double MAGMA_tally3_D_INF;

magma_tally3_int_t
magma_tally3_dnan_inf(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    const double *A, magma_tally3_int_t lda,
    magma_tally3_int_t *cnt_nan,
    magma_tally3_int_t *cnt_inf );

magma_tally3_int_t
magma_tally3_dnan_inf_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_const_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *cnt_nan,
    magma_tally3_int_t *cnt_inf );

void magma_tally3_dprint(
    magma_tally3_int_t m, magma_tally3_int_t n,
    const double *A, magma_tally3_int_t lda );

void magma_tally3_dprint_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_const_ptr dA, magma_tally3_int_t ldda );

void dpanel_to_q_tally3(
    magma_tally3_uplo_t uplo, magma_tally3_int_t ib,
    double *A, magma_tally3_int_t lda,
    double *work );

void dq_to_panel_tally3(
    magma_tally3_uplo_t uplo, magma_tally3_int_t ib,
    double *A, magma_tally3_int_t lda,
    double *work );

#ifdef __cplusplus
}
#endif

#undef REAL

#endif /* MAGMA_tally3_D_H */
