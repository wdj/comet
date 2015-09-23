/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_minproduct_z.h normal z -> d, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_minproduct_D_H
#define MAGMA_minproduct_D_H

#include "magma_minproduct_types.h"
#include "magma_minproduct_dgehrd_m.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct Auxiliary functions to get the NB used
*/
#ifdef REAL
magma_minproduct_int_t magma_minproduct_get_dlaex3_m_nb();       // defined in dlaex3_m.cpp
#endif

magma_minproduct_int_t magma_minproduct_get_dpotrf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_dgetrf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_dgetri_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_dgeqp3_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_dgeqrf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_dgeqlf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_dgehrd_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_dsytrd_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_dsytrf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_dsytrf_nopiv_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_dgelqf_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_dgebrd_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_dsygst_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_dgesvd_nb( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_dsygst_nb_m( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_get_dbulge_nb( magma_minproduct_int_t m, magma_minproduct_int_t nbthreads );
magma_minproduct_int_t magma_minproduct_get_dbulge_nb_mgpu( magma_minproduct_int_t m );
magma_minproduct_int_t magma_minproduct_dbulge_get_Vblksiz( magma_minproduct_int_t m, magma_minproduct_int_t nb, magma_minproduct_int_t nbthreads );
magma_minproduct_int_t magma_minproduct_get_dbulge_gcperf();


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct function definitions / Data on CPU (alphabetical order)
*/

#ifdef REAL
// only applicable to real [sd] precisions
void
magma_minproduct_dmove_eig(
    magma_minproduct_range_t range, magma_minproduct_int_t n, double *w,
    magma_minproduct_int_t *il, magma_minproduct_int_t *iu, double vl, double vu, magma_minproduct_int_t *m);

// defined in dlaex3.cpp
void
magma_minproduct_dvrange(
    magma_minproduct_int_t k, double *d, magma_minproduct_int_t *il, magma_minproduct_int_t *iu, double vl, double vu);

void
magma_minproduct_dirange(
    magma_minproduct_int_t k, magma_minproduct_int_t *indxq, magma_minproduct_int_t *iil, magma_minproduct_int_t *iiu, magma_minproduct_int_t il, magma_minproduct_int_t iu);
#endif

magma_minproduct_int_t
magma_minproduct_dgebrd(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double *d, double *e,
    double *tauq, double *taup,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeev(
    magma_minproduct_vec_t jobvl, magma_minproduct_vec_t jobvr, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    #ifdef COMPLEX
    double *w,
    #else
    double *wr, double *wi,
    #endif
    double *VL, magma_minproduct_int_t ldvl,
    double *VR, magma_minproduct_int_t ldvr,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgehrd(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    double *work, magma_minproduct_int_t lwork,
    magma_minproductDouble_ptr dT,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgehrd2(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgelqf(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double *A,    magma_minproduct_int_t lda,
    double *tau,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqlf(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double *A,    magma_minproduct_int_t lda,
    double *tau,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqp3(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *jpvt, double *tau,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqrf(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqrf_ooc(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqrf4(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double *A,    magma_minproduct_int_t lda,
    double *tau,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgesdd(
    magma_minproduct_vec_t jobz, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double *s,
    double *U, magma_minproduct_int_t ldu,
    double *VT, magma_minproduct_int_t ldvt,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_minproduct_int_t *iwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgesv(
    magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    double *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    double *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgesv_rbt(
    magma_minproduct_bool_t ref, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    double *A, magma_minproduct_int_t lda, 
    double *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgesvd(
    magma_minproduct_vec_t jobu, magma_minproduct_vec_t jobvt, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double *A,    magma_minproduct_int_t lda, double *s,
    double *U,    magma_minproduct_int_t ldu,
    double *VT,   magma_minproduct_int_t ldvt,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgetf2_nopiv(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgetrf(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgetrf_nopiv(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgetrf_piv(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t NB,
    double *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgetrf2(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsyevd(
    magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double *w,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsyevdx(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsyevdx_2stage(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_minproduct_int_t
magma_minproduct_dsyevr(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double vl, double vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, double abstol, magma_minproduct_int_t *m,
    double *w,
    double *Z, magma_minproduct_int_t ldz,
    magma_minproduct_int_t *isuppz,
    double *work, magma_minproduct_int_t lwork,
    double *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

// no real [sd] precisions available
magma_minproduct_int_t
magma_minproduct_dsyevx(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double vl, double vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, double abstol, magma_minproduct_int_t *m,
    double *w,
    double *Z, magma_minproduct_int_t ldz,
    double *work, magma_minproduct_int_t lwork,
    double *rwork, magma_minproduct_int_t *iwork,
    magma_minproduct_int_t *ifail,
    magma_minproduct_int_t *info);
#endif  // COMPLEX

magma_minproduct_int_t
magma_minproduct_dsygst(
    magma_minproduct_int_t itype, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsygvd(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double *B, magma_minproduct_int_t ldb,
    double *w, double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsygvdx(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n, double *A, magma_minproduct_int_t lda,
    double *B, magma_minproduct_int_t ldb,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsygvdx_2stage(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double *B, magma_minproduct_int_t ldb,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsygvr(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double *B, magma_minproduct_int_t ldb,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    double abstol, magma_minproduct_int_t *m, double *w,
    double *Z, magma_minproduct_int_t ldz,
    magma_minproduct_int_t *isuppz, double *work, magma_minproduct_int_t lwork,
    double *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsygvx(
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n, double *A, magma_minproduct_int_t lda,
    double *B, magma_minproduct_int_t ldb,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    double abstol, magma_minproduct_int_t *m, double *w,
    double *Z, magma_minproduct_int_t ldz,
    double *work, magma_minproduct_int_t lwork, double *rwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t *ifail,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsysv(magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
            double *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
            double *B, magma_minproduct_int_t ldb,
            magma_minproduct_int_t *info );

magma_minproduct_int_t
magma_minproduct_dsytrd(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double *d, double *e, double *tau,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsytrf(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsytrf_nopiv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsytrd_sb2st(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t Vblksiz,
    double *A, magma_minproduct_int_t lda,
    double *d, double *e,
    double *V, magma_minproduct_int_t ldv,
    double *TAU, magma_minproduct_int_t compT,
    double *T, magma_minproduct_int_t ldt);

magma_minproduct_int_t
magma_minproduct_dsytrd_sy2sb(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    double *work, magma_minproduct_int_t lwork,
    magma_minproductDouble_ptr dT,
    magma_minproduct_int_t *info);

#ifdef REAL
// only applicable to real [sd] precisions
magma_minproduct_int_t
magma_minproduct_dlaex0(
    magma_minproduct_int_t n, double *d, double *e,
    double *Q, magma_minproduct_int_t ldq,
    double *work, magma_minproduct_int_t *iwork,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_range_t range, double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dlaex1(
    magma_minproduct_int_t n, double *d,
    double *Q, magma_minproduct_int_t ldq,
    magma_minproduct_int_t *indxq, double rho, magma_minproduct_int_t cutpnt,
    double *work, magma_minproduct_int_t *iwork,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_range_t range, double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dlaex3(
    magma_minproduct_int_t k, magma_minproduct_int_t n, magma_minproduct_int_t n1, double *d,
    double *Q, magma_minproduct_int_t ldq,
    double rho,
    double *dlamda, double *Q2, magma_minproduct_int_t *indx,
    magma_minproduct_int_t *ctot, double *w, double *s, magma_minproduct_int_t *indxq,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_range_t range, double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *info);
#endif  // REAL

magma_minproduct_int_t
magma_minproduct_dlasyf_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    double    *hA, magma_minproduct_int_t lda,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductDouble_ptr dW, magma_minproduct_int_t lddw,
    magma_minproduct_queue_t queues[], magma_minproduct_event_t event[],
    magma_minproduct_int_t *info);

magma_minproduct_int_t
dsytrf_nopiv_cpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t ib,
    double *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsytrs_nopiv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsysv_nopiv_gpu(
    magma_minproduct_uplo_t uplo,  magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb, 
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dlahr2(
    magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dV, magma_minproduct_int_t lddv,
    double *A,  magma_minproduct_int_t lda,
    double *tau,
    double *T,  magma_minproduct_int_t ldt,
    double *Y,  magma_minproduct_int_t ldy);

magma_minproduct_int_t
magma_minproduct_dlahru(
    magma_minproduct_int_t n, magma_minproduct_int_t ihi, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    double     *A, magma_minproduct_int_t lda,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dY, magma_minproduct_int_t lddy,
    magma_minproductDouble_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductDouble_ptr dT,
    magma_minproductDouble_ptr dwork);

#ifdef REAL
magma_minproduct_int_t
magma_minproduct_dlaln2(
    magma_minproduct_int_t trans, magma_minproduct_int_t na, magma_minproduct_int_t nw,
    double smin, double ca, const double *A, magma_minproduct_int_t lda,
    double d1, double d2,   const double *B, magma_minproduct_int_t ldb,
    double wr, double wi, double *X, magma_minproduct_int_t ldx,
    double *scale, double *xnorm,
    magma_minproduct_int_t *info);
#endif

magma_minproduct_int_t
magma_minproduct_dlaqps(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t offset,
    magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    double *A,  magma_minproduct_int_t lda,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt, double *tau, double *vn1, double *vn2,
    double *auxv,
    double *F,  magma_minproduct_int_t ldf,
    magma_minproductDouble_ptr dF, magma_minproduct_int_t lddf );

#ifdef REAL
magma_minproduct_int_t
magma_minproduct_dlaqtrsd(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n,
    const double *T, magma_minproduct_int_t ldt,
    double *x,       magma_minproduct_int_t ldx,
    const double *cnorm,
    magma_minproduct_int_t *info);
#endif

magma_minproduct_int_t
magma_minproduct_dlatrd(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    double *A, magma_minproduct_int_t lda,
    double *e, double *tau,
    double *W, magma_minproduct_int_t ldw,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dW, magma_minproduct_int_t lddw);

magma_minproduct_int_t
magma_minproduct_dlatrd2(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    double *A,  magma_minproduct_int_t lda,
    double *e, double *tau,
    double *W,  magma_minproduct_int_t ldw,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dW, magma_minproduct_int_t lddw,
    magma_minproductDouble_ptr dwork, magma_minproduct_int_t ldwork);

#ifdef COMPLEX
magma_minproduct_int_t
magma_minproduct_dlatrsd(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_diag_t diag, magma_minproduct_bool_t normin,
    magma_minproduct_int_t n, const double *A, magma_minproduct_int_t lda,
    double lambda,
    double *x,
    double *scale, double *cnorm,
    magma_minproduct_int_t *info);
#endif

magma_minproduct_int_t
magma_minproduct_dlauum(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dposv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    double *A, magma_minproduct_int_t lda,
    double *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dpotrf(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dpotri(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dstedx(
    magma_minproduct_range_t range, magma_minproduct_int_t n, double vl, double vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, double *d, double *e,
    double *Z, magma_minproduct_int_t ldz,
    double *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dtrevc3(
    magma_minproduct_side_t side, magma_minproduct_vec_t howmany,
    magma_minproduct_int_t *select, magma_minproduct_int_t n,
    double *T,  magma_minproduct_int_t ldt,
    double *VL, magma_minproduct_int_t ldvl,
    double *VR, magma_minproduct_int_t ldvr,
    magma_minproduct_int_t mm, magma_minproduct_int_t *mout,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dtrevc3_mt(
    magma_minproduct_side_t side, magma_minproduct_vec_t howmany,
    magma_minproduct_int_t *select, magma_minproduct_int_t n,
    double *T,  magma_minproduct_int_t ldt,
    double *VL, magma_minproduct_int_t ldvl,
    double *VR, magma_minproduct_int_t ldvr,
    magma_minproduct_int_t mm, magma_minproduct_int_t *mout,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dtrtri(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dorghr(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    magma_minproductDouble_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dorgqr(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    magma_minproductDouble_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dorgqr2(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dormbr(
    magma_minproduct_vect_t vect, magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    double *C, magma_minproduct_int_t ldc,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dormlq(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    double *C, magma_minproduct_int_t ldc,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

// not yet implemented
//magma_minproduct_int_t magma_minproduct_dunmrq( magma_minproduct_side_t side, magma_minproduct_trans_t trans,
//                          magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
//                          double *A, magma_minproduct_int_t lda,
//                          double *tau,
//                          double *C, magma_minproduct_int_t ldc,
//                          double *work, magma_minproduct_int_t lwork,
//                          magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dormql(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    double *C, magma_minproduct_int_t ldc,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dormqr(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    double *C, magma_minproduct_int_t ldc,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dormtr(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double *A,    magma_minproduct_int_t lda,
    double *tau,
    double *C,    magma_minproduct_int_t ldc,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct function definitions / Data on CPU / Multi-GPU (alphabetical order)
*/
magma_minproduct_int_t
magma_minproduct_dgeev_m(
    magma_minproduct_vec_t jobvl, magma_minproduct_vec_t jobvr, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    #ifdef COMPLEX
    double *w,
    #else
    double *wr, double *wi,
    #endif
    double *VL, magma_minproduct_int_t ldvl,
    double *VR, magma_minproduct_int_t ldvr,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgehrd_m(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    double *work, magma_minproduct_int_t lwork,
    double *T,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgetrf_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsyevd_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double *w,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsyevdx_2stage_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsyevdx_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsygst_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t itype, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsygvd_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double *B, magma_minproduct_int_t ldb,
    double *w,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsygvdx_2stage_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double *B, magma_minproduct_int_t ldb,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsygvdx_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t itype, magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double *B, magma_minproduct_int_t ldb,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

#ifdef REAL
magma_minproduct_int_t
magma_minproduct_dlaex0_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t n, double *d, double *e,
    double *Q, magma_minproduct_int_t ldq,
    double *work, magma_minproduct_int_t *iwork,
    magma_minproduct_range_t range, double vl, double vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dlaex1_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t n, double *d,
    double *Q, magma_minproduct_int_t ldq,
    magma_minproduct_int_t *indxq, double rho, magma_minproduct_int_t cutpnt,
    double *work, magma_minproduct_int_t *iwork,
    magma_minproductDouble_ptr dwork[],
    magma_minproduct_queue_t queues[Magma_minproductMaxGPUs][2],
    magma_minproduct_range_t range, double vl, double vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dlaex3_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t k, magma_minproduct_int_t n, magma_minproduct_int_t n1, double *d,
    double *Q, magma_minproduct_int_t ldq, double rho,
    double *dlamda, double *Q2, magma_minproduct_int_t *indx,
    magma_minproduct_int_t *ctot, double *w, double *s, magma_minproduct_int_t *indxq,
    magma_minproductDouble_ptr dwork[],
    magma_minproduct_queue_t queues[Magma_minproductMaxGPUs][2],
    magma_minproduct_range_t range, double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *info);
#endif

magma_minproduct_int_t
magma_minproduct_dlahr2_m(
    magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    double *T, magma_minproduct_int_t ldt,
    double *Y, magma_minproduct_int_t ldy,
    struct dgehrd_data *data );

magma_minproduct_int_t
magma_minproduct_dlahru_m(
    magma_minproduct_int_t n, magma_minproduct_int_t ihi, magma_minproduct_int_t k, magma_minproduct_int_t nb,
    double *A, magma_minproduct_int_t lda,
    struct dgehrd_data *data );

magma_minproduct_int_t
magma_minproduct_dpotrf_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dstedx_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_range_t range, magma_minproduct_int_t n, double vl, double vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, double *d, double *e,
    double *Z, magma_minproduct_int_t ldz,
    double *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dtrsm_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transa, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n, double alpha,
    double *A, magma_minproduct_int_t lda,
    double *B, magma_minproduct_int_t ldb);

magma_minproduct_int_t
magma_minproduct_dorghr_m(
    magma_minproduct_int_t n, magma_minproduct_int_t ilo, magma_minproduct_int_t ihi,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    double *T, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dorgqr_m(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    double *T, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dormqr_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double *A,    magma_minproduct_int_t lda,
    double *tau,
    double *C,    magma_minproduct_int_t ldc,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dormtr_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double *A,    magma_minproduct_int_t lda,
    double *tau,
    double *C,    magma_minproduct_int_t ldc,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct function definitions / Data on GPU (alphabetical order)
*/
magma_minproduct_int_t
magma_minproduct_dgegqr_gpu(
    magma_minproduct_int_t ikind, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dwork, double *work,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgelqf_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *tau,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgels_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    double *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgels3_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    double *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqp3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt, double *tau,
    magma_minproductDouble_ptr dwork, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork,
    #endif
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqr2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dtau,
    magma_minproductDouble_ptr        dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqr2x_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dtau,
    magma_minproductDouble_ptr dT, magma_minproductDouble_ptr ddA,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqr2x2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dtau,
    magma_minproductDouble_ptr dT, magma_minproductDouble_ptr ddA,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqr2x3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dtau,
    magma_minproductDouble_ptr dT,
    magma_minproductDouble_ptr ddA,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqr2x4_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dtau,
    magma_minproductDouble_ptr dT, magma_minproductDouble_ptr ddA,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_queue_t queue,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqrf_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *tau,
    magma_minproductDouble_ptr dT,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqrf2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *tau,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqrf2_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dlA[], magma_minproduct_int_t ldda,
    double *tau,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqrf3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *tau,
    magma_minproductDouble_ptr dT,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqrs_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *tau,
    magma_minproductDouble_ptr dT,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    double *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgeqrs3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *tau,
    magma_minproductDouble_ptr dT,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    double *hwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgerbt_gpu(
    magma_minproduct_bool_t gen, magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda, 
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb, 
    double *U, double *V,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgessm_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t ib,
    magma_minproduct_int_t *ipiv,
    magma_minproductDouble_ptr dL1, magma_minproduct_int_t lddl1,
    magma_minproductDouble_ptr dL,  magma_minproduct_int_t lddl,
    magma_minproductDouble_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgesv_gpu(
    magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgesv_nopiv_gpu( 
    magma_minproduct_int_t n, magma_minproduct_int_t nrhs, 
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb, 
                 magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgetf2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgetrf_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgetrf_incpiv_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t ib,
    double    *hA, magma_minproduct_int_t ldha,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double    *hL, magma_minproduct_int_t ldhl,
    magma_minproductDouble_ptr dL, magma_minproduct_int_t lddl,
    magma_minproduct_int_t *ipiv,
    magma_minproductDouble_ptr dwork, magma_minproduct_int_t lddwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgetrf_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr d_lA[], magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgetrf_nopiv_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgetrf2_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t offset,
    magma_minproductDouble_ptr d_lAT[], magma_minproduct_int_t lddat, magma_minproduct_int_t *ipiv,
    magma_minproductDouble_ptr d_lAP[],
    double *W, magma_minproduct_int_t ldw,
    magma_minproduct_queue_t queues[][2],
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgetri_gpu(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductDouble_ptr dwork, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgetrs_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dgetrs_nopiv_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsyevd_gpu(
    magma_minproduct_vec_t jobz, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *w,
    double *wA,  magma_minproduct_int_t ldwa,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsyevdx_gpu(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double vl, double vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu,
    magma_minproduct_int_t *m, double *w,
    double *wA,  magma_minproduct_int_t ldwa,
    double *work, magma_minproduct_int_t lwork,
    #ifdef COMPLEX
    double *rwork, magma_minproduct_int_t lrwork,
    #endif
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

#ifdef COMPLEX
// no real [sd] precisions available
magma_minproduct_int_t
magma_minproduct_dsyevr_gpu(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double vl, double vu,
    magma_minproduct_int_t il, magma_minproduct_int_t iu, double abstol, magma_minproduct_int_t *m,
    double *w,
    magma_minproductDouble_ptr dZ, magma_minproduct_int_t lddz,
    magma_minproduct_int_t *isuppz,
    double *wA, magma_minproduct_int_t ldwa,
    double *wZ, magma_minproduct_int_t ldwz,
    double *work, magma_minproduct_int_t lwork,
    double *rwork, magma_minproduct_int_t lrwork,
    magma_minproduct_int_t *iwork, magma_minproduct_int_t liwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsyevx_gpu(
    magma_minproduct_vec_t jobz, magma_minproduct_range_t range, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double vl, double vu, magma_minproduct_int_t il, magma_minproduct_int_t iu,
    double abstol, magma_minproduct_int_t *m,
    double *w,
    magma_minproductDouble_ptr dZ, magma_minproduct_int_t lddz,
    double *wA, magma_minproduct_int_t ldwa,
    double *wZ, magma_minproduct_int_t ldwz,
    double *work, magma_minproduct_int_t lwork,
    double *rwork, magma_minproduct_int_t *iwork,
    magma_minproduct_int_t *ifail,
    magma_minproduct_int_t *info);
#endif  // COMPLEX

magma_minproduct_int_t
magma_minproduct_dsygst_gpu(
    magma_minproduct_int_t itype, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsytrd_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *d, double *e, double *tau,
    double *wA,  magma_minproduct_int_t ldwa,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsytrd_sy2sb_mgpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    double *work, magma_minproduct_int_t lwork,
    magma_minproductDouble_ptr dAmgpu[], magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dTmgpu[], magma_minproduct_int_t lddt,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t distblk,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsytrd_sy2sb_mgpu_spec(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    double *A, magma_minproduct_int_t lda,
    double *tau,
    double *work, magma_minproduct_int_t lwork,
    magma_minproductDouble_ptr dAmgpu[], magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dTmgpu[], magma_minproduct_int_t lddt,
    magma_minproduct_int_t ngpu, magma_minproduct_int_t distblk,
    magma_minproduct_queue_t queues[][20], magma_minproduct_int_t nqueue,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsytrd_mgpu(
    magma_minproduct_int_t ngpu, magma_minproduct_int_t nqueue,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double *A, magma_minproduct_int_t lda,
    double *d, double *e, double *tau,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsytrd2_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *d, double *e, double *tau,
    double *wA,  magma_minproduct_int_t ldwa,
    double *work, magma_minproduct_int_t lwork,
    magma_minproductDouble_ptr dwork, magma_minproduct_int_t ldwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dsytrf_nopiv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dlabrd_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb,
    double     *A, magma_minproduct_int_t lda,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *d, double *e, double *tauq, double *taup,
    double     *X, magma_minproduct_int_t ldx,
    magma_minproductDouble_ptr dX, magma_minproduct_int_t lddx,
    double     *Y, magma_minproduct_int_t ldy,
    magma_minproductDouble_ptr dY, magma_minproduct_int_t lddy );

magma_minproduct_int_t
magma_minproduct_dlaqps_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t offset,
    magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    magma_minproductDouble_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt, double *tau,
    double *vn1, double *vn2,
    magma_minproductDouble_ptr dauxv,
    magma_minproductDouble_ptr dF, magma_minproduct_int_t lddf);

magma_minproduct_int_t
magma_minproduct_dlaqps2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t offset,
    magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    magma_minproductDouble_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt,
    magma_minproductDouble_ptr dtau,
    magma_minproductDouble_ptr dvn1, magma_minproductDouble_ptr dvn2,
    magma_minproductDouble_ptr dauxv,
    magma_minproductDouble_ptr dF, magma_minproduct_int_t lddf);

magma_minproduct_int_t
magma_minproduct_dlaqps3_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t offset,
    magma_minproduct_int_t nb, magma_minproduct_int_t *kb,
    magma_minproductDouble_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *jpvt,
    magma_minproductDouble_ptr dtau,
    magma_minproductDouble_ptr dvn1, magma_minproductDouble_ptr dvn2,
    magma_minproductDouble_ptr dauxv,
    magma_minproductDouble_ptr dF, magma_minproduct_int_t lddf);

magma_minproduct_int_t
magma_minproduct_dlarf_gpu(
    magma_minproduct_int_t m,  magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dv, magma_minproductDouble_const_ptr dtau,
    magma_minproductDouble_ptr dC,  magma_minproduct_int_t lddc);

magma_minproduct_int_t
magma_minproduct_dlarfb_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_direct_t direct, magma_minproduct_storev_t storev,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDouble_const_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductDouble_const_ptr dT, magma_minproduct_int_t lddt,
    magma_minproductDouble_ptr dC,       magma_minproduct_int_t lddc,
    magma_minproductDouble_ptr dwork,    magma_minproduct_int_t ldwork );

magma_minproduct_int_t
magma_minproduct_dlarfb_gpu_gemm(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_direct_t direct, magma_minproduct_storev_t storev,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDouble_const_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductDouble_const_ptr dT, magma_minproduct_int_t lddt,
    magma_minproductDouble_ptr dC,       magma_minproduct_int_t lddc,
    magma_minproductDouble_ptr dwork,    magma_minproduct_int_t ldwork,
    magma_minproductDouble_ptr dworkvt,  magma_minproduct_int_t ldworkvt);

magma_minproduct_int_t
magma_minproduct_dlarfb2_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDouble_const_ptr dV, magma_minproduct_int_t lddv,
    magma_minproductDouble_const_ptr dT, magma_minproduct_int_t lddt,
    magma_minproductDouble_ptr dC,       magma_minproduct_int_t lddc,
    magma_minproductDouble_ptr dwork,    magma_minproduct_int_t ldwork );

magma_minproduct_int_t
magma_minproduct_dlatrd_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t nb0,
    double *A,  magma_minproduct_int_t lda,
    double *e, double *tau,
    double    *W,       magma_minproduct_int_t ldw,
    magma_minproductDouble_ptr dA[],    magma_minproduct_int_t ldda, magma_minproduct_int_t offset,
    magma_minproductDouble_ptr dW[],    magma_minproduct_int_t lddw,
    double    *hwork,   magma_minproduct_int_t lhwork,
    magma_minproductDouble_ptr dwork[], magma_minproduct_int_t ldwork,
    magma_minproduct_queue_t queues[] );

magma_minproduct_int_t
magma_minproduct_dlauum_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dposv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dpotf2_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dpotrf_gpu(
    magma_minproduct_uplo_t uplo,  magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dpotrf_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDouble_ptr d_lA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dpotrf_mgpu_right(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDouble_ptr d_lA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dpotrf3_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t off_i, magma_minproduct_int_t off_j, magma_minproduct_int_t nb,
    magma_minproductDouble_ptr d_lA[],  magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr d_lP[],  magma_minproduct_int_t lddp,
    double *A, magma_minproduct_int_t lda, magma_minproduct_int_t h,
    magma_minproduct_queue_t queues[][3], magma_minproduct_event_t events[][5],
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dpotri_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dpotrs_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dssssm_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m1, magma_minproduct_int_t n1,
    magma_minproduct_int_t m2, magma_minproduct_int_t n2, magma_minproduct_int_t k, magma_minproduct_int_t ib,
    magma_minproductDouble_ptr dA1, magma_minproduct_int_t ldda1,
    magma_minproductDouble_ptr dA2, magma_minproduct_int_t ldda2,
    magma_minproductDouble_ptr dL1, magma_minproduct_int_t lddl1,
    magma_minproductDouble_ptr dL2, magma_minproduct_int_t lddl2,
    magma_minproduct_int_t *ipiv,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dtrtri_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dtsqrt_gpu(
    magma_minproduct_int_t *m, magma_minproduct_int_t *n,
    double *A1, double *A2, magma_minproduct_int_t *lda,
    double *tau,
    double *work, magma_minproduct_int_t *lwork,
    magma_minproductDouble_ptr dwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dtstrf_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t ib, magma_minproduct_int_t nb,
    double    *hU, magma_minproduct_int_t ldhu,
    magma_minproductDouble_ptr dU, magma_minproduct_int_t lddu,
    double    *hA, magma_minproduct_int_t ldha,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double    *hL, magma_minproduct_int_t ldhl,
    magma_minproductDouble_ptr dL, magma_minproduct_int_t lddl,
    magma_minproduct_int_t *ipiv,
    double *hwork, magma_minproduct_int_t ldhwork,
    magma_minproductDouble_ptr dwork, magma_minproduct_int_t lddwork,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dorgqr_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *tau,
    magma_minproductDouble_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dormql2_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *tau,
    magma_minproductDouble_ptr dC, magma_minproduct_int_t lddc,
    double *wA, magma_minproduct_int_t ldwa,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dormqr_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *tau,
    magma_minproductDouble_ptr dC, magma_minproduct_int_t lddc,
    double *hwork, magma_minproduct_int_t lwork,
    magma_minproductDouble_ptr dT, magma_minproduct_int_t nb,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dormqr2_gpu(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *tau,
    magma_minproductDouble_ptr dC, magma_minproduct_int_t lddc,
    double    *wA, magma_minproduct_int_t ldwa,
    magma_minproduct_int_t *info);

magma_minproduct_int_t
magma_minproduct_dormtr_gpu(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *tau,
    magma_minproductDouble_ptr dC, magma_minproduct_int_t lddc,
    double    *wA, magma_minproduct_int_t ldwa,
    magma_minproduct_int_t *info);


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_minproduct utility function definitions
*/

extern const double MAGMA_minproduct_D_NAN;
extern const double MAGMA_minproduct_D_INF;

magma_minproduct_int_t
magma_minproduct_dnan_inf(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    const double *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *cnt_nan,
    magma_minproduct_int_t *cnt_inf );

magma_minproduct_int_t
magma_minproduct_dnan_inf_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *cnt_nan,
    magma_minproduct_int_t *cnt_inf );

void magma_minproduct_dprint(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    const double *A, magma_minproduct_int_t lda );

void magma_minproduct_dprint_gpu(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda );

void dpanel_to_q(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t ib,
    double *A, magma_minproduct_int_t lda,
    double *work );

void dq_to_panel(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t ib,
    double *A, magma_minproduct_int_t lda,
    double *work );

#ifdef __cplusplus
}
#endif

#undef REAL

#endif /* MAGMA_minproduct_D_H */
