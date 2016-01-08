/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally4_zc.h mixed zc -> ds, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_tally4_DS_H
#define MAGMA_tally4_DS_H

#include "magma_tally4_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Mixed precision */
magma_tally4_int_t
magma_tally4_dsgesv_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *ipiv,
    magma_tally4Int_ptr dipiv,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    magma_tally4Double_ptr dX, magma_tally4_int_t lddx,
    magma_tally4Double_ptr dworkd, magma_tally4Float_ptr dworks,
    magma_tally4_int_t *iter,
    magma_tally4_int_t *info );

magma_tally4_int_t
magma_tally4_dsgetrs_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Float_ptr  dA, magma_tally4_int_t ldda,
    magma_tally4Int_ptr        dipiv,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    magma_tally4Double_ptr dX, magma_tally4_int_t lddx,
    magma_tally4Float_ptr dSX,
    magma_tally4_int_t *info );

magma_tally4_int_t
magma_tally4_dsposv_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    magma_tally4Double_ptr dX, magma_tally4_int_t lddx,
    magma_tally4Double_ptr dworkd, magma_tally4Float_ptr dworks,
    magma_tally4_int_t *iter,
    magma_tally4_int_t *info );

magma_tally4_int_t
magma_tally4_dsgeqrsv_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    magma_tally4Double_ptr dX, magma_tally4_int_t lddx,
    magma_tally4_int_t *iter,
    magma_tally4_int_t *info );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally4_DS_H */
