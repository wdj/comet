/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally3_zc.h mixed zc -> ds, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_tally3_DS_H
#define MAGMA_tally3_DS_H

#include "magma_tally3_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Mixed precision */
magma_tally3_int_t
magma_tally3_dsgesv_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *ipiv,
    magma_tally3Int_ptr dipiv,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb,
    magma_tally3Double_ptr dX, magma_tally3_int_t lddx,
    magma_tally3Double_ptr dworkd, magma_tally3Float_ptr dworks,
    magma_tally3_int_t *iter,
    magma_tally3_int_t *info );

magma_tally3_int_t
magma_tally3_dsgetrs_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Float_ptr  dA, magma_tally3_int_t ldda,
    magma_tally3Int_ptr        dipiv,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb,
    magma_tally3Double_ptr dX, magma_tally3_int_t lddx,
    magma_tally3Float_ptr dSX,
    magma_tally3_int_t *info );

magma_tally3_int_t
magma_tally3_dsposv_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb,
    magma_tally3Double_ptr dX, magma_tally3_int_t lddx,
    magma_tally3Double_ptr dworkd, magma_tally3Float_ptr dworks,
    magma_tally3_int_t *iter,
    magma_tally3_int_t *info );

magma_tally3_int_t
magma_tally3_dsgeqrsv_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Double_ptr dB, magma_tally3_int_t lddb,
    magma_tally3Double_ptr dX, magma_tally3_int_t lddx,
    magma_tally3_int_t *iter,
    magma_tally3_int_t *info );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally3_DS_H */
