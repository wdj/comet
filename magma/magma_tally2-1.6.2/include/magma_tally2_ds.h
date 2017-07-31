/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally2_zc.h mixed zc -> ds, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_tally2_DS_H
#define MAGMA_tally2_DS_H

#include "magma_tally2_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Mixed precision */
magma_tally2_int_t
magma_tally2_dsgesv_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *ipiv,
    magma_tally2Int_ptr dipiv,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    magma_tally2Double_ptr dX, magma_tally2_int_t lddx,
    magma_tally2Double_ptr dworkd, magma_tally2Float_ptr dworks,
    magma_tally2_int_t *iter,
    magma_tally2_int_t *info );

magma_tally2_int_t
magma_tally2_dsgetrs_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Float_ptr  dA, magma_tally2_int_t ldda,
    magma_tally2Int_ptr        dipiv,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    magma_tally2Double_ptr dX, magma_tally2_int_t lddx,
    magma_tally2Float_ptr dSX,
    magma_tally2_int_t *info );

magma_tally2_int_t
magma_tally2_dsposv_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    magma_tally2Double_ptr dX, magma_tally2_int_t lddx,
    magma_tally2Double_ptr dworkd, magma_tally2Float_ptr dworks,
    magma_tally2_int_t *iter,
    magma_tally2_int_t *info );

magma_tally2_int_t
magma_tally2_dsgeqrsv_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    magma_tally2Double_ptr dX, magma_tally2_int_t lddx,
    magma_tally2_int_t *iter,
    magma_tally2_int_t *info );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally2_DS_H */
