/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions mixed zc -> ds
*/

#ifndef MAGMA_tally3_ZC_H
#define MAGMA_tally3_ZC_H

#include "magma_tally3_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Mixed precision */
magma_tally3_int_t
magma_tally3_zcgesv_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3_int_t *ipiv,
    magma_tally3Int_ptr dipiv,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3DoubleComplex_ptr dX, magma_tally3_int_t lddx,
    magma_tally3DoubleComplex_ptr dworkd, magma_tally3FloatComplex_ptr dworks,
    magma_tally3_int_t *iter,
    magma_tally3_int_t *info );

magma_tally3_int_t
magma_tally3_zcgetrs_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3FloatComplex_ptr  dA, magma_tally3_int_t ldda,
    magma_tally3Int_ptr        dipiv,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3DoubleComplex_ptr dX, magma_tally3_int_t lddx,
    magma_tally3FloatComplex_ptr dSX,
    magma_tally3_int_t *info );

magma_tally3_int_t
magma_tally3_zcposv_gpu(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3DoubleComplex_ptr dX, magma_tally3_int_t lddx,
    magma_tally3DoubleComplex_ptr dworkd, magma_tally3FloatComplex_ptr dworks,
    magma_tally3_int_t *iter,
    magma_tally3_int_t *info );

magma_tally3_int_t
magma_tally3_zcgeqrsv_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb,
    magma_tally3DoubleComplex_ptr dX, magma_tally3_int_t lddx,
    magma_tally3_int_t *iter,
    magma_tally3_int_t *info );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally3_ZC_H */
