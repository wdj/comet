/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions mixed zc -> ds
*/

#ifndef MAGMA_tally2_ZC_H
#define MAGMA_tally2_ZC_H

#include "magma_tally2_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Mixed precision */
magma_tally2_int_t
magma_tally2_zcgesv_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2_int_t *ipiv,
    magma_tally2Int_ptr dipiv,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex_ptr dX, magma_tally2_int_t lddx,
    magma_tally2DoubleComplex_ptr dworkd, magma_tally2FloatComplex_ptr dworks,
    magma_tally2_int_t *iter,
    magma_tally2_int_t *info );

magma_tally2_int_t
magma_tally2_zcgetrs_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex_ptr  dA, magma_tally2_int_t ldda,
    magma_tally2Int_ptr        dipiv,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex_ptr dX, magma_tally2_int_t lddx,
    magma_tally2FloatComplex_ptr dSX,
    magma_tally2_int_t *info );

magma_tally2_int_t
magma_tally2_zcposv_gpu(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex_ptr dX, magma_tally2_int_t lddx,
    magma_tally2DoubleComplex_ptr dworkd, magma_tally2FloatComplex_ptr dworks,
    magma_tally2_int_t *iter,
    magma_tally2_int_t *info );

magma_tally2_int_t
magma_tally2_zcgeqrsv_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr dB, magma_tally2_int_t lddb,
    magma_tally2DoubleComplex_ptr dX, magma_tally2_int_t lddx,
    magma_tally2_int_t *iter,
    magma_tally2_int_t *info );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally2_ZC_H */
