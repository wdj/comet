/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions mixed zc -> ds
*/

#ifndef MAGMA_tally4_ZC_H
#define MAGMA_tally4_ZC_H

#include "magma_tally4_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Mixed precision */
magma_tally4_int_t
magma_tally4_zcgesv_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4_int_t *ipiv,
    magma_tally4Int_ptr dipiv,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex_ptr dX, magma_tally4_int_t lddx,
    magma_tally4DoubleComplex_ptr dworkd, magma_tally4FloatComplex_ptr dworks,
    magma_tally4_int_t *iter,
    magma_tally4_int_t *info );

magma_tally4_int_t
magma_tally4_zcgetrs_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex_ptr  dA, magma_tally4_int_t ldda,
    magma_tally4Int_ptr        dipiv,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex_ptr dX, magma_tally4_int_t lddx,
    magma_tally4FloatComplex_ptr dSX,
    magma_tally4_int_t *info );

magma_tally4_int_t
magma_tally4_zcposv_gpu(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex_ptr dX, magma_tally4_int_t lddx,
    magma_tally4DoubleComplex_ptr dworkd, magma_tally4FloatComplex_ptr dworks,
    magma_tally4_int_t *iter,
    magma_tally4_int_t *info );

magma_tally4_int_t
magma_tally4_zcgeqrsv_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex_ptr dX, magma_tally4_int_t lddx,
    magma_tally4_int_t *iter,
    magma_tally4_int_t *info );

#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally4_ZC_H */
