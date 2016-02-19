/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#ifndef MAGMA_tally3_AUXILIARY_H
#define MAGMA_tally3_AUXILIARY_H

#include "magma_tally3_types.h"

/* ------------------------------------------------------------
 *   -- MAGMA_tally3 Auxiliary structures and functions
 * --------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

real_Double_t magma_tally3_wtime( void );
real_Double_t magma_tally3_sync_wtime( magma_tally3_queue_t queue );

size_t magma_tally3_strlcpy(char *dst, const char *src, size_t siz);

magma_tally3_int_t magma_tally3_num_gpus( void );

double magma_tally3_cabs(magma_tally3DoubleComplex x);
float  magma_tally3_cabsf(magma_tally3FloatComplex x);

magma_tally3_int_t magma_tally3_is_devptr( const void* A );

// magma_tally3 GPU-complex PCIe connection
magma_tally3_int_t magma_tally3_buildconnection_mgpu(  magma_tally3_int_t gnode[Magma_tally3MaxGPUs+2][Magma_tally3MaxGPUs+2], magma_tally3_int_t *nbcmplx, magma_tally3_int_t ngpu);

void magma_tally3_indices_1D_bcyclic( magma_tally3_int_t nb, magma_tally3_int_t ngpu, magma_tally3_int_t dev,
                               magma_tally3_int_t j0, magma_tally3_int_t j1,
                               magma_tally3_int_t* dj0, magma_tally3_int_t* dj1 );

void magma_tally3_print_environment();

void swp2pswp_tally3(magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t *ipiv, magma_tally3_int_t *newipiv);

#ifdef __cplusplus
}
#endif

#endif  // MAGMA_tally3_AUXILIARY_H
