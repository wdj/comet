/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#ifndef MAGMA_tally2_AUXILIARY_H
#define MAGMA_tally2_AUXILIARY_H

#include "magma_tally2_types.h"

/* ------------------------------------------------------------
 *   -- MAGMA_tally2 Auxiliary structures and functions
 * --------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

real_Double_t magma_tally2_wtime( void );
real_Double_t magma_tally2_sync_wtime( magma_tally2_queue_t queue );

size_t magma_tally2_strlcpy(char *dst, const char *src, size_t siz);

magma_tally2_int_t magma_tally2_num_gpus( void );

double magma_tally2_cabs(magma_tally2DoubleComplex x);
float  magma_tally2_cabsf(magma_tally2FloatComplex x);

magma_tally2_int_t magma_tally2_is_devptr( const void* A );

// magma_tally2 GPU-complex PCIe connection
magma_tally2_int_t magma_tally2_buildconnection_mgpu(  magma_tally2_int_t gnode[Magma_tally2MaxGPUs+2][Magma_tally2MaxGPUs+2], magma_tally2_int_t *nbcmplx, magma_tally2_int_t ngpu);

void magma_tally2_indices_1D_bcyclic( magma_tally2_int_t nb, magma_tally2_int_t ngpu, magma_tally2_int_t dev,
                               magma_tally2_int_t j0, magma_tally2_int_t j1,
                               magma_tally2_int_t* dj0, magma_tally2_int_t* dj1 );

void magma_tally2_print_environment();

void swp2pswp_tally2(magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t *ipiv, magma_tally2_int_t *newipiv);

#ifdef __cplusplus
}
#endif

#endif  // MAGMA_tally2_AUXILIARY_H
