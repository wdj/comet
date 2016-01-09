/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#ifndef MAGMA_tally4_AUXILIARY_H
#define MAGMA_tally4_AUXILIARY_H

#include "magma_tally4_types.h"

/* ------------------------------------------------------------
 *   -- MAGMA_tally4 Auxiliary structures and functions
 * --------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

real_Double_t magma_tally4_wtime( void );
real_Double_t magma_tally4_sync_wtime( magma_tally4_queue_t queue );

size_t magma_tally4_strlcpy(char *dst, const char *src, size_t siz);

magma_tally4_int_t magma_tally4_num_gpus( void );

double magma_tally4_cabs(magma_tally4DoubleComplex x);
float  magma_tally4_cabsf(magma_tally4FloatComplex x);

magma_tally4_int_t magma_tally4_is_devptr( const void* A );

// magma_tally4 GPU-complex PCIe connection
magma_tally4_int_t magma_tally4_buildconnection_mgpu(  magma_tally4_int_t gnode[Magma_tally4MaxGPUs+2][Magma_tally4MaxGPUs+2], magma_tally4_int_t *nbcmplx, magma_tally4_int_t ngpu);

void magma_tally4_indices_1D_bcyclic( magma_tally4_int_t nb, magma_tally4_int_t ngpu, magma_tally4_int_t dev,
                               magma_tally4_int_t j0, magma_tally4_int_t j1,
                               magma_tally4_int_t* dj0, magma_tally4_int_t* dj1 );

void magma_tally4_print_environment();

void swp2pswp_tally4(magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t *ipiv, magma_tally4_int_t *newipiv);

#ifdef __cplusplus
}
#endif

#endif  // MAGMA_tally4_AUXILIARY_H
