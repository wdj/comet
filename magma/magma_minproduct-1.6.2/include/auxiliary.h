/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#ifndef MAGMA_minproduct_AUXILIARY_H
#define MAGMA_minproduct_AUXILIARY_H

#include "magma_minproduct_types.h"

/* ------------------------------------------------------------
 *   -- MAGMA_minproduct Auxiliary structures and functions
 * --------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

real_Double_t magma_minproduct_wtime( void );
real_Double_t magma_minproduct_sync_wtime( magma_minproduct_queue_t queue );

size_t magma_minproduct_strlcpy(char *dst, const char *src, size_t siz);

magma_minproduct_int_t magma_minproduct_num_gpus( void );

double magma_minproduct_cabs(magma_minproductDoubleComplex x);
float  magma_minproduct_cabsf(magma_minproductFloatComplex x);

magma_minproduct_int_t magma_minproduct_is_devptr( const void* A );

// magma_minproduct GPU-complex PCIe connection
magma_minproduct_int_t magma_minproduct_buildconnection_mgpu(  magma_minproduct_int_t gnode[Magma_minproductMaxGPUs+2][Magma_minproductMaxGPUs+2], magma_minproduct_int_t *nbcmplx, magma_minproduct_int_t ngpu);

void magma_minproduct_indices_1D_bcyclic( magma_minproduct_int_t nb, magma_minproduct_int_t ngpu, magma_minproduct_int_t dev,
                               magma_minproduct_int_t j0, magma_minproduct_int_t j1,
                               magma_minproduct_int_t* dj0, magma_minproduct_int_t* dj1 );

void magma_minproduct_print_environment();

void swp2pswp(magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t *ipiv, magma_minproduct_int_t *newipiv);

#ifdef __cplusplus
}
#endif

#endif  // MAGMA_minproduct_AUXILIARY_H
