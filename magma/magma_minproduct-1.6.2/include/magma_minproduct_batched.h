/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#ifndef MAGMA_minproductBATCHED_H
#define MAGMA_minproductBATCHED_H
#include "magma_minproduct_types.h"

/* ------------------------------------------------------------
 * MAGMA_minproduct BATCHED functions
 * --------------------------------------------------------- */
#include "magma_minproduct_zbatched.h"
#include "magma_minproduct_cbatched.h"
#include "magma_minproduct_dbatched.h"
#include "magma_minproduct_sbatched.h"



#ifdef __cplusplus
extern "C" {
#endif


void setup_pivinfo_batched( magma_minproduct_int_t **pivinfo_array, magma_minproduct_int_t **ipiv_array, 
                              magma_minproduct_int_t m, magma_minproduct_int_t nb, 
                              magma_minproduct_int_t batchCount,  magma_minproduct_queue_t queue);


void adjust_ipiv_batched( magma_minproduct_int_t **ipiv_array, 
                       magma_minproduct_int_t m, magma_minproduct_int_t offset, 
                       magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproduct_idisplace_pointers(magma_minproduct_int_t **output_array,
               magma_minproduct_int_t **input_array, magma_minproduct_int_t lda,
               magma_minproduct_int_t row, magma_minproduct_int_t column, 
               magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void stepinit_ipiv(magma_minproduct_int_t **ipiv_array,
                 magma_minproduct_int_t pm,
                 magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void set_ipointer(magma_minproduct_int_t **output_array,
                 magma_minproduct_int_t *input,
                 magma_minproduct_int_t lda,
                 magma_minproduct_int_t row, magma_minproduct_int_t column, 
                 magma_minproduct_int_t batchSize,
                 magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

#ifdef __cplusplus
}
#endif


#endif /* MAGMA_minproductBATCHED_H */
