/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#ifndef MAGMA_tally4BATCHED_H
#define MAGMA_tally4BATCHED_H
#include "magma_tally4_types.h"

/* ------------------------------------------------------------
 * MAGMA_tally4 BATCHED functions
 * --------------------------------------------------------- */
#include "magma_tally4_zbatched.h"
#include "magma_tally4_cbatched.h"
#include "magma_tally4_dbatched.h"
#include "magma_tally4_sbatched.h"



#ifdef __cplusplus
extern "C" {
#endif


void setup_pivinfo_batched( magma_tally4_int_t **pivinfo_array, magma_tally4_int_t **ipiv_array, 
                              magma_tally4_int_t m, magma_tally4_int_t nb, 
                              magma_tally4_int_t batchCount,  magma_tally4_queue_t queue);


void adjust_ipiv_batched( magma_tally4_int_t **ipiv_array, 
                       magma_tally4_int_t m, magma_tally4_int_t offset, 
                       magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void magma_tally4_idisplace_pointers(magma_tally4_int_t **output_array,
               magma_tally4_int_t **input_array, magma_tally4_int_t lda,
               magma_tally4_int_t row, magma_tally4_int_t column, 
               magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void stepinit_ipiv(magma_tally4_int_t **ipiv_array,
                 magma_tally4_int_t pm,
                 magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void set_ipointer(magma_tally4_int_t **output_array,
                 magma_tally4_int_t *input,
                 magma_tally4_int_t lda,
                 magma_tally4_int_t row, magma_tally4_int_t column, 
                 magma_tally4_int_t batchSize,
                 magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

#ifdef __cplusplus
}
#endif


#endif /* MAGMA_tally4BATCHED_H */
