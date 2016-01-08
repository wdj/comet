/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2013
       
       @generated from zset_pointer.cu normal z -> d, Fri Jan 30 19:00:10 2015
       @author Azzam Haidar
       @author Tingxing Dong

*/

#include "common_magma_tally4.h"
///////////////////////////////////////////////////////////////////////////////////////
__global__ void kernel_set_pointer(double **output_array,
                 double *input,
                 magma_tally4_int_t lda,
                 magma_tally4_int_t row, magma_tally4_int_t column, 
                 magma_tally4_int_t batchSize)
{
     output_array[blockIdx.x] =  input + blockIdx.x * batchSize + row + column * lda;
     //printf("==> kernel_set_pointer input_array %p output_array %p  \n",input+ blockIdx.x * batchSize,output_array[blockIdx.x]);
}


extern "C"
void dset_pointer(double **output_array,
                 double *input,
                 magma_tally4_int_t lda,
                 magma_tally4_int_t row, magma_tally4_int_t column, 
                 magma_tally4_int_t batchSize,
                 magma_tally4_int_t batchCount, 
                 magma_tally4_queue_t queue)

{
/*
    convert consecutive stored variable to array stored
    for example the size  of A is N*batchCount; N is the size of A(batchSize)
    change into dA_array[0] dA_array[1],... dA_array[batchCount-1], where the size of each dA_array[i] is N
*/


    kernel_set_pointer<<<batchCount, 1, 0, queue>>>(output_array, input, lda,  row, column, batchSize);
}



__global__ void zdisplace_pointers_kernel(double **output_array,
               double **input_array, magma_tally4_int_t lda,
               magma_tally4_int_t row, magma_tally4_int_t column)
{
      double *inpt = input_array[blockIdx.x];
      output_array[blockIdx.x] = &inpt[row + column * lda];
}


extern "C"
void dset_array(double **output_array,
               double **input_array, magma_tally4_int_t lda,
               magma_tally4_int_t row, magma_tally4_int_t column, 
               magma_tally4_int_t batchCount, magma_tally4_queue_t queue)

{
/*
    compute the offset for all the matrices and save the displacment of the new pointer on output_array.
    input_array contains the pointers to the initial position.
    output_array[i] = input_array[i] + row + lda * column; 
*/
    zdisplace_pointers_kernel<<<batchCount, 1, 0, queue>>>(output_array, input_array, lda, row, column);
}

extern "C"
void magma_tally4_ddisplace_pointers(double **output_array,
               double **input_array, magma_tally4_int_t lda,
               magma_tally4_int_t row, magma_tally4_int_t column, 
               magma_tally4_int_t batchCount, magma_tally4_queue_t queue)

{
/*
    compute the offset for all the matrices and save the displacment of the new pointer on output_array.
    input_array contains the pointers to the initial position.
    output_array[i] = input_array[i] + row + lda * column; 
*/
    zdisplace_pointers_kernel<<<batchCount, 1, 0, queue>>>(output_array, input_array, lda, row, column);
}



