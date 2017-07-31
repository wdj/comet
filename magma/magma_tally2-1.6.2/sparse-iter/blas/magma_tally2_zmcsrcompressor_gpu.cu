/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt

*/
#include "common_magma_tally2sparse.h"

#define BLOCK_SIZE1 256
#define BLOCK_SIZE2 1


// copy nonzeros into new structure
__global__ void
magma_tally2_zmcsrgpu_kernel1( int num_rows,
                 magma_tally2DoubleComplex *A_val,
                 magma_tally2_index_t *A_rowptr,
                 magma_tally2_index_t *A_colind,
                 magma_tally2DoubleComplex *B_val,
                 magma_tally2_index_t *B_rowptr,
                 magma_tally2_index_t *B_colind ){

    int row = blockIdx.x*blockDim.x+threadIdx.x;
    int j;

    if(row<num_rows){
        magma_tally2DoubleComplex zero = MAGMA_tally2_Z_ZERO;
        int start = A_rowptr[ row ];
        int new_location = start;
        int end = A_rowptr[ row+1 ];
        for( j=start; j<end; j++ ){
            if( A_val[j] != zero ){
       //         B_val[new_location] = A_val[j];
       //         B_colind[new_location] = A_colind[j];
                new_location++;
            }
        }
        // this is not a correctr rowpointer! this is nn_z in this row!
        B_rowptr[ row ] = new_location-start;
    }
}


// generate a valid rowpointer
__global__ void
magma_tally2_zmcsrgpu_kernel2( int num_rows,
                 magma_tally2_index_t *B_rowptr,
                 magma_tally2_index_t *A_rowptr ){

    int idx = blockIdx.x*blockDim.x+threadIdx.x;
    int j, nnz = 0;

    if( idx == 0 ){
    A_rowptr[ 0 ] = nnz;
        for( j=0; j<num_rows; j++ ){
            nnz+=B_rowptr[ j ];
            A_rowptr[ j+1 ] = nnz;
        }
    }
}



// copy new structure into original matrix
__global__ void
magma_tally2_zmcsrgpu_kernel3( int num_rows,
                 magma_tally2DoubleComplex *B_val,
                 magma_tally2_index_t *B_rowptr,
                 magma_tally2_index_t *B_colind,
                 magma_tally2_index_t *B2_rowptr,
                 magma_tally2DoubleComplex *A_val,
                 magma_tally2_index_t *A_rowptr,
                 magma_tally2_index_t *A_colind
                                            ){

    int row = blockIdx.x*blockDim.x+threadIdx.x;
    int j, new_location;
    
    if(row<num_rows){
    new_location = A_rowptr[ row ];
        int start = B2_rowptr[ row ];
        int end = B2_rowptr[ row+1 ];
        magma_tally2DoubleComplex zero = MAGMA_tally2_Z_ZERO;
        for( j=start; j<end; j++ ){
            if( A_val[j] != zero ){
                B_val[new_location] = A_val[j];
                B_colind[new_location] = A_colind[j];
                new_location++;
            }
               // A_val[ j ] = B_val[ j ];
               // A_colind[ j ] = B_colind[ j ];
        }
    }
}


/**
    Purpose
    -------

    Removes zeros in a CSR matrix. This is a GPU implementation of the
    CSR compressor.

    Arguments
    ---------

    @param
    A           magma_tally2_z_matrix*
                input/output matrix
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_zaux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_zmcsrcompressor_gpu(
    magma_tally2_z_matrix *A,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    magma_tally2_z_matrix B={Magma_tally2_CSR}, B2={Magma_tally2_CSR};
    magma_tally2_z_matrix dA={Magma_tally2_CSR}, CSRA={Magma_tally2_CSR};
    magma_tally2_index_t *cputmp = NULL;
    
    if ( A->memory_location == Magma_tally2_DEV && A->storage_type == Magma_tally2_CSR ) {

        CHECK( magma_tally2_index_malloc( &B.drow, A->num_rows + 1 ));
        CHECK( magma_tally2_index_malloc( &B2.drow, A->num_rows + 1 ));
        
        magma_tally2_index_copyvector( (A->num_rows+1), A->drow, 1, B2.drow, 1 );

        dim3 grid1( magma_tally2_ceildiv( A->num_rows, BLOCK_SIZE1 ) );

        // copying the nonzeros into B and write in B.drow how many there are
        magma_tally2_zmcsrgpu_kernel1<<< grid1, BLOCK_SIZE1, 0, queue >>>
                ( A->num_rows, A->dval, A->drow, A->dcol, B.dval, B.drow, B.dcol );

        // correct the row pointer
        dim3 grid2( 1, 1, 1);
        magma_tally2_zmcsrgpu_kernel2<<< grid2, BLOCK_SIZE2, 0, queue >>>
                ( A->num_rows, B.drow, A->drow );
        // access the true number of nonzeros

        CHECK( magma_tally2_index_malloc_cpu( &cputmp, 1 ));

        magma_tally2_index_getvector( 1, A->row+(A->num_rows), 1, cputmp, 1 );
        A->nnz = (magma_tally2_int_t) cputmp[0];

        // reallocate with right size
        CHECK( magma_tally2_zmalloc( &B.dval, A->nnz ));
        CHECK( magma_tally2_index_malloc( &B.dcol, A->nnz ));
        
        // copy correct values back
        magma_tally2_zmcsrgpu_kernel3<<< grid1, BLOCK_SIZE1, 0, queue >>>
                ( A->num_rows, B.dval, B.drow, B.dcol, B2.drow, A->dval, A->drow, A->dcol );

        magma_tally2_free( A->dcol );
        magma_tally2_free( A->dval );

        A->dcol = B.dcol;
        A->dval = B.dval;


    }
    else {
        magma_tally2_storage_t A_storage = A->storage_type;
        magma_tally2_location_t A_location = A->memory_location;
        CHECK( magma_tally2_zmconvert( *A, &CSRA, A->storage_type, Magma_tally2_CSR, queue ));
        CHECK( magma_tally2_zmtransfer( *A, &dA, A->memory_location, Magma_tally2_DEV, queue ));

        CHECK( magma_tally2_zmcsrcompressor_gpu( &dA, queue ));

        magma_tally2_zmfree( &dA, queue );
        magma_tally2_zmfree( A, queue );
        CHECK( magma_tally2_zmtransfer( dA, &CSRA, Magma_tally2_DEV, A_location, queue ));
        CHECK( magma_tally2_zmconvert( CSRA, A, Magma_tally2_CSR, A_storage, queue ));
        magma_tally2_zmfree( &dA, queue );
        magma_tally2_zmfree( &CSRA, queue );

    }
    
cleanup:
    magma_tally2_zmfree( &dA, queue );
    magma_tally2_zmfree( &CSRA, queue );
    magma_tally2_free( B2.drow );
    magma_tally2_free( B.drow );
    return info;
}


