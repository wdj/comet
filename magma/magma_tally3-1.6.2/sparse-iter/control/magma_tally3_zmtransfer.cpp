/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt
*/
#include "common_magma_tally3sparse.h"


/**
    Purpose
    -------

    Copies a matrix from memory location src to memory location dst.


    Arguments
    ---------

    @param[in]
    A           magma_tally3_z_matrix
                sparse matrix A

    @param[out]
    B           magma_tally3_z_matrix*
                copy of A

    @param[in]
    src         magma_tally3_location_t
                original location A

    @param[in]
    dst         magma_tally3_location_t
                location of the copy of A

   
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zaux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_zmtransfer(
    magma_tally3_z_matrix A,
    magma_tally3_z_matrix *B,
    magma_tally3_location_t src,
    magma_tally3_location_t dst,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue=NULL;
    magma_tally3blasGetKernelStream( &orig_queue );

    B->val = NULL;
    B->col = NULL;
    B->row = NULL;
    B->rowidx = NULL;
    B->blockinfo = NULL;
    B->diag = NULL;
    B->dval = NULL;
    B->dcol = NULL;
    B->drow = NULL;
    B->drowidx = NULL;
    B->ddiag = NULL;

    // first case: copy matrix from host to device
    if ( src == Magma_tally3_CPU && dst == Magma_tally3_DEV ) {
        //CSR-type
        if ( A.storage_type == Magma_tally3_CSR || A.storage_type == Magma_tally3_CSC
                                        || A.storage_type == Magma_tally3_CSRD
                                        || A.storage_type == Magma_tally3_CSRL
                                        || A.storage_type == Magma_tally3_CSRU ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_DEV;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            // memory allocation
            CHECK( magma_tally3_zmalloc( &B->dval, A.nnz ));
            CHECK( magma_tally3_index_malloc( &B->drow, A.num_rows + 1 ));
            CHECK( magma_tally3_index_malloc( &B->dcol, A.nnz ));
            // data transfer
            magma_tally3_zsetvector( A.nnz, A.val, 1, B->dval, 1 );
            magma_tally3_index_setvector( A.num_rows+1, A.row, 1, B->drow, 1 );
            magma_tally3_index_setvector( A.nnz, A.col, 1, B->dcol, 1 );
        }
        //CSRCOO-type
        if ( A.storage_type == Magma_tally3_CSRCOO ) {
            // fill in information for B
            *B = A;
            B->memory_location = Magma_tally3_DEV;
            // memory allocation
            CHECK( magma_tally3_zmalloc( &B->dval, A.nnz ));
            CHECK( magma_tally3_index_malloc( &B->drow, A.num_rows + 1 ));
            CHECK( magma_tally3_index_malloc( &B->dcol, A.nnz ));
            CHECK( magma_tally3_index_malloc( &B->drowidx, A.nnz ));
            // data transfer
            magma_tally3_zsetvector( A.nnz, A.val, 1, B->dval, 1 );
            magma_tally3_index_setvector( A.num_rows+1, A.row, 1, B->drow, 1 );
            magma_tally3_index_setvector( A.nnz, A.col, 1, B->dcol, 1 );
            magma_tally3_index_setvector( A.nnz, A.rowidx, 1, B->drowidx, 1 );
        }
        //ELL/ELLPACKT-type
        if ( A.storage_type == Magma_tally3_ELLPACKT || A.storage_type == Magma_tally3_ELL ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_DEV;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            // memory allocation
            CHECK( magma_tally3_zmalloc( &B->dval, A.num_rows * A.max_nnz_row ));
            CHECK( magma_tally3_index_malloc( &B->dcol, A.num_rows * A.max_nnz_row ));
            // data transfer
            magma_tally3_zsetvector( A.num_rows * A.max_nnz_row, A.val, 1, B->dval, 1 );
            magma_tally3_index_setvector( A.num_rows * A.max_nnz_row, A.col, 1, B->dcol, 1 );
        }
        //ELLD-type
        if ( A.storage_type == Magma_tally3_ELLD ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_DEV;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            // memory allocation
            CHECK( magma_tally3_zmalloc( &B->dval, A.num_rows * A.max_nnz_row ));
            CHECK( magma_tally3_index_malloc( &B->dcol, A.num_rows * A.max_nnz_row ));
            // data transfer
            magma_tally3_zsetvector( A.num_rows * A.max_nnz_row, A.val, 1, B->dval, 1 );
            magma_tally3_index_setvector( A.num_rows * A.max_nnz_row, A.col, 1, B->dcol, 1 );
        }
        //ELLRT-type
        if ( A.storage_type == Magma_tally3_ELLRT ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_DEV;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            int threads_per_row = A.alignment;
            B->blocksize = A.blocksize;
            B->alignment = A.alignment;
            int rowlength = ( (int)((A.max_nnz_row+threads_per_row-1)
                                        /threads_per_row) ) * threads_per_row;
            // memory allocation
            CHECK( magma_tally3_zmalloc( &B->dval, A.num_rows * rowlength ));
            CHECK( magma_tally3_index_malloc( &B->dcol, A.num_rows * rowlength ));
            CHECK( magma_tally3_index_malloc( &B->drow, A.num_rows ));
            // data transfer
            magma_tally3_zsetvector( A.num_rows * rowlength, A.val, 1, B->dval, 1 );
            magma_tally3_index_setvector( A.num_rows * rowlength, A.col, 1, B->dcol, 1 );
            magma_tally3_index_setvector( A.num_rows, A.row, 1, B->drow, 1 );
        }
        //SELLP-type
        if ( A.storage_type == Magma_tally3_SELLP ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_DEV;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            B->blocksize = A.blocksize;
            B->numblocks = A.numblocks;
            B->alignment = A.alignment;
            // memory allocation
            CHECK( magma_tally3_zmalloc( &B->dval, A.nnz ));
            CHECK( magma_tally3_index_malloc( &B->dcol, A.nnz ));
            CHECK( magma_tally3_index_malloc( &B->drow, A.numblocks + 1 ));
            // data transfer
            magma_tally3_zsetvector( A.nnz, A.val, 1, B->dval, 1 );
            magma_tally3_index_setvector( A.nnz, A.col, 1, B->dcol, 1 );
            magma_tally3_index_setvector( A.numblocks+1, A.row, 1, B->drow, 1 );
        }
        //BCSR-type
        if ( A.storage_type == Magma_tally3_BCSR ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_DEV;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            B->blocksize = A.blocksize;
            B->numblocks = A.numblocks;
            magma_tally3_int_t size_b = A.blocksize;
            //magma_tally3_int_t c_blocks = ceil( (float)A.num_cols / (float)size_b );
                    // max number of blocks per row
            magma_tally3_int_t r_blocks = ceil( (float)A.num_rows / (float)size_b );
                    // max number of blocks per column
            // memory allocation
            CHECK( magma_tally3_zmalloc( &B->dval, size_b*size_b*A.numblocks ));
            CHECK( magma_tally3_index_malloc( &B->drow, r_blocks + 1 ));
            CHECK( magma_tally3_index_malloc( &B->dcol, A.numblocks ));
            // data transfer
            magma_tally3_zsetvector( size_b*size_b*A.numblocks, A.val, 1, B->dval, 1 );
            magma_tally3_index_setvector( r_blocks+1, A.row, 1, B->drow, 1 );
            magma_tally3_index_setvector( A.numblocks, A.col, 1, B->dcol, 1 );
        }
        //DENSE-type
        if ( A.storage_type == Magma_tally3_DENSE ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_DEV;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            B->ld = A.ld;
            // memory allocation
            CHECK( magma_tally3_zmalloc( &B->dval, A.num_rows*A.num_cols ));
            // data transfer
            magma_tally3_zsetvector( A.num_rows*A.num_cols, A.val, 1, B->dval, 1 );
        }
    }

    // second case: copy matrix from host to host
    if ( src == Magma_tally3_CPU && dst == Magma_tally3_CPU ) {
        //CSR-type
        if ( A.storage_type == Magma_tally3_CSR || A.storage_type == Magma_tally3_CSC
                                        || A.storage_type == Magma_tally3_CSRD
                                        || A.storage_type == Magma_tally3_CSRL
                                        || A.storage_type == Magma_tally3_CSRU ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_CPU;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            // memory allocation
            CHECK( magma_tally3_zmalloc_cpu( &B->val, A.nnz ));
            CHECK( magma_tally3_index_malloc_cpu( &B->row, A.num_rows+1 ));
            CHECK( magma_tally3_index_malloc_cpu( &B->col, A.nnz ));
            // data transfer
            for( magma_tally3_int_t i=0; i<A.nnz; i++ ) {
                B->val[i] = A.val[i];
                B->col[i] = A.col[i];
            }
            for( magma_tally3_int_t i=0; i<A.num_rows+1; i++ ) {
                B->row[i] = A.row[i];
            }
        }
        //CSRCOO-type
        if ( A.storage_type == Magma_tally3_CSRCOO ) {
            // fill in information for B
            *B = A;
            B->memory_location = Magma_tally3_CPU;
            // memory allocation
            CHECK( magma_tally3_zmalloc_cpu( &B->val, A.nnz ));
            CHECK( magma_tally3_index_malloc_cpu( &B->row, A.num_rows+1 ));
            CHECK( magma_tally3_index_malloc_cpu( &B->col, A.nnz ));
            CHECK( magma_tally3_index_malloc_cpu( &B->rowidx, A.nnz ));
            // data transfer
            for( magma_tally3_int_t i=0; i<A.nnz; i++ ) {
                B->val[i] = A.val[i];
                B->col[i] = A.col[i];
                B->rowidx[i] = A.rowidx[i];
            }
            for( magma_tally3_int_t i=0; i<A.num_rows+1; i++ ) {
                B->row[i] = A.row[i];
            }
        }
        //ELL/ELLPACKT-type
        if ( A.storage_type == Magma_tally3_ELLPACKT || A.storage_type == Magma_tally3_ELL ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_CPU;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            // memory allocation
            CHECK( magma_tally3_zmalloc_cpu( &B->val, A.num_rows*A.max_nnz_row ));
            CHECK( magma_tally3_index_malloc_cpu( &B->col, A.num_rows*A.max_nnz_row ));
            // data transfer
            for( magma_tally3_int_t i=0; i<A.num_rows*A.max_nnz_row; i++ ) {
                B->val[i] = A.val[i];
                B->col[i] = A.col[i];
            }
        }
        //ELLD-type
        if ( A.storage_type == Magma_tally3_ELLD ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_CPU;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            // memory allocation
            CHECK( magma_tally3_zmalloc_cpu( &B->val, A.num_rows*A.max_nnz_row ));
            CHECK( magma_tally3_index_malloc_cpu( &B->col, A.num_rows*A.max_nnz_row ));
            // data transfer
            for( magma_tally3_int_t i=0; i<A.num_rows*A.max_nnz_row; i++ ) {
                B->val[i] = A.val[i];
                B->col[i] = A.col[i];
            }
        }
        //ELLRT-type
        if ( A.storage_type == Magma_tally3_ELLRT ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_CPU;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            B->blocksize = A.blocksize;
            B->alignment = A.alignment;
            int threads_per_row = A.alignment;
            int rowlength = ( (int)((A.max_nnz_row+threads_per_row-1)
                                    /threads_per_row) ) * threads_per_row;
            // memory allocation
            CHECK( magma_tally3_zmalloc_cpu( &B->val, rowlength*A.num_rows ));
            CHECK( magma_tally3_index_malloc_cpu( &B->row, A.num_rows ));
            CHECK( magma_tally3_index_malloc_cpu( &B->col, rowlength*A.num_rows ));
            // data transfer
            for( magma_tally3_int_t i=0; i<A.num_rows*rowlength; i++ ) {
                B->val[i] = A.val[i];
                B->col[i] = A.col[i];
            }
            for( magma_tally3_int_t i=0; i<A.num_rows; i++ ) {
                B->row[i] = A.row[i];
            }
        }
        //SELLP-type
        if (  A.storage_type == Magma_tally3_SELLP ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_CPU;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            B->blocksize = A.blocksize;
            B->alignment = A.alignment;
            B->numblocks = A.numblocks;

            // memory allocation
            CHECK( magma_tally3_zmalloc_cpu( &B->val, A.nnz ));
            CHECK( magma_tally3_index_malloc_cpu( &B->col, A.nnz ));
            CHECK( magma_tally3_index_malloc_cpu( &B->row, A.numblocks+1 ));
            // data transfer
            for( magma_tally3_int_t i=0; i<A.nnz; i++ ) {
                B->val[i] = A.val[i];
                B->col[i] = A.col[i];
            }
            for( magma_tally3_int_t i=0; i<A.numblocks+1; i++ ) {
                B->row[i] = A.row[i];
            }
        }
        //DENSE-type
        if ( A.storage_type == Magma_tally3_DENSE ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_CPU;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            B->ld = A.ld;
            // memory allocation
            CHECK( magma_tally3_zmalloc_cpu( &B->val, A.num_rows*A.num_cols ));
            // data transfer
            for( magma_tally3_int_t i=0; i<A.num_rows*A.num_cols; i++ ) {
                B->val[i] = A.val[i];
            }
        }
    }

    // third case: copy matrix from device to host
    if ( src == Magma_tally3_DEV && dst == Magma_tally3_CPU ) {
        //CSR-type
        if ( A.storage_type == Magma_tally3_CSR || A.storage_type == Magma_tally3_CSC
                                        || A.storage_type == Magma_tally3_CSRD
                                        || A.storage_type == Magma_tally3_CSRL
                                        || A.storage_type == Magma_tally3_CSRU ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_CPU;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            // memory allocation
            CHECK( magma_tally3_zmalloc_cpu( &B->val, A.nnz ));
            CHECK( magma_tally3_index_malloc_cpu( &B->row, A.num_rows+1 ));
            CHECK( magma_tally3_index_malloc_cpu( &B->col, A.nnz ));
            // data transfer
            magma_tally3_zgetvector( A.nnz, A.dval, 1, B->val, 1 );
            magma_tally3_index_getvector( A.num_rows+1, A.drow, 1, B->row, 1 );
            magma_tally3_index_getvector( A.nnz, A.dcol, 1, B->col, 1 );
        }
        //CSRCOO-type
        if ( A.storage_type == Magma_tally3_CSRCOO ) {
            // fill in information for B
            *B = A;
            B->memory_location = Magma_tally3_CPU;
            // memory allocation
            CHECK( magma_tally3_zmalloc_cpu( &B->val, A.nnz ));
            CHECK( magma_tally3_index_malloc_cpu( &B->row, A.num_rows+1 ));
            CHECK( magma_tally3_index_malloc_cpu( &B->col, A.nnz ));
            CHECK( magma_tally3_index_malloc_cpu( &B->rowidx, A.nnz ));
            // data transfer
            magma_tally3_zgetvector( A.nnz, A.dval, 1, B->val, 1 );
            magma_tally3_index_getvector( A.num_rows+1, A.drow, 1, B->row, 1 );
            magma_tally3_index_getvector( A.nnz, A.dcol, 1, B->col, 1 );
            magma_tally3_index_getvector( A.nnz, A.drowidx, 1, B->rowidx, 1 );
        }
        //ELL/ELLPACKT-type
        if ( A.storage_type == Magma_tally3_ELL || A.storage_type == Magma_tally3_ELLPACKT ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_CPU;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            // memory allocation
            CHECK( magma_tally3_zmalloc_cpu( &B->val, A.num_rows*A.max_nnz_row ));
            CHECK( magma_tally3_index_malloc_cpu( &B->col, A.num_rows*A.max_nnz_row ));
            // data transfer
            magma_tally3_zgetvector( A.num_rows*A.max_nnz_row, A.dval, 1, B->val, 1 );
            magma_tally3_index_getvector( A.num_rows*A.max_nnz_row, A.dcol, 1, B->col, 1 );
        }
        //ELLD-type
        if (  A.storage_type == Magma_tally3_ELLD ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_CPU;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            // memory allocation
            CHECK( magma_tally3_zmalloc_cpu( &B->val, A.num_rows*A.max_nnz_row ));
            CHECK( magma_tally3_index_malloc_cpu( &B->col, A.num_rows*A.max_nnz_row ));
            // data transfer
            magma_tally3_zgetvector( A.num_rows*A.max_nnz_row, A.dval, 1, B->val, 1 );
            magma_tally3_index_getvector( A.num_rows*A.max_nnz_row, A.dcol, 1, B->col, 1 );
        }
        //ELLRT-type
        if ( A.storage_type == Magma_tally3_ELLRT ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_CPU;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            B->blocksize = A.blocksize;
            B->alignment = A.alignment;
            int threads_per_row = A.alignment;
            // memory allocation
            int rowlength = ( (int)((A.max_nnz_row+threads_per_row-1)
                                /threads_per_row) ) * threads_per_row;
            CHECK( magma_tally3_zmalloc_cpu( &B->val, rowlength*A.num_rows ));
            CHECK( magma_tally3_index_malloc_cpu( &B->row, A.num_rows ));
            CHECK( magma_tally3_index_malloc_cpu( &B->col, rowlength*A.num_rows ));
            // data transfer
            magma_tally3_zgetvector( A.num_rows*rowlength, A.dval, 1, B->val, 1 );
            magma_tally3_index_getvector( A.num_rows*rowlength, A.dcol, 1, B->col, 1 );
            magma_tally3_index_getvector( A.num_rows, A.drow, 1, B->row, 1 );
        }
        //SELLP-type
        if ( A.storage_type == Magma_tally3_SELLP ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_CPU;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            B->alignment = A.alignment;
            B->blocksize = A.blocksize;
            B->numblocks = A.numblocks;

            // memory allocation
            CHECK( magma_tally3_zmalloc_cpu( &B->val, A.nnz ));
            CHECK( magma_tally3_index_malloc_cpu( &B->col, A.nnz ));
            CHECK( magma_tally3_index_malloc_cpu( &B->row, A.numblocks+1 ));
            // data transfer
            magma_tally3_zgetvector( A.nnz, A.dval, 1, B->val, 1 );
            magma_tally3_index_getvector( A.nnz, A.dcol, 1, B->col, 1 );
            magma_tally3_index_getvector( A.numblocks+1, A.drow, 1, B->row, 1 );
        }
        //BCSR-type
        if ( A.storage_type == Magma_tally3_BCSR ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_CPU;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            B->blocksize = A.blocksize;
            B->numblocks = A.numblocks;
            magma_tally3_int_t size_b = A.blocksize;
            //magma_tally3_int_t c_blocks = ceil( (float)A.num_cols / (float)size_b );
                    // max number of blocks per row
            magma_tally3_int_t r_blocks = ceil( (float)A.num_rows / (float)size_b );
                    // max number of blocks per column
            // memory allocation
            CHECK( magma_tally3_zmalloc_cpu( &B->val, A.numblocks*A.blocksize*A.blocksize ));
            CHECK( magma_tally3_index_malloc_cpu( &B->row, r_blocks+1 ));
            CHECK( magma_tally3_index_malloc_cpu( &B->col, A.numblocks ));
            // data transfer
            magma_tally3_zgetvector( A.numblocks * A.blocksize * A.blocksize, A.dval, 1, B->val, 1 );
            magma_tally3_index_getvector( r_blocks+1, A.drow, 1, B->row, 1 );
            magma_tally3_index_getvector( A.numblocks, A.dcol, 1, B->col, 1 );
        }
        //DENSE-type
        if ( A.storage_type == Magma_tally3_DENSE ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_CPU;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            B->ld = A.ld;
            // memory allocation
            CHECK( magma_tally3_zmalloc_cpu( &B->val, A.num_rows*A.num_cols ));
            // data transfer
            magma_tally3_zgetvector( A.num_rows*A.num_cols, A.dval, 1, B->val, 1 );
        }
    }

    // fourth case: copy matrix from device to device
    if ( src == Magma_tally3_DEV && dst == Magma_tally3_DEV ) {
        //CSR-type
        if ( A.storage_type == Magma_tally3_CSR || A.storage_type == Magma_tally3_CSC
                                        || A.storage_type == Magma_tally3_CSRD
                                        || A.storage_type == Magma_tally3_CSRL
                                        || A.storage_type == Magma_tally3_CSRU ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_DEV;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            // memory allocation
            CHECK( magma_tally3_zmalloc( &B->dval, A.nnz ));
            CHECK( magma_tally3_index_malloc( &B->drow, A.num_rows + 1 ));
            CHECK( magma_tally3_index_malloc( &B->dcol, A.nnz ));
            // data transfer
            magma_tally3_zcopyvector( A.nnz, A.dval, 1, B->dval, 1 );
            magma_tally3_index_copyvector( (A.num_rows+1), A.drow, 1, B->drow, 1 );
            magma_tally3_index_copyvector( A.nnz, A.dcol, 1, B->dcol, 1 );
        }
        //CSRCOO-type
        if ( A.storage_type == Magma_tally3_CSRCOO ) {
            // fill in information for B
            *B = A;
            B->memory_location = Magma_tally3_DEV;
            // memory allocation
            CHECK( magma_tally3_zmalloc( &B->dval, A.nnz ));
            CHECK( magma_tally3_index_malloc( &B->drow, A.num_rows + 1 ));
            CHECK( magma_tally3_index_malloc( &B->dcol, A.nnz ));
            CHECK( magma_tally3_index_malloc( &B->drowidx, A.nnz ));
            // data transfer
            magma_tally3_zcopyvector( A.nnz, A.dval, 1, B->dval, 1 );
            magma_tally3_index_copyvector( (A.num_rows+1), A.drow, 1, B->drow, 1 );
            magma_tally3_index_copyvector( A.nnz, A.dcol, 1, B->dcol, 1 );
            magma_tally3_index_copyvector( A.nnz, A.drowidx, 1, B->drowidx, 1 );
        }
        //ELL/ELLPACKT-type
        if ( A.storage_type == Magma_tally3_ELL || A.storage_type == Magma_tally3_ELLPACKT ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_DEV;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            // memory allocation
            CHECK( magma_tally3_zmalloc( &B->dval, A.num_rows * A.max_nnz_row ));
            CHECK( magma_tally3_index_malloc( &B->dcol, A.num_rows * A.max_nnz_row ));
            // data transfer
            magma_tally3_zcopyvector( A.num_rows*A.max_nnz_row, A.dval, 1, B->dval, 1 );
            magma_tally3_index_copyvector( A.num_rows*A.max_nnz_row, A.dcol, 1, B->dcol, 1 );
        }
        //ELLD-type
        if ( A.storage_type == Magma_tally3_ELLD ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_DEV;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            // memory allocation
            CHECK( magma_tally3_zmalloc( &B->dval, A.num_rows * A.max_nnz_row ));
            CHECK( magma_tally3_index_malloc( &B->dcol, A.num_rows * A.max_nnz_row ));
            // data transfer
            magma_tally3_zcopyvector( A.num_rows*A.max_nnz_row, A.dval, 1, B->dval, 1 );
            magma_tally3_index_copyvector( A.num_rows*A.max_nnz_row, A.dcol, 1, B->dcol, 1 );
        }
        //ELLRT-type
        if ( A.storage_type == Magma_tally3_ELLRT ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_DEV;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            B->blocksize = A.blocksize;
            B->alignment = A.alignment;
            int threads_per_row = A.alignment;
            int rowlength = ( (int)((A.max_nnz_row+threads_per_row-1)
                                    /threads_per_row) ) * threads_per_row;
            // memory allocation
            CHECK( magma_tally3_zmalloc( &B->dval, A.num_rows * rowlength ));
            CHECK( magma_tally3_index_malloc( &B->dcol, A.num_rows * rowlength ));
            CHECK( magma_tally3_index_malloc( &B->drow, A.num_rows ));
            // data transfer
            magma_tally3_zcopyvector( A.num_rows * rowlength, A.dval, 1, B->dval, 1 );
            magma_tally3_index_copyvector( A.num_rows * rowlength, A.dcol, 1, B->dcol, 1 );
            magma_tally3_index_copyvector( A.num_rows, A.drow, 1, B->drow, 1 );
        }
        //SELLP-type
        if ( A.storage_type == Magma_tally3_SELLP ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_DEV;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            B->blocksize = A.blocksize;
            B->numblocks = A.numblocks;
            // memory allocation
            CHECK( magma_tally3_zmalloc( &B->dval, A.nnz ));
            CHECK( magma_tally3_index_malloc( &B->dcol, A.nnz ));
            CHECK( magma_tally3_index_malloc( &B->drow, A.numblocks + 1 ));
            // data transfer
            magma_tally3_zcopyvector( A.nnz, A.dval, 1, B->dval, 1 );
            magma_tally3_index_copyvector( A.nnz,         A.dcol, 1, B->dcol, 1 );
            magma_tally3_index_copyvector( A.numblocks+1, A.drow, 1, B->drow, 1 );
        }
        //BCSR-type
        if ( A.storage_type == Magma_tally3_BCSR ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_DEV;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            B->blocksize = A.blocksize;
            B->numblocks = A.numblocks;
            magma_tally3_int_t size_b = A.blocksize;
            //magma_tally3_int_t c_blocks = ceil( (float)A.num_cols / (float)size_b );
                    // max number of blocks per row
            magma_tally3_int_t r_blocks = ceil( (float)A.num_rows / (float)size_b );
                    // max number of blocks per column
            // memory allocation
            CHECK( magma_tally3_zmalloc( &B->dval, size_b*size_b*A.numblocks ));
            CHECK( magma_tally3_index_malloc( &B->drow, r_blocks + 1 ));
            CHECK( magma_tally3_index_malloc( &B->dcol, A.numblocks ));
            // data transfer
            magma_tally3_zcopyvector( size_b*size_b*A.numblocks, A.dval, 1, B->dval, 1 );
            magma_tally3_index_copyvector( (r_blocks+1), A.drow, 1, B->drow, 1 );
            magma_tally3_index_copyvector( A.numblocks, A.dcol, 1, B->dcol, 1 );
        }
        //DENSE-type
        if ( A.storage_type == Magma_tally3_DENSE ) {
            // fill in information for B
            B->storage_type = A.storage_type;
            B->diagorder_type = A.diagorder_type;
            B->memory_location = Magma_tally3_DEV;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->fill_mode = A.fill_mode;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            B->ld = A.ld;
            // memory allocation
            CHECK( magma_tally3_zmalloc( &B->dval, A.num_rows*A.num_cols ));
            // data transfer
            magma_tally3_zcopyvector( A.num_rows*A.num_cols, A.dval, 1, B->dval, 1 );
        }
    }
    
    
cleanup:
    if( info != 0 ){
        magma_tally3_zmfree( B, queue );
    }
    magma_tally3blasSetKernelStream( orig_queue );
    return info;
}
