/*
    -- MAGMA_tally4 (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @generated from magma_tally4_zmconverter.cpp normal z -> s, Tue May  5 14:02:10 2015
       @author Hartwig Anzt
*/
#include "common_magma_tally4sparse.h"


/**
    Purpose
    -------

    Helper function to compress CSR containing zero-entries.


    Arguments
    ---------

    @param[in]
    val         float**
                input val pointer to compress

    @param[in]
    row         magma_tally4_int_t**
                input row pointer to modify

    @param[in]
    col         magma_tally4_int_t**
                input col pointer to compress

    @param[in]
    valn        float**
                output val pointer

    @param[out]
    rown        magma_tally4_int_t**
                output row pointer

    @param[out]
    coln        magma_tally4_int_t**
                output col pointer

    @param[out]
    n           magma_tally4_int_t*
                number of rows in matrix

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_saux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_s_csr_compressor(
    float ** val,
    magma_tally4_index_t ** row,
    magma_tally4_index_t ** col,
    float ** valn,
    magma_tally4_index_t ** rown,
    magma_tally4_index_t ** coln,
    magma_tally4_int_t *n,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    magma_tally4_index_t i,j, nnz_new=0, (*row_nnz)=NULL, nnz_this_row;
    CHECK( magma_tally4_index_malloc_cpu( &(row_nnz), (*n) ));
    CHECK( magma_tally4_index_malloc_cpu( rown, *n+1 ));
    for( i=0; i<*n; i++ ) {
        (*rown)[i] = nnz_new;
        nnz_this_row = 0;
        for( j=(*row)[i]; j<(*row)[i+1]; j++ ) {
            if ( MAGMA_tally4_S_REAL((*val)[j]) != 0 ) {
                nnz_new++;
                nnz_this_row++;
            }
        }
        row_nnz[i] = nnz_this_row;
    }
    (*rown)[*n] = nnz_new;

    CHECK( magma_tally4_smalloc_cpu( valn, nnz_new ));
    CHECK( magma_tally4_index_malloc_cpu( coln, nnz_new ));

    nnz_new = 0;
    for( i=0; i<*n; i++ ) {
        for( j=(*row)[i]; j<(*row)[i+1]; j++ ) {
            if ( MAGMA_tally4_S_REAL((*val)[j]) != 0 ) {
                (*valn)[nnz_new]= (*val)[j];
                (*coln)[nnz_new]= (*col)[j];
                nnz_new++;
            }
        }
    }


cleanup:
    if ( info != 0 ) {
        magma_tally4_free_cpu( valn );
        magma_tally4_free_cpu( coln );
        magma_tally4_free_cpu( rown );
    }
    magma_tally4_free_cpu( row_nnz );
    row_nnz = NULL;
    return info;
}






/**
    Purpose
    -------

    Converter between different sparse storage formats.

    Arguments
    ---------

    @param[in]
    A           magma_tally4_s_matrix
                sparse matrix A

    @param[out]
    B           magma_tally4_s_matrix*
                copy of A in new format

    @param[in]
    old_format  magma_tally4_storage_t
                original storage format

    @param[in]
    new_format  magma_tally4_storage_t
                new storage format

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_saux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_smconvert(
    magma_tally4_s_matrix A,
    magma_tally4_s_matrix *B,
    magma_tally4_storage_t old_format,
    magma_tally4_storage_t new_format,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;

    magma_tally4_index_t *length=NULL;
    
    magma_tally4_s_matrix hA={Magma_tally4_CSR}, hB={Magma_tally4_CSR};
    magma_tally4_s_matrix A_d={Magma_tally4_CSR}, B_d={Magma_tally4_CSR};
    magma_tally4_index_t *row_tmp=NULL, *col_tmp=NULL;
    float *val_tmp = NULL;
    magma_tally4_index_t *row_tmp2=NULL, *col_tmp2=NULL;
    float *val_tmp2 = NULL;
    float *transpose=NULL;
    magma_tally4_index_t intnnz, *nnz_per_row=NULL;
    
    cusparseHandle_t cusparseHandle = 0;
    cusparseMatDescr_t descr = 0;
    
    // set queue for old dense routines
    magma_tally4_queue_t orig_queue=NULL;
    magma_tally4blasGetKernelStream( &orig_queue );

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

    float zero = MAGMA_tally4_S_MAKE( 0.0, 0.0 );

    // check whether matrix on CPU
    if ( A.memory_location == Magma_tally4_CPU ) {
        
        // CSR to anything
        if ( old_format == Magma_tally4_CSR ) {
            
            // CSR to CSR
            if ( new_format == Magma_tally4_CSR ) {
                // fill in information for B
                B->storage_type = Magma_tally4_CSR;
                B->memory_location = A.memory_location;
                B->fill_mode = A.fill_mode;
                B->num_rows = A.num_rows;
                B->num_cols = A.num_cols;
                B->nnz = A.nnz;
                B->max_nnz_row = A.max_nnz_row;
                B->diameter = A.diameter;
    
                CHECK( magma_tally4_smalloc_cpu( &B->val, A.nnz ));
                CHECK( magma_tally4_index_malloc_cpu( &B->row, A.num_rows+1 ));
                CHECK( magma_tally4_index_malloc_cpu( &B->col, A.nnz ));
                
                for( magma_tally4_int_t i=0; i < A.nnz; i++) {
                    B->val[i] = A.val[i];
                    B->col[i] = A.col[i];
                }
                for( magma_tally4_int_t i=0; i < A.num_rows+1; i++) {
                    B->row[i] = A.row[i];
                }
            }
            
            // CSR to CSRL
            else if ( new_format == Magma_tally4_CSRL ) {
                // fill in information for B
                B->storage_type = Magma_tally4_CSR;
                B->memory_location = A.memory_location;
                B->fill_mode = Magma_tally4_LOWER;
                B->num_rows = A.num_rows;
                B->num_cols = A.num_cols;
                B->diameter = A.diameter;
    
                magma_tally4_int_t numzeros=0;
                for( magma_tally4_int_t i=0; i < A.num_rows; i++) {
                    for( magma_tally4_int_t j=A.row[i]; j < A.row[i+1]; j++) {
                        if ( A.col[j] <= i) {
                            numzeros++;
                        }
                    }
                }
                B->nnz = numzeros;
                CHECK( magma_tally4_smalloc_cpu( &B->val, numzeros ));
                CHECK( magma_tally4_index_malloc_cpu( &B->row, A.num_rows+1 ));
                CHECK( magma_tally4_index_malloc_cpu( &B->col, numzeros ));
                
                numzeros=0;
                for( magma_tally4_int_t i=0; i < A.num_rows; i++) {
                    B->row[i]=numzeros;
                    for( magma_tally4_int_t j=A.row[i]; j < A.row[i+1]; j++) {
                        if ( A.col[j] < i) {
                            B->val[numzeros] = A.val[j];
                            B->col[numzeros] = A.col[j];
                            numzeros++;
                        }
                        else if ( A.col[j] == i &&
                                        B->diagorder_type == Magma_tally4_UNITY) {
                            B->val[numzeros] = MAGMA_tally4_S_MAKE(1.0, 0.0);
                            B->col[numzeros] = A.col[j];
                            numzeros++;
                        }
                        else if ( A.col[j] == i ) {
                            B->val[numzeros] = A.val[j];
                            B->col[numzeros] = A.col[j];
                            numzeros++;
                        }
                    }
                }
                B->row[B->num_rows] = numzeros;
            }
        
            // CSR to CSRU
            else if (  new_format == Magma_tally4_CSRU ) {
                // fill in information for B
                *B = A;
                B->fill_mode = Magma_tally4_UPPER;
                magma_tally4_int_t numzeros=0;
                for( magma_tally4_int_t i=0; i < A.num_rows; i++) {
                    for( magma_tally4_int_t j=A.row[i]; j < A.row[i+1]; j++) {
                        if ( A.col[j] >= i) {
                            numzeros++;
                        }
                    }
                }
                B->nnz = numzeros;
                CHECK( magma_tally4_smalloc_cpu( &B->val, numzeros ));
                CHECK( magma_tally4_index_malloc_cpu( &B->row, A.num_rows+1 ));
                CHECK( magma_tally4_index_malloc_cpu( &B->col, numzeros ));
                
                numzeros=0;
                for( magma_tally4_int_t i=0; i < A.num_rows; i++) {
                    B->row[i]=numzeros;
                    for( magma_tally4_int_t j=A.row[i]; j < A.row[i+1]; j++) {
                        if ( A.col[j] >= i) {
                            B->val[numzeros] = A.val[j];
                            B->col[numzeros] = A.col[j];
                            numzeros++;
                        }
                    }
                }
                B->row[B->num_rows] = numzeros;
            }
            
            // CSR to CSRD (diagonal elements first)
            else if ( new_format == Magma_tally4_CSRD ) {
                // fill in information for B
                B->storage_type = Magma_tally4_CSRD;
                B->memory_location = A.memory_location;
                B->fill_mode = A.fill_mode;
                B->num_rows = A.num_rows;
                B->num_cols = A.num_cols;
                B->nnz = A.nnz;
                B->max_nnz_row = A.max_nnz_row;
                B->diameter = A.diameter;
    
                CHECK( magma_tally4_smalloc_cpu( &B->val, A.nnz ));
                CHECK( magma_tally4_index_malloc_cpu( &B->row, A.num_rows+1 ));
                CHECK( magma_tally4_index_malloc_cpu( &B->col, A.nnz ));
                
                for(magma_tally4_int_t i=0; i < A.num_rows; i++) {
                    magma_tally4_int_t count = 1;
                    for(magma_tally4_int_t j=A.row[i]; j < A.row[i+1]; j++) {
                        if ( A.col[j] == i ) {
                            B->col[A.row[i]] = A.col[j];
                            B->val[A.row[i]] = A.val[j];
                        } else {
                            B->col[A.row[i]+count] = A.col[j];
                            B->val[A.row[i]+count] = A.val[j];
                            count++;
                        }
                    }
                }
                for( magma_tally4_int_t i=0; i < A.num_rows+1; i++) {
                    B->row[i] = A.row[i];
                }
            }
                        
            // CSR to COO
            else if ( new_format == Magma_tally4_COO ) {
                CHECK( magma_tally4_smconvert( A, B, Magma_tally4_CSR, Magma_tally4_CSR, queue ));
                B->storage_type = Magma_tally4_COO;
        
                magma_tally4_free_cpu( B->row );
                CHECK( magma_tally4_index_malloc_cpu( &B->row, A.nnz ));
                
                for(magma_tally4_int_t i=0; i < A.num_rows; i++) {
                    for(magma_tally4_int_t j=A.row[i]; j < A.row[i+1]; j++) {
                            B->row[j] = i;
                    }
                }
            }
                    
            // CSR to CSRCOO
            else if ( new_format == Magma_tally4_CSRCOO ) {
                CHECK( magma_tally4_smconvert( A, B, Magma_tally4_CSR, Magma_tally4_CSR, queue ));
                B->storage_type = Magma_tally4_CSRCOO;
    
                CHECK( magma_tally4_index_malloc_cpu( &B->rowidx, A.nnz ));
                
                for(magma_tally4_int_t i=0; i < A.num_rows; i++) {
                    for(magma_tally4_int_t j=A.row[i]; j < A.row[i+1]; j++) {
                            B->rowidx[j] = i;
                    }
                }
            }
   
            // CSR to ELLPACKT (using row-major storage)
            else if (  new_format == Magma_tally4_ELLPACKT ) {
                // fill in information for B
                B->storage_type = Magma_tally4_ELLPACKT;
                B->memory_location = A.memory_location;
                B->fill_mode = A.fill_mode;
                B->num_rows = A.num_rows;
                B->num_cols = A.num_cols;
                B->nnz = A.nnz;
                B->max_nnz_row = A.max_nnz_row;
                B->diameter = A.diameter;
                // conversion
                magma_tally4_index_t i, j, maxrowlength=0;
                CHECK( magma_tally4_index_malloc_cpu( &length, A.num_rows));
    
                for( i=0; i < A.num_rows; i++ ) {
                    length[i] = A.row[i+1]-A.row[i];
                    if (length[i] > maxrowlength)
                        maxrowlength = length[i];
                }
                //printf( "Conversion to ELLPACK with %d elements per row: ",
                                                                // maxrowlength );
                //fflush(stdout);
                CHECK( magma_tally4_smalloc_cpu( &B->val, maxrowlength*A.num_rows ));
                CHECK( magma_tally4_index_malloc_cpu( &B->col, maxrowlength*A.num_rows ));
                
                
                for( magma_tally4_int_t i=0; i<(maxrowlength*A.num_rows); i++) {
                    B->val[i] = MAGMA_tally4_S_MAKE(0., 0.);
                    B->col[i] =  -1;
                }
                for( i=0; i < A.num_rows; i++ ) {
                    magma_tally4_int_t offset = 0;
                    for( j=A.row[i]; j < A.row[i+1]; j++ ) {
                        B->val[i*maxrowlength+offset] = A.val[j];
                        B->col[i*maxrowlength+offset] = A.col[j];
                        offset++;
                    }
                }
                B->max_nnz_row = maxrowlength;
            }
            
            // CSR to ELL
            else if ( new_format == Magma_tally4_ELL ) {
                // fill in information for B
                B->storage_type = Magma_tally4_ELL;
                B->memory_location = A.memory_location;
                B->fill_mode = A.fill_mode;
                B->num_rows = A.num_rows;
                B->num_cols = A.num_cols;
                B->nnz = A.nnz;
                B->max_nnz_row = A.max_nnz_row;
                B->diameter = A.diameter;
    
                // conversion
                magma_tally4_index_t i, j, maxrowlength=0;
                CHECK( magma_tally4_index_malloc_cpu( &length, A.num_rows));

                for( i=0; i < A.num_rows; i++ ) {
                    length[i] = A.row[i+1]-A.row[i];
                    if (length[i] > maxrowlength)
                        maxrowlength = length[i];
                }
                //printf( "Conversion to ELL with %d elements per row: ",
                                                               // maxrowlength );
                //fflush(stdout);
                CHECK( magma_tally4_smalloc_cpu( &B->val, maxrowlength*A.num_rows ));
                CHECK( magma_tally4_index_malloc_cpu( &B->col, maxrowlength*A.num_rows ));
                
                for( magma_tally4_int_t i=0; i<(maxrowlength*A.num_rows); i++) {
                    B->val[i] = MAGMA_tally4_S_MAKE(0., 0.);
                    B->col[i] =  -1;
                }
    
                for( i=0; i < A.num_rows; i++ ) {
                    magma_tally4_int_t offset = 0;
                    for( j=A.row[i]; j < A.row[i+1]; j++ ) {
                        B->val[offset*A.num_rows+i] = A.val[j];
                        B->col[offset*A.num_rows+i] = A.col[j];
                        offset++;
                    }
                }
                B->max_nnz_row = maxrowlength;
                //printf( "done\n" );
            }
            
            // CSR to ELLD (ELLPACK with diagonal element first)
            else if ( new_format == Magma_tally4_ELLD ) {
                // fill in information for B
                B->storage_type = Magma_tally4_ELLD;
                B->memory_location = A.memory_location;
                B->fill_mode = A.fill_mode;
                B->num_rows = A.num_rows;
                B->num_cols = A.num_cols;
                B->nnz = A.nnz;
                B->max_nnz_row = A.max_nnz_row;
                B->diameter = A.diameter;
    
                // conversion
                magma_tally4_index_t i, j, maxrowlength=0;
                CHECK( magma_tally4_index_malloc_cpu( &length, A.num_rows));

                for( i=0; i < A.num_rows; i++ ) {
                    length[i] = A.row[i+1]-A.row[i];
                    if (length[i] > maxrowlength)
                        maxrowlength = length[i];
                }
                //printf( "Conversion to ELL with %d elements per row: ",
                                                               // maxrowlength );
                //fflush(stdout);
                CHECK( magma_tally4_smalloc_cpu( &B->val, maxrowlength*A.num_rows ));
                CHECK( magma_tally4_index_malloc_cpu( &B->col, maxrowlength*A.num_rows ));
                
                
                for( magma_tally4_int_t i=0; i<(maxrowlength*A.num_rows); i++) {
                    B->val[i] = MAGMA_tally4_S_MAKE(0., 0.);
                    B->col[i] =  -1;
                }
    
                for( i=0; i < A.num_rows; i++ ) {
                    magma_tally4_int_t offset = 1;
                    for( j=A.row[i]; j < A.row[i+1]; j++ ) {
                        if ( A.col[j] == i ) { // diagonal case
                            B->val[i*maxrowlength] = A.val[j];
                            B->col[i*maxrowlength] = A.col[j];
                        } else {
                            B->val[i*maxrowlength+offset] = A.val[j];
                            B->col[i*maxrowlength+offset] = A.col[j];
                            offset++;
                        }
                    }
                }
                B->max_nnz_row = maxrowlength;
            }
            
            // CSR to ELLRT (also ELLPACKRT)
            else if (  new_format == Magma_tally4_ELLRT ) {
                // fill in information for B
                B->storage_type = Magma_tally4_ELLRT;
                B->memory_location = A.memory_location;
                B->fill_mode = A.fill_mode;
                B->num_rows = A.num_rows;
                B->num_cols = A.num_cols;
                B->nnz = A.nnz;
                B->max_nnz_row = A.max_nnz_row;
                B->diameter = A.diameter;
    
                // conversion
                magma_tally4_index_t i, j, maxrowlength=0;
                CHECK( magma_tally4_index_malloc_cpu( &length, A.num_rows));

                for( i=0; i < A.num_rows; i++ ) {
                    length[i] = A.row[i+1]-A.row[i];
                    if (length[i] > maxrowlength)
                        maxrowlength = length[i];
                }

                //printf( "Conversion to ELLRT with %d elements per row: ",
                //                                                   maxrowlength );
    
                magma_tally4_int_t threads_per_row = B->alignment;
                magma_tally4_int_t rowlength = magma_tally4_roundup( maxrowlength, threads_per_row );
    
                CHECK( magma_tally4_smalloc_cpu( &B->val, rowlength*A.num_rows ));
                CHECK( magma_tally4_index_malloc_cpu( &B->col, rowlength*A.num_rows ));
                CHECK( magma_tally4_index_malloc_cpu( &B->row, A.num_rows ));
                
                
                for( magma_tally4_int_t i=0; i < rowlength*A.num_rows; i++) {
                    B->val[i] = MAGMA_tally4_S_MAKE(0., 0.);
                    B->col[i] =  0;
                }
    
                for( i=0; i < A.num_rows; i++ ) {
                    magma_tally4_int_t offset = 0;
                    for( j=A.row[i]; j < A.row[i+1]; j++ ) {
                        B->val[i*rowlength+offset] = A.val[j];
                        B->col[i*rowlength+offset] = A.col[j];
                        offset++;
                    }
                    B->row[i] = A.row[i+1] - A.row[i];
                }
                B->max_nnz_row = maxrowlength;
                //printf( "done\n" );
            }
            
            // CSR to SELLP
            // SELLC is SELLP using alignment 1
            // see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
            // A UNIFIED SPARSE MATRIX DATA FORMAT
            // FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
            // in SELLP we modify SELLC:
            // alignment is posible such that multiple threads can be used for SpMV
            // so the rowlength is padded (SELLP) to a multiple of the alignment
            else if ( new_format == Magma_tally4_SELLP ) {
                // fill in information for B
                B->storage_type = new_format;
                if (B->alignment > 1)
                    B->storage_type = Magma_tally4_SELLP;
                B->memory_location = A.memory_location;
                B->fill_mode = A.fill_mode;
                B->num_rows = A.num_rows;
                B->num_cols = A.num_cols;
                B->diameter = A.diameter;
                B->max_nnz_row = 0;
                magma_tally4_int_t C = B->blocksize;
                magma_tally4_int_t slices = ( A.num_rows+C-1)/(C);
                B->numblocks = slices;
                magma_tally4_int_t alignedlength, alignment = B->alignment;
                // conversion
                magma_tally4_index_t i, j, k, maxrowlength=0;
                CHECK( magma_tally4_index_malloc_cpu( &length, C));
                // B-row points to the start of each slice
                CHECK( magma_tally4_index_malloc_cpu( &B->row, slices+1 ));
                
                
                B->row[0] = 0;
                for( i=0; i < slices; i++ ) {
                    maxrowlength = 0;
                    for(j=0; j < C; j++) {
                        if (i*C+j < A.num_rows) {
                            length[j] = A.row[i*C+j+1]-A.row[i*C+j];
                        }
                        else
                            length[j]=0;
                        if (length[j] > maxrowlength) {
                            maxrowlength = length[j];
                        }
                    }
                    alignedlength = magma_tally4_roundup( maxrowlength, alignment );
                    B->row[i+1] = B->row[i] + alignedlength * C;
                    if ( alignedlength > B->max_nnz_row )
                        B->max_nnz_row = alignedlength;
                }
                B->nnz = B->row[slices];
                //printf( "Conversion to SELLC with %d slices of size %d and"
                //       " %d nonzeros.\n", slices, C, B->nnz );
    
                //fflush(stdout);
                CHECK( magma_tally4_smalloc_cpu( &B->val, B->row[slices] ));
                CHECK( magma_tally4_index_malloc_cpu( &B->col, B->row[slices] ));
                
                // zero everything
                for( i=0; i < B->row[slices]; i++ ) {
                    B->val[ i ] = MAGMA_tally4_S_MAKE(0., 0.);
                    B->col[ i ] =  0;
                }
                // fill in values
                for( i=0; i < slices; i++ ) {
                    for(j=0; j < C; j++) {
                        magma_tally4_int_t line = i*C+j;
                        magma_tally4_int_t offset = 0;
                        if ( line < A.num_rows) {
                            for( k=A.row[line]; k < A.row[line+1]; k++ ) {
                                B->val[ B->row[i] + j +offset*C ] = A.val[k];
                                B->col[ B->row[i] + j +offset*C ] = A.col[k];
                                offset++;
                            }
                        }
                    }
                }
                //B->nnz = A.nnz;
            }
            
            // CSR to DENSE
            else if ( new_format == Magma_tally4_DENSE ) {
                //printf( "Conversion to DENSE: " );
                // fill in information for B
                B->storage_type = Magma_tally4_DENSE;
                B->memory_location = A.memory_location;
                B->fill_mode = A.fill_mode;
                B->num_rows = A.num_rows;
                B->num_cols = A.num_cols;
                B->nnz = A.nnz;
                B->max_nnz_row = A.max_nnz_row;
                B->diameter = A.diameter;
    
                // conversion
                CHECK( magma_tally4_smalloc_cpu( &B->val, A.num_rows*A.num_cols ));
                
                for( magma_tally4_int_t i=0; i<(A.num_rows)*(A.num_cols); i++) {
                    B->val[i] = MAGMA_tally4_S_MAKE(0., 0.);
                }
    
                for(magma_tally4_int_t i=0; i < A.num_rows; i++ ) {
                    for(magma_tally4_int_t j=A.row[i]; j < A.row[i+1]; j++ )
                        B->val[i * (A.num_cols) + A.col[j] ] = A.val[ j ];
                }
    
                //printf( "done\n" );
            }
            
            // CSR to BCSR
            else if ( new_format == Magma_tally4_BCSR ) {
                CHECK( magma_tally4_smtransfer(A, &A_d, Magma_tally4_CPU, Magma_tally4_DEV, queue ) );              
                B_d.blocksize = B->blocksize;                                                  
                CHECK( magma_tally4_smconvert(A_d, &B_d, Magma_tally4_CSR, Magma_tally4_BCSR, queue ) );            
                CHECK( magma_tally4_smtransfer(B_d, B, Magma_tally4_DEV, Magma_tally4_CPU, queue ) );               
            }
            else {
                printf("error: format not supported.\n");
                info = MAGMA_tally4_ERR_NOT_SUPPORTED;
            }
        }
        // anything to CSR
        else if ( new_format == Magma_tally4_CSR ) {
            // CSRU/CSRCSCU to CSR
            if ( old_format == Magma_tally4_CSRU ) {
                CHECK( magma_tally4_smconvert( A, B, Magma_tally4_CSR, Magma_tally4_CSR, queue ));
            }
                    
            // CSRD to CSR (diagonal elements first)
            else if ( old_format == Magma_tally4_CSRD ) {
                // fill in information for B
                B->storage_type = Magma_tally4_CSR;
                B->memory_location = A.memory_location;
                B->fill_mode = A.fill_mode;
                B->num_rows = A.num_rows;
                B->num_cols = A.num_cols;
                B->nnz = A.nnz;
                B->max_nnz_row = A.max_nnz_row;
                B->diameter = A.diameter;

    
                CHECK( magma_tally4_smalloc_cpu( &B->val, A.nnz ));
                CHECK( magma_tally4_index_malloc_cpu( &B->row, A.num_rows+1 ));
                CHECK( magma_tally4_index_malloc_cpu( &B->col, A.nnz ));
                
                for(magma_tally4_int_t i=0; i < A.num_rows; i++) {
                    float diagval = A.val[A.row[i]];
                    magma_tally4_index_t diagcol = A.col[A.row[i]];
                    magma_tally4_int_t smaller = 0;
                    for( magma_tally4_int_t k=A.row[i]; k < A.row[i+1]; k++ ) {
                        if ( (A.col[k] < diagcol) )
                            smaller++;
                    }
                    for( magma_tally4_int_t k=A.row[i]; k < A.row[i]+smaller; k++ ) {
                        B->col[k] = A.col[k+1];
                        B->val[k] = A.val[k+1];
                    }
                    B->col[A.row[i]+smaller] = diagcol;
                    B->val[A.row[i]+smaller] = diagval;
                    for( magma_tally4_int_t k=A.row[i]+smaller+1; k < A.row[i+1]; k++ ) {
                        B->col[k] = A.col[k];
                        B->val[k] = A.val[k];
                    }
                }
                for( magma_tally4_int_t i=0; i < A.num_rows+1; i++) {
                    B->row[i] = A.row[i];
                }
            }
            
            // CSRCOO to CSR
            else if ( old_format == Magma_tally4_CSRCOO ) {
                CHECK( magma_tally4_smconvert( A, B, Magma_tally4_CSR, Magma_tally4_CSR, queue ));
            }
            
            // CSRCSC to CSR
            else if ( old_format == Magma_tally4_COO ) {
               // A.storage_type = Magma_tally4_CSR;
              //  magma_tally4_smconvert( A, B, Magma_tally4_CSR, Magma_tally4_CSR, queue );
            }
                
            // ELL/ELLPACK to CSR
            else if ( old_format == Magma_tally4_ELLPACKT ) {
                //printf( "Conversion to CSR: " );
                // fill in information for B
                B->storage_type = Magma_tally4_CSR;
                B->memory_location = A.memory_location;
                B->fill_mode = A.fill_mode;
                B->num_rows = A.num_rows;
                B->num_cols = A.num_cols;
                B->nnz = A.nnz;
                B->max_nnz_row = A.max_nnz_row;
                B->diameter = A.diameter;
        
                // conversion

                CHECK( magma_tally4_index_malloc_cpu( &row_tmp, A.num_rows+1 ));
                //fill the row-pointer
                for( magma_tally4_int_t i=0; i < A.num_rows+1; i++ )
                    row_tmp[i] = i*A.max_nnz_row;
                //now use AA_ELL, IA_ELL, row_tmp as CSR with some zeros.
                //The CSR compressor removes these
                CHECK( magma_tally4_s_csr_compressor(&A.val, &row_tmp, &A.col,
                           &B->val, &B->row, &B->col, &B->num_rows, queue ));
                B->nnz = B->row[B->num_rows];
            }
            
            // ELL (column-major) to CSR
            else if ( old_format == Magma_tally4_ELL ) {
                //printf( "Conversion to CSR: " );
                //fflush(stdout);
                // fill in information for B
                B->storage_type = Magma_tally4_CSR;
                B->memory_location = A.memory_location;
                B->fill_mode = A.fill_mode;
                B->num_rows = A.num_rows;
                B->num_cols = A.num_cols;
                B->nnz = A.nnz;
                B->max_nnz_row = A.max_nnz_row;
                B->diameter = A.diameter;
    
                // conversion
                CHECK( magma_tally4_smalloc_cpu( &val_tmp, A.num_rows*A.max_nnz_row ));
                CHECK( magma_tally4_index_malloc_cpu( &row_tmp, A.num_rows+1 ));
                CHECK( magma_tally4_index_malloc_cpu( &col_tmp, A.num_rows*A.max_nnz_row ));

                //fill the row-pointer
                for( magma_tally4_int_t i=0; i < A.num_rows+1; i++ )
                    row_tmp[i] = i*A.max_nnz_row;
                //transform RowMajor to ColMajor
                for( magma_tally4_int_t j=0; j < A.max_nnz_row; j++ ) {
                    for( magma_tally4_int_t i=0; i < A.num_rows; i++ ) {
                        col_tmp[i*A.max_nnz_row+j] = A.col[j*A.num_rows+i];
                        val_tmp[i*A.max_nnz_row+j] = A.val[j*A.num_rows+i];
                    }
                }
                //now use AA_ELL, IA_ELL, row_tmp as CSR with some zeros.
                //The CSR compressor removes these
                CHECK( magma_tally4_s_csr_compressor(&val_tmp, &row_tmp, &col_tmp,
                           &B->val, &B->row, &B->col, &B->num_rows, queue ));
    
                B->nnz = B->row[B->num_rows];
            }
            
            // ELLD (ELLPACK with diagonal element first) to CSR
            else if ( old_format == Magma_tally4_ELLD ) {
                //printf( "Conversion to CSR: " );
                //fflush(stdout);
                // fill in information for B
                B->storage_type = Magma_tally4_CSR;
                B->memory_location = A.memory_location;
                B->fill_mode = A.fill_mode;
                B->num_rows = A.num_rows;
                B->num_cols = A.num_cols;
                B->nnz = A.nnz;
                B->max_nnz_row = A.max_nnz_row;
                B->diameter = A.diameter;
    
                // conversion
                CHECK( magma_tally4_index_malloc_cpu( &row_tmp, A.num_rows+1 ));
                //fill the row-pointer
                for( magma_tally4_int_t i=0; i < A.num_rows+1; i++ )
                    row_tmp[i] = i*A.max_nnz_row;
                // sort the diagonal element into the right place
                CHECK( magma_tally4_smalloc_cpu( &val_tmp2, A.num_rows*A.max_nnz_row ));
                CHECK( magma_tally4_index_malloc_cpu( &col_tmp2, A.num_rows*A.max_nnz_row ));

                for( magma_tally4_int_t j=0; j < A.num_rows; j++ ) {
                    magma_tally4_index_t diagcol = A.col[j*A.max_nnz_row];
                    magma_tally4_int_t smaller = 0;
                    for( magma_tally4_int_t i=1; i < A.max_nnz_row; i++ ) {
                        if ( (A.col[j*A.max_nnz_row+i] < diagcol)
                             && (A.val[j*A.max_nnz_row+i] !=  zero) )
                            smaller++;
                    }
                    for( magma_tally4_int_t i=0; i < smaller; i++ ) {
                        col_tmp2[j*A.max_nnz_row+i] = A.col[j*A.max_nnz_row+i+1];
                        val_tmp2[j*A.max_nnz_row+i] = A.val[j*A.max_nnz_row+i+1];
                    }
                    col_tmp2[j*A.max_nnz_row+smaller] = A.col[j*A.max_nnz_row];
                    val_tmp2[j*A.max_nnz_row+smaller] = A.val[j*A.max_nnz_row];
                    for( magma_tally4_int_t i=smaller+1; i < A.max_nnz_row; i++ ) {
                        col_tmp2[j*A.max_nnz_row+i] = A.col[j*A.max_nnz_row+i];
                        val_tmp2[j*A.max_nnz_row+i] = A.val[j*A.max_nnz_row+i];
                    }
                }
    
                //now use AA_ELL, IA_ELL, row_tmp as CSR with some zeros.
                //The CSR compressor removes these
                CHECK( magma_tally4_s_csr_compressor(&val_tmp2, &row_tmp, &col_tmp2,
                           &B->val, &B->row, &B->col, &B->num_rows, queue ));
                B->nnz = B->row[B->num_rows];
            }
            
            // ELLRT to CSR
            else if ( old_format == Magma_tally4_ELLRT ) {
                //printf( "Conversion to CSR: " );
                //fflush(stdout);
                // fill in information for B
                B->storage_type = Magma_tally4_CSR;
                B->memory_location = A.memory_location;
                B->fill_mode = A.fill_mode;
                B->num_rows = A.num_rows;
                B->num_cols = A.num_cols;
                B->nnz = A.nnz;
                B->max_nnz_row = A.max_nnz_row;
                B->diameter = A.diameter;
    
                magma_tally4_int_t threads_per_row = A.alignment;
                magma_tally4_int_t rowlength = magma_tally4_roundup( A.max_nnz_row, threads_per_row );
                // conversion
                magma_tally4_index_t *row_tmp;
                CHECK( magma_tally4_index_malloc_cpu( &row_tmp, A.num_rows+1 ));
                //fill the row-pointer
                for( magma_tally4_int_t i=0; i < A.num_rows+1; i++ )
                    row_tmp[i] = i*rowlength;
                //now use AA_ELL, IA_ELL, row_tmp as CSR with some zeros.
                //The CSR compressor removes these
                CHECK( magma_tally4_s_csr_compressor(&A.val, &row_tmp, &A.col,
                       &B->val, &B->row, &B->col, &B->num_rows, queue ));
                B->nnz = B->row[B->num_rows];
                //printf( "done\n" );
            }
            
            // SELLP to CSR
            else if ( old_format == Magma_tally4_SELLP ) {
                // printf( "Conversion to CSR: " );
                // fill in information for B
                B->storage_type = Magma_tally4_CSR;
                B->memory_location = A.memory_location;
                B->fill_mode = A.fill_mode;
                B->num_rows = A.num_rows;
                B->num_cols = A.num_cols;
                B->nnz = A.nnz;
                B->max_nnz_row = A.max_nnz_row;
                B->diameter = A.diameter;
                magma_tally4_int_t C = A.blocksize;
                magma_tally4_int_t slices = A.numblocks;
                B->blocksize = A.blocksize;
                B->numblocks = A.numblocks;
                // conversion
                CHECK( magma_tally4_smalloc_cpu( &val_tmp, A.max_nnz_row*(A.num_rows+C) ));
                CHECK( magma_tally4_index_malloc_cpu( &row_tmp, A.num_rows+C ));
                CHECK( magma_tally4_index_malloc_cpu( &col_tmp, A.max_nnz_row*(A.num_rows+C) ));
                // zero everything
                for(magma_tally4_int_t i=0; i < A.max_nnz_row*(A.num_rows+C); i++ ) {
                    val_tmp[ i ] = MAGMA_tally4_S_MAKE(0., 0.);
                    col_tmp[ i ] =  0;
                }
    
                //fill the row-pointer
                for( magma_tally4_int_t i=0; i < A.num_rows+1; i++ ) {
                    row_tmp[i] = A.max_nnz_row*i;
                }
    
                //transform RowMajor to ColMajor
                for( magma_tally4_int_t k=0; k < slices; k++) {
                    magma_tally4_int_t blockinfo = (A.row[k+1]-A.row[k])/A.blocksize;
                    for( magma_tally4_int_t j=0; j < C; j++ ) {
                        for( magma_tally4_int_t i=0; i < blockinfo; i++ ) {
                            col_tmp[ (k*C+j)*A.max_nnz_row+i ] =
                                                    A.col[A.row[k]+i*C+j];
                            val_tmp[ (k*C+j)*A.max_nnz_row+i ] =
                                                    A.val[A.row[k]+i*C+j];
                        }
                    }
                }
    
                //now use AA_ELL, IA_ELL, row_tmp as CSR with some zeros.
                //The CSR compressor removes these
    
                CHECK( magma_tally4_s_csr_compressor(&val_tmp, &row_tmp, &col_tmp,
                           &B->val, &B->row, &B->col, &B->num_rows, queue ));
                B->nnz = B->row[B->num_rows];
                //printf( "done\n" );
            }
            
            // DENSE to CSR
            else if ( old_format == Magma_tally4_DENSE ) {
                //printf( "Conversion to CSR: " );
                // fill in information for B
                B->storage_type = Magma_tally4_CSR;
                B->memory_location = A.memory_location;
                B->fill_mode = A.fill_mode;
                B->num_rows = A.num_rows;
                B->num_cols = A.num_cols;
                B->nnz = A.nnz;
                B->max_nnz_row = A.max_nnz_row;
                B->diameter = A.diameter;
    
                // conversion
    
                B->nnz=0;
                for( magma_tally4_int_t i=0; i<(A.num_rows)*(A.num_cols); i++ ) {
                    if ( MAGMA_tally4_S_REAL(A.val[i]) != 0.0 )
                        (B->nnz)++;
                }
                CHECK( magma_tally4_smalloc_cpu( &B->val, B->nnz));
                CHECK( magma_tally4_index_malloc_cpu( &B->row, B->num_rows+1 ));
                CHECK( magma_tally4_index_malloc_cpu( &B->col, B->nnz ));
                
                magma_tally4_int_t i = 0;
                magma_tally4_int_t j = 0;
                magma_tally4_int_t k = 0;
    
                for(i=0; i<(A.num_rows)*(A.num_cols); i++)
                {
                    if ( i%(B->num_cols) == 0 )
                    {
                        (B->row)[k] = j;
                        k++;
                    }
                    if ( MAGMA_tally4_S_REAL(A.val[i]) != 0 )
                    {
                        (B->val)[j] = A.val[i];
                        (B->col)[j] = i%(B->num_cols);
                        j++;
                    }
                }
                (B->row)[B->num_rows]=B->nnz;
    
                //printf( "done\n" );
            }
            
            // BCSR to CSR
            else if ( old_format == Magma_tally4_BCSR ) {
                CHECK( magma_tally4_smtransfer(A, &A_d, Magma_tally4_CPU, Magma_tally4_DEV, queue ) );
                CHECK( magma_tally4_smconvert(A_d, &B_d, Magma_tally4_BCSR, Magma_tally4_CSR, queue ) );
                CHECK( magma_tally4_smtransfer(B_d, B, Magma_tally4_DEV, Magma_tally4_CPU, queue ) );
            }
            
            else {
                printf("error: format not supported.\n");
                magma_tally4blasSetKernelStream( orig_queue );
                info = MAGMA_tally4_ERR_NOT_SUPPORTED;
            }
        }
        else {
            printf("error: conversion not supported.\n");
            magma_tally4blasSetKernelStream( orig_queue );
            info = MAGMA_tally4_ERR_NOT_SUPPORTED;
        }
    } // end CPU case
    else if ( A.memory_location == Magma_tally4_DEV ) {
        // CSR to CSR
        if ( old_format == Magma_tally4_CSR && new_format == Magma_tally4_CSR ) {
            CHECK( magma_tally4_smtransfer( A, B, Magma_tally4_DEV, Magma_tally4_DEV, queue ));
        }
        // CSR to DENSE
        if ( old_format == Magma_tally4_CSR && new_format == Magma_tally4_DENSE ) {
            // fill in information for B
            B->storage_type = Magma_tally4_DENSE;
            B->memory_location = A.memory_location;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;

            // CUSPARSE context //
            CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
            CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
            CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));
            // end CUSPARSE context //

            CHECK( magma_tally4_smalloc( &B->dval, A.num_rows*A.num_cols ));
            
            
            // conversion using CUSPARSE
            cusparseScsr2dense( cusparseHandle, A.num_rows, A.num_cols,
                                descr, A.dval, A.drow, A.dcol,
                                B->dval, A.num_rows );
        }
        // DENSE to CSR
        else if ( old_format == Magma_tally4_DENSE && new_format == Magma_tally4_CSR ) {
            // fill in information for B
            B->storage_type = Magma_tally4_CSR;
            B->memory_location = A.memory_location;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;

            // CUSPARSE context //
            CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
            CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
            CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));
            // end CUSPARSE context //


            intnnz = B->nnz;
            CHECK( magma_tally4_index_malloc( &nnz_per_row, A.num_rows ));
            //magma_tally4_sprint_gpu( A.num_rows, 1, nnz_per_row, A.num_rows )
            cusparseSnnz( cusparseHandle, CUSPARSE_DIRECTION_COLUMN,
                          A.num_rows, A.num_cols,
                          descr,
                          A.dval, A.num_rows, nnz_per_row, &intnnz );

            CHECK( magma_tally4_smalloc( &B->dval, B->nnz ));
            CHECK( magma_tally4_index_malloc( &B->drow, B->num_rows+1 ));
            CHECK( magma_tally4_index_malloc( &B->dcol, B->nnz ));
            
            
            // conversion using CUSPARSE
            cusparseSdense2csr( cusparseHandle, A.num_rows, A.num_cols,
                                descr,
                                A.dval, A.num_rows, nnz_per_row,
                                B->dval, B->drow, B->dcol );
        }
        // CSR to BCSR
        else if ( old_format == Magma_tally4_CSR && new_format == Magma_tally4_BCSR ) {
            //printf( "Conversion to BCSR: " );
            // fill in information for B
            B->storage_type = Magma_tally4_BCSR;
            B->memory_location = A.memory_location;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;
            magma_tally4_int_t size_b = B->blocksize;

            // CUSPARSE context //
            CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
            CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
            CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));
            // end CUSPARSE context //

            magma_tally4_index_t base, nnzb;
            magma_tally4_int_t mb = magma_tally4_ceildiv( A.num_rows, size_b );
            // nnzTotalDevHostPtr points to host memory
            magma_tally4_index_t *nnzTotalDevHostPtr = &nnzb;

            CHECK( magma_tally4_index_malloc( &B->drow, mb+1 ));
            cusparseXcsr2bsrNnz( cusparseHandle, CUSPARSE_DIRECTION_COLUMN,
                                 A.num_rows, A.num_cols, descr,
                                 A.drow, A.dcol, size_b,
                                 descr, B->drow, nnzTotalDevHostPtr );
            
            
            if (NULL != nnzTotalDevHostPtr) {
                nnzb = *nnzTotalDevHostPtr;
            } else {
                magma_tally4_index_getvector( 1, B->row+mb, 1, &nnzb, 1 );
                magma_tally4_index_getvector( 1, B->row, 1, &base, 1 );
                nnzb -= base;
            }
            B->numblocks = nnzb; // number of blocks

            CHECK( magma_tally4_smalloc( &B->dval, nnzb*size_b*size_b ));
            CHECK( magma_tally4_index_malloc( &B->dcol, nnzb ));
            
            
            // conversion using CUSPARSE
            cusparseScsr2bsr( cusparseHandle, CUSPARSE_DIRECTION_ROW,
                              A.num_rows, A.num_cols, descr,
                              A.dval, A.drow, A.dcol,
                              size_b, descr,
                              B->dval, B->drow, B->dcol);
        }
        // BCSR to CSR
        else if ( old_format == Magma_tally4_BCSR && new_format == Magma_tally4_CSR ) {
            //printf( "Conversion to CSR: " );
            // fill in information for B
            B->storage_type = Magma_tally4_CSR;
            B->memory_location = A.memory_location;
            B->diameter = A.diameter;

            magma_tally4_int_t size_b = A.blocksize;

            // CUSPARSE context //
            CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
            CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
            CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));
            // end CUSPARSE context //

            magma_tally4_int_t mb = magma_tally4_ceildiv( A.num_rows, size_b );
            magma_tally4_int_t nb = magma_tally4_ceildiv( A.num_cols, size_b );
            magma_tally4_int_t nnzb = A.numblocks; // number of blocks
            B->nnz  = nnzb * size_b * size_b; // number of elements
            B->num_rows = mb * size_b;
            B->num_cols = nb * size_b;

            CHECK( magma_tally4_smalloc( &B->dval, B->nnz ));
            CHECK( magma_tally4_index_malloc( &B->drow, B->num_rows+1 ));
            CHECK( magma_tally4_index_malloc( &B->dcol, B->nnz ));
            
            
            // conversion using CUSPARSE
            cusparseSbsr2csr( cusparseHandle, CUSPARSE_DIRECTION_ROW,
                              mb, nb, descr, A.dval, A.drow, A.dcol,
                              size_b, descr,
                              B->dval, B->drow, B->dcol );
        }
        // CSR to CSC
        else if ( old_format == Magma_tally4_CSR && new_format == Magma_tally4_CSC ) {
            // fill in information for B
            B->storage_type = Magma_tally4_CSC;
            B->memory_location = A.memory_location;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;

            // CUSPARSE context //
            CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
            CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
            CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));
            // end CUSPARSE context //

            CHECK( magma_tally4_smalloc( &B->dval, B->nnz ));
            CHECK( magma_tally4_index_malloc( &B->drow, B->nnz ));
            CHECK( magma_tally4_index_malloc( &B->dcol, B->num_cols+1 ));
            
            
            // conversion using CUSPARSE
            cusparseScsr2csc(cusparseHandle, A.num_rows, A.num_cols, A.nnz,
                             A.dval, A.drow, A.dcol,
                             B->dval, B->drow, B->dcol,
                             CUSPARSE_ACTION_NUMERIC,
                             CUSPARSE_INDEX_BASE_ZERO);
        }
        // CSC to CSR
        else if ( old_format == Magma_tally4_CSC && new_format == Magma_tally4_CSR ) {
            // fill in information for B
            B->storage_type = Magma_tally4_CSR;
            B->memory_location = A.memory_location;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;

            // CUSPARSE context //
            CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
            CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
            CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));
            // end CUSPARSE context //

            CHECK( magma_tally4_smalloc( &B->dval, B->nnz ));
            CHECK( magma_tally4_index_malloc( &B->drow, B->num_rows+1 ));
            CHECK( magma_tally4_index_malloc( &B->dcol, B->nnz ));
            
            
            // conversion using CUSPARSE
            cusparseScsr2csc(cusparseHandle, A.num_rows, A.num_cols, A.nnz,
                             A.dval, A.dcol, A.drow,
                             B->dval, B->dcol, B->drow,
                             CUSPARSE_ACTION_NUMERIC,
                             CUSPARSE_INDEX_BASE_ZERO);
        }
        // CSR to COO
        else if ( old_format == Magma_tally4_CSR && new_format == Magma_tally4_COO ) {
            // fill in information for B
            B->storage_type = Magma_tally4_COO;
            B->memory_location = A.memory_location;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;

            // CUSPARSE context //
            CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
            CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
            CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));
            // end CUSPARSE context //

            CHECK( magma_tally4_smalloc( &B->dval, B->nnz ));
            CHECK( magma_tally4_index_malloc( &B->drow, B->nnz ));
            CHECK( magma_tally4_index_malloc( &B->dcol, B->nnz ));
            
            
            magma_tally4_scopyvector( A.nnz, A.dval, 1, B->dval, 1 );
            magma_tally4_index_copyvector( A.nnz, A.dcol, 1, B->dcol, 1 );

            // conversion using CUSPARSE
            cusparseXcsr2coo( cusparseHandle, A.drow,
                              A.nnz, A.num_rows, B->drow,
                              CUSPARSE_INDEX_BASE_ZERO );
        }
        // COO to CSR
        else if ( old_format == Magma_tally4_COO && new_format == Magma_tally4_CSR ) {
            // fill in information for B
            B->storage_type = Magma_tally4_CSR;
            B->memory_location = A.memory_location;
            B->num_rows = A.num_rows;
            B->num_cols = A.num_cols;
            B->nnz = A.nnz;
            B->max_nnz_row = A.max_nnz_row;
            B->diameter = A.diameter;

            // CUSPARSE context //

            CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
            CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
            CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));
            // end CUSPARSE context //

            CHECK( magma_tally4_smalloc( &B->dval, B->nnz ));
            CHECK( magma_tally4_index_malloc( &B->drow, B->nnz ));
            CHECK( magma_tally4_index_malloc( &B->dcol, B->nnz ));
            
            
            magma_tally4_scopyvector( A.nnz, A.val, 1, B->val, 1 );
            magma_tally4_index_copyvector( A.nnz, A.col, 1, B->col, 1 );

            // conversion using CUSPARSE
            cusparseXcoo2csr( cusparseHandle, A.drow,
                              A.nnz, A.num_rows, B->drow,
                              CUSPARSE_INDEX_BASE_ZERO );
        }
        else {
            printf("warning: format not supported on GPU. "
            "Conversion handled by CPU.\n");
            CHECK( magma_tally4_smtransfer( A, &hA, A.memory_location, Magma_tally4_CPU, queue ));
            CHECK( magma_tally4_smconvert( hA, &hB, old_format, new_format, queue ));
            CHECK( magma_tally4_smtransfer( hB, B, Magma_tally4_CPU, A.memory_location, queue ));
        }
    }
    
cleanup:    
    cusparseDestroyMatDescr(descr); 
    cusparseDestroy(cusparseHandle);
    descr = NULL;
    cusparseHandle = NULL;
    magma_tally4_free( nnz_per_row );
    magma_tally4_free_cpu( row_tmp );
    magma_tally4_free_cpu( col_tmp );
    magma_tally4_free_cpu( val_tmp );
    magma_tally4_free_cpu( row_tmp2 );
    magma_tally4_free_cpu( col_tmp2 );
    magma_tally4_free_cpu( val_tmp2 );
    row_tmp = NULL;
    col_tmp = NULL;
    val_tmp = NULL;
    row_tmp2 = NULL;
    col_tmp2 = NULL;
    val_tmp2 = NULL;   
    magma_tally4_free( transpose );
    magma_tally4_free_cpu( length );
    length = NULL;
    magma_tally4_smfree( &hA, queue );
    magma_tally4_smfree( &hB, queue );
    magma_tally4_smfree( &A_d, queue );
    magma_tally4_smfree( &B_d, queue );
    if ( info != 0 ) {
        magma_tally4_smfree( B, queue );
    }
    magma_tally4blasSetKernelStream( orig_queue );
    return info;
}