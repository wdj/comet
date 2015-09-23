/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt
       @author Mark Gates

*/
#include "common_magma_minproductsparse.h"

/**
    Purpose
    -------
    Transposes a matrix stored in CSR format on the CPU host.


    Arguments
    ---------
    @param[in]
    n_rows      magma_minproduct_int_t
                number of rows in input matrix

    @param[in]
    n_cols      magma_minproduct_int_t
                number of columns in input matrix

    @param[in]
    nnz         magma_minproduct_int_t
                number of nonzeros in input matrix

    @param[in]
    values      magma_minproductDoubleComplex*
                value array of input matrix

    @param[in]
    rowptr      magma_minproduct_index_t*
                row pointer of input matrix

    @param[in]
    colind      magma_minproduct_index_t*
                column indices of input matrix

    @param[in]
    new_n_rows  magma_minproduct_index_t*
                number of rows in transposed matrix

    @param[in]
    new_n_cols  magma_minproduct_index_t*
                number of columns in transposed matrix

    @param[in]
    new_nnz     magma_minproduct_index_t*
                number of nonzeros in transposed matrix

    @param[in]
    new_values  magma_minproductDoubleComplex**
                value array of transposed matrix

    @param[in]
    new_rowptr  magma_minproduct_index_t**
                row pointer of transposed matrix

    @param[in]
    new_colind  magma_minproduct_index_t**
                column indices of transposed matrix

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zaux
    ********************************************************************/

extern "C" magma_minproduct_int_t
z_transpose_csr(
    magma_minproduct_int_t n_rows,
    magma_minproduct_int_t n_cols,
    magma_minproduct_int_t nnz,
    magma_minproductDoubleComplex *values,
    magma_minproduct_index_t *rowptr,
    magma_minproduct_index_t *colind,
    magma_minproduct_int_t *new_n_rows,
    magma_minproduct_int_t *new_n_cols,
    magma_minproduct_int_t *new_nnz,
    magma_minproductDoubleComplex **new_values,
    magma_minproduct_index_t **new_rowptr,
    magma_minproduct_index_t **new_colind,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // easier to keep names straight if we convert CSR to CSC,
    // which is the same as tranposing CSR.
    magma_minproductDoubleComplex *csc_values=NULL;
    magma_minproduct_index_t *csc_colptr=NULL, *csc_rowind=NULL;
    
    // i, j are actual row & col indices (0 <= i < nrows, 0 <= j < ncols).
    // k is index into col and values (0 <= k < nnz).
    magma_minproduct_int_t i, j, k, total, tmp;
    
    CHECK( magma_minproduct_zmalloc_cpu( &csc_values, nnz ) );
    CHECK( magma_minproduct_index_malloc_cpu( &csc_colptr, n_cols + 1 ) );
    CHECK( magma_minproduct_index_malloc_cpu( &csc_rowind, nnz ) );
    
    // example matrix
    // [ x x 0 x ]
    // [ x 0 x x ]
    // [ x x 0 0 ]
    // rowptr = [ 0 3 6, 8 ]
    // colind = [ 0 1 3 ; 0 2 3 ; 0 1 ]
    
    // sum up nnz in each original column
    // colptr = [ 3 2 1 2, X ]
    for( j=0; j < n_cols; ++j ) {
        csc_colptr[ j ] = 0;
    }
    for( k=0; k < nnz; ++k ) {
        csc_colptr[ colind[k] ]++;
    }
    
    // running sum to convert to new colptr
    // colptr = [ 0 3 5 6, 8 ]
    total = 0;
    for( j=0; j < n_cols; ++j ) {
        tmp = csc_colptr[ j ];
        csc_colptr[ j ] = total;
        total += tmp;
    }
    csc_colptr[ n_cols ] = total;
    assert( total == nnz );
    
    // copy row indices and values
    // this increments colptr until it effectively shifts left one
    // colptr = [ 3 5 6 8, 8 ]
    // rowind = [ 0 1 2 ; 0 2 ; 1 ; 0 1 ]
    for( i=0; i < n_rows; ++i ) {
        for( k=rowptr[i]; k < rowptr[i+1]; ++k ) {
            j = colind[k];
            csc_rowind[ csc_colptr[ j ] ] = i;
            csc_values[ csc_colptr[ j ] ] = values[k];
            csc_colptr[ j ]++;
        }
    }
    assert( csc_colptr[ n_cols-1 ] == nnz );
    
    // shift colptr right one
    // colptr = [ 0 3 5 6, 8 ]
    for( j=n_cols-1; j > 0; --j ) {
        csc_colptr[j] = csc_colptr[j-1];
    }
    csc_colptr[0] = 0;

    // save into output variables
    *new_n_rows = n_cols;
    *new_n_cols = n_rows;
    *new_nnz    = nnz;
    *new_values = csc_values;
    *new_rowptr = csc_colptr;
    *new_colind = csc_rowind;
    
cleanup:
    magma_minproduct_free_cpu( csc_values );
    magma_minproduct_free_cpu( csc_colptr );
    magma_minproduct_free_cpu( csc_rowind );
    return info;
}


/**
    Purpose
    -------

    Interface to cuSPARSE transpose.

    Arguments
    ---------

    @param[in]
    A           magma_minproduct_z_matrix
                input matrix (CSR)

    @param[out]
    B           magma_minproduct_z_matrix*
                output matrix (CSR)
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zaux
    ********************************************************************/
    
    
extern "C" magma_minproduct_int_t
magma_minproduct_zmtranspose(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix *B,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    CHECK( magma_minproduct_z_cucsrtranspose( A, B, queue ));
    
cleanup:
    return info;
}


/**
    Purpose
    -------

    Helper function to transpose CSR matrix.
    Using the CUSPARSE CSR2CSC function.


    Arguments
    ---------

    @param[in]
    A           magma_minproduct_z_matrix
                input matrix (CSR)

    @param[out]
    B           magma_minproduct_z_matrix*
                output matrix (CSR)
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zaux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_z_cucsrtranspose(
    magma_minproduct_z_matrix A,
    magma_minproduct_z_matrix *B,
    magma_minproduct_queue_t queue )
{
    // for symmetric matrices: convert to csc using cusparse
    
    magma_minproduct_int_t info = 0;
    cusparseHandle_t handle=NULL;
    cusparseMatDescr_t descrA=NULL;
    cusparseMatDescr_t descrB=NULL;
    
    
    magma_minproduct_z_matrix ACSR={Magma_minproduct_CSR}, BCSR={Magma_minproduct_CSR};
    magma_minproduct_z_matrix A_d={Magma_minproduct_CSR}, B_d={Magma_minproduct_CSR};

    if( A.storage_type == Magma_minproduct_CSR && A.memory_location == Magma_minproduct_DEV ) {
                  
        // fill in information for B
        B->storage_type    = A.storage_type;
        B->diagorder_type  = A.diagorder_type;
        B->memory_location = Magma_minproduct_DEV;
        B->num_rows        = A.num_cols;  // transposed
        B->num_cols        = A.num_rows;  // transposed
        B->nnz             = A.nnz;
        
        if ( A.fill_mode == Magma_minproduct_FULL ) {
            B->fill_mode = Magma_minproduct_FULL;
        }
        else if ( A.fill_mode == Magma_minproduct_LOWER ) {
            B->fill_mode = Magma_minproduct_UPPER;
        }
        else if ( A.fill_mode == Magma_minproduct_UPPER ) {
            B->fill_mode = Magma_minproduct_LOWER;
        }
        
        B->dval = NULL;
        B->drow = NULL;
        B->dcol = NULL;
        
        // memory allocation
        CHECK( magma_minproduct_zmalloc( &B->dval, B->nnz ));
        CHECK( magma_minproduct_index_malloc( &B->drow, B->num_rows + 1 ));
        CHECK( magma_minproduct_index_malloc( &B->dcol, B->nnz ));
        
        // CUSPARSE context //
        CHECK_CUSPARSE( cusparseCreate( &handle ));
        CHECK_CUSPARSE( cusparseSetStream( handle, queue ));
        CHECK_CUSPARSE( cusparseCreateMatDescr( &descrA ));
        CHECK_CUSPARSE( cusparseCreateMatDescr( &descrB ));
        CHECK_CUSPARSE( cusparseSetMatType( descrA, CUSPARSE_MATRIX_TYPE_GENERAL ));
        CHECK_CUSPARSE( cusparseSetMatType( descrB, CUSPARSE_MATRIX_TYPE_GENERAL ));
        CHECK_CUSPARSE( cusparseSetMatIndexBase( descrA, CUSPARSE_INDEX_BASE_ZERO ));
        CHECK_CUSPARSE( cusparseSetMatIndexBase( descrB, CUSPARSE_INDEX_BASE_ZERO ));
        CHECK_CUSPARSE(
        cusparseZcsr2csc( handle, A.num_rows, A.num_cols, A.nnz,
                          A.dval, A.drow, A.dcol, B->dval, B->dcol, B->drow,
                          CUSPARSE_ACTION_NUMERIC,
                          CUSPARSE_INDEX_BASE_ZERO) );
        
    }else if( A.storage_type == Magma_minproduct_CSR && A.memory_location == Magma_minproduct_CPU ){
        CHECK( magma_minproduct_zmtransfer( A, &A_d, A.memory_location, Magma_minproduct_DEV, queue ));
        CHECK( magma_minproduct_z_cucsrtranspose( A_d, &B_d, queue ));
        CHECK( magma_minproduct_zmtransfer( B_d, B, Magma_minproduct_DEV, A.memory_location, queue ));
                
    }else {
        CHECK( magma_minproduct_zmconvert( A, &ACSR, A.storage_type, Magma_minproduct_CSR, queue ));
        CHECK( magma_minproduct_z_cucsrtranspose( ACSR, &BCSR, queue ));
        CHECK( magma_minproduct_zmconvert( BCSR, B, Magma_minproduct_CSR, A.storage_type, queue ));
    }
cleanup:
    cusparseDestroyMatDescr( descrA );
    cusparseDestroyMatDescr( descrB );
    cusparseDestroy( handle );
    magma_minproduct_zmfree( &A_d, queue );
    magma_minproduct_zmfree( &B_d, queue );
    magma_minproduct_zmfree( &ACSR, queue );
    magma_minproduct_zmfree( &BCSR, queue );
    if( info != 0 ){
        magma_minproduct_zmfree( B, queue );
    }
    return info;
}



