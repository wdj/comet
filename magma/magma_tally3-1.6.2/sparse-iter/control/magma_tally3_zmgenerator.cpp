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

    Generate a symmetric n x n CSR matrix for a stencil.

    Arguments
    ---------

    @param[in]
    n           magma_tally3_int_t
                number of rows

    @param[in]
    offdiags    magma_tally3_int_t
                number of offdiagonals

    @param[in]
    diag_offset magma_tally3_int_t*
                array containing the offsets

                                                (length offsets+1)
    @param[in]
    diag_vals   magma_tally3DoubleComplex*
                array containing the values

                                                (length offsets+1)
    @param[out]
    A           magma_tally3_z_matrix*
                matrix to generate
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zaux
    ********************************************************************/

extern "C"
magma_tally3_int_t
magma_tally3_zmgenerator(
    magma_tally3_int_t n,
    magma_tally3_int_t offdiags,
    magma_tally3_index_t *diag_offset,
    magma_tally3DoubleComplex *diag_vals,
    magma_tally3_z_matrix *A,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    magma_tally3_z_matrix B={Magma_tally3_CSR};
    
    B.val = NULL;
    B.col = NULL;
    B.row = NULL;
    B.rowidx = NULL;
    B.blockinfo = NULL;
    B.diag = NULL;
    B.dval = NULL;
    B.dcol = NULL;
    B.drow = NULL;
    B.drowidx = NULL;
    B.ddiag = NULL;
    
    B.num_rows = n;
    B.num_cols = n;
    B.fill_mode = Magma_tally3_FULL;
    B.memory_location = Magma_tally3_CPU;
    B.storage_type = Magma_tally3_ELLPACKT;
    B.max_nnz_row = (2*offdiags+1);

    CHECK( magma_tally3_zmalloc_cpu( &B.val, B.max_nnz_row*n ));
    CHECK( magma_tally3_index_malloc_cpu( &B.col, B.max_nnz_row*n ));
    
    for( int i=0; i<n; i++ ) { // stride over rows
        // stride over the number of nonzeros in each row
        // left of diagonal
        for( int j=0; j<offdiags; j++ ) {
            B.val[ i*B.max_nnz_row + j ] = diag_vals[ offdiags - j ];
            B.col[ i*B.max_nnz_row + j ] = -1 * diag_offset[ offdiags-j ] + i;
        }
        // elements on the diagonal
        B.val[ i*B.max_nnz_row + offdiags ] = diag_vals[ 0 ];
        B.col[ i*B.max_nnz_row + offdiags ] = i;
        // right of diagonal
        for( int j=0; j<offdiags; j++ ) {
            B.val[ i*B.max_nnz_row + j + offdiags +1 ] = diag_vals[ j+1 ];
            B.col[ i*B.max_nnz_row + j + offdiags +1 ] = diag_offset[ j+1 ] + i;
        }
    }

    // set invalid entries to zero
    for( int i=0; i<n; i++ ) { // stride over rows
        for( int j=0; j<B.max_nnz_row; j++ ) { // nonzeros in every row
            if ( (B.col[i*B.max_nnz_row + j] < 0) ||
                    (B.col[i*B.max_nnz_row + j] >= n) ) {
                B.val[ i*B.max_nnz_row + j ] = MAGMA_tally3_Z_MAKE( 0.0, 0.0 );
            }
        }

    }

    B.nnz = 0;

    for( int i=0; i<n; i++ ) { // stride over rows
        for( int j=0; j<B.max_nnz_row; j++ ) { // nonzeros in every row
            if ( MAGMA_tally3_Z_REAL( B.val[i*B.max_nnz_row + j]) != 0.0 )
                B.nnz++;
        }

    }

    // converting it to CSR will remove the invalit entries completely
    CHECK( magma_tally3_zmconvert( B, A, Magma_tally3_ELLPACKT, Magma_tally3_CSR, queue ));

cleanup:
    if( info != 0 ){
        magma_tally3_zmfree( &B, queue );
    }
    return info;
}




/**
    Purpose
    -------

    Generate a 27-point stencil for a 3D FD discretization.

    Arguments
    ---------

    @param[in]
    n           magma_tally3_int_t
                number of rows

    @param[out]
    A           magma_tally3_z_matrix*
                matrix to generate
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zaux
    ********************************************************************/

extern "C"
magma_tally3_int_t
magma_tally3_zm_27stencil(
    magma_tally3_int_t n,
    magma_tally3_z_matrix *A,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    magma_tally3_int_t i,j,k;
    magma_tally3_z_matrix hA={Magma_tally3_CSR};

    
    // generate matrix of desired structure and size (3d 27-point stencil)
    magma_tally3_int_t nn = n*n*n;
    magma_tally3_int_t offdiags = 13;
    magma_tally3_index_t *diag_offset=NULL;
    magma_tally3DoubleComplex *diag_vals=NULL;
    CHECK( magma_tally3_zmalloc_cpu( &diag_vals, offdiags+1 ));
    CHECK( magma_tally3_index_malloc_cpu( &diag_offset, offdiags+1 ));

    diag_offset[0] = 0;
    diag_offset[1] = 1;
    diag_offset[2] = n-1;
    diag_offset[3] = n;
    diag_offset[4] = n+1;
    diag_offset[5] = n*n-n-1;
    diag_offset[6] = n*n-n;
    diag_offset[7] = n*n-n+1;
    diag_offset[8] = n*n-1;
    diag_offset[9] = n*n;
    diag_offset[10] = n*n+1;
    diag_offset[11] = n*n+n-1;
    diag_offset[12] = n*n+n;
    diag_offset[13] = n*n+n+1;

    diag_vals[0] = MAGMA_tally3_Z_MAKE( 26.0, 0.0 );
    diag_vals[1] = MAGMA_tally3_Z_MAKE( -1.0, 0.0 );
    diag_vals[2] = MAGMA_tally3_Z_MAKE( -1.0, 0.0 );
    diag_vals[3] = MAGMA_tally3_Z_MAKE( -1.0, 0.0 );
    diag_vals[4] = MAGMA_tally3_Z_MAKE( -1.0, 0.0 );
    diag_vals[5] = MAGMA_tally3_Z_MAKE( -1.0, 0.0 );
    diag_vals[6] = MAGMA_tally3_Z_MAKE( -1.0, 0.0 );
    diag_vals[7] = MAGMA_tally3_Z_MAKE( -1.0, 0.0 );
    diag_vals[8] = MAGMA_tally3_Z_MAKE( -1.0, 0.0 );
    diag_vals[9] = MAGMA_tally3_Z_MAKE( -1.0, 0.0 );
    diag_vals[10] = MAGMA_tally3_Z_MAKE( -1.0, 0.0 );
    diag_vals[11] = MAGMA_tally3_Z_MAKE( -1.0, 0.0 );
    diag_vals[12] = MAGMA_tally3_Z_MAKE( -1.0, 0.0 );
    diag_vals[13] = MAGMA_tally3_Z_MAKE( -1.0, 0.0 );
    CHECK( magma_tally3_zmgenerator( nn, offdiags, diag_offset, diag_vals, &hA, queue ));


    // now set some entries to zero (boundary...)
    for( i=0; i<n*n; i++ ) {
    for( j=0; j<n; j++ ) {
        magma_tally3_index_t row = i*n+j;
        for( k=hA.row[row]; k<hA.row[row+1]; k++) {

            if ((hA.col[k] == row-1 ||
                hA.col[k] == row-n-1 ||
                hA.col[k] == row+n-1 ||
                hA.col[k] == row-n*n+n-1 ||
                hA.col[k] == row+n*n-n-1 ||
                hA.col[k] == row-n*n-1 ||
                hA.col[k] == row+n*n-1 ||
                hA.col[k] == row-n*n-n-1 ||
                hA.col[k] == row+n*n+n-1 ) && (row+1)%n == 1 )
                    
                    hA.val[k] = MAGMA_tally3_Z_MAKE( 0.0, 0.0 );

            if ((hA.col[k] == row+1 ||
                hA.col[k] == row-n+1 ||
                hA.col[k] == row+n+1 ||
                hA.col[k] == row-n*n+n+1 ||
                hA.col[k] == row+n*n-n+1 ||
                hA.col[k] == row-n*n+1 ||
                hA.col[k] == row+n*n+1 ||
                hA.col[k] == row-n*n-n+1 ||
                hA.col[k] == row+n*n+n+1 ) && (row)%n ==n-1 )
                    
                    hA.val[k] = MAGMA_tally3_Z_MAKE( 0.0, 0.0 );
        }
        
    }
    }
    CHECK( magma_tally3_zmconvert( hA, A, Magma_tally3_CSR, Magma_tally3_CSR, queue ));

cleanup:
    magma_tally3_free_cpu( diag_vals );
    magma_tally3_free_cpu( diag_offset );
    magma_tally3_zmfree( &hA, queue );
    return info;
}



/**
    Purpose
    -------

    Generate a 5-point stencil for a 2D FD discretization.

    Arguments
    ---------

    @param[in]
    n           magma_tally3_int_t
                number of rows

    @param[out]
    A           magma_tally3_z_matrix*
                matrix to generate
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zaux
    ********************************************************************/

extern "C"
magma_tally3_int_t
magma_tally3_zm_5stencil(
    magma_tally3_int_t n,
    magma_tally3_z_matrix *A,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    magma_tally3_int_t i,j,k;
    magma_tally3_z_matrix hA={Magma_tally3_CSR};
    
    // generate matrix of desired structure and size (2d 5-point stencil)
    magma_tally3_int_t nn = n*n;
    magma_tally3_int_t offdiags = 2;
    magma_tally3_index_t *diag_offset=NULL;
    magma_tally3DoubleComplex *diag_vals=NULL;
    CHECK( magma_tally3_zmalloc_cpu( &diag_vals, offdiags+1 ));
    CHECK( magma_tally3_index_malloc_cpu( &diag_offset, offdiags+1 ));

    diag_offset[0] = 0;
    diag_offset[1] = 1;
    diag_offset[2] = n;
    
    #define COMPLEX
    
    #ifdef COMPLEX
        // complex case
        diag_vals[0] = MAGMA_tally3_Z_MAKE( 4.0, 4.0 );
        diag_vals[1] = MAGMA_tally3_Z_MAKE( -1.0, -1.0 );
        diag_vals[2] = MAGMA_tally3_Z_MAKE( -1.0, -1.0 );
        
    #else
        // real case
        diag_vals[0] = MAGMA_tally3_Z_MAKE( 4.0, 0.0 );
        diag_vals[1] = MAGMA_tally3_Z_MAKE( -1.0, 0.0 );
        diag_vals[2] = MAGMA_tally3_Z_MAKE( -1.0, 0.0 );
    #endif
    CHECK( magma_tally3_zmgenerator( nn, offdiags, diag_offset, diag_vals, &hA, queue ));

    // now set some entries to zero (boundary...)
    for( i=0; i<n; i++ ) {
    for( j=0; j<n; j++ ) {
        magma_tally3_index_t row = i*n+j;
        for( k=hA.row[row]; k<hA.row[row+1]; k++) {

            if ((hA.col[k] == row-1 ) && (row+1)%n == 1 )
                    
                    hA.val[k] = MAGMA_tally3_Z_MAKE( 0.0, 0.0 );

            if ((hA.col[k] == row+1 ) && (row)%n ==n-1 )
                    
                    hA.val[k] = MAGMA_tally3_Z_MAKE( 0.0, 0.0 );
        }
        
    }
    }

    CHECK( magma_tally3_zmconvert( hA, A, Magma_tally3_CSR, Magma_tally3_CSR, queue ));
    
cleanup:
    magma_tally3_free_cpu( diag_vals );
    magma_tally3_free_cpu( diag_offset );
    magma_tally3_zmfree( &hA, queue );
    return info;
}