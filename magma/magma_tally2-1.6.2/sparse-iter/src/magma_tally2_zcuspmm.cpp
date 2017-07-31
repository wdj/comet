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

#define RTOLERANCE     lapackf77_dlamch( "E" )
#define ATOLERANCE     lapackf77_dlamch( "E" )


/**
    Purpose
    -------

    This is an interface to the cuSPARSE routine csrmm computing the product
    of two sparse matrices stored in csr format.


    Arguments
    ---------

    @param[in]
    A           magma_tally2_z_matrix
                input matrix

    @param[in]
    B           magma_tally2_z_matrix
                input matrix

    @param[out]
    AB          magma_tally2_z_matrix*
                output matrix AB = A * B

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_zblas
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_zcuspmm(
    magma_tally2_z_matrix A, magma_tally2_z_matrix B,
    magma_tally2_z_matrix *AB,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    
    magma_tally2_z_matrix C={Magma_tally2_CSR};
    C.num_rows = A.num_rows;
    C.num_cols = B.num_cols;
    C.storage_type = A.storage_type;
    C.memory_location = A.memory_location;
    C.fill_mode = Magma_tally2_FULL;
    
    C.val = NULL;
    C.col = NULL;
    C.row = NULL;
    C.rowidx = NULL;
    C.blockinfo = NULL;
    C.diag = NULL;
    C.dval = NULL;
    C.dcol = NULL;
    C.drow = NULL;
    C.drowidx = NULL;
    C.ddiag = NULL;
    
    magma_tally2_index_t base_t, nnz_t, baseC;
    
    cusparseHandle_t handle=NULL;
    cusparseMatDescr_t descrA=NULL;
    cusparseMatDescr_t descrB=NULL;
    cusparseMatDescr_t descrC=NULL;
    
    if (    A.memory_location == Magma_tally2_DEV
        && B.memory_location == Magma_tally2_DEV
        && ( A.storage_type == Magma_tally2_CSR ||
             A.storage_type == Magma_tally2_CSRCOO )
        && ( B.storage_type == Magma_tally2_CSR ||
             B.storage_type == Magma_tally2_CSRCOO ) ) {
            
        // CUSPARSE context /
        CHECK_CUSPARSE( cusparseCreate( &handle ));
        CHECK_CUSPARSE( cusparseSetStream( handle, queue ));
        CHECK_CUSPARSE( cusparseCreateMatDescr( &descrA ));
        CHECK_CUSPARSE( cusparseCreateMatDescr( &descrB ));
        CHECK_CUSPARSE( cusparseCreateMatDescr( &descrC ));
        CHECK_CUSPARSE( cusparseSetMatType( descrA, CUSPARSE_MATRIX_TYPE_GENERAL ));
        CHECK_CUSPARSE( cusparseSetMatType( descrB, CUSPARSE_MATRIX_TYPE_GENERAL ));
        CHECK_CUSPARSE( cusparseSetMatType( descrC, CUSPARSE_MATRIX_TYPE_GENERAL ));
        CHECK_CUSPARSE( cusparseSetMatIndexBase( descrA, CUSPARSE_INDEX_BASE_ZERO ));
        CHECK_CUSPARSE( cusparseSetMatIndexBase( descrB, CUSPARSE_INDEX_BASE_ZERO ));
        CHECK_CUSPARSE( cusparseSetMatIndexBase( descrC, CUSPARSE_INDEX_BASE_ZERO ));

        // nnzTotalDevHostPtr points to host memory
        magma_tally2_index_t *nnzTotalDevHostPtr = (magma_tally2_index_t*) &C.nnz;
        CHECK_CUSPARSE( cusparseSetPointerMode( handle, CUSPARSE_POINTER_MODE_HOST ));
        CHECK( magma_tally2_index_malloc( &C.drow, (A.num_rows + 1) ));
        CHECK_CUSPARSE( cusparseXcsrgemmNnz( handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                                    A.num_rows, B.num_cols, A.num_cols,
                                    descrA, A.nnz, A.drow, A.dcol,
                                    descrB, B.nnz, B.drow, B.dcol,
                                    descrC, C.drow, nnzTotalDevHostPtr ));
        if (NULL != nnzTotalDevHostPtr) {
            C.nnz = *nnzTotalDevHostPtr;
        } else {
            // workaround as nnz and base C are magma_tally2_int_t
            magma_tally2_index_getvector( 1, C.drow+C.num_rows, 1, &nnz_t, 1 );
            magma_tally2_index_getvector( 1, C.drow,   1, &base_t,    1 );
            C.nnz = (magma_tally2_int_t) nnz_t;
            baseC = (magma_tally2_int_t) base_t;
            C.nnz -= baseC;
        }
        CHECK( magma_tally2_index_malloc( &C.dcol, C.nnz ));
        CHECK( magma_tally2_zmalloc( &C.dval, C.nnz ));
        CHECK_CUSPARSE( cusparseZcsrgemm( handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                        A.num_rows, B.num_cols, A.num_cols,
                        descrA, A.nnz,
                        A.dval, A.drow, A.dcol,
                        descrB, B.nnz,
                        B.dval, B.drow, B.dcol,
                        descrC,
                        C.dval, C.drow, C.dcol ));
        // end CUSPARSE context //

        CHECK( magma_tally2_zmtransfer( C, AB, Magma_tally2_DEV, Magma_tally2_DEV, queue ));
    }
    else {
        info = MAGMA_tally2_ERR_NOT_SUPPORTED; 
    }
    
cleanup:
    cusparseDestroyMatDescr( descrA );
    cusparseDestroyMatDescr( descrB );
    cusparseDestroyMatDescr( descrC );
    cusparseDestroy( handle );
    magma_tally2_zmfree( &C, queue );
    return info;
}





