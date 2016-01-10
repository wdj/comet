/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally4_zcuspaxpy.cpp normal z -> c, Sun May  3 11:22:59 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally4sparse.h"

#define RTOLERANCE     lapackf77_slamch( "E" )
#define ATOLERANCE     lapackf77_slamch( "E" )


/**
    Purpose
    -------

    This is an interface to the cuSPARSE routine csrgeam computing the sum
    of two sparse matrices stored in csr format:

        C = alpha * A + beta * B


    Arguments
    ---------

    @param[in]
    alpha       magma_tally4FloatComplex*
                scalar

    @param[in]
    A           magma_tally4_c_matrix
                input matrix

    @param[in]
    beta        magma_tally4FloatComplex*
                scalar

    @param[in]
    B           magma_tally4_c_matrix
                input matrix

    @param[out]
    AB          magma_tally4_c_matrix*
                output matrix AB = alpha * A + beta * B

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_cblas
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_ccuspaxpy(
    magma_tally4FloatComplex *alpha, magma_tally4_c_matrix A,
    magma_tally4FloatComplex *beta, magma_tally4_c_matrix B,
    magma_tally4_c_matrix *AB,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    magma_tally4_c_matrix C={Magma_tally4_CSR};
    C.num_rows = A.num_rows;
    C.num_cols = A.num_cols;
    C.storage_type = A.storage_type;
    C.memory_location = A.memory_location;
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
   
    magma_tally4_index_t base_t, nnz_t, baseC;
    
    cusparseHandle_t handle=NULL;
    cusparseMatDescr_t descrA=NULL;
    cusparseMatDescr_t descrB=NULL;
    cusparseMatDescr_t descrC=NULL;
                             
    if (    A.memory_location == Magma_tally4_DEV
        && B.memory_location == Magma_tally4_DEV
        && ( A.storage_type == Magma_tally4_CSR ||
             A.storage_type == Magma_tally4_CSRCOO )
        && ( B.storage_type == Magma_tally4_CSR ||
             B.storage_type == Magma_tally4_CSRCOO ) ) {

        // CUSPARSE context //

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
        magma_tally4_index_t *nnzTotalDevHostPtr = (magma_tally4_index_t*) &C.nnz;
        CHECK_CUSPARSE( cusparseSetPointerMode( handle, CUSPARSE_POINTER_MODE_HOST ));
        CHECK( magma_tally4_index_malloc( &C.drow, (A.num_rows + 1) ));
        CHECK_CUSPARSE( cusparseXcsrgeamNnz( handle, A.num_rows, A.num_cols,
                    descrA, A.nnz, A.drow, A.dcol,
                    descrB, B.nnz, B.drow, B.dcol,
                    descrC, C.row, nnzTotalDevHostPtr ));

        if (NULL != nnzTotalDevHostPtr) {
            C.nnz = *nnzTotalDevHostPtr;
        } else {
            // workaround as nnz and base C are magma_tally4_int_t
            magma_tally4_index_getvector( 1, C.drow+C.num_rows, 1, &nnz_t, 1 );
            magma_tally4_index_getvector( 1, C.drow,   1, &base_t,    1 );
            C.nnz = (magma_tally4_int_t) nnz_t;
            baseC = (magma_tally4_int_t) base_t;
            C.nnz -= baseC;
        }
        CHECK( magma_tally4_index_malloc( &C.dcol, C.nnz ));
        CHECK( magma_tally4_cmalloc( &C.dval, C.nnz ));
        CHECK_CUSPARSE( cusparseCcsrgeam( handle, A.num_rows, A.num_cols,
                alpha,
                descrA, A.nnz,
                A.dval, A.drow, A.dcol,
                beta,
                descrB, B.nnz,
                B.dval, B.drow, B.dcol,
                descrC,
                C.dval, C.drow, C.dcol ));
        // end CUSPARSE context //

        CHECK( magma_tally4_cmtransfer( C, AB, Magma_tally4_DEV, Magma_tally4_DEV, queue ));
    }
    else {
        info = MAGMA_tally4_ERR_NOT_SUPPORTED; 
    }
    
cleanup:
    cusparseDestroyMatDescr( descrA );
    cusparseDestroyMatDescr( descrB );
    cusparseDestroyMatDescr( descrC );
    cusparseDestroy( handle );
    magma_tally4_cmfree( &C, queue );
    return info;
}




