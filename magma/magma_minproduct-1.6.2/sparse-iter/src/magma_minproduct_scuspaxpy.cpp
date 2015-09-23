/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_minproduct_zcuspaxpy.cpp normal z -> s, Sun May  3 11:22:59 2015
       @author Hartwig Anzt

*/
#include "common_magma_minproductsparse.h"

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
    alpha       float*
                scalar

    @param[in]
    A           magma_minproduct_s_matrix
                input matrix

    @param[in]
    beta        float*
                scalar

    @param[in]
    B           magma_minproduct_s_matrix
                input matrix

    @param[out]
    AB          magma_minproduct_s_matrix*
                output matrix AB = alpha * A + beta * B

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_sblas
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_scuspaxpy(
    float *alpha, magma_minproduct_s_matrix A,
    float *beta, magma_minproduct_s_matrix B,
    magma_minproduct_s_matrix *AB,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_s_matrix C={Magma_minproduct_CSR};
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
   
    magma_minproduct_index_t base_t, nnz_t, baseC;
    
    cusparseHandle_t handle=NULL;
    cusparseMatDescr_t descrA=NULL;
    cusparseMatDescr_t descrB=NULL;
    cusparseMatDescr_t descrC=NULL;
                             
    if (    A.memory_location == Magma_minproduct_DEV
        && B.memory_location == Magma_minproduct_DEV
        && ( A.storage_type == Magma_minproduct_CSR ||
             A.storage_type == Magma_minproduct_CSRCOO )
        && ( B.storage_type == Magma_minproduct_CSR ||
             B.storage_type == Magma_minproduct_CSRCOO ) ) {

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
        magma_minproduct_index_t *nnzTotalDevHostPtr = (magma_minproduct_index_t*) &C.nnz;
        CHECK_CUSPARSE( cusparseSetPointerMode( handle, CUSPARSE_POINTER_MODE_HOST ));
        CHECK( magma_minproduct_index_malloc( &C.drow, (A.num_rows + 1) ));
        CHECK_CUSPARSE( cusparseXcsrgeamNnz( handle, A.num_rows, A.num_cols,
                    descrA, A.nnz, A.drow, A.dcol,
                    descrB, B.nnz, B.drow, B.dcol,
                    descrC, C.row, nnzTotalDevHostPtr ));

        if (NULL != nnzTotalDevHostPtr) {
            C.nnz = *nnzTotalDevHostPtr;
        } else {
            // workaround as nnz and base C are magma_minproduct_int_t
            magma_minproduct_index_getvector( 1, C.drow+C.num_rows, 1, &nnz_t, 1 );
            magma_minproduct_index_getvector( 1, C.drow,   1, &base_t,    1 );
            C.nnz = (magma_minproduct_int_t) nnz_t;
            baseC = (magma_minproduct_int_t) base_t;
            C.nnz -= baseC;
        }
        CHECK( magma_minproduct_index_malloc( &C.dcol, C.nnz ));
        CHECK( magma_minproduct_smalloc( &C.dval, C.nnz ));
        CHECK_CUSPARSE( cusparseScsrgeam( handle, A.num_rows, A.num_cols,
                alpha,
                descrA, A.nnz,
                A.dval, A.drow, A.dcol,
                beta,
                descrB, B.nnz,
                B.dval, B.drow, B.dcol,
                descrC,
                C.dval, C.drow, C.dcol ));
        // end CUSPARSE context //

        CHECK( magma_minproduct_smtransfer( C, AB, Magma_minproduct_DEV, Magma_minproduct_DEV, queue ));
    }
    else {
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED; 
    }
    
cleanup:
    cusparseDestroyMatDescr( descrA );
    cusparseDestroyMatDescr( descrB );
    cusparseDestroyMatDescr( descrC );
    cusparseDestroy( handle );
    magma_minproduct_smfree( &C, queue );
    return info;
}





