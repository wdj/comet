/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @generated from zcumilu.cpp normal z -> d, Sun May  3 11:22:59 2015
*/
#include "common_magma_tally4sparse.h"

#define PRECISION_d


/**
    Purpose
    -------

    Prepares the ILU preconditioner via the cuSPARSE.

    Arguments
    ---------

    @param[in]
    A           magma_tally4_d_matrix
                input matrix A

    @param[in,out]
    precond     magma_tally4_d_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_dgepr
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_dcumilusetup(
    magma_tally4_d_matrix A,
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrA=NULL;
    cusparseMatDescr_t descrL=NULL;
    cusparseMatDescr_t descrU=NULL;
    
    //magma_tally4_dprint_matrix(A, queue );
    // copy matrix into preconditioner parameter
    magma_tally4_d_matrix hA={Magma_tally4_CSR}, hACSR={Magma_tally4_CSR};
    magma_tally4_d_matrix hL={Magma_tally4_CSR}, hU={Magma_tally4_CSR};
    CHECK( magma_tally4_dmtransfer( A, &hA, A.memory_location, Magma_tally4_CPU, queue ));
    CHECK( magma_tally4_dmconvert( hA, &hACSR, hA.storage_type, Magma_tally4_CSR, queue ));

        // in case using fill-in
    if( precond->levels > 0 ){
        magma_tally4_d_matrix hAL={Magma_tally4_CSR}, hAUt={Magma_tally4_CSR};
        CHECK( magma_tally4_dsymbilu( &hACSR, precond->levels, &hAL, &hAUt,  queue ));
        magma_tally4_dmfree(&hAL, queue);
        magma_tally4_dmfree(&hAUt, queue);
    }

    CHECK( magma_tally4_dmtransfer(hACSR, &(precond->M), Magma_tally4_CPU, Magma_tally4_DEV, queue ));

    magma_tally4_dmfree( &hA, queue );
    magma_tally4_dmfree( &hACSR, queue );

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrA ));
    CHECK_CUSPARSE( cusparseSetMatType( descrA, CUSPARSE_MATRIX_TYPE_GENERAL ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrA, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrA, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &(precond->cuinfo) ));
    CHECK_CUSPARSE( cusparseDcsrsm_analysis( cusparseHandle,
                CUSPARSE_OPERATION_NON_TRANSPOSE,
                precond->M.num_rows, precond->M.nnz, descrA,
                precond->M.dval, precond->M.drow, precond->M.dcol,
                precond->cuinfo ));
    CHECK_CUSPARSE( cusparseDcsrilu0( cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                      precond->M.num_rows, descrA,
                      precond->M.dval,
                      precond->M.drow,
                      precond->M.dcol,
                      precond->cuinfo ));

    CHECK( magma_tally4_dmtransfer( precond->M, &hA, Magma_tally4_DEV, Magma_tally4_CPU, queue ));

    hL.diagorder_type = Magma_tally4_UNITY;
    CHECK( magma_tally4_dmconvert( hA, &hL , Magma_tally4_CSR, Magma_tally4_CSRL, queue ));
    hU.diagorder_type = Magma_tally4_VALUE;
    CHECK( magma_tally4_dmconvert( hA, &hU , Magma_tally4_CSR, Magma_tally4_CSRU, queue ));
    CHECK( magma_tally4_dmtransfer( hL, &(precond->L), Magma_tally4_CPU, Magma_tally4_DEV, queue ));
    CHECK( magma_tally4_dmtransfer( hU, &(precond->U), Magma_tally4_CPU, Magma_tally4_DEV, queue ));


    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrL ));
    CHECK_CUSPARSE( cusparseSetMatType( descrL, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrL, CUSPARSE_DIAG_TYPE_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrL, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrL, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseDcsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->L.num_rows,
        precond->L.nnz, descrL,
        precond->L.dval, precond->L.drow, precond->L.dcol, precond->cuinfoL ));

    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_UPPER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoU ));
    CHECK_CUSPARSE( cusparseDcsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->U.num_rows,
        precond->U.nnz, descrU,
        precond->U.dval, precond->U.drow, precond->U.dcol, precond->cuinfoU ));


    if( precond->maxiter < 50 ){
        //prepare for iterative solves
        
        // extract the diagonal of L into precond->d
        CHECK( magma_tally4_djacobisetup_diagscal( precond->L, &precond->d, queue ));
        CHECK( magma_tally4_dvinit( &precond->work1, Magma_tally4_DEV, hA.num_rows, 1, MAGMA_tally4_D_ZERO, queue ));
        
        // extract the diagonal of U into precond->d2
        CHECK( magma_tally4_djacobisetup_diagscal( precond->U, &precond->d2, queue ));
        CHECK( magma_tally4_dvinit( &precond->work2, Magma_tally4_DEV, hA.num_rows, 1, MAGMA_tally4_D_ZERO, queue ));
    }

    
cleanup:
    cusparseDestroySolveAnalysisInfo( precond->cuinfo );
    cusparseDestroyMatDescr( descrA );
    cusparseDestroyMatDescr( descrL );
    cusparseDestroyMatDescr( descrU );
    cusparseDestroy( cusparseHandle );
    magma_tally4_dmfree( &hA, queue );
    magma_tally4_dmfree( &hACSR, queue );
    magma_tally4_dmfree(&hA, queue );
    magma_tally4_dmfree(&hL, queue );
    magma_tally4_dmfree(&hU, queue );

    return info;
}



/**
    Purpose
    -------

    Prepares the ILU triangular solves via cuSPARSE using an ILU factorization
    matrix stored either in precond->M or on the device as
    precond->L and precond->U.

    Arguments
    ---------

    @param[in,out]
    precond     magma_tally4_d_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_dgepr
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_dcumilugeneratesolverinfo(
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrL=NULL;
    cusparseMatDescr_t descrU=NULL;
    
    magma_tally4_d_matrix hA={Magma_tally4_CSR}, hL={Magma_tally4_CSR}, hU={Magma_tally4_CSR};
    
    if (precond->L.memory_location != Magma_tally4_DEV ){
        
        CHECK( magma_tally4_dmtransfer( precond->M, &hA,
            precond->M.memory_location, Magma_tally4_CPU, queue ));

        hL.diagorder_type = Magma_tally4_UNITY;
        CHECK( magma_tally4_dmconvert( hA, &hL , Magma_tally4_CSR, Magma_tally4_CSRL, queue ));
        hU.diagorder_type = Magma_tally4_VALUE;
        CHECK( magma_tally4_dmconvert( hA, &hU , Magma_tally4_CSR, Magma_tally4_CSRU, queue ));
        CHECK( magma_tally4_dmtransfer( hL, &(precond->L), Magma_tally4_CPU, Magma_tally4_DEV, queue ));
        CHECK( magma_tally4_dmtransfer( hU, &(precond->U), Magma_tally4_CPU, Magma_tally4_DEV, queue ));
        
        magma_tally4_dmfree(&hA, queue );
        magma_tally4_dmfree(&hL, queue );
        magma_tally4_dmfree(&hU, queue );
        
    }
    
    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));


    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrL ));
    CHECK_CUSPARSE( cusparseSetMatType( descrL, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrL, CUSPARSE_DIAG_TYPE_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrL, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrL, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseDcsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->L.num_rows,
        precond->L.nnz, descrL,
        precond->L.dval, precond->L.drow, precond->L.dcol, precond->cuinfoL ));


    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_UPPER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoU ));
    CHECK_CUSPARSE( cusparseDcsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->U.num_rows,
        precond->U.nnz, descrU,
        precond->U.dval, precond->U.drow, precond->U.dcol, precond->cuinfoU ));

    
    if( precond->maxiter < 50 ){
        //prepare for iterative solves

        // extract the diagonal of L into precond->d
        CHECK( magma_tally4_djacobisetup_diagscal( precond->L, &precond->d, queue ));
        CHECK( magma_tally4_dvinit( &precond->work1, Magma_tally4_DEV, precond->U.num_rows, 1, MAGMA_tally4_D_ZERO, queue ));
        
        // extract the diagonal of U into precond->d2
        CHECK( magma_tally4_djacobisetup_diagscal( precond->U, &precond->d2, queue ));
        CHECK( magma_tally4_dvinit( &precond->work2, Magma_tally4_DEV, precond->U.num_rows, 1, MAGMA_tally4_D_ZERO, queue ));
    }
    
cleanup:
    cusparseDestroyMatDescr( descrL );
    cusparseDestroyMatDescr( descrU );
    cusparseDestroy( cusparseHandle );
     
    return info;
}








/**
    Purpose
    -------

    Performs the left triangular solves using the ILU preconditioner.

    Arguments
    ---------

    @param[in]
    b           magma_tally4_d_matrix
                RHS

    @param[in,out]
    x           magma_tally4_d_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_tally4_d_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_dgepr
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_dapplycumilu_l(
    magma_tally4_d_matrix b,
    magma_tally4_d_matrix *x,
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrL=NULL;
    
    double one = MAGMA_tally4_D_MAKE( 1.0, 0.0);

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrL ));
    CHECK_CUSPARSE( cusparseSetMatType( descrL, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrL, CUSPARSE_DIAG_TYPE_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrL, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrL, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseDcsrsm_solve( cusparseHandle,
                            CUSPARSE_OPERATION_NON_TRANSPOSE,
                            precond->L.num_rows,
                            b.num_rows*b.num_cols/precond->L.num_rows,
                            &one,
                            descrL,
                            precond->L.dval,
                            precond->L.drow,
                            precond->L.dcol,
                            precond->cuinfoL,
                            b.dval,
                            precond->L.num_rows,
                            x->dval,
                            precond->L.num_rows ));
    
    magma_tally4_device_sync();

cleanup:
    cusparseDestroyMatDescr( descrL );
    cusparseDestroy( cusparseHandle );
    return info;
}


/**
    Purpose
    -------

    Performs the right triangular solves using the ILU preconditioner.

    Arguments
    ---------

    @param[in]
    b           magma_tally4_d_matrix
                RHS

    @param[in,out]
    x           magma_tally4_d_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_tally4_d_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_dgepr
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_dapplycumilu_r(
    magma_tally4_d_matrix b,
    magma_tally4_d_matrix *x,
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrU=NULL;
    
    double one = MAGMA_tally4_D_MAKE( 1.0, 0.0);

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_UPPER ));
    CHECK_CUSPARSE( cusparseDcsrsm_solve( cusparseHandle,
                            CUSPARSE_OPERATION_NON_TRANSPOSE,
                            precond->U.num_rows,
                            b.num_rows*b.num_cols/precond->U.num_rows,
                            &one,
                            descrU,
                            precond->U.dval,
                            precond->U.drow,
                            precond->U.dcol,
                            precond->cuinfoU,
                            b.dval,
                            precond->U.num_rows,
                            x->dval,
                            precond->U.num_rows ));
    
    magma_tally4_device_sync();

cleanup:
    cusparseDestroyMatDescr( descrU );
    cusparseDestroy( cusparseHandle );
    return info; 
}




/**
    Purpose
    -------

    Prepares the IC preconditioner via cuSPARSE.

    Arguments
    ---------

    @param[in]
    A           magma_tally4_d_matrix
                input matrix A

    @param[in,out]
    precond     magma_tally4_d_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_dhepr
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_dcumiccsetup(
    magma_tally4_d_matrix A,
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrA=NULL;
    cusparseMatDescr_t descrL=NULL;
    cusparseMatDescr_t descrU=NULL;
    
    magma_tally4_d_matrix hA={Magma_tally4_CSR}, hACSR={Magma_tally4_CSR}, U={Magma_tally4_CSR};
    CHECK( magma_tally4_dmtransfer( A, &hA, A.memory_location, Magma_tally4_CPU, queue ));
    U.diagorder_type = Magma_tally4_VALUE;
    CHECK( magma_tally4_dmconvert( hA, &hACSR, hA.storage_type, Magma_tally4_CSR, queue ));

    // in case using fill-in
    if( precond->levels > 0 ){
            magma_tally4_d_matrix hAL={Magma_tally4_CSR}, hAUt={Magma_tally4_CSR};
            CHECK( magma_tally4_dsymbilu( &hACSR, precond->levels, &hAL, &hAUt,  queue ));
            magma_tally4_dmfree(&hAL, queue);
            magma_tally4_dmfree(&hAUt, queue);
    }

    CHECK( magma_tally4_dmconvert( hACSR, &U, Magma_tally4_CSR, Magma_tally4_CSRL, queue ));
    magma_tally4_dmfree( &hACSR, queue );
    CHECK( magma_tally4_dmtransfer(U, &(precond->M), Magma_tally4_CPU, Magma_tally4_DEV, queue ));

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrA ));
    CHECK_CUSPARSE( cusparseSetMatType( descrA, CUSPARSE_MATRIX_TYPE_SYMMETRIC ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrA, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrA, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrA, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &(precond->cuinfo) ));
    CHECK_CUSPARSE( cusparseDcsrsm_analysis( cusparseHandle,
                CUSPARSE_OPERATION_NON_TRANSPOSE,
                precond->M.num_rows, precond->M.nnz, descrA,
                precond->M.dval, precond->M.drow, precond->M.dcol,
                precond->cuinfo ));
    CHECK_CUSPARSE( cusparseDcsric0( cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                      precond->M.num_rows, descrA,
                      precond->M.dval,
                      precond->M.drow,
                      precond->M.dcol,
                      precond->cuinfo ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrL ));
    CHECK_CUSPARSE( cusparseSetMatType( descrL, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrL, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrL, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrL, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseDcsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->M.num_rows,
        precond->M.nnz, descrL,
        precond->M.dval, precond->M.drow, precond->M.dcol, precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoU ));
    CHECK_CUSPARSE( cusparseDcsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_TRANSPOSE, precond->M.num_rows,
        precond->M.nnz, descrU,
        precond->M.dval, precond->M.drow, precond->M.dcol, precond->cuinfoU ));
    
    if( precond->maxiter < 50 ){
        //prepare for iterative solves
        
        // copy the matrix to precond->L and (transposed) to precond->U
        CHECK( magma_tally4_dmtransfer(precond->M, &(precond->L), Magma_tally4_DEV, Magma_tally4_DEV, queue ));
        CHECK( magma_tally4_dmtranspose( precond->L, &(precond->U), queue ));
        
        // extract the diagonal of L into precond->d
        CHECK( magma_tally4_djacobisetup_diagscal( precond->L, &precond->d, queue ));
        CHECK( magma_tally4_dvinit( &precond->work1, Magma_tally4_DEV, hA.num_rows, 1, MAGMA_tally4_D_ZERO, queue ));
        
        // extract the diagonal of U into precond->d2
        CHECK( magma_tally4_djacobisetup_diagscal( precond->U, &precond->d2, queue ));
        CHECK( magma_tally4_dvinit( &precond->work2, Magma_tally4_DEV, hA.num_rows, 1, MAGMA_tally4_D_ZERO, queue ));
    }



/*
    // to enable also the block-asynchronous iteration for the triangular solves
    CHECK( magma_tally4_dmtransfer( precond->M, &hA, Magma_tally4_DEV, Magma_tally4_CPU, queue ));
    hA.storage_type = Magma_tally4_CSR;

    magma_tally4_d_matrix hD, hR, hAt

    CHECK( magma_tally4_dcsrsplit( 256, hA, &hD, &hR, queue ));

    CHECK( magma_tally4_dmtransfer( hD, &precond->LD, Magma_tally4_CPU, Magma_tally4_DEV, queue ));
    CHECK( magma_tally4_dmtransfer( hR, &precond->L, Magma_tally4_CPU, Magma_tally4_DEV, queue ));

    magma_tally4_dmfree(&hD, queue );
    magma_tally4_dmfree(&hR, queue );

    CHECK( magma_tally4_d_cucsrtranspose(   hA, &hAt, queue ));

    CHECK( magma_tally4_dcsrsplit( 256, hAt, &hD, &hR, queue ));

    CHECK( magma_tally4_dmtransfer( hD, &precond->UD, Magma_tally4_CPU, Magma_tally4_DEV, queue ));
    CHECK( magma_tally4_dmtransfer( hR, &precond->U, Magma_tally4_CPU, Magma_tally4_DEV, queue ));
    
    magma_tally4_dmfree(&hD, queue );
    magma_tally4_dmfree(&hR, queue );
    magma_tally4_dmfree(&hA, queue );
    magma_tally4_dmfree(&hAt, queue );
*/

cleanup:
    cusparseDestroySolveAnalysisInfo( precond->cuinfo );
    cusparseDestroyMatDescr( descrL );
    cusparseDestroyMatDescr( descrU );
    cusparseDestroyMatDescr( descrA );
    cusparseDestroy( cusparseHandle );
    magma_tally4_dmfree(&U, queue );
    magma_tally4_dmfree(&hA, queue );
    
    return info;
}


/**
    Purpose
    -------

    Prepares the IC preconditioner solverinfo via cuSPARSE for a triangular
    matrix present on the device in precond->M.

    Arguments
    ---------
    
    @param[in,out]
    precond     magma_tally4_d_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_dhepr
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_dcumicgeneratesolverinfo(
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrL=NULL;
    cusparseMatDescr_t descrU=NULL;
    
    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrL ));
    CHECK_CUSPARSE( cusparseSetMatType( descrL, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrL, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrL, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrL, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseDcsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->M.num_rows,
        precond->M.nnz, descrL,
        precond->M.dval, precond->M.drow, precond->M.dcol, precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoU ));
    CHECK_CUSPARSE( cusparseDcsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_TRANSPOSE, precond->M.num_rows,
        precond->M.nnz, descrU,
        precond->M.dval, precond->M.drow, precond->M.dcol, precond->cuinfoU ));


/*
    // to enable also the block-asynchronous iteration for the triangular solves
    CHECK( magma_tally4_dmtransfer( precond->M, &hA, Magma_tally4_DEV, Magma_tally4_CPU, queue ));
    hA.storage_type = Magma_tally4_CSR;

    CHECK( magma_tally4_dcsrsplit( 256, hA, &hD, &hR, queue ));

    CHECK( magma_tally4_dmtransfer( hD, &precond->LD, Magma_tally4_CPU, Magma_tally4_DEV, queue ));
    CHECK( magma_tally4_dmtransfer( hR, &precond->L, Magma_tally4_CPU, Magma_tally4_DEV, queue ));

    magma_tally4_dmfree(&hD, queue );
    magma_tally4_dmfree(&hR, queue );

    CHECK( magma_tally4_d_cucsrtranspose(   hA, &hAt, queue ));

    CHECK( magma_tally4_dcsrsplit( 256, hAt, &hD, &hR, queue ));

    CHECK( magma_tally4_dmtransfer( hD, &precond->UD, Magma_tally4_CPU, Magma_tally4_DEV, queue ));
    CHECK( magma_tally4_dmtransfer( hR, &precond->U, Magma_tally4_CPU, Magma_tally4_DEV, queue ));
    
    magma_tally4_dmfree(&hD, queue );
    magma_tally4_dmfree(&hR, queue );
    magma_tally4_dmfree(&hA, queue );
    magma_tally4_dmfree(&hAt, queue );
*/

cleanup:
    cusparseDestroyMatDescr( descrL );
    cusparseDestroyMatDescr( descrU );
    cusparseDestroy( cusparseHandle );
    return info;
}



/**
    Purpose
    -------

    Performs the left triangular solves using the ICC preconditioner.

    Arguments
    ---------

    @param[in]
    b           magma_tally4_d_matrix
                RHS

    @param[in,out]
    x           magma_tally4_d_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_tally4_d_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_dhepr
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_dapplycumicc_l(
    magma_tally4_d_matrix b,
    magma_tally4_d_matrix *x,
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrL=NULL;
    
    double one = MAGMA_tally4_D_MAKE( 1.0, 0.0);

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrL ));
    CHECK_CUSPARSE( cusparseSetMatType( descrL, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrL, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrL, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrL, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseDcsrsm_solve( cusparseHandle,
                            CUSPARSE_OPERATION_NON_TRANSPOSE,
                            precond->M.num_rows,
                            b.num_rows*b.num_cols/precond->M.num_rows,
                            &one,
                            descrL,
                            precond->M.dval,
                            precond->M.drow,
                            precond->M.dcol,
                            precond->cuinfoL,
                            b.dval,
                            precond->M.num_rows,
                            x->dval,
                            precond->M.num_rows ));
    
    magma_tally4_device_sync();

cleanup:
    cusparseDestroyMatDescr( descrL );
    cusparseDestroy( cusparseHandle );
    return info; 
}




/**
    Purpose
    -------

    Performs the right triangular solves using the ICC preconditioner.

    Arguments
    ---------

    @param[in]
    b           magma_tally4_d_matrix
                RHS

    @param[in,out]
    x           magma_tally4_d_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_tally4_d_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_dhepr
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_dapplycumicc_r(
    magma_tally4_d_matrix b,
    magma_tally4_d_matrix *x,
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrU=NULL;
    
    double one = MAGMA_tally4_D_MAKE( 1.0, 0.0);

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseDcsrsm_solve( cusparseHandle,
                            CUSPARSE_OPERATION_TRANSPOSE,
                            precond->M.num_rows,
                            b.num_rows*b.num_cols/precond->M.num_rows,
                            &one,
                            descrU,
                            precond->M.dval,
                            precond->M.drow,
                            precond->M.dcol,
                            precond->cuinfoU,
                            b.dval,
                            precond->M.num_rows,
                            x->dval,
                            precond->M.num_rows ));
    
    magma_tally4_device_sync();

cleanup:
    cusparseDestroyMatDescr( descrU );
    cusparseDestroy( cusparseHandle );
    return info; 
}









