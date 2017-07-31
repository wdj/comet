/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @generated from zcumilu.cpp normal z -> s, Sun May  3 11:22:59 2015
*/
#include "common_magma_tally2sparse.h"

#define PRECISION_s


/**
    Purpose
    -------

    Prepares the ILU preconditioner via the cuSPARSE.

    Arguments
    ---------

    @param[in]
    A           magma_tally2_s_matrix
                input matrix A

    @param[in,out]
    precond     magma_tally2_s_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_sgepr
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_scumilusetup(
    magma_tally2_s_matrix A,
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrA=NULL;
    cusparseMatDescr_t descrL=NULL;
    cusparseMatDescr_t descrU=NULL;
    
    //magma_tally2_sprint_matrix(A, queue );
    // copy matrix into preconditioner parameter
    magma_tally2_s_matrix hA={Magma_tally2_CSR}, hACSR={Magma_tally2_CSR};
    magma_tally2_s_matrix hL={Magma_tally2_CSR}, hU={Magma_tally2_CSR};
    CHECK( magma_tally2_smtransfer( A, &hA, A.memory_location, Magma_tally2_CPU, queue ));
    CHECK( magma_tally2_smconvert( hA, &hACSR, hA.storage_type, Magma_tally2_CSR, queue ));

        // in case using fill-in
    if( precond->levels > 0 ){
        magma_tally2_s_matrix hAL={Magma_tally2_CSR}, hAUt={Magma_tally2_CSR};
        CHECK( magma_tally2_ssymbilu( &hACSR, precond->levels, &hAL, &hAUt,  queue ));
        magma_tally2_smfree(&hAL, queue);
        magma_tally2_smfree(&hAUt, queue);
    }

    CHECK( magma_tally2_smtransfer(hACSR, &(precond->M), Magma_tally2_CPU, Magma_tally2_DEV, queue ));

    magma_tally2_smfree( &hA, queue );
    magma_tally2_smfree( &hACSR, queue );

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrA ));
    CHECK_CUSPARSE( cusparseSetMatType( descrA, CUSPARSE_MATRIX_TYPE_GENERAL ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrA, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrA, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &(precond->cuinfo) ));
    CHECK_CUSPARSE( cusparseScsrsm_analysis( cusparseHandle,
                CUSPARSE_OPERATION_NON_TRANSPOSE,
                precond->M.num_rows, precond->M.nnz, descrA,
                precond->M.dval, precond->M.drow, precond->M.dcol,
                precond->cuinfo ));
    CHECK_CUSPARSE( cusparseScsrilu0( cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                      precond->M.num_rows, descrA,
                      precond->M.dval,
                      precond->M.drow,
                      precond->M.dcol,
                      precond->cuinfo ));

    CHECK( magma_tally2_smtransfer( precond->M, &hA, Magma_tally2_DEV, Magma_tally2_CPU, queue ));

    hL.diagorder_type = Magma_tally2_UNITY;
    CHECK( magma_tally2_smconvert( hA, &hL , Magma_tally2_CSR, Magma_tally2_CSRL, queue ));
    hU.diagorder_type = Magma_tally2_VALUE;
    CHECK( magma_tally2_smconvert( hA, &hU , Magma_tally2_CSR, Magma_tally2_CSRU, queue ));
    CHECK( magma_tally2_smtransfer( hL, &(precond->L), Magma_tally2_CPU, Magma_tally2_DEV, queue ));
    CHECK( magma_tally2_smtransfer( hU, &(precond->U), Magma_tally2_CPU, Magma_tally2_DEV, queue ));


    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrL ));
    CHECK_CUSPARSE( cusparseSetMatType( descrL, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrL, CUSPARSE_DIAG_TYPE_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrL, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrL, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseScsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->L.num_rows,
        precond->L.nnz, descrL,
        precond->L.dval, precond->L.drow, precond->L.dcol, precond->cuinfoL ));

    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_UPPER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoU ));
    CHECK_CUSPARSE( cusparseScsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->U.num_rows,
        precond->U.nnz, descrU,
        precond->U.dval, precond->U.drow, precond->U.dcol, precond->cuinfoU ));


    if( precond->maxiter < 50 ){
        //prepare for iterative solves
        
        // extract the diagonal of L into precond->d
        CHECK( magma_tally2_sjacobisetup_diagscal( precond->L, &precond->d, queue ));
        CHECK( magma_tally2_svinit( &precond->work1, Magma_tally2_DEV, hA.num_rows, 1, MAGMA_tally2_S_ZERO, queue ));
        
        // extract the diagonal of U into precond->d2
        CHECK( magma_tally2_sjacobisetup_diagscal( precond->U, &precond->d2, queue ));
        CHECK( magma_tally2_svinit( &precond->work2, Magma_tally2_DEV, hA.num_rows, 1, MAGMA_tally2_S_ZERO, queue ));
    }

    
cleanup:
    cusparseDestroySolveAnalysisInfo( precond->cuinfo );
    cusparseDestroyMatDescr( descrA );
    cusparseDestroyMatDescr( descrL );
    cusparseDestroyMatDescr( descrU );
    cusparseDestroy( cusparseHandle );
    magma_tally2_smfree( &hA, queue );
    magma_tally2_smfree( &hACSR, queue );
    magma_tally2_smfree(&hA, queue );
    magma_tally2_smfree(&hL, queue );
    magma_tally2_smfree(&hU, queue );

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
    precond     magma_tally2_s_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_sgepr
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_scumilugeneratesolverinfo(
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrL=NULL;
    cusparseMatDescr_t descrU=NULL;
    
    magma_tally2_s_matrix hA={Magma_tally2_CSR}, hL={Magma_tally2_CSR}, hU={Magma_tally2_CSR};
    
    if (precond->L.memory_location != Magma_tally2_DEV ){
        
        CHECK( magma_tally2_smtransfer( precond->M, &hA,
            precond->M.memory_location, Magma_tally2_CPU, queue ));

        hL.diagorder_type = Magma_tally2_UNITY;
        CHECK( magma_tally2_smconvert( hA, &hL , Magma_tally2_CSR, Magma_tally2_CSRL, queue ));
        hU.diagorder_type = Magma_tally2_VALUE;
        CHECK( magma_tally2_smconvert( hA, &hU , Magma_tally2_CSR, Magma_tally2_CSRU, queue ));
        CHECK( magma_tally2_smtransfer( hL, &(precond->L), Magma_tally2_CPU, Magma_tally2_DEV, queue ));
        CHECK( magma_tally2_smtransfer( hU, &(precond->U), Magma_tally2_CPU, Magma_tally2_DEV, queue ));
        
        magma_tally2_smfree(&hA, queue );
        magma_tally2_smfree(&hL, queue );
        magma_tally2_smfree(&hU, queue );
        
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
    CHECK_CUSPARSE( cusparseScsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->L.num_rows,
        precond->L.nnz, descrL,
        precond->L.dval, precond->L.drow, precond->L.dcol, precond->cuinfoL ));


    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_UPPER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoU ));
    CHECK_CUSPARSE( cusparseScsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->U.num_rows,
        precond->U.nnz, descrU,
        precond->U.dval, precond->U.drow, precond->U.dcol, precond->cuinfoU ));

    
    if( precond->maxiter < 50 ){
        //prepare for iterative solves

        // extract the diagonal of L into precond->d
        CHECK( magma_tally2_sjacobisetup_diagscal( precond->L, &precond->d, queue ));
        CHECK( magma_tally2_svinit( &precond->work1, Magma_tally2_DEV, precond->U.num_rows, 1, MAGMA_tally2_S_ZERO, queue ));
        
        // extract the diagonal of U into precond->d2
        CHECK( magma_tally2_sjacobisetup_diagscal( precond->U, &precond->d2, queue ));
        CHECK( magma_tally2_svinit( &precond->work2, Magma_tally2_DEV, precond->U.num_rows, 1, MAGMA_tally2_S_ZERO, queue ));
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
    b           magma_tally2_s_matrix
                RHS

    @param[in,out]
    x           magma_tally2_s_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_tally2_s_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_sgepr
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_sapplycumilu_l(
    magma_tally2_s_matrix b,
    magma_tally2_s_matrix *x,
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrL=NULL;
    
    float one = MAGMA_tally2_S_MAKE( 1.0, 0.0);

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrL ));
    CHECK_CUSPARSE( cusparseSetMatType( descrL, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrL, CUSPARSE_DIAG_TYPE_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrL, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrL, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseScsrsm_solve( cusparseHandle,
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
    
    magma_tally2_device_sync();

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
    b           magma_tally2_s_matrix
                RHS

    @param[in,out]
    x           magma_tally2_s_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_tally2_s_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_sgepr
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_sapplycumilu_r(
    magma_tally2_s_matrix b,
    magma_tally2_s_matrix *x,
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrU=NULL;
    
    float one = MAGMA_tally2_S_MAKE( 1.0, 0.0);

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_UPPER ));
    CHECK_CUSPARSE( cusparseScsrsm_solve( cusparseHandle,
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
    
    magma_tally2_device_sync();

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
    A           magma_tally2_s_matrix
                input matrix A

    @param[in,out]
    precond     magma_tally2_s_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_shepr
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_scumiccsetup(
    magma_tally2_s_matrix A,
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrA=NULL;
    cusparseMatDescr_t descrL=NULL;
    cusparseMatDescr_t descrU=NULL;
    
    magma_tally2_s_matrix hA={Magma_tally2_CSR}, hACSR={Magma_tally2_CSR}, U={Magma_tally2_CSR};
    CHECK( magma_tally2_smtransfer( A, &hA, A.memory_location, Magma_tally2_CPU, queue ));
    U.diagorder_type = Magma_tally2_VALUE;
    CHECK( magma_tally2_smconvert( hA, &hACSR, hA.storage_type, Magma_tally2_CSR, queue ));

    // in case using fill-in
    if( precond->levels > 0 ){
            magma_tally2_s_matrix hAL={Magma_tally2_CSR}, hAUt={Magma_tally2_CSR};
            CHECK( magma_tally2_ssymbilu( &hACSR, precond->levels, &hAL, &hAUt,  queue ));
            magma_tally2_smfree(&hAL, queue);
            magma_tally2_smfree(&hAUt, queue);
    }

    CHECK( magma_tally2_smconvert( hACSR, &U, Magma_tally2_CSR, Magma_tally2_CSRL, queue ));
    magma_tally2_smfree( &hACSR, queue );
    CHECK( magma_tally2_smtransfer(U, &(precond->M), Magma_tally2_CPU, Magma_tally2_DEV, queue ));

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrA ));
    CHECK_CUSPARSE( cusparseSetMatType( descrA, CUSPARSE_MATRIX_TYPE_SYMMETRIC ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrA, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrA, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrA, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &(precond->cuinfo) ));
    CHECK_CUSPARSE( cusparseScsrsm_analysis( cusparseHandle,
                CUSPARSE_OPERATION_NON_TRANSPOSE,
                precond->M.num_rows, precond->M.nnz, descrA,
                precond->M.dval, precond->M.drow, precond->M.dcol,
                precond->cuinfo ));
    CHECK_CUSPARSE( cusparseScsric0( cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
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
    CHECK_CUSPARSE( cusparseScsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->M.num_rows,
        precond->M.nnz, descrL,
        precond->M.dval, precond->M.drow, precond->M.dcol, precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoU ));
    CHECK_CUSPARSE( cusparseScsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_TRANSPOSE, precond->M.num_rows,
        precond->M.nnz, descrU,
        precond->M.dval, precond->M.drow, precond->M.dcol, precond->cuinfoU ));
    
    if( precond->maxiter < 50 ){
        //prepare for iterative solves
        
        // copy the matrix to precond->L and (transposed) to precond->U
        CHECK( magma_tally2_smtransfer(precond->M, &(precond->L), Magma_tally2_DEV, Magma_tally2_DEV, queue ));
        CHECK( magma_tally2_smtranspose( precond->L, &(precond->U), queue ));
        
        // extract the diagonal of L into precond->d
        CHECK( magma_tally2_sjacobisetup_diagscal( precond->L, &precond->d, queue ));
        CHECK( magma_tally2_svinit( &precond->work1, Magma_tally2_DEV, hA.num_rows, 1, MAGMA_tally2_S_ZERO, queue ));
        
        // extract the diagonal of U into precond->d2
        CHECK( magma_tally2_sjacobisetup_diagscal( precond->U, &precond->d2, queue ));
        CHECK( magma_tally2_svinit( &precond->work2, Magma_tally2_DEV, hA.num_rows, 1, MAGMA_tally2_S_ZERO, queue ));
    }



/*
    // to enable also the block-asynchronous iteration for the triangular solves
    CHECK( magma_tally2_smtransfer( precond->M, &hA, Magma_tally2_DEV, Magma_tally2_CPU, queue ));
    hA.storage_type = Magma_tally2_CSR;

    magma_tally2_s_matrix hD, hR, hAt

    CHECK( magma_tally2_scsrsplit( 256, hA, &hD, &hR, queue ));

    CHECK( magma_tally2_smtransfer( hD, &precond->LD, Magma_tally2_CPU, Magma_tally2_DEV, queue ));
    CHECK( magma_tally2_smtransfer( hR, &precond->L, Magma_tally2_CPU, Magma_tally2_DEV, queue ));

    magma_tally2_smfree(&hD, queue );
    magma_tally2_smfree(&hR, queue );

    CHECK( magma_tally2_s_cucsrtranspose(   hA, &hAt, queue ));

    CHECK( magma_tally2_scsrsplit( 256, hAt, &hD, &hR, queue ));

    CHECK( magma_tally2_smtransfer( hD, &precond->UD, Magma_tally2_CPU, Magma_tally2_DEV, queue ));
    CHECK( magma_tally2_smtransfer( hR, &precond->U, Magma_tally2_CPU, Magma_tally2_DEV, queue ));
    
    magma_tally2_smfree(&hD, queue );
    magma_tally2_smfree(&hR, queue );
    magma_tally2_smfree(&hA, queue );
    magma_tally2_smfree(&hAt, queue );
*/

cleanup:
    cusparseDestroySolveAnalysisInfo( precond->cuinfo );
    cusparseDestroyMatDescr( descrL );
    cusparseDestroyMatDescr( descrU );
    cusparseDestroyMatDescr( descrA );
    cusparseDestroy( cusparseHandle );
    magma_tally2_smfree(&U, queue );
    magma_tally2_smfree(&hA, queue );
    
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
    precond     magma_tally2_s_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_shepr
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_scumicgeneratesolverinfo(
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
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
    CHECK_CUSPARSE( cusparseScsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->M.num_rows,
        precond->M.nnz, descrL,
        precond->M.dval, precond->M.drow, precond->M.dcol, precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoU ));
    CHECK_CUSPARSE( cusparseScsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_TRANSPOSE, precond->M.num_rows,
        precond->M.nnz, descrU,
        precond->M.dval, precond->M.drow, precond->M.dcol, precond->cuinfoU ));


/*
    // to enable also the block-asynchronous iteration for the triangular solves
    CHECK( magma_tally2_smtransfer( precond->M, &hA, Magma_tally2_DEV, Magma_tally2_CPU, queue ));
    hA.storage_type = Magma_tally2_CSR;

    CHECK( magma_tally2_scsrsplit( 256, hA, &hD, &hR, queue ));

    CHECK( magma_tally2_smtransfer( hD, &precond->LD, Magma_tally2_CPU, Magma_tally2_DEV, queue ));
    CHECK( magma_tally2_smtransfer( hR, &precond->L, Magma_tally2_CPU, Magma_tally2_DEV, queue ));

    magma_tally2_smfree(&hD, queue );
    magma_tally2_smfree(&hR, queue );

    CHECK( magma_tally2_s_cucsrtranspose(   hA, &hAt, queue ));

    CHECK( magma_tally2_scsrsplit( 256, hAt, &hD, &hR, queue ));

    CHECK( magma_tally2_smtransfer( hD, &precond->UD, Magma_tally2_CPU, Magma_tally2_DEV, queue ));
    CHECK( magma_tally2_smtransfer( hR, &precond->U, Magma_tally2_CPU, Magma_tally2_DEV, queue ));
    
    magma_tally2_smfree(&hD, queue );
    magma_tally2_smfree(&hR, queue );
    magma_tally2_smfree(&hA, queue );
    magma_tally2_smfree(&hAt, queue );
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
    b           magma_tally2_s_matrix
                RHS

    @param[in,out]
    x           magma_tally2_s_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_tally2_s_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_shepr
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_sapplycumicc_l(
    magma_tally2_s_matrix b,
    magma_tally2_s_matrix *x,
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrL=NULL;
    
    float one = MAGMA_tally2_S_MAKE( 1.0, 0.0);

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrL ));
    CHECK_CUSPARSE( cusparseSetMatType( descrL, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrL, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrL, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrL, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseScsrsm_solve( cusparseHandle,
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
    
    magma_tally2_device_sync();

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
    b           magma_tally2_s_matrix
                RHS

    @param[in,out]
    x           magma_tally2_s_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_tally2_s_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_shepr
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_sapplycumicc_r(
    magma_tally2_s_matrix b,
    magma_tally2_s_matrix *x,
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrU=NULL;
    
    float one = MAGMA_tally2_S_MAKE( 1.0, 0.0);

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseScsrsm_solve( cusparseHandle,
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
    
    magma_tally2_device_sync();

cleanup:
    cusparseDestroyMatDescr( descrU );
    cusparseDestroy( cusparseHandle );
    return info; 
}









