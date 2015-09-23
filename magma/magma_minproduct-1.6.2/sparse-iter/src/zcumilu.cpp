/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Hartwig Anzt

       @precisions normal z -> s d c
*/
#include "common_magma_minproductsparse.h"

#define PRECISION_z


/**
    Purpose
    -------

    Prepares the ILU preconditioner via the cuSPARSE.

    Arguments
    ---------

    @param[in]
    A           magma_minproduct_z_matrix
                input matrix A

    @param[in,out]
    precond     magma_minproduct_z_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zgepr
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zcumilusetup(
    magma_minproduct_z_matrix A,
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrA=NULL;
    cusparseMatDescr_t descrL=NULL;
    cusparseMatDescr_t descrU=NULL;
    
    //magma_minproduct_zprint_matrix(A, queue );
    // copy matrix into preconditioner parameter
    magma_minproduct_z_matrix hA={Magma_minproduct_CSR}, hACSR={Magma_minproduct_CSR};
    magma_minproduct_z_matrix hL={Magma_minproduct_CSR}, hU={Magma_minproduct_CSR};
    CHECK( magma_minproduct_zmtransfer( A, &hA, A.memory_location, Magma_minproduct_CPU, queue ));
    CHECK( magma_minproduct_zmconvert( hA, &hACSR, hA.storage_type, Magma_minproduct_CSR, queue ));

        // in case using fill-in
    if( precond->levels > 0 ){
        magma_minproduct_z_matrix hAL={Magma_minproduct_CSR}, hAUt={Magma_minproduct_CSR};
        CHECK( magma_minproduct_zsymbilu( &hACSR, precond->levels, &hAL, &hAUt,  queue ));
        magma_minproduct_zmfree(&hAL, queue);
        magma_minproduct_zmfree(&hAUt, queue);
    }

    CHECK( magma_minproduct_zmtransfer(hACSR, &(precond->M), Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));

    magma_minproduct_zmfree( &hA, queue );
    magma_minproduct_zmfree( &hACSR, queue );

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrA ));
    CHECK_CUSPARSE( cusparseSetMatType( descrA, CUSPARSE_MATRIX_TYPE_GENERAL ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrA, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrA, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &(precond->cuinfo) ));
    CHECK_CUSPARSE( cusparseZcsrsm_analysis( cusparseHandle,
                CUSPARSE_OPERATION_NON_TRANSPOSE,
                precond->M.num_rows, precond->M.nnz, descrA,
                precond->M.dval, precond->M.drow, precond->M.dcol,
                precond->cuinfo ));
    CHECK_CUSPARSE( cusparseZcsrilu0( cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                      precond->M.num_rows, descrA,
                      precond->M.dval,
                      precond->M.drow,
                      precond->M.dcol,
                      precond->cuinfo ));

    CHECK( magma_minproduct_zmtransfer( precond->M, &hA, Magma_minproduct_DEV, Magma_minproduct_CPU, queue ));

    hL.diagorder_type = Magma_minproduct_UNITY;
    CHECK( magma_minproduct_zmconvert( hA, &hL , Magma_minproduct_CSR, Magma_minproduct_CSRL, queue ));
    hU.diagorder_type = Magma_minproduct_VALUE;
    CHECK( magma_minproduct_zmconvert( hA, &hU , Magma_minproduct_CSR, Magma_minproduct_CSRU, queue ));
    CHECK( magma_minproduct_zmtransfer( hL, &(precond->L), Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
    CHECK( magma_minproduct_zmtransfer( hU, &(precond->U), Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));


    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrL ));
    CHECK_CUSPARSE( cusparseSetMatType( descrL, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrL, CUSPARSE_DIAG_TYPE_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrL, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrL, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseZcsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->L.num_rows,
        precond->L.nnz, descrL,
        precond->L.dval, precond->L.drow, precond->L.dcol, precond->cuinfoL ));

    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_UPPER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoU ));
    CHECK_CUSPARSE( cusparseZcsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->U.num_rows,
        precond->U.nnz, descrU,
        precond->U.dval, precond->U.drow, precond->U.dcol, precond->cuinfoU ));


    if( precond->maxiter < 50 ){
        //prepare for iterative solves
        
        // extract the diagonal of L into precond->d
        CHECK( magma_minproduct_zjacobisetup_diagscal( precond->L, &precond->d, queue ));
        CHECK( magma_minproduct_zvinit( &precond->work1, Magma_minproduct_DEV, hA.num_rows, 1, MAGMA_minproduct_Z_ZERO, queue ));
        
        // extract the diagonal of U into precond->d2
        CHECK( magma_minproduct_zjacobisetup_diagscal( precond->U, &precond->d2, queue ));
        CHECK( magma_minproduct_zvinit( &precond->work2, Magma_minproduct_DEV, hA.num_rows, 1, MAGMA_minproduct_Z_ZERO, queue ));
    }

    
cleanup:
    cusparseDestroySolveAnalysisInfo( precond->cuinfo );
    cusparseDestroyMatDescr( descrA );
    cusparseDestroyMatDescr( descrL );
    cusparseDestroyMatDescr( descrU );
    cusparseDestroy( cusparseHandle );
    magma_minproduct_zmfree( &hA, queue );
    magma_minproduct_zmfree( &hACSR, queue );
    magma_minproduct_zmfree(&hA, queue );
    magma_minproduct_zmfree(&hL, queue );
    magma_minproduct_zmfree(&hU, queue );

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
    precond     magma_minproduct_z_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zgepr
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zcumilugeneratesolverinfo(
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrL=NULL;
    cusparseMatDescr_t descrU=NULL;
    
    magma_minproduct_z_matrix hA={Magma_minproduct_CSR}, hL={Magma_minproduct_CSR}, hU={Magma_minproduct_CSR};
    
    if (precond->L.memory_location != Magma_minproduct_DEV ){
        
        CHECK( magma_minproduct_zmtransfer( precond->M, &hA,
            precond->M.memory_location, Magma_minproduct_CPU, queue ));

        hL.diagorder_type = Magma_minproduct_UNITY;
        CHECK( magma_minproduct_zmconvert( hA, &hL , Magma_minproduct_CSR, Magma_minproduct_CSRL, queue ));
        hU.diagorder_type = Magma_minproduct_VALUE;
        CHECK( magma_minproduct_zmconvert( hA, &hU , Magma_minproduct_CSR, Magma_minproduct_CSRU, queue ));
        CHECK( magma_minproduct_zmtransfer( hL, &(precond->L), Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
        CHECK( magma_minproduct_zmtransfer( hU, &(precond->U), Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
        
        magma_minproduct_zmfree(&hA, queue );
        magma_minproduct_zmfree(&hL, queue );
        magma_minproduct_zmfree(&hU, queue );
        
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
    CHECK_CUSPARSE( cusparseZcsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->L.num_rows,
        precond->L.nnz, descrL,
        precond->L.dval, precond->L.drow, precond->L.dcol, precond->cuinfoL ));


    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_UPPER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoU ));
    CHECK_CUSPARSE( cusparseZcsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->U.num_rows,
        precond->U.nnz, descrU,
        precond->U.dval, precond->U.drow, precond->U.dcol, precond->cuinfoU ));

    
    if( precond->maxiter < 50 ){
        //prepare for iterative solves

        // extract the diagonal of L into precond->d
        CHECK( magma_minproduct_zjacobisetup_diagscal( precond->L, &precond->d, queue ));
        CHECK( magma_minproduct_zvinit( &precond->work1, Magma_minproduct_DEV, precond->U.num_rows, 1, MAGMA_minproduct_Z_ZERO, queue ));
        
        // extract the diagonal of U into precond->d2
        CHECK( magma_minproduct_zjacobisetup_diagscal( precond->U, &precond->d2, queue ));
        CHECK( magma_minproduct_zvinit( &precond->work2, Magma_minproduct_DEV, precond->U.num_rows, 1, MAGMA_minproduct_Z_ZERO, queue ));
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
    b           magma_minproduct_z_matrix
                RHS

    @param[in,out]
    x           magma_minproduct_z_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_minproduct_z_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zgepr
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zapplycumilu_l(
    magma_minproduct_z_matrix b,
    magma_minproduct_z_matrix *x,
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrL=NULL;
    
    magma_minproductDoubleComplex one = MAGMA_minproduct_Z_MAKE( 1.0, 0.0);

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrL ));
    CHECK_CUSPARSE( cusparseSetMatType( descrL, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrL, CUSPARSE_DIAG_TYPE_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrL, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrL, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseZcsrsm_solve( cusparseHandle,
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
    
    magma_minproduct_device_sync();

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
    b           magma_minproduct_z_matrix
                RHS

    @param[in,out]
    x           magma_minproduct_z_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_minproduct_z_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zgepr
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zapplycumilu_r(
    magma_minproduct_z_matrix b,
    magma_minproduct_z_matrix *x,
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrU=NULL;
    
    magma_minproductDoubleComplex one = MAGMA_minproduct_Z_MAKE( 1.0, 0.0);

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_UPPER ));
    CHECK_CUSPARSE( cusparseZcsrsm_solve( cusparseHandle,
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
    
    magma_minproduct_device_sync();

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
    A           magma_minproduct_z_matrix
                input matrix A

    @param[in,out]
    precond     magma_minproduct_z_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zhepr
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zcumiccsetup(
    magma_minproduct_z_matrix A,
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrA=NULL;
    cusparseMatDescr_t descrL=NULL;
    cusparseMatDescr_t descrU=NULL;
    
    magma_minproduct_z_matrix hA={Magma_minproduct_CSR}, hACSR={Magma_minproduct_CSR}, U={Magma_minproduct_CSR};
    CHECK( magma_minproduct_zmtransfer( A, &hA, A.memory_location, Magma_minproduct_CPU, queue ));
    U.diagorder_type = Magma_minproduct_VALUE;
    CHECK( magma_minproduct_zmconvert( hA, &hACSR, hA.storage_type, Magma_minproduct_CSR, queue ));

    // in case using fill-in
    if( precond->levels > 0 ){
            magma_minproduct_z_matrix hAL={Magma_minproduct_CSR}, hAUt={Magma_minproduct_CSR};
            CHECK( magma_minproduct_zsymbilu( &hACSR, precond->levels, &hAL, &hAUt,  queue ));
            magma_minproduct_zmfree(&hAL, queue);
            magma_minproduct_zmfree(&hAUt, queue);
    }

    CHECK( magma_minproduct_zmconvert( hACSR, &U, Magma_minproduct_CSR, Magma_minproduct_CSRL, queue ));
    magma_minproduct_zmfree( &hACSR, queue );
    CHECK( magma_minproduct_zmtransfer(U, &(precond->M), Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrA ));
    CHECK_CUSPARSE( cusparseSetMatType( descrA, CUSPARSE_MATRIX_TYPE_SYMMETRIC ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrA, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrA, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrA, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &(precond->cuinfo) ));
    CHECK_CUSPARSE( cusparseZcsrsm_analysis( cusparseHandle,
                CUSPARSE_OPERATION_NON_TRANSPOSE,
                precond->M.num_rows, precond->M.nnz, descrA,
                precond->M.dval, precond->M.drow, precond->M.dcol,
                precond->cuinfo ));
    CHECK_CUSPARSE( cusparseZcsric0( cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
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
    CHECK_CUSPARSE( cusparseZcsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->M.num_rows,
        precond->M.nnz, descrL,
        precond->M.dval, precond->M.drow, precond->M.dcol, precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoU ));
    CHECK_CUSPARSE( cusparseZcsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_TRANSPOSE, precond->M.num_rows,
        precond->M.nnz, descrU,
        precond->M.dval, precond->M.drow, precond->M.dcol, precond->cuinfoU ));
    
    if( precond->maxiter < 50 ){
        //prepare for iterative solves
        
        // copy the matrix to precond->L and (transposed) to precond->U
        CHECK( magma_minproduct_zmtransfer(precond->M, &(precond->L), Magma_minproduct_DEV, Magma_minproduct_DEV, queue ));
        CHECK( magma_minproduct_zmtranspose( precond->L, &(precond->U), queue ));
        
        // extract the diagonal of L into precond->d
        CHECK( magma_minproduct_zjacobisetup_diagscal( precond->L, &precond->d, queue ));
        CHECK( magma_minproduct_zvinit( &precond->work1, Magma_minproduct_DEV, hA.num_rows, 1, MAGMA_minproduct_Z_ZERO, queue ));
        
        // extract the diagonal of U into precond->d2
        CHECK( magma_minproduct_zjacobisetup_diagscal( precond->U, &precond->d2, queue ));
        CHECK( magma_minproduct_zvinit( &precond->work2, Magma_minproduct_DEV, hA.num_rows, 1, MAGMA_minproduct_Z_ZERO, queue ));
    }



/*
    // to enable also the block-asynchronous iteration for the triangular solves
    CHECK( magma_minproduct_zmtransfer( precond->M, &hA, Magma_minproduct_DEV, Magma_minproduct_CPU, queue ));
    hA.storage_type = Magma_minproduct_CSR;

    magma_minproduct_z_matrix hD, hR, hAt

    CHECK( magma_minproduct_zcsrsplit( 256, hA, &hD, &hR, queue ));

    CHECK( magma_minproduct_zmtransfer( hD, &precond->LD, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
    CHECK( magma_minproduct_zmtransfer( hR, &precond->L, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));

    magma_minproduct_zmfree(&hD, queue );
    magma_minproduct_zmfree(&hR, queue );

    CHECK( magma_minproduct_z_cucsrtranspose(   hA, &hAt, queue ));

    CHECK( magma_minproduct_zcsrsplit( 256, hAt, &hD, &hR, queue ));

    CHECK( magma_minproduct_zmtransfer( hD, &precond->UD, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
    CHECK( magma_minproduct_zmtransfer( hR, &precond->U, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
    
    magma_minproduct_zmfree(&hD, queue );
    magma_minproduct_zmfree(&hR, queue );
    magma_minproduct_zmfree(&hA, queue );
    magma_minproduct_zmfree(&hAt, queue );
*/

cleanup:
    cusparseDestroySolveAnalysisInfo( precond->cuinfo );
    cusparseDestroyMatDescr( descrL );
    cusparseDestroyMatDescr( descrU );
    cusparseDestroyMatDescr( descrA );
    cusparseDestroy( cusparseHandle );
    magma_minproduct_zmfree(&U, queue );
    magma_minproduct_zmfree(&hA, queue );
    
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
    precond     magma_minproduct_z_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zhepr
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zcumicgeneratesolverinfo(
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
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
    CHECK_CUSPARSE( cusparseZcsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->M.num_rows,
        precond->M.nnz, descrL,
        precond->M.dval, precond->M.drow, precond->M.dcol, precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoU ));
    CHECK_CUSPARSE( cusparseZcsrsm_analysis( cusparseHandle,
        CUSPARSE_OPERATION_TRANSPOSE, precond->M.num_rows,
        precond->M.nnz, descrU,
        precond->M.dval, precond->M.drow, precond->M.dcol, precond->cuinfoU ));


/*
    // to enable also the block-asynchronous iteration for the triangular solves
    CHECK( magma_minproduct_zmtransfer( precond->M, &hA, Magma_minproduct_DEV, Magma_minproduct_CPU, queue ));
    hA.storage_type = Magma_minproduct_CSR;

    CHECK( magma_minproduct_zcsrsplit( 256, hA, &hD, &hR, queue ));

    CHECK( magma_minproduct_zmtransfer( hD, &precond->LD, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
    CHECK( magma_minproduct_zmtransfer( hR, &precond->L, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));

    magma_minproduct_zmfree(&hD, queue );
    magma_minproduct_zmfree(&hR, queue );

    CHECK( magma_minproduct_z_cucsrtranspose(   hA, &hAt, queue ));

    CHECK( magma_minproduct_zcsrsplit( 256, hAt, &hD, &hR, queue ));

    CHECK( magma_minproduct_zmtransfer( hD, &precond->UD, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
    CHECK( magma_minproduct_zmtransfer( hR, &precond->U, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
    
    magma_minproduct_zmfree(&hD, queue );
    magma_minproduct_zmfree(&hR, queue );
    magma_minproduct_zmfree(&hA, queue );
    magma_minproduct_zmfree(&hAt, queue );
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
    b           magma_minproduct_z_matrix
                RHS

    @param[in,out]
    x           magma_minproduct_z_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_minproduct_z_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zhepr
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zapplycumicc_l(
    magma_minproduct_z_matrix b,
    magma_minproduct_z_matrix *x,
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrL=NULL;
    
    magma_minproductDoubleComplex one = MAGMA_minproduct_Z_MAKE( 1.0, 0.0);

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrL ));
    CHECK_CUSPARSE( cusparseSetMatType( descrL, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrL, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrL, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrL, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseZcsrsm_solve( cusparseHandle,
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
    
    magma_minproduct_device_sync();

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
    b           magma_minproduct_z_matrix
                RHS

    @param[in,out]
    x           magma_minproduct_z_matrix*
                vector to precondition

    @param[in,out]
    precond     magma_minproduct_z_preconditioner*
                preconditioner parameters
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zhepr
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zapplycumicc_r(
    magma_minproduct_z_matrix b,
    magma_minproduct_z_matrix *x,
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrU=NULL;
    
    magma_minproductDoubleComplex one = MAGMA_minproduct_Z_MAKE( 1.0, 0.0);

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseZcsrsm_solve( cusparseHandle,
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
    
    magma_minproduct_device_sync();

cleanup:
    cusparseDestroyMatDescr( descrU );
    cusparseDestroy( cusparseHandle );
    return info; 
}









