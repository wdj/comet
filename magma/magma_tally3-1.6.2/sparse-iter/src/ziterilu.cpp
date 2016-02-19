/*
    -- MAGMA_tally3 (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @author Hartwig Anzt

       @precisions normal z -> s d c
*/
#include "common_magma_tally3sparse.h"

#define PRECISION_z


/**
    Purpose
    -------

    Prepares the ILU preconditioner via the iterative ILU iteration.

    Arguments
    ---------

    @param[in]
    A           magma_tally3_z_matrix
                input matrix A

    @param[in][out]
    precond     magma_tally3_z_preconditioner*
                preconditioner parameters
                
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zgepr
    ********************************************************************/
extern "C"
magma_tally3_int_t
magma_tally3_ziterilusetup(
    magma_tally3_z_matrix A,
    magma_tally3_z_matrix b,
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;

    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrL=NULL;
    cusparseMatDescr_t descrU=NULL;

    magma_tally3_z_matrix hAh={Magma_tally3_CSR}, hA={Magma_tally3_CSR}, hL={Magma_tally3_CSR}, hU={Magma_tally3_CSR},
    hAcopy={Magma_tally3_CSR}, hAL={Magma_tally3_CSR}, hAU={Magma_tally3_CSR}, hAUt={Magma_tally3_CSR},
    hUT={Magma_tally3_CSR}, hAtmp={Magma_tally3_CSR}, hACSRCOO={Magma_tally3_CSR}, dAinitguess={Magma_tally3_CSR},
    dL={Magma_tally3_CSR}, dU={Magma_tally3_CSR}, DL={Magma_tally3_CSR}, RL={Magma_tally3_CSR}, DU={Magma_tally3_CSR}, RU={Magma_tally3_CSR};

    // copy original matrix as CSRCOO to device
    CHECK( magma_tally3_zmtransfer(A, &hAh, A.memory_location, Magma_tally3_CPU, queue ));
    CHECK( magma_tally3_zmconvert( hAh, &hA, hAh.storage_type, Magma_tally3_CSR , queue ));
    magma_tally3_zmfree(&hAh, queue );

    CHECK( magma_tally3_zmtransfer( hA, &hAcopy, Magma_tally3_CPU, Magma_tally3_CPU , queue ));

    // in case using fill-in
    CHECK( magma_tally3_zsymbilu( &hAcopy, precond->levels, &hAL, &hAUt,  queue ));
    // add a unit diagonal to L for the algorithm
    CHECK( magma_tally3_zmLdiagadd( &hAL , queue ));
    // transpose U for the algorithm
    CHECK( magma_tally3_z_cucsrtranspose(  hAUt, &hAU , queue ));
    magma_tally3_zmfree( &hAUt , queue );

    // ---------------- initial guess ------------------- //
    CHECK( magma_tally3_zmconvert( hAcopy, &hACSRCOO, Magma_tally3_CSR, Magma_tally3_CSRCOO , queue ));
    CHECK( magma_tally3_zmtransfer( hACSRCOO, &dAinitguess, Magma_tally3_CPU, Magma_tally3_DEV , queue ));
    magma_tally3_zmfree(&hACSRCOO, queue );
    magma_tally3_zmfree(&hAcopy, queue );

    // transfer the factor L and U
    CHECK( magma_tally3_zmtransfer( hAL, &dL, Magma_tally3_CPU, Magma_tally3_DEV , queue ));
    CHECK( magma_tally3_zmtransfer( hAU, &dU, Magma_tally3_CPU, Magma_tally3_DEV , queue ));
    magma_tally3_zmfree(&hAL, queue );
    magma_tally3_zmfree(&hAU, queue );

    for(int i=0; i<precond->sweeps; i++){
        CHECK( magma_tally3_ziterilu_csr( dAinitguess, dL, dU , queue ));
    }

    CHECK( magma_tally3_zmtransfer( dL, &hL, Magma_tally3_DEV, Magma_tally3_CPU , queue ));
    CHECK( magma_tally3_zmtransfer( dU, &hU, Magma_tally3_DEV, Magma_tally3_CPU , queue ));
    CHECK( magma_tally3_z_cucsrtranspose(  hU, &hUT , queue ));

    magma_tally3_zmfree(&dL, queue );
    magma_tally3_zmfree(&dU, queue );
    magma_tally3_zmfree(&hU, queue );
    CHECK( magma_tally3_zmlumerge( hL, hUT, &hAtmp, queue ));

    magma_tally3_zmfree(&hL, queue );
    magma_tally3_zmfree(&hUT, queue );

    CHECK( magma_tally3_zmtransfer( hAtmp, &precond->M, Magma_tally3_CPU, Magma_tally3_DEV , queue ));

    hAL.diagorder_type = Magma_tally3_UNITY;
    CHECK( magma_tally3_zmconvert(hAtmp, &hAL, Magma_tally3_CSR, Magma_tally3_CSRL, queue ));
    hAL.storage_type = Magma_tally3_CSR;
    CHECK( magma_tally3_zmconvert(hAtmp, &hAU, Magma_tally3_CSR, Magma_tally3_CSRU, queue ));
    hAU.storage_type = Magma_tally3_CSR;

    magma_tally3_zmfree(&hAtmp, queue );

    CHECK( magma_tally3_zcsrsplit( 256, hAL, &DL, &RL , queue ));
    CHECK( magma_tally3_zcsrsplit( 256, hAU, &DU, &RU , queue ));

    CHECK( magma_tally3_zmtransfer( DL, &precond->LD, Magma_tally3_CPU, Magma_tally3_DEV , queue ));
    CHECK( magma_tally3_zmtransfer( DU, &precond->UD, Magma_tally3_CPU, Magma_tally3_DEV , queue ));

    // for cusparse uncomment this
    CHECK( magma_tally3_zmtransfer( hAL, &precond->L, Magma_tally3_CPU, Magma_tally3_DEV , queue ));
    CHECK( magma_tally3_zmtransfer( hAU, &precond->U, Magma_tally3_CPU, Magma_tally3_DEV , queue ));
    
/*

    //-- for ba-solve uncomment this

    if( RL.nnz != 0 )
        CHECK( magma_tally3_zmtransfer( RL, &precond->L, Magma_tally3_CPU, Magma_tally3_DEV , queue ));
    else{
        precond->L.nnz = 0;
        precond->L.val = NULL;
        precond->L.col = NULL;
        precond->L.row = NULL;
        precond->L.blockinfo = NULL;
    }

    if( RU.nnz != 0 )
        CHECK( magma_tally3_zmtransfer( RU, &precond->U, Magma_tally3_CPU, Magma_tally3_DEV , queue ));
    else{
        precond->U.nnz = 0;
        precond->L.val = NULL;
        precond->L.col = NULL;
        precond->L.row = NULL;
        precond->L.blockinfo = NULL;
    }

    //-- for ba-solve uncomment this
*/

        // extract the diagonal of L into precond->d
    CHECK( magma_tally3_zjacobisetup_diagscal( precond->L, &precond->d, queue ));
    CHECK( magma_tally3_zvinit( &precond->work1, Magma_tally3_DEV, hA.num_rows, 1, MAGMA_tally3_Z_ZERO, queue ));
    
    // extract the diagonal of U into precond->d2
    CHECK( magma_tally3_zjacobisetup_diagscal( precond->U, &precond->d2, queue ));
    CHECK( magma_tally3_zvinit( &precond->work2, Magma_tally3_DEV, hA.num_rows, 1, MAGMA_tally3_Z_ZERO, queue ));

    magma_tally3_zmfree(&hAL, queue );
    magma_tally3_zmfree(&hAU, queue );
    magma_tally3_zmfree(&DL, queue );
    magma_tally3_zmfree(&RL, queue );
    magma_tally3_zmfree(&DU, queue );
    magma_tally3_zmfree(&RU, queue );

    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrL ));
    CHECK_CUSPARSE( cusparseSetMatType( descrL, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrL, CUSPARSE_DIAG_TYPE_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrL, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrL, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseZcsrsv_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->L.num_rows,
        precond->L.nnz, descrL,
        precond->L.val, precond->L.row, precond->L.col, precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_UPPER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoU ));
    CHECK_CUSPARSE( cusparseZcsrsv_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->U.num_rows,
        precond->U.nnz, descrU,
        precond->U.val, precond->U.row, precond->U.col, precond->cuinfoU ));


cleanup:
    cusparseDestroy( cusparseHandle );
    cusparseDestroyMatDescr( descrL );
    cusparseDestroyMatDescr( descrU );
    cusparseHandle=NULL;
    descrL=NULL;
    descrU=NULL;
    magma_tally3_zmfree( &hAh, queue );
    magma_tally3_zmfree( &hA, queue );
    magma_tally3_zmfree( &hL, queue );
    magma_tally3_zmfree( &hU, queue );
    magma_tally3_zmfree( &hAcopy, queue );
    magma_tally3_zmfree( &hAL, queue );
    magma_tally3_zmfree( &hAU, queue );
    magma_tally3_zmfree( &hAUt, queue );
    magma_tally3_zmfree( &hUT, queue );
    magma_tally3_zmfree( &hAtmp, queue );
    magma_tally3_zmfree( &hACSRCOO, queue );
    magma_tally3_zmfree( &dAinitguess, queue );
    magma_tally3_zmfree( &dL, queue );
    magma_tally3_zmfree( &dU, queue );
    magma_tally3_zmfree( &DL, queue );
    magma_tally3_zmfree( &DU, queue );
    magma_tally3_zmfree( &RL, queue );
    magma_tally3_zmfree( &RU, queue );

    return info;
}






/**
    Purpose
    -------

    Prepares the IC preconditioner via the iterative IC iteration.

    Arguments
    ---------

    @param[in]
    A           magma_tally3_z_matrix
                input matrix A

    @param[in][out]
    precond     magma_tally3_z_preconditioner*
                preconditioner parameters
                
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zhepr
    ********************************************************************/
extern "C"
magma_tally3_int_t
magma_tally3_zitericsetup(
    magma_tally3_z_matrix A,
    magma_tally3_z_matrix b,
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    cusparseHandle_t cusparseHandle=NULL;
    cusparseMatDescr_t descrL=NULL;
    cusparseMatDescr_t descrU=NULL;

    magma_tally3_z_matrix hAh={Magma_tally3_CSR}, hA={Magma_tally3_CSR}, hAtmp={Magma_tally3_CSR},
    hAL={Magma_tally3_CSR}, hAUt={Magma_tally3_CSR}, hALt={Magma_tally3_CSR}, hM={Magma_tally3_CSR},
    hACSRCOO={Magma_tally3_CSR}, dAinitguess={Magma_tally3_CSR}, dL={Magma_tally3_CSR};
    magma_tally3_z_matrix d_h={Magma_tally3_CSR};


    // copy original matrix as CSRCOO to device
    CHECK( magma_tally3_zmtransfer(A, &hAh, A.memory_location, Magma_tally3_CPU, queue ));
    CHECK( magma_tally3_zmconvert( hAh, &hA, hAh.storage_type, Magma_tally3_CSR , queue ));
    magma_tally3_zmfree(&hAh, queue );

    // in case using fill-in
    CHECK( magma_tally3_zsymbilu( &hA, precond->levels, &hAL, &hAUt , queue ));

    // need only lower triangular
    magma_tally3_zmfree(&hAUt, queue );
    magma_tally3_zmfree(&hAL, queue );
    CHECK( magma_tally3_zmconvert( hA, &hAtmp, Magma_tally3_CSR, Magma_tally3_CSRL , queue ));
    magma_tally3_zmfree(&hA, queue );

    // ---------------- initial guess ------------------- //
    CHECK( magma_tally3_zmconvert( hAtmp, &hACSRCOO, Magma_tally3_CSR, Magma_tally3_CSRCOO , queue ));
    //int blocksize = 1;
    //magma_tally3_zmreorder( hACSRCOO, n, blocksize, blocksize, blocksize, &hAinitguess , queue );
    CHECK( magma_tally3_zmtransfer( hACSRCOO, &dAinitguess, Magma_tally3_CPU, Magma_tally3_DEV , queue ));
    magma_tally3_zmfree(&hACSRCOO, queue );
    CHECK( magma_tally3_zmtransfer( hAtmp, &dL, Magma_tally3_CPU, Magma_tally3_DEV , queue ));
    magma_tally3_zmfree(&hAtmp, queue );

    for(int i=0; i<precond->sweeps; i++){
        CHECK( magma_tally3_ziteric_csr( dAinitguess, dL , queue ));
    }
    CHECK( magma_tally3_zmtransfer( dL, &hAL, Magma_tally3_DEV, Magma_tally3_CPU , queue ));
    magma_tally3_zmfree(&dL, queue );
    magma_tally3_zmfree(&dAinitguess, queue );


    // for CUSPARSE
    CHECK( magma_tally3_zmtransfer( hAL, &precond->M, Magma_tally3_CPU, Magma_tally3_DEV , queue ));

    // Jacobi setup
    CHECK( magma_tally3_zjacobisetup_matrix( precond->M, &precond->L, &precond->d , queue ));

    // for Jacobi, we also need U
    CHECK( magma_tally3_z_cucsrtranspose(   hAL, &hALt , queue ));
    CHECK( magma_tally3_zjacobisetup_matrix( hALt, &hM, &d_h , queue ));

    CHECK( magma_tally3_zmtransfer( hM, &precond->U, Magma_tally3_CPU, Magma_tally3_DEV , queue ));

    magma_tally3_zmfree(&hM, queue );

    magma_tally3_zmfree(&d_h, queue );


        // copy the matrix to precond->L and (transposed) to precond->U
    CHECK( magma_tally3_zmtransfer(precond->M, &(precond->L), Magma_tally3_DEV, Magma_tally3_DEV, queue ));
    CHECK( magma_tally3_zmtranspose( precond->L, &(precond->U), queue ));

    // extract the diagonal of L into precond->d
    CHECK( magma_tally3_zjacobisetup_diagscal( precond->L, &precond->d, queue ));
    CHECK( magma_tally3_zvinit( &precond->work1, Magma_tally3_DEV, hAL.num_rows, 1, MAGMA_tally3_Z_ZERO, queue ));

    // extract the diagonal of U into precond->d2
    CHECK( magma_tally3_zjacobisetup_diagscal( precond->U, &precond->d2, queue ));
    CHECK( magma_tally3_zvinit( &precond->work2, Magma_tally3_DEV, hAL.num_rows, 1, MAGMA_tally3_Z_ZERO, queue ));


    magma_tally3_zmfree(&hAL, queue );
    magma_tally3_zmfree(&hALt, queue );


    // CUSPARSE context //
    CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrL ));
    CHECK_CUSPARSE( cusparseSetMatType( descrL, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrL, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrL, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrL, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseZcsrsv_analysis( cusparseHandle,
        CUSPARSE_OPERATION_NON_TRANSPOSE, precond->M.num_rows,
        precond->M.nnz, descrL,
        precond->M.val, precond->M.row, precond->M.col, precond->cuinfoL ));
    CHECK_CUSPARSE( cusparseCreateMatDescr( &descrU ));
    CHECK_CUSPARSE( cusparseSetMatType( descrU, CUSPARSE_MATRIX_TYPE_TRIANGULAR ));
    CHECK_CUSPARSE( cusparseSetMatDiagType( descrU, CUSPARSE_DIAG_TYPE_NON_UNIT ));
    CHECK_CUSPARSE( cusparseSetMatIndexBase( descrU, CUSPARSE_INDEX_BASE_ZERO ));
    CHECK_CUSPARSE( cusparseSetMatFillMode( descrU, CUSPARSE_FILL_MODE_LOWER ));
    CHECK_CUSPARSE( cusparseCreateSolveAnalysisInfo( &precond->cuinfoU ));
    CHECK_CUSPARSE( cusparseZcsrsv_analysis( cusparseHandle,
        CUSPARSE_OPERATION_TRANSPOSE, precond->M.num_rows,
        precond->M.nnz, descrU,
        precond->M.val, precond->M.row, precond->M.col, precond->cuinfoU ));

    
    cleanup:
    cusparseDestroy( cusparseHandle );
    cusparseDestroyMatDescr( descrL );
    cusparseDestroyMatDescr( descrU );
    cusparseHandle=NULL;
    descrL=NULL;
    descrU=NULL;    
    magma_tally3_zmfree( &hAh, queue );
    magma_tally3_zmfree( &hA, queue );
    magma_tally3_zmfree( &hAtmp, queue );
    magma_tally3_zmfree( &hAL, queue );
    magma_tally3_zmfree( &hAUt, queue );
    magma_tally3_zmfree( &hALt, queue );
    magma_tally3_zmfree( &hM, queue );
    magma_tally3_zmfree( &hACSRCOO, queue );
    magma_tally3_zmfree( &dAinitguess, queue );
    magma_tally3_zmfree( &dL, queue );
    magma_tally3_zmfree( &d_h, queue );
    
    return info;
}


/**
    Purpose
    -------

    Updates an existing preconditioner via additional iterative IC sweeps for
    previous factorization initial guess (PFIG).
    See  Anzt et al., Parallel Computing, 2015.

    Arguments
    ---------

    @param[in]
    precond         magma_tally3_z_preconditioner*
                    preconditioner parameters

    @param[in]
    magma_tally3_int_t     number of updates
    
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.
                
    @ingroup magma_tally3sparse_zhepr
    ********************************************************************/
extern "C"
magma_tally3_int_t
magma_tally3_zitericupdate(
    magma_tally3_z_matrix A,
    magma_tally3_z_preconditioner *precond,
    magma_tally3_int_t updates,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;

    magma_tally3_z_matrix hALt={Magma_tally3_CSR};
    magma_tally3_z_matrix d_h={Magma_tally3_CSR};


    // copy original matrix as CSRCOO to device

    for(int i=0; i<updates; i++){
        CHECK( magma_tally3_ziteric_csr( A, precond->M , queue ));
    }
    //magma_tally3_zmtransfer( precond->M, &precond->M, Magma_tally3_DEV, Magma_tally3_DEV , queue );
    magma_tally3_zmfree(&precond->L, queue );
    magma_tally3_zmfree(&precond->U, queue );
    magma_tally3_zmfree( &precond->d , queue );


    // Jacobi setup
    CHECK( magma_tally3_zjacobisetup_matrix( precond->M, &precond->L, &precond->d , queue ));

    // for Jacobi, we also need U
    CHECK( magma_tally3_z_cucsrtranspose(   precond->M, &hALt , queue ));
    CHECK( magma_tally3_zjacobisetup_matrix( hALt, &precond->U, &d_h , queue ));


cleanup:
    magma_tally3_zmfree(&d_h, queue );
    magma_tally3_zmfree(&hALt, queue );
    
    return info;
}


/**
    Purpose
    -------

    Performs the left triangular solves using the IC preconditioner via Jacobi.

    Arguments
    ---------

    @param[in]
    b           magma_tally3_z_matrix
                RHS

    @param[out]
    x           magma_tally3_z_matrix*
                vector to precondition

    @param[in]
    precond     magma_tally3_z_preconditioner*
                preconditioner parameters
                
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zgepr
    ********************************************************************/
extern "C"
magma_tally3_int_t
magma_tally3_zapplyiteric_l(
    magma_tally3_z_matrix b,
    magma_tally3_z_matrix *x,
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    magma_tally3_int_t dofs = precond->L.num_rows;
    magma_tally3_z_solver_par jacobiiter_par;
    jacobiiter_par.maxiter = precond->maxiter;

    // compute c = D^{-1}b and copy c as initial guess to x
    CHECK( magma_tally3_zjacobisetup_vector_gpu( dofs, b, precond->d,
                                                precond->work1, x, queue ));
    // Jacobi iterator
    CHECK( magma_tally3_zjacobiiter_precond( precond->L, x, &jacobiiter_par, precond , queue ));

cleanup:
    return info;
}


/**
    Purpose
    -------

    Performs the right triangular solves using the IC preconditioner via Jacobi.

    Arguments
    ---------

    @param[in]
    b           magma_tally3_z_matrix
                RHS

    @param[out]
    x           magma_tally3_z_matrix*
                vector to precondition

    @param[in]
    precond     magma_tally3_z_preconditioner*
                preconditioner parameters
                
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zgepr
    ********************************************************************/
extern "C"
magma_tally3_int_t
magma_tally3_zapplyiteric_r(
    magma_tally3_z_matrix b,
    magma_tally3_z_matrix *x,
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;

    magma_tally3_int_t dofs = precond->U.num_rows;
    magma_tally3_z_solver_par jacobiiter_par;
    jacobiiter_par.maxiter = precond->maxiter;

    // compute c = D^{-1}b and copy c as initial guess to x
    CHECK( magma_tally3_zjacobisetup_vector_gpu( dofs, b, precond->d,
                                                precond->work1, x, queue ));

    // Jacobi iterator
    CHECK( magma_tally3_zjacobiiter_precond( precond->U, x, &jacobiiter_par, precond , queue ));
    
cleanup:
    return info;
}
