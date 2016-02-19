/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally3_z_blaswrapper.cpp normal z -> c, Sun May  3 11:22:57 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally3sparse.h"


#include "common_magma_tally3sparse.h"
#include "magma_tally3blas.h"
#include "magma_tally3sparse_types.h"




/**
    Purpose
    -------

    For a given input matrix A and vectors x, y and scalars alpha, beta
    the wrapper determines the suitable SpMV computing
              y = alpha * A * x + beta * y.
    Arguments
    ---------

    @param[in]
    alpha       magma_tally3FloatComplex
                scalar alpha

    @param[in]
    A           magma_tally3_c_matrix
                sparse matrix A

    @param[in]
    x           magma_tally3_c_matrix
                input vector x
                
    @param[in]
    beta        magma_tally3FloatComplex
                scalar beta
    @param[out]
    y           magma_tally3_c_matrix
                output vector y
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_c
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_c_spmv(
    magma_tally3FloatComplex alpha,
    magma_tally3_c_matrix A,
    magma_tally3_c_matrix x,
    magma_tally3FloatComplex beta,
    magma_tally3_c_matrix y,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue=NULL;
    magma_tally3blasGetKernelStream( &orig_queue );

    magma_tally3_c_matrix x2={Magma_tally3_CSR};

    cusparseHandle_t cusparseHandle = 0;
    cusparseMatDescr_t descr = 0;
    // make sure RHS is a dense matrix
    if ( x.storage_type != Magma_tally3_DENSE ) {
         printf("error: only dense vectors are supported.\n");
         info = MAGMA_tally3_ERR_NOT_SUPPORTED;
         goto cleanup;
    }

    if ( A.memory_location != x.memory_location ||
                            x.memory_location != y.memory_location ) {
        printf("error: linear algebra objects are not located in same memory!\n");
        printf("memory locations are: %d   %d   %d\n",
                        A.memory_location, x.memory_location, y.memory_location );
        info = MAGMA_tally3_ERR_INVALID_PTR;
        goto cleanup;
    }

    // DEV case
    if ( A.memory_location == Magma_tally3_DEV ) {
        if ( A.num_cols == x.num_rows && x.num_cols == 1 ) {

             if ( A.storage_type == Magma_tally3_CSR
                            || A.storage_type == Magma_tally3_CSRL
                            || A.storage_type == Magma_tally3_CSRU ) {
                 //printf("using CSR kernel for SpMV: ");
                 //magma_tally3_cgecsrmv( Magma_tally3NoTrans, A.num_rows, A.num_cols, alpha,
                 //                A.dval, A.drow, A.dcol, x.dval, beta, y.dval );
                 //printf("done.\n");
                CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
                CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
                CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));

                CHECK_CUSPARSE( cusparseSetMatType( descr, CUSPARSE_MATRIX_TYPE_GENERAL ));
                CHECK_CUSPARSE( cusparseSetMatIndexBase( descr, CUSPARSE_INDEX_BASE_ZERO ));

                cusparseCcsrmv( cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                            A.num_rows, A.num_cols, A.nnz, &alpha, descr,
                            A.dval, A.drow, A.dcol, x.dval, &beta, y.dval );
             }
             else if ( A.storage_type == Magma_tally3_ELL ) {
                 //printf("using ELLPACKT kernel for SpMV: ");
                 CHECK( magma_tally3_cgeelltmv( Magma_tally3NoTrans, A.num_rows, A.num_cols,
                    A.max_nnz_row, alpha, A.dval, A.dcol, x.dval, beta,
                    y.dval, queue ));
                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_tally3_ELLPACKT ) {
                 //printf("using ELL kernel for SpMV: ");
                 CHECK( magma_tally3_cgeellmv( Magma_tally3NoTrans, A.num_rows, A.num_cols,
                    A.max_nnz_row, alpha, A.dval, A.dcol, x.dval, beta,
                    y.dval, queue ));
                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_tally3_ELLRT ) {
                 //printf("using ELLRT kernel for SpMV: ");
                 CHECK( magma_tally3_cgeellrtmv( Magma_tally3NoTrans, A.num_rows, A.num_cols,
                            A.max_nnz_row, alpha, A.dval, A.dcol, A.drow, x.dval,
                         beta, y.dval, A.alignment, A.blocksize, queue ));
                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_tally3_SELLP ) {
                 //printf("using SELLP kernel for SpMV: ");
                 CHECK( magma_tally3_cgesellpmv( Magma_tally3NoTrans, A.num_rows, A.num_cols,
                    A.blocksize, A.numblocks, A.alignment,
                    alpha, A.dval, A.dcol, A.drow, x.dval, beta, y.dval, queue ));

                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_tally3_DENSE ) {
                 //printf("using DENSE kernel for SpMV: ");
                 magma_tally3blas_cgemv( Magma_tally3NoTrans, A.num_rows, A.num_cols, alpha,
                                A.dval, A.num_rows, x.dval, 1, beta,  y.dval,
                                1 );
                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_tally3_SPMVFUNCTION ) {
                 //printf("using DENSE kernel for SpMV: ");
                 CHECK( magma_tally3_ccustomspmv( alpha, x, beta, y, queue ));
                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_tally3_BCSR ) {
                 //printf("using CUSPARSE BCSR kernel for SpMV: ");
                // CUSPARSE context //
                cusparseDirection_t dirA = CUSPARSE_DIRECTION_ROW;
                int mb = (A.num_rows + A.blocksize-1)/A.blocksize;
                int nb = (A.num_cols + A.blocksize-1)/A.blocksize;
                CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
                CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
                CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));
                cusparseCbsrmv( cusparseHandle, dirA,
                    CUSPARSE_OPERATION_NON_TRANSPOSE, mb, nb, A.numblocks,
                    &alpha, descr, A.dval, A.drow, A.dcol, A.blocksize, x.dval,
                    &beta, y.dval );

             }
             else {
                 printf("error: format not supported.\n");
                 info = MAGMA_tally3_ERR_NOT_SUPPORTED; 
             }
        }
        else if ( A.num_cols < x.num_rows || x.num_cols > 1 ) {
            magma_tally3_int_t num_vecs = x.num_rows / A.num_cols * x.num_cols;
            if ( A.storage_type == Magma_tally3_CSR ) {

                CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
                CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
                CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));
                CHECK_CUSPARSE( cusparseSetMatType( descr, CUSPARSE_MATRIX_TYPE_GENERAL ));
                CHECK_CUSPARSE( cusparseSetMatIndexBase( descr, CUSPARSE_INDEX_BASE_ZERO ));

                if ( x.major == Magma_tally3ColMajor) {
                    cusparseCcsrmm(cusparseHandle,
                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                    A.num_rows,   num_vecs, A.num_cols, A.nnz,
                    &alpha, descr, A.dval, A.drow, A.dcol,
                    x.dval, A.num_cols, &beta, y.dval, A.num_cols);
                } else if ( x.major == Magma_tally3RowMajor) {
                    cusparseCcsrmm2(cusparseHandle,
                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                    CUSPARSE_OPERATION_TRANSPOSE,
                    A.num_rows,   num_vecs, A.num_cols, A.nnz,
                    &alpha, descr, A.dval, A.drow, A.dcol,
                    x.dval, A.num_cols, &beta, y.dval, A.num_cols);
                }

             } else if ( A.storage_type == Magma_tally3_SELLP ) {
                if ( x.major == Magma_tally3RowMajor) {
                 CHECK( magma_tally3_cmgesellpmv( Magma_tally3NoTrans, A.num_rows, A.num_cols,
                    num_vecs, A.blocksize, A.numblocks, A.alignment,
                    alpha, A.dval, A.dcol, A.drow, x.dval, beta, y.dval, queue ));
                }
                else if ( x.major == Magma_tally3ColMajor) {
                    // transpose first to row major
                    CHECK( magma_tally3_cvtranspose( x, &x2, queue ));
                    CHECK( magma_tally3_cmgesellpmv( Magma_tally3NoTrans, A.num_rows, A.num_cols,
                    num_vecs, A.blocksize, A.numblocks, A.alignment,
                    alpha, A.dval, A.dcol, A.drow, x2.dval, beta, y.dval, queue ));
                }
             }
             /*if ( A.storage_type == Magma_tally3_DENSE ) {
                 //printf("using DENSE kernel for SpMV: ");
                 magma_tally3blas_cmgemv( Magma_tally3NoTrans, A.num_rows, A.num_cols,
                            num_vecs, alpha, A.dval, A.num_rows, x.dval, 1,
                            beta,  y.dval, 1 );
                 //printf("done.\n");
             }*/
             else {
                 printf("error: format not supported.\n");
                 info = MAGMA_tally3_ERR_NOT_SUPPORTED;
             }
        }
         
         
    }
    // CPU case missing!
    else {
        printf("error: CPU not yet supported.\n");
        info = MAGMA_tally3_ERR_NOT_SUPPORTED;
    }

cleanup:
    cusparseDestroyMatDescr( descr );
    cusparseDestroy( cusparseHandle );
    cusparseHandle = 0;
    descr = 0;
    magma_tally3_cmfree(&x2, queue );
    
    magma_tally3blasSetKernelStream( orig_queue );
    return info;


}





/**
    Purpose
    -------

    For a given input matrix A and vectors x, y and scalars alpha, beta
    the wrapper determines the suitable SpMV computing
              y = alpha * ( A - lambda I ) * x + beta * y.
    Arguments
    ---------

    @param
    alpha       magma_tally3FloatComplex
                scalar alpha

    @param
    A           magma_tally3_c_matrix
                sparse matrix A

    @param
    lambda      magma_tally3FloatComplex
                scalar lambda

    @param
    x           magma_tally3_c_matrix
                input vector x

    @param
    beta        magma_tally3FloatComplex
                scalar beta
                
    @param
    offset      magma_tally3_int_t
                in case not the main diagonal is scaled
                
    @param
    blocksize   magma_tally3_int_t
                in case of processing multiple vectors
                
    @param
    add_rows    magma_tally3_int_t*
                in case the matrixpowerskernel is used
                
    @param
    y           magma_tally3_c_matrix
                output vector y
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_caux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_c_spmv_shift(
    magma_tally3FloatComplex alpha,
    magma_tally3_c_matrix A,
    magma_tally3FloatComplex lambda,
    magma_tally3_c_matrix x,
    magma_tally3FloatComplex beta,
    magma_tally3_int_t offset,
    magma_tally3_int_t blocksize,
    magma_tally3_index_t *add_rows,
    magma_tally3_c_matrix y,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;

    // make sure RHS is a dense matrix
    if ( x.storage_type != Magma_tally3_DENSE ) {
         printf("error: only dense vectors are supported.\n");
         info = MAGMA_tally3_ERR_NOT_SUPPORTED;
         goto cleanup;
    }


    if ( A.memory_location != x.memory_location
                || x.memory_location != y.memory_location ) {
        printf("error: linear algebra objects are not located in same memory!\n");
        printf("memory locations are: %d   %d   %d\n",
                    A.memory_location, x.memory_location, y.memory_location );
        info = MAGMA_tally3_ERR_INVALID_PTR;
        goto cleanup;
    }

    // DEV case
    if ( A.memory_location == Magma_tally3_DEV ) {
         if ( A.storage_type == Magma_tally3_CSR ) {
             //printf("using CSR kernel for SpMV: ");
             CHECK( magma_tally3_cgecsrmv_shift( Magma_tally3NoTrans, A.num_rows, A.num_cols,
                alpha, lambda, A.dval, A.drow, A.dcol, x.dval, beta, offset,
                blocksize, add_rows, y.dval, queue ));
             //printf("done.\n");
         }
         else if ( A.storage_type == Magma_tally3_ELLPACKT ) {
             //printf("using ELLPACKT kernel for SpMV: ");
             CHECK( magma_tally3_cgeellmv_shift( Magma_tally3NoTrans, A.num_rows, A.num_cols,
                A.max_nnz_row, alpha, lambda, A.dval, A.dcol, x.dval, beta, offset,
                blocksize, add_rows, y.dval, queue ));
             //printf("done.\n");
         }
         else if ( A.storage_type == Magma_tally3_ELL ) {
             //printf("using ELL kernel for SpMV: ");
             CHECK( magma_tally3_cgeelltmv_shift( Magma_tally3NoTrans, A.num_rows, A.num_cols,
                A.max_nnz_row, alpha, lambda, A.dval, A.dcol, x.dval, beta, offset,
                blocksize, add_rows, y.dval, queue ));
             //printf("done.\n");
         }
         else {
             printf("error: format not supported.\n");
             info = MAGMA_tally3_ERR_NOT_SUPPORTED;
         }
    }
    // CPU case missing!
    else {
        printf("error: CPU not yet supported.\n");
        info = MAGMA_tally3_ERR_NOT_SUPPORTED;
    }
cleanup:
    return info;
}



/**
    Purpose
    -------

    For a given input matrix A and B and scalar alpha,
    the wrapper determines the suitable SpMV computing
              C = alpha * A * B.
    Arguments
    ---------

    @param[in]
    alpha       magma_tally3FloatComplex
                scalar alpha

    @param[in]
    A           magma_tally3_c_matrix
                sparse matrix A
                
    @param[in]
    B           magma_tally3_c_matrix
                sparse matrix C
                
    @param[out]
    C           magma_tally3_c_matrix *
                outpur sparse matrix C

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_c
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_c_spmm(
    magma_tally3FloatComplex alpha,
    magma_tally3_c_matrix A,
    magma_tally3_c_matrix B,
    magma_tally3_c_matrix *C,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    // set queue for old dense routines
    //magma_tally3_queue_t orig_queue=NULL;
    //magma_tally3blasGetKernelStream( &orig_queue );

    if ( A.memory_location != B.memory_location ) {
        printf("error: linear algebra objects are not located in same memory!\n");
        printf("memory locations are: %d   %d\n",
                        A.memory_location, B.memory_location );
        info = MAGMA_tally3_ERR_INVALID_PTR;
        goto cleanup;
    }

    // DEV case
    if ( A.memory_location == Magma_tally3_DEV ) {
        if ( A.num_cols == B.num_rows ) {

             if ( A.storage_type == Magma_tally3_CSR
                            || A.storage_type == Magma_tally3_CSRL
                            || A.storage_type == Magma_tally3_CSRU
                            || A.storage_type == Magma_tally3_CSRCOO ) {
                    
                CHECK( magma_tally3_ccuspmm( A, B, C, queue ));
                
             }
             else {
                 printf("error: format not supported.\n");
                 // magma_tally3blasSetKernelStream( orig_queue );
                 info = MAGMA_tally3_ERR_NOT_SUPPORTED;
             }
        }
         
    }
    // CPU case missing!
    else {
        printf("error: CPU not yet supported.\n");
        // magma_tally3blasSetKernelStream( orig_queue );
        info = MAGMA_tally3_ERR_NOT_SUPPORTED; // TODO change to goto cleanup?
    }
    
cleanup:
    // magma_tally3blasSetKernelStream( orig_queue );
    return info;
}





