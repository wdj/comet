/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally2_z_blaswrapper.cpp normal z -> d, Sun May  3 11:22:57 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally2sparse.h"


#include "common_magma_tally2sparse.h"
#include "magma_tally2blas.h"
#include "magma_tally2sparse_types.h"




/**
    Purpose
    -------

    For a given input matrix A and vectors x, y and scalars alpha, beta
    the wrapper determines the suitable SpMV computing
              y = alpha * A * x + beta * y.
    Arguments
    ---------

    @param[in]
    alpha       double
                scalar alpha

    @param[in]
    A           magma_tally2_d_matrix
                sparse matrix A

    @param[in]
    x           magma_tally2_d_matrix
                input vector x
                
    @param[in]
    beta        double
                scalar beta
    @param[out]
    y           magma_tally2_d_matrix
                output vector y
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_d
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_d_spmv(
    double alpha,
    magma_tally2_d_matrix A,
    magma_tally2_d_matrix x,
    double beta,
    magma_tally2_d_matrix y,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    // set queue for old dense routines
    magma_tally2_queue_t orig_queue=NULL;
    magma_tally2blasGetKernelStream( &orig_queue );

    magma_tally2_d_matrix x2={Magma_tally2_CSR};

    cusparseHandle_t cusparseHandle = 0;
    cusparseMatDescr_t descr = 0;
    // make sure RHS is a dense matrix
    if ( x.storage_type != Magma_tally2_DENSE ) {
         printf("error: only dense vectors are supported.\n");
         info = MAGMA_tally2_ERR_NOT_SUPPORTED;
         goto cleanup;
    }

    if ( A.memory_location != x.memory_location ||
                            x.memory_location != y.memory_location ) {
        printf("error: linear algebra objects are not located in same memory!\n");
        printf("memory locations are: %d   %d   %d\n",
                        A.memory_location, x.memory_location, y.memory_location );
        info = MAGMA_tally2_ERR_INVALID_PTR;
        goto cleanup;
    }

    // DEV case
    if ( A.memory_location == Magma_tally2_DEV ) {
        if ( A.num_cols == x.num_rows && x.num_cols == 1 ) {

             if ( A.storage_type == Magma_tally2_CSR
                            || A.storage_type == Magma_tally2_CSRL
                            || A.storage_type == Magma_tally2_CSRU ) {
                 //printf("using CSR kernel for SpMV: ");
                 //magma_tally2_dgecsrmv( Magma_tally2NoTrans, A.num_rows, A.num_cols, alpha,
                 //                A.dval, A.drow, A.dcol, x.dval, beta, y.dval );
                 //printf("done.\n");
                CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
                CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
                CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));

                CHECK_CUSPARSE( cusparseSetMatType( descr, CUSPARSE_MATRIX_TYPE_GENERAL ));
                CHECK_CUSPARSE( cusparseSetMatIndexBase( descr, CUSPARSE_INDEX_BASE_ZERO ));

                cusparseDcsrmv( cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                            A.num_rows, A.num_cols, A.nnz, &alpha, descr,
                            A.dval, A.drow, A.dcol, x.dval, &beta, y.dval );
             }
             else if ( A.storage_type == Magma_tally2_ELL ) {
                 //printf("using ELLPACKT kernel for SpMV: ");
                 CHECK( magma_tally2_dgeelltmv( Magma_tally2NoTrans, A.num_rows, A.num_cols,
                    A.max_nnz_row, alpha, A.dval, A.dcol, x.dval, beta,
                    y.dval, queue ));
                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_tally2_ELLPACKT ) {
                 //printf("using ELL kernel for SpMV: ");
                 CHECK( magma_tally2_dgeellmv( Magma_tally2NoTrans, A.num_rows, A.num_cols,
                    A.max_nnz_row, alpha, A.dval, A.dcol, x.dval, beta,
                    y.dval, queue ));
                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_tally2_ELLRT ) {
                 //printf("using ELLRT kernel for SpMV: ");
                 CHECK( magma_tally2_dgeellrtmv( Magma_tally2NoTrans, A.num_rows, A.num_cols,
                            A.max_nnz_row, alpha, A.dval, A.dcol, A.drow, x.dval,
                         beta, y.dval, A.alignment, A.blocksize, queue ));
                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_tally2_SELLP ) {
                 //printf("using SELLP kernel for SpMV: ");
                 CHECK( magma_tally2_dgesellpmv( Magma_tally2NoTrans, A.num_rows, A.num_cols,
                    A.blocksize, A.numblocks, A.alignment,
                    alpha, A.dval, A.dcol, A.drow, x.dval, beta, y.dval, queue ));

                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_tally2_DENSE ) {
                 //printf("using DENSE kernel for SpMV: ");
                 magma_tally2blas_dgemv( Magma_tally2NoTrans, A.num_rows, A.num_cols, alpha,
                                A.dval, A.num_rows, x.dval, 1, beta,  y.dval,
                                1 );
                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_tally2_SPMVFUNCTION ) {
                 //printf("using DENSE kernel for SpMV: ");
                 CHECK( magma_tally2_dcustomspmv( alpha, x, beta, y, queue ));
                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_tally2_BCSR ) {
                 //printf("using CUSPARSE BCSR kernel for SpMV: ");
                // CUSPARSE context //
                cusparseDirection_t dirA = CUSPARSE_DIRECTION_ROW;
                int mb = (A.num_rows + A.blocksize-1)/A.blocksize;
                int nb = (A.num_cols + A.blocksize-1)/A.blocksize;
                CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
                CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
                CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));
                cusparseDbsrmv( cusparseHandle, dirA,
                    CUSPARSE_OPERATION_NON_TRANSPOSE, mb, nb, A.numblocks,
                    &alpha, descr, A.dval, A.drow, A.dcol, A.blocksize, x.dval,
                    &beta, y.dval );

             }
             else {
                 printf("error: format not supported.\n");
                 info = MAGMA_tally2_ERR_NOT_SUPPORTED; 
             }
        }
        else if ( A.num_cols < x.num_rows || x.num_cols > 1 ) {
            magma_tally2_int_t num_vecs = x.num_rows / A.num_cols * x.num_cols;
            if ( A.storage_type == Magma_tally2_CSR ) {

                CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
                CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
                CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));
                CHECK_CUSPARSE( cusparseSetMatType( descr, CUSPARSE_MATRIX_TYPE_GENERAL ));
                CHECK_CUSPARSE( cusparseSetMatIndexBase( descr, CUSPARSE_INDEX_BASE_ZERO ));

                if ( x.major == Magma_tally2ColMajor) {
                    cusparseDcsrmm(cusparseHandle,
                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                    A.num_rows,   num_vecs, A.num_cols, A.nnz,
                    &alpha, descr, A.dval, A.drow, A.dcol,
                    x.dval, A.num_cols, &beta, y.dval, A.num_cols);
                } else if ( x.major == Magma_tally2RowMajor) {
                    cusparseDcsrmm2(cusparseHandle,
                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                    CUSPARSE_OPERATION_TRANSPOSE,
                    A.num_rows,   num_vecs, A.num_cols, A.nnz,
                    &alpha, descr, A.dval, A.drow, A.dcol,
                    x.dval, A.num_cols, &beta, y.dval, A.num_cols);
                }

             } else if ( A.storage_type == Magma_tally2_SELLP ) {
                if ( x.major == Magma_tally2RowMajor) {
                 CHECK( magma_tally2_dmgesellpmv( Magma_tally2NoTrans, A.num_rows, A.num_cols,
                    num_vecs, A.blocksize, A.numblocks, A.alignment,
                    alpha, A.dval, A.dcol, A.drow, x.dval, beta, y.dval, queue ));
                }
                else if ( x.major == Magma_tally2ColMajor) {
                    // transpose first to row major
                    CHECK( magma_tally2_dvtranspose( x, &x2, queue ));
                    CHECK( magma_tally2_dmgesellpmv( Magma_tally2NoTrans, A.num_rows, A.num_cols,
                    num_vecs, A.blocksize, A.numblocks, A.alignment,
                    alpha, A.dval, A.dcol, A.drow, x2.dval, beta, y.dval, queue ));
                }
             }
             /*if ( A.storage_type == Magma_tally2_DENSE ) {
                 //printf("using DENSE kernel for SpMV: ");
                 magma_tally2blas_dmgemv( Magma_tally2NoTrans, A.num_rows, A.num_cols,
                            num_vecs, alpha, A.dval, A.num_rows, x.dval, 1,
                            beta,  y.dval, 1 );
                 //printf("done.\n");
             }*/
             else {
                 printf("error: format not supported.\n");
                 info = MAGMA_tally2_ERR_NOT_SUPPORTED;
             }
        }
         
         
    }
    // CPU case missing!
    else {
        printf("error: CPU not yet supported.\n");
        info = MAGMA_tally2_ERR_NOT_SUPPORTED;
    }

cleanup:
    cusparseDestroyMatDescr( descr );
    cusparseDestroy( cusparseHandle );
    cusparseHandle = 0;
    descr = 0;
    magma_tally2_dmfree(&x2, queue );
    
    magma_tally2blasSetKernelStream( orig_queue );
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
    alpha       double
                scalar alpha

    @param
    A           magma_tally2_d_matrix
                sparse matrix A

    @param
    lambda      double
                scalar lambda

    @param
    x           magma_tally2_d_matrix
                input vector x

    @param
    beta        double
                scalar beta
                
    @param
    offset      magma_tally2_int_t
                in case not the main diagonal is scaled
                
    @param
    blocksize   magma_tally2_int_t
                in case of processing multiple vectors
                
    @param
    add_rows    magma_tally2_int_t*
                in case the matrixpowerskernel is used
                
    @param
    y           magma_tally2_d_matrix
                output vector y
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_daux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_d_spmv_shift(
    double alpha,
    magma_tally2_d_matrix A,
    double lambda,
    magma_tally2_d_matrix x,
    double beta,
    magma_tally2_int_t offset,
    magma_tally2_int_t blocksize,
    magma_tally2_index_t *add_rows,
    magma_tally2_d_matrix y,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;

    // make sure RHS is a dense matrix
    if ( x.storage_type != Magma_tally2_DENSE ) {
         printf("error: only dense vectors are supported.\n");
         info = MAGMA_tally2_ERR_NOT_SUPPORTED;
         goto cleanup;
    }


    if ( A.memory_location != x.memory_location
                || x.memory_location != y.memory_location ) {
        printf("error: linear algebra objects are not located in same memory!\n");
        printf("memory locations are: %d   %d   %d\n",
                    A.memory_location, x.memory_location, y.memory_location );
        info = MAGMA_tally2_ERR_INVALID_PTR;
        goto cleanup;
    }

    // DEV case
    if ( A.memory_location == Magma_tally2_DEV ) {
         if ( A.storage_type == Magma_tally2_CSR ) {
             //printf("using CSR kernel for SpMV: ");
             CHECK( magma_tally2_dgecsrmv_shift( Magma_tally2NoTrans, A.num_rows, A.num_cols,
                alpha, lambda, A.dval, A.drow, A.dcol, x.dval, beta, offset,
                blocksize, add_rows, y.dval, queue ));
             //printf("done.\n");
         }
         else if ( A.storage_type == Magma_tally2_ELLPACKT ) {
             //printf("using ELLPACKT kernel for SpMV: ");
             CHECK( magma_tally2_dgeellmv_shift( Magma_tally2NoTrans, A.num_rows, A.num_cols,
                A.max_nnz_row, alpha, lambda, A.dval, A.dcol, x.dval, beta, offset,
                blocksize, add_rows, y.dval, queue ));
             //printf("done.\n");
         }
         else if ( A.storage_type == Magma_tally2_ELL ) {
             //printf("using ELL kernel for SpMV: ");
             CHECK( magma_tally2_dgeelltmv_shift( Magma_tally2NoTrans, A.num_rows, A.num_cols,
                A.max_nnz_row, alpha, lambda, A.dval, A.dcol, x.dval, beta, offset,
                blocksize, add_rows, y.dval, queue ));
             //printf("done.\n");
         }
         else {
             printf("error: format not supported.\n");
             info = MAGMA_tally2_ERR_NOT_SUPPORTED;
         }
    }
    // CPU case missing!
    else {
        printf("error: CPU not yet supported.\n");
        info = MAGMA_tally2_ERR_NOT_SUPPORTED;
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
    alpha       double
                scalar alpha

    @param[in]
    A           magma_tally2_d_matrix
                sparse matrix A
                
    @param[in]
    B           magma_tally2_d_matrix
                sparse matrix C
                
    @param[out]
    C           magma_tally2_d_matrix *
                outpur sparse matrix C

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_d
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_d_spmm(
    double alpha,
    magma_tally2_d_matrix A,
    magma_tally2_d_matrix B,
    magma_tally2_d_matrix *C,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    // set queue for old dense routines
    //magma_tally2_queue_t orig_queue=NULL;
    //magma_tally2blasGetKernelStream( &orig_queue );

    if ( A.memory_location != B.memory_location ) {
        printf("error: linear algebra objects are not located in same memory!\n");
        printf("memory locations are: %d   %d\n",
                        A.memory_location, B.memory_location );
        info = MAGMA_tally2_ERR_INVALID_PTR;
        goto cleanup;
    }

    // DEV case
    if ( A.memory_location == Magma_tally2_DEV ) {
        if ( A.num_cols == B.num_rows ) {

             if ( A.storage_type == Magma_tally2_CSR
                            || A.storage_type == Magma_tally2_CSRL
                            || A.storage_type == Magma_tally2_CSRU
                            || A.storage_type == Magma_tally2_CSRCOO ) {
                    
                CHECK( magma_tally2_dcuspmm( A, B, C, queue ));
                
             }
             else {
                 printf("error: format not supported.\n");
                 // magma_tally2blasSetKernelStream( orig_queue );
                 info = MAGMA_tally2_ERR_NOT_SUPPORTED;
             }
        }
         
    }
    // CPU case missing!
    else {
        printf("error: CPU not yet supported.\n");
        // magma_tally2blasSetKernelStream( orig_queue );
        info = MAGMA_tally2_ERR_NOT_SUPPORTED; // TODO change to goto cleanup?
    }
    
cleanup:
    // magma_tally2blasSetKernelStream( orig_queue );
    return info;
}





