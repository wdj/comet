/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_minproduct_z_blaswrapper.cpp normal z -> c, Sun May  3 11:22:57 2015
       @author Hartwig Anzt

*/
#include "common_magma_minproductsparse.h"


#include "common_magma_minproductsparse.h"
#include "magma_minproductblas.h"
#include "magma_minproductsparse_types.h"




/**
    Purpose
    -------

    For a given input matrix A and vectors x, y and scalars alpha, beta
    the wrapper determines the suitable SpMV computing
              y = alpha * A * x + beta * y.
    Arguments
    ---------

    @param[in]
    alpha       magma_minproductFloatComplex
                scalar alpha

    @param[in]
    A           magma_minproduct_c_matrix
                sparse matrix A

    @param[in]
    x           magma_minproduct_c_matrix
                input vector x
                
    @param[in]
    beta        magma_minproductFloatComplex
                scalar beta
    @param[out]
    y           magma_minproduct_c_matrix
                output vector y
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_c
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_c_spmv(
    magma_minproductFloatComplex alpha,
    magma_minproduct_c_matrix A,
    magma_minproduct_c_matrix x,
    magma_minproductFloatComplex beta,
    magma_minproduct_c_matrix y,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );

    magma_minproduct_c_matrix x2={Magma_minproduct_CSR};

    cusparseHandle_t cusparseHandle = 0;
    cusparseMatDescr_t descr = 0;
    // make sure RHS is a dense matrix
    if ( x.storage_type != Magma_minproduct_DENSE ) {
         printf("error: only dense vectors are supported.\n");
         info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
         goto cleanup;
    }

    if ( A.memory_location != x.memory_location ||
                            x.memory_location != y.memory_location ) {
        printf("error: linear algebra objects are not located in same memory!\n");
        printf("memory locations are: %d   %d   %d\n",
                        A.memory_location, x.memory_location, y.memory_location );
        info = MAGMA_minproduct_ERR_INVALID_PTR;
        goto cleanup;
    }

    // DEV case
    if ( A.memory_location == Magma_minproduct_DEV ) {
        if ( A.num_cols == x.num_rows && x.num_cols == 1 ) {

             if ( A.storage_type == Magma_minproduct_CSR
                            || A.storage_type == Magma_minproduct_CSRL
                            || A.storage_type == Magma_minproduct_CSRU ) {
                 //printf("using CSR kernel for SpMV: ");
                 //magma_minproduct_cgecsrmv( Magma_minproductNoTrans, A.num_rows, A.num_cols, alpha,
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
             else if ( A.storage_type == Magma_minproduct_ELL ) {
                 //printf("using ELLPACKT kernel for SpMV: ");
                 CHECK( magma_minproduct_cgeelltmv( Magma_minproductNoTrans, A.num_rows, A.num_cols,
                    A.max_nnz_row, alpha, A.dval, A.dcol, x.dval, beta,
                    y.dval, queue ));
                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_minproduct_ELLPACKT ) {
                 //printf("using ELL kernel for SpMV: ");
                 CHECK( magma_minproduct_cgeellmv( Magma_minproductNoTrans, A.num_rows, A.num_cols,
                    A.max_nnz_row, alpha, A.dval, A.dcol, x.dval, beta,
                    y.dval, queue ));
                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_minproduct_ELLRT ) {
                 //printf("using ELLRT kernel for SpMV: ");
                 CHECK( magma_minproduct_cgeellrtmv( Magma_minproductNoTrans, A.num_rows, A.num_cols,
                            A.max_nnz_row, alpha, A.dval, A.dcol, A.drow, x.dval,
                         beta, y.dval, A.alignment, A.blocksize, queue ));
                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_minproduct_SELLP ) {
                 //printf("using SELLP kernel for SpMV: ");
                 CHECK( magma_minproduct_cgesellpmv( Magma_minproductNoTrans, A.num_rows, A.num_cols,
                    A.blocksize, A.numblocks, A.alignment,
                    alpha, A.dval, A.dcol, A.drow, x.dval, beta, y.dval, queue ));

                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_minproduct_DENSE ) {
                 //printf("using DENSE kernel for SpMV: ");
                 magma_minproductblas_cgemv( Magma_minproductNoTrans, A.num_rows, A.num_cols, alpha,
                                A.dval, A.num_rows, x.dval, 1, beta,  y.dval,
                                1 );
                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_minproduct_SPMVFUNCTION ) {
                 //printf("using DENSE kernel for SpMV: ");
                 CHECK( magma_minproduct_ccustomspmv( alpha, x, beta, y, queue ));
                 //printf("done.\n");
             }
             else if ( A.storage_type == Magma_minproduct_BCSR ) {
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
                 info = MAGMA_minproduct_ERR_NOT_SUPPORTED; 
             }
        }
        else if ( A.num_cols < x.num_rows || x.num_cols > 1 ) {
            magma_minproduct_int_t num_vecs = x.num_rows / A.num_cols * x.num_cols;
            if ( A.storage_type == Magma_minproduct_CSR ) {

                CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
                CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
                CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));
                CHECK_CUSPARSE( cusparseSetMatType( descr, CUSPARSE_MATRIX_TYPE_GENERAL ));
                CHECK_CUSPARSE( cusparseSetMatIndexBase( descr, CUSPARSE_INDEX_BASE_ZERO ));

                if ( x.major == Magma_minproductColMajor) {
                    cusparseCcsrmm(cusparseHandle,
                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                    A.num_rows,   num_vecs, A.num_cols, A.nnz,
                    &alpha, descr, A.dval, A.drow, A.dcol,
                    x.dval, A.num_cols, &beta, y.dval, A.num_cols);
                } else if ( x.major == Magma_minproductRowMajor) {
                    cusparseCcsrmm2(cusparseHandle,
                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                    CUSPARSE_OPERATION_TRANSPOSE,
                    A.num_rows,   num_vecs, A.num_cols, A.nnz,
                    &alpha, descr, A.dval, A.drow, A.dcol,
                    x.dval, A.num_cols, &beta, y.dval, A.num_cols);
                }

             } else if ( A.storage_type == Magma_minproduct_SELLP ) {
                if ( x.major == Magma_minproductRowMajor) {
                 CHECK( magma_minproduct_cmgesellpmv( Magma_minproductNoTrans, A.num_rows, A.num_cols,
                    num_vecs, A.blocksize, A.numblocks, A.alignment,
                    alpha, A.dval, A.dcol, A.drow, x.dval, beta, y.dval, queue ));
                }
                else if ( x.major == Magma_minproductColMajor) {
                    // transpose first to row major
                    CHECK( magma_minproduct_cvtranspose( x, &x2, queue ));
                    CHECK( magma_minproduct_cmgesellpmv( Magma_minproductNoTrans, A.num_rows, A.num_cols,
                    num_vecs, A.blocksize, A.numblocks, A.alignment,
                    alpha, A.dval, A.dcol, A.drow, x2.dval, beta, y.dval, queue ));
                }
             }
             /*if ( A.storage_type == Magma_minproduct_DENSE ) {
                 //printf("using DENSE kernel for SpMV: ");
                 magma_minproductblas_cmgemv( Magma_minproductNoTrans, A.num_rows, A.num_cols,
                            num_vecs, alpha, A.dval, A.num_rows, x.dval, 1,
                            beta,  y.dval, 1 );
                 //printf("done.\n");
             }*/
             else {
                 printf("error: format not supported.\n");
                 info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
             }
        }
         
         
    }
    // CPU case missing!
    else {
        printf("error: CPU not yet supported.\n");
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
    }

cleanup:
    cusparseDestroyMatDescr( descr );
    cusparseDestroy( cusparseHandle );
    cusparseHandle = 0;
    descr = 0;
    magma_minproduct_cmfree(&x2, queue );
    
    magma_minproductblasSetKernelStream( orig_queue );
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
    alpha       magma_minproductFloatComplex
                scalar alpha

    @param
    A           magma_minproduct_c_matrix
                sparse matrix A

    @param
    lambda      magma_minproductFloatComplex
                scalar lambda

    @param
    x           magma_minproduct_c_matrix
                input vector x

    @param
    beta        magma_minproductFloatComplex
                scalar beta
                
    @param
    offset      magma_minproduct_int_t
                in case not the main diagonal is scaled
                
    @param
    blocksize   magma_minproduct_int_t
                in case of processing multiple vectors
                
    @param
    add_rows    magma_minproduct_int_t*
                in case the matrixpowerskernel is used
                
    @param
    y           magma_minproduct_c_matrix
                output vector y
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_caux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_c_spmv_shift(
    magma_minproductFloatComplex alpha,
    magma_minproduct_c_matrix A,
    magma_minproductFloatComplex lambda,
    magma_minproduct_c_matrix x,
    magma_minproductFloatComplex beta,
    magma_minproduct_int_t offset,
    magma_minproduct_int_t blocksize,
    magma_minproduct_index_t *add_rows,
    magma_minproduct_c_matrix y,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;

    // make sure RHS is a dense matrix
    if ( x.storage_type != Magma_minproduct_DENSE ) {
         printf("error: only dense vectors are supported.\n");
         info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
         goto cleanup;
    }


    if ( A.memory_location != x.memory_location
                || x.memory_location != y.memory_location ) {
        printf("error: linear algebra objects are not located in same memory!\n");
        printf("memory locations are: %d   %d   %d\n",
                    A.memory_location, x.memory_location, y.memory_location );
        info = MAGMA_minproduct_ERR_INVALID_PTR;
        goto cleanup;
    }

    // DEV case
    if ( A.memory_location == Magma_minproduct_DEV ) {
         if ( A.storage_type == Magma_minproduct_CSR ) {
             //printf("using CSR kernel for SpMV: ");
             CHECK( magma_minproduct_cgecsrmv_shift( Magma_minproductNoTrans, A.num_rows, A.num_cols,
                alpha, lambda, A.dval, A.drow, A.dcol, x.dval, beta, offset,
                blocksize, add_rows, y.dval, queue ));
             //printf("done.\n");
         }
         else if ( A.storage_type == Magma_minproduct_ELLPACKT ) {
             //printf("using ELLPACKT kernel for SpMV: ");
             CHECK( magma_minproduct_cgeellmv_shift( Magma_minproductNoTrans, A.num_rows, A.num_cols,
                A.max_nnz_row, alpha, lambda, A.dval, A.dcol, x.dval, beta, offset,
                blocksize, add_rows, y.dval, queue ));
             //printf("done.\n");
         }
         else if ( A.storage_type == Magma_minproduct_ELL ) {
             //printf("using ELL kernel for SpMV: ");
             CHECK( magma_minproduct_cgeelltmv_shift( Magma_minproductNoTrans, A.num_rows, A.num_cols,
                A.max_nnz_row, alpha, lambda, A.dval, A.dcol, x.dval, beta, offset,
                blocksize, add_rows, y.dval, queue ));
             //printf("done.\n");
         }
         else {
             printf("error: format not supported.\n");
             info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
         }
    }
    // CPU case missing!
    else {
        printf("error: CPU not yet supported.\n");
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
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
    alpha       magma_minproductFloatComplex
                scalar alpha

    @param[in]
    A           magma_minproduct_c_matrix
                sparse matrix A
                
    @param[in]
    B           magma_minproduct_c_matrix
                sparse matrix C
                
    @param[out]
    C           magma_minproduct_c_matrix *
                outpur sparse matrix C

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_c
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_c_spmm(
    magma_minproductFloatComplex alpha,
    magma_minproduct_c_matrix A,
    magma_minproduct_c_matrix B,
    magma_minproduct_c_matrix *C,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // set queue for old dense routines
    //magma_minproduct_queue_t orig_queue=NULL;
    //magma_minproductblasGetKernelStream( &orig_queue );

    if ( A.memory_location != B.memory_location ) {
        printf("error: linear algebra objects are not located in same memory!\n");
        printf("memory locations are: %d   %d\n",
                        A.memory_location, B.memory_location );
        info = MAGMA_minproduct_ERR_INVALID_PTR;
        goto cleanup;
    }

    // DEV case
    if ( A.memory_location == Magma_minproduct_DEV ) {
        if ( A.num_cols == B.num_rows ) {

             if ( A.storage_type == Magma_minproduct_CSR
                            || A.storage_type == Magma_minproduct_CSRL
                            || A.storage_type == Magma_minproduct_CSRU
                            || A.storage_type == Magma_minproduct_CSRCOO ) {
                    
                CHECK( magma_minproduct_ccuspmm( A, B, C, queue ));
                
             }
             else {
                 printf("error: format not supported.\n");
                 // magma_minproductblasSetKernelStream( orig_queue );
                 info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
             }
        }
         
    }
    // CPU case missing!
    else {
        printf("error: CPU not yet supported.\n");
        // magma_minproductblasSetKernelStream( orig_queue );
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED; // TODO change to goto cleanup?
    }
    
cleanup:
    // magma_minproductblasSetKernelStream( orig_queue );
    return info;
}





