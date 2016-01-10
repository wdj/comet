/*
    -- MAGMA_tally4 (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @generated from testing_zspmv.cpp normal z -> s, Tue May  5 14:04:40 2015
       @author Hartwig Anzt
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef MAGMA_tally4_WITH_MKL
    #include <mkl_spblas.h>

    #define PRECISION_s
    #if defined(PRECISION_z)
    #define MKL_ADDR(a) (float*)(a)
    #elif defined(PRECISION_c)
    #define MKL_ADDR(a) (MKL_Complex8*)(a)
    #else
    #define MKL_ADDR(a) (a)
    #endif
#endif

// includes, project
#include "flops.h"
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"
#include "testings.h"
#include "common_magma_tally4sparse.h"


/* ////////////////////////////////////////////////////////////////////////////
   -- testing sparse matrix vector product
*/
int main(  int argc, char** argv )
{
    magma_tally4_int_t info = 0;
    TESTING_INIT();
    magma_tally4_queue_t queue=NULL;
    magma_tally4_queue_create( /*devices[ opts->device ],*/ &queue );

    magma_tally4_s_matrix hA={Magma_tally4_CSR}, hA_SELLP={Magma_tally4_CSR}, hA_ELL={Magma_tally4_CSR}, 
    dA={Magma_tally4_CSR}, dA_SELLP={Magma_tally4_CSR}, dA_ELL={Magma_tally4_CSR};
    
    magma_tally4_s_matrix hx={Magma_tally4_CSR}, hy={Magma_tally4_CSR}, dx={Magma_tally4_CSR}, 
    dy={Magma_tally4_CSR}, hrefvec={Magma_tally4_CSR}, hcheck={Magma_tally4_CSR};
            
    hA_SELLP.blocksize = 8;
    hA_SELLP.alignment = 8;
    real_Double_t start, end, res;

    float c_one  = MAGMA_tally4_S_MAKE(1.0, 0.0);
    float c_zero = MAGMA_tally4_S_MAKE(0.0, 0.0);
    
    cusparseMatDescr_t descrA=NULL;
    cusparseHandle_t cusparseHandle = NULL;
    cusparseHybMat_t hybA=NULL;
    cusparseMatDescr_t descr = NULL;
    
    #ifdef MAGMA_tally4_WITH_MKL
        magma_tally4_int_t *pntre=NULL;
    #endif
    
    magma_tally4_int_t i, j;
    for( i = 1; i < argc; ++i ) {
        if ( strcmp("--blocksize", argv[i]) == 0 ) {
            hA_SELLP.blocksize = atoi( argv[++i] );
        } else if ( strcmp("--alignment", argv[i]) == 0 ) {
            hA_SELLP.alignment = atoi( argv[++i] );
        } else
            break;
    }
    printf( "\n#    usage: ./run_sspmv"
        " [ --blocksize %d --alignment %d (for SELLP) ]"
        " matrices \n\n", (int) hA_SELLP.blocksize, (int) hA_SELLP.alignment );

    while( i < argc ) {
        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_tally4_int_t laplace_size = atoi( argv[i] );
            CHECK( magma_tally4_sm_5stencil(  laplace_size, &hA, queue ));
        } else {                        // file-matrix test
            CHECK( magma_tally4_s_csr_mtx( &hA,  argv[i], queue ));
        }

        printf( "\n# matrix info: %d-by-%d with %d nonzeros\n\n",
                            (int) hA.num_rows,(int) hA.num_cols,(int) hA.nnz );

        real_Double_t FLOPS = 2.0*hA.nnz/1e9;

        // init CPU vectors
        CHECK( magma_tally4_svinit( &hx, Magma_tally4_CPU, hA.num_rows, 1, c_zero, queue ));
        CHECK( magma_tally4_svinit( &hy, Magma_tally4_CPU, hA.num_rows, 1, c_zero, queue ));

        // init DEV vectors
        CHECK( magma_tally4_svinit( &dx, Magma_tally4_DEV, hA.num_rows, 1, c_one, queue ));
        CHECK( magma_tally4_svinit( &dy, Magma_tally4_DEV, hA.num_rows, 1, c_zero, queue ));
       

        #ifdef MAGMA_tally4_WITH_MKL
            // calling MKL with CSR
            CHECK( magma_tally4_imalloc_cpu( &pntre, hA.num_rows + 1 ) );
            pntre[0] = 0;
            for (j=0; j<hA.num_rows; j++ ) {
                pntre[j] = hA.row[j+1];
            }
             MKL_INT num_rows = hA.num_rows;
             MKL_INT num_cols = hA.num_cols;
             MKL_INT nnz = hA.nnz;

            MKL_INT *col;
            TESTING_MALLOC_CPU( col, MKL_INT, nnz );
            for( magma_tally4_int_t t=0; t < hA.nnz; ++t ) {
                col[ t ] = hA.col[ t ];
            }
            MKL_INT *row;
            TESTING_MALLOC_CPU( row, MKL_INT, num_rows );
            for( magma_tally4_int_t t=0; t < hA.num_rows; ++t ) {
                row[ t ] = hA.col[ t ];
            }
    
            start = magma_tally4_wtime();
            for (j=0; j<10; j++ ) {
                mkl_scsrmv( "N", &num_rows, &num_cols,
                            MKL_ADDR(&c_one), "GFNC", MKL_ADDR(hA.val),
                            col, row, pntre,
                                                    MKL_ADDR(hx.val),
                            MKL_ADDR(&c_zero),        MKL_ADDR(hy.val) );
            }
            end = magma_tally4_wtime();
            printf( "\n > MKL  : %.2e seconds %.2e GFLOP/s    (CSR).\n",
                                            (end-start)/10, FLOPS*10/(end-start) );

            TESTING_FREE_CPU( row );
            TESTING_FREE_CPU( col );
            row = NULL;
            col = NULL;
        #endif // MAGMA_tally4_WITH_MKL

        // copy matrix to GPU
        CHECK( magma_tally4_smtransfer( hA, &dA, Magma_tally4_CPU, Magma_tally4_DEV, queue ));
        // SpMV on GPU (CSR) -- this is the reference!
        start = magma_tally4_sync_wtime( queue );
        for (j=0; j<10; j++)
            CHECK( magma_tally4_s_spmv( c_one, dA, dx, c_zero, dy, queue ));
        end = magma_tally4_sync_wtime( queue );
        printf( " > MAGMA_tally4: %.2e seconds %.2e GFLOP/s    (standard CSR).\n",
                                        (end-start)/10, FLOPS*10/(end-start) );
        magma_tally4_smfree(&dA, queue );
        CHECK( magma_tally4_smtransfer( dy, &hrefvec , Magma_tally4_DEV, Magma_tally4_CPU, queue ));

        // convert to ELL and copy to GPU
        CHECK( magma_tally4_smconvert(  hA, &hA_ELL, Magma_tally4_CSR, Magma_tally4_ELL, queue ));
        CHECK( magma_tally4_smtransfer( hA_ELL, &dA_ELL, Magma_tally4_CPU, Magma_tally4_DEV, queue ));
        magma_tally4_smfree(&hA_ELL, queue );
        magma_tally4_smfree( &dy, queue );
        CHECK( magma_tally4_svinit( &dy, Magma_tally4_DEV, hA.num_rows, 1, c_zero, queue ));
        // SpMV on GPU (ELL)
        start = magma_tally4_sync_wtime( queue );
        for (j=0; j<10; j++)
            CHECK( magma_tally4_s_spmv( c_one, dA_ELL, dx, c_zero, dy, queue ));
        end = magma_tally4_sync_wtime( queue );
        printf( " > MAGMA_tally4: %.2e seconds %.2e GFLOP/s    (standard ELL).\n",
                                        (end-start)/10, FLOPS*10/(end-start) );
        magma_tally4_smfree(&dA_ELL, queue );
        CHECK( magma_tally4_smtransfer( dy, &hcheck , Magma_tally4_DEV, Magma_tally4_CPU, queue ));
        res = 0.0;
        for(magma_tally4_int_t k=0; k<hA.num_rows; k++ )
            res=res + MAGMA_tally4_S_REAL(hcheck.val[k]) - MAGMA_tally4_S_REAL(hrefvec.val[k]);
        if ( res < .000001 )
            printf("# tester spmv ELL:  ok\n");
        else
            printf("# tester spmv ELL:  failed\n");
        magma_tally4_smfree( &hcheck, queue );

        // convert to SELLP and copy to GPU
        CHECK( magma_tally4_smconvert(  hA, &hA_SELLP, Magma_tally4_CSR, Magma_tally4_SELLP, queue ));
        CHECK( magma_tally4_smtransfer( hA_SELLP, &dA_SELLP, Magma_tally4_CPU, Magma_tally4_DEV, queue ));
        magma_tally4_smfree(&hA_SELLP, queue );
        magma_tally4_smfree( &dy, queue );
        CHECK( magma_tally4_svinit( &dy, Magma_tally4_DEV, hA.num_rows, 1, c_zero, queue ));
        // SpMV on GPU (SELLP)
        start = magma_tally4_sync_wtime( queue );
        for (j=0; j<10; j++)
            CHECK( magma_tally4_s_spmv( c_one, dA_SELLP, dx, c_zero, dy, queue ));
        end = magma_tally4_sync_wtime( queue );
        printf( " > MAGMA_tally4: %.2e seconds %.2e GFLOP/s    (SELLP).\n",
                                        (end-start)/10, FLOPS*10/(end-start) );

        CHECK( magma_tally4_smtransfer( dy, &hcheck , Magma_tally4_DEV, Magma_tally4_CPU, queue ));
        res = 0.0;
        for(magma_tally4_int_t k=0; k<hA.num_rows; k++ )
            res=res + MAGMA_tally4_S_REAL(hcheck.val[k]) - MAGMA_tally4_S_REAL(hrefvec.val[k]);
        printf("# |x-y|_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# tester spmv SELL-P:  ok\n");
        else
            printf("# tester spmv SELL-P:  failed\n");
        magma_tally4_smfree( &hcheck, queue );

        magma_tally4_smfree(&dA_SELLP, queue );


        // SpMV on GPU (CUSPARSE - CSR)
        // CUSPARSE context //

        CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
        CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
        CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));

        CHECK_CUSPARSE( cusparseSetMatType( descr, CUSPARSE_MATRIX_TYPE_GENERAL ));
        CHECK_CUSPARSE( cusparseSetMatIndexBase( descr, CUSPARSE_INDEX_BASE_ZERO ));
        float alpha = c_one;
        float beta = c_zero;
        magma_tally4_smfree( &dy, queue );
        CHECK( magma_tally4_svinit( &dy, Magma_tally4_DEV, hA.num_rows, 1, c_zero, queue ));

        // copy matrix to GPU
        CHECK( magma_tally4_smtransfer( hA, &dA, Magma_tally4_CPU, Magma_tally4_DEV, queue ));

        start = magma_tally4_sync_wtime( queue );
        for (j=0; j<10; j++)
            CHECK_CUSPARSE(
            cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                        hA.num_rows, hA.num_cols, hA.nnz, &alpha, descr,
                        dA.dval, dA.drow, dA.dcol, dx.dval, &beta, dy.dval) );
        end = magma_tally4_sync_wtime( queue );
        printf( " > CUSPARSE: %.2e seconds %.2e GFLOP/s    (CSR).\n",
                                        (end-start)/10, FLOPS*10/(end-start) );
        CHECK_CUSPARSE( cusparseCreateMatDescr( &descrA ));
        cusparseCreateHybMat( &hybA );
        CHECK( magma_tally4_smtransfer( dy, &hcheck , Magma_tally4_DEV, Magma_tally4_CPU, queue ));
        res = 0.0;
        for(magma_tally4_int_t k=0; k<hA.num_rows; k++ )
            res=res + MAGMA_tally4_S_REAL(hcheck.val[k]) - MAGMA_tally4_S_REAL(hrefvec.val[k]);
        printf("# |x-y|_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# tester spmv cuSPARSE CSR:  ok\n");
        else
            printf("# tester spmv cuSPARSE CSR:  failed\n");
        magma_tally4_smfree( &hcheck, queue );
        magma_tally4_smfree( &dy, queue );
        CHECK( magma_tally4_svinit( &dy, Magma_tally4_DEV, hA.num_rows, 1, c_zero, queue ));
       
        cusparseScsr2hyb(cusparseHandle,  hA.num_rows, hA.num_cols,
                        descrA, dA.dval, dA.drow, dA.dcol,
                        hybA, 0, CUSPARSE_HYB_PARTITION_AUTO);

        start = magma_tally4_sync_wtime( queue );
        for (j=0; j<10; j++)
            CHECK_CUSPARSE(
            cusparseShybmv( cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
               &alpha, descrA, hybA,
               dx.dval, &beta, dy.dval) );
        end = magma_tally4_sync_wtime( queue );
        printf( " > CUSPARSE: %.2e seconds %.2e GFLOP/s    (HYB).\n",
                                        (end-start)/10, FLOPS*10/(end-start) );

        CHECK( magma_tally4_smtransfer( dy, &hcheck , Magma_tally4_DEV, Magma_tally4_CPU, queue ));
        res = 0.0;
        for(magma_tally4_int_t k=0; k<hA.num_rows; k++ )
            res=res + MAGMA_tally4_S_REAL(hcheck.val[k]) - MAGMA_tally4_S_REAL(hrefvec.val[k]);
        printf("# |x-y|_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# tester spmv cuSPARSE HYB:  ok\n");
        else
            printf("# tester spmv cuSPARSE HYB:  failed\n");
        magma_tally4_smfree( &hcheck, queue );

        cusparseDestroyMatDescr( descrA );
        cusparseDestroyHybMat( hybA );
        cusparseDestroy( cusparseHandle ); 
        descrA=NULL;
        cusparseHandle = NULL;
        hybA=NULL;
        descr = NULL;

        // free CPU memory
        magma_tally4_smfree(&hA, queue );
        magma_tally4_smfree(&hx, queue );
        magma_tally4_smfree(&hy, queue );
        magma_tally4_smfree(&hrefvec, queue );
        // free GPU memory
        magma_tally4_smfree(&dA, queue );
        magma_tally4_smfree(&dx, queue );
        magma_tally4_smfree(&dy, queue );
        
        printf("\n\n");

        i++;
    }
    
cleanup:
    #ifdef MAGMA_tally4_WITH_MKL
        magma_tally4_free_cpu(pntre);
    #endif    
    cusparseDestroyMatDescr( descrA );
    cusparseDestroyHybMat( hybA );
    cusparseDestroy( cusparseHandle ); 
    magma_tally4_smfree(&hA, queue );
    magma_tally4_smfree(&dA, queue );
    magma_tally4_smfree(&hA_ELL, queue );
    magma_tally4_smfree(&dA_ELL, queue );
    magma_tally4_smfree(&hA_SELLP, queue );
    magma_tally4_smfree(&dA_SELLP, queue );
    magma_tally4_smfree( &hcheck, queue );
    magma_tally4_smfree(&hx, queue );
    magma_tally4_smfree(&hy, queue );
    magma_tally4_smfree(&hrefvec, queue );
    magma_tally4_smfree(&dx, queue );
    magma_tally4_smfree(&dy, queue );
    
    magma_tally4_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}