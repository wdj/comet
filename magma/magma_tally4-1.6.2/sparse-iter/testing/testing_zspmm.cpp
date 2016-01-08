/*
    -- MAGMA_tally4 (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @precisions normal z -> s d c
       @author Hartwig Anzt
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <cusparse_v2.h>
#include <cuda_profiler_api.h>

#ifdef MAGMA_tally4_WITH_MKL
    #include "mkl_spblas.h"
    
    #define PRECISION_z
    #if defined(PRECISION_z)
    #define MKL_ADDR(a) (MKL_Complex16*)(a)
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
    
    magma_tally4_z_matrix hA={Magma_tally4_CSR}, hA_SELLP={Magma_tally4_CSR}, hA_ELL={Magma_tally4_CSR}, 
    dA={Magma_tally4_CSR}, dA_SELLP={Magma_tally4_CSR}, dA_ELL={Magma_tally4_CSR};
    
    magma_tally4_z_matrix hx={Magma_tally4_CSR}, hy={Magma_tally4_CSR}, dx={Magma_tally4_CSR}, 
    dy={Magma_tally4_CSR}, hrefvec={Magma_tally4_CSR}, hcheck={Magma_tally4_CSR};
        
    hA_SELLP.blocksize = 8;
    hA_SELLP.alignment = 8;
    real_Double_t start, end, res;
    #ifdef MAGMA_tally4_WITH_MKL
        magma_tally4_int_t *pntre=NULL;
    #endif
    cusparseHandle_t cusparseHandle = NULL;
    cusparseMatDescr_t descr = NULL;

    magma_tally4DoubleComplex c_one  = MAGMA_tally4_Z_MAKE(1.0, 0.0);
    magma_tally4DoubleComplex c_zero = MAGMA_tally4_Z_MAKE(0.0, 0.0);
    
    magma_tally4_int_t i, j;
    for( i = 1; i < argc; ++i ) {
        if ( strcmp("--blocksize", argv[i]) == 0 ) {
            hA_SELLP.blocksize = atoi( argv[++i] );
        } else if ( strcmp("--alignment", argv[i]) == 0 ) {
            hA_SELLP.alignment = atoi( argv[++i] );
        } else
            break;
    }
    printf( "\n#    usage: ./run_zspmm"
        " [ --blocksize %d --alignment %d (for SELLP) ]"
        " matrices \n\n", (int) hA_SELLP.blocksize, (int) hA_SELLP.alignment );

    while( i < argc ) {
        if ( strcmp("LAPLACE2D", argv[i]) == 0 && i+1 < argc ) {   // Laplace test
            i++;
            magma_tally4_int_t laplace_size = atoi( argv[i] );
            CHECK( magma_tally4_zm_5stencil(  laplace_size, &hA, queue ));
        } else {                        // file-matrix test
            CHECK( magma_tally4_z_csr_mtx( &hA,  argv[i], queue ));
        }

        printf( "# matrix info: %d-by-%d with %d nonzeros\n",
                            (int) hA.num_rows,(int) hA.num_cols,(int) hA.nnz );

        real_Double_t FLOPS = 2.0*hA.nnz/1e9;



        // m - number of rows for the sparse matrix
        // n - number of vectors to be multiplied in the SpMM product
        magma_tally4_int_t m, n;

        m = hA.num_rows;
        n = 48;

        // init CPU vectors
        CHECK( magma_tally4_zvinit( &hx, Magma_tally4_CPU, m, n, c_one, queue ));
        CHECK( magma_tally4_zvinit( &hy, Magma_tally4_CPU, m, n, c_zero, queue ));

        // init DEV vectors
        CHECK( magma_tally4_zvinit( &dx, Magma_tally4_DEV, m, n, c_one, queue ));
        CHECK( magma_tally4_zvinit( &dy, Magma_tally4_DEV, m, n, c_zero, queue ));


        // calling MKL with CSR
        #ifdef MAGMA_tally4_WITH_MKL
            magma_tally4_int_t *pntre=NULL;
            CHECK( magma_tally4_imalloc_cpu( &pntre, m + 1 ) );
            pntre[0] = 0;
            for (j=0; j < m; j++ ) {
                pntre[j] = hA.row[j+1];
            }

            MKL_INT num_rows = hA.num_rows;
            MKL_INT num_cols = hA.num_cols;
            MKL_INT nnz = hA.nnz;
            MKL_INT num_vecs = n;

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

            // === Call MKL with consecutive SpMVs, using mkl_zcsrmv ===
            // warmp up
            mkl_zcsrmv( "N", &num_rows, &num_cols,
                        MKL_ADDR(&c_one), "GFNC", MKL_ADDR(hA.val), col, row, pntre,
                                                MKL_ADDR(hx.val),
                        MKL_ADDR(&c_zero),        MKL_ADDR(hy.val) );
    
            start = magma_tally4_wtime();
            for (j=0; j<10; j++ )
                mkl_zcsrmv( "N", &num_rows, &num_cols,
                        MKL_ADDR(&c_one), "GFNC", MKL_ADDR(hA.val), col, row, pntre,
                                                MKL_ADDR(hx.val),
                        MKL_ADDR(&c_zero),        MKL_ADDR(hy.val) );
            end = magma_tally4_wtime();
            printf( "\n > MKL SpMVs : %.2e seconds %.2e GFLOP/s    (CSR).\n",
                                            (end-start)/10, FLOPS*10/(end-start) );
    
            // === Call MKL with blocked SpMVs, using mkl_zcsrmm ===
            char transa = 'n';
            MKL_INT ldb = n, ldc=n;
            char matdescra[6] = {'g', 'l', 'n', 'c', 'x', 'x'};
    
            // warm up
            mkl_zcsrmm( &transa, &num_rows, &num_vecs, &num_cols, MKL_ADDR(&c_one), matdescra,
                      MKL_ADDR(hA.val), col, row, pntre,
                      MKL_ADDR(hx.val), &ldb,
                      MKL_ADDR(&c_zero),
                      MKL_ADDR(hy.val), &ldc );
    
            start = magma_tally4_wtime();
            for (j=0; j<10; j++ )
                mkl_zcsrmm( &transa, &num_rows, &num_vecs, &num_cols, MKL_ADDR(&c_one), matdescra,
                          MKL_ADDR(hA.val), col, row, pntre,
                          MKL_ADDR(hx.val), &ldb,
                          MKL_ADDR(&c_zero),
                          MKL_ADDR(hy.val), &ldc );
            end = magma_tally4_wtime();
            printf( "\n > MKL SpMM  : %.2e seconds %.2e GFLOP/s    (CSR).\n",
                    (end-start)/10, FLOPS*10.*n/(end-start) );

            TESTING_FREE_CPU( row );
            TESTING_FREE_CPU( col );
            row = NULL;
            col = NULL;

        #endif // MAGMA_tally4_WITH_MKL

        // copy matrix to GPU
        CHECK( magma_tally4_zmtransfer( hA, &dA, Magma_tally4_CPU, Magma_tally4_DEV, queue ));
        // SpMV on GPU (CSR)
        start = magma_tally4_sync_wtime( queue );
        for (j=0; j<10; j++)
            CHECK( magma_tally4_z_spmv( c_one, dA, dx, c_zero, dy, queue ));
        end = magma_tally4_sync_wtime( queue );
        printf( " > MAGMA_tally4: %.2e seconds %.2e GFLOP/s    (standard CSR).\n",
                                        (end-start)/10, FLOPS*10.*n/(end-start) );

        CHECK( magma_tally4_zmtransfer( dy, &hrefvec , Magma_tally4_DEV, Magma_tally4_CPU, queue ));
        magma_tally4_zmfree(&dA, queue );


        // convert to SELLP and copy to GPU
        CHECK( magma_tally4_zmconvert(  hA, &hA_SELLP, Magma_tally4_CSR, Magma_tally4_SELLP, queue ));
        CHECK( magma_tally4_zmtransfer( hA_SELLP, &dA_SELLP, Magma_tally4_CPU, Magma_tally4_DEV, queue ));
        magma_tally4_zmfree(&hA_SELLP, queue );
        magma_tally4_zmfree( &dy, queue );
        CHECK( magma_tally4_zvinit( &dy, Magma_tally4_DEV, dx.num_rows, dx.num_cols, c_zero, queue ));
        // SpMV on GPU (SELLP)
        start = magma_tally4_sync_wtime( queue );
        for (j=0; j<10; j++)
            CHECK( magma_tally4_z_spmv( c_one, dA_SELLP, dx, c_zero, dy, queue ));
        end = magma_tally4_sync_wtime( queue );
        printf( " > MAGMA_tally4: %.2e seconds %.2e GFLOP/s    (SELLP).\n",
                                        (end-start)/10, FLOPS*10.*n/(end-start) );

        CHECK( magma_tally4_zmtransfer( dy, &hcheck , Magma_tally4_DEV, Magma_tally4_CPU, queue ));
        res = 0.0;
        for(magma_tally4_int_t k=0; k<hA.num_rows; k++ )
            res=res + MAGMA_tally4_Z_REAL(hcheck.val[k]) - MAGMA_tally4_Z_REAL(hrefvec.val[k]);
        printf("# |x-y|_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# tester spmm SELL-P:  ok\n");
        else
            printf("# tester spmm SELL-P:  failed\n");
        magma_tally4_zmfree( &hcheck, queue );
        magma_tally4_zmfree(&dA_SELLP, queue );



        // SpMV on GPU (CUSPARSE - CSR)
        // CUSPARSE context //
        magma_tally4_zmfree( &dy, queue );
        CHECK( magma_tally4_zvinit( &dy, Magma_tally4_DEV, dx.num_rows, dx.num_cols, c_zero, queue ));
        //#ifdef PRECISION_d
        start = magma_tally4_sync_wtime( queue );
        CHECK_CUSPARSE( cusparseCreate( &cusparseHandle ));
        CHECK_CUSPARSE( cusparseSetStream( cusparseHandle, queue ));
        CHECK_CUSPARSE( cusparseCreateMatDescr( &descr ));
        CHECK_CUSPARSE( cusparseSetMatType( descr, CUSPARSE_MATRIX_TYPE_GENERAL ));
        CHECK_CUSPARSE( cusparseSetMatIndexBase( descr, CUSPARSE_INDEX_BASE_ZERO ));
        magma_tally4DoubleComplex alpha = c_one;
        magma_tally4DoubleComplex beta = c_zero;

        // copy matrix to GPU
        CHECK( magma_tally4_zmtransfer( hA, &dA, Magma_tally4_CPU, Magma_tally4_DEV, queue) );

        for (j=0; j<10; j++)
        cusparseZcsrmm(cusparseHandle,
            CUSPARSE_OPERATION_NON_TRANSPOSE,
                    dA.num_rows,   n, dA.num_cols, dA.nnz,
                    &alpha, descr, dA.dval, dA.drow, dA.dcol,
                    dx.dval, dA.num_cols, &beta, dy.dval, dA.num_cols);
        end = magma_tally4_sync_wtime( queue );
        printf( " > CUSPARSE: %.2e seconds %.2e GFLOP/s    (CSR).\n",
                                        (end-start)/10, FLOPS*10*n/(end-start) );

        CHECK( magma_tally4_zmtransfer( dy, &hcheck , Magma_tally4_DEV, Magma_tally4_CPU, queue ));
        res = 0.0;
        for(magma_tally4_int_t k=0; k<hA.num_rows; k++ )
            res=res + MAGMA_tally4_Z_REAL(hcheck.val[k]) - MAGMA_tally4_Z_REAL(hrefvec.val[k]);
        printf("# |x-y|_F = %8.2e\n", res);
        if ( res < .000001 )
            printf("# tester spmm cuSPARSE:  ok\n");
        else
            printf("# tester spmm cuSPARSE:  failed\n");
        magma_tally4_zmfree( &hcheck, queue );

        cusparseDestroyMatDescr( descr ); 
        cusparseDestroy( cusparseHandle );
        descr = NULL;
        cusparseHandle = NULL;
        //#endif

        printf("\n\n");


        // free CPU memory
        magma_tally4_zmfree(&hA, queue );
        magma_tally4_zmfree(&hx, queue );
        magma_tally4_zmfree(&hy, queue );
        magma_tally4_zmfree(&hrefvec, queue );
        // free GPU memory
        magma_tally4_zmfree(&dx, queue );
        magma_tally4_zmfree(&dy, queue );
        magma_tally4_zmfree(&dA, queue);

        i++;
    }

cleanup:
    #ifdef MAGMA_tally4_WITH_MKL
        magma_tally4_free_cpu(pntre);
    #endif
    cusparseDestroyMatDescr( descr ); 
    cusparseDestroy( cusparseHandle );
    magma_tally4_zmfree(&hA, queue );
    magma_tally4_zmfree(&dA, queue );
    magma_tally4_zmfree(&hA_ELL, queue );
    magma_tally4_zmfree(&dA_ELL, queue );
    magma_tally4_zmfree(&hA_SELLP, queue );
    magma_tally4_zmfree(&dA_SELLP, queue );
    
    magma_tally4_queue_destroy( queue );
    TESTING_FINALIZE();
    return info;
}
