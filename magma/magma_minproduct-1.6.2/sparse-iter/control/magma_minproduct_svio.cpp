/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_minproduct_zvio.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/
#include "common_magma_minproductsparse.h"


/**
    Purpose
    -------

    Visualizes part of a vector of type magma_minproduct_s_matrix.
    With input vector x , offset, visulen, the entries
    offset - (offset +  visulen) of x are visualized.

    Arguments
    ---------

    @param[in]
    x           magma_minproduct_s_matrix
                vector to visualize

    @param[in]
    offset      magma_minproduct_int_t
                start inex of visualization

    @param[in]
    visulen     magma_minproduct_int_t
                number of entries to visualize

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_saux
    ********************************************************************/

extern "C"
magma_minproduct_int_t
magma_minproduct_sprint_vector(
    magma_minproduct_s_matrix x,
    magma_minproduct_int_t offset,
    magma_minproduct_int_t  visulen,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_s_matrix y={Magma_minproduct_CSR};
    
    //**************************************************************
    #define REAL
    float c_zero = MAGMA_minproduct_S_ZERO;
    
    #ifdef COMPLEX
    #define magma_minproduct_sprintval( tmp )       {                                  \
        if ( MAGMA_minproduct_S_EQUAL( tmp, c_zero )) {                                \
            printf( "   0.              \n" );                                \
        }                                                                   \
        else {                                                              \
            printf( " %8.4f+%8.4fi\n",                                        \
                    MAGMA_minproduct_S_REAL( tmp ), MAGMA_minproduct_S_IMAG( tmp ));              \
        }                                                                   \
    }
    #else
    #define magma_minproduct_sprintval( tmp )       {                                  \
        if ( MAGMA_minproduct_S_EQUAL( tmp, c_zero )) {                                \
            printf( "   0.    \n" );                                          \
        }                                                                   \
        else {                                                              \
            printf( " %8.4f\n", MAGMA_minproduct_S_REAL( tmp ));                         \
        }                                                                   \
    }
    #endif
    //**************************************************************
    
    printf("visualize entries %d - %d of vector ",
                    (int) offset, (int) (offset + visulen) );
    fflush(stdout);
    if ( x.memory_location == Magma_minproduct_CPU ) {
        printf("located on CPU:\n");
        for( magma_minproduct_int_t i=offset; i<offset + visulen; i++ )
            magma_minproduct_sprintval(x.val[i]);
    }
    else if ( x.memory_location == Magma_minproduct_DEV ) {
        printf("located on DEV:\n");
        CHECK( magma_minproduct_smtransfer( x, &y, Magma_minproduct_DEV, Magma_minproduct_CPU, queue ));
        for( magma_minproduct_int_t i=offset; i<offset +  visulen; i++ )
            magma_minproduct_sprintval(y.val[i]);


    }

cleanup:
    magma_minproduct_free_cpu(y.val);
    return info;
}





/**
    Purpose
    -------

    Reads in a float vector of length "length".

    Arguments
    ---------

    @param[out]
    x           magma_minproduct_s_matrix *
                vector to read in

    @param[in]
    length      magma_minproduct_int_t
                length of vector
    @param[in]
    filename    char*
                file where vector is stored
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_saux
    ********************************************************************/

extern "C"
magma_minproduct_int_t
magma_minproduct_svread(
    magma_minproduct_s_matrix *x,
    magma_minproduct_int_t length,
    char * filename,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_int_t nnz=0, i=0;
    FILE *fid;
    
    x->memory_location = Magma_minproduct_CPU;
    x->storage_type = Magma_minproduct_DENSE;
    x->num_rows = length;
    x->num_cols = 1;
    x->major = Magma_minproductColMajor;
    CHECK( magma_minproduct_smalloc_cpu( &x->val, length ));
    
    fid = fopen(filename, "r");
    
    while( i<length )  // eof() is 'true' at the end of data
    {
        float VAL1;

        float VAL;
        #define REAL
        
        #ifdef COMPLEX
            float VAL2;
            fscanf(fid, " %f %f \n", &VAL1, &VAL2);
            VAL = MAGMA_minproduct_S_MAKE(VAL1, VAL2);
        #else
            fscanf(fid, " %f \n", &VAL1);
            VAL = MAGMA_minproduct_S_MAKE(VAL1, 0.0);
        #endif
        
        if ( VAL != MAGMA_minproduct_S_ZERO )
            nnz++;
        x->val[i] = VAL;
        i++;
    }
    fclose(fid);
    
    x->nnz = nnz;
    
cleanup:
    return info;
}




/**
    Purpose
    -------

    Reads in a sparse vector-block stored in COO format.

    Arguments
    ---------

    @param[out]
    x           magma_minproduct_s_matrix *
                vector to read in

    @param[in]
    filename    char*
                file where vector is stored
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_saux
    ********************************************************************/

extern "C"
magma_minproduct_int_t
magma_minproduct_svspread(
    magma_minproduct_s_matrix *x,
    const char * filename,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_s_matrix A={Magma_minproduct_CSR}, B={Magma_minproduct_CSR};
    magma_minproduct_int_t entry=0;
     //   char *vfilename[] = {"/mnt/sparse_matrices/mtx/rail_79841_B.mtx"};
    CHECK( magma_minproduct_s_csr_mtx( &A,  filename, queue  ));
    CHECK( magma_minproduct_smconvert( A, &B, Magma_minproduct_CSR, Magma_minproduct_DENSE, queue ));
    CHECK( magma_minproduct_svinit( x, Magma_minproduct_CPU, A.num_cols, A.num_rows, MAGMA_minproduct_S_ZERO, queue ));
    x->major = Magma_minproductRowMajor;
    for(magma_minproduct_int_t i=0; i<A.num_cols; i++) {
        for(magma_minproduct_int_t j=0; j<A.num_rows; j++) {
            x->val[i*A.num_rows+j] = B.val[ i+j*A.num_cols ];
            entry++;
        }
    }
    x->num_rows = A.num_rows;
    x->num_cols = A.num_cols;
    
cleanup:
    magma_minproduct_smfree( &A, queue );
    magma_minproduct_smfree( &B, queue );
    return info;
}


