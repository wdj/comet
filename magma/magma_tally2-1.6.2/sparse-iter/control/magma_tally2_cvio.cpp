/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally2_zvio.cpp normal z -> c, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/
#include "common_magma_tally2sparse.h"


/**
    Purpose
    -------

    Visualizes part of a vector of type magma_tally2_c_matrix.
    With input vector x , offset, visulen, the entries
    offset - (offset +  visulen) of x are visualized.

    Arguments
    ---------

    @param[in]
    x           magma_tally2_c_matrix
                vector to visualize

    @param[in]
    offset      magma_tally2_int_t
                start inex of visualization

    @param[in]
    visulen     magma_tally2_int_t
                number of entries to visualize

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_caux
    ********************************************************************/

extern "C"
magma_tally2_int_t
magma_tally2_cprint_vector(
    magma_tally2_c_matrix x,
    magma_tally2_int_t offset,
    magma_tally2_int_t  visulen,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    magma_tally2_c_matrix y={Magma_tally2_CSR};
    
    //**************************************************************
    #define COMPLEX
    magma_tally2FloatComplex c_zero = MAGMA_tally2_C_ZERO;
    
    #ifdef COMPLEX
    #define magma_tally2_cprintval( tmp )       {                                  \
        if ( MAGMA_tally2_C_EQUAL( tmp, c_zero )) {                                \
            printf( "   0.              \n" );                                \
        }                                                                   \
        else {                                                              \
            printf( " %8.4f+%8.4fi\n",                                        \
                    MAGMA_tally2_C_REAL( tmp ), MAGMA_tally2_C_IMAG( tmp ));              \
        }                                                                   \
    }
    #else
    #define magma_tally2_cprintval( tmp )       {                                  \
        if ( MAGMA_tally2_C_EQUAL( tmp, c_zero )) {                                \
            printf( "   0.    \n" );                                          \
        }                                                                   \
        else {                                                              \
            printf( " %8.4f\n", MAGMA_tally2_C_REAL( tmp ));                         \
        }                                                                   \
    }
    #endif
    //**************************************************************
    
    printf("visualize entries %d - %d of vector ",
                    (int) offset, (int) (offset + visulen) );
    fflush(stdout);
    if ( x.memory_location == Magma_tally2_CPU ) {
        printf("located on CPU:\n");
        for( magma_tally2_int_t i=offset; i<offset + visulen; i++ )
            magma_tally2_cprintval(x.val[i]);
    }
    else if ( x.memory_location == Magma_tally2_DEV ) {
        printf("located on DEV:\n");
        CHECK( magma_tally2_cmtransfer( x, &y, Magma_tally2_DEV, Magma_tally2_CPU, queue ));
        for( magma_tally2_int_t i=offset; i<offset +  visulen; i++ )
            magma_tally2_cprintval(y.val[i]);


    }

cleanup:
    magma_tally2_free_cpu(y.val);
    return info;
}





/**
    Purpose
    -------

    Reads in a float vector of length "length".

    Arguments
    ---------

    @param[out]
    x           magma_tally2_c_matrix *
                vector to read in

    @param[in]
    length      magma_tally2_int_t
                length of vector
    @param[in]
    filename    char*
                file where vector is stored
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_caux
    ********************************************************************/

extern "C"
magma_tally2_int_t
magma_tally2_cvread(
    magma_tally2_c_matrix *x,
    magma_tally2_int_t length,
    char * filename,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    magma_tally2_int_t nnz=0, i=0;
    FILE *fid;
    
    x->memory_location = Magma_tally2_CPU;
    x->storage_type = Magma_tally2_DENSE;
    x->num_rows = length;
    x->num_cols = 1;
    x->major = Magma_tally2ColMajor;
    CHECK( magma_tally2_cmalloc_cpu( &x->val, length ));
    
    fid = fopen(filename, "r");
    
    while( i<length )  // eof() is 'true' at the end of data
    {
        float VAL1;

        magma_tally2FloatComplex VAL;
        #define COMPLEX
        
        #ifdef COMPLEX
            float VAL2;
            fscanf(fid, " %f %f \n", &VAL1, &VAL2);
            VAL = MAGMA_tally2_C_MAKE(VAL1, VAL2);
        #else
            fscanf(fid, " %f \n", &VAL1);
            VAL = MAGMA_tally2_C_MAKE(VAL1, 0.0);
        #endif
        
        if ( VAL != MAGMA_tally2_C_ZERO )
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
    x           magma_tally2_c_matrix *
                vector to read in

    @param[in]
    filename    char*
                file where vector is stored
    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_caux
    ********************************************************************/

extern "C"
magma_tally2_int_t
magma_tally2_cvspread(
    magma_tally2_c_matrix *x,
    const char * filename,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    magma_tally2_c_matrix A={Magma_tally2_CSR}, B={Magma_tally2_CSR};
    magma_tally2_int_t entry=0;
     //   char *vfilename[] = {"/mnt/sparse_matrices/mtx/rail_79841_B.mtx"};
    CHECK( magma_tally2_c_csr_mtx( &A,  filename, queue  ));
    CHECK( magma_tally2_cmconvert( A, &B, Magma_tally2_CSR, Magma_tally2_DENSE, queue ));
    CHECK( magma_tally2_cvinit( x, Magma_tally2_CPU, A.num_cols, A.num_rows, MAGMA_tally2_C_ZERO, queue ));
    x->major = Magma_tally2RowMajor;
    for(magma_tally2_int_t i=0; i<A.num_cols; i++) {
        for(magma_tally2_int_t j=0; j<A.num_rows; j++) {
            x->val[i*A.num_rows+j] = B.val[ i+j*A.num_cols ];
            entry++;
        }
    }
    x->num_rows = A.num_rows;
    x->num_cols = A.num_cols;
    
cleanup:
    magma_tally2_cmfree( &A, queue );
    magma_tally2_cmfree( &B, queue );
    return info;
}


