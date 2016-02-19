/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt
*/
#include "common_magma_tally3sparse.h"


/**
    Purpose
    -------

    Visualizes part of a vector of type magma_tally3_z_matrix.
    With input vector x , offset, visulen, the entries
    offset - (offset +  visulen) of x are visualized.

    Arguments
    ---------

    @param[in]
    x           magma_tally3_z_matrix
                vector to visualize

    @param[in]
    offset      magma_tally3_int_t
                start inex of visualization

    @param[in]
    visulen     magma_tally3_int_t
                number of entries to visualize

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zaux
    ********************************************************************/

extern "C"
magma_tally3_int_t
magma_tally3_zprint_vector(
    magma_tally3_z_matrix x,
    magma_tally3_int_t offset,
    magma_tally3_int_t  visulen,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    magma_tally3_z_matrix y={Magma_tally3_CSR};
    
    //**************************************************************
    #define COMPLEX
    magma_tally3DoubleComplex c_zero = MAGMA_tally3_Z_ZERO;
    
    #ifdef COMPLEX
    #define magma_tally3_zprintval( tmp )       {                                  \
        if ( MAGMA_tally3_Z_EQUAL( tmp, c_zero )) {                                \
            printf( "   0.              \n" );                                \
        }                                                                   \
        else {                                                              \
            printf( " %8.4f+%8.4fi\n",                                        \
                    MAGMA_tally3_Z_REAL( tmp ), MAGMA_tally3_Z_IMAG( tmp ));              \
        }                                                                   \
    }
    #else
    #define magma_tally3_zprintval( tmp )       {                                  \
        if ( MAGMA_tally3_Z_EQUAL( tmp, c_zero )) {                                \
            printf( "   0.    \n" );                                          \
        }                                                                   \
        else {                                                              \
            printf( " %8.4f\n", MAGMA_tally3_Z_REAL( tmp ));                         \
        }                                                                   \
    }
    #endif
    //**************************************************************
    
    printf("visualize entries %d - %d of vector ",
                    (int) offset, (int) (offset + visulen) );
    fflush(stdout);
    if ( x.memory_location == Magma_tally3_CPU ) {
        printf("located on CPU:\n");
        for( magma_tally3_int_t i=offset; i<offset + visulen; i++ )
            magma_tally3_zprintval(x.val[i]);
    }
    else if ( x.memory_location == Magma_tally3_DEV ) {
        printf("located on DEV:\n");
        CHECK( magma_tally3_zmtransfer( x, &y, Magma_tally3_DEV, Magma_tally3_CPU, queue ));
        for( magma_tally3_int_t i=offset; i<offset +  visulen; i++ )
            magma_tally3_zprintval(y.val[i]);


    }

cleanup:
    magma_tally3_free_cpu(y.val);
    return info;
}





/**
    Purpose
    -------

    Reads in a double vector of length "length".

    Arguments
    ---------

    @param[out]
    x           magma_tally3_z_matrix *
                vector to read in

    @param[in]
    length      magma_tally3_int_t
                length of vector
    @param[in]
    filename    char*
                file where vector is stored
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zaux
    ********************************************************************/

extern "C"
magma_tally3_int_t
magma_tally3_zvread(
    magma_tally3_z_matrix *x,
    magma_tally3_int_t length,
    char * filename,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    magma_tally3_int_t nnz=0, i=0;
    FILE *fid;
    
    x->memory_location = Magma_tally3_CPU;
    x->storage_type = Magma_tally3_DENSE;
    x->num_rows = length;
    x->num_cols = 1;
    x->major = Magma_tally3ColMajor;
    CHECK( magma_tally3_zmalloc_cpu( &x->val, length ));
    
    fid = fopen(filename, "r");
    
    while( i<length )  // eof() is 'true' at the end of data
    {
        double VAL1;

        magma_tally3DoubleComplex VAL;
        #define COMPLEX
        
        #ifdef COMPLEX
            double VAL2;
            fscanf(fid, " %lf %lf \n", &VAL1, &VAL2);
            VAL = MAGMA_tally3_Z_MAKE(VAL1, VAL2);
        #else
            fscanf(fid, " %lf \n", &VAL1);
            VAL = MAGMA_tally3_Z_MAKE(VAL1, 0.0);
        #endif
        
        if ( VAL != MAGMA_tally3_Z_ZERO )
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
    x           magma_tally3_z_matrix *
                vector to read in

    @param[in]
    filename    char*
                file where vector is stored
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zaux
    ********************************************************************/

extern "C"
magma_tally3_int_t
magma_tally3_zvspread(
    magma_tally3_z_matrix *x,
    const char * filename,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    magma_tally3_z_matrix A={Magma_tally3_CSR}, B={Magma_tally3_CSR};
    magma_tally3_int_t entry=0;
     //   char *vfilename[] = {"/mnt/sparse_matrices/mtx/rail_79841_B.mtx"};
    CHECK( magma_tally3_z_csr_mtx( &A,  filename, queue  ));
    CHECK( magma_tally3_zmconvert( A, &B, Magma_tally3_CSR, Magma_tally3_DENSE, queue ));
    CHECK( magma_tally3_zvinit( x, Magma_tally3_CPU, A.num_cols, A.num_rows, MAGMA_tally3_Z_ZERO, queue ));
    x->major = Magma_tally3RowMajor;
    for(magma_tally3_int_t i=0; i<A.num_cols; i++) {
        for(magma_tally3_int_t j=0; j<A.num_rows; j++) {
            x->val[i*A.num_rows+j] = B.val[ i+j*A.num_cols ];
            entry++;
        }
    }
    x->num_rows = A.num_rows;
    x->num_cols = A.num_cols;
    
cleanup:
    magma_tally3_zmfree( &A, queue );
    magma_tally3_zmfree( &B, queue );
    return info;
}


