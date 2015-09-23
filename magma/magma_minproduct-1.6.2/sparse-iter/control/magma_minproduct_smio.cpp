/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_minproduct_zmio.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from
//  the IO functions provided by MatrixMarket

#include "common_magma_minproductsparse.h"
#include "mmio.h"



/**
    Purpose
    -------

    Reads in a matrix stored in coo format from a Matrix Market (.mtx)
    file and converts it into CSR format. It duplicates the off-diagonal
    entries in the symmetric case.

    Arguments
    ---------
    
    @param[out]
    type        magma_minproduct_storage_t*
                storage type of matrix
                
    @param[out]
    location    magma_minproduct_location_t*
                location of matrix
                
    @param[out]
    n_row       magma_minproduct_int_t*
                number of rows in matrix
                
    @param[out]
    n_col       magma_minproduct_int_t*
                number of columns in matrix
                
    @param[out]
    nnz         magma_minproduct_int_t*
                number of nonzeros in matrix
                
    @param[out]
    val         float**
                value array of CSR output

    @param[out]
    row         magma_minproduct_index_t**
                row pointer of CSR output

    @param[out]
    col         magma_minproduct_index_t**
                column indices of CSR output

    @param[in]
    filename    const char*
                filname of the mtx matrix
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_saux
    ********************************************************************/

extern "C"
magma_minproduct_int_t read_s_csr_from_mtx(
    magma_minproduct_storage_t *type,
    magma_minproduct_location_t *location,
    magma_minproduct_int_t* n_row,
    magma_minproduct_int_t* n_col,
    magma_minproduct_int_t* nnz,
    float **val,
    magma_minproduct_index_t **row,
    magma_minproduct_index_t **col,
    const char *filename,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    FILE *fid;
    MM_typecode matcode;
    
    magma_minproduct_index_t *coo_col=NULL, *coo_row=NULL;
    float *coo_val=NULL;
    magma_minproduct_index_t *new_col=NULL, *new_row=NULL;
    float *new_val=NULL;
    
    fid = fopen(filename, "r");
    
    if (fid == NULL) {
        printf("# Unable to open file %s\n", filename);
        info = MAGMA_minproduct_ERR_NOT_FOUND;
        goto cleanup;
    }
    
    if (mm_read_banner(fid, &matcode) != 0) {
        printf("#Could not process lMatrix Market banner.\n");
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
        goto cleanup;
    }
    
    if (!mm_is_valid(matcode)) {
        printf("#Invalid lMatrix Market file.\n");
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
        goto cleanup;
    }
    
    if ( ! ( (mm_is_real(matcode) || mm_is_integer(matcode)
           || mm_is_pattern(matcode) || mm_is_real(matcode) )
             && mm_is_coordinate(matcode)
             && mm_is_sparse(matcode) ) )
    {
        printf("#Sorry, this application does not support ");
        printf("#Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        printf("#Only real-valued or pattern coordinate matrices are supported\n");
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
        goto cleanup;
    }
    
    magma_minproduct_index_t num_rows, num_cols, num_nonzeros;
    if (mm_read_mtx_crd_size(fid, &num_rows, &num_cols, &num_nonzeros) != 0){
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
        goto cleanup;
    }
    
    *type     = Magma_minproduct_CSR;
    *location = Magma_minproduct_CPU;
    *n_row    = num_rows;
    *n_col    = num_cols;
    *nnz      = num_nonzeros;


    CHECK( magma_minproduct_index_malloc_cpu( &coo_col, *nnz ) );
    CHECK( magma_minproduct_index_malloc_cpu( &coo_row, *nnz ) );
    CHECK( magma_minproduct_smalloc_cpu( &coo_val, *nnz ) );

    if (mm_is_real(matcode) || mm_is_integer(matcode)) {
        for(magma_minproduct_int_t i = 0; i < *nnz; ++i) {
            magma_minproduct_index_t ROW, COL;
            float VAL;  // always read in a float and convert later if necessary
            
            fscanf(fid, " %d %d %f \n", &ROW, &COL, &VAL);
            
            coo_row[i] = ROW - 1;
            coo_col[i] = COL - 1;
            coo_val[i] = MAGMA_minproduct_S_MAKE( VAL, 0.);
        }
    } else if (mm_is_pattern(matcode) ) {
        for(magma_minproduct_int_t i = 0; i < *nnz; ++i) {
            magma_minproduct_index_t ROW, COL;
            
            fscanf(fid, " %d %d \n", &ROW, &COL );
            
            coo_row[i] = ROW - 1;
            coo_col[i] = COL - 1;
            coo_val[i] = MAGMA_minproduct_S_MAKE( 1.0, 0.);
        }
    } else if (mm_is_real(matcode) ){
       for(magma_minproduct_int_t i = 0; i < *nnz; ++i) {
            magma_minproduct_index_t ROW, COL;
            float VAL, VALC;  // always read in a float and convert later if necessary
            
            fscanf(fid, " %d %d %f %f\n", &ROW, &COL, &VAL, &VALC);
            
            coo_row[i] = ROW - 1;
            coo_col[i] = COL - 1;
            coo_val[i] = MAGMA_minproduct_S_MAKE( VAL, VALC);
        }
        // printf(" ...successfully read real matrix... ");
    } else {
        printf("Unrecognized data type\n");
        info = MAGMA_minproduct_ERR_NOT_FOUND;
    }
    
    fclose(fid);
    printf(" done\n");
    
    if (mm_is_symmetric(matcode)) { // duplicate off diagonal entries
        printf("detected symmetric case\n");
        magma_minproduct_index_t off_diagonals = 0;
        for(magma_minproduct_int_t i = 0; i < *nnz; ++i) {
            if (coo_row[i] != coo_col[i])
                ++off_diagonals;
        }
        
        magma_minproduct_index_t true_nonzeros = 2*off_diagonals + (*nnz - off_diagonals);
        
        printf("total number of nonzeros: %d\n", (int) *nnz);

        
        CHECK( magma_minproduct_index_malloc_cpu( &new_col, true_nonzeros ) );
        CHECK( magma_minproduct_index_malloc_cpu( &new_row, true_nonzeros ) );
        CHECK( magma_minproduct_smalloc_cpu( &new_val, true_nonzeros ) );
    
        magma_minproduct_index_t ptr = 0;
        for(magma_minproduct_int_t i = 0; i < *nnz; ++i) {
            if (coo_row[i] != coo_col[i]) {
                new_row[ptr] = coo_row[i];
                new_col[ptr] = coo_col[i];
                new_val[ptr] = coo_val[i];
                ptr++;
                new_col[ptr] = coo_row[i];
                new_row[ptr] = coo_col[i];
                new_val[ptr] = coo_val[i];
                ptr++;
            } else {
                new_row[ptr] = coo_row[i];
                new_col[ptr] = coo_col[i];
                new_val[ptr] = coo_val[i];
                ptr++;
            }
        }
        
        magma_minproduct_free_cpu(coo_row); 
        magma_minproduct_free_cpu(coo_col); 
        magma_minproduct_free_cpu(coo_val);

        coo_row = new_row;
        coo_col = new_col;
        coo_val = new_val;
        
        *nnz = true_nonzeros;
    } // end symmetric case
    
    float tv;
    magma_minproduct_index_t ti;
    
    // If matrix is not in standard format, sorting is necessary
    /*
    printf( "Sorting the cols....\n" );
    // bubble sort (by cols)
    for (int i=0; i < *nnz-1; ++i) {
        for (int j=0; j < *nnz-i-1; ++j) {
            if (coo_col[j] > coo_col[j+1] ) {
                ti = coo_col[j];
                coo_col[j] = coo_col[j+1];
                coo_col[j+1] = ti;

                ti = coo_row[j];
                coo_row[j] = coo_row[j+1];
                coo_row[j+1] = ti;

                tv = coo_val[j];
                coo_val[j] = coo_val[j+1];
                coo_val[j+1] = tv;
            }
        }
    }

    printf( "Sorting the rows....\n" );
    // bubble sort (by rows)
    for (int i=0; i < *nnz-1; ++i) {
        for (int j=0; j < *nnz-i-1; ++j) {
            if ( coo_row[j] > coo_row[j+1] ) {
                ti = coo_col[j];
                coo_col[j] = coo_col[j+1];
                coo_col[j+1] = ti;

                ti = coo_row[j];
                coo_row[j] = coo_row[j+1];
                coo_row[j+1] = ti;

                tv = coo_val[j];
                coo_val[j] = coo_val[j+1];
                coo_val[j+1] = tv;
            }
        }
    }
    printf( "Sorting: done\n" );
    
    */
    CHECK( magma_minproduct_smalloc_cpu( val, *nnz ) );
    
    CHECK( magma_minproduct_index_malloc_cpu( col, *nnz ) );
    CHECK( magma_minproduct_index_malloc_cpu( row, (*n_row+1) ) );
    CHECK( magma_minproduct_smalloc_cpu( val, *nnz ) );

    // original code from  Nathan Bell and Michael Garland
    // the output CSR structure is NOT sorted!

    for (magma_minproduct_index_t i = 0; i < num_rows; i++)
        (*row)[i] = 0;
    
    for (magma_minproduct_index_t i = 0; i < *nnz; i++)
        (*row)[coo_row[i]]++;
    
    // cumsum the nnz per row to get Bp[]
    for(magma_minproduct_int_t i = 0, cumsum = 0; i < num_rows; i++) {
        magma_minproduct_index_t temp = (*row)[i];
        (*row)[i] = cumsum;
        cumsum += temp;
    }
    (*row)[num_rows] = *nnz;
    
    // write Aj,Ax into Bj,Bx
    for(magma_minproduct_int_t i = 0; i < *nnz; i++) {
        magma_minproduct_index_t row_  = coo_row[i];
        magma_minproduct_index_t dest = (*row)[row_];
        
        (*col)[dest] = coo_col[i];
        
        (*val)[dest] = coo_val[i];
        
        (*row)[row_]++;
    }
    
    for(int i = 0, last = 0; i <= num_rows; i++) {
        int temp  = (*row)[i];
        (*row)[i] = last;
        last      = temp;
    }
    
    (*row)[*n_row] = *nnz;

    for (magma_minproduct_index_t k=0; k < *n_row; ++k) {
        for (magma_minproduct_index_t i=(*row)[k]; i < (*row)[k+1]-1; ++i) {
            for (magma_minproduct_index_t j=(*row)[k]; j < (*row)[k+1]-1; ++j) {
                if ( (*col)[j] > (*col)[j+1] ) {
                    ti = (*col)[j];
                    (*col)[j] = (*col)[j+1];
                    (*col)[j+1] = ti;
    
                    tv = (*val)[j];
                    (*val)[j] = (*val)[j+1];
                    (*val)[j+1] = tv;
                }
            }
        }
    }

cleanup:
    magma_minproduct_free_cpu(coo_row);
    magma_minproduct_free_cpu(coo_col);
    magma_minproduct_free_cpu(coo_val);
    return info;
}


extern "C" magma_minproduct_int_t
magma_minproduct_swrite_csrtomtx(
    magma_minproduct_s_matrix B,
    const char *filename,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    // TODO: why does this hard code Magma_minproductColMajor?
    CHECK( magma_minproduct_swrite_csr_mtx( B, Magma_minproductColMajor, filename, queue ));
cleanup:
    return info;
}


/**
    Purpose
    -------

    Writes a CSR matrix to a file using Matrix Market format.

    Arguments
    ---------

    @param[in]
    n_row       magma_minproduct_int_t
                number of rows in matrix
                
    @param[in]
    n_col       magma_minproduct_int_t
                number of columns in matrix
                
    @param[in]
    nnz         magma_minproduct_int_t
                number of nonzeros in matrix
                
    @param[in]
    val         float**
                value array of CSR
                TODO: why are these ** pointers? Wouldn't * pointers work?

    @param[in]
    row         magma_minproduct_index_t**
                row pointer of CSR

    @param[in]
    col         magma_minproduct_index_t**
                column indices of CSR

    @param[in]
    MajorType   magma_minproduct_index_t
                Row or Column sort
                default: 0 = RowMajor, 1 = ColMajor
                TODO: use named constants (e.g., Magma_minproductRowMajor), not numbers.

    @param[in]
    filename    const char*
                output-filname of the mtx matrix
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_saux
    ********************************************************************/

extern "C"
magma_minproduct_int_t
magma_minproduct_swrite_csr_mtx(
    magma_minproduct_s_matrix A,
    magma_minproduct_order_t MajorType,
    const char *filename,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    FILE *fp;
    magma_minproduct_s_matrix B = {Magma_minproduct_CSR};
    
    if ( MajorType == Magma_minproductColMajor ) {
        // to obtain ColMajor output we transpose the matrix
        // and flip the row and col pointer in the output
        
        CHECK( magma_minproduct_s_cucsrtranspose( A, &B, queue ));
        magma_minproduct_sprint_matrix( A, queue );
        
        
        magma_minproduct_sprint_matrix( B, queue );
        
        // TODO avoid duplicating this code below.
        printf("# Writing sparse matrix to file (%s):", filename);
        fflush(stdout);
        
        fp = fopen (filename, "w");
        if (  fp == NULL ){
            printf("error writing matrix: file exists or missing write permission\n");
            info = -1;
            goto cleanup;
        }
            
        #define REAL
        printf("check1\n");
        #ifdef COMPLEX
        // real case
        fprintf( fp, "%%%%MatrixMarket matrix coordinate real general\n" );
        fprintf( fp, "%d %d %d\n",B.num_cols, B.num_rows, B.nnz);
        
        // TODO what's the difference between i (or i+1) and rowindex?
        magma_minproduct_index_t i=0, j=0, rowindex=1;
        
        for(i=0; i < B.num_cols; i++) {
            magma_minproduct_index_t rowtemp1 = B.row[i];
            magma_minproduct_index_t rowtemp2 = B.row[i+1];
            for(j=0; j < rowtemp2 - rowtemp1; j++) {
                fprintf( fp, "%d %d %.16g %.16g\n",
                    ((B.col)[rowtemp1+j]+1), rowindex,
                    MAGMA_minproduct_S_REAL((B.val)[rowtemp1+j]),
                    MAGMA_minproduct_S_IMAG((B.val)[rowtemp1+j]) );

            }
            rowindex++;
        }
        #else
        // real case
        fprintf( fp, "%%%%MatrixMarket matrix coordinate real general\n" );
        fprintf( fp, "%d %d %d\n",B.num_cols, B.num_rows, B.nnz);
        
        // TODO what's the difference between i (or i+1) and rowindex?
        magma_minproduct_index_t i=0, j=0, rowindex=1;
                
        for(i=0; i < B.num_cols; i++) {
            magma_minproduct_index_t rowtemp1 = B.row[i];
            magma_minproduct_index_t rowtemp2 = B.row[i+1];
            for(j=0; j < rowtemp2 - rowtemp1; j++) {
                fprintf( fp, "%d %d %.16g\n",
                    ((B.col)[rowtemp1+j]+1), rowindex,
                    MAGMA_minproduct_S_REAL((B.val)[rowtemp1+j]) );

            }
            rowindex++;
        }
        #endif
       
        if(fclose (fp)!=0)
            printf("error: writing matrix failed\n");
        else
            printf(" done\n");
    }
    else {
        printf("# Writing sparse matrix to file (%s):", filename);
        fflush(stdout);
        
        fp = fopen (filename, "w");
        if (  fp == NULL ){
            printf("error writing matrix: file exists or missing write permission\n");
            info = -1;
            goto cleanup;
        }
             
            
        #define REAL
        printf("check1\n");
        #ifdef COMPLEX
        // real case
        fprintf( fp, "%%%%MatrixMarket matrix coordinate real general\n" );
        fprintf( fp, "%d %d %d\n",A.num_cols, A.num_rows, A.nnz);
        
        // TODO what's the difference between i (or i+1) and rowindex?
        magma_minproduct_index_t i=0, j=0, rowindex=1;
        
        for(i=0; i < A.num_cols; i++) {
            magma_minproduct_index_t rowtemp1 = A.row[i];
            magma_minproduct_index_t rowtemp2 = A.row[i+1];
            for(j=0; j < rowtemp2 - rowtemp1; j++) {
                fprintf( fp, "%d %d %.16g %.16g\n",
                    ((A.col)[rowtemp1+j]+1), rowindex,
                    MAGMA_minproduct_S_REAL((A.val)[rowtemp1+j]),
                    MAGMA_minproduct_S_IMAG((A.val)[rowtemp1+j]) );

            }
            rowindex++;
        }
        #else
        // real case
        fprintf( fp, "%%%%MatrixMarket matrix coordinate real general\n" );
        fprintf( fp, "%d %d %d\n",A.num_cols, A.num_rows, A.nnz);
        
        // TODO what's the difference between i (or i+1) and rowindex?
        magma_minproduct_index_t i=0, j=0, rowindex=1;
                
        for(i=0; i < B.num_cols; i++) {
            magma_minproduct_index_t rowtemp1 = A.row[i];
            magma_minproduct_index_t rowtemp2 = A.row[i+1];
            for(j=0; j < rowtemp2 - rowtemp1; j++) {
                fprintf( fp, "%d %d %.16g\n",
                    ((A.col)[rowtemp1+j]+1), rowindex,
                    MAGMA_minproduct_S_REAL((A.val)[rowtemp1+j]));

            }
            rowindex++;
        }
        #endif

        if(fclose (fp)!=0)
            printf("error: writing matrix failed\n");
        else
            printf(" done\n");
    }
cleanup:
    return info;
}


/**
    Purpose
    -------

    Prints a CSR matrix in Matrix Market format.

    Arguments
    ---------

    @param[in]
    n_row       magma_minproduct_int_t*
                number of rows in matrix
                
    @param[in]
    n_col       magma_minproduct_int_t*
                number of columns in matrix
                
    @param[in]
    nnz         magma_minproduct_int_t*
                number of nonzeros in matrix
                
    @param[in]
    val         float**
                value array of CSR

    @param[in]
    row         magma_minproduct_index_t**
                row pointer of CSR

    @param[in]
    col         magma_minproduct_index_t**
                column indices of CSR

    @param[in]
    MajorType   magma_minproduct_index_t
                Row or Column sort
                default: 0 = RowMajor, 1 = ColMajor
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_saux
    ********************************************************************/

extern "C"
magma_minproduct_int_t
magma_minproduct_sprint_csr_mtx(
    magma_minproduct_int_t n_row,
    magma_minproduct_int_t n_col,
    magma_minproduct_int_t nnz,
    float **val,
    magma_minproduct_index_t **row,
    magma_minproduct_index_t **col,
    magma_minproduct_order_t MajorType,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
        
    if ( MajorType == Magma_minproductColMajor ) {
        // to obtain ColMajor output we transpose the matrix
        // and flip the row and col pointer in the output
        float *new_val=NULL;
        magma_minproduct_index_t *new_row;
        magma_minproduct_index_t *new_col;
        magma_minproduct_int_t new_n_row;
        magma_minproduct_int_t new_n_col;
        magma_minproduct_int_t new_nnz;
        
        CHECK( s_transpose_csr( n_row, n_col, nnz, *val, *row, *col,
            &new_n_row, &new_n_col, &new_nnz, &new_val, &new_row, &new_col, queue) );
       
 
            
        #define REAL
        
        #ifdef COMPLEX
        // real case
        printf( "%%%%MatrixMarket matrix coordinate real general\n" );
        printf( "%d %d %d\n",new_n_col, new_n_row, new_nnz);
        
        // TODO what's the difference between i (or i+1) and rowindex?
        magma_minproduct_index_t i=0, j=0, rowindex=1;
        
        for(i=0; i < n_col; i++) {
            magma_minproduct_index_t rowtemp1 = (new_row)[i];
            magma_minproduct_index_t rowtemp2 = (new_row)[i+1];
            for(j=0; j < rowtemp2 - rowtemp1; j++) {
                printf( "%d %d %.6e %.6e\n",
                    ((new_col)[rowtemp1+j]+1), rowindex,
                    MAGMA_minproduct_S_REAL((new_val)[rowtemp1+j]),
                    MAGMA_minproduct_S_IMAG((new_val)[rowtemp1+j]) );

            }
            rowindex++;
        }
        
        #else
        // real case
        printf( "%%%%MatrixMarket matrix coordinate real general\n" );
        printf( "%d %d %d\n",new_n_col, new_n_row, new_nnz);
        
        // TODO what's the difference between i (or i+1) and rowindex?
        magma_minproduct_index_t i=0, j=0, rowindex=1;
        
        for(i=0; i < n_col; i++) {
            magma_minproduct_index_t rowtemp1 = (new_row)[i];
            magma_minproduct_index_t rowtemp2 = (new_row)[i+1];
            for(j=0; j < rowtemp2 - rowtemp1; j++) {
                printf( "%d %d %.6e\n",
                    ((new_col)[rowtemp1+j]+1), rowindex,
                    MAGMA_minproduct_S_REAL((new_val)[rowtemp1+j]) );
            }
            rowindex++;
        }
        #endif
       
        
    }
    else {

            
        #define REAL
        
        #ifdef COMPLEX
        // real case
        printf( "%%%%MatrixMarket matrix coordinate real general\n" );
        printf( "%d %d %d\n",n_col, n_row, nnz);
        
        // TODO what's the difference between i (or i+1) and rowindex?
        magma_minproduct_index_t i=0, j=0, rowindex=1;
        
        for(i=0; i < n_col; i++) {
            magma_minproduct_index_t rowtemp1 = (*row)[i];
            magma_minproduct_index_t rowtemp2 = (*row)[i+1];
            for(j=0; j < rowtemp2 - rowtemp1; j++) {
                printf( "%d %d %.6e %.6e\n",
                    rowindex, ((*col)[rowtemp1+j]+1),
                    MAGMA_minproduct_S_REAL((*val)[rowtemp1+j]),
                    MAGMA_minproduct_S_IMAG((*val)[rowtemp1+j]) );

            }
            rowindex++;
        }
        
        #else
        // real case
        printf( "%%%%MatrixMarket matrix coordinate real general\n" );
        printf( "%d %d %d\n",n_col, n_row, nnz);
        
        // TODO what's the difference between i (or i+1) and rowindex?
        magma_minproduct_index_t i=0, j=0, rowindex=1;
        
        for(i=0; i < n_col; i++) {
            magma_minproduct_index_t rowtemp1 = (*row)[i];
            magma_minproduct_index_t rowtemp2 = (*row)[i+1];
            for(j=0; j < rowtemp2 - rowtemp1; j++) {
                printf( "%d %d %.6e\n",
                    rowindex, ((*col)[rowtemp1+j]+1),
                    MAGMA_minproduct_S_REAL((*val)[rowtemp1+j]) );
            }
            rowindex++;
        }
        #endif
        
    }

cleanup:
    return info;
}


/**
    Purpose
    -------

    Prints a CSR matrix in CSR format.

    Arguments
    ---------
    
    @param[in]
    n_row       magma_minproduct_int_t*
                number of rows in matrix
                
    @param[in]
    n_col       magma_minproduct_int_t*
                number of columns in matrix
                
    @param[in]
    nnz         magma_minproduct_int_t*
                number of nonzeros in matrix
                
    @param[in]
    val         float**
                value array of CSR

    @param[in]
    row         magma_minproduct_index_t**
                row pointer of CSR

    @param[in]
    col         magma_minproduct_index_t**
                column indices of CSR

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_saux
    ********************************************************************/

extern "C"
magma_minproduct_int_t
magma_minproduct_sprint_csr(
    magma_minproduct_int_t n_row,
    magma_minproduct_int_t n_col,
    magma_minproduct_int_t nnz,
    float **val,
    magma_minproduct_index_t **row,
    magma_minproduct_index_t **col,
    magma_minproduct_queue_t queue )
{
    printf( "Matrix in CSR format (row col val)\n" );
    printf( " %d %d %d\n", n_row, n_col, nnz );
     
    magma_minproduct_index_t info = 0, i=0, j=0;

    for(i=0; i < n_col; i++) {
        magma_minproduct_index_t rowtemp1 = (*row)[i];
        magma_minproduct_index_t rowtemp2 = (*row)[i+1];
        for(j=0; j < rowtemp2 - rowtemp1; j++) {
                printf(" %d %d %.2f\n", (rowtemp1+1), (*col)[rowtemp1+j]+1,
                    MAGMA_minproduct_S_REAL((*val)[rowtemp1+j]) );
        }
    }
    
    return info;
}


/**
    Purpose
    -------

    Prints a sparse matrix in CSR format.

    Arguments
    ---------

    @param[in]
    A           magma_minproduct_s_matrix
                sparse matrix in Magma_minproduct_CSR format
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_saux
    ********************************************************************/

extern "C"
magma_minproduct_int_t
magma_minproduct_sprint_matrix(
    magma_minproduct_s_matrix A,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    //**************************************************************
    #define REAL
    
    #ifdef COMPLEX
    #define magma_minproduct_sprintval( tmp )       {                                  \
        if ( MAGMA_minproduct_S_EQUAL( tmp, c_zero )) {                                \
            printf( "   0.              " );                                \
        }                                                                   \
        else {                                                              \
            printf( " %8.4f+%8.4fi",                                        \
                    MAGMA_minproduct_S_REAL( tmp ), MAGMA_minproduct_S_IMAG( tmp ));              \
        }                                                                   \
    }
    #else
    #define magma_minproduct_sprintval( tmp )       {                                  \
        if ( MAGMA_minproduct_S_EQUAL( tmp, c_zero )) {                                \
            printf( "   0.    " );                                          \
        }                                                                   \
        else {                                                              \
            printf( " %8.4f", MAGMA_minproduct_S_REAL( tmp ));                         \
        }                                                                   \
    }
    #endif
    //**************************************************************
    
    magma_minproduct_index_t i, j, k;
    float c_zero = MAGMA_minproduct_S_ZERO;
    magma_minproduct_s_matrix C={Magma_minproduct_CSR};

    if ( A.memory_location == Magma_minproduct_CPU ) {
        printf("visualizing matrix of size %d %d with %d nonzeros:\n",
            (int) A.num_rows, (int) A.num_cols, (int) A.nnz);
        
        if ( A.storage_type == Magma_minproduct_DENSE ) {
            for( i=0; i < (A.num_rows); i++ ) {
                for( j=0; j < A.num_cols; j++ ) {
                    magma_minproduct_sprintval( A.val[i*(A.num_cols)+j] );
                }
                printf( "\n" );
            }
        }
        else if ( A.storage_type == Magma_minproduct_CSR ) {
            // visualize only small matrices like dense
            if ( A.num_rows < 11 && A.num_cols < 11 ) {
                CHECK( magma_minproduct_smconvert( A, &C, A.storage_type, Magma_minproduct_DENSE, queue ));
                CHECK( magma_minproduct_sprint_matrix(  C, queue ));
                magma_minproduct_smfree( &C, queue );
            }
            // otherwise visualize only coners
            else {
                // 4 beginning and 4 last elements of first four rows
                for( i=0; i < 4; i++ ) {
                    // upper left corner
                    for( j=0; j < 4; j++ ) {
                        float tmp = MAGMA_minproduct_S_ZERO;
                        magma_minproduct_index_t rbound = min( A.row[i]+4, A.row[i+1]);
                        magma_minproduct_index_t lbound = max( A.row[i], A.row[i]);
                        for( k=lbound; k < rbound; k++ ) {
                            if ( A.col[k] == j ) {
                                tmp = A.val[k];
                            }
                        }
                        magma_minproduct_sprintval( tmp );
                    }
                    if ( i == 0 ) {
                        printf( "    . . .    " );
                    } else {
                        printf( "             " );
                    }
                    // upper right corner
                    for( j=A.num_rows-4; j < A.num_rows; j++ ) {
                        float tmp = MAGMA_minproduct_S_ZERO;
                        magma_minproduct_index_t rbound = min( A.row[i+1], A.row[i+1]);
                        magma_minproduct_index_t lbound = max( A.row[i+1]-4, A.row[i]);
                        for( k=lbound; k < rbound; k++ ) {
                            if ( A.col[k] == j ) {
                                tmp = A.val[k];
                                                                
                            }
                        }
                        magma_minproduct_sprintval( tmp );
                    }
                    printf( "\n");
                }
                printf( "     .                     .         .         .\n"
                        "     .                         .         .         .\n"
                        "     .                             .         .         .\n"
                        "     .                                 .         .         .\n" );
                for( i=A.num_rows-4; i < A.num_rows; i++ ) {
                    // lower left corner
                    for( j=0; j < 4; j++ ) {
                        float tmp = MAGMA_minproduct_S_ZERO;
                        magma_minproduct_index_t rbound = min( A.row[i]+4, A.row[i+1]);
                        magma_minproduct_index_t lbound = max( A.row[i], A.row[i]);
                        for( k=lbound; k < rbound; k++ ) {
                            if ( A.col[k] == j ) {
                                tmp = A.val[k];
                            }
                        }
                        magma_minproduct_sprintval( tmp );
                    }
                    printf( "             ");
                    // lower right corner
                    for( j=A.num_rows-4; j < A.num_rows; j++ ) {
                        float tmp = MAGMA_minproduct_S_ZERO;
                        magma_minproduct_index_t rbound = min( A.row[i+1], A.row[i+1]);
                        magma_minproduct_index_t lbound = max( A.row[i+1]-4, A.row[i]);
                        for( k=lbound; k < rbound; k++ ) {
                            if ( A.col[k] == j ) {
                                tmp = A.val[k];
                            }
                        }
                        magma_minproduct_sprintval( tmp );
                    }
                    printf( "\n");
                }
            }
        }
        else {
            CHECK( magma_minproduct_smconvert( A, &C, A.storage_type, Magma_minproduct_CSR, queue ));
            CHECK( magma_minproduct_sprint_matrix(  C, queue ));
        }
    }
    else {
        magma_minproduct_s_matrix C={Magma_minproduct_CSR};
        CHECK( magma_minproduct_smtransfer( A, &C, A.memory_location, Magma_minproduct_CPU, queue ));
        CHECK( magma_minproduct_sprint_matrix(  C, queue ));
    }

cleanup:
    magma_minproduct_smfree( &C, queue );
    return info;
}


/**
    Purpose
    -------

    Reads in a matrix stored in coo format from a Matrix Market (.mtx)
    file and converts it into CSR format. It duplicates the off-diagonal
    entries in the symmetric case.

    Arguments
    ---------

    @param[out]
    A           magma_minproduct_s_matrix*
                matrix in magma_minproduct sparse matrix format

    @param[in]
    filename    const char*
                filname of the mtx matrix
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_saux
    ********************************************************************/

extern "C"
magma_minproduct_int_t
magma_minproduct_s_csr_mtx(
    magma_minproduct_s_matrix *A,
    const char *filename,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;

    int csr_compressor = 0;       // checks for zeros in original file
    
    magma_minproduct_s_matrix B={Magma_minproduct_CSR};

    magma_minproduct_index_t *coo_col = NULL;
    magma_minproduct_index_t *coo_row = NULL;
    float *coo_val = NULL;
    float *new_val = NULL;
    magma_minproduct_index_t* new_row = NULL;
    magma_minproduct_index_t* new_col = NULL;

    FILE *fid;
    MM_typecode matcode;
    fid = fopen(filename, "r");
    
    if (fid == NULL) {
        printf("#Unable to open file %s\n", filename);
        info = MAGMA_minproduct_ERR_NOT_FOUND;
        goto cleanup;
    }
    
    if (mm_read_banner(fid, &matcode) != 0) {
        printf("#Could not process lMatrix Market banner: %s.\n", matcode);
        info = MAGMA_minproduct_ERR_NOT_FOUND;
        goto cleanup;
    }
    
    if (!mm_is_valid(matcode)) {
        printf("#Invalid lMatrix Market file.\n");
        info = MAGMA_minproduct_ERR_NOT_FOUND;
        goto cleanup;
    }
    
    if ( ! ( (mm_is_real(matcode) || mm_is_integer(matcode)
           || mm_is_pattern(matcode) || mm_is_real(matcode) )
             && mm_is_coordinate(matcode)
             && mm_is_sparse(matcode) ) )
    {
        printf("#Sorry, this application does not support ");
        printf("#Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        printf("#Only real-valued or pattern coordinate matrices are supported\n");
        info = MAGMA_minproduct_ERR_NOT_FOUND;
        goto cleanup;
    }

    magma_minproduct_index_t num_rows, num_cols, num_nonzeros;
    if (mm_read_mtx_crd_size(fid, &num_rows, &num_cols, &num_nonzeros) != 0)
        info = MAGMA_minproduct_ERR_UNKNOWN;
    
    A->storage_type    = Magma_minproduct_CSR;
    A->memory_location = Magma_minproduct_CPU;
    A->num_rows        = num_rows;
    A->num_cols        = num_cols;
    A->nnz             = num_nonzeros;
    A->fill_mode       = Magma_minproduct_FULL;
    


    CHECK( magma_minproduct_index_malloc_cpu( &coo_col, A->nnz ) );
    CHECK( magma_minproduct_index_malloc_cpu( &coo_row, A->nnz ) );
    CHECK( magma_minproduct_smalloc_cpu( &coo_val, A->nnz ) );

    printf("# Reading sparse matrix from file (%s):", filename);
    fflush(stdout);
    if (mm_is_real(matcode) || mm_is_integer(matcode)) {
        for(magma_minproduct_int_t i = 0; i < A->nnz; ++i) {
            magma_minproduct_index_t ROW, COL;
            float VAL;  // always read in a float and convert later if necessary
            
            fscanf(fid, " %d %d %f \n", &ROW, &COL, &VAL);
            if ( VAL == 0 )
                csr_compressor = 1;
            coo_row[i] = ROW - 1;
            coo_col[i] = COL - 1;
            coo_val[i] = MAGMA_minproduct_S_MAKE( VAL, 0.);
        }
    } else if (mm_is_pattern(matcode) ) {
        for(magma_minproduct_int_t i = 0; i < A->nnz; ++i) {
            magma_minproduct_index_t ROW, COL;
            
            fscanf(fid, " %d %d \n", &ROW, &COL );
            
            coo_row[i] = ROW - 1;
            coo_col[i] = COL - 1;
            coo_val[i] = MAGMA_minproduct_S_MAKE( 1.0, 0.);
        }
    } else if (mm_is_real(matcode) ){
       for(magma_minproduct_int_t i = 0; i < A->nnz; ++i) {
            magma_minproduct_index_t ROW, COL;
            float VAL, VALC;  // always read in a float and convert later if necessary
            
            fscanf(fid, " %d %d %f %f\n", &ROW, &COL, &VAL, &VALC);
            
            coo_row[i] = ROW - 1;
            coo_col[i] = COL - 1;
            coo_val[i] = MAGMA_minproduct_S_MAKE( VAL, VALC);
        }
        // printf(" ...successfully read real matrix... ");
    } else {
        printf("Unrecognized data type\n");
        info = MAGMA_minproduct_ERR_NOT_FOUND;
    }
    fclose(fid);
    printf(" done.\n");
        
    A->sym = Magma_minproduct_GENERAL;

    if (mm_is_symmetric(matcode)) { // duplicate off diagonal entries
        A->sym = Magma_minproduct_SYMMETRIC;
        //printf("detected symmetric case\n");
        magma_minproduct_index_t off_diagonals = 0;
        for(magma_minproduct_int_t i = 0; i < A->nnz; ++i) {
            if (coo_row[i] != coo_col[i])
                ++off_diagonals;
        }
        magma_minproduct_index_t true_nonzeros = 2 * off_diagonals + (A->nnz - off_diagonals);
         


        CHECK( magma_minproduct_smalloc_cpu( &new_val, true_nonzeros ));
        CHECK( magma_minproduct_index_malloc_cpu( &new_row, true_nonzeros ));
        CHECK( magma_minproduct_index_malloc_cpu( &new_col, true_nonzeros ));
        
        magma_minproduct_index_t ptr = 0;
        for(magma_minproduct_int_t i = 0; i < A->nnz; ++i) {
            if (coo_row[i] != coo_col[i]) {
                new_row[ptr] = coo_row[i];
                new_col[ptr] = coo_col[i];
                new_val[ptr] = coo_val[i];
                ptr++;
                new_col[ptr] = coo_row[i];
                new_row[ptr] = coo_col[i];
                new_val[ptr] = coo_val[i];
                ptr++;
            } else {
                new_row[ptr] = coo_row[i];
                new_col[ptr] = coo_col[i];
                new_val[ptr] = coo_val[i];
                ptr++;
            }
        }
        
        magma_minproduct_free_cpu(coo_row);
        magma_minproduct_free_cpu(coo_col);
        magma_minproduct_free_cpu(coo_val);

        coo_row = new_row;
        coo_col = new_col;
        coo_val = new_val;
        A->nnz = true_nonzeros;
        //printf("total number of nonzeros: %d\n", A->nnz);
    } // end symmetric case
    
    float tv;
    magma_minproduct_index_t ti;
    
    // If matrix is not in standard format, sorting is necessary
    /*
    printf( "Sorting the cols....\n" );
    // bubble sort (by cols)
    for (int i=0; i < A->nnz-1; ++i) {
        for (int j=0; j < A->nnz-i-1; ++j) {
            if (coo_col[j] > coo_col[j+1] ) {
                ti = coo_col[j];
                coo_col[j] = coo_col[j+1];
                coo_col[j+1] = ti;
                
                ti = coo_row[j];
                coo_row[j] = coo_row[j+1];
                coo_row[j+1] = ti;
                
                tv = coo_val[j];
                coo_val[j] = coo_val[j+1];
                coo_val[j+1] = tv;
            }
        }
    }

    printf( "Sorting the rows....\n" );
    // bubble sort (by rows)
    for (int i=0; i < A->nnz-1; ++i) {
        for (int j=0; j < A->nnz-i-1; ++j) {
            if ( coo_row[j] > coo_row[j+1] ) {
                ti = coo_col[j];
                coo_col[j] = coo_col[j+1];
                coo_col[j+1] = ti;
                
                ti = coo_row[j];
                coo_row[j] = coo_row[j+1];
                coo_row[j+1] = ti;
                
                tv = coo_val[j];
                coo_val[j] = coo_val[j+1];
                coo_val[j+1] = tv;
            }
        }
    }
    printf( "Sorting: done\n" );
    
    */
    CHECK( magma_minproduct_smalloc_cpu( &A->val, A->nnz ));
    CHECK( magma_minproduct_index_malloc_cpu( &A->col, A->nnz ));
    CHECK( magma_minproduct_index_malloc_cpu( &A->row, A->num_rows+1 ));
    
    // original code from  Nathan Bell and Michael Garland
    // the output CSR structure is NOT sorted!

    for (magma_minproduct_index_t i = 0; i < num_rows; i++)
        (A->row)[i] = 0;
    
    for (magma_minproduct_index_t i = 0; i < A->nnz; i++)
        (A->row)[coo_row[i]]++;
        
    // cumsum the nnz per row to get Bp[]
    for(magma_minproduct_int_t i = 0, cumsum = 0; i < num_rows; i++) {
        magma_minproduct_index_t temp = (A->row)[i];
        (A->row)[i] = cumsum;
        cumsum += temp;
    }
    (A->row)[num_rows] = A->nnz;
    
    // write Aj,Ax into Bj,Bx
    for(magma_minproduct_int_t i = 0; i < A->nnz; i++) {
        magma_minproduct_index_t row_  = coo_row[i];
        magma_minproduct_index_t dest = (A->row)[row_];
        
        (A->col)[dest] = coo_col[i];
        
        (A->val)[dest] = coo_val[i];
        
        (A->row)[row_]++;
    }    
    magma_minproduct_free_cpu(coo_row);
    magma_minproduct_free_cpu(coo_col);
    magma_minproduct_free_cpu(coo_val);
    coo_row = NULL;
    coo_col = NULL;
    coo_val = NULL;

    for(int i = 0, last = 0; i <= num_rows; i++) {
        int temp    = (A->row)[i];
        (A->row)[i] = last;
        last        = temp;
    }
    
    (A->row)[A->num_rows]=A->nnz;
    
    for (magma_minproduct_index_t k=0; k < A->num_rows; ++k) {
        for (magma_minproduct_index_t i=(A->row)[k]; i < (A->row)[k+1]-1; ++i) {
            for (magma_minproduct_index_t j=(A->row)[k]; j < (A->row)[k+1]-1; ++j) {
                if ( (A->col)[j] > (A->col)[j+1] ) {
                    ti            = (A->col)[j];
                    (A->col)[j]   = (A->col)[j+1];
                    (A->col)[j+1] = ti;
                    
                    tv            = (A->val)[j];
                    (A->val)[j]   = (A->val)[j+1];
                    (A->val)[j+1] = tv;
                }
            }
        }
    }
    if ( csr_compressor > 0) { // run the CSR compressor to remove zeros
        //printf("removing zeros: ");
        CHECK( magma_minproduct_smtransfer( *A, &B, Magma_minproduct_CPU, Magma_minproduct_CPU, queue ));
        CHECK( magma_minproduct_s_csr_compressor(
            &(A->val), &(A->row), &(A->col),
            &B.val, &B.row, &B.col, &B.num_rows, queue ));
        B.nnz = B.row[num_rows];
        //printf(" remaining nonzeros:%d ", B.nnz);
        magma_minproduct_free_cpu( A->val );
        magma_minproduct_free_cpu( A->row );
        magma_minproduct_free_cpu( A->col );
        CHECK( magma_minproduct_smtransfer( B, A, Magma_minproduct_CPU, Magma_minproduct_CPU, queue ));
        //printf("done.\n");

    }
cleanup:
    magma_minproduct_smfree( &B, queue );
    magma_minproduct_free_cpu(coo_row);
    magma_minproduct_free_cpu(coo_col);
    magma_minproduct_free_cpu(coo_val);
    return info;
}


/**
    Purpose
    -------

    Reads in a SYMMETRIC matrix stored in coo format from a Matrix Market (.mtx)
    file and converts it into CSR format. It does not duplicate the off-diagonal
    entries!

    Arguments
    ---------

    @param[out]
    A           magma_minproduct_s_matrix*
                matrix in magma_minproduct sparse matrix format

    @param[in]
    filename    const char*
                filname of the mtx matrix
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_saux
    ********************************************************************/

extern "C"
magma_minproduct_int_t
magma_minproduct_s_csr_mtxsymm(
    magma_minproduct_s_matrix *A,
    const char *filename,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    magma_minproduct_s_matrix B={Magma_minproduct_CSR};
        
    int csr_compressor = 0;       // checks for zeros in original file
    
    magma_minproduct_index_t *coo_col=NULL, *coo_row=NULL;
    float *coo_val=NULL;

    FILE *fid;
    MM_typecode matcode;
      
    fid = fopen(filename, "r");
    
    if (fid == NULL) {
        printf("#Unable to open file %s\n", filename);
        info = MAGMA_minproduct_ERR_NOT_FOUND;
        goto cleanup;
    }
    
    if (mm_read_banner(fid, &matcode) != 0) {
        printf("#Could not process lMatrix Market banner.\n");
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
        goto cleanup;
    }
    
    if (!mm_is_valid(matcode)) {
        printf("#Invalid lMatrix Market file.\n");
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
        goto cleanup;
    }
    
    if ( ! ( (mm_is_real(matcode) || mm_is_integer(matcode)
           || mm_is_pattern(matcode) || mm_is_real(matcode) )
             && mm_is_coordinate(matcode)
             && mm_is_sparse(matcode) ) )
    {
        printf("#Sorry, this application does not support ");
        printf("#Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        printf("#Only real-valued or pattern coordinate matrices are supported\n");
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
        goto cleanup;
    }
    
    magma_minproduct_index_t num_rows, num_cols, num_nonzeros;
    if (mm_read_mtx_crd_size(fid, &num_rows, &num_cols, &num_nonzeros) != 0)
        info = MAGMA_minproduct_ERR_NOT_FOUND;
    
    A->storage_type    = Magma_minproduct_CSR;
    A->memory_location = Magma_minproduct_CPU;
    A->num_rows        = num_rows;
    A->num_cols        = num_cols;
    A->nnz             = num_nonzeros;
    A->fill_mode       = Magma_minproduct_FULL;
  
    
    
    CHECK( magma_minproduct_index_malloc_cpu( &coo_col, A->nnz ) );
    CHECK( magma_minproduct_index_malloc_cpu( &coo_row, A->nnz ) );
    CHECK( magma_minproduct_smalloc_cpu( &coo_val, A->nnz ) );
    
    printf("# Reading sparse matrix from file (%s):", filename);
    fflush(stdout);

    if (mm_is_real(matcode) || mm_is_integer(matcode)) {
        for(magma_minproduct_int_t i = 0; i < A->nnz; ++i) {
            magma_minproduct_index_t ROW, COL;
            float VAL;  // always read in a float and convert later if necessary
            
            fscanf(fid, " %d %d %f \n", &ROW, &COL, &VAL);
            if ( VAL == 0 )
                csr_compressor = 1;
            coo_row[i] = ROW - 1;
            coo_col[i] = COL - 1;
            coo_val[i] = MAGMA_minproduct_S_MAKE( VAL, 0.);
        }
    } else if (mm_is_pattern(matcode) ) {
        for(magma_minproduct_int_t i = 0; i < A->nnz; ++i) {
            magma_minproduct_index_t ROW, COL;
            
            fscanf(fid, " %d %d \n", &ROW, &COL);
            
            coo_row[i] = ROW - 1;
            coo_col[i] = COL - 1;
            coo_val[i] = MAGMA_minproduct_S_MAKE( 1.0, 0.);
        }
    } else if (mm_is_real(matcode) ){
       for(magma_minproduct_int_t i = 0; i < A->nnz; ++i) {
            magma_minproduct_index_t ROW, COL;
            float VAL, VALC;  // always read in a float and convert later if necessary
            
            fscanf(fid, " %d %d %f %f\n", &ROW, &COL, &VAL, &VALC);
            
            coo_row[i] = ROW - 1;
            coo_col[i] = COL - 1;
            coo_val[i] = MAGMA_minproduct_S_MAKE( VAL, VALC);
        }
        // printf(" ...successfully read real matrix... ");
    } else {
        printf("Unrecognized data type\n");
        info = MAGMA_minproduct_ERR_NOT_FOUND;
    }
    
    fclose(fid);
    printf(" done\n");
    
    A->sym = Magma_minproduct_GENERAL;

    if (mm_is_symmetric(matcode)) { // do not duplicate off diagonal entries!
        A->sym = Magma_minproduct_SYMMETRIC;
    } // end symmetric case
    
    float tv;
    magma_minproduct_index_t ti;
    
    // If matrix is not in standard format, sorting is necessary
    /*
    printf( "Sorting the cols....\n" );
    // bubble sort (by cols)
    for (int i=0; i < A->nnz-1; ++i) {
        for (int j=0; j < A->nnz-i-1; ++j) {
            if (coo_col[j] > coo_col[j+1] ) {
                ti = coo_col[j];
                coo_col[j] = coo_col[j+1];
                coo_col[j+1] = ti;
                
                ti = coo_row[j];
                coo_row[j] = coo_row[j+1];
                coo_row[j+1] = ti;
                
                tv = coo_val[j];
                coo_val[j] = coo_val[j+1];
                coo_val[j+1] = tv;
            }
        }
    }

    printf( "Sorting the rows....\n" );
    // bubble sort (by rows)
    for (int i=0; i < A->nnz-1; ++i) {
        for (int j=0; j < A->nnz-i-1; ++j) {
            if ( coo_row[j] > coo_row[j+1] ) {
                ti = coo_col[j];
                coo_col[j] = coo_col[j+1];
                coo_col[j+1] = ti;
                
                ti = coo_row[j];
                coo_row[j] = coo_row[j+1];
                coo_row[j+1] = ti;
                
                tv = coo_val[j];
                coo_val[j] = coo_val[j+1];
                coo_val[j+1] = tv;
            }
        }
    }
    printf( "Sorting: done\n" );
    
    */

    CHECK( magma_minproduct_index_malloc_cpu( &A->col, A->nnz ) );
    CHECK( magma_minproduct_index_malloc_cpu( &A->row, A->num_rows+1 ) );
    CHECK( magma_minproduct_smalloc_cpu( &A->val, A->nnz ) );

    // original code from  Nathan Bell and Michael Garland
    // the output CSR structure is NOT sorted!

    for (magma_minproduct_index_t i = 0; i < num_rows; i++)
        (A->row)[i] = 0;
    
    for (magma_minproduct_index_t i = 0; i < A->nnz; i++)
        (A->row)[coo_row[i]]++;
    
    // cumsum the nnz per row to get Bp[]
    for(magma_minproduct_int_t i = 0, cumsum = 0; i < num_rows; i++) {
        magma_minproduct_index_t temp = (A->row)[i];
        (A->row)[i] = cumsum;
        cumsum += temp;
    }
    (A->row)[num_rows] = A->nnz;
    
    // write Aj,Ax into Bj,Bx
    for(magma_minproduct_int_t i = 0; i < A->nnz; i++) {
        magma_minproduct_index_t row_  = coo_row[i];
        magma_minproduct_index_t dest = (A->row)[row_];
        
        (A->col)[dest] = coo_col[i];
        
        (A->val)[dest] = coo_val[i];
        
        (A->row)[row_]++;
    }
    magma_minproduct_free_cpu(coo_row);
    magma_minproduct_free_cpu(coo_col);
    magma_minproduct_free_cpu(coo_val);
    coo_row = NULL;
    coo_col = NULL;
    coo_val = NULL;
        
    for(int i = 0, last = 0; i <= num_rows; i++) {
        int temp    = (A->row)[i];
        (A->row)[i] = last;
        last        = temp;
    }
    
    (A->row)[A->num_rows]=A->nnz;
       
    for (magma_minproduct_index_t k=0; k < A->num_rows; ++k) {
        for (magma_minproduct_index_t i=(A->row)[k]; i < (A->row)[k+1]-1; ++i) {
            for (magma_minproduct_index_t j=(A->row)[k]; j < (A->row)[k+1]-1; ++j) {
                if ( (A->col)[j] > (A->col)[j+1] ) {
                    ti            = (A->col)[j];
                    (A->col)[j]   = (A->col)[j+1];
                    (A->col)[j+1] = ti;
                    
                    tv            = (A->val)[j];
                    (A->val)[j]   = (A->val)[j+1];
                    (A->val)[j+1] = tv;
                }
            }
        }
    }
    if ( csr_compressor > 0) { // run the CSR compressor to remove zeros
        //printf("removing zeros: ");
        CHECK( magma_minproduct_smtransfer( *A, &B, Magma_minproduct_CPU, Magma_minproduct_CPU, queue ));
        CHECK( magma_minproduct_s_csr_compressor(
            &(A->val), &(A->row), &(A->col),
            &B.val, &B.row, &B.col, &B.num_rows, queue ));
        B.nnz = B.row[num_rows];
        //printf(" remaining nonzeros:%d ", B.nnz);
        magma_minproduct_free_cpu( A->val );
        magma_minproduct_free_cpu( A->row );
        magma_minproduct_free_cpu( A->col );
        CHECK( magma_minproduct_smtransfer( B, A, Magma_minproduct_CPU, Magma_minproduct_CPU, queue ));

        //printf("done.\n");
    }
cleanup:
    magma_minproduct_smfree( &B, queue );
    magma_minproduct_free_cpu(coo_row);
    magma_minproduct_free_cpu(coo_col);
    magma_minproduct_free_cpu(coo_val);
    return info;
}
