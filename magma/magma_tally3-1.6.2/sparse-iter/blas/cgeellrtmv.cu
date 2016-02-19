/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zgeellrtmv.cu normal z -> c, Sun May  3 11:22:58 2015

*/

#include "common_magma_tally3.h"

//F. Vázquez, G. Ortega, J.J. Fernández, E.M. Garzón, Almeria University
__global__ void 
cgeellrtmv_kernel_32( 
    int num_rows, 
    int num_cols,
    magma_tally3FloatComplex alpha, 
    magma_tally3FloatComplex * dval, 
    magma_tally3_index_t * dcolind,
    magma_tally3_index_t * drowlength,
    magma_tally3FloatComplex * dx,
    magma_tally3FloatComplex beta, 
    magma_tally3FloatComplex * dy,
    int T,
    int alignment )
{
int idx = blockIdx.y * gridDim.x * blockDim.x + 
          blockDim.x * blockIdx.x + threadIdx.x ; // global thread index
int idb = threadIdx.x ;  // local thread index
int idp = idb%T;  // number of threads assigned to one row
int i = idx/T;  // row index

extern __shared__ magma_tally3FloatComplex shared[];

    if(i < num_rows ){
        magma_tally3FloatComplex dot = MAGMA_tally3_C_MAKE(0.0, 0.0);
        int max_ = magma_tally3_ceildiv( drowlength[i], T );  
            // number of elements each thread handles

        for ( int k = 0; k < max_ ; k++ ){

            // original code in paper (not working for me)
            //magma_tally3FloatComplex val = dval[ k*(T*alignment)+(i*T)+idp ];  
            //int col = dcolind [ k*(T*alignment)+(i*T)+idp ];    

            // new code (working for me)        
            magma_tally3FloatComplex val = dval[ k*(T)+(i*alignment)+idp ];
            int col = dcolind [ k*(T)+(i*alignment)+idp ];

            dot += val * dx[ col ];
        }
        shared[idb]  = dot;
        if( idp < 16 ){
            shared[idb]+=shared[idb+16];
            if( idp < 8 ) shared[idb]+=shared[idb+8];
            if( idp < 4 ) shared[idb]+=shared[idb+4];
            if( idp < 2 ) shared[idb]+=shared[idb+2];
            if( idp == 0 ) {
                dy[i] = (shared[idb]+shared[idb+1])*alpha + beta*dy [i];
            }

        }
    }

}

//F. Vázquez, G. Ortega, J.J. Fernández, E.M. Garzón, Almeria University
__global__ void 
cgeellrtmv_kernel_16( 
    int num_rows, 
    int num_cols,
    magma_tally3FloatComplex alpha, 
    magma_tally3FloatComplex * dval, 
    magma_tally3_index_t * dcolind,
    magma_tally3_index_t * drowlength,
    magma_tally3FloatComplex * dx,
    magma_tally3FloatComplex beta, 
    magma_tally3FloatComplex * dy,
    int T,
    int alignment )
{
int idx = blockIdx.y * gridDim.x * blockDim.x + 
          blockDim.x * blockIdx.x + threadIdx.x ; // global thread index
int idb = threadIdx.x ;  // local thread index
int idp = idb%T;  // number of threads assigned to one row
int i = idx/T;  // row index

extern __shared__ magma_tally3FloatComplex shared[];

    if(i < num_rows ){
        magma_tally3FloatComplex dot = MAGMA_tally3_C_MAKE(0.0, 0.0);
        int max_ = magma_tally3_ceildiv( drowlength[i], T );  
            // number of elements each thread handles

        for ( int k = 0; k < max_ ; k++ ){

            // original code in paper (not working for me)
            //magma_tally3FloatComplex val = dval[ k*(T*alignment)+(i*T)+idp ];  
            //int col = dcolind [ k*(T*alignment)+(i*T)+idp ];    

            // new code (working for me)        
            magma_tally3FloatComplex val = dval[ k*(T)+(i*alignment)+idp ];
            int col = dcolind [ k*(T)+(i*alignment)+idp ];

            dot += val * dx[ col ];
        }
        shared[idb]  = dot;
        if( idp < 8 ){
            shared[idb]+=shared[idb+8];
            if( idp < 4 ) shared[idb]+=shared[idb+4];
            if( idp < 2 ) shared[idb]+=shared[idb+2];
            if( idp == 0 ) {
                dy[i] = (shared[idb]+shared[idb+1])*alpha + beta*dy [i];
            }

        }
    }

}

//F. Vázquez, G. Ortega, J.J. Fernández, E.M. Garzón, Almeria University
__global__ void 
cgeellrtmv_kernel_8( 
    int num_rows, 
    int num_cols,
    magma_tally3FloatComplex alpha, 
    magma_tally3FloatComplex * dval, 
    magma_tally3_index_t * dcolind,
    magma_tally3_index_t * drowlength,
    magma_tally3FloatComplex * dx,
    magma_tally3FloatComplex beta, 
    magma_tally3FloatComplex * dy,
    int T,
    int alignment )
{
int idx = blockIdx.y * gridDim.x * blockDim.x + 
          blockDim.x * blockIdx.x + threadIdx.x ; // global thread index
int idb = threadIdx.x ;  // local thread index
int idp = idb%T;  // number of threads assigned to one row
int i = idx/T;  // row index

extern __shared__ magma_tally3FloatComplex shared[];

    if(i < num_rows ){
        magma_tally3FloatComplex dot = MAGMA_tally3_C_MAKE(0.0, 0.0);
        int max_ = magma_tally3_ceildiv( drowlength[i], T );  
            // number of elements each thread handles

        for ( int k = 0; k < max_ ; k++ ){

            // original code in paper (not working for me)
            //magma_tally3FloatComplex val = dval[ k*(T*alignment)+(i*T)+idp ];  
            //int col = dcolind [ k*(T*alignment)+(i*T)+idp ];    

            // new code (working for me)        
            magma_tally3FloatComplex val = dval[ k*(T)+(i*alignment)+idp ];
            int col = dcolind [ k*(T)+(i*alignment)+idp ];

            dot += val * dx[ col ];
        }
        shared[idb]  = dot;
        if( idp < 4 ){
            shared[idb]+=shared[idb+4];
            if( idp < 2 ) shared[idb]+=shared[idb+2];
            if( idp == 0 ) {
                dy[i] = (shared[idb]+shared[idb+1])*alpha + beta*dy [i];
            }

        }
    }

}



/**
    Purpose
    -------
    
    This routine computes y = alpha *  A *  x + beta * y on the GPU.
    Input format is ELLRT. The ideas are taken from 
    "Improving the performance of the sparse matrix
    vector product with GPUs", (CIT 2010), 
    and modified to provide correct values.

    
    Arguments
    ---------

    @param[in]
    transA      magma_tally3_trans_t
                transposition parameter for A
    @param[in]
    m           magma_tally3_int_t
                number of rows 

    @param[in]
    n           magma_tally3_int_t
                number of columns

    @param[in]
    nnz_per_row magma_tally3_int_t
                max number of nonzeros in a row

    @param[in]
    alpha       magma_tally3FloatComplex
                scalar alpha

    @param[in]
    dval        magma_tally3FloatComplex_ptr
                val array

    @param[in]
    dcolind     magma_tally3Index_ptr
                col indices  

    @param[in]
    drowlength  magma_tally3Index_ptr
                number of elements in each row

    @param[in]
    dx          magma_tally3FloatComplex_ptr
                input vector x

    @param[in]
    beta        magma_tally3FloatComplex
                scalar beta

    @param[out]
    dy          magma_tally3FloatComplex_ptr
                output vector y

    @param[in]
    blocksize   magma_tally3_int_t
                threads per block

    @param[in]
    alignment   magma_tally3_int_t
                threads assigned to each row

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_cblas
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_cgeellrtmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowlength,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dy,
    magma_tally3_int_t alignment,
    magma_tally3_int_t blocksize,
    magma_tally3_queue_t queue )
{
    int num_blocks = magma_tally3_ceildiv( m, blocksize );

    magma_tally3_int_t num_threads = alignment*blocksize;
    magma_tally3_int_t threads = alignment*blocksize;

    int real_row_length = magma_tally3_roundup( nnz_per_row, alignment );

    magma_tally3_int_t arch = magma_tally3_getdevice_arch();
    if ( arch < 200 && num_threads > 256 )
        printf("error: too much shared memory requested.\n");

    int dimgrid1 = (int) sqrt( (float) num_blocks );
    int dimgrid2 = magma_tally3_ceildiv( num_blocks, dimgrid1 );
    dim3 grid( dimgrid1, dimgrid2, 1);

    int Ms = alignment * blocksize * sizeof( magma_tally3FloatComplex );
    // printf("launch kernel: %dx%d %d %d\n", grid.x, grid.y, num_threads , Ms);

    if ( alignment == 32 ) {
        cgeellrtmv_kernel_32<<< grid, threads , Ms, queue >>>
                 ( m, n, alpha, dval, dcolind, drowlength, dx, beta, dy, 
                                                 alignment, real_row_length );
    }
    else if ( alignment == 16 ) {
        cgeellrtmv_kernel_16<<< grid, threads , Ms, queue >>>
                 ( m, n, alpha, dval, dcolind, drowlength, dx, beta, dy, 
                                                 alignment, real_row_length );
    }
    else if ( alignment == 8 ) {
        cgeellrtmv_kernel_8<<< grid, threads , Ms, queue >>>
                 ( m, n, alpha, dval, dcolind, drowlength, dx, beta, dy, 
                                                 alignment, real_row_length );
    }
    else {
        printf("error: alignment %d not supported.\n", alignment);
        return MAGMA_tally3_ERR_NOT_SUPPORTED;
    }



   return MAGMA_tally3_SUCCESS;
}


