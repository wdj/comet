/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zbajac_csr.cu normal z -> c, Sun May  3 11:22:57 2015

*/
#include "common_magma_tally4sparse.h"

#define PRECISION_c
#define BLOCKSIZE 256


__global__ void
magma_tally4_cbajac_csr_ls_kernel(int localiters, int n, 
                            magma_tally4FloatComplex * valD, 
                            magma_tally4_index_t * rowD, 
                            magma_tally4_index_t * colD, 
                            magma_tally4FloatComplex * valR, 
                            magma_tally4_index_t * rowR,
                            magma_tally4_index_t * colR, 
                            const magma_tally4FloatComplex *  __restrict__ b,                            
                            magma_tally4FloatComplex * x ){

    int inddiag =  blockIdx.x*blockDim.x;
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    int i, j, start, end;   


    if(index<n){
    
        start=rowR[index];
        end  =rowR[index+1];

        magma_tally4FloatComplex zero = MAGMA_tally4_C_MAKE(0.0, 0.0);
        magma_tally4FloatComplex bl, tmp = zero, v = zero; 

#if (__CUDA_ARCH__ >= 350) && (defined(PRECISION_d) || defined(PRECISION_s))
        bl = __ldg( b+index );
#else
        bl = b[index];
#endif

        #pragma unroll
        for( i=start; i<end; i++ )
             v += valR[i] * x[ colR[i] ];

        start=rowD[index];
        end  =rowD[index+1];

        #pragma unroll
        for( i=start; i<end; i++ )
            tmp += valD[i] * x[ colD[i] ];

        v =  bl - v;

        /* add more local iterations */           
        __shared__ magma_tally4FloatComplex local_x[ BLOCKSIZE ];
        local_x[threadIdx.x] = x[index] + ( v - tmp) / (valD[start]);
        __syncthreads();

        #pragma unroll
        for( j=0; j<localiters; j++ )
        {
            tmp = zero;
            #pragma unroll
            for( i=start; i<end; i++ )
                tmp += valD[i] * local_x[ colD[i] - inddiag];

            local_x[threadIdx.x] +=  ( v - tmp) / (valD[start]);
        }
        x[index] = local_x[threadIdx.x];
    }
}



__global__ void
magma_tally4_cbajac_csr_kernel(    
    int n, 
    magma_tally4FloatComplex * valD, 
    magma_tally4_index_t * rowD, 
    magma_tally4_index_t * colD, 
    magma_tally4FloatComplex * valR, 
    magma_tally4_index_t * rowR,
    magma_tally4_index_t * colR, 
    magma_tally4FloatComplex * b,                                
    magma_tally4FloatComplex * x ){

    int index = blockIdx.x*blockDim.x+threadIdx.x;
    int i, start, end;   

    if(index<n){
        
        magma_tally4FloatComplex zero = MAGMA_tally4_C_MAKE(0.0, 0.0);
        magma_tally4FloatComplex bl, tmp = zero, v = zero; 

#if (__CUDA_ARCH__ >= 350) && (defined(PRECISION_d) || defined(PRECISION_s))
        bl = __ldg( b+index );
#else
        bl = b[index];
#endif

        start=rowR[index];
        end  =rowR[index+1];

        #pragma unroll
        for( i=start; i<end; i++ )
             v += valR[i] * x[ colR[i] ];

        v =  bl - v;

        start=rowD[index];
        end  =rowD[index+1];

        #pragma unroll
        for( i=start; i<end; i++ )
            tmp += valD[i] * x[ colD[i] ];

        x[index] = x[index] + ( v - tmp ) / (valD[start]); 
    }
}









/**
    Purpose
    -------
    
    This routine is a block-asynchronous Jacobi iteration performing s
    local Jacobi-updates within the block. Input format is two CSR matrices,
    one containing the diagonal blocks, one containing the rest.

    Arguments
    ---------

    @param[in]
    localiters  magma_tally4_int_t
                number of local Jacobi-like updates

    @param[in]
    D           magma_tally4_c_matrix
                input matrix with diagonal blocks

    @param[in]
    R           magma_tally4_c_matrix
                input matrix with non-diagonal parts

    @param[in]
    b           magma_tally4_c_matrix
                RHS

    @param[in]
    x           magma_tally4_c_matrix*
                iterate/solution

    
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_cgegpuk
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_cbajac_csr(
    magma_tally4_int_t localiters,
    magma_tally4_c_matrix D,
    magma_tally4_c_matrix R,
    magma_tally4_c_matrix b,
    magma_tally4_c_matrix *x,
    magma_tally4_queue_t queue )
{
    int blocksize1 = BLOCKSIZE;
    int blocksize2 = 1;

    int dimgrid1 = magma_tally4_ceildiv(  D.num_rows, blocksize1 );
    int dimgrid2 = 1;
    int dimgrid3 = 1;

    dim3 grid( dimgrid1, dimgrid2, dimgrid3 );
    dim3 block( blocksize1, blocksize2, 1 );
    if ( R.nnz > 0 ) { 
        if ( localiters == 1 )
        magma_tally4_cbajac_csr_kernel<<< grid, block, 0, queue >>>
            ( D.num_rows, D.dval, D.drow, D.dcol, 
                            R.dval, R.drow, R.dcol, b.dval, x->dval );
        else
            magma_tally4_cbajac_csr_ls_kernel<<< grid, block, 0, queue >>>
            ( localiters, D.num_rows, D.dval, D.drow, D.dcol, 
                            R.dval, R.drow, R.dcol, b.dval, x->dval );
    }
    else {
        printf("error: all elements in diagonal block.\n");
    }

    return MAGMA_tally4_SUCCESS;
}


