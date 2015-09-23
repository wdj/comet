/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zbajac_csr.cu normal z -> s, Sun May  3 11:22:58 2015

*/
#include "common_magma_minproductsparse.h"

#define PRECISION_s
#define BLOCKSIZE 256


__global__ void
magma_minproduct_sbajac_csr_ls_kernel(int localiters, int n, 
                            float * valD, 
                            magma_minproduct_index_t * rowD, 
                            magma_minproduct_index_t * colD, 
                            float * valR, 
                            magma_minproduct_index_t * rowR,
                            magma_minproduct_index_t * colR, 
                            const float *  __restrict__ b,                            
                            float * x ){

    int inddiag =  blockIdx.x*blockDim.x;
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    int i, j, start, end;   


    if(index<n){
    
        start=rowR[index];
        end  =rowR[index+1];

        float zero = MAGMA_minproduct_S_MAKE(0.0, 0.0);
        float bl, tmp = zero, v = zero; 

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
        __shared__ float local_x[ BLOCKSIZE ];
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
magma_minproduct_sbajac_csr_kernel(    
    int n, 
    float * valD, 
    magma_minproduct_index_t * rowD, 
    magma_minproduct_index_t * colD, 
    float * valR, 
    magma_minproduct_index_t * rowR,
    magma_minproduct_index_t * colR, 
    float * b,                                
    float * x ){

    int index = blockIdx.x*blockDim.x+threadIdx.x;
    int i, start, end;   

    if(index<n){
        
        float zero = MAGMA_minproduct_S_MAKE(0.0, 0.0);
        float bl, tmp = zero, v = zero; 

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
    localiters  magma_minproduct_int_t
                number of local Jacobi-like updates

    @param[in]
    D           magma_minproduct_s_matrix
                input matrix with diagonal blocks

    @param[in]
    R           magma_minproduct_s_matrix
                input matrix with non-diagonal parts

    @param[in]
    b           magma_minproduct_s_matrix
                RHS

    @param[in]
    x           magma_minproduct_s_matrix*
                iterate/solution

    
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_sgegpuk
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_sbajac_csr(
    magma_minproduct_int_t localiters,
    magma_minproduct_s_matrix D,
    magma_minproduct_s_matrix R,
    magma_minproduct_s_matrix b,
    magma_minproduct_s_matrix *x,
    magma_minproduct_queue_t queue )
{
    int blocksize1 = BLOCKSIZE;
    int blocksize2 = 1;

    int dimgrid1 = magma_minproduct_ceildiv(  D.num_rows, blocksize1 );
    int dimgrid2 = 1;
    int dimgrid3 = 1;

    dim3 grid( dimgrid1, dimgrid2, dimgrid3 );
    dim3 block( blocksize1, blocksize2, 1 );
    if ( R.nnz > 0 ) { 
        if ( localiters == 1 )
        magma_minproduct_sbajac_csr_kernel<<< grid, block, 0, queue >>>
            ( D.num_rows, D.dval, D.drow, D.dcol, 
                            R.dval, R.drow, R.dcol, b.dval, x->dval );
        else
            magma_minproduct_sbajac_csr_ls_kernel<<< grid, block, 0, queue >>>
            ( localiters, D.num_rows, D.dval, D.drow, D.dcol, 
                            R.dval, R.drow, R.dcol, b.dval, x->dval );
    }
    else {
        printf("error: all elements in diagonal block.\n");
    }

    return MAGMA_minproduct_SUCCESS;
}



