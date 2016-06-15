/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zmergebicgstab2.cu normal z -> c, Sun May  3 11:22:58 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally3sparse.h"

#define BLOCK_SIZE 256

#define PRECISION_c


// These routines merge multiple kernels from cmergebicgstab into one
// This is the code used for the ASHES2014 paper
// "Accelerating Krylov Subspace Solvers on Graphics Processing Units".
// notice that only CSR format is supported so far.


// accelerated reduction for one vector
__global__ void
magma_tally3_creduce_kernel_spmv1(    
    int Gs,
    int n, 
    magma_tally3FloatComplex * vtmp,
    magma_tally3FloatComplex * vtmp2 )
{

    extern __shared__ magma_tally3FloatComplex temp[];    
    int Idx = threadIdx.x;
    int blockSize = 128;
    int gridSize = blockSize  * 2 * gridDim.x; 
    temp[Idx] = MAGMA_tally3_C_MAKE( 0.0, 0.0);
    int i = blockIdx.x * ( blockSize * 2 ) + Idx;   
    while (i < Gs ) {
        temp[ Idx  ] += vtmp[ i ]; 
        temp[ Idx  ] += ( i + blockSize < Gs ) ? vtmp[ i + blockSize ] 
                                                : MAGMA_tally3_C_MAKE( 0.0, 0.0); 
        i += gridSize;
    }
    __syncthreads();
    if ( Idx < 64 ){
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            temp[ Idx ] += temp[ Idx + 32 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 8 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 4 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 2 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 1 ];__syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    if ( Idx == 0 ){
        vtmp2[ blockIdx.x ] = temp[ 0 ];
    }
}


__global__ void
magma_tally3_cbicgmerge_spmv1_kernel(  
    int n,
    magma_tally3FloatComplex * dval, 
    magma_tally3_index_t * drowptr, 
    magma_tally3_index_t * dcolind,
    magma_tally3FloatComplex * p,
    magma_tally3FloatComplex * r,
    magma_tally3FloatComplex * v,
    magma_tally3FloatComplex * vtmp)
{

    extern __shared__ magma_tally3FloatComplex temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;
    int j;

    if( i<n ){
        magma_tally3FloatComplex dot = MAGMA_tally3_C_ZERO;
        int start = drowptr[ i ];
        int end = drowptr[ i+1 ];
        for( j=start; j<end; j++)
            dot += dval[ j ] * p[ dcolind[j] ];
        v[ i ] =  dot;
    }

    __syncthreads(); 

    temp[ Idx ] = ( i < n ) ? v[ i ] * r[ i ] : MAGMA_tally3_C_MAKE( 0.0, 0.0);
    __syncthreads();
    if ( Idx < 128 ){
        temp[ Idx ] += temp[ Idx + 128 ];
    }
    __syncthreads();
    if ( Idx < 64 ){
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            temp[ Idx ] += temp[ Idx + 32 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 8 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 4 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 2 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 1 ];__syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif

    if ( Idx == 0 ){
            vtmp[ blockIdx.x ] = temp[ 0 ];
    }
}

__global__ void
magma_tally3_cbicgstab_alphakernel(  
                    magma_tally3FloatComplex * skp ){
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if( i==0 ){
        magma_tally3FloatComplex tmp = skp[0];
        skp[0] = skp[4]/tmp;
    }
}

/**
    Purpose
    -------

    Merges the first SpmV using CSR with the dot product 
    and the computation of alpha

    Arguments
    ---------

    @param[in]
    A           magma_tally3_c_matrix
                system matrix

    @param[in]
    d1          magma_tally3FloatComplex_ptr
                temporary vector

    @param[in]
    d2          magma_tally3FloatComplex_ptr
                temporary vector

    @param[in]
    dp          magma_tally3FloatComplex_ptr
                input vector p

    @param[in]
    dr          magma_tally3FloatComplex_ptr
                input vector r

    @param[in]
    dv          magma_tally3FloatComplex_ptr
                output vector v

    @param[in/out]
    skp         magma_tally3FloatComplex_ptr
                array for parameters ( skp[0]=alpha )

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_cgegpuk
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_cbicgmerge_spmv1(
    magma_tally3_c_matrix A,
    magma_tally3FloatComplex_ptr d1,
    magma_tally3FloatComplex_ptr d2,
    magma_tally3FloatComplex_ptr dp,
    magma_tally3FloatComplex_ptr dr,
    magma_tally3FloatComplex_ptr dv,
    magma_tally3FloatComplex_ptr skp,
    magma_tally3_queue_t queue )
{
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue;
    magma_tally3blasGetKernelStream( &orig_queue );

    int n = A.num_rows;
    int local_block_size=256;
    dim3 Bs( local_block_size );
    dim3 Gs( magma_tally3_ceildiv( n, local_block_size ) );
    dim3 Gs_next;
    int Ms =  local_block_size * sizeof( magma_tally3FloatComplex ); 
    magma_tally3FloatComplex_ptr aux1 = d1, aux2 = d2;
    int b = 1;        

    if ( A.storage_type == Magma_tally3_CSR)
        magma_tally3_cbicgmerge_spmv1_kernel<<<Gs, Bs, Ms>>>
                    ( n, A.dval, A.drow, A.dcol, dp, dr, dv, d1 );
    else
        printf("error: only CSR format supported.\n");

    while( Gs.x > 1 ) {
        Gs_next.x = ( Gs.x+Bs.x-1 )/ Bs.x ;
        if ( Gs_next.x == 1 ) Gs_next.x = 2;
        magma_tally3_creduce_kernel_spmv1<<< Gs_next.x/2, Bs.x/2, Ms/2 >>> 
                            ( Gs.x, n, aux1, aux2 );
        Gs_next.x = Gs_next.x /2;
        Gs.x = Gs_next.x;
        b = 1 - b;
        if ( b ) { aux1 = d1; aux2 = d2; }
        else   { aux2 = d1; aux1 = d2; }
    }


    magma_tally3_ccopyvector( 1, aux1, 1, skp, 1 );
    dim3 Bs2( 2 );
    dim3 Gs2( 1 );
    magma_tally3_cbicgstab_alphakernel<<<Gs2, Bs2, 0>>>( skp );

   magma_tally3blasSetKernelStream( orig_queue );
   return MAGMA_tally3_SUCCESS;
}

/* -------------------------------------------------------------------------- */

// accelerated block reduction for multiple vectors
__global__ void
magma_tally3_creduce_kernel_spmv2( 
    int Gs,
    int n, 
    magma_tally3FloatComplex * vtmp,
    magma_tally3FloatComplex * vtmp2 )
{

    extern __shared__ magma_tally3FloatComplex temp[];    
    int Idx = threadIdx.x;
    int blockSize = 128;
    int gridSize = blockSize  * 2 * gridDim.x; 
    int j;

    for( j=0; j<2; j++){
        int i = blockIdx.x * ( blockSize * 2 ) + Idx;   
        temp[Idx+j*(blockSize)] = MAGMA_tally3_C_MAKE( 0.0, 0.0);
        while (i < Gs ) {
            temp[ Idx+j*(blockSize)  ] += vtmp[ i+j*n ]; 
            temp[ Idx+j*(blockSize)  ] += 
                ( i + (blockSize) < Gs ) ? vtmp[ i+j*n + (blockSize) ] 
                : MAGMA_tally3_C_MAKE( 0.0, 0.0); 
            i += gridSize;
        }
    }
    __syncthreads();
    if ( Idx < 64 ){
        for( j=0; j<2; j++){
            temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 64 ];
        }
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            for( j=0; j<2; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 32 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 16 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 8 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 4 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 2 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 1 ];
                __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            for( j=0; j<2; j++){
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 32 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 16 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 8 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 4 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 2 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 1 ];
            }
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            for( j=0; j<2; j++){
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 32 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 16 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 8 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 4 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 2 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 1 ];
            }
        }
    #endif
    if ( Idx == 0 ){
        for( j=0; j<2; j++){
            vtmp2[ blockIdx.x+j*n ] = temp[ j*(blockSize) ];
        }
    }
}

__global__ void
magma_tally3_cbicgmerge_spmv2_kernel(  
    int n,
    magma_tally3FloatComplex * dval, 
    magma_tally3_index_t * drowptr, 
    magma_tally3_index_t * dcolind,
    magma_tally3FloatComplex * s,
    magma_tally3FloatComplex * t,
    magma_tally3FloatComplex * vtmp )
{

    extern __shared__ magma_tally3FloatComplex temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;
    int j;

    if( i<n ){
        magma_tally3FloatComplex dot = MAGMA_tally3_C_ZERO;
        int start = drowptr[ i ];
        int end = drowptr[ i+1 ];
        for( j=start; j<end; j++)
            dot += dval[ j ] * s[ dcolind[j] ];
        t[ i ] =  dot;
    }

    __syncthreads(); 

    // 2 vectors 
    if (i<n){
            magma_tally3FloatComplex tmp2 = t[i];
            temp[Idx] = s[i] * tmp2;
            temp[Idx+blockDim.x] = tmp2 * tmp2;
    }
    else{
        for( j=0; j<2; j++)
            temp[Idx+j*blockDim.x] =MAGMA_tally3_C_MAKE( 0.0, 0.0);
    }
    __syncthreads();
    if ( Idx < 128 ){
        for( j=0; j<2; j++){
            temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 128 ];
        }
    }
    __syncthreads();
    if ( Idx < 64 ){
        for( j=0; j<2; j++){
            temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 64 ];
        }
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 32 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 16 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 8 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 4 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 2 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 1 ];
                __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            for( j=0; j<2; j++){
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 32 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 16 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 8 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 4 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 2 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 1 ];
            }
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            for( j=0; j<2; j++){
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 32 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 16 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 8 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 4 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 2 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 1 ];
            }
        }
    #endif
    if ( Idx == 0 ){
        for( j=0; j<2; j++){
            vtmp[ blockIdx.x+j*n ] = temp[ j*blockDim.x ];
        }
    }
}

__global__ void
magma_tally3_cbicgstab_omegakernel(  
                    magma_tally3FloatComplex * skp ){
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if( i==0 ){
        skp[2] = skp[6]/skp[7];
        skp[3] = skp[4];
    }
}

/**
    Purpose
    -------

    Merges the second SpmV using CSR with the dot product 
    and the computation of omega

    Arguments
    ---------

    @param[in]
    A           magma_tally3_c_matrix
                input matrix 

    @param[in]
    d1          magma_tally3FloatComplex_ptr
                temporary vector

    @param[in]
    d2          magma_tally3FloatComplex_ptr
                temporary vector

    @param[in]
    ds          magma_tally3FloatComplex_ptr
                input vector s

    @param[in]
    dt          magma_tally3FloatComplex_ptr
                output vector t

    @param[in/out]
    skp         magma_tally3FloatComplex_ptr
                array for parameters

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_cgegpuk
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_cbicgmerge_spmv2(
    magma_tally3_c_matrix A,
    magma_tally3FloatComplex_ptr d1,
    magma_tally3FloatComplex_ptr d2,
    magma_tally3FloatComplex_ptr ds,
    magma_tally3FloatComplex_ptr dt,
    magma_tally3FloatComplex_ptr skp,
    magma_tally3_queue_t queue )
{
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue;
    magma_tally3blasGetKernelStream( &orig_queue );

    int n = A.num_rows;
    int local_block_size=256;
    dim3 Bs( local_block_size );
    dim3 Gs( magma_tally3_ceildiv( n, local_block_size ) );
    dim3 Gs_next;
    int Ms =  2*local_block_size * sizeof( magma_tally3FloatComplex ); 
    magma_tally3FloatComplex_ptr aux1 = d1, aux2 = d2;
    int b = 1;        
    if ( A.storage_type == Magma_tally3_CSR)
        magma_tally3_cbicgmerge_spmv2_kernel<<<Gs, Bs, Ms>>>
                    ( n, A.dval, A.drow, A.dcol, ds, dt, d1 );
    else
        printf("error: only CSR format supported.\n");

    while( Gs.x > 1 ) {
        Gs_next.x = ( Gs.x+Bs.x-1 )/ Bs.x ;
        if ( Gs_next.x == 1 ) Gs_next.x = 2;
        magma_tally3_creduce_kernel_spmv2<<< Gs_next.x/2, Bs.x/2, Ms/2 >>> 
                    ( Gs.x, n, aux1, aux2 );
        Gs_next.x = Gs_next.x /2;
        Gs.x = Gs_next.x;
        b = 1 - b;
        if ( b ) { aux1 = d1; aux2 = d2; }
        else   { aux2 = d1; aux1 = d2; }
    }


    magma_tally3_ccopyvector( 1, aux1, 1, skp+6, 1 );
    magma_tally3_ccopyvector( 1, aux1+n, 1, skp+7, 1 );
    dim3 Bs2( 2 );
    dim3 Gs2( 1 );
    magma_tally3_cbicgstab_omegakernel<<<Gs2, Bs2, 0>>>( skp );

   magma_tally3blasSetKernelStream( orig_queue );
   return MAGMA_tally3_SUCCESS;
}

/* -------------------------------------------------------------------------- */

__global__ void
magma_tally3_cbicgmerge_xrbeta_kernel(  
    int n, 
    magma_tally3FloatComplex * rr,
    magma_tally3FloatComplex * r,
    magma_tally3FloatComplex * p,
    magma_tally3FloatComplex * s,
    magma_tally3FloatComplex * t,
    magma_tally3FloatComplex * x, 
    magma_tally3FloatComplex * skp,
    magma_tally3FloatComplex * vtmp )
{

    extern __shared__ magma_tally3FloatComplex temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;
    int j;

    magma_tally3FloatComplex alpha=skp[0];
    magma_tally3FloatComplex omega=skp[2];

    if( i<n ){
        magma_tally3FloatComplex sl;
        sl = s[i];
        x[i] = x[i] + alpha * p[i] + omega * sl;
        r[i] = sl - omega * t[i];
    }

    __syncthreads(); 

    // 2 vectors 
    if (i<n){
            magma_tally3FloatComplex tmp2 = r[i];
            temp[Idx] = rr[i] * tmp2;
            temp[Idx+blockDim.x] = tmp2 * tmp2;
    }
    else{
        for( j=0; j<2; j++)
            temp[Idx+j*blockDim.x] =MAGMA_tally3_C_MAKE( 0.0, 0.0);
    }
    __syncthreads();
    if ( Idx < 128 ){
        for( j=0; j<2; j++){
            temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 128 ];
        }
    }
    __syncthreads();
    if ( Idx < 64 ){
        for( j=0; j<2; j++){
            temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 64 ];
        }
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 32 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 16 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 8 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 4 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 2 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 1 ];
                __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            for( j=0; j<2; j++){
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 32 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 16 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 8 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 4 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 2 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 1 ];
            }
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            for( j=0; j<2; j++){
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 32 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 16 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 8 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 4 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 2 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 1 ];
            }
        }
    #endif
    if ( Idx == 0 ){
        for( j=0; j<2; j++){
            vtmp[ blockIdx.x+j*n ] = temp[ j*blockDim.x ];
        }
    }
}

__global__ void
magma_tally3_cbicgstab_betakernel(  
    magma_tally3FloatComplex * skp )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if( i==0 ){
        magma_tally3FloatComplex tmp1 = skp[4]/skp[3];
        magma_tally3FloatComplex tmp2 = skp[0] / skp[2];
        skp[1] =  tmp1*tmp2;
    }
}

/**
    Purpose
    -------

    Merges the second SpmV using CSR with the dot product 
    and the computation of omega

    Arguments
    ---------

    @param[in]
    n           int
                dimension n

    @param[in]
    d1          magma_tally3FloatComplex_ptr
                temporary vector

    @param[in]
    d2          magma_tally3FloatComplex_ptr
                temporary vector

    @param[in]
    rr          magma_tally3FloatComplex_ptr
                input vector rr

    @param[in]
    r           magma_tally3FloatComplex_ptr
                input/output vector r

    @param[in]
    p           magma_tally3FloatComplex_ptr
                input vector p

    @param[in]
    s           magma_tally3FloatComplex_ptr
                input vector s

    @param[in]
    t           magma_tally3FloatComplex_ptr
                input vector t

    @param[out]
    x           magma_tally3FloatComplex_ptr
                output vector x

    @param[in]
    skp         magma_tally3FloatComplex_ptr
                array for parameters

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_cgegpuk
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_cbicgmerge_xrbeta(
    int n,
    magma_tally3FloatComplex_ptr d1,
    magma_tally3FloatComplex_ptr d2,
    magma_tally3FloatComplex_ptr rr,
    magma_tally3FloatComplex_ptr r,
    magma_tally3FloatComplex_ptr p,
    magma_tally3FloatComplex_ptr s,
    magma_tally3FloatComplex_ptr t,
    magma_tally3FloatComplex_ptr x, 
    magma_tally3FloatComplex_ptr skp,
    magma_tally3_queue_t queue )
{
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue;
    magma_tally3blasGetKernelStream( &orig_queue );

    int local_block_size=256;
    dim3 Bs( local_block_size );
    dim3 Gs( magma_tally3_ceildiv( n, local_block_size ) );
    dim3 Gs_next;
    int Ms =  2*local_block_size * sizeof( magma_tally3FloatComplex ); 
    magma_tally3FloatComplex_ptr aux1 = d1, aux2 = d2;
    int b = 1;        
    magma_tally3_cbicgmerge_xrbeta_kernel<<<Gs, Bs, Ms>>>
                    ( n, rr, r, p, s, t, x, skp, d1);  

    while( Gs.x > 1 ) {
        Gs_next.x = ( Gs.x+Bs.x-1 )/ Bs.x ;
        if ( Gs_next.x == 1 ) Gs_next.x = 2;
        magma_tally3_creduce_kernel_spmv2<<< Gs_next.x/2, Bs.x/2, Ms/2 >>> 
                            ( Gs.x, n, aux1, aux2 );
        Gs_next.x = Gs_next.x /2;
        Gs.x = Gs_next.x;
        b = 1 - b;
        if ( b ) { aux1 = d1; aux2 = d2; }
        else   { aux2 = d1; aux1 = d2; }
    }


    magma_tally3_ccopyvector( 1, aux1, 1, skp+4, 1 );
    magma_tally3_ccopyvector( 1, aux1+n, 1, skp+5, 1 );
    dim3 Bs2( 2 );
    dim3 Gs2( 1 );
    magma_tally3_cbicgstab_betakernel<<<Gs2, Bs2, 0>>>( skp );

   magma_tally3blasSetKernelStream( orig_queue );
   return MAGMA_tally3_SUCCESS;
}

/* -------------------------------------------------------------------------- */
