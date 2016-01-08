/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zpipelinedgmres.cu normal z -> s, Sun May  3 11:22:58 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally4.h"

#define REAL

#define BLOCK_SIZE 512


template< int n >
__device__ void sum_reduce( /*int n,*/ int i, float* x )
{
    __syncthreads();
    if ( n > 1024 ) { if ( i < 1024 && i + 1024 < n ) { x[i] += x[i+1024]; }  
        __syncthreads(); }
    if ( n >  512 ) { if ( i <  512 && i +  512 < n ) { x[i] += x[i+ 512]; }  
        __syncthreads(); }
    if ( n >  256 ) { if ( i <  256 && i +  256 < n ) { x[i] += x[i+ 256]; }  
        __syncthreads(); }
    if ( n >  128 ) { if ( i <  128 && i +  128 < n ) { x[i] += x[i+ 128]; }  
        __syncthreads(); }
    if ( n >   64 ) { if ( i <   64 && i +   64 < n ) { x[i] += x[i+  64]; }  
        __syncthreads(); }
    if ( n >   32 ) { if ( i <   32 && i +   32 < n ) { x[i] += x[i+  32]; }  
        __syncthreads(); }
    // probably don't need __syncthreads for < 16 threads
    // because of implicit warp level synchronization.
    if ( n >   16 ) { if ( i <   16 && i +   16 < n ) { x[i] += x[i+  16]; }  
        __syncthreads(); }
    if ( n >    8 ) { if ( i <    8 && i +    8 < n ) { x[i] += x[i+   8]; }  
        __syncthreads(); }
    if ( n >    4 ) { if ( i <    4 && i +    4 < n ) { x[i] += x[i+   4]; }  
        __syncthreads(); }
    if ( n >    2 ) { if ( i <    2 && i +    2 < n ) { x[i] += x[i+   2]; }  
        __syncthreads(); }
    if ( n >    1 ) { if ( i <    1 && i +    1 < n ) { x[i] += x[i+   1]; }  
        __syncthreads(); }
}

__global__ void
magma_tally4_spipelined_correction( 
    int n,  
    int k,
    float * skp, 
    float * r,
    float * v )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float zz= 0.0, tmp= 0.0;

    extern __shared__ float temp[];    
    
    temp[ i ] = ( i < k ) ? skp[ i ] * skp[ i ] : MAGMA_tally4_S_MAKE( 0.0, 0.0);
    __syncthreads();
     if (i < 64) { temp[ i ] += temp[ i + 64 ]; } __syncthreads(); 
     if( i < 32 ){
        temp[ i ] += temp[ i + 32 ];__syncthreads();    
        temp[ i ] += temp[ i + 16 ];__syncthreads(); 
        temp[ i ] += temp[ i +  8 ];__syncthreads(); 
        temp[ i ] += temp[ i +  4 ];__syncthreads(); 
        temp[ i ] += temp[ i +  2 ];__syncthreads(); 
        temp[ i ] += temp[ i +  1 ];__syncthreads();      
    }
    if( i == 0 ){
        tmp = MAGMA_tally4_S_REAL( temp[ i ] );
        zz = MAGMA_tally4_S_REAL( skp[(k)] );
        skp[k] = MAGMA_tally4_S_MAKE( sqrt(zz-tmp),0.0 );
    }
}

__global__ void
magma_tally4_spipelined_copyscale( 
    int n,  
    int k,
    float * skp, 
    float * r,
    float * v )
{

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    float rr=skp[k];

    if( i<n ){
        v[i] =  r[i] * 1.0 / rr;

    }
}

//----------------------------------------------------------------------------//

__global__ void
magma_tally4_spipelinedsnrm2_kernel( 
    int m, 
    float * da, 
    int ldda, 
    float * dxnorm )
{
    const int i = threadIdx.x;
    magma_tally4Float_ptr dx = da + blockIdx.x * ldda;

    __shared__ float sum[ 512 ];
    float re, lsum;

    // get norm of dx
    lsum = 0;
    for( int j = i; j < m; j += 512 ) {
        #ifdef REAL
            re = dx[j];
            lsum += re*re;
        #else
            re = MAGMA_tally4_S_REAL( dx[j] );
            float im = MAGMA_tally4_S_IMAG( dx[j] );
            lsum += re*re + im*im;
        #endif
    }
    sum[i] = lsum;
    sum_reduce< 512 >( i, sum );

    if (i==0)
        dxnorm[blockIdx.x] = MAGMA_tally4_S_MAKE( sqrt(sum[0]), 0.0 );
}

//----------------------------------------------------------------------------//

__global__ void
magma_tally4_spipelinesscale( 
    int n, 
    float * r, 
    float * drnorm )
{

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if( i<n ){
        r[i] =  r[i] * 1.0 / drnorm[0];
    }
}

/**
    Purpose
    -------

    Computes the correction term of the pipelined GMRES according to P. Ghysels 
    and scales and copies the new search direction
    
    Returns the vector v = r/ ( skp[k] - (sum_i=1^k skp[i]^2) ) .

    Arguments
    ---------

    @param[in]
    n           int
                length of v_i

    @param[in]
    k           int
                # skp entries v_i^T * r ( without r )

    @param[in]
    r           magma_tally4Float_ptr 
                vector of length n

    @param[in]
    v           magma_tally4Float_ptr 
                vector of length n
                
    @param[in]  
    skp         magma_tally4Float_ptr 
                array of parameters

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_saux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_scopyscale(
    int n, 
    int k,
    magma_tally4Float_ptr r, 
    magma_tally4Float_ptr v,
    magma_tally4Float_ptr skp,
    magma_tally4_queue_t queue )
{
    dim3 Bs( BLOCK_SIZE );
    dim3 Gs( magma_tally4_ceildiv( k, BLOCK_SIZE ) );
    unsigned int Ms =   Bs.x * sizeof( float ); 

    dim3 Gs2( magma_tally4_ceildiv( n, BLOCK_SIZE ) );


    magma_tally4_spipelined_correction<<<Gs, Bs, Ms, queue >>>
                                            ( n, k, skp, r, v );
    magma_tally4_spipelined_copyscale<<<Gs2, Bs, 0, queue >>>
                                            ( n, k, skp, r, v );

    return MAGMA_tally4_SUCCESS;
}


extern "C" magma_tally4_int_t
magma_tally4_snrm2scale(
    int m, 
    magma_tally4Float_ptr r, 
    int lddr, 
    magma_tally4Float_ptr drnorm,
    magma_tally4_queue_t queue )
{
    dim3  blocks( 1 );
    dim3 threads( 512 );
    magma_tally4_spipelinedsnrm2_kernel<<< blocks, threads, 0, queue >>>
                                ( m, r, lddr, drnorm );

    dim3 Bs( BLOCK_SIZE );
    dim3 Gs2( magma_tally4_ceildiv( m, BLOCK_SIZE ) );
    magma_tally4_spipelinesscale<<<Gs2, Bs, 0, queue >>>( m, r, drnorm );

    return MAGMA_tally4_SUCCESS;
}

