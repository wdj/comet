/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s
       @author Hartwig Anzt

*/
#include "common_magma_tally3.h"

#define COMPLEX

#define BLOCK_SIZE 512


template< int n >
__device__ void sum_reduce( /*int n,*/ int i, double* x )
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
magma_tally3_zpipelined_correction( 
    int n,  
    int k,
    magma_tally3DoubleComplex * skp, 
    magma_tally3DoubleComplex * r,
    magma_tally3DoubleComplex * v )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    double zz= 0.0, tmp= 0.0;

    extern __shared__ magma_tally3DoubleComplex temp[];    
    
    temp[ i ] = ( i < k ) ? skp[ i ] * skp[ i ] : MAGMA_tally3_Z_MAKE( 0.0, 0.0);
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
        tmp = MAGMA_tally3_Z_REAL( temp[ i ] );
        zz = MAGMA_tally3_Z_REAL( skp[(k)] );
        skp[k] = MAGMA_tally3_Z_MAKE( sqrt(zz-tmp),0.0 );
    }
}

__global__ void
magma_tally3_zpipelined_copyscale( 
    int n,  
    int k,
    magma_tally3DoubleComplex * skp, 
    magma_tally3DoubleComplex * r,
    magma_tally3DoubleComplex * v )
{

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    magma_tally3DoubleComplex rr=skp[k];

    if( i<n ){
        v[i] =  r[i] * 1.0 / rr;

    }
}

//----------------------------------------------------------------------------//

__global__ void
magma_tally3_zpipelineddznrm2_kernel( 
    int m, 
    magma_tally3DoubleComplex * da, 
    int ldda, 
    magma_tally3DoubleComplex * dxnorm )
{
    const int i = threadIdx.x;
    magma_tally3DoubleComplex_ptr dx = da + blockIdx.x * ldda;

    __shared__ double sum[ 512 ];
    double re, lsum;

    // get norm of dx
    lsum = 0;
    for( int j = i; j < m; j += 512 ) {
        #ifdef REAL
            re = dx[j];
            lsum += re*re;
        #else
            re = MAGMA_tally3_Z_REAL( dx[j] );
            double im = MAGMA_tally3_Z_IMAG( dx[j] );
            lsum += re*re + im*im;
        #endif
    }
    sum[i] = lsum;
    sum_reduce< 512 >( i, sum );

    if (i==0)
        dxnorm[blockIdx.x] = MAGMA_tally3_Z_MAKE( sqrt(sum[0]), 0.0 );
}

//----------------------------------------------------------------------------//

__global__ void
magma_tally3_zpipelinedscale( 
    int n, 
    magma_tally3DoubleComplex * r, 
    magma_tally3DoubleComplex * drnorm )
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
    r           magma_tally3DoubleComplex_ptr 
                vector of length n

    @param[in]
    v           magma_tally3DoubleComplex_ptr 
                vector of length n
                
    @param[in]  
    skp         magma_tally3DoubleComplex_ptr 
                array of parameters

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zaux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_zcopyscale(
    int n, 
    int k,
    magma_tally3DoubleComplex_ptr r, 
    magma_tally3DoubleComplex_ptr v,
    magma_tally3DoubleComplex_ptr skp,
    magma_tally3_queue_t queue )
{
    dim3 Bs( BLOCK_SIZE );
    dim3 Gs( magma_tally3_ceildiv( k, BLOCK_SIZE ) );
    unsigned int Ms =   Bs.x * sizeof( magma_tally3DoubleComplex ); 

    dim3 Gs2( magma_tally3_ceildiv( n, BLOCK_SIZE ) );


    magma_tally3_zpipelined_correction<<<Gs, Bs, Ms, queue >>>
                                            ( n, k, skp, r, v );
    magma_tally3_zpipelined_copyscale<<<Gs2, Bs, 0, queue >>>
                                            ( n, k, skp, r, v );

    return MAGMA_tally3_SUCCESS;
}


extern "C" magma_tally3_int_t
magma_tally3_dznrm2scale(
    int m, 
    magma_tally3DoubleComplex_ptr r, 
    int lddr, 
    magma_tally3DoubleComplex_ptr drnorm,
    magma_tally3_queue_t queue )
{
    dim3  blocks( 1 );
    dim3 threads( 512 );
    magma_tally3_zpipelineddznrm2_kernel<<< blocks, threads, 0, queue >>>
                                ( m, r, lddr, drnorm );

    dim3 Bs( BLOCK_SIZE );
    dim3 Gs2( magma_tally3_ceildiv( m, BLOCK_SIZE ) );
    magma_tally3_zpipelinedscale<<<Gs2, Bs, 0, queue >>>( m, r, drnorm );

    return MAGMA_tally3_SUCCESS;
}
