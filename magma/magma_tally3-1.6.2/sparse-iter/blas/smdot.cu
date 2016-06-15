/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zmdot.cu normal z -> s, Sun May  3 11:22:58 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally3.h"

#define BLOCK_SIZE 256

#define PRECISION_s


// initialize arrays with zero
__global__ void
magma_tally3_sgpumemzero(  
    float * d, 
    int n, 
    int k )
{

   int i = blockIdx.x * blockDim.x + threadIdx.x;

   if( i < n ){
    for( int j=0; j<k; j++)
      d[ i+j*n ] = MAGMA_tally3_S_MAKE( 0.0, 0.0 );
    }
}

// dot product
__global__ void
magma_tally3_sdot_kernel( 
    int Gs,
    int n, 
    float * v,
    float * r,
    float * vtmp)
{

    extern __shared__ float temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;

    temp[ Idx ] = ( i < n ) ? v[ i ] * r[ i ] : MAGMA_tally3_S_MAKE( 0.0, 0.0);
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

// dot product for multiple vectors
__global__ void
magma_tally3_sblockdot_kernel( 
    int Gs,
    int n, 
    int k,
    float * v,
    float * r,
    float * vtmp)
{

    extern __shared__ float temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;
    int j;

    // k vectors v(i)
    if (i<n){
        for( j=0; j<k; j++)
            temp[Idx+j*blockDim.x] = v[i+j*n] * r[i];
    }
    else{
        for( j=0; j<k; j++)
            temp[Idx+j*blockDim.x] =MAGMA_tally3_S_MAKE( 0.0, 0.0);
    }
    __syncthreads();
    if ( Idx < 128 ){
        for( j=0; j<k; j++){
            temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 128 ];
        }
    }
    __syncthreads();
    if ( Idx < 64 ){
        for( j=0; j<k; j++){
            temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 64 ];
        }
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            for( j=0; j<k; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 32 ];
                __syncthreads();
            for( j=0; j<k; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 16 ];
                __syncthreads();
            for( j=0; j<k; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 8 ];
                __syncthreads();
            for( j=0; j<k; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 4 ];
                __syncthreads();
            for( j=0; j<k; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 2 ];
                __syncthreads();
            for( j=0; j<k; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 1 ];
                __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            for( j=0; j<k; j++){
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
            for( j=0; j<k; j++){
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
        for( j=0; j<k; j++){
            vtmp[ blockIdx.x+j*n ] = temp[ j*blockDim.x ];
        }
    }
}

// block reduction for multiple vectors
__global__ void
magma_tally3_sblockreduce_kernel( 
    int Gs,
    int n, 
    int k,
    float * vtmp,
    float * vtmp2 )
{

    extern __shared__ float temp[];    
    int Idx = threadIdx.x;
    int i = blockIdx.x * blockDim.x + Idx;  
    int j;
    for( j=0; j<k; j++){
        temp[ Idx+j*blockDim.x ] =  ( i < n ) ? vtmp[ i+j*n ] 
                                        : MAGMA_tally3_S_MAKE( 0.0, 0.0);
    }
    __syncthreads();
    if ( Idx < 128 ){
        for( j=0; j<k; j++){
            temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 128 ];
        }
    }
    __syncthreads();
    if ( Idx < 64 ){
        for( j=0; j<k; j++){
            temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 64 ];
        }
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            for( j=0; j<k; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 32 ];
                __syncthreads();
            for( j=0; j<k; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 16 ];
                __syncthreads();
            for( j=0; j<k; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 8 ];
                __syncthreads();
            for( j=0; j<k; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 4 ];
                __syncthreads();
            for( j=0; j<k; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 2 ];
                __syncthreads();
            for( j=0; j<k; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 1 ];
                __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            for( j=0; j<k; j++){
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
            for( j=0; j<k; j++){
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
        for( j=0; j<k; j++){
            vtmp2[ blockIdx.x+j*n ] = temp[ j*blockDim.x ];
        }
    }
}

// accelerated reduction for one vector
__global__ void
magma_tally3_sreduce_kernel_fast( int Gs,
                           int n, 
                           float * vtmp,
                           float * vtmp2 ){

    extern __shared__ float temp[];    
    int Idx = threadIdx.x;
    int blockSize = 128;
    int gridSize = blockSize  * 2 * gridDim.x; 
    temp[Idx] = MAGMA_tally3_S_MAKE( 0.0, 0.0);
    int i = blockIdx.x * ( blockSize * 2 ) + Idx;   
    while (i < Gs ) {
        temp[ Idx  ] += vtmp[ i ]; 
        temp[ Idx  ] += ( i + blockSize < Gs ) ? vtmp[ i + blockSize ] 
                                                : MAGMA_tally3_S_MAKE( 0.0, 0.0); 
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

// accelerated block reduction for multiple vectors
__global__ void
magma_tally3_sblockreduce_kernel_fast( 
    int Gs,
    int n, 
    int k,
    float * vtmp,
    float * vtmp2 )
{

    extern __shared__ float temp[];    
    int Idx = threadIdx.x;
    int blockSize = 128;
    int gridSize = blockSize  * 2 * gridDim.x; 
    int j;

    for( j=0; j<k; j++){
        int i = blockIdx.x * ( blockSize * 2 ) + Idx;   
        temp[Idx+j*(blockSize)] = MAGMA_tally3_S_MAKE( 0.0, 0.0);
        while (i < Gs ) {
            temp[ Idx+j*(blockSize)  ] += vtmp[ i+j*n ]; 
            temp[ Idx+j*(blockSize)  ] += 
                ( i + (blockSize) < Gs ) ? vtmp[ i+j*n + (blockSize) ] 
                                                : MAGMA_tally3_S_MAKE( 0.0, 0.0); 
            i += gridSize;
        }
    }
    __syncthreads();
    if ( Idx < 64 ){
        for( j=0; j<k; j++){
            temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 64 ];
        }
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            for( j=0; j<k; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 32 ];
                __syncthreads();
            for( j=0; j<k; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 16 ];
                __syncthreads();
            for( j=0; j<k; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 8 ];
                __syncthreads();
            for( j=0; j<k; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 4 ];
                __syncthreads();
            for( j=0; j<k; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 2 ];
                __syncthreads();
            for( j=0; j<k; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 1 ];
                __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            for( j=0; j<k; j++){
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
            for( j=0; j<k; j++){
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
        for( j=0; j<k; j++){
            vtmp2[ blockIdx.x+j*n ] = temp[ j*(blockSize) ];
        }
    }
}

/**
    Purpose
    -------

    Computes the scalar product of a set of vectors v_i such that

    skp = ( <v_0,r>, <v_1,r>, .. )

    Returns the vector skp.

    Arguments
    ---------

    @param[in]
    n           int
                length of v_i and r

    @param[in]
    k           int
                # vectors v_i

    @param[in]
    v           magma_tally3Float_ptr 
                v = (v_0 .. v_i.. v_k)

    @param[in]
    r           magma_tally3Float_ptr 
                r

    @param[in]
    d1          magma_tally3Float_ptr 
                workspace

    @param[in]
    d2          magma_tally3Float_ptr 
                workspace

    @param[out]
    skp         magma_tally3Float_ptr 
                vector[k] of scalar products (<v_i,r>...)

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_sblas
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_smdotc(
    int n, 
    int k, 
    magma_tally3Float_ptr v, 
    magma_tally3Float_ptr r,
    magma_tally3Float_ptr d1,
    magma_tally3Float_ptr d2,
    magma_tally3Float_ptr skp,
    magma_tally3_queue_t queue )
{
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue;
    magma_tally3blasGetKernelStream( &orig_queue );

    int local_block_size=256;
    dim3 Bs( local_block_size );
    dim3 Gs( magma_tally3_ceildiv( n, local_block_size ) );
    dim3 Gs_next;
    int Ms =  (k)* (local_block_size) * sizeof( float ); // k vecs 
    magma_tally3Float_ptr aux1 = d1, aux2 = d2;
    int b = 1;        

    if (k>1) {
        magma_tally3_sblockdot_kernel<<<Gs, Bs, Ms>>>( Gs.x, n, k, v, r, d1 );
    }
    else {
        magma_tally3_sdot_kernel<<<Gs, Bs, Ms>>>( Gs.x, n, v, r, d1 );
    }
/*
    // not necessary to zero GPU mem
    magma_tally3_sgpumemzero<<<Gs, Bs, 0>>>( d1, n*k,1 );
    magma_tally3_sgpumemzero<<<Gs, Bs, 0>>>( d2, n*k,1 );
    //magma_tally3blas_slaset( Magma_tally3UpperLower, n, k, d1, n );
    //magma_tally3blas_slaset( Magma_tally3UpperLower, n, k, d2, n );
    while( Gs.x > 1 ) {
        Gs_next.x = ( Gs.x+Bs.x-1 )/ Bs.x ;
        magma_tally3_sblockreduce_kernel<<< Gs_next.x, Bs.x, Ms >>> 
                                        ( Gs.x, n, k, aux1, aux2 );
        Gs.x = Gs_next.x;
        b = 1 - b;
        if ( b ) { aux1 = d1; aux2 = d2; }
        else   { aux2 = d1; aux1 = d2; }
    }
    for( int j=0; j<k; j++) {
            magma_tally3_scopyvector( 1, aux1+j*n, 1, skp+j, 1 );
    }
*/
   
    if ( k>1) {
        while( Gs.x > 1 ) {
            Gs_next.x = ( Gs.x+Bs.x-1 )/ Bs.x ;
            if ( Gs_next.x == 1 ) Gs_next.x = 2;
            magma_tally3_sblockreduce_kernel_fast<<< Gs_next.x/2, Bs.x/2, Ms/2 >>> 
                        ( Gs.x, n, k, aux1, aux2 );
            Gs_next.x = Gs_next.x /2;
            Gs.x = Gs_next.x;
            b = 1 - b;
            if ( b ) { aux1 = d1; aux2 = d2; }
            else   { aux2 = d1; aux1 = d2; }
        }
    }
    else {
        while( Gs.x > 1 ) {
            Gs_next.x = ( Gs.x+Bs.x-1 )/ Bs.x ;
            if ( Gs_next.x == 1 ) Gs_next.x = 2;
            magma_tally3_sreduce_kernel_fast<<< Gs_next.x/2, Bs.x/2, Ms/2 >>> 
                        ( Gs.x, n, aux1, aux2 );
            Gs_next.x = Gs_next.x /2;
            Gs.x = Gs_next.x;
            b = 1 - b;
            if ( b ) { aux1 = d1; aux2 = d2; }
            else   { aux2 = d1; aux1 = d2; }
        }
    }


    for( int j=0; j<k; j++) {
            magma_tally3_scopyvector( 1, aux1+j*n, 1, skp+j, 1 );
    }

   magma_tally3blasSetKernelStream( orig_queue );
   return MAGMA_tally3_SUCCESS;
}

/**
    Purpose
    -------

    This is an extension of the merged dot product above by chunking
    the set of vectors v_i such that the data always fits into cache.
    It is equivalent to a matrix vecor product Vr where V
    contains few rows and many columns. The computation is the same:

    skp = ( <v_0,r>, <v_1,r>, .. )

    Returns the vector skp.

    Arguments
    ---------

    @param[in]
    n           int
                length of v_i and r

    @param[in]
    k           int
                # vectors v_i

    @param[in]
    v           magma_tally3Float_ptr 
                v = (v_0 .. v_i.. v_k)

    @param[in]
    r           magma_tally3Float_ptr 
                r

    @param[in]
    d1          magma_tally3Float_ptr 
                workspace

    @param[in]
    d2          magma_tally3Float_ptr 
                workspace

    @param[out]
    skp         magma_tally3Float_ptr 
                vector[k] of scalar products (<v_i,r>...)

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_s
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_sgemvmdot(
    int n, 
    int k, 
    magma_tally3Float_ptr v, 
    magma_tally3Float_ptr r,
    magma_tally3Float_ptr d1,
    magma_tally3Float_ptr d2,
    magma_tally3Float_ptr skp,
    magma_tally3_queue_t queue )
{
    int rows_left = k;
    int offset = 0;
    int chunk_size = 4;
    // process in chunks of 10 - has to be adapted to hardware and precision
    while( rows_left > (chunk_size) ) {
        magma_tally3_smdotc( n, chunk_size, v+offset*n, r, d1, d2, skp+offset, queue );
        offset = offset + chunk_size;
        rows_left = rows_left-chunk_size;

    }
    // process rest
    magma_tally3_smdotc( n, rows_left, v+offset*n, r, d1, d2, skp+offset, queue ); 


   return MAGMA_tally3_SUCCESS;
}


