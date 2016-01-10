/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s
       @author Hartwig Anzt

*/
#include "common_magma_tally4.h"

#define BLOCK_SIZE 512

#define PRECISION_z


// These routines merge multiple kernels from zmergebicgstab into one
// The difference to zmergedbicgstab2 is that the SpMV is not merged into the
// kernes. This results in higher flexibility at the price of lower performance.

/* -------------------------------------------------------------------------- */

__global__ void
magma_tally4_zbicgmerge1_kernel(  
    int n, 
    magma_tally4DoubleComplex * skp,
    magma_tally4DoubleComplex * v, 
    magma_tally4DoubleComplex * r, 
    magma_tally4DoubleComplex * p )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    magma_tally4DoubleComplex beta=skp[1];
    magma_tally4DoubleComplex omega=skp[2];
    if( i<n ){
        p[i] =  r[i] + beta * ( p[i] - omega * v[i] );

    }

}

/**
    Purpose
    -------

    Mergels multiple operations into one kernel:

    p = beta*p
    p = p-omega*beta*v
    p = p+r
    
    -> p = r + beta * ( p - omega * v ) 

    Arguments
    ---------

    @param[in]
    n           int
                dimension n

    @param[in]
    skp         magma_tally4DoubleComplex_ptr 
                set of scalar parameters

    @param[in]
    v           magma_tally4DoubleComplex_ptr 
                input v

    @param[in]
    r           magma_tally4DoubleComplex_ptr 
                input r

    @param[in/out]
    p           magma_tally4DoubleComplex_ptr 
                input/output p

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zgegpuk
    ********************************************************************/

extern "C" int
magma_tally4_zbicgmerge1(  
    int n, 
    magma_tally4DoubleComplex_ptr skp,
    magma_tally4DoubleComplex_ptr v, 
    magma_tally4DoubleComplex_ptr r, 
    magma_tally4DoubleComplex_ptr p ){

    
    dim3 Bs( BLOCK_SIZE );
    dim3 Gs( magma_tally4_ceildiv( n, BLOCK_SIZE ) );
    magma_tally4_zbicgmerge1_kernel<<<Gs, Bs, 0>>>( n, skp, v, r, p );

   return MAGMA_tally4_SUCCESS;
}

/* -------------------------------------------------------------------------- */

__global__ void
magma_tally4_zbicgmerge2_kernel(  
    int n, 
    magma_tally4DoubleComplex * skp, 
    magma_tally4DoubleComplex * r,
    magma_tally4DoubleComplex * v, 
    magma_tally4DoubleComplex * s )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    magma_tally4DoubleComplex alpha=skp[0];
    if( i<n ){
        s[i] =  r[i] - alpha * v[i] ;
    }

}

/**
    Purpose
    -------

    Mergels multiple operations into one kernel:

    s=r
    s=s-alpha*v
        
    -> s = r - alpha * v

    Arguments
    ---------

    @param[in]
    n           int
                dimension n

    @param[in]
    skp         magma_tally4DoubleComplex_ptr 
                set of scalar parameters

    @param[in]
    r           magma_tally4DoubleComplex_ptr 
                input r

    @param[in]
    v           magma_tally4DoubleComplex_ptr 
                input v

    @param[s]
    s           magma_tally4DoubleComplex_ptr 
                output s

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zgegpuk
    ********************************************************************/

extern "C" int
magma_tally4_zbicgmerge2(  
    int n, 
    magma_tally4DoubleComplex_ptr skp, 
    magma_tally4DoubleComplex_ptr r,
    magma_tally4DoubleComplex_ptr v, 
    magma_tally4DoubleComplex_ptr s )
{

    
    dim3 Bs( BLOCK_SIZE );
    dim3 Gs( magma_tally4_ceildiv( n, BLOCK_SIZE ) );

    magma_tally4_zbicgmerge2_kernel<<<Gs, Bs, 0>>>( n, skp, r, v, s );

   return MAGMA_tally4_SUCCESS;
}

/* -------------------------------------------------------------------------- */

__global__ void
magma_tally4_zbicgmerge3_kernel(  
    int n, 
    magma_tally4DoubleComplex * skp, 
    magma_tally4DoubleComplex * p,
    magma_tally4DoubleComplex * se,
    magma_tally4DoubleComplex * t,
    magma_tally4DoubleComplex * x, 
    magma_tally4DoubleComplex * r )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    magma_tally4DoubleComplex alpha=skp[0];
    magma_tally4DoubleComplex omega=skp[2];
    if( i<n ){
        magma_tally4DoubleComplex s;
        s = se[i];
        x[i] = x[i] + alpha * p[i] + omega * s;
        r[i] = s - omega * t[i];
    }

}

/**
    Purpose
    -------

    Mergels multiple operations into one kernel:

    x=x+alpha*p
    x=x+omega*s
    r=s
    r=r-omega*t
        
    -> x = x + alpha * p + omega * s
    -> r = s - omega * t

    Arguments
    ---------

    @param[in]
    n           int
                dimension n

    @param[in]
    skp         magma_tally4DoubleComplex_ptr 
                set of scalar parameters

    @param[in]
    p           magma_tally4DoubleComplex_ptr 
                input p

    @param[in]
    s           magma_tally4DoubleComplex_ptr 
                input s

    @param[in]
    t           magma_tally4DoubleComplex_ptr 
                input t

    @param[in/out]
    x           magma_tally4DoubleComplex_ptr 
                input/output x

    @param[in/out]
    r           magma_tally4DoubleComplex_ptr 
                input/output r

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zgegpuk
    ********************************************************************/

extern "C" int
magma_tally4_zbicgmerge3(  
    int n, 
    magma_tally4DoubleComplex_ptr skp,
    magma_tally4DoubleComplex_ptr p,
    magma_tally4DoubleComplex_ptr s,
    magma_tally4DoubleComplex_ptr t,
    magma_tally4DoubleComplex_ptr x, 
    magma_tally4DoubleComplex_ptr r )
{

    
    dim3 Bs( BLOCK_SIZE );
    dim3 Gs( magma_tally4_ceildiv( n, BLOCK_SIZE ) );
    magma_tally4_zbicgmerge3_kernel<<<Gs, Bs, 0>>>( n, skp, p, s, t, x, r );

   return MAGMA_tally4_SUCCESS;
}

/* -------------------------------------------------------------------------- */

__global__ void
magma_tally4_zbicgmerge4_kernel_1(  
    magma_tally4DoubleComplex * skp )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if( i==0 ){
        magma_tally4DoubleComplex tmp = skp[0];
        skp[0] = skp[4]/tmp;
    }
}

__global__ void
magma_tally4_zbicgmerge4_kernel_2(  
    magma_tally4DoubleComplex * skp )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if( i==0 ){
        skp[2] = skp[6]/skp[7];
        skp[3] = skp[4];
    }
}

__global__ void
magma_tally4_zbicgmerge4_kernel_3(  
    magma_tally4DoubleComplex * skp )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if( i==0 ){
        magma_tally4DoubleComplex tmp1 = skp[4]/skp[3];
        magma_tally4DoubleComplex tmp2 = skp[0] / skp[2];
        skp[1] =  tmp1*tmp2;
        //skp[1] =  skp[4]/skp[3] * skp[0] / skp[2];

    }
}

/**
    Purpose
    -------

    Performs some parameter operations for the BiCGSTAB with scalars on GPU.

    Arguments
    ---------

    @param[in]
    type        int
                kernel type

    @param[in/out]
    skp         magma_tally4DoubleComplex_ptr 
                vector with parameters

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zgegpuk
    ********************************************************************/

extern "C" int
magma_tally4_zbicgmerge4(  
    int type, 
    magma_tally4DoubleComplex_ptr skp )
{

    dim3 Bs( 1 );
    dim3 Gs( 1 );
    if( type == 1 )
        magma_tally4_zbicgmerge4_kernel_1<<<Gs, Bs, 0>>>( skp );
    else if( type == 2 )
        magma_tally4_zbicgmerge4_kernel_2<<<Gs, Bs, 0>>>( skp );
    else if( type == 3 )
        magma_tally4_zbicgmerge4_kernel_3<<<Gs, Bs, 0>>>( skp );
    else
        printf("error: no kernel called\n");

   return MAGMA_tally4_SUCCESS;
}
