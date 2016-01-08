/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c


       @author Adrien REMY
*/
#include "common_magma_tally4.h"
#include "zgerbt.h"


#define block_height  32
#define block_width  4
#define block_length 256
#define NB 64
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
    Purpose
    -------
    ZPRBT_MVT compute B = UTB to randomize B
    
    Arguments
    ---------
    @param[in]
    n       INTEGER
            The number of values of db.  n >= 0.

    @param[in]
    du     COMPLEX_16 array, dimension (n,2)
            The 2*n vector representing the random butterfly matrix V
    
    @param[in,out]
    db     COMPLEX_16 array, dimension (n)
            The n vector db computed by ZGESV_NOPIV_GPU
            On exit db = du*db
    
    @param[in]
    queue   magma_tally4_queue_t
            Queue to execute in.
    ********************************************************************/
extern "C" void
magma_tally4blas_zprbt_mtv_q(
    magma_tally4_int_t n, 
    magma_tally4DoubleComplex *du, magma_tally4DoubleComplex *db,
    magma_tally4_queue_t queue)
{
    /*

     */
    magma_tally4_int_t threads = block_length;
    magma_tally4_int_t grid = n/(4*block_length) + ((n%(4*block_length))!=0);

    magma_tally4blas_zapply_transpose_vector_kernel<<< grid, threads, 0, queue >>>(n/2, du, n, db, 0);
    magma_tally4blas_zapply_transpose_vector_kernel<<< grid, threads, 0, queue >>>(n/2, du, n+n/2, db, n/2);

    threads = block_length;
    grid = n/(2*block_length) + ((n%(2*block_length))!=0);
    magma_tally4blas_zapply_transpose_vector_kernel<<< grid, threads, 0, queue >>>(n, du, 0, db, 0);
}

/**
    @see magma_tally4blas_zprbt_mtv_q
    ********************************************************************/
extern "C" void
magma_tally4blas_zprbt_mtv(
    magma_tally4_int_t n, 
    magma_tally4DoubleComplex *du, magma_tally4DoubleComplex *db)
{
    magma_tally4blas_zprbt_mtv_q(n, du, db, magma_tally4_stream);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
    Purpose
    -------
    ZPRBT_MV compute B = VB to obtain the non randomized solution
    
    Arguments
    ---------
    @param[in]
    n       INTEGER
            The number of values of db.  n >= 0.
    
    @param[in,out]
    db      COMPLEX_16 array, dimension (n)
            The n vector db computed by ZGESV_NOPIV_GPU
            On exit db = dv*db
    
    @param[in]
    dv      COMPLEX_16 array, dimension (n,2)
            The 2*n vector representing the random butterfly matrix V
    
    @param[in]
    queue   magma_tally4_queue_t
            Queue to execute in.
    ********************************************************************/
extern "C" void
magma_tally4blas_zprbt_mv_q(
    magma_tally4_int_t n, 
    magma_tally4DoubleComplex *dv, magma_tally4DoubleComplex *db,
    magma_tally4_queue_t queue)
{

    magma_tally4_int_t threads = block_length;
    magma_tally4_int_t grid = n/(2*block_length) + ((n%(2*block_length))!=0);

    magma_tally4blas_zapply_vector_kernel<<< grid, threads, 0, queue >>>(n, dv, 0, db, 0);


    threads = block_length;
    grid = n/(4*block_length) + ((n%(4*block_length))!=0);

    magma_tally4blas_zapply_vector_kernel<<< grid, threads, 0, queue >>>(n/2, dv, n, db, 0);
    magma_tally4blas_zapply_vector_kernel<<< grid, threads, 0, queue >>>(n/2, dv, n+n/2, db, n/2);
}

/**
    @see magma_tally4blas_zprbt_mtv_q
    ********************************************************************/
extern "C" void
magma_tally4blas_zprbt_mv(
    magma_tally4_int_t n, 
    magma_tally4DoubleComplex *dv, magma_tally4DoubleComplex *db)
{
    magma_tally4blas_zprbt_mv_q(n, dv, db, magma_tally4_stream);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
    Purpose
    -------
    ZPRBT randomize a square general matrix using partial randomized transformation
    
    Arguments
    ---------
    @param[in]
    n       INTEGER
            The number of columns and rows of the matrix dA.  n >= 0.
    
    @param[in,out]
    dA      COMPLEX_16 array, dimension (n,ldda)
            The n-by-n matrix dA
            On exit dA = duT*dA*d_V
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDA >= max(1,n).
    
    @param[in]
    du      COMPLEX_16 array, dimension (n,2)
            The 2*n vector representing the random butterfly matrix U
    
    @param[in]
    dv      COMPLEX_16 array, dimension (n,2)
            The 2*n vector representing the random butterfly matrix V
    
    @param[in]
    queue   magma_tally4_queue_t
            Queue to execute in.

    ********************************************************************/
extern "C" void 
magma_tally4blas_zprbt_q(
    magma_tally4_int_t n, 
    magma_tally4DoubleComplex *dA, magma_tally4_int_t ldda, 
    magma_tally4DoubleComplex *du, magma_tally4DoubleComplex *dv,
    magma_tally4_queue_t queue)
{
    du += ldda;
    dv += ldda;

    dim3 threads(block_height, block_width);
    dim3 grid(n/(4*block_height) + ((n%(4*block_height))!=0), 
            n/(4*block_width)  + ((n%(4*block_width))!=0));

    magma_tally4blas_zelementary_multiplication_kernel<<< grid, threads, 0, queue >>>(n/2, dA,            0, ldda, du,   0, dv,   0);
    magma_tally4blas_zelementary_multiplication_kernel<<< grid, threads, 0, queue >>>(n/2, dA,     ldda*n/2, ldda, du,   0, dv, n/2);
    magma_tally4blas_zelementary_multiplication_kernel<<< grid, threads, 0, queue >>>(n/2, dA,          n/2, ldda, du, n/2, dv,   0);
    magma_tally4blas_zelementary_multiplication_kernel<<< grid, threads, 0, queue >>>(n/2, dA, ldda*n/2+n/2, ldda, du, n/2, dv, n/2);

    dim3 threads2(block_height, block_width);
    dim3 grid2(n/(2*block_height) + ((n%(2*block_height))!=0), 
            n/(2*block_width)  + ((n%(2*block_width))!=0));
    magma_tally4blas_zelementary_multiplication_kernel<<< grid2, threads2, 0, queue >>>(n, dA, 0, ldda, du, -ldda, dv, -ldda);
}


/**
    @see magma_tally4blas_zprbt_q
    ********************************************************************/
extern "C" void 
magma_tally4blas_zprbt(
    magma_tally4_int_t n, 
    magma_tally4DoubleComplex *dA, magma_tally4_int_t ldda, 
    magma_tally4DoubleComplex *du, magma_tally4DoubleComplex *dv)
{
    magma_tally4blas_zprbt_q(n, dA, ldda, du, dv, magma_tally4_stream);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// adds   x += r  --and--
// copies r = b
// each thread does one index, x[i] and r[i]
__global__ void
zaxpycp2_kernel(
    int m, magma_tally4DoubleComplex *r, magma_tally4DoubleComplex *x,
    const magma_tally4DoubleComplex *b)
{
    const int i = threadIdx.x + blockIdx.x*NB;
    if ( i < m ) {
        x[i] = MAGMA_tally4_Z_ADD( x[i], r[i] );
        r[i] = b[i];
    }
}


// ----------------------------------------------------------------------
// adds   x += r  --and--
// copies r = b
extern "C" void
magma_tally4blas_zaxpycp2_q(
    magma_tally4_int_t m, magma_tally4DoubleComplex *r, magma_tally4DoubleComplex *x,
    const magma_tally4DoubleComplex *b,
    magma_tally4_queue_t queue )
{
    dim3 threads( NB );
    dim3 grid( (m + NB - 1)/NB );
    zaxpycp2_kernel <<< grid, threads, 0, queue >>> ( m, r, x, b );
}


extern "C" void
magma_tally4blas_zaxpycp2(
    magma_tally4_int_t m, magma_tally4DoubleComplex *r, magma_tally4DoubleComplex *x,
    const magma_tally4DoubleComplex *b)
{
    magma_tally4blas_zaxpycp2_q( m, r, x, b, magma_tally4_stream );
}
