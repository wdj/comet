/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgerbt.cu normal z -> d, Fri Jan 30 19:00:08 2015


       @author Adrien REMY
*/
#include "common_magma_tally3.h"
#include "dgerbt.h"


#define block_height  32
#define block_width  4
#define block_length 256
#define NB 64
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
    Purpose
    -------
    DPRBT_MVT compute B = UTB to randomize B
    
    Arguments
    ---------
    @param[in]
    n       INTEGER
            The number of values of db.  n >= 0.

    @param[in]
    du     DOUBLE_PRECISION array, dimension (n,2)
            The 2*n vector representing the random butterfly matrix V
    
    @param[in,out]
    db     DOUBLE_PRECISION array, dimension (n)
            The n vector db computed by DGESV_NOPIV_GPU
            On exit db = du*db
    
    @param[in]
    queue   magma_tally3_queue_t
            Queue to execute in.
    ********************************************************************/
extern "C" void
magma_tally3blas_dprbt_mtv_q(
    magma_tally3_int_t n, 
    double *du, double *db,
    magma_tally3_queue_t queue)
{
    /*

     */
    magma_tally3_int_t threads = block_length;
    magma_tally3_int_t grid = n/(4*block_length) + ((n%(4*block_length))!=0);

    magma_tally3blas_dapply_transpose_vector_kernel<<< grid, threads, 0, queue >>>(n/2, du, n, db, 0);
    magma_tally3blas_dapply_transpose_vector_kernel<<< grid, threads, 0, queue >>>(n/2, du, n+n/2, db, n/2);

    threads = block_length;
    grid = n/(2*block_length) + ((n%(2*block_length))!=0);
    magma_tally3blas_dapply_transpose_vector_kernel<<< grid, threads, 0, queue >>>(n, du, 0, db, 0);
}

/**
    @see magma_tally3blas_dprbt_mtv_q
    ********************************************************************/
extern "C" void
magma_tally3blas_dprbt_mtv(
    magma_tally3_int_t n, 
    double *du, double *db)
{
    magma_tally3blas_dprbt_mtv_q(n, du, db, magma_tally3_stream);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
    Purpose
    -------
    DPRBT_MV compute B = VB to obtain the non randomized solution
    
    Arguments
    ---------
    @param[in]
    n       INTEGER
            The number of values of db.  n >= 0.
    
    @param[in,out]
    db      DOUBLE_PRECISION array, dimension (n)
            The n vector db computed by DGESV_NOPIV_GPU
            On exit db = dv*db
    
    @param[in]
    dv      DOUBLE_PRECISION array, dimension (n,2)
            The 2*n vector representing the random butterfly matrix V
    
    @param[in]
    queue   magma_tally3_queue_t
            Queue to execute in.
    ********************************************************************/
extern "C" void
magma_tally3blas_dprbt_mv_q(
    magma_tally3_int_t n, 
    double *dv, double *db,
    magma_tally3_queue_t queue)
{

    magma_tally3_int_t threads = block_length;
    magma_tally3_int_t grid = n/(2*block_length) + ((n%(2*block_length))!=0);

    magma_tally3blas_dapply_vector_kernel<<< grid, threads, 0, queue >>>(n, dv, 0, db, 0);


    threads = block_length;
    grid = n/(4*block_length) + ((n%(4*block_length))!=0);

    magma_tally3blas_dapply_vector_kernel<<< grid, threads, 0, queue >>>(n/2, dv, n, db, 0);
    magma_tally3blas_dapply_vector_kernel<<< grid, threads, 0, queue >>>(n/2, dv, n+n/2, db, n/2);
}

/**
    @see magma_tally3blas_dprbt_mtv_q
    ********************************************************************/
extern "C" void
magma_tally3blas_dprbt_mv(
    magma_tally3_int_t n, 
    double *dv, double *db)
{
    magma_tally3blas_dprbt_mv_q(n, dv, db, magma_tally3_stream);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
    Purpose
    -------
    DPRBT randomize a square general matrix using partial randomized transformation
    
    Arguments
    ---------
    @param[in]
    n       INTEGER
            The number of columns and rows of the matrix dA.  n >= 0.
    
    @param[in,out]
    dA      DOUBLE_PRECISION array, dimension (n,ldda)
            The n-by-n matrix dA
            On exit dA = duT*dA*d_V
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDA >= max(1,n).
    
    @param[in]
    du      DOUBLE_PRECISION array, dimension (n,2)
            The 2*n vector representing the random butterfly matrix U
    
    @param[in]
    dv      DOUBLE_PRECISION array, dimension (n,2)
            The 2*n vector representing the random butterfly matrix V
    
    @param[in]
    queue   magma_tally3_queue_t
            Queue to execute in.

    ********************************************************************/
extern "C" void 
magma_tally3blas_dprbt_q(
    magma_tally3_int_t n, 
    double *dA, magma_tally3_int_t ldda, 
    double *du, double *dv,
    magma_tally3_queue_t queue)
{
    du += ldda;
    dv += ldda;

    dim3 threads(block_height, block_width);
    dim3 grid(n/(4*block_height) + ((n%(4*block_height))!=0), 
            n/(4*block_width)  + ((n%(4*block_width))!=0));

    magma_tally3blas_delementary_multiplication_kernel<<< grid, threads, 0, queue >>>(n/2, dA,            0, ldda, du,   0, dv,   0);
    magma_tally3blas_delementary_multiplication_kernel<<< grid, threads, 0, queue >>>(n/2, dA,     ldda*n/2, ldda, du,   0, dv, n/2);
    magma_tally3blas_delementary_multiplication_kernel<<< grid, threads, 0, queue >>>(n/2, dA,          n/2, ldda, du, n/2, dv,   0);
    magma_tally3blas_delementary_multiplication_kernel<<< grid, threads, 0, queue >>>(n/2, dA, ldda*n/2+n/2, ldda, du, n/2, dv, n/2);

    dim3 threads2(block_height, block_width);
    dim3 grid2(n/(2*block_height) + ((n%(2*block_height))!=0), 
            n/(2*block_width)  + ((n%(2*block_width))!=0));
    magma_tally3blas_delementary_multiplication_kernel<<< grid2, threads2, 0, queue >>>(n, dA, 0, ldda, du, -ldda, dv, -ldda);
}


/**
    @see magma_tally3blas_dprbt_q
    ********************************************************************/
extern "C" void 
magma_tally3blas_dprbt(
    magma_tally3_int_t n, 
    double *dA, magma_tally3_int_t ldda, 
    double *du, double *dv)
{
    magma_tally3blas_dprbt_q(n, dA, ldda, du, dv, magma_tally3_stream);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// adds   x += r  --and--
// copies r = b
// each thread does one index, x[i] and r[i]
__global__ void
daxpycp2_kernel(
    int m, double *r, double *x,
    const double *b)
{
    const int i = threadIdx.x + blockIdx.x*NB;
    if ( i < m ) {
        x[i] = MAGMA_tally3_D_ADD( x[i], r[i] );
        r[i] = b[i];
    }
}


// ----------------------------------------------------------------------
// adds   x += r  --and--
// copies r = b
extern "C" void
magma_tally3blas_daxpycp2_q(
    magma_tally3_int_t m, double *r, double *x,
    const double *b,
    magma_tally3_queue_t queue )
{
    dim3 threads( NB );
    dim3 grid( (m + NB - 1)/NB );
    daxpycp2_kernel <<< grid, threads, 0, queue >>> ( m, r, x, b );
}


extern "C" void
magma_tally3blas_daxpycp2(
    magma_tally3_int_t m, double *r, double *x,
    const double *b)
{
    magma_tally3blas_daxpycp2_q( m, r, x, b, magma_tally3_stream );
}
