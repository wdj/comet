/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zjacobisetup.cu normal z -> c, Sun May  3 11:22:58 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally4sparse.h"

#define BLOCK_SIZE 128


#define PRECISION_c

__global__ void 
cvjacobisetup_gpu(  int num_rows, 
                    int num_vecs,
                    magma_tally4FloatComplex *b, 
                    magma_tally4FloatComplex *d, 
                    magma_tally4FloatComplex *c,
                    magma_tally4FloatComplex *x){

    int row = blockDim.x * blockIdx.x + threadIdx.x ;

    if(row < num_rows ){
        for( int i=0; i<num_vecs; i++ ){
            c[row+i*num_rows] = b[row+i*num_rows] / d[row];
            x[row+i*num_rows] = c[row+i*num_rows];
        }
    }
}





/**
    Purpose
    -------

    Prepares the Jacobi Iteration according to
       x^(k+1) = D^(-1) * b - D^(-1) * (L+U) * x^k
       x^(k+1) =      c     -       M        * x^k.

    Returns the vector c. It calls a GPU kernel

    Arguments
    ---------

    @param[in]
    num_rows    magma_tally4_int_t
                number of rows
                
    @param[in]
    b           magma_tally4_c_matrix
                RHS b

    @param[in]
    d           magma_tally4_c_matrix
                vector with diagonal entries

    @param[out]
    c           magma_tally4_c_matrix*
                c = D^(-1) * b

    @param[out]
    x           magma_tally4_c_matrix*
                iteration vector
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_cgegpuk
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_cjacobisetup_vector_gpu(
    int num_rows, 
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix d, 
    magma_tally4_c_matrix c,
    magma_tally4_c_matrix *x,
    magma_tally4_queue_t queue )
{
    dim3 grid( magma_tally4_ceildiv( num_rows, BLOCK_SIZE ) );
    int num_vecs = b.num_rows / num_rows;
    magma_tally4_int_t threads = BLOCK_SIZE;
    cvjacobisetup_gpu<<< grid, threads, 0 >>>
                ( num_rows, num_vecs, b.dval, d.dval, c.dval, x->val );

    return MAGMA_tally4_SUCCESS;
}






__global__ void 
cjacobidiagscal_kernel(  int num_rows,
                         int num_vecs, 
                    magma_tally4FloatComplex *b, 
                    magma_tally4FloatComplex *d, 
                    magma_tally4FloatComplex *c){

    int row = blockDim.x * blockIdx.x + threadIdx.x ;

    if(row < num_rows ){
        for( int i=0; i<num_vecs; i++)
            c[row+i*num_rows] = b[row+i*num_rows] * d[row];
    }
}





/**
    Purpose
    -------

    Prepares the Jacobi Iteration according to
       x^(k+1) = D^(-1) * b - D^(-1) * (L+U) * x^k
       x^(k+1) =      c     -       M        * x^k.

    Returns the vector c. It calls a GPU kernel

    Arguments
    ---------

    @param[in]
    num_rows    magma_tally4_int_t
                number of rows
                
    @param[in]
    b           magma_tally4_c_matrix
                RHS b

    @param[in]
    d           magma_tally4_c_matrix
                vector with diagonal entries

    @param[out]
    c           magma_tally4_c_matrix*
                c = D^(-1) * b
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_c
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_cjacobi_diagscal(
    int num_rows, 
    magma_tally4_c_matrix d, 
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *c,
    magma_tally4_queue_t queue )
{
    dim3 grid( magma_tally4_ceildiv( num_rows, BLOCK_SIZE ));
    int num_vecs = b.num_rows*b.num_cols/num_rows;
    magma_tally4_int_t threads = BLOCK_SIZE;
    cjacobidiagscal_kernel<<< grid, threads, 0 >>>( num_rows, num_vecs, b.dval, d.dval, c->val );

    return MAGMA_tally4_SUCCESS;
}













__global__ void 
cjacobiupdate_kernel(  int num_rows,
                       int num_cols, 
                    magma_tally4FloatComplex *t, 
                    magma_tally4FloatComplex *b, 
                    magma_tally4FloatComplex *d, 
                    magma_tally4FloatComplex *x){

    int row = blockDim.x * blockIdx.x + threadIdx.x ;

    if(row < num_rows ){
        for( int i=0; i<num_cols; i++)
            x[row+i*num_rows] += (b[row+i*num_rows]-t[row+i*num_rows]) * d[row];
    }
}


/**
    Purpose
    -------

    Updates the iteration vector x for the Jacobi iteration
    according to
        x=x+d.*(b-t)
    where d is the diagonal of the system matrix A and t=Ax.

    Arguments
    ---------

    @param[in]
    num_rows    magma_tally4_int_t
                number of rows
                
    @param[in]
    num_cols    magma_tally4_int_t
                number of cols
                
    @param[in]
    t           magma_tally4_c_matrix
                t = A*x
                
    @param[in]
    b           magma_tally4_c_matrix
                RHS b
                
    @param[in]
    d           magma_tally4_c_matrix
                vector with diagonal entries

    @param[out]
    x           magma_tally4_c_matrix*
                iteration vector
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_c
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_cjacobiupdate(
    magma_tally4_c_matrix t, 
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix d, 
    magma_tally4_c_matrix *x,
    magma_tally4_queue_t queue )
{

    dim3 grid( magma_tally4_ceildiv( t.num_rows, BLOCK_SIZE ));
    magma_tally4_int_t threads = BLOCK_SIZE;
    cjacobiupdate_kernel<<< grid, threads, 0 >>>( t.num_rows, t.num_cols, t.dval, b.dval, d.dval, x->dval );

    return MAGMA_tally4_SUCCESS;
}










__global__ void 
cjacobispmvupdate_kernel(  
    int num_rows,
    int num_cols, 
    magma_tally4FloatComplex * dval, 
    magma_tally4_index_t * drowptr, 
    magma_tally4_index_t * dcolind,
    magma_tally4FloatComplex *t, 
    magma_tally4FloatComplex *b, 
    magma_tally4FloatComplex *d, 
    magma_tally4FloatComplex *x ){



    int row = blockDim.x * blockIdx.x + threadIdx.x ;
    int j;

    if(row<num_rows){
        magma_tally4FloatComplex dot = MAGMA_tally4_C_ZERO;
        int start = drowptr[ row ];
        int end = drowptr[ row+1 ];
        for( int i=0; i<num_cols; i++){
            for( j=start; j<end; j++){
                dot += dval[ j ] * x[ dcolind[j]+i*num_rows ];
            }
            x[row+i*num_rows] += (b[row+i*num_rows]-dot) * d[row];
        }
    }
}





/**
    Purpose
    -------

    Updates the iteration vector x for the Jacobi iteration
    according to
        x=x+d.*(b-Ax)


    Arguments
    ---------

    @param[in]
    maxiter     magma_tally4_int_t
                number of Jacobi iterations   
                
    @param[in]
    A           magma_tally4_c_matrix
                system matrix
                
    @param[in]
    t           magma_tally4_c_matrix
                workspace
                
    @param[in]
    b           magma_tally4_c_matrix
                RHS b
                
    @param[in]
    d           magma_tally4_c_matrix
                vector with diagonal entries

    @param[out]
    x           magma_tally4_c_matrix*
                iteration vector
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_c
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_cjacobispmvupdate(
    magma_tally4_int_t maxiter,
    magma_tally4_c_matrix A,
    magma_tally4_c_matrix t, 
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix d, 
    magma_tally4_c_matrix *x,
    magma_tally4_queue_t queue )
{

    // local variables
    magma_tally4FloatComplex c_zero = MAGMA_tally4_C_ZERO, c_one = MAGMA_tally4_C_ONE;
    dim3 grid( magma_tally4_ceildiv( t.num_rows, BLOCK_SIZE ));
    magma_tally4_int_t threads = BLOCK_SIZE;

    for( magma_tally4_int_t i=0; i<maxiter; i++ ) {
        // distinct routines imply synchronization
        // magma_tally4_c_spmv( c_one, A, *x, c_zero, t, queue );                // t =  A * x
        // cjacobiupdate_kernel<<< grid, threads, 0 >>>( t.num_rows, t.num_cols, t.dval, b.dval, d.dval, x->dval );
        // merged in one implies asynchronous update
        cjacobispmvupdate_kernel<<< grid, threads, 0 >>>
            ( t.num_rows, t.num_cols, A.dval, A.drow, A.dcol, t.dval, b.dval, d.dval, x->dval );

    }

    return MAGMA_tally4_SUCCESS;
}









