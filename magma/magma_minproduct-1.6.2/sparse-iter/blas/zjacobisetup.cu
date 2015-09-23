/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s
       @author Hartwig Anzt

*/
#include "common_magma_minproductsparse.h"

#define BLOCK_SIZE 128


#define PRECISION_z

__global__ void 
zvjacobisetup_gpu(  int num_rows, 
                    int num_vecs,
                    magma_minproductDoubleComplex *b, 
                    magma_minproductDoubleComplex *d, 
                    magma_minproductDoubleComplex *c,
                    magma_minproductDoubleComplex *x){

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
    num_rows    magma_minproduct_int_t
                number of rows
                
    @param[in]
    b           magma_minproduct_z_matrix
                RHS b

    @param[in]
    d           magma_minproduct_z_matrix
                vector with diagonal entries

    @param[out]
    c           magma_minproduct_z_matrix*
                c = D^(-1) * b

    @param[out]
    x           magma_minproduct_z_matrix*
                iteration vector
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zgegpuk
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zjacobisetup_vector_gpu(
    int num_rows, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix d, 
    magma_minproduct_z_matrix c,
    magma_minproduct_z_matrix *x,
    magma_minproduct_queue_t queue )
{
    dim3 grid( magma_minproduct_ceildiv( num_rows, BLOCK_SIZE ) );
    int num_vecs = b.num_rows / num_rows;
    magma_minproduct_int_t threads = BLOCK_SIZE;
    zvjacobisetup_gpu<<< grid, threads, 0 >>>
                ( num_rows, num_vecs, b.dval, d.dval, c.dval, x->val );

    return MAGMA_minproduct_SUCCESS;
}






__global__ void 
zjacobidiagscal_kernel(  int num_rows,
                         int num_vecs, 
                    magma_minproductDoubleComplex *b, 
                    magma_minproductDoubleComplex *d, 
                    magma_minproductDoubleComplex *c){

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
    num_rows    magma_minproduct_int_t
                number of rows
                
    @param[in]
    b           magma_minproduct_z_matrix
                RHS b

    @param[in]
    d           magma_minproduct_z_matrix
                vector with diagonal entries

    @param[out]
    c           magma_minproduct_z_matrix*
                c = D^(-1) * b
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_z
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zjacobi_diagscal(
    int num_rows, 
    magma_minproduct_z_matrix d, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *c,
    magma_minproduct_queue_t queue )
{
    dim3 grid( magma_minproduct_ceildiv( num_rows, BLOCK_SIZE ));
    int num_vecs = b.num_rows*b.num_cols/num_rows;
    magma_minproduct_int_t threads = BLOCK_SIZE;
    zjacobidiagscal_kernel<<< grid, threads, 0 >>>( num_rows, num_vecs, b.dval, d.dval, c->val );

    return MAGMA_minproduct_SUCCESS;
}













__global__ void 
zjacobiupdate_kernel(  int num_rows,
                       int num_cols, 
                    magma_minproductDoubleComplex *t, 
                    magma_minproductDoubleComplex *b, 
                    magma_minproductDoubleComplex *d, 
                    magma_minproductDoubleComplex *x){

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
    num_rows    magma_minproduct_int_t
                number of rows
                
    @param[in]
    num_cols    magma_minproduct_int_t
                number of cols
                
    @param[in]
    t           magma_minproduct_z_matrix
                t = A*x
                
    @param[in]
    b           magma_minproduct_z_matrix
                RHS b
                
    @param[in]
    d           magma_minproduct_z_matrix
                vector with diagonal entries

    @param[out]
    x           magma_minproduct_z_matrix*
                iteration vector
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_z
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zjacobiupdate(
    magma_minproduct_z_matrix t, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix d, 
    magma_minproduct_z_matrix *x,
    magma_minproduct_queue_t queue )
{

    dim3 grid( magma_minproduct_ceildiv( t.num_rows, BLOCK_SIZE ));
    magma_minproduct_int_t threads = BLOCK_SIZE;
    zjacobiupdate_kernel<<< grid, threads, 0 >>>( t.num_rows, t.num_cols, t.dval, b.dval, d.dval, x->dval );

    return MAGMA_minproduct_SUCCESS;
}










__global__ void 
zjacobispmvupdate_kernel(  
    int num_rows,
    int num_cols, 
    magma_minproductDoubleComplex * dval, 
    magma_minproduct_index_t * drowptr, 
    magma_minproduct_index_t * dcolind,
    magma_minproductDoubleComplex *t, 
    magma_minproductDoubleComplex *b, 
    magma_minproductDoubleComplex *d, 
    magma_minproductDoubleComplex *x ){



    int row = blockDim.x * blockIdx.x + threadIdx.x ;
    int j;

    if(row<num_rows){
        magma_minproductDoubleComplex dot = MAGMA_minproduct_Z_ZERO;
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
    maxiter     magma_minproduct_int_t
                number of Jacobi iterations   
                
    @param[in]
    A           magma_minproduct_z_matrix
                system matrix
                
    @param[in]
    t           magma_minproduct_z_matrix
                workspace
                
    @param[in]
    b           magma_minproduct_z_matrix
                RHS b
                
    @param[in]
    d           magma_minproduct_z_matrix
                vector with diagonal entries

    @param[out]
    x           magma_minproduct_z_matrix*
                iteration vector
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_z
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zjacobispmvupdate(
    magma_minproduct_int_t maxiter,
    magma_minproduct_z_matrix A,
    magma_minproduct_z_matrix t, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix d, 
    magma_minproduct_z_matrix *x,
    magma_minproduct_queue_t queue )
{

    // local variables
    magma_minproductDoubleComplex c_zero = MAGMA_minproduct_Z_ZERO, c_one = MAGMA_minproduct_Z_ONE;
    dim3 grid( magma_minproduct_ceildiv( t.num_rows, BLOCK_SIZE ));
    magma_minproduct_int_t threads = BLOCK_SIZE;

    for( magma_minproduct_int_t i=0; i<maxiter; i++ ) {
        // distinct routines imply synchronization
        // magma_minproduct_z_spmv( c_one, A, *x, c_zero, t, queue );                // t =  A * x
        // zjacobiupdate_kernel<<< grid, threads, 0 >>>( t.num_rows, t.num_cols, t.dval, b.dval, d.dval, x->dval );
        // merged in one implies asynchronous update
        zjacobispmvupdate_kernel<<< grid, threads, 0 >>>
            ( t.num_rows, t.num_cols, A.dval, A.drow, A.dcol, t.dval, b.dval, d.dval, x->dval );

    }

    return MAGMA_minproduct_SUCCESS;
}









