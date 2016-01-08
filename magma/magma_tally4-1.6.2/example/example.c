// This is a simple standalone example. See README.txt

#include <stdio.h>
#include <stdlib.h>

#include "cublas_v2.h"     // if you need CUBLAS, include before magma_tally4.h
#include "magma_tally4.h"
#include "magma_tally4_lapack.h"  // if you need BLAS & LAPACK

#include "zfill.h"         // code to fill matrix; replace with your application code


// ------------------------------------------------------------
// Solve A * X = B, where A and X are stored in CPU host memory.
// Internally, MAGMA_tally4 transfers data to the GPU device
// and uses a hybrid CPU + GPU algorithm.
void cpu_interface( magma_tally4_int_t n, magma_tally4_int_t nrhs )
{
    magma_tally4DoubleComplex *A=NULL, *X=NULL;
    magma_tally4_int_t *ipiv=NULL;
    magma_tally4_int_t lda  = n;
    magma_tally4_int_t ldx  = lda;
    magma_tally4_int_t info = 0;
    
    // magma_tally4 malloc_cpu routines are type-safe and align to memory boundaries,
    // but you can use malloc or new if you prefer.
    magma_tally4_zmalloc_cpu( &A, lda*n );
    magma_tally4_zmalloc_cpu( &X, ldx*nrhs );
    magma_tally4_imalloc_cpu( &ipiv, n );
    if ( A == NULL || X == NULL || ipiv == NULL ) {
        fprintf( stderr, "malloc failed\n" );
        goto cleanup;
    }
    
    zfill_matrix( n, n, A, lda );
    zfill_rhs( n, nrhs, X, ldx );
    
    magma_tally4_zgesv( n, 1, A, lda, ipiv, X, lda, &info );
    if ( info != 0 ) {
        fprintf( stderr, "magma_tally4_zgesv failed with info=%d\n", info );
    }
    
    // TODO: use result in X
    
cleanup:
    magma_tally4_free_cpu( A );
    magma_tally4_free_cpu( X );
    magma_tally4_free_cpu( ipiv );
}


// ------------------------------------------------------------
// Solve dA * dX = dB, where dA and dX are stored in GPU device memory.
// Internally, MAGMA_tally4 uses a hybrid CPU + GPU algorithm.
void gpu_interface( magma_tally4_int_t n, magma_tally4_int_t nrhs )
{
    magma_tally4DoubleComplex *dA=NULL, *dX=NULL;
    magma_tally4_int_t *ipiv=NULL;
    magma_tally4_int_t ldda = ((n+31)/32)*32;  // round up to multiple of 32 for best GPU performance
    magma_tally4_int_t lddx = ldda;
    magma_tally4_int_t info = 0;
    
    // magma_tally4 malloc (GPU) routines are type-safe,
    // but you can use cudaMalloc if you prefer.
    magma_tally4_zmalloc( &dA, ldda*n );
    magma_tally4_zmalloc( &dX, lddx*nrhs );
    magma_tally4_imalloc_cpu( &ipiv, n );  // ipiv always on CPU
    if ( dA == NULL || dX == NULL || ipiv == NULL ) {
        fprintf( stderr, "malloc failed\n" );
        goto cleanup;
    }
    
    zfill_matrix_gpu( n, n, dA, ldda );
    zfill_rhs_gpu( n, nrhs, dX, lddx );
    
    magma_tally4_zgesv_gpu( n, 1, dA, ldda, ipiv, dX, ldda, &info );
    if ( info != 0 ) {
        fprintf( stderr, "magma_tally4_zgesv_gpu failed with info=%d\n", info );
    }
    
    // TODO: use result in dX
    
cleanup:
    magma_tally4_free( dA );
    magma_tally4_free( dX );
    magma_tally4_free_cpu( ipiv );
}


// ------------------------------------------------------------
int main( int argc, char** argv )
{
    magma_tally4_init();
    
    magma_tally4_int_t n = 1000;
    magma_tally4_int_t nrhs = 1;
    
    printf( "using MAGMA_tally4 CPU interface\n" );
    cpu_interface( n, nrhs );

    printf( "using MAGMA_tally4 GPU interface\n" );
    gpu_interface( n, nrhs );
    
    magma_tally4_finalize();
    return 0;
}
