/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
       @author Azzam Haidar

*/
#include "common_magma_minproduct.h"
#include "magma_minproduct_templates.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE 512

#define BLOCK_SIZEx  32
#define BLOCK_SIZEy  16


//==============================================================================
//==============================================================================

__global__
void magma_minproduct_zlarf_kernel( int m, const magma_minproductDoubleComplex *dv, const magma_minproductDoubleComplex *dtau,
                         magma_minproductDoubleComplex *dc, int lddc )
{
    if ( !MAGMA_minproduct_Z_EQUAL(*dtau, MAGMA_minproduct_Z_ZERO) ) {
        const int tx = threadIdx.x;
        dc = dc + blockIdx.x * lddc;

        __shared__ magma_minproductDoubleComplex sum[ BLOCK_SIZE ];
        magma_minproductDoubleComplex tmp;

        /* perform  w := v**H * C  */
        if (tx==0)
            tmp = dc[0]; //since V[0] should be one
        else
            tmp = MAGMA_minproduct_Z_ZERO;
        for( int j = tx+1; j < m; j += BLOCK_SIZE ){
            tmp += MAGMA_minproduct_Z_MUL( MAGMA_minproduct_Z_CNJG( dv[j] ), dc[j] );
        }
        sum[tx] = tmp;
        magma_minproduct_sum_reduce< BLOCK_SIZE >( tx, sum );

        /*  C := C - v * w  */
        __syncthreads();
        tmp = - MAGMA_minproduct_Z_CNJG(*dtau) * sum[0];
        for( int j = m-tx-1; j>0 ; j -= BLOCK_SIZE )
             dc[j] += tmp * dv[j];

        if(tx==0) dc[0] += tmp;
    }
}

//==============================================================================
//==============================================================================

__global__
void magma_minproduct_zlarf_smkernel( int m, int n, magma_minproductDoubleComplex *dv, magma_minproductDoubleComplex *dtau,
                           magma_minproductDoubleComplex *dc, int lddc )
{
    if ( ! MAGMA_minproduct_Z_EQUAL(*dtau, MAGMA_minproduct_Z_ZERO) ) {
        const int i = threadIdx.x, col= threadIdx.y;

        for( int k = col; k < n; k += BLOCK_SIZEy ) {
            dc = dc + k * lddc;
    
            __shared__ magma_minproductDoubleComplex sum[ BLOCK_SIZEx ][ BLOCK_SIZEy + 1];
            magma_minproductDoubleComplex lsum;
    
            /*  w := v**H * C  */
            lsum = MAGMA_minproduct_Z_ZERO;
            for( int j = i; j < m; j += BLOCK_SIZEx ){
                if (j==0)
                   lsum += MAGMA_minproduct_Z_MUL( MAGMA_minproduct_Z_ONE, dc[j] );
                else
                   lsum += MAGMA_minproduct_Z_MUL( MAGMA_minproduct_Z_CNJG( dv[j] ), dc[j] );
            }
            sum[i][col] = lsum;
            magma_minproduct_sum_reduce_2d< BLOCK_SIZEx, BLOCK_SIZEy+1 >( i, col, sum );
    
            /*  C := C - v * w  */
            __syncthreads();
            magma_minproductDoubleComplex z__1 = - MAGMA_minproduct_Z_CNJG(*dtau) * sum[0][col];
            for( int j = m-i-1; j>=0 ; j -= BLOCK_SIZEx ) {
                 if (j==0)
                    dc[j] += z__1;
                 else
                    dc[j] += z__1 * dv[j];
            }
        }
    }
}

//==============================================================================

/*
    Apply a complex elementary reflector H to a complex M-by-N
    matrix C from the left. H is represented in the form
          H = I - tau * v * v**H
    where tau is a complex scalar and v is a complex vector.
    If tau = 0, then H is taken to be the unit matrix.

    To apply H**H (the conjugate transpose of H), supply conjg(tau)
    instead tau.

    This routine uses only one SM (block).
 */
extern "C" void
magma_minproduct_zlarf_sm(magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproductDoubleComplex *dv, magma_minproductDoubleComplex *dtau,
               magma_minproductDoubleComplex *dc, magma_minproduct_int_t lddc)
{
    dim3  blocks( 1 );
    dim3 threads( BLOCK_SIZEx, BLOCK_SIZEy );

    magma_minproduct_zlarf_smkernel<<< blocks, threads, 0, magma_minproduct_stream >>>( m, n, dv, dtau, dc, lddc );
}
//==============================================================================
/*
    Apply a complex elementary reflector H to a complex M-by-N
    matrix C from the left. H is represented in the form
          H = I - tau * v * v**H
    where tau is a complex scalar and v is a complex vector.
    If tau = 0, then H is taken to be the unit matrix.

    To apply H**H (the conjugate transpose of H), supply conjg(tau) 
    instead tau.

 */

extern "C" magma_minproduct_int_t
magma_minproduct_zlarf_gpu(
    magma_minproduct_int_t m,  magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr dv,
    magma_minproductDoubleComplex_const_ptr dtau,
    magma_minproductDoubleComplex_ptr dC,  magma_minproduct_int_t lddc)
{
    dim3 grid( n, 1, 1 );
    dim3 threads( BLOCK_SIZE );
    if ( n > 0 ) {
        magma_minproduct_zlarf_kernel<<< grid, threads, 0, magma_minproduct_stream >>>( m, dv, dtau, dC, lddc);
    }

    // The computation can be done on 1 SM with the following routine.
    // magma_minproduct_zlarf_sm(m, n, dv, dtau, dc, lddc);

    return MAGMA_minproduct_SUCCESS;
}

//==============================================================================
