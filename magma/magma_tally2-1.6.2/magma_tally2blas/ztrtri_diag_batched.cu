/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> c d s

       @author Peng Du
       @author Tingxing Dong
       @author Mark Gates
       @author Azzam Haidar
       
       File named ztrtri_diag.cu to avoid name conflict with src/ztrtri.o
       in the library. The actual kernels are in ztrtri_lower.cu and ztrtri_upper.cu
*/

#include "common_magma_tally2.h"
#include "ztrtri.h"


/**
    Inverts the NB x NB diagonal blocks of a triangular matrix.
    This routine is used in ztrsm.
    
    Same as ztrtri_diag, but adds queue argument.
    
    @ingroup magma_tally2_zblas3
    ********************************************************************/
/**
    Purpose
    -------
    ztrtri_diag inverts the NB x NB diagonal blocks of A.

    Arguments
    ----------
    @param[in]
    uplo    magma_tally2_uplo_t.
            On entry, uplo specifies whether the matrix A is an upper or
            lower triangular matrix as follows:
      -     = Magma_tally2Upper:  A is an upper triangular matrix.
      -     = Magma_tally2Lower:  A is a  lower triangular matrix.

    @param[in]
    diag    magma_tally2_diag_t.
            On entry, diag specifies whether or not A is unit triangular
            as follows:
      -     = Magma_tally2Unit:     A is assumed to be unit triangular.
      -     = Magma_tally2NonUnit:  A is not assumed to be unit triangular.

    @param[in]
    n       INTEGER.
            On entry, n specifies the order of the matrix A. N >= 0.

    @param[in]
    dA_array      COMPLEX_16 array of dimension ( ldda, n )
            The triangular matrix A.
    \n
            If UPLO = 'U', the leading N-by-N upper triangular part of A
            contains the upper triangular matrix, and the strictly lower
            triangular part of A is not referenced.
    \n
            If UPLO = 'L', the leading N-by-N lower triangular part of A
            contains the lower triangular matrix, and the strictly upper
            triangular part of A is not referenced.
    \n
            If DIAG = 'U', the diagonal elements of A are also not referenced
            and are assumed to be 1.

    @param[in]
    ldda    INTEGER.
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[out]
    dinvA_array COMPLEX_16 array of dimension (NB, ((n+NB-1)/NB)*NB),
            where NB = 128.
            On exit, contains inverses of the NB-by-NB diagonal blocks of A.

    @param[in]
    queue   magma_tally2_queue_t
            Queue to execute in.

    @ingroup magma_tally2_zblas3
    ********************************************************************/
extern "C" void
magma_tally2blas_ztrtri_diag_batched(
    magma_tally2_uplo_t uplo, magma_tally2_diag_t diag, magma_tally2_int_t n,
    magma_tally2DoubleComplex const * const *dA_array, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex **dinvA_array, 
    magma_tally2_int_t resetozero, magma_tally2_int_t batchCount, magma_tally2_queue_t queue)
{
    magma_tally2_int_t info = 0;
    if (uplo != Magma_tally2Lower && uplo != Magma_tally2Upper)
        info = -1;
    else if (diag != Magma_tally2NonUnit && diag != Magma_tally2Unit)
        info = -2;
    else if (n < 0)
        info = -3;
    else if (ldda < n)
        info = -5;

    if (info != 0) {
        magma_tally2_xerbla( __func__, -(info) );
        return;  //info
    }
    
    int nblocks = (n + IB - 1)/IB;


    if(resetozero)
    { 
        magma_tally2blas_zlaset_batched(Magma_tally2Full, ((n+NB-1)/NB)*NB, NB, MAGMA_tally2_Z_ZERO, MAGMA_tally2_Z_ZERO, dinvA_array, ((n+NB-1)/NB)*NB, batchCount, queue);
       //magma_tally2blas_zmemset_batched( dinvA_array, ((n+NB-1)/NB)*NB*NB, batchCount, queue);
    }
    // if someone want to use cudamemset he need to set the whole vectors 
    // of initial size otherwise it is a bug and thus need to have dinvA_length 
    // in input parameter and has been tested and was slower.
    //was not the largest size computed by the high API getrf_batched then it is bug and need to use magma_tally2blas_zlaset_batched


    if ( uplo == Magma_tally2Lower ) {
        // invert diagonal IB x IB inner blocks
        dim3 diaggrid( nblocks, 1, batchCount );  // emulate 3D grid
        ztrtri_diag_lower_kernel_batched<<< diaggrid, IB, 0, queue >>>( diag, n, dA_array, ldda, dinvA_array );

        // build up NB x NB blocks (assuming IB=16 here):
        // use   16 x 16  blocks to build  32 x 32  blocks,  1 x (1 x npages) grid,  4 x 4 threads;
        // then  32 x 32  blocks to build  64 x 64  blocks,  1 x (2 x npages) grid,  8 x 4 threads;
        // then  64 x 64  blocks to build 128 x 128 blocks,  1 x (4 x npages) grid, 16 x 4 threads;
        // then 128 x 128 blocks to build 256 x 256 blocks,  2 x (8 x npages) grid, 16 x 4 threads.
        for( int jb=IB; jb < NB; jb *= 2 ) {
            int kb = jb*2;
            int npages = (n + kb - 1)/kb;
            dim3 threads( (jb <= 32 ? jb/4 : 16), 4 );
            dim3 grid( jb/(threads.x*threads.y), npages*(jb/16), batchCount );  // emulate 3D grid: NX * (NY*npages), for CUDA ARCH 1.x
            
            //printf( "n %d, jb %d, grid %d x %d (%d x %d)\n", n, jb, grid.x, grid.y, grid.y / npages, npages );
            switch (jb) {
                case 16:
                    triple_zgemm16_part1_lower_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    triple_zgemm16_part2_lower_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    break;
                case 32:
                    triple_zgemm32_part1_lower_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    triple_zgemm32_part2_lower_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    break;
                case 64:
                    triple_zgemm64_part1_lower_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    triple_zgemm64_part2_lower_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    break;
                default:
                    triple_zgemm_above64_part1_lower_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    triple_zgemm_above64_part2_lower_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    triple_zgemm_above64_part3_lower_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    break;
            }
            if ( kb >= n ) break;
        }
    }
    else {
        dim3 diaggrid( nblocks, 1, batchCount );  // emulate 3D grid
        ztrtri_diag_upper_kernel_batched<<< diaggrid, IB, 0, queue >>>( diag, n, dA_array, ldda, dinvA_array );

        // update the inverse up to the size of IB
        for( int jb=IB; jb < NB; jb*=2 ) {
            int kb = jb*2;
            int npages = (n + kb - 1)/kb;
            dim3 threads( (jb <= 32 ? jb/4 : 16), 4 );
            dim3 grid( jb/(threads.x*threads.y), npages*(jb/16), batchCount );  // emulate 3D grid: NX * (NY*npages), for CUDA ARCH 1.x
            
            switch (jb) {
                case 16:
                    triple_zgemm16_part1_upper_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    triple_zgemm16_part2_upper_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    break;
                case 32:
                    triple_zgemm32_part1_upper_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    triple_zgemm32_part2_upper_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    break;
                case 64:
                    triple_zgemm64_part1_upper_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    triple_zgemm64_part2_upper_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    break;
                default:
                    triple_zgemm_above64_part1_upper_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    triple_zgemm_above64_part2_upper_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    triple_zgemm_above64_part3_upper_kernel_batched<<< grid, threads, 0, queue >>>( n, dA_array, ldda, dinvA_array, jb, npages );
                    break;
            }
            if ( kb >= n ) break;
        }
    }
}


