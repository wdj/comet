/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal d -> s
*/
#include "common_magma_tally2.h"
#include "commonblas_d.h"

/*
 * daxpy computes c += alpha*b, where b and c are 16-element vectors.
 */
static __device__ void daxpy(
    double alpha,
    const double* __restrict__ b,
    double*       __restrict__ c )
{
    c[0]  += alpha * b[0];
    c[1]  += alpha * b[1];
    c[2]  += alpha * b[2];
    c[3]  += alpha * b[3];
    c[4]  += alpha * b[4];
    c[5]  += alpha * b[5];
    c[6]  += alpha * b[6];
    c[7]  += alpha * b[7];
    c[8]  += alpha * b[8];
    c[9]  += alpha * b[9];
    c[10] += alpha * b[10];
    c[11] += alpha * b[11];
    c[12] += alpha * b[12];
    c[13] += alpha * b[13];
    c[14] += alpha * b[14];
    c[15] += alpha * b[15];
}


/**
    Purpose:
    --------
    This routine computes
        C = alpha * A^T*B^T + beta * C

    B is put into shared memory
    Parameters Used:
        blk_M=64 blk_N=16 blk_K=16 nthd_x=16 nthd_y=4

    This code should run for any matrix size.
    This kernel outperforms cuda-2.2 when m, n, k >= 512

    @ingroup magma_tally2_dblas3
    ********************************************************************/
__global__ void
dgemm_kernel_T_T_64_16_16_16_4(
    double*       __restrict__ C,
    const double* __restrict__ A,
    const double* __restrict__ B,
    int m, int n, int k,
    int lda, int ldb, int ldc,
    double alpha, double beta )
{
    __shared__ double Bb[16][17];

    const int tx = threadIdx.x;
    const int ty = threadIdx.y;

    int ibx = blockIdx.x * 64;
    int iby = blockIdx.y * 16;

    const int idt = ty * 16 + tx;

    /*
        Taking care of invalid memory access in dimension M
    */
    if ( ibx + idt >= m )
        A += ibx + 0;
    else
        A += ibx + idt;

    C += __mul24(ibx + idt, ldc) + iby;
    B += tx + __mul24(iby, ldb);

    /*
        These variables guide the threads to avoid invalid memory
        accesses in dimension N
        Simply it's the stopping criterion.
        or you can say that access index wraps around to a valid
        memory location.
    */
    int s1=0, s2=4*ldb, s3=8*ldb, s4=12*ldb;
    
    if ( iby+ty    >= n ) { s1=1;  s2=0*ldb;  s3=0*ldb;  s4=0*ldb; } else
    if ( iby+ty+4  >= n ) { s1=0;  s2=0*ldb;  s3=0*ldb;  s4=0*ldb; } else
    if ( iby+ty+8  >= n ) { s1=0;  s2=4*ldb;  s3=0*ldb;  s4=0*ldb; } else
    if ( iby+ty+12 >= n ) { s1=0;  s2=4*ldb;  s3=8*ldb;  s4=0*ldb; }

    if ( s1 == 0 )
        B += __mul24(ty, ldb);
    else
        s1 = 0;

    const double *Bend = B + k - k % 16;
    double Cb[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    if ( k > 15 ) {
        do {
            double Ab[4] = { A[0], A[lda], A[2*lda], A[3*lda] };

            Bb[tx][ty+0 ] = B[s1];
            Bb[tx][ty+4 ] = B[s2];
            Bb[tx][ty+8 ] = B[s3];
            Bb[tx][ty+12] = B[s4];
            
            __syncthreads();
      
            A += 4 * lda;
            daxpy( Ab[0], &Bb[0][0], Cb );  Ab[0] = A[0*lda];
            daxpy( Ab[1], &Bb[1][0], Cb );  Ab[1] = A[1*lda];
            daxpy( Ab[2], &Bb[2][0], Cb );  Ab[2] = A[2*lda];
            daxpy( Ab[3], &Bb[3][0], Cb );  Ab[3] = A[3*lda];
            
            A += 4 * lda;
            daxpy( Ab[0], &Bb[4][0], Cb );  Ab[0] = A[0*lda];
            daxpy( Ab[1], &Bb[5][0], Cb );  Ab[1] = A[1*lda];
            daxpy( Ab[2], &Bb[6][0], Cb );  Ab[2] = A[2*lda];
            daxpy( Ab[3], &Bb[7][0], Cb );  Ab[3] = A[3*lda];
            
            A += 4 * lda;
            daxpy( Ab[0], &Bb[8][0],  Cb );  Ab[0] = A[0*lda];
            daxpy( Ab[1], &Bb[9][0],  Cb );  Ab[1] = A[1*lda];
            daxpy( Ab[2], &Bb[10][0], Cb );  Ab[2] = A[2*lda];
            daxpy( Ab[3], &Bb[11][0], Cb );  Ab[3] = A[3*lda];

            A += 4 * lda;
            daxpy( Ab[0], &Bb[12][0], Cb );
            daxpy( Ab[1], &Bb[13][0], Cb );
            daxpy( Ab[2], &Bb[14][0], Cb );
            daxpy( Ab[3], &Bb[15][0], Cb );

            B += 16;

            __syncthreads();
        } while (B < Bend);
    }
    
    /*
        Common sub expression elimination.
    */
    ibx = ibx + idt - m;

    /*
        remembering k dimension
    */
    ldb = m = k;

    /*
        k changed to support the generic case and reuse valuable registers
    */
    k = k % 16;

    m -= k;

    /*
        Here we are taking care of k % dim_k portions
    */
    if ( k != 0 ) {
        /*
            Avoid Invalid Memory access in dimension K
            If some thread enters this if ( ) block first access to B
            should be valid as K isn't divisible by blk_K
            Note that dimension N has been taken care of by s1, s2, s3, s4
            But depending upon K and thread index tx, some memory access
            may be still invalid, so take care of them now by
            setting s1, s2, s3, s4 = 0
            B might have been advanced in the previous loop, take care
            of that, this is about right bottom corner.
        */
        if ( m + tx >= ldb ) {
            s1 = s2 = s3 = s4 = 0;
            B -= tx;
        }

        Bb[tx][ty+0 ] = B[s1];
        Bb[tx][ty+4 ] = B[s2];
        Bb[tx][ty+8 ] = B[s3];
        Bb[tx][ty+12] = B[s4];
        __syncthreads();

        for(int i=0; i < k; i++) {
            daxpy( A[0], &Bb[i+0][0], Cb );
            A += lda;
        }
    }

    /*
        Now taking care of dimension M, N that doesnt fit into blocks
    */
    if ( (iby + 16) >= n ) {
        lda = n - iby;
    }
    else {
        lda = 16;
    }
    if ( ibx >= 0 )
        lda = 0;
    else
        lda = lda;

    switch(lda) {
        case 16:
                C[0]  = alpha * Cb[0]  + beta * C[0];
                C[1]  = alpha * Cb[1]  + beta * C[1];
                C[2]  = alpha * Cb[2]  + beta * C[2];
                C[3]  = alpha * Cb[3]  + beta * C[3];
                C[4]  = alpha * Cb[4]  + beta * C[4];
                C[5]  = alpha * Cb[5]  + beta * C[5];
                C[6]  = alpha * Cb[6]  + beta * C[6];
                C[7]  = alpha * Cb[7]  + beta * C[7];
                C[8]  = alpha * Cb[8]  + beta * C[8];
                C[9]  = alpha * Cb[9]  + beta * C[9];
                C[10] = alpha * Cb[10] + beta * C[10];
                C[11] = alpha * Cb[11] + beta * C[11];
                C[12] = alpha * Cb[12] + beta * C[12];
                C[13] = alpha * Cb[13] + beta * C[13];
                C[14] = alpha * Cb[14] + beta * C[14];
                C[15] = alpha * Cb[15] + beta * C[15];
                break;
        case 15:
                C[0]  = alpha * Cb[0]  + beta * C[0];
                C[1]  = alpha * Cb[1]  + beta * C[1];
                C[2]  = alpha * Cb[2]  + beta * C[2];
                C[3]  = alpha * Cb[3]  + beta * C[3];
                C[4]  = alpha * Cb[4]  + beta * C[4];
                C[5]  = alpha * Cb[5]  + beta * C[5];
                C[6]  = alpha * Cb[6]  + beta * C[6];
                C[7]  = alpha * Cb[7]  + beta * C[7];
                C[8]  = alpha * Cb[8]  + beta * C[8];
                C[9]  = alpha * Cb[9]  + beta * C[9];
                C[10] = alpha * Cb[10] + beta * C[10];
                C[11] = alpha * Cb[11] + beta * C[11];
                C[12] = alpha * Cb[12] + beta * C[12];
                C[13] = alpha * Cb[13] + beta * C[13];
                C[14] = alpha * Cb[14] + beta * C[14];
                break;
        case 14:
                C[0]  = alpha * Cb[0]  + beta * C[0];
                C[1]  = alpha * Cb[1]  + beta * C[1];
                C[2]  = alpha * Cb[2]  + beta * C[2];
                C[3]  = alpha * Cb[3]  + beta * C[3];
                C[4]  = alpha * Cb[4]  + beta * C[4];
                C[5]  = alpha * Cb[5]  + beta * C[5];
                C[6]  = alpha * Cb[6]  + beta * C[6];
                C[7]  = alpha * Cb[7]  + beta * C[7];
                C[8]  = alpha * Cb[8]  + beta * C[8];
                C[9]  = alpha * Cb[9]  + beta * C[9];
                C[10] = alpha * Cb[10] + beta * C[10];
                C[11] = alpha * Cb[11] + beta * C[11];
                C[12] = alpha * Cb[12] + beta * C[12];
                C[13] = alpha * Cb[13] + beta * C[13];
                break;
        case 13:
                C[0]  = alpha * Cb[0]  + beta * C[0];
                C[1]  = alpha * Cb[1]  + beta * C[1];
                C[2]  = alpha * Cb[2]  + beta * C[2];
                C[3]  = alpha * Cb[3]  + beta * C[3];
                C[4]  = alpha * Cb[4]  + beta * C[4];
                C[5]  = alpha * Cb[5]  + beta * C[5];
                C[6]  = alpha * Cb[6]  + beta * C[6];
                C[7]  = alpha * Cb[7]  + beta * C[7];
                C[8]  = alpha * Cb[8]  + beta * C[8];
                C[9]  = alpha * Cb[9]  + beta * C[9];
                C[10] = alpha * Cb[10] + beta * C[10];
                C[11] = alpha * Cb[11] + beta * C[11];
                C[12] = alpha * Cb[12] + beta * C[12];
                break;
        case 12:
                C[0]  = alpha * Cb[0]  + beta * C[0];
                C[1]  = alpha * Cb[1]  + beta * C[1];
                C[2]  = alpha * Cb[2]  + beta * C[2];
                C[3]  = alpha * Cb[3]  + beta * C[3];
                C[4]  = alpha * Cb[4]  + beta * C[4];
                C[5]  = alpha * Cb[5]  + beta * C[5];
                C[6]  = alpha * Cb[6]  + beta * C[6];
                C[7]  = alpha * Cb[7]  + beta * C[7];
                C[8]  = alpha * Cb[8]  + beta * C[8];
                C[9]  = alpha * Cb[9]  + beta * C[9];
                C[10] = alpha * Cb[10] + beta * C[10];
                C[11] = alpha * Cb[11] + beta * C[11];
                break;
        case 11:
                C[0]  = alpha * Cb[0]  + beta * C[0];
                C[1]  = alpha * Cb[1]  + beta * C[1];
                C[2]  = alpha * Cb[2]  + beta * C[2];
                C[3]  = alpha * Cb[3]  + beta * C[3];
                C[4]  = alpha * Cb[4]  + beta * C[4];
                C[5]  = alpha * Cb[5]  + beta * C[5];
                C[6]  = alpha * Cb[6]  + beta * C[6];
                C[7]  = alpha * Cb[7]  + beta * C[7];
                C[8]  = alpha * Cb[8]  + beta * C[8];
                C[9]  = alpha * Cb[9]  + beta * C[9];
                C[10] = alpha * Cb[10] + beta * C[10];
                break;
        case 10:
                C[0]  = alpha * Cb[0]  + beta * C[0];
                C[1]  = alpha * Cb[1]  + beta * C[1];
                C[2]  = alpha * Cb[2]  + beta * C[2];
                C[3]  = alpha * Cb[3]  + beta * C[3];
                C[4]  = alpha * Cb[4]  + beta * C[4];
                C[5]  = alpha * Cb[5]  + beta * C[5];
                C[6]  = alpha * Cb[6]  + beta * C[6];
                C[7]  = alpha * Cb[7]  + beta * C[7];
                C[8]  = alpha * Cb[8]  + beta * C[8];
                C[9]  = alpha * Cb[9]  + beta * C[9];
                break;
        case 9:
                C[0]  = alpha * Cb[0]  + beta * C[0];
                C[1]  = alpha * Cb[1]  + beta * C[1];
                C[2]  = alpha * Cb[2]  + beta * C[2];
                C[3]  = alpha * Cb[3]  + beta * C[3];
                C[4]  = alpha * Cb[4]  + beta * C[4];
                C[5]  = alpha * Cb[5]  + beta * C[5];
                C[6]  = alpha * Cb[6]  + beta * C[6];
                C[7]  = alpha * Cb[7]  + beta * C[7];
                C[8]  = alpha * Cb[8]  + beta * C[8];
                break;
        case 8:
                C[0]  = alpha * Cb[0]  + beta * C[0];
                C[1]  = alpha * Cb[1]  + beta * C[1];
                C[2]  = alpha * Cb[2]  + beta * C[2];
                C[3]  = alpha * Cb[3]  + beta * C[3];
                C[4]  = alpha * Cb[4]  + beta * C[4];
                C[5]  = alpha * Cb[5]  + beta * C[5];
                C[6]  = alpha * Cb[6]  + beta * C[6];
                C[7]  = alpha * Cb[7]  + beta * C[7];
                break;
        case 7:
                C[0]  = alpha * Cb[0]  + beta * C[0];
                C[1]  = alpha * Cb[1]  + beta * C[1];
                C[2]  = alpha * Cb[2]  + beta * C[2];
                C[3]  = alpha * Cb[3]  + beta * C[3];
                C[4]  = alpha * Cb[4]  + beta * C[4];
                C[5]  = alpha * Cb[5]  + beta * C[5];
                C[6]  = alpha * Cb[6]  + beta * C[6];
                break;
        case 6:
                C[0]  = alpha * Cb[0]  + beta * C[0];
                C[1]  = alpha * Cb[1]  + beta * C[1];
                C[2]  = alpha * Cb[2]  + beta * C[2];
                C[3]  = alpha * Cb[3]  + beta * C[3];
                C[4]  = alpha * Cb[4]  + beta * C[4];
                C[5]  = alpha * Cb[5]  + beta * C[5];
                break;
        case 5:
                C[0]  = alpha * Cb[0]  + beta * C[0];
                C[1]  = alpha * Cb[1]  + beta * C[1];
                C[2]  = alpha * Cb[2]  + beta * C[2];
                C[3]  = alpha * Cb[3]  + beta * C[3];
                C[4]  = alpha * Cb[4]  + beta * C[4];
                break;
        case 4:
                C[0]  = alpha * Cb[0]  + beta * C[0];
                C[1]  = alpha * Cb[1]  + beta * C[1];
                C[2]  = alpha * Cb[2]  + beta * C[2];
                C[3]  = alpha * Cb[3]  + beta * C[3];
                break;
        case 3:
                C[0]  = alpha * Cb[0]  + beta * C[0];
                C[1]  = alpha * Cb[1]  + beta * C[1];
                C[2]  = alpha * Cb[2]  + beta * C[2];
                break;
        case 2:
                C[0]  = alpha * Cb[0]  + beta * C[0];
                C[1]  = alpha * Cb[1]  + beta * C[1];
                break;
        case 1:
                C[0]  = alpha * Cb[0]  + beta * C[0];
                break;
        case 0:
                break;
    }
}


extern "C" void
magma_tally2blas_dgemm_T_T_64_16_16_16_4(
    double *C, const double *A, const double *B,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2_int_t lda, magma_tally2_int_t ldb, magma_tally2_int_t ldc,
    double alpha, double beta )
{
    dim3 threads( 16, 4 );
    dim3 grid( (m - 1)/64 + 1, (n - 1)/16 + 1 );
    dgemm_kernel_T_T_64_16_16_16_4<<< grid, threads, 0, magma_tally2_stream >>>
        ( C, A, B, m, n, k, lda, ldb, ldc, alpha, beta );
}