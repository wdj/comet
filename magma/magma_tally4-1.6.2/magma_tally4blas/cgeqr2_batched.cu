/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @author Azzam Haidar
       @author Tingxing Dong

       @generated from zgeqr2_batched.cu normal z -> c, Fri Jan 30 19:00:10 2015
*/

#include "common_magma_tally4.h"
#include "magma_tally4_templates.h"
#include "batched_kernel_param.h"

#define BLOCK_SIZE 256
#define PRECISION_c

#define dA(a_1,a_2) (dA  + (a_1) + (a_2)*(local_lda))

#define COMPLEX

//==============================================================================
static __device__ void
clarfg_device(
    int n,
    magma_tally4FloatComplex* dalpha, magma_tally4FloatComplex* dx, int incx,
    magma_tally4FloatComplex* dtau,  float* swork, float* sscale, magma_tally4FloatComplex* scale)
{

    const int tx = threadIdx.x;

    magma_tally4FloatComplex tmp;
    
    // find max of [dalpha, dx], to use as scaling to avoid unnecesary under- and overflow    

    if ( tx == 0 ) {
        tmp = *dalpha;
        #ifdef COMPLEX
        swork[tx] = max( fabs(real(tmp)), fabs(imag(tmp)) );
        #else
        swork[tx] = fabs(tmp);
        #endif
    }
    else {
        swork[tx] = 0;
    }
    if(tx<BLOCK_SIZE)
    {
        for( int j = tx; j < n-1; j += BLOCK_SIZE ) {
            tmp = dx[j*incx];
            #ifdef COMPLEX
            swork[tx] = max( swork[tx], max( fabs(real(tmp)), fabs(imag(tmp)) ));
            #else
            swork[tx] = max( swork[tx], fabs(tmp) );
            #endif
         }
    }

    magma_tally4_max_reduce<BLOCK_SIZE>( tx, swork );

    if ( tx == 0 )
        *sscale = swork[0];
    __syncthreads();
    
    // sum norm^2 of dx/sscale
    // dx has length n-1
    if(tx<BLOCK_SIZE) swork[tx] = 0;
    if ( *sscale > 0 ) {
        if(tx<BLOCK_SIZE)
        {
            for( int j = tx; j < n-1; j += BLOCK_SIZE ) {
                tmp = dx[j*incx] / *sscale;
                swork[tx] += real(tmp)*real(tmp) + imag(tmp)*imag(tmp);
            }
        }
        magma_tally4_sum_reduce<BLOCK_SIZE>( tx, swork );

    }
    
    if ( tx == 0 ) {
        magma_tally4FloatComplex alpha = *dalpha;

        if ( swork[0] == 0 && imag(alpha) == 0 ) {
            // H = I
            *dtau = MAGMA_tally4_C_ZERO;
        }
        else {
            // beta = norm( [dalpha, dx] )
            float beta;
            tmp  = alpha / *sscale;
            beta = *sscale * sqrt( real(tmp)*real(tmp) + imag(tmp)*imag(tmp) + swork[0] );
            beta = -copysign( beta, real(alpha) );
            // todo: deal with badly scaled vectors (see lapack's larfg)
            *dtau   = MAGMA_tally4_C_MAKE( (beta - real(alpha)) / beta, -imag(alpha) / beta );
            *dalpha = MAGMA_tally4_C_MAKE( beta, 0 );
            *scale = 1 / (alpha - beta);
        }
    }
    
    // scale x (if norm was not 0)
    __syncthreads();
    if ( swork[0] != 0 ) {
        if(tx<BLOCK_SIZE)
        {
            for( int j = tx; j < n-1; j += BLOCK_SIZE ) {
                dx[j*incx] *= *scale;
            }
        }
    }

}




//==============================================================================

static __device__
void clarfx_device( int m, int n,  magma_tally4FloatComplex *v, magma_tally4FloatComplex *tau,
                         magma_tally4FloatComplex *dc, magma_tally4_int_t ldc, magma_tally4FloatComplex* sum)
{


    if(n <=0) return ;
    if (MAGMA_tally4_C_EQUAL(*tau, MAGMA_tally4_C_ZERO) )  return; // check singularity

    const int tx = threadIdx.x;

    magma_tally4FloatComplex lsum;
       
    for(int k=0;k<n;k++)
    {
        /* perform  w := v' * C  */
        if(tx<BLOCK_SIZE)
        {
            if (tx==0)
                lsum = dc[0+ldc*k]; //since V[0] should be one
            else
                lsum = MAGMA_tally4_C_ZERO;
            for( int j = tx+1; j < m; j += BLOCK_SIZE ){
                lsum += MAGMA_tally4_C_MUL( MAGMA_tally4_C_CNJG( v[j] ), dc[j+ldc*k] );
            }

            sum[tx] = lsum;
        }

        magma_tally4_sum_reduce< BLOCK_SIZE >( tx, sum );
        __syncthreads();

        magma_tally4FloatComplex z__1 = - MAGMA_tally4_C_CNJG(*tau) * sum[0];
        /*  C := C - v * w  */
        if(tx<BLOCK_SIZE)
        {    
           for( int j = tx+1; j<m ; j += BLOCK_SIZE )
                 dc[j+ldc*k] += z__1 * v[j];
        }
        if(tx==0) dc[0+ldc*k] += z__1;

        __syncthreads();


    } 
}

//==============================================================================

extern __shared__ magma_tally4FloatComplex shared_data[];


__global__
void cgeqr2_sm_kernel_batched( int m, int n, magma_tally4FloatComplex** dA_array, magma_tally4_int_t lda,
                               magma_tally4FloatComplex **dtau_array)
{

    magma_tally4FloatComplex* dA = dA_array[blockIdx.z];
    magma_tally4FloatComplex* dtau = dtau_array[blockIdx.z];

    magma_tally4FloatComplex *sdata = (magma_tally4FloatComplex*)shared_data;

    const int tx = threadIdx.x;

    __shared__ magma_tally4FloatComplex scale;
    __shared__ magma_tally4FloatComplex sum[ BLOCK_SIZE ];

    __shared__ float swork[ BLOCK_SIZE ];
    __shared__ float sscale;
    
    //load data from global to shared memory
    for(int s=0;s<n;s++)
    {
        for( int j = tx; j < m; j += BLOCK_SIZE )
        {
            sdata[j + s * m] = dA[j + s * lda] ;
        }
    }


    __syncthreads();
 
    for(int s=0; s<min(m,n); s++)
    {

       //lapack clarfg, compute the norm, scale and generate the householder vector   

       clarfg_device(m-s, &(sdata[s+s*m]), &(sdata[s+1+s*m]), 1, dtau+s, swork, &sscale, &scale); 
       __syncthreads();

       
       //update the trailing matix with the householder
       clarfx_device(m-s, n-(s+1), &(sdata[s+s*m]), dtau+s,&(sdata[s+(s+1)*m]), m, sum);

    }// end of s

    //copy back to global memory
    for(int s=0;s<n;s++)
    {
        for( int j = tx; j < m; j += BLOCK_SIZE )
        {
            dA[j + s * lda] = sdata[j + s * m];
        }
    }

}




//==============================================================================



static __device__
void cgeqr2_device( magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4FloatComplex* dA, magma_tally4_int_t lda,
                               magma_tally4FloatComplex *dtau, 
                               magma_tally4FloatComplex *sdata,
                               magma_tally4FloatComplex *sum,
                               float *swork,
                               magma_tally4FloatComplex *scale,
                               float *sscale)
{

    const int tx = threadIdx.x;


    for(int s=0; s<min(m,n); s++)
    {
       //load one vector in shared memory: sdata
       for( int j = tx; j < m-s; j += BLOCK_SIZE )
       {
           sdata[j] = dA[s + j + s * lda] ;
       }

       __syncthreads();

       //if(tx== 0) printf("m-s=%d",m-s);
       //lapack clarfg, compute the norm, scale and generate the householder vector   
       clarfg_device(m-s, sdata, &(sdata[1]), 1, dtau+s, swork, sscale, scale); 

       __syncthreads();

       //update the trailing matix with the householder
       clarfx_device(m-s, n-(s+1), sdata, dtau+s, &(dA[s+(s+1)*lda]), lda, sum);

       for( int j = tx; j < m-s; j += BLOCK_SIZE )
       {
           dA[s + j + s * lda] = sdata[j];
       }

       __syncthreads();

    }// end of s

}



//==============================================================================

__global__
void cgeqr2_kernel_batched( int m, int n, magma_tally4FloatComplex** dA_array, magma_tally4_int_t lda,
                               magma_tally4FloatComplex **dtau_array)
{

    magma_tally4FloatComplex* dA = dA_array[blockIdx.z];
    magma_tally4FloatComplex* dtau = dtau_array[blockIdx.z];

    magma_tally4FloatComplex *sdata = (magma_tally4FloatComplex*)shared_data;


    __shared__ magma_tally4FloatComplex scale;
    __shared__ magma_tally4FloatComplex sum[ BLOCK_SIZE ];

    __shared__ float swork[ BLOCK_SIZE ];
    __shared__ float sscale;

    cgeqr2_device(m, n, dA, lda, dtau, sdata, sum, swork, &scale, &sscale); 
 
}




//==============================================================================


/**
    Purpose
    -------
    CGEQR2 computes a QR factorization of a complex m by n matrix A:
    A = Q * R.

    This expert routine requires two more arguments than the standard
    cgeqr2, namely, dT and ddA, explained below. The storage for A is
    also not as in the LAPACK's cgeqr2 routine (see below).

    The first is used to output the triangular
    n x n factor T of the block reflector used in the factorization.
    The second holds the diagonal nxn blocks of A, i.e., the diagonal
    submatrices of R.

    This version implements the right-looking QR with non-blocking.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in,out]
    dA      COMPLEX array, dimension (LDA,N)
            On entry, the m by n matrix A.
            On exit, the unitary matrix Q as a
            product of elementary reflectors (see Further Details).
    \n
            the elements on and above the diagonal of the array
            contain the min(m,n) by n upper trapezoidal matrix R (R is
            upper triangular if m >= n); the elements below the diagonal,
            with the array TAU, represent the unitary matrix Q as a
            product of elementary reflectors (see Further Details).

    @param[in]
    lda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    @param[out]
    dtau    COMPLEX array, dimension (min(M,N))
            The scalar factors of the elementary reflectors (see Further
            Details).

    @param[out]
    dT      COMPLEX array, dimension N x N.
            Stores the triangular N x N factor T of the block reflector
            used in the factorization. The lower triangular part is 0.


    @param
    dwork   (workspace) COMPLEX array, dimension (N) * ( sizeof(float) + sizeof(magma_tally4FloatComplex)) 

    @param[out]
    info    INTEGER
      -     = 0: successful exit
      -     < 0: if INFO = -i, the i-th argument had an illegal value

    Further Details
    ---------------
    The matrix Q is represented as a product of elementary reflectors

       Q = H(1) H(2) . . . H(k), where k = min(m,n).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a complex scalar, and v is a complex vector with
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
    and tau in TAU(i).

    @ingroup magma_tally4_cgeqrf_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_cgeqr2_batched(magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4FloatComplex **dA_array,
                  magma_tally4_int_t lda, magma_tally4FloatComplex **dtau_array,
                  magma_tally4_int_t *info_array, magma_tally4_int_t batchCount, magma_tally4_queue_t queue)
{
    
    magma_tally4_int_t k;

    /* Check arguments */
    magma_tally4_int_t arginfo = 0;
    if (m < 0)
        arginfo = -1;
    else if (n < 0)
        arginfo = -2;
    else if (lda < max(1,m))
        arginfo = -4;

    if (arginfo != 0) {
        magma_tally4_xerbla( __func__, -(arginfo) );
        return arginfo;
    }


    k = min(m,n);

    dim3 blocks(1, 1, batchCount);
    dim3 threads(BLOCK_SIZE);

    if(sizeof(magma_tally4FloatComplex)*(m*k) <= 128 /*sizeof(magma_tally4FloatComplex) * 128 * k*/) // there are some static shared memory besides of dynamic ones 
    {   
        //load panel in shared memory and factorize it and copy back to gloabl memory
        //intend for small panel to avoid overfill of shared memory.
        //this kernel is composed of device routine and thus clean
        cgeqr2_sm_kernel_batched<<< blocks, threads, sizeof(magma_tally4FloatComplex)*(m*k), queue >>>
                                      (m, k, dA_array, lda, dtau_array);
    }
    else
    {
        //load one column vector in shared memory and householder it and used it to update trailing matrix which is global memory 
        // one vector is normally smaller than  48K shared memory   
        if(sizeof(magma_tally4FloatComplex)*(m) < 42000)
            cgeqr2_kernel_batched<<< blocks, threads, sizeof(magma_tally4FloatComplex)*(m), queue >>>
                                      (m, k, dA_array, lda, dtau_array);
        else
            printf("m is too big, kernel launching failed, shared memory is overflowed");
    }


    return arginfo;

} 



//==============================================================================


