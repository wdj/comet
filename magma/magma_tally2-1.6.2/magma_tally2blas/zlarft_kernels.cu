/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @precisions normal z -> s d c
       @author Azzam Haidar
*/

#include "common_magma_tally2.h"
#include "magma_tally2_templates.h"
#define zgemv_bs 32
#define BLOCK_SIZE 512

#define use_gemm_larft

extern __shared__ magma_tally2DoubleComplex shared_data[];


//===================================================================================================
static __device__
void zlarft_gemvcolwise_device( int m, magma_tally2DoubleComplex *v, magma_tally2DoubleComplex *tau,
                         magma_tally2DoubleComplex *c, int ldc, magma_tally2DoubleComplex *T, int ldt, int step )
{

    const int thblk =  blockIdx.x;
    if (thblk > step)
        return;
    /* if blockIdx.x<step step performs the z = V(tx:n,tx)' * V(tx:n,1:tx-1) used for computing T:*/

    if ( !MAGMA_tally2_Z_EQUAL(*tau, MAGMA_tally2_Z_ZERO) ) {
        if(thblk<step){    
            const int tx = threadIdx.x;
            magma_tally2DoubleComplex *dc = c + blockIdx.x * ldc;
           
            __shared__ magma_tally2DoubleComplex sum[ BLOCK_SIZE ];
            magma_tally2DoubleComplex tmp;
           
            /* perform  {T_i}^H := V(:,i)' * V(:,1:i-1)  */
            if (tx==0)
                tmp = dc[0]; //since V[0] should be one
            else
                tmp = MAGMA_tally2_Z_ZERO;
            for( int j = tx+1; j < m; j += BLOCK_SIZE ){
                tmp +=  MAGMA_tally2_Z_CNJG( v[j] ) * dc[j];
            }
            sum[tx] = tmp;
            magma_tally2_sum_reduce< BLOCK_SIZE >( tx, sum );
            #if defined (use_gemm_larft)
            *(T+thblk) = MAGMA_tally2_Z_CNJG(sum[0]);
            #else
            tmp = - MAGMA_tally2_Z_CNJG(*tau) * sum[0]; 
            *(T+thblk) = MAGMA_tally2_Z_CNJG(tmp); // T = - tau(tx) * V(tx:n,1:tx-1)' * V(tx:n,tx) = tmp'
            //*(T+thblk) = - MAGMA_tally2_Z_CNJG(sum[0]) * (*tau); // T = - tau(tx) * V(tx:n,1:tx-1)' * V(tx:n,tx) = tmp'
            #endif
        }
        else{
            #if defined (use_gemm_larft)
            *(T+thblk) = MAGMA_tally2_Z_ONE;
            #else
            *(T+thblk) = *tau;
            #endif
        }
    }// in case tau is zero put the corresponding column of T to zero
    else 
    {
        *(T+thblk) = MAGMA_tally2_Z_ZERO;
    }
}
//===================================================================================================
__global__
void zlarft_gemvcolwise_kernel( int m, magma_tally2DoubleComplex *v, int ldv, magma_tally2DoubleComplex *tau,
                          magma_tally2DoubleComplex *T, int ldt, int step )
{
    zlarft_gemvcolwise_device(m, v+step+step*ldv, tau+step, v+step, ldv, T+step*ldt, ldt, step);
}
//===================================================================================================
__global__
void zlarft_gemvcolwise_kernel_batched( int m, magma_tally2DoubleComplex **v_array, int ldv, magma_tally2DoubleComplex **tau_array,
                          magma_tally2DoubleComplex **T_array, int ldt, int step )
{
    int batchid = blockIdx.z;
    zlarft_gemvcolwise_device(m, v_array[batchid]+step+step*ldv, tau_array[batchid]+step, v_array[batchid]+step, ldv, T_array[batchid]+step*ldt, ldt, step);
}
//===================================================================================================
extern "C" 
void magma_tally2blas_zlarft_gemvcolwise(
    magma_tally2_int_t m,  magma_tally2_int_t step,
    magma_tally2DoubleComplex *v, magma_tally2_int_t ldv, 
    magma_tally2DoubleComplex *T,  magma_tally2_int_t ldt,
    magma_tally2DoubleComplex *tau)
{
    dim3 grid( step+1, 1, 1 );
    dim3 threads( BLOCK_SIZE );
    zlarft_gemvcolwise_kernel<<< grid, threads, 0, magma_tally2_stream >>>( m, v, ldv, tau, T, ldt, step);

}
//===================================================================================================
extern "C" 
void magma_tally2blas_zlarft_gemvcolwise_batched(
    magma_tally2_int_t m,  magma_tally2_int_t step,
    magma_tally2DoubleComplex **v_array, magma_tally2_int_t ldv, 
    magma_tally2DoubleComplex **T_array,  magma_tally2_int_t ldt,
    magma_tally2DoubleComplex **tau_array, magma_tally2_int_t batchCount, magma_tally2_queue_t queue )
{
    dim3 grid( step+1, 1, batchCount );
    dim3 threads( BLOCK_SIZE );
    zlarft_gemvcolwise_kernel_batched<<< grid, threads, 0, queue >>>( m, v_array, ldv, tau_array, T_array, ldt, step);

}
//===================================================================================================




//===================================================================================================
// zgemv(y=alpha*A*x) interface: T/W=tau*v*x, 
static __device__ void
zlarft_gemvrowwise_device(
    int m, int i,
    magma_tally2DoubleComplex *tau, 
    magma_tally2DoubleComplex *v_ptr, int ldv, 
    magma_tally2DoubleComplex *x_ptr, int incx,
    magma_tally2DoubleComplex *T_ptr, int ldt,
    magma_tally2DoubleComplex *W, magma_tally2DoubleComplex* sdata)
{
    int tx = threadIdx.x; 
    int ty = threadIdx.y; 


    if(tx ==0 && ty == 0)
    {
        T_ptr[0] = *tau;
    } 

    if(i <= 0) return;
    
    magma_tally2DoubleComplex res = MAGMA_tally2_Z_ZERO;

    v_ptr += ldv * ty;
            

   
    if(tx < zgemv_bs)
    {
        for(int s=tx; s<m; s+= zgemv_bs)
        {
            res += MAGMA_tally2_Z_CNJG (v_ptr[s]) * x_ptr[s*incx];
        }
    
        sdata[ty * zgemv_bs + tx] = res;
    }
    __syncthreads();

    magma_tally2_sum_reduce<zgemv_bs>(tx, &(sdata[ty*zgemv_bs+0]));

    #if defined (use_gemm_larft)
    if(tx == 0)
    {
            W[ty] = -sdata[ty * zgemv_bs + 0];
    } 
    #else
    if(tx == 0)
    {
            W[ty] = -sdata[ty * zgemv_bs + 0] * (*tau) ;
    }
    #endif 
}




//T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
//T(i,i) = tau(i)
//===================================================================================================
 __global__ void
zlarft_gemvrowwise_kernel(
    int m, int i, 
    magma_tally2DoubleComplex *tau, 
    magma_tally2DoubleComplex *v, int ldv, 
    magma_tally2DoubleComplex *T, int ldt)
{

    magma_tally2DoubleComplex *W =  T +i*ldt;

    magma_tally2DoubleComplex *sdata = (magma_tally2DoubleComplex*)shared_data;

    zlarft_gemvrowwise_device(m, i, tau+i, v+i, ldv,  v+i+i*ldv, 1,  
                           T+i+i*ldt , ldt, W, sdata);
}

//===================================================================================================
__global__ void
zlarft_gemvrowwise_kernel_batched(
    int m, int i,
    magma_tally2DoubleComplex **tau_array, 
    magma_tally2DoubleComplex **v_array, int ldv, 
    magma_tally2DoubleComplex **T_array, int ldt)
{

    int batchid = blockIdx.z;

    magma_tally2DoubleComplex *W =  T_array[batchid] +i*ldt;

    magma_tally2DoubleComplex *sdata = (magma_tally2DoubleComplex*)shared_data;

    zlarft_gemvrowwise_device(m, i, tau_array[batchid]+i, v_array[batchid]+i, ldv,  v_array[batchid]+i+i*ldv, 1,  
                           T_array[batchid] +i+i*ldt , ldt, W, sdata);
}

//===================================================================================================
extern "C"
void magma_tally2blas_zlarft_gemvrowwise(
    magma_tally2_int_t m, magma_tally2_int_t i, 
    magma_tally2DoubleComplex *tau, 
    magma_tally2DoubleComplex *v, magma_tally2_int_t ldv, 
    magma_tally2DoubleComplex *T, magma_tally2_int_t ldt,
    magma_tally2DoubleComplex *W)
{

    dim3 grid(1);


    dim3 threads(zgemv_bs, max(i,1), 1);


    zlarft_gemvrowwise_kernel <<< grid, threads, sizeof(magma_tally2DoubleComplex)*zgemv_bs*(i+1), magma_tally2_stream>>>(m, i, tau, v, ldv, T, ldt);
}
//===================================================================================================
extern "C"
void magma_tally2blas_zlarft_gemvrowwise_batched(
    magma_tally2_int_t m, magma_tally2_int_t i, 
    magma_tally2DoubleComplex **tau_array, 
    magma_tally2DoubleComplex **v_array, magma_tally2_int_t ldv, 
    magma_tally2DoubleComplex **T_array, magma_tally2_int_t ldt,
    magma_tally2_int_t batchCount, magma_tally2_queue_t queue)
{

    dim3 grid(1, 1, batchCount);
    dim3 threads(zgemv_bs, max(i,1), 1);

    /*  zgemvrowwise used a bigger shared memory and has more data reuse and performs better
    */
    zlarft_gemvrowwise_kernel_batched <<< grid, threads, sizeof(magma_tally2DoubleComplex)*zgemv_bs*(i+1), queue>>>(m, i,  tau_array, v_array, ldv, T_array, ldt);
}
//===================================================================================================
   


//===================================================================================================
/*
   loop_inside
*/
static __device__ void
zlarft_gemv_loop_inside_device(
    int n, int k, 
    magma_tally2DoubleComplex *tau, 
    magma_tally2DoubleComplex *v, int ldv, 
    magma_tally2DoubleComplex *T, int ldt)
{
    int tx = threadIdx.x; 
    int ty = threadIdx.y; 
    
    int incx = 1;
    magma_tally2DoubleComplex *sdata = (magma_tally2DoubleComplex*)shared_data;

    magma_tally2DoubleComplex res;

    // write the first elment
    if(tx ==0 && ty == 0)
    {
        T[0] = tau[0];
    } 
 
    for(int i=1; i<k;i++)
    {

        int m = n-i; 

        magma_tally2DoubleComplex *v_ptr = v;

        v_ptr += i;

        magma_tally2DoubleComplex *x_ptr = v_ptr + i * ldv;
            
        res = MAGMA_tally2_Z_ZERO;
            
        if(tx < zgemv_bs && ty < i)
        {
            v_ptr += ldv * ty;

            for(int s=tx; s<m; s+= zgemv_bs)
            {
                res += MAGMA_tally2_Z_CNJG (v_ptr[s]) * x_ptr[s*incx];
            }
    
            sdata[ty * zgemv_bs + tx] = res;
        }
        __syncthreads();

        magma_tally2_sum_reduce<zgemv_bs>(tx, &(sdata[ty*zgemv_bs+0]));
        

       __syncthreads();
       #if defined (use_gemm_larft)
       if(tx < i && ty == 0)
       {
            T[i* ldt + tx] = sdata[tx * zgemv_bs + 0];  
       } 
       // not needed since it is overwritten in trmv
       /*
       if(tx == i && ty == 0)
       {
           T[i * ldt + i] = tau[i];
       }
       */
       #else
       if(tx < i && ty == 0)
       {
           T[i* ldt + tx] = -sdata[tx * zgemv_bs + 0] * (tau[i]) ;  
       } 
      
       if(tx == i && ty == 0)
       {
           T[i * ldt + i] = tau[i];
       }
       #endif
     
       v_ptr -= i;

    }// end of loop k
}
//===================================================================================================
__global__ void
zlarft_gemv_loop_inside_kernel(
    int n, int k, 
    magma_tally2DoubleComplex *tau, 
    magma_tally2DoubleComplex *v, int ldv, 
    magma_tally2DoubleComplex *T, int ldt)
{
    zlarft_gemv_loop_inside_device(n, k, tau, v, ldv, T, ldt);
}
//===================================================================================================
__global__ void
zlarft_gemv_loop_inside_kernel_batched(
    int n, int k, 
    magma_tally2DoubleComplex **tau_array, 
    magma_tally2DoubleComplex **v_array, int ldv, 
    magma_tally2DoubleComplex **T_array, int ldt)
{
    int batchid = blockIdx.z;
    zlarft_gemv_loop_inside_device(n, k, tau_array[batchid], v_array[batchid], ldv, T_array[batchid], ldt);
}
//===================================================================================================
//===================================================================================================
//===================================================================================================
extern "C"
void magma_tally2blas_zlarft_gemv_loop_inside(
    int n, int k, 
    magma_tally2DoubleComplex *tau, 
    magma_tally2DoubleComplex *v, int ldv, 
    magma_tally2DoubleComplex *T, int ldt)
{

    dim3 grid(1);
    dim3 threads(zgemv_bs, max(k,1), 1);
    zlarft_gemv_loop_inside_kernel<<<grid, threads, sizeof(magma_tally2DoubleComplex) * (zgemv_bs*(k+1)), magma_tally2_stream>>>(n, k, tau, v, ldv, T, ldt); 
}
//===================================================================================================
extern "C"
void magma_tally2blas_zlarft_gemv_loop_inside_batched(
    int n, int k, 
    magma_tally2DoubleComplex **tau_array, 
    magma_tally2DoubleComplex **v_array, int ldv, 
    magma_tally2DoubleComplex **T_array, int ldt, magma_tally2_int_t batchCount, magma_tally2_queue_t queue)
{

    dim3 grid(1, 1, batchCount);
    dim3 threads(zgemv_bs, max(k,1), 1);
    zlarft_gemv_loop_inside_kernel_batched<<<grid, threads, sizeof(magma_tally2DoubleComplex) * (zgemv_bs*(k+1)), queue>>>(n, k, tau_array, v_array, ldv, T_array, ldt); 
}
//===================================================================================================





//===================================================================================================
static  __device__ void 
zlarft_ztrmv_sm32x32_device(
    int n, int k, magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *Tin, int ldtin,  magma_tally2DoubleComplex *Tout, int ldtout )
{
    int tx = threadIdx.x; 
    magma_tally2DoubleComplex *sdata = (magma_tally2DoubleComplex*)shared_data;
    magma_tally2DoubleComplex res;

    // this routine apply a sequence of trmv to update k column of the triangular
    // T starting at n-k to n where T is of size n by n and where the first n-k 
    // columns of T are supposed updated previously.
    // So the routine load all of T nxn to the shared memory 
    // and apply the sequence of trmv.
    // to update a certain column i, threads go in horizontal fashion where
    // every thread read one row and do it gemv(dot) to generate 
    // one element of the column of T then move to the next column

    // read T into shared
    for(int s=0; s<n-k; s++)
    {
        sdata[tx + s*n] = Tin[tx + s * ldtin];
    }
    
#if defined(use_gemm_larft)
    for(int s=n-k; s<n; s++)
    {
        if(tx == s)
            sdata[tx + s*n] = tau[s];
        else
            sdata[tx + s*n] = -tau[s] * Tin[tx + s * ldtin];
    }
#else
    for(int s=n-k; s<n; s++)
    {
        sdata[tx + s*n] = Tin[tx + s * ldtin];
    }
#endif

    // perform trmv
    for(int i=n-k; i<n;i++)
    {
       __syncthreads();  
       res = MAGMA_tally2_Z_ZERO;
       if(tx < i)
       {
           for(int j=tx; j<i; j++)
           {
               res += sdata[tx + j * n] * sdata[j+ i * n];      
           }
       }       
       __syncthreads();  
       if(tx < i)
       {
           sdata[tx + i * n] = res;
       }
    } 

    __syncthreads();  
    // write back the updated block of k column of T
    for(int s=n-k; s<n; s++)
    {
       Tout[tx + s * ldtout] = sdata[tx + s*n];
    }

}
//===================================================================================================
__global__ void 
zlarft_ztrmv_sm32x32_kernel(
    int n, int k, magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *Tin, int ldtin,  magma_tally2DoubleComplex *Tout, int ldtout )
{
    zlarft_ztrmv_sm32x32_device( n, k, tau, Tin, ldtin, Tout, ldtout);
}
//===================================================================================================
__global__ void 
zlarft_ztrmv_sm32x32_kernel_batched(
    int n, int k, magma_tally2DoubleComplex **tau_array,
    magma_tally2DoubleComplex **Tin_array, int ldtin,  magma_tally2DoubleComplex **Tout_array, int ldtout )
{
    int batchId = blockIdx.z;
    zlarft_ztrmv_sm32x32_device( n, k, tau_array[batchId], Tin_array[batchId], ldtin, Tout_array[batchId], ldtout);
}
//===================================================================================================
//===================================================================================================
extern "C"
void magma_tally2blas_zlarft_ztrmv_sm32x32(
    magma_tally2_int_t m, magma_tally2_int_t n, 
    magma_tally2DoubleComplex *tau, 
    magma_tally2DoubleComplex *Tin, magma_tally2_int_t ldtin, 
    magma_tally2DoubleComplex *Tout, magma_tally2_int_t ldtout)
{

    dim3 grid(1);
    dim3 threads(max(m,1), 1, 1);
    zlarft_ztrmv_sm32x32_kernel <<< grid, threads, sizeof(magma_tally2DoubleComplex)*(m*m), magma_tally2_stream >>> (m, n,  tau, Tin, ldtin, Tout, ldtout);
}
//===================================================================================================
extern "C"
void magma_tally2blas_zlarft_ztrmv_sm32x32_batched(
    magma_tally2_int_t m, magma_tally2_int_t n, 
    magma_tally2DoubleComplex **tau_array, 
    magma_tally2DoubleComplex **Tin_array, magma_tally2_int_t ldtin, 
    magma_tally2DoubleComplex **Tout_array, magma_tally2_int_t ldtout,
    magma_tally2_int_t batchCount, magma_tally2_queue_t queue)
{

    dim3 grid(1, 1, batchCount);
    dim3 threads(max(m,1), 1, 1);
    zlarft_ztrmv_sm32x32_kernel_batched <<< grid, threads, sizeof(magma_tally2DoubleComplex)*(m*m), queue >>> (m, n,  tau_array, Tin_array, ldtin, Tout_array, ldtout);
}
//===================================================================================================




//===================================================================================================
//===================================================================================================
static __device__ void 
zlarft_recztrmv_sm32x32_device(
    int m, int n, magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *Trec, int ldtrec, magma_tally2DoubleComplex *Ttri, int ldttri)
{
    int tx = threadIdx.x; 
    magma_tally2DoubleComplex *sdata = (magma_tally2DoubleComplex*)shared_data;
    magma_tally2DoubleComplex res;

    // to update a certain column i, threads go in horizontal fashion where
    // every thread read one row and do it gemv(dot) to generate 
    // one element of the column of T then move to the next column

    // read T into shared
    for(int s=0; s<n; s++)
    {
        sdata[tx + s*n] = Trec[tx + s * ldtrec];
    }
    __syncthreads();  
    
    // perform sequence of n-1 gemv
    for(int i=0; i<n;i++)
    {
       res = MAGMA_tally2_Z_ZERO;
       for(int j=0; j<i; j++)
       {
           res += sdata[tx + j * n] * Ttri[j+ i * ldttri];      
       }
       __syncthreads();   // a enlever
       sdata[tx + i * n] = -tau[i] * (sdata[tx + i * n] + res);
       __syncthreads();  
    } 

    // write back the updated block of k column of T  multiplying by -tau
    for(int s=0; s<n; s++)
    {
       Trec[tx + s * ldtrec] = sdata[tx + s*n];
    }

}

//===================================================================================================
__global__ void 
zlarft_recztrmv_sm32x32_kernel(
    int m, int n, magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *Trec, int ldtrec, magma_tally2DoubleComplex *Ttri, int ldttri)
{
    zlarft_recztrmv_sm32x32_device(m, n, tau, Trec, ldtrec, Ttri, ldttri);
}
//===================================================================================================
__global__ void 
zlarft_recztrmv_sm32x32_kernel_batched(
    int m, int n, magma_tally2DoubleComplex **tau_array,
    magma_tally2DoubleComplex **Trec_array, int ldtrec, magma_tally2DoubleComplex **Ttri_array, int ldttri)
{
    int batchId = blockIdx.z;
    zlarft_recztrmv_sm32x32_device(m, n, tau_array[batchId], Trec_array[batchId], ldtrec, Ttri_array[batchId], ldttri);
}
//===================================================================================================
extern "C"
void magma_tally2blas_zlarft_recztrmv_sm32x32(
    magma_tally2_int_t m, magma_tally2_int_t n, 
    magma_tally2DoubleComplex *tau, 
    magma_tally2DoubleComplex *Trec, magma_tally2_int_t ldtrec, 
    magma_tally2DoubleComplex *Ttri, magma_tally2_int_t ldttri)
{

    dim3 grid(1);
    dim3 threads(max(m,1), 1, 1);
    zlarft_recztrmv_sm32x32_kernel <<< grid, threads, sizeof(magma_tally2DoubleComplex)*(m*n), magma_tally2_stream >>> (m, n,  tau, Trec, ldtrec, Ttri, ldttri);
}
//===================================================================================================
extern "C"
void magma_tally2blas_zlarft_recztrmv_sm32x32_batched(
    magma_tally2_int_t m, magma_tally2_int_t n, 
    magma_tally2DoubleComplex **tau_array, 
    magma_tally2DoubleComplex **Trec_array, magma_tally2_int_t ldtrec, 
    magma_tally2DoubleComplex **Ttri_array, magma_tally2_int_t ldttri,
    magma_tally2_int_t batchCount, magma_tally2_queue_t queue)
{

    dim3 grid(1, 1, batchCount);
    dim3 threads(max(m,1), 1, 1);
    zlarft_recztrmv_sm32x32_kernel_batched <<< grid, threads, sizeof(magma_tally2DoubleComplex)*(m*n), queue >>> (m, n,  tau_array, Trec_array, ldtrec, Ttri_array, ldttri);
}
//===================================================================================================


