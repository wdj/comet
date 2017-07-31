/*
   -- MAGMA_tally2 (version 1.6.1) --
   Univ. of Tennessee, Knoxville
   Univ. of California, Berkeley
   Univ. of Colorado, Denver
   @date January 2015

   @author Azzam Haidar
   @author Adrien Remy

   @precisions normal z -> s d c
 */
#include "common_magma_tally2.h"
#include "batched_kernel_param.h"

///////////////////////////////////////////////////////////////////////////////////////
/**
    Purpose
    -------
    ZGETRF computes an LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.

    The factorization has the form
        A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 3 BLAS version of the algorithm.
    
    If the current stream is NULL, this version replaces it with a new
    stream to overlap computation with communication.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in,out]
    dA      COMPLEX_16 array on the GPU, dimension (LDDA,N).
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    @param[in]
    ldda     INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

    @param[out]
    ipiv    INTEGER array, dimension (min(M,N))
            The pivot indices; for 1 <= i <= min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
                  or another error occured, such as memory allocation failed.
      -     > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                  has been completed, but the factor U is exactly
                  singular, and division by zero will occur if it is used
                  to solve a system of equations.

    @ingroup magma_tally2_zgesv_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_zgetrf_nopiv_batched(
        magma_tally2_int_t m, magma_tally2_int_t n,
        magma_tally2DoubleComplex **dA_array, 
        magma_tally2_int_t lda,
        magma_tally2_int_t *info_array, 
        magma_tally2_int_t batchCount, magma_tally2_queue_t queue)
{
#define A(i_, j_)  (A + (i_) + (j_)*lda)   
    magma_tally2_int_t min_mn = min(m, n);
    cudaMemset(info_array, 0, batchCount*sizeof(magma_tally2_int_t));

    /* Check arguments */
    magma_tally2_int_t arginfo = 0;
    if (m < 0)
        arginfo = -1;
    else if (n < 0)
        arginfo = -2;
    else if (lda < max(1,m))
        arginfo = -4;

    if (arginfo != 0) {
        magma_tally2_xerbla( __func__, -(arginfo) );
        return arginfo;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0)
        if(min_mn == 0 ) return arginfo;

    if( m >  2048 || n > 2048 ){
        printf("=========================================================================================\n");
        printf("   WARNING batched routines are designed for small sizes it might be better to use the\n   Native/Hybrid classical routines if you want performance\n");
        printf("=========================================================================================\n");
    }


    magma_tally2DoubleComplex neg_one = MAGMA_tally2_Z_NEG_ONE;
    magma_tally2DoubleComplex one  = MAGMA_tally2_Z_ONE;
    magma_tally2_int_t ib, i, k, pm;
    magma_tally2_int_t nb = BATRF_NB;
    magma_tally2_int_t gemm_crossover = nb > 32 ? 127 : 160;
    // magma_tally2_int_t gemm_crossover = 0;// use only stream gemm

#if defined(USE_CUOPT)    
    cublasHandle_t myhandle;
    cublasCreate_v2(&myhandle);
#else
    cublasHandle_t myhandle=NULL;
#endif
 

    magma_tally2DoubleComplex **dA_displ   = NULL;
    magma_tally2DoubleComplex **dW0_displ  = NULL;
    magma_tally2DoubleComplex **dW1_displ  = NULL;
    magma_tally2DoubleComplex **dW2_displ  = NULL;
    magma_tally2DoubleComplex **dW3_displ  = NULL;
    magma_tally2DoubleComplex **dW4_displ  = NULL;
    magma_tally2DoubleComplex **dinvA_array = NULL;
    magma_tally2DoubleComplex **dwork_array = NULL;

    magma_tally2_malloc((void**)&dA_displ,   batchCount * sizeof(*dA_displ));
    magma_tally2_malloc((void**)&dW0_displ,  batchCount * sizeof(*dW0_displ));
    magma_tally2_malloc((void**)&dW1_displ,  batchCount * sizeof(*dW1_displ));
    magma_tally2_malloc((void**)&dW2_displ,  batchCount * sizeof(*dW2_displ));
    magma_tally2_malloc((void**)&dW3_displ,  batchCount * sizeof(*dW3_displ));
    magma_tally2_malloc((void**)&dW4_displ,  batchCount * sizeof(*dW4_displ));
    magma_tally2_malloc((void**)&dinvA_array, batchCount * sizeof(*dinvA_array));
    magma_tally2_malloc((void**)&dwork_array, batchCount * sizeof(*dwork_array));

    magma_tally2_int_t invA_msize = ((n+TRI_NB-1)/TRI_NB)*TRI_NB*TRI_NB;
    magma_tally2_int_t dwork_msize = max(m,n)*nb;
    magma_tally2DoubleComplex* dinvA      = NULL;
    magma_tally2DoubleComplex* dwork      = NULL;// dinvA and dwork are workspace in ztrsm
    magma_tally2DoubleComplex **cpuAarray = NULL;
    magma_tally2_zmalloc( &dinvA, invA_msize * batchCount);
    magma_tally2_zmalloc( &dwork, dwork_msize * batchCount );
    magma_tally2_malloc_cpu((void**) &cpuAarray, batchCount*sizeof(magma_tally2DoubleComplex*));
   /* check allocation */
    if ( dA_displ  == NULL || dW0_displ == NULL || dW1_displ   == NULL || dW2_displ   == NULL || 
         dW3_displ == NULL || dW4_displ == NULL || dinvA_array == NULL || dwork_array == NULL || 
         dinvA     == NULL || dwork     == NULL || cpuAarray   == NULL ) {
        magma_tally2_free(dA_displ);
        magma_tally2_free(dW0_displ);
        magma_tally2_free(dW1_displ);
        magma_tally2_free(dW2_displ);
        magma_tally2_free(dW3_displ);
        magma_tally2_free(dW4_displ);
        magma_tally2_free(dinvA_array);
        magma_tally2_free(dwork_array);
        magma_tally2_free( dinvA );
        magma_tally2_free( dwork );
        free(cpuAarray);
        magma_tally2_int_t info = MAGMA_tally2_ERR_DEVICE_ALLOC;
        magma_tally2_xerbla( __func__, -(info) );
        return info;
    }

    magma_tally2blas_zlaset_q(Magma_tally2Full, invA_msize, batchCount, MAGMA_tally2_Z_ZERO, MAGMA_tally2_Z_ZERO, dinvA, invA_msize, queue);
    magma_tally2blas_zlaset_q(Magma_tally2Full, dwork_msize, batchCount, MAGMA_tally2_Z_ZERO, MAGMA_tally2_Z_ZERO, dwork, dwork_msize, queue);
    zset_pointer(dwork_array, dwork, n, 0, 0, dwork_msize, batchCount, queue);
    zset_pointer(dinvA_array, dinvA, TRI_NB, 0, 0, invA_msize, batchCount, queue);

    // printf(" I am in zgetrfbatched\n");
    magma_tally2_queue_t cstream;
    magma_tally2blasGetKernelStream(&cstream);
    magma_tally2_int_t streamid;
    const magma_tally2_int_t nbstreams=32;
    magma_tally2_queue_t stream[nbstreams];
    for(i=0; i<nbstreams; i++){
        magma_tally2_queue_create( &stream[i] );
    }
    magma_tally2_getvector( batchCount, sizeof(magma_tally2DoubleComplex*), dA_array, 1, cpuAarray, 1);


    //printf(" I am after malloc\n");



    for(i = 0; i < min_mn; i+=nb) 
    {
        magma_tally2blasSetKernelStream(NULL);

        ib = min(nb, min_mn-i);
        pm = m-i;


        magma_tally2_zdisplace_pointers(dA_displ, dA_array, lda, i, i, batchCount, queue);
        zset_pointer(dwork_array, dwork, nb, 0, 0, dwork_msize, batchCount, queue);
#if 0
        /* buggy: TODO */
        arginfo = magma_tally2_zgetrf_panel_nopiv_batched(
                pm, ib,
                dA_displ, lda,
                dwork_array, nb, 
                dinvA_array, invA_msize, 
                dW0_displ, dW1_displ, dW2_displ, 
                dW3_displ, dW4_displ,
                info_array, i,
                batchCount, myhandle, queue); 
 
#else
        arginfo = magma_tally2_zgetrf_recpanel_nopiv_batched(
                pm, ib, 32,
                dA_displ, lda,
                dwork_array, nb, 
                dinvA_array, invA_msize, 
                dW0_displ, dW1_displ, dW2_displ, 
                dW3_displ, dW4_displ,
                info_array, i,
                batchCount, myhandle, queue);   
#endif

        if(arginfo != 0 ) goto fin;

#define RUN_ALL
#ifdef RUN_ALL

        if( (i + ib) < n)
        {
            // swap right side and trsm     
            //magma_tally2_zdisplace_pointers(dA_displ, dA_array, lda, i, i+ib, batchCount);
            zset_pointer(dwork_array, dwork, nb, 0, 0, dwork_msize, batchCount, queue); // I don't think it is needed Azzam

            magma_tally2_zdisplace_pointers(dA_displ, dA_array, lda, i, i, batchCount, queue);
            magma_tally2_zdisplace_pointers(dW0_displ, dA_array, lda, i, i+ib, batchCount, queue);
            magma_tally2blas_ztrsm_work_batched(Magma_tally2Left, Magma_tally2Lower, Magma_tally2NoTrans, Magma_tally2Unit, 1,
                    ib, n-i-ib,
                    MAGMA_tally2_Z_ONE,
                    dA_displ,    lda, // dA
                    dW0_displ,   lda, // dB
                    dwork_array,  nb, // dX
                    dinvA_array,  invA_msize, 
                    dW1_displ,   dW2_displ, 
                    dW3_displ,   dW4_displ,
                    1, batchCount, queue);

            if( (i + ib) < m)
            {    
                // if gemm size is >160 use a streamed classical cublas gemm since it is faster
                // the batched is faster only when M=N<=160 for K40c
                //-------------------------------------------
                //          USE STREAM  GEMM
                //-------------------------------------------
                if( (m-i-ib) > gemm_crossover  && (n-i-ib) > gemm_crossover)   
                { 
                    //printf("caling streamed dgemm %d %d %d \n", m-i-ib, n-i-ib, ib);

                    // since it use different stream I need to wait the TRSM and swap.
                    // But since the code use the NULL stream everywhere, 
                    // so I don't need it, because the NULL stream do the sync by itself
                    // magma_tally2_queue_sync(NULL);
                    for(k=0; k<batchCount; k++)
                    {
                        streamid = k%nbstreams;                                       
                        magma_tally2blasSetKernelStream(stream[streamid]);
                        magma_tally2_zgemm(Magma_tally2NoTrans, Magma_tally2NoTrans, 
                                m-i-ib, n-i-ib, ib,
                                neg_one, cpuAarray[k] + (i+ib)+i*lda, lda, 
                                         cpuAarray[k] + i+(i+ib)*lda, lda,
                                one,     cpuAarray[k] + (i+ib)+(i+ib)*lda, lda);
                    }
                    // need to synchronise to be sure that zgetf2 do not start before
                    // finishing the update at least of the next panel
                    // BUT no need for it as soon as the other portion of the code 
                    // use the NULL stream which do the sync by itself 
                    //magma_tally2_device_sync(); 
                }
                //-------------------------------------------
                //          USE BATCHED GEMM
                //-------------------------------------------
                else
                {
                    magma_tally2_zdisplace_pointers(dA_displ, dA_array,  lda, i+ib,    i, batchCount, queue);
                    magma_tally2_zdisplace_pointers(dW1_displ, dA_array, lda,    i, i+ib, batchCount, queue);
                    magma_tally2_zdisplace_pointers(dW2_displ, dA_array, lda, i+ib, i+ib, batchCount, queue);
                    //printf("caling batched dgemm %d %d %d \n", m-i-ib, n-i-ib, ib);
                    magma_tally2blas_zgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, m-i-ib, n-i-ib, ib, 
                            neg_one, dA_displ, lda, 
                            dW1_displ, lda, 
                            one,  dW2_displ, lda, 
                            batchCount, queue);
                } // end of batched/stream gemm
            } // end of  if( (i + ib) < m) 
        } // end of if( (i + ib) < n)
#endif
    }// end of for

fin:
    magma_tally2_queue_sync(NULL);

    for(i=0; i<nbstreams; i++){
        magma_tally2_queue_destroy( stream[i] );
    }
    magma_tally2blasSetKernelStream(cstream);


#if defined(USE_CUOPT)    
    cublasDestroy_v2(myhandle);
#endif


    magma_tally2_free(dA_displ);
    magma_tally2_free(dW0_displ);
    magma_tally2_free(dW1_displ);
    magma_tally2_free(dW2_displ);
    magma_tally2_free(dW3_displ);
    magma_tally2_free(dW4_displ);
    magma_tally2_free(dinvA_array);
    magma_tally2_free(dwork_array);
    magma_tally2_free( dinvA );
    magma_tally2_free( dwork );
    free(cpuAarray);

    return arginfo;

}

