/*
   -- MAGMA_tally4 (version 1.6.1) --
   Univ. of Tennessee, Knoxville
   Univ. of California, Berkeley
   Univ. of Colorado, Denver
   @date January 2015

   @author Azzam Haidar
   @author Tingxing Dong

   @generated from zgetrf_batched.cpp normal z -> c, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_tally4.h"
#include "batched_kernel_param.h"
///////////////////////////////////////////////////////////////////////////////////////
/**
    Purpose
    -------
    CGETRF computes an LU factorization of a general M-by-N matrix A
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
    dA      COMPLEX array on the GPU, dimension (LDDA,N).
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

    @ingroup magma_tally4_cgesv_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_cgetrf_batched(
        magma_tally4_int_t m, magma_tally4_int_t n,
        magma_tally4FloatComplex **dA_array, 
        magma_tally4_int_t ldda,
        magma_tally4_int_t **ipiv_array, 
        magma_tally4_int_t *info_array, 
        magma_tally4_int_t batchCount, magma_tally4_queue_t queue)
{
#define A(i_, j_)  (A + (i_) + (j_)*ldda)   

    magma_tally4_int_t min_mn = min(m, n);
    cudaMemset(info_array, 0, batchCount*sizeof(magma_tally4_int_t));

    /* Check arguments */
    magma_tally4_int_t arginfo = 0;
    if (m < 0)
        arginfo = -1;
    else if (n < 0)
        arginfo = -2;
    else if (ldda < max(1,m))
        arginfo = -4;

    if (arginfo != 0) {
        magma_tally4_xerbla( __func__, -(arginfo) );
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


//#define ENABLE_TIMER3

#if defined(ENABLE_TIMER3)
    real_Double_t   tall=0.0, tloop=0., talloc=0., tdalloc=0.;
    tall   = magma_tally4_sync_wtime(0);
    talloc = magma_tally4_sync_wtime(0);
#endif

    magma_tally4FloatComplex neg_one = MAGMA_tally4_C_NEG_ONE;
    magma_tally4FloatComplex one  = MAGMA_tally4_C_ONE;
    magma_tally4_int_t ib, i, k, pm;
    magma_tally4_int_t nb = BATRF_NB;
    magma_tally4_int_t gemm_crossover = nb > 32 ? 127 : 160;
    // magma_tally4_int_t gemm_crossover = n;// use only stream gemm

#if defined(USE_CUOPT)    
    cublasHandle_t myhandle;
    cublasCreate_v2(&myhandle);
#else
    cublasHandle_t myhandle=NULL;
#endif

    magma_tally4_int_t     **dipiv_displ   = NULL;
    magma_tally4FloatComplex **dA_displ   = NULL;
    magma_tally4FloatComplex **dW0_displ  = NULL;
    magma_tally4FloatComplex **dW1_displ  = NULL;
    magma_tally4FloatComplex **dW2_displ  = NULL;
    magma_tally4FloatComplex **dW3_displ  = NULL;
    magma_tally4FloatComplex **dW4_displ  = NULL;
    magma_tally4FloatComplex **dinvA_array = NULL;
    magma_tally4FloatComplex **dwork_array = NULL;


    magma_tally4_malloc((void**)&dipiv_displ,   batchCount * sizeof(*dipiv_displ));
    magma_tally4_malloc((void**)&dA_displ,   batchCount * sizeof(*dA_displ));
    magma_tally4_malloc((void**)&dW0_displ,  batchCount * sizeof(*dW0_displ));
    magma_tally4_malloc((void**)&dW1_displ,  batchCount * sizeof(*dW1_displ));
    magma_tally4_malloc((void**)&dW2_displ,  batchCount * sizeof(*dW2_displ));
    magma_tally4_malloc((void**)&dW3_displ,  batchCount * sizeof(*dW3_displ));
    magma_tally4_malloc((void**)&dW4_displ,  batchCount * sizeof(*dW4_displ));
    magma_tally4_malloc((void**)&dinvA_array, batchCount * sizeof(*dinvA_array));
    magma_tally4_malloc((void**)&dwork_array, batchCount * sizeof(*dwork_array));


    magma_tally4_int_t invA_msize = ((n+TRI_NB-1)/TRI_NB)*TRI_NB*TRI_NB;
    magma_tally4_int_t dwork_msize = n*nb;
    magma_tally4_int_t **pivinfo_array    = NULL;
    magma_tally4_int_t *pivinfo           = NULL; 
    magma_tally4FloatComplex* dinvA      = NULL;
    magma_tally4FloatComplex* dwork      = NULL;// dinvA and dwork are workspace in ctrsm
    magma_tally4FloatComplex **cpuAarray = NULL;
    magma_tally4_cmalloc( &dinvA, invA_msize * batchCount);
    magma_tally4_cmalloc( &dwork, dwork_msize * batchCount );
    magma_tally4_malloc((void**)&pivinfo_array, batchCount * sizeof(*pivinfo_array));
    magma_tally4_malloc((void**)&pivinfo, batchCount * m * sizeof(magma_tally4_int_t));
    magma_tally4_malloc_cpu((void**) &cpuAarray, batchCount*sizeof(magma_tally4FloatComplex*));

   /* check allocation */
    if ( dA_displ  == NULL || dW0_displ == NULL || dW1_displ   == NULL || dW2_displ   == NULL || 
         dW3_displ == NULL || dW4_displ == NULL || dinvA_array == NULL || dwork_array == NULL || 
         dinvA     == NULL || dwork     == NULL || cpuAarray   == NULL || 
         dipiv_displ == NULL || pivinfo_array == NULL || pivinfo == NULL) {
        magma_tally4_free(dA_displ);
        magma_tally4_free(dW0_displ);
        magma_tally4_free(dW1_displ);
        magma_tally4_free(dW2_displ);
        magma_tally4_free(dW3_displ);
        magma_tally4_free(dW4_displ);
        magma_tally4_free(dinvA_array);
        magma_tally4_free(dwork_array);
        magma_tally4_free( dinvA );
        magma_tally4_free( dwork );
        free(cpuAarray);
        magma_tally4_free(dipiv_displ);
        magma_tally4_free(pivinfo_array);
        magma_tally4_free(pivinfo);
        magma_tally4_int_t info = MAGMA_tally4_ERR_DEVICE_ALLOC;
        magma_tally4_xerbla( __func__, -(info) );
        return info;
    }


    magma_tally4blas_claset_q(Magma_tally4Full, invA_msize, batchCount, MAGMA_tally4_C_ZERO, MAGMA_tally4_C_ZERO, dinvA, invA_msize, queue);
    magma_tally4blas_claset_q(Magma_tally4Full, dwork_msize, batchCount, MAGMA_tally4_C_ZERO, MAGMA_tally4_C_ZERO, dwork, dwork_msize, queue);
    cset_pointer(dwork_array, dwork, 1, 0, 0, dwork_msize, batchCount, queue);
    cset_pointer(dinvA_array, dinvA, TRI_NB, 0, 0, invA_msize, batchCount, queue);
    set_ipointer(pivinfo_array, pivinfo, 1, 0, 0, m, batchCount, queue);


    // printf(" I am in cgetrfbatched\n");
    magma_tally4_queue_t cstream;
    magma_tally4blasGetKernelStream(&cstream);
    magma_tally4_int_t streamid;
    const magma_tally4_int_t nbstreams=32;
    magma_tally4_queue_t stream[nbstreams];
    for(i=0; i<nbstreams; i++){
        magma_tally4_queue_create( &stream[i] );
    }
    magma_tally4_getvector( batchCount, sizeof(magma_tally4FloatComplex*), dA_array, 1, cpuAarray, 1);



#if defined(ENABLE_TIMER3)
    printf(" I am after malloc\n");
    talloc = magma_tally4_sync_wtime(0) - talloc;
    tloop  = magma_tally4_sync_wtime(0);
#endif


    for(i = 0; i < min_mn; i+=nb) 
    {
        magma_tally4blasSetKernelStream(NULL);

        ib = min(nb, min_mn-i);
        pm = m-i;
        magma_tally4_idisplace_pointers(dipiv_displ, ipiv_array, ldda, i, 0, batchCount, queue);
        magma_tally4_cdisplace_pointers(dA_displ, dA_array, ldda, i, i, batchCount, queue);
        //===============================================
        //  panel factorization
        //===============================================
#if 0
        arginfo = magma_tally4_cgetf2_batched(
                pm, ib,
                dA_displ, ldda,
                dW1_displ, dW2_displ, dW3_displ,
                dipiv_displ, 
                info_array, i, batchCount, myhandle);   
#else
        arginfo = magma_tally4_cgetrf_recpanel_batched(
                pm, ib, 16,
                dA_displ, ldda,
                dipiv_displ, pivinfo_array,
                dwork_array, nb, 
                dinvA_array, invA_msize, 
                dW0_displ, dW1_displ, dW2_displ, 
                dW3_displ, dW4_displ,
                info_array, i, 
                batchCount, myhandle, queue);   
#endif
        if(arginfo != 0 ) goto fin;
        //===============================================
        // end of panel
        //===============================================

#define RUN_ALL
#ifdef RUN_ALL
        // setup pivinfo before adjusting ipiv
        setup_pivinfo_batched(pivinfo_array, dipiv_displ, pm, ib, batchCount, queue);
        adjust_ipiv_batched(dipiv_displ, ib, i, batchCount, queue);

        // stepinit_ipiv(pivinfo_array, pm, batchCount);// for debug and check swap, it create an ipiv


#if 0
        claswp_batched( i, dA_displ, ldda,
                i, i+ib,
                dipiv_displ, pivinfo_array, batchCount);
#else
        magma_tally4_cdisplace_pointers(dA_displ, dA_array, ldda, i, 0, batchCount, queue);
        magma_tally4_cdisplace_pointers(dW0_displ, dA_array, ldda, i, 0, batchCount, queue);
        magma_tally4_claswp_rowparallel_batched( i, dA_displ, ldda,
                dW0_displ, ldda,
                i, i+ib,
                pivinfo_array, batchCount, queue);

#endif

        if( (i + ib) < n)
        {
            // swap right side and trsm     
            magma_tally4_cdisplace_pointers(dA_displ, dA_array, ldda, i, i+ib, batchCount, queue);
            cset_pointer(dwork_array, dwork, nb, 0, 0, dwork_msize, batchCount, queue); // I don't think it is needed Azzam
            magma_tally4_claswp_rowparallel_batched( n-(i+ib), dA_displ, ldda,
                    dwork_array, nb,
                    i, i+ib,
                    pivinfo_array, batchCount, queue);


            magma_tally4_cdisplace_pointers(dA_displ, dA_array, ldda, i, i, batchCount, queue);
            magma_tally4_cdisplace_pointers(dW0_displ, dA_array, ldda, i, i+ib, batchCount, queue);
            magma_tally4blas_ctrsm_outofplace_batched(Magma_tally4Left, Magma_tally4Lower, Magma_tally4NoTrans, Magma_tally4Unit, 1,
                    ib, n-i-ib,
                    MAGMA_tally4_C_ONE,
                    dA_displ,    ldda, // dA
                    dwork_array,  nb, // dB
                    dW0_displ,   ldda, // dX
                    dinvA_array,  invA_msize, 
                    dW1_displ,   dW2_displ, 
                    dW3_displ,   dW4_displ,
                    0, batchCount, queue);


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
                    //magma_tally4_queue_sync(NULL); 
                    //
                    for(k=0; k<batchCount; k++)
                    {
                        streamid = k%nbstreams;                                       
                        magma_tally4blasSetKernelStream(stream[streamid]);
                        magma_tally4_cgemm(Magma_tally4NoTrans, Magma_tally4NoTrans, 
                                m-i-ib, n-i-ib, ib,
                                neg_one, cpuAarray[k] + (i+ib)+i*ldda, ldda, 
                                         cpuAarray[k] + i+(i+ib)*ldda, ldda,
                                one,     cpuAarray[k] + (i+ib)+(i+ib)*ldda, ldda);
                    }
                    // need to synchronise to be sure that cgetf2 do not start before
                    // finishing the update at least of the next panel
                    // BUT no need for it as soon as the other portion of the code 
                    // use the NULL stream which do the sync by itself 
                    //magma_tally4_device_sync(); 
                }
                //-------------------------------------------
                //          USE BATCHED GEMM
                //-------------------------------------------
                else
                {
                    magma_tally4_cdisplace_pointers(dA_displ, dA_array,  ldda, i+ib,    i, batchCount, queue);
                    magma_tally4_cdisplace_pointers(dW1_displ, dA_array, ldda,    i, i+ib, batchCount, queue);
                    magma_tally4_cdisplace_pointers(dW2_displ, dA_array, ldda, i+ib, i+ib, batchCount, queue);
                    //printf("caling batched dgemm %d %d %d \n", m-i-ib, n-i-ib, ib);
                    magma_tally4blas_cgemm_batched( Magma_tally4NoTrans, Magma_tally4NoTrans, m-i-ib, n-i-ib, ib, 
                            neg_one, dA_displ, ldda, 
                            dW1_displ, ldda, 
                            one,  dW2_displ, ldda, 
                            batchCount, queue);
                } // end of batched/stream gemm
            } // end of  if( (i + ib) < m) 
        } // end of if( (i + ib) < n)
#endif
    }// end of for

fin:
    magma_tally4_queue_sync(NULL);

#if defined(ENABLE_TIMER3)
    tloop   = magma_tally4_sync_wtime(0) - tloop;
    tdalloc = magma_tally4_sync_wtime(0);

#endif

    for(i=0; i<nbstreams; i++){
        magma_tally4_queue_destroy( stream[i] );
    }
    magma_tally4blasSetKernelStream(cstream);


#if defined(USE_CUOPT)    
    cublasDestroy_v2(myhandle);
#endif

    magma_tally4_free(dA_displ);
    magma_tally4_free(dW0_displ);
    magma_tally4_free(dW1_displ);
    magma_tally4_free(dW2_displ);
    magma_tally4_free(dW3_displ);
    magma_tally4_free(dW4_displ);
    magma_tally4_free(dinvA_array);
    magma_tally4_free(dwork_array);
    magma_tally4_free( dinvA );
    magma_tally4_free( dwork );
    free(cpuAarray);
    magma_tally4_free(dipiv_displ);
    magma_tally4_free(pivinfo_array);
    magma_tally4_free(pivinfo);

#if defined(ENABLE_TIMER3)
    tdalloc = magma_tally4_sync_wtime(0) - tdalloc;
    tall = magma_tally4_sync_wtime(0) - tall;
    printf("here is the timing from inside cgetrf_batched talloc: %10.5f  tloop: %10.5f tdalloc: %10.5f tall: %10.5f sum: %10.5f\n", talloc, tloop, tdalloc, tall, talloc+tloop+tdalloc );
#endif
    
    return arginfo;

}





