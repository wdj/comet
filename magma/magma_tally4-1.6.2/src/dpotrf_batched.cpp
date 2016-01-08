/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Azzam Haidar
       @author Tingxing Dong

       @generated from zpotrf_batched.cpp normal z -> d, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_tally4.h"
#include "batched_kernel_param.h"
///////////////////////////////////////////////////////////////////////////////////////
/**
    Purpose
    -------
    DPOTRF computes the Cholesky factorization of a real symmetric
    positive definite matrix dA.

    The factorization has the form
        dA = U**H * U,   if UPLO = Magma_tally4Upper, or
        dA = L  * L**H,  if UPLO = Magma_tally4Lower,
    where U is an upper triangular matrix and L is lower triangular.

    This is the block version of the algorithm, calling Level 3 BLAS.
    If the current stream is NULL, this version replaces it with a new
    stream to overlap computation with communication.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally4_uplo_t
      -     = Magma_tally4Upper:  Upper triangle of dA is stored;
      -     = Magma_tally4Lower:  Lower triangle of dA is stored.

    @param[in]
    n       INTEGER
            The order of the matrix dA.  N >= 0.

    @param[in,out]
    dA      DOUBLE_PRECISION array on the GPU, dimension (LDDA,N)
            On entry, the symmetric matrix dA.  If UPLO = Magma_tally4Upper, the leading
            N-by-N upper triangular part of dA contains the upper
            triangular part of the matrix dA, and the strictly lower
            triangular part of dA is not referenced.  If UPLO = Magma_tally4Lower, the
            leading N-by-N lower triangular part of dA contains the lower
            triangular part of the matrix dA, and the strictly upper
            triangular part of dA is not referenced.
    \n
            On exit, if INFO = 0, the factor U or L from the Cholesky
            factorization dA = U**H * U or dA = L * L**H.

    @param[in]
    ldda     INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,N).
            To benefit from coalescent memory accesses LDDA must be
            divisible by 16.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     > 0:  if INFO = i, the leading minor of order i is not
                  positive definite, and the factorization could not be
                  completed.

    @ingroup magma_tally4_dposv_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_dpotrf_batched(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double **dA_array, magma_tally4_int_t ldda,
    magma_tally4_int_t *info_array,  magma_tally4_int_t batchCount, magma_tally4_queue_t queue)
{
#define A(i_, j_)  (A + (i_) + (j_)*ldda)   
    double d_alpha = -1.0;
    double d_beta  = 1.0;
    cudaMemset(info_array, 0, batchCount*sizeof(magma_tally4_int_t));

    magma_tally4_int_t arginfo = 0;
    if ( uplo != Magma_tally4Upper && uplo != Magma_tally4Lower) {
        arginfo = -1;
    } else if (n < 0) {
        arginfo = -2;
    } else if (ldda < max(1,n)) {
        arginfo = -4;
    }

    if (arginfo != 0) {
        magma_tally4_xerbla( __func__, -(arginfo) );
        return arginfo;
    }

    // Quick return if possible
    if (n == 0) {
        return arginfo;
    }

    if( n > 2048 ){
        printf("=========================================================================================\n");
        printf("   WARNING batched routines are designed for small sizes it might be better to use the\n   Native/Hybrid classical routines if you want performance\n");
        printf("=========================================================================================\n");
    }


    magma_tally4_int_t j, k, ib;
    magma_tally4_int_t nb = POTRF_NB;
    magma_tally4_int_t gemm_crossover = 127;//nb > 32 ? 127 : 160;

#if defined(USE_CUOPT)    
    cublasHandle_t myhandle;
    cublasCreate_v2(&myhandle);
#else
    cublasHandle_t myhandle=NULL;
#endif

    double **dA_displ   = NULL;
    double **dW0_displ  = NULL;
    double **dW1_displ  = NULL;
    double **dW2_displ  = NULL;
    double **dW3_displ  = NULL;
    double **dW4_displ  = NULL;
    double **dinvA_array = NULL;
    double **dwork_array = NULL;

    magma_tally4_malloc((void**)&dA_displ,   batchCount * sizeof(*dA_displ));
    magma_tally4_malloc((void**)&dW0_displ,  batchCount * sizeof(*dW0_displ));
    magma_tally4_malloc((void**)&dW1_displ,  batchCount * sizeof(*dW1_displ));
    magma_tally4_malloc((void**)&dW2_displ,  batchCount * sizeof(*dW2_displ));
    magma_tally4_malloc((void**)&dW3_displ,  batchCount * sizeof(*dW3_displ));
    magma_tally4_malloc((void**)&dW4_displ,  batchCount * sizeof(*dW4_displ));
    magma_tally4_malloc((void**)&dinvA_array, batchCount * sizeof(*dinvA_array));
    magma_tally4_malloc((void**)&dwork_array,    batchCount * sizeof(*dwork_array));

    magma_tally4_int_t invA_msize = ((n+TRI_NB-1)/TRI_NB)*TRI_NB*TRI_NB;
    magma_tally4_int_t dwork_msize = n*nb;
    double* dinvA      = NULL;
    double* dwork      = NULL;// dinvA and dwork are workspace in dtrsm
    double **cpuAarray = NULL;
    magma_tally4_dmalloc( &dinvA, invA_msize * batchCount);
    magma_tally4_dmalloc( &dwork, dwork_msize * batchCount );
    magma_tally4_malloc_cpu((void**) &cpuAarray, batchCount*sizeof(double*));
   /* check allocation */
    if ( dA_displ  == NULL || dW0_displ == NULL || dW1_displ   == NULL || dW2_displ   == NULL || 
         dW3_displ == NULL || dW4_displ == NULL || dinvA_array == NULL || dwork_array == NULL || 
         dinvA     == NULL || dwork     == NULL || cpuAarray   == NULL ) {
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
        magma_tally4_int_t info = MAGMA_tally4_ERR_DEVICE_ALLOC;
        magma_tally4_xerbla( __func__, -(info) );
        return info;
    }

    magma_tally4blas_dlaset_q(Magma_tally4Full, invA_msize, batchCount, MAGMA_tally4_D_ZERO, MAGMA_tally4_D_ZERO, dinvA, invA_msize, queue);
    magma_tally4blas_dlaset_q(Magma_tally4Full, dwork_msize, batchCount, MAGMA_tally4_D_ZERO, MAGMA_tally4_D_ZERO, dwork, dwork_msize, queue);
    dset_pointer(dwork_array, dwork, 1, 0, 0, dwork_msize, batchCount, queue);
    dset_pointer(dinvA_array, dinvA, TRI_NB, 0, 0, invA_msize, batchCount, queue);


    magma_tally4_queue_t cstream;
    magma_tally4blasGetKernelStream(&cstream);
    magma_tally4_int_t streamid;
    const magma_tally4_int_t nbstreams=32;
    magma_tally4_queue_t stream[nbstreams];
    for(k=0; k<nbstreams; k++){
        magma_tally4_queue_create( &stream[k] );
    }
    magma_tally4_getvector( batchCount, sizeof(double*), dA_array, 1, cpuAarray, 1);

    magma_tally4blasSetKernelStream(NULL);

    if (uplo == Magma_tally4Upper) {
        printf("Upper side is unavailable \n");
        goto fin;
    }
    else {
        for(j = 0; j < n; j+=nb) {
            ib = min(nb, n-j);
#if 1
            //===============================================
            //  panel factorization
            //===============================================
            magma_tally4_ddisplace_pointers(dA_displ, dA_array, ldda, j, j, batchCount, queue);
            dset_pointer(dwork_array, dwork, 1, 0, 0, dwork_msize, batchCount, queue);
            dset_pointer(dinvA_array, dinvA, TRI_NB, 0, 0, invA_msize, batchCount, queue);


            #if 0
            arginfo = magma_tally4_dpotrf_panel_batched(
                               uplo, n-j, ib,
                               dA_displ, ldda,
                               dwork_array, dwork_msize,
                               dinvA_array, invA_msize,
                               dW0_displ, dW1_displ, dW2_displ,
                               dW3_displ, dW4_displ,
                               info_array, j, batchCount, myhandle);
            #else
            //arginfo = magma_tally4_dpotrf_rectile_batched(
            arginfo = magma_tally4_dpotrf_recpanel_batched(
                               uplo, n-j, ib, 32,
                               dA_displ, ldda,
                               dwork_array, dwork_msize,
                               dinvA_array, invA_msize,
                               dW0_displ, dW1_displ, dW2_displ,
                               dW3_displ, dW4_displ, 
                               info_array, j, batchCount, myhandle, queue);
            #endif
            if(arginfo != 0 ) goto fin;
            //===============================================
            // end of panel
            //===============================================
#endif            
#if 1
            //real_Double_t gpu_time;
            //gpu_time = magma_tally4_sync_wtime(NULL);
            if( (n-j-ib) > 0){
                if( (n-j-ib) > gemm_crossover)   
                { 
                    //-------------------------------------------
                    //          USE STREAM  HERK
                    //-------------------------------------------
                    // since it use different stream I need to wait the panel.
                    // But since the code use the NULL stream everywhere, 
                    // so I don't need it, because the NULL stream do the sync by itself
                    //magma_tally4_queue_sync(NULL); 
                    /* you must know the matrix layout inorder to do it */  
                    for(k=0; k<batchCount; k++)
                    {
                        streamid = k%nbstreams;                                       
                        magma_tally4blasSetKernelStream(stream[streamid]);
                        // call herk, class dsyrk must call cpu pointer 
                        magma_tally4_dsyrk(Magma_tally4Lower, Magma_tally4NoTrans, n-j-ib, ib, 
                            d_alpha, 
                            (const double*) cpuAarray[k] + j+ib+j*ldda, ldda, 
                            d_beta,
                            cpuAarray[k] + j+ib+(j+ib)*ldda, ldda);

                     }
                     // need to synchronise to be sure that panel do not start before
                     // finishing the update at least of the next panel
                     // BUT no need for it as soon as the other portion of the code 
                     // use the NULL stream which do the sync by itself 
                     //magma_tally4_device_sync(); 
                     magma_tally4blasSetKernelStream(NULL);
                }
                else
                {
                    //-------------------------------------------
                    //          USE BATCHED GEMM(which is a HERK in fact, since it only access the lower part)
                    //-------------------------------------------
                    magma_tally4_ddisplace_pointers(dA_displ, dA_array, ldda, j+ib, j, batchCount, queue);
                    magma_tally4_ddisplace_pointers(dW1_displ, dA_array, ldda, j+ib, j+ib, batchCount, queue);
                    magma_tally4blas_dsyrk_batched(uplo, Magma_tally4NoTrans, n-j-ib, ib,
                                          d_alpha, dA_displ, ldda, 
                                          d_beta,  dW1_displ, ldda, 
                                          batchCount, queue);
                }
            } 
            //gpu_time = magma_tally4_sync_wtime(NULL) - gpu_time;
            //real_Double_t flops = (n-j-ib) * (n-j-ib) * ib / 1e9 * batchCount;
            //real_Double_t gpu_perf = flops / gpu_time;
            //printf("Rows= %d, Colum=%d, herk time = %7.2fms, Gflops= %7.2f\n", n-j-ib, ib, gpu_time*1000, gpu_perf);
#endif
        }
    }

fin:
    magma_tally4_queue_sync(NULL);
    for(k=0; k<nbstreams; k++){
        magma_tally4_queue_destroy( stream[k] );
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

    return arginfo;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


