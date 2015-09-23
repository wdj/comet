/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Azzam Haidar
       @author Ichi Yamazaki

       @generated from zherk_mgpu.cpp normal z -> d, Fri Jan 30 19:00:10 2015

*/
#include "common_magma_minproduct.h"
#include "trace.h"

/**
    Purpose
    -------
    This dsyrk_mgpu is internal routine used by dpotrf_mgpu_right.
    it has specific assumption on the block diagonal.
    
    @ingroup magma_minproduct_dblas3
    ********************************************************************/

extern "C" void
magma_minproduct_dsyrk_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t nb, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    double beta,
    magma_minproductDouble_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t nqueue, magma_minproduct_queue_t queues[][10])
{
#define dB(id, i, j)  (dB[(id)]+(j)*lddb + (i)+b_offset)
#define dC(id, i, j)  (dC[(id)]+(j)*lddc + (i))
#define STREAM_ID(i) (nqueue > 1 ? 1+((i)/nb)%(nqueue-1) : 0)

    magma_minproduct_int_t i, id, ib, ii, kk, n1;
    double z_alpha = MAGMA_minproduct_D_MAKE(alpha,0.0);
    double z_beta  = MAGMA_minproduct_D_MAKE(beta, 0.0);

    magma_minproduct_device_t orig_dev;
    magma_minproduct_getdevice( &orig_dev );
    magma_minproduct_queue_t orig_stream;
    magma_minproductblasGetKernelStream( &orig_stream );
    
    /* diagonal update */
    for( i=0; i < n; i += nb ) {
        id = ((i+c_offset)/nb)%ngpu;
        kk = STREAM_ID( i+c_offset );

        ib = min(nb, n-i);
        ii = nb*((i+c_offset)/(nb*ngpu));

        /* dsyr2k on diagonal block */
        magma_minproduct_setdevice(id);
        magma_minproductblasSetKernelStream( queues[id][kk] );
        trace_gpu_start( id, kk, "syr2k", "syr2k" );
        magma_minproduct_dsyrk(uplo, trans, ib, k,
                    alpha,  dB(id, i,          0 ), lddb,
                     beta,  dC(id, i+c_offset, ii), lddc);
        trace_gpu_end( id, kk );
    }

    /* off-diagonal update */
    if (uplo == Magma_minproductUpper) {
        for( i=nb; i < n; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = STREAM_ID( i+c_offset );

            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));

            magma_minproduct_setdevice(id);
            magma_minproductblasSetKernelStream( queues[id][kk] );
            magma_minproduct_dgemm(Magma_minproductNoTrans, Magma_minproductConjTrans, i, ib, k,
                        z_alpha, dB(id, 0, 0 ), lddb,
                                 dB(id, i, 0 ), lddb,
                        z_beta,  dC(id, 0, ii), lddc);
        }
    }
    else {
        for( i=0; i < n-nb; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = STREAM_ID( i+c_offset );

            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));
            n1 = n-i-ib;

            /* dgemm on off-diagonal blocks */
            magma_minproduct_setdevice(id);
            magma_minproductblasSetKernelStream( queues[id][kk] );
            trace_gpu_start( id, kk, "gemm_up", "gemm_up" );
            magma_minproduct_dgemm(Magma_minproductNoTrans, Magma_minproductConjTrans, n1, ib, k,
                        z_alpha, dB(id, i+ib,           0 ), lddb,
                                 dB(id,  i,             0 ), lddb,
                        z_beta,  dC(id,  i+c_offset+ib, ii), lddc);
            trace_gpu_end( id, kk );
        }
    }

    // TODO why not sync?
    //for( id=0; id < ngpu; id++ ) {
    //    magma_minproduct_setdevice(id);
    //    //for( kk=0; kk < nqueue; kk++ )
    //    //    magma_minproduct_queue_sync( queues[id][kk] );
    //}
    magma_minproduct_setdevice( orig_dev );
    magma_minproductblasSetKernelStream( orig_stream );
}
#undef dB
#undef dC
#undef STREAM_ID

// ----------------------------------------------------------------------
extern "C" void
magma_minproduct_dsyrk_mgpu2(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t nb, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_ptr dB[], magma_minproduct_int_t lddb, magma_minproduct_int_t b_offset,
    double beta,
    magma_minproductDouble_ptr dC[], magma_minproduct_int_t lddc, magma_minproduct_int_t c_offset,
    magma_minproduct_int_t nqueue, magma_minproduct_queue_t queues[][10])
{
#define dB(id, i, j)  (dB[(id)]+(j)*lddb + (i)+b_offset)
#define dC(id, i, j)  (dC[(id)]+(j)*lddc + (i))
#define STREAM_ID(i) (nqueue > 1 ? 1+((i)/nb)%(nqueue-1) : 0)

    magma_minproduct_int_t i, id, ib, ii, kk, n1;
    double z_alpha = MAGMA_minproduct_D_MAKE(alpha,0.0);
    double z_beta  = MAGMA_minproduct_D_MAKE(beta, 0.0);

    magma_minproduct_device_t orig_dev;
    magma_minproduct_getdevice( &orig_dev );
    magma_minproduct_queue_t orig_stream;
    magma_minproductblasGetKernelStream( &orig_stream );
    
    /* diagonal update */
    for( i=0; i < n; i += nb ) {
        id = ((i+c_offset)/nb)%ngpu;
        kk = STREAM_ID( i+c_offset );

        ib = min(nb, n-i);
        ii = nb*((i+c_offset)/(nb*ngpu));
    }

    if (uplo == Magma_minproductUpper) {
        for( i=0; i < n; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = STREAM_ID( i+c_offset );

            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));
            n1 = i+ib;

            magma_minproduct_setdevice(id);
            magma_minproductblasSetKernelStream( queues[id][kk] );

            /* dgemm on diag and off-diagonal blocks */
            magma_minproduct_dgemm(Magma_minproductNoTrans, Magma_minproductConjTrans, n1, ib, k,
                        z_alpha, dB(id, 0, 0 ), lddb,
                                 dB(id, i, 0 ), lddb,
                        z_beta,  dC(id, 0, ii), lddc);
        }
    }
    else {
        for( i=0; i < n; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = STREAM_ID( i+c_offset );

            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));
            n1 = n-i;

            magma_minproduct_setdevice(id);
            magma_minproductblasSetKernelStream( queues[id][kk] );
            trace_gpu_start( id, kk, "gemm_up", "gemm_up" );
            /* dgemm on diag and off-diagonal blocks */
            magma_minproduct_dgemm(Magma_minproductNoTrans, Magma_minproductConjTrans, n1, ib, k,
                        z_alpha, dB(id, i,           0), lddb,
                                 dB(id, i,           0), lddb,
                        z_beta,  dC(id, i+c_offset, ii), lddc);
            trace_gpu_end( id, kk );
        }
    }

    // TODO: why not sync?
    //for( id=0; id < ngpu; id++ ) {
    //    magma_minproduct_setdevice(id);
    //    //for( kk=0; kk < nqueue; kk++ )
    //    //    magma_minproduct_queue_sync( queues[id][kk] );
    //}
    magma_minproduct_setdevice( orig_dev );
    magma_minproductblasSetKernelStream( orig_stream );
}

#undef dB
#undef dC
#undef STREAM_ID

