/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Azzam Haidar
       @author Ichi Yamazaki

       @generated from zherk_mgpu.cpp normal z -> s, Fri Jan 30 19:00:10 2015

*/
#include "common_magma_tally3.h"
#include "trace.h"

/**
    Purpose
    -------
    This ssyrk_mgpu is internal routine used by spotrf_mgpu_right.
    it has specific assumption on the block diagonal.
    
    @ingroup magma_tally3_sblas3
    ********************************************************************/

extern "C" void
magma_tally3_ssyrk_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t nb, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_ptr dB[], magma_tally3_int_t lddb, magma_tally3_int_t b_offset,
    float beta,
    magma_tally3Float_ptr dC[], magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3_int_t nqueue, magma_tally3_queue_t queues[][10])
{
#define dB(id, i, j)  (dB[(id)]+(j)*lddb + (i)+b_offset)
#define dC(id, i, j)  (dC[(id)]+(j)*lddc + (i))
#define STREAM_ID(i) (nqueue > 1 ? 1+((i)/nb)%(nqueue-1) : 0)

    magma_tally3_int_t i, id, ib, ii, kk, n1;
    float z_alpha = MAGMA_tally3_S_MAKE(alpha,0.0);
    float z_beta  = MAGMA_tally3_S_MAKE(beta, 0.0);

    magma_tally3_device_t orig_dev;
    magma_tally3_getdevice( &orig_dev );
    magma_tally3_queue_t orig_stream;
    magma_tally3blasGetKernelStream( &orig_stream );
    
    /* diagonal update */
    for( i=0; i < n; i += nb ) {
        id = ((i+c_offset)/nb)%ngpu;
        kk = STREAM_ID( i+c_offset );

        ib = min(nb, n-i);
        ii = nb*((i+c_offset)/(nb*ngpu));

        /* ssyr2k on diagonal block */
        magma_tally3_setdevice(id);
        magma_tally3blasSetKernelStream( queues[id][kk] );
        trace_gpu_start( id, kk, "syr2k", "syr2k" );
        magma_tally3_ssyrk(uplo, trans, ib, k,
                    alpha,  dB(id, i,          0 ), lddb,
                     beta,  dC(id, i+c_offset, ii), lddc);
        trace_gpu_end( id, kk );
    }

    /* off-diagonal update */
    if (uplo == Magma_tally3Upper) {
        for( i=nb; i < n; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = STREAM_ID( i+c_offset );

            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));

            magma_tally3_setdevice(id);
            magma_tally3blasSetKernelStream( queues[id][kk] );
            magma_tally3_sgemm(Magma_tally3NoTrans, Magma_tally3ConjTrans, i, ib, k,
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

            /* sgemm on off-diagonal blocks */
            magma_tally3_setdevice(id);
            magma_tally3blasSetKernelStream( queues[id][kk] );
            trace_gpu_start( id, kk, "gemm_up", "gemm_up" );
            magma_tally3_sgemm(Magma_tally3NoTrans, Magma_tally3ConjTrans, n1, ib, k,
                        z_alpha, dB(id, i+ib,           0 ), lddb,
                                 dB(id,  i,             0 ), lddb,
                        z_beta,  dC(id,  i+c_offset+ib, ii), lddc);
            trace_gpu_end( id, kk );
        }
    }

    // TODO why not sync?
    //for( id=0; id < ngpu; id++ ) {
    //    magma_tally3_setdevice(id);
    //    //for( kk=0; kk < nqueue; kk++ )
    //    //    magma_tally3_queue_sync( queues[id][kk] );
    //}
    magma_tally3_setdevice( orig_dev );
    magma_tally3blasSetKernelStream( orig_stream );
}
#undef dB
#undef dC
#undef STREAM_ID

// ----------------------------------------------------------------------
extern "C" void
magma_tally3_ssyrk_mgpu2(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t nb, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_ptr dB[], magma_tally3_int_t lddb, magma_tally3_int_t b_offset,
    float beta,
    magma_tally3Float_ptr dC[], magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3_int_t nqueue, magma_tally3_queue_t queues[][10])
{
#define dB(id, i, j)  (dB[(id)]+(j)*lddb + (i)+b_offset)
#define dC(id, i, j)  (dC[(id)]+(j)*lddc + (i))
#define STREAM_ID(i) (nqueue > 1 ? 1+((i)/nb)%(nqueue-1) : 0)

    magma_tally3_int_t i, id, ib, ii, kk, n1;
    float z_alpha = MAGMA_tally3_S_MAKE(alpha,0.0);
    float z_beta  = MAGMA_tally3_S_MAKE(beta, 0.0);

    magma_tally3_device_t orig_dev;
    magma_tally3_getdevice( &orig_dev );
    magma_tally3_queue_t orig_stream;
    magma_tally3blasGetKernelStream( &orig_stream );
    
    /* diagonal update */
    for( i=0; i < n; i += nb ) {
        id = ((i+c_offset)/nb)%ngpu;
        kk = STREAM_ID( i+c_offset );

        ib = min(nb, n-i);
        ii = nb*((i+c_offset)/(nb*ngpu));
    }

    if (uplo == Magma_tally3Upper) {
        for( i=0; i < n; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = STREAM_ID( i+c_offset );

            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));
            n1 = i+ib;

            magma_tally3_setdevice(id);
            magma_tally3blasSetKernelStream( queues[id][kk] );

            /* sgemm on diag and off-diagonal blocks */
            magma_tally3_sgemm(Magma_tally3NoTrans, Magma_tally3ConjTrans, n1, ib, k,
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

            magma_tally3_setdevice(id);
            magma_tally3blasSetKernelStream( queues[id][kk] );
            trace_gpu_start( id, kk, "gemm_up", "gemm_up" );
            /* sgemm on diag and off-diagonal blocks */
            magma_tally3_sgemm(Magma_tally3NoTrans, Magma_tally3ConjTrans, n1, ib, k,
                        z_alpha, dB(id, i,           0), lddb,
                                 dB(id, i,           0), lddb,
                        z_beta,  dC(id, i+c_offset, ii), lddc);
            trace_gpu_end( id, kk );
        }
    }

    // TODO: why not sync?
    //for( id=0; id < ngpu; id++ ) {
    //    magma_tally3_setdevice(id);
    //    //for( kk=0; kk < nqueue; kk++ )
    //    //    magma_tally3_queue_sync( queues[id][kk] );
    //}
    magma_tally3_setdevice( orig_dev );
    magma_tally3blasSetKernelStream( orig_stream );
}

#undef dB
#undef dC
#undef STREAM_ID

