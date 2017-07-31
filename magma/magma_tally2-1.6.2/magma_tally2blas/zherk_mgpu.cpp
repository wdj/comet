/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Azzam Haidar
       @author Ichi Yamazaki

       @precisions normal z -> s d c

*/
#include "common_magma_tally2.h"
#include "trace.h"

/**
    Purpose
    -------
    This zherk_mgpu is internal routine used by zpotrf_mgpu_right.
    it has specific assumption on the block diagonal.
    
    @ingroup magma_tally2_zblas3
    ********************************************************************/

extern "C" void
magma_tally2_zherk_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t nb, magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2DoubleComplex_ptr dB[], magma_tally2_int_t lddb, magma_tally2_int_t b_offset,
    double beta,
    magma_tally2DoubleComplex_ptr dC[], magma_tally2_int_t lddc, magma_tally2_int_t c_offset,
    magma_tally2_int_t nqueue, magma_tally2_queue_t queues[][10])
{
#define dB(id, i, j)  (dB[(id)]+(j)*lddb + (i)+b_offset)
#define dC(id, i, j)  (dC[(id)]+(j)*lddc + (i))
#define STREAM_ID(i) (nqueue > 1 ? 1+((i)/nb)%(nqueue-1) : 0)

    magma_tally2_int_t i, id, ib, ii, kk, n1;
    magma_tally2DoubleComplex z_alpha = MAGMA_tally2_Z_MAKE(alpha,0.0);
    magma_tally2DoubleComplex z_beta  = MAGMA_tally2_Z_MAKE(beta, 0.0);

    magma_tally2_device_t orig_dev;
    magma_tally2_getdevice( &orig_dev );
    magma_tally2_queue_t orig_stream;
    magma_tally2blasGetKernelStream( &orig_stream );
    
    /* diagonal update */
    for( i=0; i < n; i += nb ) {
        id = ((i+c_offset)/nb)%ngpu;
        kk = STREAM_ID( i+c_offset );

        ib = min(nb, n-i);
        ii = nb*((i+c_offset)/(nb*ngpu));

        /* zher2k on diagonal block */
        magma_tally2_setdevice(id);
        magma_tally2blasSetKernelStream( queues[id][kk] );
        trace_gpu_start( id, kk, "syr2k", "syr2k" );
        magma_tally2_zherk(uplo, trans, ib, k,
                    alpha,  dB(id, i,          0 ), lddb,
                     beta,  dC(id, i+c_offset, ii), lddc);
        trace_gpu_end( id, kk );
    }

    /* off-diagonal update */
    if (uplo == Magma_tally2Upper) {
        for( i=nb; i < n; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = STREAM_ID( i+c_offset );

            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));

            magma_tally2_setdevice(id);
            magma_tally2blasSetKernelStream( queues[id][kk] );
            magma_tally2_zgemm(Magma_tally2NoTrans, Magma_tally2ConjTrans, i, ib, k,
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

            /* zgemm on off-diagonal blocks */
            magma_tally2_setdevice(id);
            magma_tally2blasSetKernelStream( queues[id][kk] );
            trace_gpu_start( id, kk, "gemm_up", "gemm_up" );
            magma_tally2_zgemm(Magma_tally2NoTrans, Magma_tally2ConjTrans, n1, ib, k,
                        z_alpha, dB(id, i+ib,           0 ), lddb,
                                 dB(id,  i,             0 ), lddb,
                        z_beta,  dC(id,  i+c_offset+ib, ii), lddc);
            trace_gpu_end( id, kk );
        }
    }

    // TODO why not sync?
    //for( id=0; id < ngpu; id++ ) {
    //    magma_tally2_setdevice(id);
    //    //for( kk=0; kk < nqueue; kk++ )
    //    //    magma_tally2_queue_sync( queues[id][kk] );
    //}
    magma_tally2_setdevice( orig_dev );
    magma_tally2blasSetKernelStream( orig_stream );
}
#undef dB
#undef dC
#undef STREAM_ID

// ----------------------------------------------------------------------
extern "C" void
magma_tally2_zherk_mgpu2(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t nb, magma_tally2_int_t n, magma_tally2_int_t k,
    double alpha,
    magma_tally2DoubleComplex_ptr dB[], magma_tally2_int_t lddb, magma_tally2_int_t b_offset,
    double beta,
    magma_tally2DoubleComplex_ptr dC[], magma_tally2_int_t lddc, magma_tally2_int_t c_offset,
    magma_tally2_int_t nqueue, magma_tally2_queue_t queues[][10])
{
#define dB(id, i, j)  (dB[(id)]+(j)*lddb + (i)+b_offset)
#define dC(id, i, j)  (dC[(id)]+(j)*lddc + (i))
#define STREAM_ID(i) (nqueue > 1 ? 1+((i)/nb)%(nqueue-1) : 0)

    magma_tally2_int_t i, id, ib, ii, kk, n1;
    magma_tally2DoubleComplex z_alpha = MAGMA_tally2_Z_MAKE(alpha,0.0);
    magma_tally2DoubleComplex z_beta  = MAGMA_tally2_Z_MAKE(beta, 0.0);

    magma_tally2_device_t orig_dev;
    magma_tally2_getdevice( &orig_dev );
    magma_tally2_queue_t orig_stream;
    magma_tally2blasGetKernelStream( &orig_stream );
    
    /* diagonal update */
    for( i=0; i < n; i += nb ) {
        id = ((i+c_offset)/nb)%ngpu;
        kk = STREAM_ID( i+c_offset );

        ib = min(nb, n-i);
        ii = nb*((i+c_offset)/(nb*ngpu));
    }

    if (uplo == Magma_tally2Upper) {
        for( i=0; i < n; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = STREAM_ID( i+c_offset );

            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));
            n1 = i+ib;

            magma_tally2_setdevice(id);
            magma_tally2blasSetKernelStream( queues[id][kk] );

            /* zgemm on diag and off-diagonal blocks */
            magma_tally2_zgemm(Magma_tally2NoTrans, Magma_tally2ConjTrans, n1, ib, k,
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

            magma_tally2_setdevice(id);
            magma_tally2blasSetKernelStream( queues[id][kk] );
            trace_gpu_start( id, kk, "gemm_up", "gemm_up" );
            /* zgemm on diag and off-diagonal blocks */
            magma_tally2_zgemm(Magma_tally2NoTrans, Magma_tally2ConjTrans, n1, ib, k,
                        z_alpha, dB(id, i,           0), lddb,
                                 dB(id, i,           0), lddb,
                        z_beta,  dC(id, i+c_offset, ii), lddc);
            trace_gpu_end( id, kk );
        }
    }

    // TODO: why not sync?
    //for( id=0; id < ngpu; id++ ) {
    //    magma_tally2_setdevice(id);
    //    //for( kk=0; kk < nqueue; kk++ )
    //    //    magma_tally2_queue_sync( queues[id][kk] );
    //}
    magma_tally2_setdevice( orig_dev );
    magma_tally2blasSetKernelStream( orig_stream );
}

#undef dB
#undef dC
#undef STREAM_ID

