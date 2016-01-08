/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Azzam Haidar
       @author Ichi Yamazaki

       @precisions normal z -> s d c

*/
#include "common_magma_tally4.h"
#include "trace.h"

/**
    Purpose
    -------
    This zherk_mgpu is internal routine used by zpotrf_mgpu_right.
    it has specific assumption on the block diagonal.
    
    @ingroup magma_tally4_zblas3
    ********************************************************************/

extern "C" void
magma_tally4_zherk_mgpu(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t nb, magma_tally4_int_t n, magma_tally4_int_t k,
    double alpha,
    magma_tally4DoubleComplex_ptr dB[], magma_tally4_int_t lddb, magma_tally4_int_t b_offset,
    double beta,
    magma_tally4DoubleComplex_ptr dC[], magma_tally4_int_t lddc, magma_tally4_int_t c_offset,
    magma_tally4_int_t nqueue, magma_tally4_queue_t queues[][10])
{
#define dB(id, i, j)  (dB[(id)]+(j)*lddb + (i)+b_offset)
#define dC(id, i, j)  (dC[(id)]+(j)*lddc + (i))
#define STREAM_ID(i) (nqueue > 1 ? 1+((i)/nb)%(nqueue-1) : 0)

    magma_tally4_int_t i, id, ib, ii, kk, n1;
    magma_tally4DoubleComplex z_alpha = MAGMA_tally4_Z_MAKE(alpha,0.0);
    magma_tally4DoubleComplex z_beta  = MAGMA_tally4_Z_MAKE(beta, 0.0);

    magma_tally4_device_t orig_dev;
    magma_tally4_getdevice( &orig_dev );
    magma_tally4_queue_t orig_stream;
    magma_tally4blasGetKernelStream( &orig_stream );
    
    /* diagonal update */
    for( i=0; i < n; i += nb ) {
        id = ((i+c_offset)/nb)%ngpu;
        kk = STREAM_ID( i+c_offset );

        ib = min(nb, n-i);
        ii = nb*((i+c_offset)/(nb*ngpu));

        /* zher2k on diagonal block */
        magma_tally4_setdevice(id);
        magma_tally4blasSetKernelStream( queues[id][kk] );
        trace_gpu_start( id, kk, "syr2k", "syr2k" );
        magma_tally4_zherk(uplo, trans, ib, k,
                    alpha,  dB(id, i,          0 ), lddb,
                     beta,  dC(id, i+c_offset, ii), lddc);
        trace_gpu_end( id, kk );
    }

    /* off-diagonal update */
    if (uplo == Magma_tally4Upper) {
        for( i=nb; i < n; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = STREAM_ID( i+c_offset );

            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));

            magma_tally4_setdevice(id);
            magma_tally4blasSetKernelStream( queues[id][kk] );
            magma_tally4_zgemm(Magma_tally4NoTrans, Magma_tally4ConjTrans, i, ib, k,
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
            magma_tally4_setdevice(id);
            magma_tally4blasSetKernelStream( queues[id][kk] );
            trace_gpu_start( id, kk, "gemm_up", "gemm_up" );
            magma_tally4_zgemm(Magma_tally4NoTrans, Magma_tally4ConjTrans, n1, ib, k,
                        z_alpha, dB(id, i+ib,           0 ), lddb,
                                 dB(id,  i,             0 ), lddb,
                        z_beta,  dC(id,  i+c_offset+ib, ii), lddc);
            trace_gpu_end( id, kk );
        }
    }

    // TODO why not sync?
    //for( id=0; id < ngpu; id++ ) {
    //    magma_tally4_setdevice(id);
    //    //for( kk=0; kk < nqueue; kk++ )
    //    //    magma_tally4_queue_sync( queues[id][kk] );
    //}
    magma_tally4_setdevice( orig_dev );
    magma_tally4blasSetKernelStream( orig_stream );
}
#undef dB
#undef dC
#undef STREAM_ID

// ----------------------------------------------------------------------
extern "C" void
magma_tally4_zherk_mgpu2(
    magma_tally4_int_t ngpu,
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t nb, magma_tally4_int_t n, magma_tally4_int_t k,
    double alpha,
    magma_tally4DoubleComplex_ptr dB[], magma_tally4_int_t lddb, magma_tally4_int_t b_offset,
    double beta,
    magma_tally4DoubleComplex_ptr dC[], magma_tally4_int_t lddc, magma_tally4_int_t c_offset,
    magma_tally4_int_t nqueue, magma_tally4_queue_t queues[][10])
{
#define dB(id, i, j)  (dB[(id)]+(j)*lddb + (i)+b_offset)
#define dC(id, i, j)  (dC[(id)]+(j)*lddc + (i))
#define STREAM_ID(i) (nqueue > 1 ? 1+((i)/nb)%(nqueue-1) : 0)

    magma_tally4_int_t i, id, ib, ii, kk, n1;
    magma_tally4DoubleComplex z_alpha = MAGMA_tally4_Z_MAKE(alpha,0.0);
    magma_tally4DoubleComplex z_beta  = MAGMA_tally4_Z_MAKE(beta, 0.0);

    magma_tally4_device_t orig_dev;
    magma_tally4_getdevice( &orig_dev );
    magma_tally4_queue_t orig_stream;
    magma_tally4blasGetKernelStream( &orig_stream );
    
    /* diagonal update */
    for( i=0; i < n; i += nb ) {
        id = ((i+c_offset)/nb)%ngpu;
        kk = STREAM_ID( i+c_offset );

        ib = min(nb, n-i);
        ii = nb*((i+c_offset)/(nb*ngpu));
    }

    if (uplo == Magma_tally4Upper) {
        for( i=0; i < n; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = STREAM_ID( i+c_offset );

            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));
            n1 = i+ib;

            magma_tally4_setdevice(id);
            magma_tally4blasSetKernelStream( queues[id][kk] );

            /* zgemm on diag and off-diagonal blocks */
            magma_tally4_zgemm(Magma_tally4NoTrans, Magma_tally4ConjTrans, n1, ib, k,
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

            magma_tally4_setdevice(id);
            magma_tally4blasSetKernelStream( queues[id][kk] );
            trace_gpu_start( id, kk, "gemm_up", "gemm_up" );
            /* zgemm on diag and off-diagonal blocks */
            magma_tally4_zgemm(Magma_tally4NoTrans, Magma_tally4ConjTrans, n1, ib, k,
                        z_alpha, dB(id, i,           0), lddb,
                                 dB(id, i,           0), lddb,
                        z_beta,  dC(id, i+c_offset, ii), lddc);
            trace_gpu_end( id, kk );
        }
    }

    // TODO: why not sync?
    //for( id=0; id < ngpu; id++ ) {
    //    magma_tally4_setdevice(id);
    //    //for( kk=0; kk < nqueue; kk++ )
    //    //    magma_tally4_queue_sync( queues[id][kk] );
    //}
    magma_tally4_setdevice( orig_dev );
    magma_tally4blasSetKernelStream( orig_stream );
}

#undef dB
#undef dC
#undef STREAM_ID

