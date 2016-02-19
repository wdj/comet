/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zhemm_mgpu.cpp normal z -> s, Fri Jan 30 19:00:10 2015
       @author Mark Gates
       @author Azzam Haidar
       
       This still has poor performance. Work in progress.
*/
#include "common_magma_tally3.h"
#include "magma_tally3_bulge.h"
//#include "trace.h"

extern "C"
void magma_tally3blas_ssymm_mgpu_com(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_ptr dA[],    magma_tally3_int_t ldda,  magma_tally3_int_t offset,
    magma_tally3Float_ptr dB[],    magma_tally3_int_t lddb,
    float beta,
    magma_tally3Float_ptr dC[],    magma_tally3_int_t lddc,
    magma_tally3Float_ptr dwork[], magma_tally3_int_t dworksiz,
    float *C,          magma_tally3_int_t ldc,
    float *work[],     magma_tally3_int_t worksiz,  // TODO unused
    magma_tally3_int_t ngpu, magma_tally3_int_t nb, 
    magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue, 
    magma_tally3_event_t redevents[][Magma_tally3MaxGPUs*Magma_tally3MaxGPUs+10], magma_tally3_int_t nbevents, 
    magma_tally3_int_t gnode[Magma_tally3MaxGPUs][Magma_tally3MaxGPUs+2], magma_tally3_int_t nbcmplx )
{
    #define dA(dev, i, j) (dA[dev] + (i) + (j)*ldda)
    #define dB(dev, i, j) (dB[dev] + (i) + (j)*lddb)
    #define dC(dev, i, j) (dC[dev] + (i) + (j)*lddc)
    #define dwork(dev, i, j) (dwork[dev] + (i) + (j)*lddwork)
    #define C(i, j) (C + (i) + (j)*ldc)
    //printf("####################################################\n");
    //printf("                      start ssymm                   \n");
    //printf("####################################################\n");
   
    if ( side != Magma_tally3Left || uplo != Magma_tally3Lower ) {
        fprintf( stderr, "%s: only Left Lower implemented\n", __func__ );
    }
    
    assert( ldda >= m );
    assert( lddb >= m );
    assert( lddc >= m );
    assert( nqueue >= ngpu );
    assert( nbevents >= ngpu*ngpu );
   
    
    float c_one  = MAGMA_tally3_S_ONE;

    magma_tally3Float_ptr dwork2[Magma_tally3MaxGPUs];


    magma_tally3_int_t maxgsize    = n*m;
    magma_tally3_int_t lddwork = lddc;
    magma_tally3_int_t ldwork  = m;
    for( magma_tally3_int_t dev = 0; dev < ngpu; ++dev ) {
        dwork2[dev] = dwork[dev]+n*lddwork;  // size of dwork2 is maxgsize*ngpu
    }
    assert( dworksiz >= (n*lddwork+maxgsize*ngpu) );
    assert( worksiz  >= (n*ldwork) );

        
    magma_tally3_device_t cdev;
    magma_tally3_getdevice( &cdev );
    magma_tally3_queue_t cstream;
    magma_tally3blasGetKernelStream(&cstream);


    magma_tally3_int_t dev, devperm, myblk, mycolsize, myblkoffst;
    magma_tally3_int_t gmaster;
    magma_tally3_int_t masterdev, lcdev, lccolsize, myngpu;

    magma_tally3_int_t stdev       = (offset/nb)%ngpu;  
    magma_tally3_int_t blockoffset = offset % nb;  
    magma_tally3_int_t fstblksiz   = 0;
    if(blockoffset>0){
        fstblksiz   = min(m, (nb - blockoffset));
    }
    //magma_tally3_int_t nbblk       = magma_tally3_ceildiv(m, nb);
    magma_tally3_int_t nbblk       = magma_tally3_ceildiv((m+blockoffset), nb);
    magma_tally3_int_t remm        = m- fstblksiz;
    magma_tally3_int_t nbblkoffst  = offset/nb;


    magma_tally3_int_t nblstblks = -1;
    magma_tally3_int_t devlstblk = -1;
    magma_tally3_int_t lstblksiz = remm%nb;
    if(lstblksiz>0){
        nblstblks = nbblk%ngpu;
        devlstblk = (nblstblks-1+ngpu)%ngpu;
    }

    magma_tally3_int_t nbcmplxactive =  0;
    magma_tally3_int_t cmplxisactive[Magma_tally3MaxGPUs];
    magma_tally3_int_t gpuisactive[Magma_tally3MaxGPUs];
    memset(gpuisactive, 0, Magma_tally3MaxGPUs*sizeof(magma_tally3_int_t));
    memset(cmplxisactive, 0, Magma_tally3MaxGPUs*sizeof(magma_tally3_int_t));


    for( magma_tally3_int_t dev = 0; dev < ngpu; ++dev ) {
        magma_tally3_setdevice( dev );
        magma_tally3blasSetKernelStream( queues[ dev ][ 0 ] );
        cudaMemset(dwork(dev,0,0), 0, (lddwork)*(n)*sizeof(float) );
        // put all dC on all dev to 0 except the one which
        // hold i==0 because this one has to multiply by beta.
        if(dev!=stdev){
           cudaMemset(dC(dev,0,0), 0, (lddc)*(n)*sizeof(float) );
        }
    }

    magma_tally3_int_t newoffset = offset;
    // 1. symmetrize
    if(blockoffset>0){
        newoffset  = offset+fstblksiz; // newoffset is adjusted over nb
        magma_tally3_int_t myblkoffst = (nbblkoffst/ngpu)+(nbblkoffst%ngpu > stdev?1:0);
        //printf("STDEV %d  voici offset %d remm %d   myblockoffset %d    siz %d \n", stdev, offset, remm, myblkoffst, fstblksiz);
        magma_tally3_setdevice( stdev );
        magma_tally3blasSetKernelStream( queues[ stdev ][ 0 ] );
        magma_tally3blas_ssymmetrize_tiles(  Magma_tally3Lower,  fstblksiz,  dA(stdev, offset, myblkoffst*nb+blockoffset),  ldda,  1,  ngpu*nb,  nb  );         
    }

    for( magma_tally3_int_t dev = 0; dev < ngpu; ++dev ) {
        magma_tally3_int_t newstdev      = (newoffset/nb)%ngpu;
        magma_tally3_int_t nbblk = remm/nb; // number of block of size nb. if m%nb>0 then a last block exist and is of size ib=m%nb
        magma_tally3_int_t myblk = (nbblk/ngpu) + (nbblk%ngpu > ((dev-newstdev+ngpu)%ngpu) ?  1:0 );
        magma_tally3_int_t devperm   = (dev-newstdev+ngpu)%ngpu;
        magma_tally3_int_t nbblkoffst = newoffset/nb;
        magma_tally3_int_t myblkoffst = (nbblkoffst/ngpu)+(nbblkoffst%ngpu > dev?1:0);
        //printf("dev %d  devperm %d   newoffset %d  rowoff %d    coloff %d    myblk %d  \n", dev, devperm, newoffset, newoffset+devperm*nb, myblkoffst*nb, myblk);
        magma_tally3_setdevice( dev );
        magma_tally3blasSetKernelStream( queues[ dev ][ 0 ] );
        magma_tally3blas_ssymmetrize_tiles(  Magma_tally3Lower,  nb,  dA(dev, newoffset+devperm*nb, myblkoffst*nb),  ldda,  myblk,  ngpu*nb,  nb  );
        if(remm%nb>0){
            magma_tally3_int_t nblstblks = (nbblk+1)%ngpu;
            magma_tally3_int_t devlstblk = (nblstblks-1+ngpu)%ngpu;
            //printf("==> siz %d devperm %d,    devlstblk %d,    newoffset+nbblk*nb %d,   myblkoffst*nb+ myblk*nb %d\n", remm % nb, devperm, devlstblk, newoffset+nbblk*nb, myblkoffst*nb+ myblk*nb);
            if(devperm==devlstblk)
                magma_tally3blas_ssymmetrize(  Magma_tally3Lower,  remm % nb,  dA(dev, newoffset+nbblk*nb, myblkoffst*nb+ myblk*nb),  ldda );  // last partial tile
        }
    }


    

/*
    magma_tally3_int_t siz = m+offset;
    float *R;
    magma_tally3_smalloc_cpu( &R, siz*siz );
    // collecte back A
    magma_tally3blas_sgetmatrix_1D_bcyclic( siz, siz, dA, ldda, R, siz, ngpu, nb );
    magma_tally3_setdevice( 0 );
    magma_tally3blasSetKernelStream( queues[ dev ][ 0 ] );
    //magma_tally3_sgetmatrix( siz, siz, dA[0], ldda, R, siz );
    FILE *trace_file;
    trace_file = fopen("AJETE/Aafter", "w");
    for (int j = 0; j < siz ; j++) 
          for (int i = 0; i < siz ; i++) 
                         fprintf(trace_file, "%10d%10d%40.30e\n", i+1, j+1, R[j*siz+i]);
    fclose(trace_file);
return;
*/
    

    // ROW GEMM transpose a row and make a gemm with a block
    // if only 1 GPU used the ROW GEMM is integrated with the 
    // COL GEMM (better accuracy observed) and better perf
    if(ngpu>1){
        for( magma_tally3_int_t i = fstblksiz; i < m; i += nb ) {
            magma_tally3_int_t ib     = min( nb, m-i );      // block size
            magma_tally3_int_t ioff   = i + offset;          // start global index in parent matrix
            //magma_tally3_int_t dev    = (ioff / nb) % ngpu;
            magma_tally3_int_t nbblkoffst = offset/nb;
            magma_tally3_int_t nbblk      = magma_tally3_ceildiv(i, nb);
            for( magma_tally3_int_t dev = 0; dev < ngpu; ++dev ) {


                magma_tally3_int_t myblk = (nbblk/ngpu) + (nbblk%ngpu > ((dev-stdev+ngpu)%ngpu) ?  1:0 );
                magma_tally3_int_t myblkoffst = (nbblkoffst/ngpu)+(nbblkoffst%ngpu > dev?1:0);

                magma_tally3_int_t myrowsize = myblk * nb;
                magma_tally3_int_t coloffset = myblkoffst*nb;
                if(dev==stdev) {
                    myrowsize = myrowsize -blockoffset;
                    coloffset = myblkoffst*nb+blockoffset;
                }
                //printf("ROW GEMM: voici i %d   ib %d    ioff %d   nbblkoffst %d stdev %d  dev %d myblk %d  myblkoffset %d  coloffset %d  rowsize %d\n", i, ib, ioff, nbblkoffst, stdev, dev, myblk, myblkoffst, coloffset, myrowsize);
                if(myrowsize>0){
                    magma_tally3_setdevice( dev );
                    magma_tally3blasSetKernelStream( queues[ dev ][ 1 ] );    
                    magma_tally3_sgemm( Magma_tally3ConjTrans, Magma_tally3NoTrans, myrowsize, n, ib,
                                 alpha, dA(dev,ioff,coloffset), ldda,
                                        dB(dev,i,0),    lddb,
                                 c_one, dwork(dev,0,0), lddwork );
                }
            }
        }
        for( magma_tally3_int_t dev = 0; dev < ngpu; ++dev ) {
            magma_tally3_setdevice( dev );
            magma_tally3_event_record(redevents[dev][1], queues[dev][1]);
        }
    }
    

    // COL GEMM
    // blockoffset is offset within first block; for subsequent blocks it is 0
    if(blockoffset>0){
        magma_tally3_int_t ib     = min( nb-blockoffset, m );  // block size
        magma_tally3_int_t iblock = (offset / nb) / ngpu;          // local block id
        magma_tally3_int_t di     = iblock*nb+blockoffset;       // local index in parent matrix
        magma_tally3_setdevice( stdev );
        magma_tally3blasSetKernelStream( queues[ stdev ][ 0 ] );        
        //printf("DEV %d COL GEMM first   ioff %d  di %d   m %d   n %d   ib %d \n", stdev, offset, di, m, n, ib);
        magma_tally3_sgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, m, n, ib,
                        alpha, dA(stdev,offset,di), ldda,
                               dB(stdev,0,0),     lddb,
                        beta,  dC(stdev,0,0),     lddc );
    }
   


    // COL GEMM
    for( magma_tally3_int_t i = fstblksiz; i < m; i += nb ) {
        magma_tally3_int_t ib     = min( nb, m-i );      // block size
        magma_tally3_int_t ioff   = i + offset;          // start global index in parent matrix
        magma_tally3_int_t iblock = (ioff / nb) / ngpu;  // local block id
        magma_tally3_int_t dev    = (ioff / nb) % ngpu;
        magma_tally3_int_t di     = iblock*nb;           // local index in parent matrix
        
        //printf("DEV %d COL GEMM i %d      ioff %d  di %d m-i %d    n %d   ib %d \n", dev, i, ioff, di, m-i, n, ib);
        
        magma_tally3_setdevice( dev );
        magma_tally3blasSetKernelStream( queues[ dev ][ 0 ] );
        if(i==0){
           magma_tally3_sgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, m-i, n, ib,
                        alpha, dA(dev,ioff,di), ldda,
                               dB(dev,i,0),     lddb,
                        beta,  dC(dev,i,0),     lddc );
        }else{
           magma_tally3_sgemm( Magma_tally3NoTrans, Magma_tally3NoTrans, m-i, n, ib,
                        alpha, dA(dev,ioff,di), ldda,
                               dB(dev,i,0),        lddb,
                        c_one, dC(dev,i,0),     lddc );
        }
        magma_tally3_event_record(redevents[dev][0], queues[dev][0]);
        // if only 1 GPU is used, do the ROW GEMM
        if(ngpu==1){
            // NOTE THAT because the COL gemm write dC below the diagonal (i) 
            // and the ROW GEMM write dC from 0 to diag-1, so they could 
            // run in parallel on different queues. 
            // 
            // NO NO NO because
            // it might happen that col finished i and strated i+1 while row still at i    
            // magma_tally3blasSetKernelStream( queues[ dev ][ 0 ] );
            magma_tally3_sgemm( Magma_tally3ConjTrans, Magma_tally3NoTrans, i, n, ib,
                         alpha, dA(dev,ioff,offset), ldda,
                                dB(dev,i,0),    lddb,
                         c_one, dC(dev,0,0),    lddc );
        }
    }


    
    if(ngpu>1){
        for( magma_tally3_int_t dev = 0; dev < ngpu; ++dev ) {
            magma_tally3_int_t nbblk    = magma_tally3_ceildiv((m+blockoffset), nb);
            magma_tally3_int_t nbblkrow = nbblk-1; 
            magma_tally3_int_t devperm  = (dev-stdev+ngpu)%ngpu;
            magma_tally3_int_t myblk = (nbblkrow/ngpu) + (nbblkrow%ngpu > devperm ?  1:0 );
            magma_tally3_int_t myrowsize = myblk * nb;
             if(dev==stdev) {
                myrowsize = myrowsize - blockoffset;
            }
      
            //printf("blockoffset %d nbblkrow %d devperm %d  DEV %d RECEIVING myblk %d  myrowsize %d\n", blockoffset, nbblkrow, devperm, dev, myblk, myrowsize);
            if(myrowsize>0){
                magma_tally3_setdevice( dev );
                magma_tally3blasSetKernelStream( queues[ dev ][ 0 ] );
                magma_tally3_queue_wait_event(queues[ dev ][ 0 ], redevents[dev][1]);
                //magma_tally3_queue_sync( queues[ dev ][ 1 ] );
                // for each dev add the computed ROW block each on its placment with dC
                for( magma_tally3_int_t blki = 0; blki < myblk; ++blki){
                    magma_tally3_int_t gbblki = (blki*ngpu + devperm)*nb - blockoffset;
                    magma_tally3_int_t lcblki = blki*nb;
                    magma_tally3_int_t ib     = nb;// min(nb, m-gbblki);
                    if(dev==stdev){
                        lcblki = blki*nb-blockoffset;
                        if(blki==0){
                            gbblki = 0;
                            lcblki = 0;
                            ib     = nb-blockoffset;
                        }
                    }
                    magma_tally3blas_sgeadd(ib, n, c_one, 
                                    &dwork[dev][lcblki], lddwork, 
                                    &dC[dev][gbblki]   , lddc   );
                }
                magma_tally3_event_record(redevents[dev][0], queues[dev][0]);                
            }
        }
    }




    // ===========================================================
    //             COMMUNICATION ALL_REDUCE_SUM 
    // ===========================================================
    if(ngpu==1){
        return;
    }
    // INITIALIZE COMM
    for( magma_tally3_int_t cmplxid = 0; cmplxid < nbcmplx; ++cmplxid ) {
        masterdev     = -1;
        gnode[cmplxid][Magma_tally3MaxGPUs+1] = -1;
        myngpu = gnode[cmplxid][Magma_tally3MaxGPUs];
        for( magma_tally3_int_t idev = 0; idev < myngpu; ++idev ) {
            dev         = gnode[cmplxid][idev];
            devperm     = (dev-stdev+ngpu)%ngpu;
            myblk       = (nbblk/ngpu) + (nbblk%ngpu > devperm ?  1:0 );
            mycolsize   = myblk*nb;
            myblkoffst  = nb*((nbblkoffst/ngpu)+(nbblkoffst%ngpu > dev?1:0));            
            if(dev==stdev){
                mycolsize  -=  blockoffset;
                myblkoffst +=  blockoffset;     // local index in parent matrix
            }
            if((devperm==devlstblk)&&(lstblksiz>0)){
                mycolsize -=  (nb-(remm%nb));
            }
            mycolsize = min(mycolsize, m);
            if(mycolsize>0){
                gpuisactive[dev] = mycolsize;
                if(masterdev==-1) {
                    masterdev     = dev;
                    nbcmplxactive = nbcmplxactive +1;
                    cmplxisactive[cmplxid] = 1;
                    gnode[cmplxid][Magma_tally3MaxGPUs+1] = masterdev;
                }
            }
        }
    }
/*
    for( magma_tally3_int_t dev = 0; dev < ngpu; ++dev ) {
        magma_tally3_setdevice( dev );
        magma_tally3_device_sync();
    }
*/
    //*******************************
    //  each GPU send its result
    //  to its master. The master make
    //  the addition and then send to 
    //  to the masters of other real
    //  and receive from the masters of 
    //  other real make the addition 
    //  and broadcast locally the final 
    //  result.
    //*******************************
    //printf("=======================================================================\n");
    //printf("                     sending to my master                             \n");
    //printf("=======================================================================\n");
    for( magma_tally3_int_t cmplxid = 0; cmplxid < nbcmplx; ++cmplxid ) {
        myngpu    = gnode[cmplxid][Magma_tally3MaxGPUs];
        masterdev = gnode[cmplxid][Magma_tally3MaxGPUs+1];
        //check if real is active
        if(masterdev!=-1){ 
            for( magma_tally3_int_t idev = 0; idev < myngpu; ++idev ) {
                dev         = gnode[cmplxid][idev];
                mycolsize   = gpuisactive[dev];
                if(mycolsize>0){
                    // I am an active GPU. if I am not the master, then send my result to my master.
                    // store result on dwork[masterdev][dev*maxgsize]
                    if(dev!=masterdev){
                        magma_tally3_setdevice( dev );        
                        //printf("             GPU %d sending to my master %d\n", dev, masterdev);
                        // wait the geadd of my ROW and COL GEMM is done
                        magma_tally3_queue_wait_event(queues[ dev ][ 0 ], redevents[dev][0]);
                        // sending to the master of my real
                        magma_tally3_scopymatrix_async(
                            m, n,
                            &dC[dev][0], lddc,
                            &dwork2[masterdev][maxgsize*dev], m, queues[dev][0] );
                        magma_tally3_event_record(redevents[dev][masterdev], queues[dev][0]);
                    } // end I am not the masterdev
                }// end if mycolsize>0
            }// for idev
        }// end of if masterdev!=-1 maening real is active
    }// for cmplxid
/*
    for( magma_tally3_int_t dev = 0; dev < ngpu; ++dev ) {
        magma_tally3_setdevice( dev );
        magma_tally3_device_sync();
    }
*/

    //printf("=======================================================================\n");
    //printf(" each master do addition of local result and broadcast to other masters \n");
    //printf("=======================================================================\n");
    for( magma_tally3_int_t cmplxid = 0; cmplxid < nbcmplx; ++cmplxid ) {
        myngpu    = gnode[cmplxid][Magma_tally3MaxGPUs];
        masterdev = gnode[cmplxid][Magma_tally3MaxGPUs+1];
        //check if real is active
        if(masterdev!=-1){ 
            magma_tally3_setdevice( masterdev ); 
            // addition is done on stream 0 sequentially
            magma_tally3blasSetKernelStream( queues[ masterdev ][ 0 ] );
            // wait the geadd of my ROW and COL GEMM is done
            magma_tally3_queue_wait_event(queues[ masterdev ][ 0 ], redevents[masterdev][0]);
            // ========================================
            //     local addition
            // ========================================
            for( magma_tally3_int_t l = 0; l < myngpu; ++l ) {
                lcdev         = gnode[cmplxid][l];
                lccolsize     = gpuisactive[lcdev];
                if((lcdev!=masterdev)&&(lccolsize>0)){
                    //printf("             master %d receiving from %d and adding \n", masterdev, lcdev);
                    // this is an active GPU of my real. 
                    // wait I received what he send it to me and then do addition.
                    magma_tally3_queue_wait_event(queues[ masterdev ][ 0 ], redevents[lcdev][masterdev]);
                    magma_tally3blas_sgeadd(m, n, c_one, 
                                    &dwork2[masterdev][maxgsize*lcdev], m, 
                                    &dC[masterdev][0]   , lddc   );
                }
            }// for l=1:myngpu
            // because addition is done sequentially on stream 0, 
            // I have to record this to be able to synch using it 
            magma_tally3_event_record(redevents[masterdev][masterdev], queues[masterdev][0]);
            // ========================================
            //
            // ========================================
            //      send to other masters
            // ========================================
            for( magma_tally3_int_t k = 0; k < nbcmplx; ++k ) {
                if(k!=cmplxid){
                    gmaster = gnode[k][Magma_tally3MaxGPUs+1];
                    if(gmaster!=-1){ //real is active
                         //Master has to  wait until finish the local addition then send using gmaster stream.
                         //use stream 0 to make it sequential or stream gmaster to make it parallel.
                         //Now both re the same.
                        //printf("             master %d from cmplx %d sending to other master %d on cmplx %d \n", masterdev, cmplxid, gmaster, k);
                        magma_tally3_queue_wait_event(queues[ masterdev ][ gmaster ], redevents[masterdev][masterdev]);
                        magma_tally3_scopymatrix_async(
                            m, n,
                            &dC[masterdev][0], lddc,
                            &dwork2[gmaster][maxgsize*masterdev], m, queues[masterdev][gmaster] );
                        magma_tally3_event_record(redevents[masterdev][gmaster], queues[masterdev][gmaster]);
                        magma_tally3_event_record(redevents[masterdev][masterdev], queues[masterdev][gmaster]);
                      } // end of gmaster!=-1
                } // end of k!=cmplxid
            }// for k = 0: nbcmplx
            // ========================================
        }// end of if masterdev!=-1 maening real is active
    }// for cmplxid
/*
    for( magma_tally3_int_t dev = 0; dev < ngpu; ++dev ) {
        magma_tally3_setdevice( dev );
        magma_tally3_device_sync();
    }
*/
    //printf("=======================================================================\n");
    //printf(" each master wait receiving other masters results, do the addition and broadcast locally \n");
    //printf("=======================================================================\n");
    for( magma_tally3_int_t cmplxid = 0; cmplxid < nbcmplx; ++cmplxid ) {
        myngpu    = gnode[cmplxid][Magma_tally3MaxGPUs];
        masterdev = gnode[cmplxid][Magma_tally3MaxGPUs+1];
        //check if real is active
        if(masterdev!=-1){ 
            magma_tally3_setdevice( masterdev ); 
            // addition is done on stream 0 sequentially
            magma_tally3blasSetKernelStream( queues[ masterdev ][ 0 ] );
            // master has to wait until finishing all the send to other masters.
            magma_tally3_queue_wait_event(queues[ masterdev ][ 0 ], redevents[masterdev][masterdev]);
            // ========================================
            //  addition of results from other masters
            // ========================================
            for( magma_tally3_int_t k = 0; k < nbcmplx; ++k ) {
                if(k!=cmplxid){
                    gmaster = gnode[k][Magma_tally3MaxGPUs+1];
                    if(gmaster!=-1){ //real is active
                        //Master has to  wait until receiving from gmaster, then do addition using stream 0
                        //printf("             master %d from cmplx %d receiving from other master %d on cmplx %d and adding \n", masterdev, cmplxid, gmaster, k);
                        magma_tally3_queue_wait_event(queues[ masterdev ][ 0 ], redevents[gmaster][masterdev]);
                        magma_tally3blas_sgeadd(m, n, c_one, 
                                        &dwork2[masterdev][maxgsize*gmaster], m, 
                                        &dC[masterdev][0]   , lddc   );
                    } // end of gmaster!=-1
                } // end of k!=cmplxid
            }// for k = 0: nbcmplx
            // because addition is done sequentially on stream 0, 
            // I have to record this to be able to synch using it 
            magma_tally3_event_record(redevents[masterdev][masterdev], queues[masterdev][0]);
            // ========================================
            // ========================================
            //     local broadcast of final results
            // ========================================
            for( magma_tally3_int_t l = 0; l < myngpu; ++l ) {
                lcdev         = gnode[cmplxid][l];
                lccolsize     = gpuisactive[lcdev];
                if((lcdev!=masterdev)&&(lccolsize>0)){
                    // this is an active GPU of my real. 
                    // wait the previous addition is done maening stream 0 is finished and broadcast sequentially for now.
                    // to make it parallel put stream lcdev instead of stream 0
                    //printf("             master %d broadcasting local to %d  \n", masterdev, lcdev);
                    magma_tally3_queue_wait_event(queues[ masterdev ][ 0 ], redevents[masterdev][masterdev]);
                    magma_tally3_scopymatrix_async(
                        m, n,
                        &dC[masterdev][0], lddc,
                        &dC[lcdev][0],     lddc, queues[masterdev][0] );
                    magma_tally3_event_record(redevents[masterdev][lcdev], queues[masterdev][0]);
                }
            }// for l=1:myngpu
            // ========================================
        }// end of if masterdev!=-1 maening real is active
    }// for cmplxid
/*
    for( magma_tally3_int_t dev = 0; dev < ngpu; ++dev ) {
        magma_tally3_setdevice( dev );
        magma_tally3_device_sync();
    }
*/


    for( magma_tally3_int_t cmplxid = 0; cmplxid < nbcmplx; ++cmplxid ) {
        myngpu    = gnode[cmplxid][Magma_tally3MaxGPUs];
        masterdev = gnode[cmplxid][Magma_tally3MaxGPUs+1];
        //check if real is active
        if(masterdev!=-1){ 
            for( magma_tally3_int_t l = 0; l < myngpu; ++l ) {
                lcdev         = gnode[cmplxid][l];
                lccolsize     = gpuisactive[lcdev];
                if(lccolsize>0){
                    magma_tally3_setdevice( lcdev );
                    magma_tally3_queue_wait_event(queues[ lcdev ][ 0 ], redevents[lcdev][0]);
                    magma_tally3_queue_wait_event(queues[ lcdev ][ 0 ], redevents[masterdev][lcdev]);
                }
            }// for l=1:myngpu
        }// end of if masterdev!=-1 maening real is active
    }// for cmplxid


 
   //printf("****************************************************\n");
   //printf("                      finish ssymm                   \n");
   //printf("****************************************************\n");

    magma_tally3_setdevice( cdev );
    magma_tally3blasSetKernelStream( cstream );

}
