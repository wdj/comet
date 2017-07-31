/*
 * Copyright (c) 2011      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *
 *     @author Azzam Haidar
 *     @author Stan Tomov
 *     @author Raffaele Solca
 *
 *     @precisions normal z -> s d c
 *
 */

#include "common_magma_tally2.h"
#include "magma_tally2_bulge.h"
#include "magma_tally2_zbulgeinc.h"


////////////////////////////////////////////////////////////////////////////////////////////////////

#define dE(dev,i,j)  (dE[dev]+(i) + ldde*(j))
#define V(j)     (V+(j))
#define T(j)     (T+(j))
/***************************************************************************
 *  Parallel apply Q2 from bulgechasing symetric matrices - static scheduling
 *  Lower case is treated
 **/
    /*
     * side == magma_tally2Left:
     *     meaning apply E = Q*E = (q_1*q_2*.....*q_n) * E ==> so
     *     traverse Vs in reverse order (forward) from q_n to q_1 Also
     *     E is splitten by block of col over cores because each apply
     *     consist in a block of row (horizontal block)
     */
    /*
     * side == magma_tally2Right:
     *     meaning apply E = E*Q = E * (q_1*q_2*.....*q_n) ==> so
     *     traverse Vs in normal order (forward) from q_1 to q_n Also
     *     E is splitten by block of row over core because each apply
     *     consist in a block of col (vertical block)
     */
/***************************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_zbulge_applyQ_v2_m(
    magma_tally2_int_t ngpu,
    magma_tally2_side_t side,
    magma_tally2_int_t NE, magma_tally2_int_t N,
    magma_tally2_int_t NB, magma_tally2_int_t Vblksiz,
    magma_tally2DoubleComplex *E, magma_tally2_int_t lde,
    magma_tally2DoubleComplex *V, magma_tally2_int_t ldv,
    magma_tally2DoubleComplex *T, magma_tally2_int_t ldt,
    magma_tally2_int_t *info)
{
    //%===========================
    //%   local variables
    //%===========================
    magma_tally2_int_t Vm, Vn, mt, nt;
    magma_tally2_int_t myrow, mycol, blkj, blki;
    magma_tally2_int_t blkid,vpos,tpos;
    magma_tally2_int_t firstrow, nbcolinvolvd;
    magma_tally2_int_t versionL  = 113;
    magma_tally2_int_t versionR  = 92;
    magma_tally2_int_t Vchunksiz = 10;
    *info=0;

    /* Quick return */
    if ( NE == 0 ) {
        return *info;
    }
    if ( N == 0 ) {
        return *info;
    }
    if ( NB == 0 ) {
        return *info;
    }
    /* ==========================================
     * some infos for developer
     * Initialisation and checking nb of cores
     * ==========================================*/
    /* we have 2 algo for left (113 114) and 2 algo for right (91 92)
     * which correspond to versionL versionR.
     * They are very similar (detail explained in tech report and matlab code)
     * however version 114 and 92 improve locality.
     * while version 113 is used in case WNATZ=1 (construct Q2) which allow
     * the construction to be done in an optimized way taking into
     * consideration that the matrix is Identity so making less flops.
     *
    */

    // Initialize streaming and events
    magma_tally2_device_sync();
    magma_tally2_device_t orig_dev;
    magma_tally2_getdevice( &orig_dev );
    magma_tally2_queue_t orig_stream;
    magma_tally2blasGetKernelStream( &orig_stream );

    magma_tally2_int_t  nbevents =2, nstream=2;
    magma_tally2_queue_t streams[Magma_tally2MaxGPUs][20];
    magma_tally2_event_t  myevent[Magma_tally2MaxGPUs][20];
    for( magma_tally2_int_t dev = 0; dev < ngpu; ++dev ) {
        magma_tally2_setdevice( dev );
        for( magma_tally2_int_t i = 0; i < nstream; ++i ) {
            magma_tally2_queue_create( &streams[dev][i] );
        }
        for( magma_tally2_int_t i = 0; i < nbevents; ++i ) {
            cudaEventCreateWithFlags(&myevent[dev][i],cudaEventDisableTiming);
        }
    }



    // Azzam 21/11/2012
    // NOTE THAT dwork was of size 2*NE*Vblksiz+...
    // but I am thinking why not modifing it to NE*Vblksiz+...
    // BUT NO because the 2* is used because of making 2 streams working and so
    // they might be using dwork in parallel
    magma_tally2DoubleComplex *dE[Magma_tally2MaxGPUs];
    magma_tally2DoubleComplex *dwork[Magma_tally2MaxGPUs], *dwork0[Magma_tally2MaxGPUs], *dwork1[Magma_tally2MaxGPUs];
    //magma_tally2DoubleComplex *dwvt[Magma_tally2MaxGPUs];
    magma_tally2DoubleComplex *dwvt0[Magma_tally2MaxGPUs], *dwvt1[Magma_tally2MaxGPUs];
    magma_tally2DoubleComplex *dT0[Magma_tally2MaxGPUs], *dV0[Magma_tally2MaxGPUs], *dT1[Magma_tally2MaxGPUs], *dV1[Magma_tally2MaxGPUs];
    magma_tally2_int_t dev;
    magma_tally2_int_t ldde = N;
    magma_tally2_int_t lddv = ldv;
    magma_tally2_int_t lddt = ldt;
    magma_tally2_int_t ne_loc = magma_tally2_ceildiv(NE, ngpu);
    if (ne_loc < 256)
       ne_loc=256;
    magma_tally2_int_t dwVTsiz  = lddv*Vblksiz;  // lddv*lddv + lddv*NE; // lddv*Vblksiz;
    magma_tally2_int_t dworksiz = ne_loc*Vblksiz;  // lddv*Vblksiz;      // NE*Vblksiz;

    ngpu = min(ngpu, magma_tally2_ceildiv(NE,ne_loc)); // Don't use GPU that will not have data.
    // copy dE to GPUs
    for (dev=0; dev < ngpu; ++dev) {
        magma_tally2_setdevice( dev );
        if (MAGMA_tally2_SUCCESS != magma_tally2_zmalloc( &dE[dev], ldde * ne_loc)) {
            printf ("!!!!  magma_tally2_zbulge_applyQ magma_tally2_alloc failed for: dE\n" );
            exit(-1);
        }
        if (MAGMA_tally2_SUCCESS != magma_tally2_zmalloc( &dwork[dev], 2*dworksiz + 2*dwVTsiz + 2*Vchunksiz* (Vblksiz* (lddv+lddt)) )) {
            printf ("!!!!  magma_tally2_zbulge_applyQ magma_tally2_alloc failed for: dwork\n" );
            exit(-1);
        }

        dwork0[dev] = dwork[dev];               // size = dworksiz;
        dwork1[dev] = dwork0[dev] + dworksiz;   // size = dworksiz;
        dwvt0[dev]  = dwork[dev] + 2*dworksiz;  // size = dwVTsiz;
        dwvt1[dev]  = dwvt0[dev] + dwVTsiz;     // size = dwVTsiz;
        dV0[dev]    = dwork[dev] + 2*dworksiz + 2*dwVTsiz;
        dT0[dev]    = dV0[dev]   + Vchunksiz*Vblksiz*lddv;
        dV1[dev]    = dT0[dev]   + Vchunksiz*Vblksiz*lddt;
        dT1[dev]    = dV1[dev]   + Vchunksiz*Vblksiz*lddv;

        magma_tally2_int_t ie_loc = min(ne_loc, NE - ne_loc*dev);
        magma_tally2_zsetmatrix_async( N, ie_loc, E+lde*ne_loc*dev, lde, dE(dev, 0, 0), ldde, streams[dev][1] );
    }


    // make overlapped copy
    magma_tally2_int_t ncpy = 0;
    magma_tally2_int_t copyed=0, copyst=0;
    magma_tally2_int_t blkcnt,nothing, mysiz, flip, vld,tld, locpos;
    findVTsiz(N, NB, Vblksiz, &blkcnt, &nothing);
    
    flip = 0;


    /* SIDE LEFT  meaning apply E = Q*E = (q_1*q_2*.....*q_n) * E ==> so traverse Vs in reverse order (forward) from q_n to q_1
     *            Also E is splitten by row meaning each apply consist in a block of row (horizontal block) */
    /* SIDE RIGHT meaning apply E = E*Q = E * (q_1*q_2*.....*q_n) ==> so tarverse Vs in normal  order (forward) from q_1 to q_n
     *            Also E is splitten by col meaning each apply consist in a block of col (vertical block) */

    #ifdef ENABLE_DEBUG
    printf("  APPLY Q_v22_m GPU with NGPU %d  N %d, NE %d,  NB %d, Vblksiz %d, versionL %d versionR %d  SIDE %c \n",
          ngpu, N, NE, NB, Vblksiz, versionL, versionR, side);
    #endif


    /*
     * Magma_tally2maLeft
     */
    if (side == Magma_tally2Left) {
        /*
         * Version 113:
         * loop over the block_col (nt) and for each find the
         * number of tiles (mt) in this block_col. then loop over mt, find
         * the size of the V's(Vm,Vn) and apply it to the corresponding
         * portion of E.
         */
        if ( versionL == 113 ) {
            nt  = magma_tally2_ceildiv((N-1),Vblksiz);
            for (blkj=nt-1; blkj >= 0; blkj--) {
                /* the index of the first row on the top of block (blkj) */
                firstrow = blkj * Vblksiz + 1;
                /*find the number of tile for this block */
                if ( blkj == nt-1 )
                    mt = magma_tally2_ceildiv( N -  firstrow,    NB);
                else
                    mt = magma_tally2_ceildiv( N - (firstrow+1), NB);
                /*loop over the tiles find the size of the Vs and apply it */
                for (blki=mt; blki > 0; blki--) {
                    /*calculate the size of each losange of Vs= (Vm,Vn)*/
                    myrow     = firstrow + (mt-blki)*NB;
                    mycol     = blkj*Vblksiz;
                    Vm = min( NB+Vblksiz-1, N-myrow);
                    if ( ( blkj == nt-1 ) && ( blki == mt ) ) {
                        Vn = min (Vblksiz, Vm);
                    } else {
                        Vn = min (Vblksiz, Vm-1);
                    }
                    /*calculate the pointer to the Vs and the Ts.
                     * Note that Vs and Ts have special storage done
                     * by the bulgechasing function*/
                    //printf("voici blkj %d blki %d  Vm %d  Vn %d mycol %d vpos %d \n",blkj,blki,Vm, Vn,mycol,vpos);
                    magma_tally2_bulge_findpos113(N, NB, Vblksiz, mycol, myrow, &blkid);
               
                    // COPY Vchunksiz Vs and Vchunksiz Ts to GPU and store it in dV0/dV1 and dT0/dT1
                    if (ncpy == 0) {
                        // flip = 1 for this.
                        copyst = 0;                             // meaning that copy will start copying from blkid =copyst
                        copyed = min(copyst+Vchunksiz, blkcnt); // meaning that copy will end copying at blkid =copyed-1==> next copy had to start at copyed
                        mysiz  = copyed-copyst;                 // the size of the chunk to be copied
                        if (mysiz > 0) {
                            ncpy = 1;
                            flip = 1;
                            vpos = copyst*Vblksiz*ldv;
                            tpos = copyst*Vblksiz*ldt;
                            vld  = mysiz * ldv;
                            tld  = mysiz * ldt;
                            for( dev = 0; dev < ngpu; ++dev ) {
                                magma_tally2_setdevice( dev );
                                magma_tally2blasSetKernelStream( streams[ dev ][ 1 ] );
                                magma_tally2_zsetmatrix_async(vld, Vblksiz, V(vpos), vld, dV1[dev], vld, streams[dev][1]);
                                magma_tally2_zsetmatrix_async(tld, Vblksiz, T(tpos), tld, dT1[dev], tld, streams[dev][1]);
                            }
                            //printf("doing the first copy   of mysiz %2d copyst %2d copyed %2d vpos %8d tpos %8d into dV1 dT1\n",mysiz,copyst,copyed,vpos,tpos);
                        }
                    }
                   
                    if (blkid == copyst) {
                        flip   = ncpy % 2;
                        copyst = copyed;                             // meaning that copy will start copying from blkid =copyst
                        copyed = min(copyst+Vchunksiz, blkcnt); // meaning that copy will end copying at blkid =copyed-1==> next copy had to start at copyed
                        mysiz  = copyed-copyst;                 // the size of the chunk to be copied
                        //printf(" get to copy blkid %d blkid+(2*Vchunksiz) %d copyst %d copyed %d\n",blkid,blkid+(Vchunksiz),copyst,copyed);
                        if (mysiz > 0) {
                            ncpy = ncpy + 1;
                            vpos = copyst*Vblksiz*ldv;
                            tpos = copyst*Vblksiz*ldt;
                            vld  = mysiz * ldv;
                            tld  = mysiz * ldt;
                            if (flip == 0) { // now I am working on dV0 so copy the next and put it on dV1
                                //printf("doing overlapping copy of mysiz %2d copyst %2d copyed %2d vpos %8d tpos %8d into dV1 dT1\n",mysiz,copyst,copyed,vpos,tpos);
                                for( dev = 0; dev < ngpu; ++dev ) {
                                    magma_tally2_setdevice( dev );
                                    magma_tally2blasSetKernelStream( streams[ dev ][ 1 ] );
                                    magma_tally2_zsetmatrix_async(vld, Vblksiz, V(vpos), vld, dV1[dev], vld, streams[dev][1]);
                                    magma_tally2_zsetmatrix_async(tld, Vblksiz, T(tpos), tld, dT1[dev], tld, streams[dev][1]);
                                }
                            } else { // now I am working on dV1 so copy the next and put it on dV0
                                //printf("doing overlapping copy of mysiz %2d copyst %2d copyed %2d vpos %8d tpos %8d into dV0 dT0\n",mysiz,copyst,copyed,vpos,tpos);
                                for( dev = 0; dev < ngpu; ++dev ) {
                                    magma_tally2_setdevice( dev );
                                    magma_tally2blasSetKernelStream( streams[ dev ][ 0 ] );
                                    magma_tally2_zsetmatrix_async(vld, Vblksiz, V(vpos), vld, dV0[dev], vld, streams[dev][0]);
                                    magma_tally2_zsetmatrix_async(tld, Vblksiz, T(tpos), tld, dT0[dev], tld, streams[dev][0]);
                                }
                            }
                        }
                    }

                    if ((Vm > 0) && (Vn > 0)) {
                        locpos = blkid%Vchunksiz;
                        magma_tally2_int_t lcvpos   = locpos*Vblksiz*lddv;
                        magma_tally2_int_t lctpos   = locpos*Vblksiz*lddt;
                        //printf("voici blkj %d blki %d  Vm %d  Vn %d mycol %d locvpos %5d loctpos %5d  blkid %2d  using data in dV%1d dT%1d \n",blkj,blki,Vm, Vn,mycol,lcvpos,lctpos, blkid,flip,flip);
                        if (flip == 0) {
                            for( dev = 0; dev < ngpu; ++dev ) {
                                magma_tally2_int_t ie_loc   = min(ne_loc, NE - ne_loc*dev);
                                magma_tally2_int_t nr_bl = magma_tally2_ceildiv(ie_loc,10000);       //nr of blocks
                                magma_tally2_int_t sz_bl = magma_tally2_ceildiv(ie_loc,nr_bl*64)*64; //maximum size of blocks (to have blocks of around the same size and multiple of 64)
                                magma_tally2_int_t ib;                                        //size of current block
                                magma_tally2_setdevice( dev );
                                magma_tally2blasSetKernelStream(streams[dev][0]);
                                magma_tally2_queue_wait_event( streams[dev][0], myevent[dev][1] );
                                for (magma_tally2_int_t i=0; i < ie_loc; i += sz_bl) {
                                    ib = min(sz_bl, ie_loc-i);
                                    //magma_tally2_zlarfb_gpu( Magma_tally2Left, Magma_tally2NoTrans, Magma_tally2Forward, Magma_tally2Columnwise, Vm, ib, Vn, dV0[dev]+lcvpos, lddv, dT0[dev]+lctpos, lddt, dE(dev,myrow,i), ldde, dwork0[dev], ib);
                                    magma_tally2_zlarfb_gpu_gemm( Magma_tally2Left, Magma_tally2NoTrans, Magma_tally2Forward, Magma_tally2Columnwise, Vm, ib, Vn, dV0[dev]+lcvpos, lddv, dT0[dev]+lctpos, lddt, dE(dev,myrow,i), ldde, dwork0[dev], ib, dwvt0[dev], Vm);
                                }
                                magma_tally2_event_record( myevent[dev][0], streams[dev][0] );
                            }
                        } else {
                            for( dev = 0; dev < ngpu; ++dev ) {
                                magma_tally2_int_t ie_loc   = min(ne_loc, NE - ne_loc*dev);
                                magma_tally2_int_t nr_bl = magma_tally2_ceildiv(ie_loc,10000);       //nr of blocks
                                magma_tally2_int_t sz_bl = magma_tally2_ceildiv(ie_loc,nr_bl*64)*64; //maximum size of blocks (to have blocks of around the same size and multiple of 64)
                                magma_tally2_int_t ib;                                        //size of current block
                                magma_tally2_setdevice( dev );
                                magma_tally2blasSetKernelStream(streams[dev][1]);
                                magma_tally2_queue_wait_event( streams[dev][1], myevent[dev][0] );
                                for (magma_tally2_int_t i=0; i < ie_loc; i += sz_bl) {
                                    ib = min(sz_bl, ie_loc-i);
                                    //magma_tally2_zlarfb_gpu( Magma_tally2Left, Magma_tally2NoTrans, Magma_tally2Forward, Magma_tally2Columnwise, Vm, ib, Vn, dV1[dev]+lcvpos, lddv, dT1[dev]+lctpos, lddt, dE(dev,myrow,i), ldde, dwork1[dev], ib);
                                    magma_tally2_zlarfb_gpu_gemm( Magma_tally2Left, Magma_tally2NoTrans, Magma_tally2Forward, Magma_tally2Columnwise, Vm, ib, Vn, dV1[dev]+lcvpos, lddv, dT1[dev]+lctpos, lddt, dE(dev,myrow,i), ldde, dwork1[dev], ib, dwvt1[dev], Vm);
                                }
                                magma_tally2_event_record( myevent[dev][1], streams[dev][1] );
                            }
                        }
                    }  // end for (Vm &Vn) > 0
                } // end for blki
            } // end for blkj
        } // end if version=113
        /*
         * Version 114:
         * loop over the block_row (mt) and for each find diagonally the
         * number of tiles (nt) in this block_row. then loop over nt, find
         * the size of the V's(Vm,Vn) and apply it to the corresponding
         * portion of E.
         */
        else {
            printf("versionL 114 not implemented in zbulge_applyQ_v2_m\n");
            exit(-1);
            mt    = magma_tally2_ceildiv((N-1),NB);
            for (blki = mt; blki > 0; blki--) {
                /* nbcolinvolvd = number of column corresponding to this block_row (blki) */
                nbcolinvolvd = min(N-1, blki*NB);
                /*find the number of tile for this block (diagonal row of tiles) */
                nt = magma_tally2_ceildiv(nbcolinvolvd,Vblksiz);
                /*loop over the tiles find the size of the Vs and apply it */
                for (blkj = nt-1; blkj >= 0; blkj--) {
                    /* the index of the first row of the first col meaning
                     * the block on the top left (blki) */
                    firstrow = (mt-blki)*NB+1;
                    /*calculate the size of each losange of Vs= (Vm,Vn)*/
                    myrow    = firstrow + blkj*Vblksiz;
                    mycol    = blkj*Vblksiz;
                    Vm = min( NB+Vblksiz-1, N-myrow);
                    if ( ( blkj == nt-1 ) && ( blki == mt ) ) {
                        Vn = min (Vblksiz, Vm);
                    } else {
                        Vn = min (Vblksiz, Vm-1);
                    }

                    if ((Vm > 0) && (Vn > 0)) {
                        /*calculate the pointer to the Vs and the Ts.
                         * Note that Vs and Ts have special storage done
                         * by the bulgechasing function*/
                        /*
                        magma_tally2_bulge_findVTpos(N, NB, Vblksiz, mycol, myrow, ldv, ldt, &vpos, &tpos);
                        magma_tally2_zsetmatrix_async(Vm, Vn, V(vpos), ldv, dV0, lddv, NULL);
                        magma_tally2_zsetmatrix_async(Vn,  Vn, T(tpos), ldt, dT0, lddt, NULL);
                        //printf("voici blki %d  rownbm %d mycol %d  coled %d  blkid %d vpos %d  tpos %d\n", blki, rownbm, mycol, coled, blkid, vpos, tpos);
                        for (magma_tally2_int_t i=0; i < NE; i += sz_bl) {
                            ib = min(sz_bl, NE-i);
                            magma_tally2_zlarfb_gpu( Magma_tally2Left, Magma_tally2NoTrans, Magma_tally2Forward, Magma_tally2Columnwise, Vm, ib, Vn, dV0, lddv, dT0, lddt, dE(myrow,i), ldde, dwork, NE);
                        }
                        */
                    } // end for (Vm &Vn) > 0
                } // end for blkj
            } // end for blki
        } // end version 114
    } // end LEFT
    /*
     * Magma_tally2Right
     */
    else {
        printf("Side 'R' not implemented in zbulge_applyQ_v2_m\n");
        exit(-1);
        /*
         * Version 91:
         */
        if ( versionR == 91 ) {
            nt  = magma_tally2_ceildiv((N-1),Vblksiz);
            for (blkj=0; blkj < nt; blkj++) {
                /* the index of the first myrow on the top of block (blkj) */
                firstrow = blkj * Vblksiz + 1;
                /*find the number of tile for this block */
                if ( blkj == nt-1 )
                    mt = magma_tally2_ceildiv( N -  firstrow,    NB);
                else
                    mt = magma_tally2_ceildiv( N - (firstrow+1), NB);
                /*loop over the tiles find the size of the Vs and apply it */
                for (blki=1; blki <= mt; blki++) {
                    /*calculate the size of each losange of Vs= (Vm,Vn)*/
                    myrow  = firstrow + (mt-blki)*NB;
                    Vm = min( NB+Vblksiz-1, N-myrow);
                    if ( (blkj == nt-1) && (blki == mt) ) {
                        Vn = min (Vblksiz, Vm);
                    } else {
                        Vn = min (Vblksiz, Vm-1);
                    }
                    mycol     = blkj*Vblksiz;
                    if ((Vm > 0) && (Vn > 0)) {
                        /*calculate the pointer to the Vs and the Ts.
                         * Note that Vs and Ts have special storage done
                         * by the bulgechasing function*/
                        /*
                        magma_tally2_bulge_findVTpos(N, NB, Vblksiz, mycol, myrow, ldv, ldt, &vpos, &tpos);
                        magma_tally2_zsetmatrix_async(Vm, Vn, V(vpos), ldv, dV0, lddv, NULL);
                        magma_tally2_zsetmatrix_async(Vn,  Vn, T(tpos), ldt, dT0, lddt, NULL);
                        magma_tally2_zlarfb_gpu( Magma_tally2Right, Magma_tally2NoTrans, Magma_tally2Forward, Magma_tally2Columnwise, NE, Vm, Vn, dV0, lddv, dT0, lddt, dE(0, myrow), ldde, dwork, NE);
                        */
                    } // end for (Vm &Vn) > 0
                } // end for blki
            } // end fo blkj
        } // end of version 91
        /*
         * Version 92:
         */
        else {
            mt    = magma_tally2_ceildiv((N-1),NB);
            for (blki = 1; blki <= mt; blki++) {
                /* nbcolinvolvd = number of column corresponding to this block_row (blki) */
                nbcolinvolvd = min(N-1, blki*NB);
                /*find the number of tile for this block (diagonal row of tiles) */
                nt = magma_tally2_ceildiv(nbcolinvolvd,Vblksiz);
                /*loop over the tiles find the size of the Vs and apply it */
                for (blkj = 0; blkj < nt; blkj++) {
                    /* the index of the first row of the first col meaning
                     * the block on the top left (blki) */
                    firstrow = (mt-blki)*NB+1;
                    /*calculate the size of each losange of Vs= (Vm,Vn)*/
                    myrow    = firstrow + blkj*Vblksiz;
                    mycol    = blkj*Vblksiz;
                    Vm = min( NB+Vblksiz-1, N-myrow);
                    if ( ( blkj == nt-1 ) && ( blki == mt ) ) {
                        Vn = min (Vblksiz, Vm);
                    } else {
                        Vn = min (Vblksiz, Vm-1);
                    }
                    if ((Vm > 0) && (Vn > 0)) {
                        /*calculate the pointer to the Vs and the Ts.
                         * Note that Vs and Ts have special storage done
                         * by the bulgechasing function*/
                        /*
                        magma_tally2_bulge_findVTpos(N, NB, Vblksiz, mycol, myrow, ldv, ldt, &vpos, &tpos);
                        magma_tally2_zsetmatrix_async(Vm, Vn, V(vpos), ldv, dV0, lddv, NULL);
                        magma_tally2_zsetmatrix_async(Vn,  Vn, T(tpos), ldt, dT0, lddt, NULL);
                        magma_tally2_zlarfb_gpu( Magma_tally2Right, Magma_tally2NoTrans, Magma_tally2Forward, Magma_tally2Columnwise, NE, Vm, Vn, dV0, lddv, dT0, lddt, dE(0, myrow), ldde, dwork, NE);
                        */
                    } // end for (Vm &Vn) > 0
                } //end for blkj
            } // end for blki
        } //end of version 92
    } // end RIGHT

    // copy back the dE form each GPU
    for( dev = 0; dev < ngpu; ++dev ) {
        magma_tally2_setdevice( dev );
        magma_tally2blasSetKernelStream(streams[dev][0]);
        magma_tally2_queue_wait_event( streams[dev][0], myevent[dev][1] );
        magma_tally2_queue_wait_event( streams[dev][0], myevent[dev][0] );
        magma_tally2_int_t ie_loc   = min(ne_loc, NE - ne_loc*dev);
        magma_tally2_zgetmatrix_async( N, ie_loc, dE(dev, 0, 0), ldde, E+lde*ne_loc*dev, lde, streams[dev][0] );
        magma_tally2_event_record( myevent[dev][0], streams[dev][0] );
    }


    for( magma_tally2_int_t dev = 0; dev < ngpu; ++dev ) {
        magma_tally2_setdevice( dev );
        magma_tally2blasSetKernelStream(streams[dev][0]);
        magma_tally2_queue_wait_event( streams[dev][0], myevent[dev][0] );
        magma_tally2_device_sync(); // no need for synchronize
        magma_tally2_free(dwork[dev]);
        magma_tally2_free(dE[dev]);
        for( magma_tally2_int_t i = 0; i < nbevents; ++i ) {
            magma_tally2_event_destroy( myevent[dev][i] );
        }
        for( magma_tally2_int_t i = 0; i < nstream; ++i ) {
            magma_tally2_queue_destroy( streams[dev][i] );
        }
    }

    magma_tally2_setdevice( orig_dev );
    magma_tally2blasSetKernelStream( orig_stream );
    return *info;
}
#undef V
#undef T
#undef dE
////////////////////////////////////////////////////////////////////////////////////////////////////
