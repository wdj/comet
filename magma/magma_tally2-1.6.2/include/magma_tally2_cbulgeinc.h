/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from magma_tally2_zbulgeinc.h normal z -> c, Fri Jan 30 19:00:05 2015
*/

#ifndef MAGMA_tally2_CBULGEINC_H
#define MAGMA_tally2_CBULGEINC_H

#include "magma_tally2_types.h"

#ifdef __cplusplus
extern "C" {
#endif


/***************************************************************************//**
 *  Configuration
 **/

 // maximum contexts
#define MAX_THREADS_BLG         256

void findVTpos(magma_tally2_int_t N, magma_tally2_int_t NB, magma_tally2_int_t Vblksiz, magma_tally2_int_t sweep, magma_tally2_int_t st, magma_tally2_int_t *Vpos, magma_tally2_int_t *TAUpos, magma_tally2_int_t *Tpos, magma_tally2_int_t *myblkid);
void findVTsiz(magma_tally2_int_t N, magma_tally2_int_t NB, magma_tally2_int_t Vblksiz, magma_tally2_int_t *blkcnt, magma_tally2_int_t *LDV);
magma_tally2_int_t plasma_ceildiv(magma_tally2_int_t a, magma_tally2_int_t b);


/*
extern volatile magma_tally2_int_t barrier_in[MAX_THREADS_BLG];
extern volatile magma_tally2_int_t barrier_out[MAX_THREADS_BLG];
extern volatile magma_tally2_int_t *ss_prog;
*/

 /***************************************************************************//**
 *  Static scheduler
 **/
/*
#define ssched_init(nbtiles) \
{ \
        volatile int   prog_ol[2*nbtiles+10];\
                 int   iamdone[MAX_THREADS_BLG]; \
                 int   thread_num[MAX_THREADS_BLG];\
        pthread_t      thread_id[MAX_THREADS_BLG];\
        pthread_attr_t thread_attr;\
}
*/
////////////////////////////////////////////////////////////////////////////////////////////////////



 struct gbstrct_blg {
    magma_tally2FloatComplex *dQ1;
    magma_tally2FloatComplex *dT1;
    magma_tally2FloatComplex *dT2;
    magma_tally2FloatComplex *dV2;
    magma_tally2FloatComplex *dE;
    magma_tally2FloatComplex *T;
    magma_tally2FloatComplex *A;
    magma_tally2FloatComplex *V;
    magma_tally2FloatComplex *TAU;
    magma_tally2FloatComplex *E;
    magma_tally2FloatComplex *E_CPU;
    int cores_num;
    int locores_num;
    int overlapQ1;
    int usemulticpu;
    int NB;
    int NBTILES;
    int N;
    int NE;
    int N_CPU;
    int N_GPU;
    int LDA;
    int LDE;
    int BAND;
    int grsiz;
    int Vblksiz;
    int WANTZ;
    magma_tally2_side_t SIDE;
    real_Double_t *timeblg;
    real_Double_t *timeaplQ;
    volatile int *ss_prog;
} ;

// declare globals here; defined in chetrd_bhe2trc.cpp
extern struct gbstrct_blg core_in_all;


#ifdef __cplusplus
}
#endif

#endif /* MAGMA_tally2_CBULGEINC_H */
