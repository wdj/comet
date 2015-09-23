/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
*/

#ifndef MAGMA_minproduct_ZBULGEINC_H
#define MAGMA_minproduct_ZBULGEINC_H

#include "magma_minproduct_types.h"

#ifdef __cplusplus
extern "C" {
#endif


/***************************************************************************//**
 *  Configuration
 **/

 // maximum contexts
#define MAX_THREADS_BLG         256

void findVTpos(magma_minproduct_int_t N, magma_minproduct_int_t NB, magma_minproduct_int_t Vblksiz, magma_minproduct_int_t sweep, magma_minproduct_int_t st, magma_minproduct_int_t *Vpos, magma_minproduct_int_t *TAUpos, magma_minproduct_int_t *Tpos, magma_minproduct_int_t *myblkid);
void findVTsiz(magma_minproduct_int_t N, magma_minproduct_int_t NB, magma_minproduct_int_t Vblksiz, magma_minproduct_int_t *blkcnt, magma_minproduct_int_t *LDV);
magma_minproduct_int_t plasma_ceildiv(magma_minproduct_int_t a, magma_minproduct_int_t b);


/*
extern volatile magma_minproduct_int_t barrier_in[MAX_THREADS_BLG];
extern volatile magma_minproduct_int_t barrier_out[MAX_THREADS_BLG];
extern volatile magma_minproduct_int_t *ss_prog;
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
    magma_minproductDoubleComplex *dQ1;
    magma_minproductDoubleComplex *dT1;
    magma_minproductDoubleComplex *dT2;
    magma_minproductDoubleComplex *dV2;
    magma_minproductDoubleComplex *dE;
    magma_minproductDoubleComplex *T;
    magma_minproductDoubleComplex *A;
    magma_minproductDoubleComplex *V;
    magma_minproductDoubleComplex *TAU;
    magma_minproductDoubleComplex *E;
    magma_minproductDoubleComplex *E_CPU;
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
    magma_minproduct_side_t SIDE;
    real_Double_t *timeblg;
    real_Double_t *timeaplQ;
    volatile int *ss_prog;
} ;

// declare globals here; defined in zhetrd_bhe2trc.cpp
extern struct gbstrct_blg core_in_all;


#ifdef __cplusplus
}
#endif

#endif /* MAGMA_minproduct_ZBULGEINC_H */
