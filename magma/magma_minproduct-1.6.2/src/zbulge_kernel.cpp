/*
 * Copyright (c) 2011      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *
 *     @author Azzam Haidar
 *     @author Stan Tomov
 *
 *     @precisions normal z -> s d c
 *
 */

#include "common_magma_minproduct.h"
#include "magma_minproduct_zbulgeinc.h"

#define PRECISION_z
 
#ifdef __cplusplus
extern "C" {
#endif

void magma_minproduct_ztrdtype1cbHLsym_withQ(
    magma_minproduct_int_t N, magma_minproduct_int_t NB,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t LDA,
    magma_minproductDoubleComplex *V, magma_minproductDoubleComplex *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz);

void magma_minproduct_ztrdtype2cbHLsym_withQ(
    magma_minproduct_int_t N, magma_minproduct_int_t NB,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t LDA,
    magma_minproductDoubleComplex *V, magma_minproductDoubleComplex *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz);
   
void magma_minproduct_ztrdtype3cbHLsym_withQ(
    magma_minproduct_int_t N, magma_minproduct_int_t NB,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t LDA,
    magma_minproductDoubleComplex *V, magma_minproductDoubleComplex *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz);

void magma_minproduct_zlarfxsym(
    magma_minproduct_int_t N,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t LDA,
    magma_minproductDoubleComplex *V, magma_minproductDoubleComplex *TAU);

#ifdef __cplusplus
}
#endif


////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C" void
magma_minproduct_zlarfxsym(
    magma_minproduct_int_t N,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t LDA,
    magma_minproductDoubleComplex *V, magma_minproductDoubleComplex *TAU)
{
    magma_minproduct_int_t IONE=1;
    magma_minproductDoubleComplex dtmp;
    magma_minproductDoubleComplex Z_ZERO =  MAGMA_minproduct_Z_ZERO;
    //magma_minproductDoubleComplex Z_ONE  =  MAGMA_minproduct_Z_ONE;
    magma_minproductDoubleComplex Z_MONE =  MAGMA_minproduct_Z_NEG_ONE;
    magma_minproductDoubleComplex Z_HALF =  MAGMA_minproduct_Z_HALF;
    //magma_minproductDoubleComplex WORK[N];
    magma_minproductDoubleComplex *WORK;
    magma_minproduct_zmalloc_cpu( &WORK, N );
    
    /* apply left and right on A(st:ed,st:ed)*/
    //magma_minproduct_zlarfxsym(len,A(st,st),LDX,V(st),TAU(st));
    /* X = AVtau */
    blasf77_zhemv("L",&N, TAU, A, &LDA, V, &IONE, &Z_ZERO, WORK, &IONE);
    /* je calcul dtmp= X'*V */
    dtmp = magma_minproduct_cblas_zdotc(N, WORK, IONE, V, IONE);
    /* je calcul 1/2 X'*V*t = 1/2*dtmp*tau  */
    dtmp = -dtmp * Z_HALF * (*TAU);
    /* je calcul W=X-1/2VX'Vt = X - dtmp*V */
    /*
    for (j = 0; j < N; j++)
        WORK[j] = WORK[j] + (dtmp*V[j]); */
    blasf77_zaxpy(&N, &dtmp, V, &IONE, WORK, &IONE);
    /* performs the symmetric rank 2 operation A := alpha*x*y' + alpha*y*x' + A */
    blasf77_zher2("L",&N,&Z_MONE,WORK,&IONE,V,&IONE,A,&LDA);
    
    magma_minproduct_free_cpu(WORK);
}
////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
//                  TYPE 1-BAND Householder
///////////////////////////////////////////////////////////
//// add -1 because of C
#define A(m,n)   &(A[((m)-(n)) + LDA*((n)-1)])
#define V(m)     &(V[(m)])
#define TAU(m)   &(TAU[(m)])
extern "C" void
magma_minproduct_ztrdtype1cbHLsym_withQ(
    magma_minproduct_int_t N, magma_minproduct_int_t NB,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t LDA,
    magma_minproductDoubleComplex *V, magma_minproductDoubleComplex *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz)
{
    //magma_minproduct_int_t    J1, J2, J3, i, j;
    magma_minproduct_int_t    len, LDX;
    magma_minproduct_int_t    IONE=1;
    magma_minproduct_int_t    blkid, vpos, taupos, tpos;
    //magma_minproductDoubleComplex conjtmp;
    magma_minproductDoubleComplex Z_ONE  =  MAGMA_minproduct_Z_ONE;
    magma_minproductDoubleComplex *WORK;
    magma_minproduct_zmalloc_cpu( &WORK, N );
    
    
    findVTpos(N,NB,Vblksiz,sweep-1,st-1, &vpos, &taupos, &tpos, &blkid);
    //printf("voici vpos %d taupos %d  tpos %d  blkid %d \n", vpos, taupos, tpos, blkid);
    LDX     = LDA-1;
    len     = ed-st+1;
    *V(vpos)  = Z_ONE;
    memcpy(V(vpos+1), A(st+1, st-1), (len-1)*sizeof(magma_minproductDoubleComplex));
    memset(A(st+1, st-1), 0, (len-1)*sizeof(magma_minproductDoubleComplex));
    /* Eliminate the col  at st-1 */
    lapackf77_zlarfg( &len, A(st, st-1), V(vpos+1), &IONE, TAU(taupos) );
    /* apply left and right on A(st:ed,st:ed)*/
    magma_minproduct_zlarfxsym(len,A(st,st),LDX,V(vpos),TAU(taupos));
    //conjtmp = MAGMA_minproduct_Z_CNJG(*TAU(taupos));
    //lapackf77_zlarfy("L", &len, V(vpos), &IONE, &conjtmp, A(st,st), &LDX, WORK); //&(MAGMA_minproduct_Z_CNJG(*TAU(taupos)))
    magma_minproduct_free_cpu(WORK);
}
#undef A
#undef V
#undef TAU
///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
//                  TYPE 1-LPK Householder
///////////////////////////////////////////////////////////
//// add -1 because of C
#define A(m,n)   &(A[((m)-(n)) + LDA*((n)-1)])
#define V(m)     &(V[(m)])
#define TAU(m)   &(TAU[(m)])
extern "C" void
magma_minproduct_ztrdtype2cbHLsym_withQ(
    magma_minproduct_int_t N, magma_minproduct_int_t NB,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t LDA,
    magma_minproductDoubleComplex *V, magma_minproductDoubleComplex *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz)
{
    magma_minproduct_int_t    J1, J2, len, lem, LDX;
    //magma_minproduct_int_t    i, j;
    magma_minproduct_int_t    IONE=1;
    magma_minproduct_int_t    blkid, vpos, taupos, tpos;
    magma_minproductDoubleComplex conjtmp;
    magma_minproductDoubleComplex Z_ONE  =  MAGMA_minproduct_Z_ONE;
    //magma_minproductDoubleComplex WORK[NB];
    magma_minproductDoubleComplex *WORK;
    magma_minproduct_zmalloc_cpu( &WORK, NB );
    
    
    findVTpos(N,NB,Vblksiz,sweep-1,st-1, &vpos, &taupos, &tpos, &blkid);
    LDX    = LDA-1;
    J1     = ed+1;
    J2     = min(ed+NB,N);
    len    = ed-st+1;
    lem    = J2-J1+1;
    if (lem > 0) {
        /* apply remaining right commming from the top block */
        lapackf77_zlarfx("R", &lem, &len, V(vpos), TAU(taupos), A(J1, st), &LDX, WORK);
    }
    if (lem > 1) {
        findVTpos(N,NB,Vblksiz,sweep-1,J1-1, &vpos, &taupos, &tpos, &blkid);
        /* remove the first column of the created bulge */
        *V(vpos)  = Z_ONE;
        memcpy(V(vpos+1), A(J1+1, st), (lem-1)*sizeof(magma_minproductDoubleComplex));
        memset(A(J1+1, st),0,(lem-1)*sizeof(magma_minproductDoubleComplex));
        /* Eliminate the col at st */
        lapackf77_zlarfg( &lem, A(J1, st), V(vpos+1), &IONE, TAU(taupos) );
        /* apply left on A(J1:J2,st+1:ed) */
        len = len-1; /* because we start at col st+1 instead of st. col st is the col that has been revomved; */
        conjtmp = MAGMA_minproduct_Z_CNJG(*TAU(taupos));
        lapackf77_zlarfx("L", &lem, &len, V(vpos),  &conjtmp, A(J1, st+1), &LDX, WORK);
    }
    magma_minproduct_free_cpu(WORK);
}
#undef A
#undef V
#undef TAU
///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
//                  TYPE 1-LPK Householder
///////////////////////////////////////////////////////////
//// add -1 because of C
#define A(m,n)   &(A[((m)-(n)) + LDA*((n)-1)])
#define V(m)     &(V[(m)])
#define TAU(m)   &(TAU[(m)])
extern "C" void
magma_minproduct_ztrdtype3cbHLsym_withQ(
    magma_minproduct_int_t N, magma_minproduct_int_t NB,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t LDA,
    magma_minproductDoubleComplex *V, magma_minproductDoubleComplex *TAU,
    magma_minproduct_int_t st, magma_minproduct_int_t ed, magma_minproduct_int_t sweep, magma_minproduct_int_t Vblksiz)
{
    //magma_minproduct_int_t    J1, J2, J3, i, j;
    magma_minproduct_int_t    len, LDX;
    //magma_minproduct_int_t    IONE=1;
    magma_minproduct_int_t    blkid, vpos, taupos, tpos;
    //magma_minproductDoubleComplex conjtmp;
    magma_minproductDoubleComplex *WORK;
    magma_minproduct_zmalloc_cpu( &WORK, N );
    
    
    findVTpos(N,NB,Vblksiz,sweep-1,st-1, &vpos, &taupos, &tpos, &blkid);
    LDX    = LDA-1;
    len    = ed-st+1;
    
    /* apply left and right on A(st:ed,st:ed)*/
    magma_minproduct_zlarfxsym(len,A(st,st),LDX,V(vpos),TAU(taupos));
    //conjtmp = MAGMA_minproduct_Z_CNJG(*TAU(taupos));
    //lapackf77_zlarfy("L", &len, V(vpos), &IONE,  &(MAGMA_minproduct_Z_CNJG(*TAU(taupos))), A(st,st), &LDX, WORK);
    magma_minproduct_free_cpu(WORK);
}
#undef A
#undef V
#undef TAU
///////////////////////////////////////////////////////////
