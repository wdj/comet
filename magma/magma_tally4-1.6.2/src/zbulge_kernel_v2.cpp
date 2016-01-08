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

#include "common_magma_tally4.h"
#include "magma_tally4_bulge.h"

#define PRECISION_z

inline static void
magma_tally4_zlarfxsym_v2(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *V, magma_tally4DoubleComplex *TAU,
    magma_tally4DoubleComplex *work);
///////////////////////////////////////////////////////////

inline static void
magma_tally4_zlarfxsym_v2(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *V, magma_tally4DoubleComplex *TAU,
    magma_tally4DoubleComplex *work)
{
/*
    WORK (workspace) double complex array, dimension N
*/

    magma_tally4_int_t ione = 1;
    magma_tally4DoubleComplex dtmp;
    magma_tally4DoubleComplex c_zero   =  MAGMA_tally4_Z_ZERO;
    magma_tally4DoubleComplex c_neg_one=  MAGMA_tally4_Z_NEG_ONE;
    magma_tally4DoubleComplex c_half   =  MAGMA_tally4_Z_HALF;

    /* X = AVtau */
    blasf77_zhemv("L",&n, TAU, A, &lda, V, &ione, &c_zero, work, &ione);

    /* compute dtmp= X'*V */
    dtmp = magma_tally4_cblas_zdotc(n, work, ione, V, ione);

    /* compute 1/2 X'*V*t = 1/2*dtmp*tau  */
    dtmp = -dtmp * c_half * (*TAU);

    /* compute W=X-1/2VX'Vt = X - dtmp*V */
    blasf77_zaxpy(&n, &dtmp, V, &ione, work, &ione);

    /* performs the symmetric rank 2 operation A := alpha*x*y' + alpha*y*x' + A */
    blasf77_zher2("L", &n, &c_neg_one, work, &ione, V, &ione, A, &lda);
}

///////////////////////////////////////////////////////////
//                  TYPE 1-BAND Householder
///////////////////////////////////////////////////////////
//// add -1 because of C
#define A(i,j)   &(A[((i)-(j)) + lda*((j)-1)])
#define V(i)     &(V[(i)])
#define TAU(i)   &(TAU[(i)])
extern "C" void
magma_tally4_ztrdtype1cbHLsym_withQ_v2(
    magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *V, magma_tally4_int_t ldv,
    magma_tally4DoubleComplex *TAU,
    magma_tally4_int_t st, magma_tally4_int_t ed,
    magma_tally4_int_t sweep, magma_tally4_int_t Vblksiz,
    magma_tally4DoubleComplex *work)
{
/*
    WORK (workspace) double complex array, dimension N
*/

    magma_tally4_int_t ione = 1;
    magma_tally4_int_t vpos, taupos, len, len2;

    magma_tally4DoubleComplex c_one    =  MAGMA_tally4_Z_ONE;

    magma_tally4_bulge_findVTAUpos(n, nb, Vblksiz, sweep-1, st-1, ldv, &vpos, &taupos);
    //printf("voici vpos %d taupos %d  tpos %d  blkid %d \n", vpos, taupos, tpos, blkid);

    len     = ed-st+1;
    *V(vpos)  = c_one;

    len2 = len-1;
    blasf77_zcopy( &len2, A(st+1, st-1), &ione, V(vpos+1), &ione );
    //memcpy(V(vpos+1), A(st+1, st-1), (len-1)*sizeof(magma_tally4DoubleComplex));
    memset(A(st+1, st-1), 0, (len-1)*sizeof(magma_tally4DoubleComplex));

    /* Eliminate the col  at st-1 */
    lapackf77_zlarfg( &len, A(st, st-1), V(vpos+1), &ione, TAU(taupos) );

    /* apply left and right on A(st:ed,st:ed)*/
    magma_tally4_zlarfxsym_v2(len, A(st,st), lda-1, V(vpos), TAU(taupos), work);
}
#undef A
#undef V
#undef TAU

///////////////////////////////////////////////////////////
//                  TYPE 1-LPK Householder
///////////////////////////////////////////////////////////
//// add -1 because of C
#define A(i,j)   &(A[((i)-(j)) + lda*((j)-1)])
#define V(i)     &(V[(i)])
#define TAU(i)   &(TAU[(i)])
extern "C" void
magma_tally4_ztrdtype2cbHLsym_withQ_v2(
    magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *V, magma_tally4_int_t ldv,
    magma_tally4DoubleComplex *TAU,
    magma_tally4_int_t st, magma_tally4_int_t ed,
    magma_tally4_int_t sweep, magma_tally4_int_t Vblksiz,
    magma_tally4DoubleComplex *work)
{
    /*
     WORK (workspace) double complex array, dimension NB
    */

    magma_tally4_int_t ione = 1;
    magma_tally4_int_t vpos, taupos;

    magma_tally4DoubleComplex conjtmp;

    magma_tally4DoubleComplex c_one = MAGMA_tally4_Z_ONE;

    magma_tally4_int_t ldx = lda-1;
    magma_tally4_int_t len = ed - st + 1;
    magma_tally4_int_t lem = min(ed+nb, n) - ed;
    magma_tally4_int_t lem2;
    
    if (lem > 0) {
        magma_tally4_bulge_findVTAUpos(n, nb, Vblksiz, sweep-1, st-1, ldv, &vpos, &taupos);
        /* apply remaining right coming from the top block */
        lapackf77_zlarfx("R", &lem, &len, V(vpos), TAU(taupos), A(ed+1, st), &ldx, work);
    }
    if (lem > 1) {
        magma_tally4_bulge_findVTAUpos(n, nb, Vblksiz, sweep-1, ed, ldv, &vpos, &taupos);

        /* remove the first column of the created bulge */
        *V(vpos)  = c_one;
        //memcpy(V(vpos+1), A(ed+2, st), (lem-1)*sizeof(magma_tally4DoubleComplex));
        lem2 = lem-1;
        blasf77_zcopy( &lem2, A(ed+2, st), &ione, V(vpos+1), &ione );
        memset(A(ed+2, st),0,(lem-1)*sizeof(magma_tally4DoubleComplex));

        /* Eliminate the col at st */
        lapackf77_zlarfg( &lem, A(ed+1, st), V(vpos+1), &ione, TAU(taupos) );

        /* apply left on A(J1:J2,st+1:ed) */
        len = len-1; /* because we start at col st+1 instead of st. col st is the col that has been removed; */
        conjtmp = MAGMA_tally4_Z_CNJG(*TAU(taupos));
        lapackf77_zlarfx("L", &lem, &len, V(vpos),  &conjtmp, A(ed+1, st+1), &ldx, work);
    }
}
#undef A
#undef V
#undef TAU

///////////////////////////////////////////////////////////
//                  TYPE 1-LPK Householder
///////////////////////////////////////////////////////////
//// add -1 because of C
#define A(i,j)   &(A[((i)-(j)) + lda*((j)-1)])
#define V(i)     &(V[(i)])
#define TAU(i)   &(TAU[(i)])
extern "C" void
magma_tally4_ztrdtype3cbHLsym_withQ_v2(
    magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *V, magma_tally4_int_t ldv,
    magma_tally4DoubleComplex *TAU,
    magma_tally4_int_t st, magma_tally4_int_t ed,
    magma_tally4_int_t sweep, magma_tally4_int_t Vblksiz,
    magma_tally4DoubleComplex *work)
{
    /*
     WORK (workspace) double complex array, dimension N
     */

    magma_tally4_int_t vpos, taupos;

    magma_tally4_bulge_findVTAUpos(n, nb, Vblksiz, sweep-1, st-1, ldv, &vpos, &taupos);

    magma_tally4_int_t len = ed-st+1;

    /* apply left and right on A(st:ed,st:ed)*/
    magma_tally4_zlarfxsym_v2(len, A(st,st), lda-1, V(vpos), TAU(taupos), work);
}
#undef A
#undef V
#undef TAU
///////////////////////////////////////////////////////////
