/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#ifndef MAGMA_tally4_BULGE_H
#define MAGMA_tally4_BULGE_H

#include "magma_tally4_types.h"

#ifdef __cplusplus
extern "C" {
#endif

    magma_tally4_int_t magma_tally4_bulge_get_nb(magma_tally4_int_t n);

    void cmp_vals(int n, double *wr1, double *wr2, double *nrmI, double *nrm1, double *nrm2);

    void magma_tally4_bulge_findVTAUpos(magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t Vblksiz, magma_tally4_int_t sweep, magma_tally4_int_t st, magma_tally4_int_t ldv,
                                 magma_tally4_int_t *Vpos, magma_tally4_int_t *TAUpos);

    void magma_tally4_bulge_findVTpos(magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t Vblksiz, magma_tally4_int_t sweep, magma_tally4_int_t st, magma_tally4_int_t ldv, magma_tally4_int_t ldt,
                               magma_tally4_int_t *Vpos, magma_tally4_int_t *Tpos);

    void magma_tally4_bulge_findVTAUTpos(magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t Vblksiz, magma_tally4_int_t sweep, magma_tally4_int_t st, magma_tally4_int_t ldv, magma_tally4_int_t ldt,
                                  magma_tally4_int_t *Vpos, magma_tally4_int_t *TAUpos, magma_tally4_int_t *Tpos, magma_tally4_int_t *blkid);

    void magma_tally4_bulge_findpos(magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t Vblksiz, magma_tally4_int_t sweep, magma_tally4_int_t st, magma_tally4_int_t *myblkid);
    void magma_tally4_bulge_findpos113(magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t Vblksiz, magma_tally4_int_t sweep, magma_tally4_int_t st, magma_tally4_int_t *myblkid);

    magma_tally4_int_t magma_tally4_bulge_get_blkcnt(magma_tally4_int_t n, magma_tally4_int_t nb, magma_tally4_int_t Vblksiz);

#ifdef __cplusplus
}
#endif

#endif  // MAGMA_tally4_BULGE_H
