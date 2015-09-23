/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#ifndef MAGMA_minproduct_BULGE_H
#define MAGMA_minproduct_BULGE_H

#include "magma_minproduct_types.h"

#ifdef __cplusplus
extern "C" {
#endif

    magma_minproduct_int_t magma_minproduct_bulge_get_nb(magma_minproduct_int_t n);

    void cmp_vals(int n, double *wr1, double *wr2, double *nrmI, double *nrm1, double *nrm2);

    void magma_minproduct_bulge_findVTAUpos(magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t Vblksiz, magma_minproduct_int_t sweep, magma_minproduct_int_t st, magma_minproduct_int_t ldv,
                                 magma_minproduct_int_t *Vpos, magma_minproduct_int_t *TAUpos);

    void magma_minproduct_bulge_findVTpos(magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t Vblksiz, magma_minproduct_int_t sweep, magma_minproduct_int_t st, magma_minproduct_int_t ldv, magma_minproduct_int_t ldt,
                               magma_minproduct_int_t *Vpos, magma_minproduct_int_t *Tpos);

    void magma_minproduct_bulge_findVTAUTpos(magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t Vblksiz, magma_minproduct_int_t sweep, magma_minproduct_int_t st, magma_minproduct_int_t ldv, magma_minproduct_int_t ldt,
                                  magma_minproduct_int_t *Vpos, magma_minproduct_int_t *TAUpos, magma_minproduct_int_t *Tpos, magma_minproduct_int_t *blkid);

    void magma_minproduct_bulge_findpos(magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t Vblksiz, magma_minproduct_int_t sweep, magma_minproduct_int_t st, magma_minproduct_int_t *myblkid);
    void magma_minproduct_bulge_findpos113(magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t Vblksiz, magma_minproduct_int_t sweep, magma_minproduct_int_t st, magma_minproduct_int_t *myblkid);

    magma_minproduct_int_t magma_minproduct_bulge_get_blkcnt(magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t Vblksiz);

#ifdef __cplusplus
}
#endif

#endif  // MAGMA_minproduct_BULGE_H
