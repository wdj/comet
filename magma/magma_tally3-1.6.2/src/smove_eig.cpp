/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Azzam Haidar

       @generated from dmove_eig.cpp normal d -> s, Fri Jan 30 19:00:17 2015

*/
#include "common_magma_tally3.h"

extern "C" void
magma_tally3_smove_eig(
    magma_tally3_range_t range, magma_tally3_int_t n, float *w, magma_tally3_int_t *il,
    magma_tally3_int_t *iu, float vl, float vu, magma_tally3_int_t *m)
{
    magma_tally3_int_t valeig, indeig, i;

    valeig = (range == Magma_tally3RangeV);
    indeig = (range == Magma_tally3RangeI);

    if (indeig) {
        *m = *iu - *il + 1;
        if (*il > 1)
            for (i = 0; i < *m; ++i)
                w[i] = w[*il - 1 + i];
    }
    else if (valeig) {
        *il=1;
        *iu=n;
        for (i = 0; i < n; ++i) {
            if (w[i] > vu) {
                *iu = i;
                break;
            }
            else if (w[i] < vl)
                ++*il;
            else if (*il > 1)
                w[i-*il+1]=w[i];
        }
        *m = *iu - *il + 1;
    }
    else {
        *il = 1;
        *iu = n;
        *m = n;
    }

    return;
}
