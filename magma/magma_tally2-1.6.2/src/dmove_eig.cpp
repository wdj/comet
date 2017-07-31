/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Azzam Haidar

       @precisions normal d -> s

*/
#include "common_magma_tally2.h"

extern "C" void
magma_tally2_dmove_eig(
    magma_tally2_range_t range, magma_tally2_int_t n, double *w, magma_tally2_int_t *il,
    magma_tally2_int_t *iu, double vl, double vu, magma_tally2_int_t *m)
{
    magma_tally2_int_t valeig, indeig, i;

    valeig = (range == Magma_tally2RangeV);
    indeig = (range == Magma_tally2RangeI);

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
