#include "magma_tally2.h"
#include "magma_tally2_mangling.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#define magma_tally2f_init FORTRAN_NAME( magma_tally2f_init, MAGMA_tally2F_INIT )
void magma_tally2f_init( void )
{
    magma_tally2_init();
}

#define magma_tally2f_finalize FORTRAN_NAME( magma_tally2f_finalize, MAGMA_tally2F_FINALIZE )
void magma_tally2f_finalize( void )
{
    magma_tally2_finalize();
}

#ifdef __cplusplus
}
#endif
