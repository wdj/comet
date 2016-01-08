#include "magma_tally4.h"
#include "magma_tally4_mangling.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#define magma_tally4f_init FORTRAN_NAME( magma_tally4f_init, MAGMA_tally4F_INIT )
void magma_tally4f_init( void )
{
    magma_tally4_init();
}

#define magma_tally4f_finalize FORTRAN_NAME( magma_tally4f_finalize, MAGMA_tally4F_FINALIZE )
void magma_tally4f_finalize( void )
{
    magma_tally4_finalize();
}

#ifdef __cplusplus
}
#endif
