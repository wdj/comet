#include "magma_minproduct.h"
#include "magma_minproduct_mangling.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#define magma_minproductf_init FORTRAN_NAME( magma_minproductf_init, MAGMA_minproductF_INIT )
void magma_minproductf_init( void )
{
    magma_minproduct_init();
}

#define magma_minproductf_finalize FORTRAN_NAME( magma_minproductf_finalize, MAGMA_minproductF_FINALIZE )
void magma_minproductf_finalize( void )
{
    magma_minproduct_finalize();
}

#ifdef __cplusplus
}
#endif
