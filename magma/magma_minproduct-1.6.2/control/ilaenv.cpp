/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Mark Gates
*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    Wrapper around LAPACK ilaenv,
    this takes name and opts as C strings,
    then calls magma_minproductf77_ilaenv with them as character arrays with length, which
    then calls LAPACK ilaenv with them as character strings.
    It also converts from C to Fortran calling conventions (pointers).
    ********************************************************************/
magma_minproduct_int_t magma_minproduct_ilaenv(
    magma_minproduct_int_t ispec, const char* name, const char* opts,
    magma_minproduct_int_t n1, magma_minproduct_int_t n2, magma_minproduct_int_t n3, magma_minproduct_int_t n4 )
{
    magma_minproduct_int_t namelen = strlen( name );
    magma_minproduct_int_t optslen = strlen( opts );
    return magma_minproductf77_ilaenv( &ispec, name, &namelen, opts, &optslen, &n1, &n2, &n3, &n4 );
}
