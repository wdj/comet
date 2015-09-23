/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    magma_minproduct_xerbla is an error handler for the MAGMA_minproduct routines.
    It is called by a MAGMA_minproduct routine if an input parameter has an
    invalid value. It prints an error message.

    Installers may consider modifying it to
    call system-specific exception-handling facilities.

    Arguments
    ---------
    @param[in]
    srname  CHAR*
            The name of the subroutine that called XERBLA.
            In C/C++ it is convenient to use "__func__".

    @param[in]
    minfo   INTEGER
            Note minfo's sign is opposite info's normal sign.
            
            Normally:
            - minfo > 0: The position of the invalid parameter
                         in the parameter list of the calling routine.
            
            These conditions are also reported, but normally code should not
            call xerbla for these runtime errors:
            - minfo <  0:          Function-specific error.
            - minfo >= -MAGMA_minproduct_ERR: Pre-defined MAGMA_minproduct error, such as malloc failure.
            - minfo == 0:          No error.

    @ingroup magma_minproduct_util
    ********************************************************************/
extern "C"
void magma_minproduct_xerbla(const char *srname , magma_minproduct_int_t minfo)
{
    // the first 3 cases are unusual for calling xerbla;
    // normally runtime errors are passed back in info.
    if ( minfo < 0 ) {
        fprintf( stderr, "Error in %s, function-specific error (info = %d)\n",
                 srname, (int) -minfo );
    }
    else if ( minfo == 0 ) {
        fprintf( stderr, "No error, why is %s calling xerbla? (info = %d)\n",
                 srname, (int) -minfo );
    }
    else if ( minfo >= -MAGMA_minproduct_ERR ) {
        fprintf( stderr, "Error in %s, %s (info = %d)\n",
                 srname, magma_minproduct_strerror(-minfo), (int) -minfo );
    }
    else {
        // this is the normal case for calling xerbla;
        // invalid parameter values are usually logic errors, not runtime errors.
        fprintf( stderr, "On entry to %s, parameter %d had an illegal value (info = %d)\n",
                 srname, (int) minfo, (int) -minfo );
    }
}
