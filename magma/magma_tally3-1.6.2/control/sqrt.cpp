/*
    -- MAGMA_tally3 (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date
       
       @author Mark Gates
*/

#include <complex>

#include "magma_tally3_operators.h"

/**
    @return Complex square root of x.
    
    @param[in]
    x       COMPLEX_16
    
    @ingroup magma_tally3_zaux0
    ********************************************************************/
magma_tally3DoubleComplex magma_tally3_zsqrt( magma_tally3DoubleComplex x )
{
    std::complex<double> y = std::sqrt( std::complex<double>( real(x), imag(x) ));
    return MAGMA_tally3_Z_MAKE( real(y), imag(y) );
}


/**
    @return Complex square root of x.
    
    @param[in]
    x       COMPLEX
    
    @ingroup magma_tally3_caux0
    ********************************************************************/
magma_tally3FloatComplex magma_tally3_csqrt( magma_tally3FloatComplex x )
{
    std::complex<float> y = std::sqrt( std::complex<float>( real(x), imag(x) ));
    return MAGMA_tally3_C_MAKE( real(y), imag(y) );
}
