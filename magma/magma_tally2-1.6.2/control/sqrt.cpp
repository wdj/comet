/*
    -- MAGMA_tally2 (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date
       
       @author Mark Gates
*/

#include <complex>

#include "magma_tally2_operators.h"

/**
    @return Complex square root of x.
    
    @param[in]
    x       COMPLEX_16
    
    @ingroup magma_tally2_zaux0
    ********************************************************************/
magma_tally2DoubleComplex magma_tally2_zsqrt( magma_tally2DoubleComplex x )
{
    std::complex<double> y = std::sqrt( std::complex<double>( real(x), imag(x) ));
    return MAGMA_tally2_Z_MAKE( real(y), imag(y) );
}


/**
    @return Complex square root of x.
    
    @param[in]
    x       COMPLEX
    
    @ingroup magma_tally2_caux0
    ********************************************************************/
magma_tally2FloatComplex magma_tally2_csqrt( magma_tally2FloatComplex x )
{
    std::complex<float> y = std::sqrt( std::complex<float>( real(x), imag(x) ));
    return MAGMA_tally2_C_MAKE( real(y), imag(y) );
}
