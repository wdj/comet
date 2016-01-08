/*
    -- MAGMA_tally4 (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date
       
       @author Mark Gates
*/

#include <complex>

#include "magma_tally4_operators.h"

/**
    @return Complex square root of x.
    
    @param[in]
    x       COMPLEX_16
    
    @ingroup magma_tally4_zaux0
    ********************************************************************/
magma_tally4DoubleComplex magma_tally4_zsqrt( magma_tally4DoubleComplex x )
{
    std::complex<double> y = std::sqrt( std::complex<double>( real(x), imag(x) ));
    return MAGMA_tally4_Z_MAKE( real(y), imag(y) );
}


/**
    @return Complex square root of x.
    
    @param[in]
    x       COMPLEX
    
    @ingroup magma_tally4_caux0
    ********************************************************************/
magma_tally4FloatComplex magma_tally4_csqrt( magma_tally4FloatComplex x )
{
    std::complex<float> y = std::sqrt( std::complex<float>( real(x), imag(x) ));
    return MAGMA_tally4_C_MAKE( real(y), imag(y) );
}
