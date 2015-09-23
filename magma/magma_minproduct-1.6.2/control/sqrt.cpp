/*
    -- MAGMA_minproduct (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date
       
       @author Mark Gates
*/

#include <complex>

#include "magma_minproduct_operators.h"

/**
    @return Complex square root of x.
    
    @param[in]
    x       COMPLEX_16
    
    @ingroup magma_minproduct_zaux0
    ********************************************************************/
magma_minproductDoubleComplex magma_minproduct_zsqrt( magma_minproductDoubleComplex x )
{
    std::complex<double> y = std::sqrt( std::complex<double>( real(x), imag(x) ));
    return MAGMA_minproduct_Z_MAKE( real(y), imag(y) );
}


/**
    @return Complex square root of x.
    
    @param[in]
    x       COMPLEX
    
    @ingroup magma_minproduct_caux0
    ********************************************************************/
magma_minproductFloatComplex magma_minproduct_csqrt( magma_minproductFloatComplex x )
{
    std::complex<float> y = std::sqrt( std::complex<float>( real(x), imag(x) ));
    return MAGMA_minproduct_C_MAKE( real(y), imag(y) );
}
