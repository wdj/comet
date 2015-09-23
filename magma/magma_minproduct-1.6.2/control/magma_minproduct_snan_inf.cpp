/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
       @generated from magma_minproduct_znan_inf.cpp normal z -> s, Fri Jan 30 19:00:20 2015

*/
#include <limits>

#include "common_magma_minproduct.h"

#define REAL


const float MAGMA_minproduct_S_NAN
    = MAGMA_minproduct_S_MAKE( std::numeric_limits<float>::quiet_NaN(),
                    std::numeric_limits<float>::quiet_NaN() );

const float MAGMA_minproduct_S_INF
    = MAGMA_minproduct_S_MAKE( std::numeric_limits<float>::infinity(),
                    std::numeric_limits<float>::infinity() );

/** @return true if either real(x) or imag(x) is NAN. */
inline bool magma_minproduct_s_isnan( float x )
{
#ifdef COMPLEX
    return isnan( MAGMA_minproduct_S_REAL( x )) ||
           isnan( MAGMA_minproduct_S_IMAG( x ));
#else
    return isnan( x );
#endif
}


/** @return true if either real(x) or imag(x) is INF. */
inline bool magma_minproduct_s_isinf( float x )
{
#ifdef COMPLEX
    return isinf( MAGMA_minproduct_S_REAL( x )) ||
           isinf( MAGMA_minproduct_S_IMAG( x ));
#else
    return isinf( x );
#endif
}


/**
    Purpose
    -------

    magma_minproduct_snan_inf checks a matrix that is located on the CPU host
    for NAN (not-a-number) and INF (infinity) values.
    
    NAN is created by 0/0 and similar.
    INF is created by x/0 and similar, where x != 0.

    Arguments
    ---------
    @param[in]
    uplo    magma_minproduct_uplo_t
            Specifies what part of the matrix A to check.
      -     = Magma_minproductUpper:  Upper triangular part of A
      -     = Magma_minproductLower:  Lower triangular part of A
      -     = Magma_minproductFull:   All of A

    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in]
    A       REAL array, dimension (LDA,N), on the CPU host.
            The M-by-N matrix to be printed.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    @param[out]
    cnt_nan INTEGER*
            If non-NULL, on exit contains the number of NAN values in A.

    @param[out]
    cnt_inf INTEGER*
            If non-NULL, on exit contains the number of INF values in A.
            
    @return
      -     >= 0:  Returns number of NAN + number of INF values.
      -     <  0:  If it returns -i, the i-th argument had an illegal value,
                   or another error occured, such as memory allocation failed.

    @ingroup magma_minproduct_saux2
    ********************************************************************/
extern "C"
magma_minproduct_int_t magma_minproduct_snan_inf(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    const float *A, magma_minproduct_int_t lda,
    magma_minproduct_int_t *cnt_nan,
    magma_minproduct_int_t *cnt_inf )
{
    #define A(i,j) (A + (i) + (j)*lda)
    
    magma_minproduct_int_t info = 0;
    if ( uplo != Magma_minproductLower && uplo != Magma_minproductUpper && uplo != Magma_minproductFull )
        info = -1;
    else if ( m < 0 )
        info = -2;
    else if ( n < 0 )
        info = -3;
    else if ( lda < max(1,m) )
        info = -5;
    
    if (info != 0) {
        magma_minproduct_xerbla( __func__, -(info) );
        return info;
    }
    
    int c_nan = 0;
    int c_inf = 0;
    
    if ( uplo == Magma_minproductLower ) {
        for( int j = 0; j < n; ++j ) {
            for( int i = j; i < m; ++i ) {  // i >= j
                if      ( magma_minproduct_s_isnan( *A(i,j) )) { c_nan++; }
                else if ( magma_minproduct_s_isinf( *A(i,j) )) { c_inf++; }
            }
        }
    }
    else if ( uplo == Magma_minproductUpper ) {
        for( int j = 0; j < n; ++j ) {
            for( int i = 0; i < m && i <= j; ++i ) {  // i <= j
                if      ( magma_minproduct_s_isnan( *A(i,j) )) { c_nan++; }
                else if ( magma_minproduct_s_isinf( *A(i,j) )) { c_inf++; }
            }
        }
    }
    else if ( uplo == Magma_minproductFull ) {
        for( int j = 0; j < n; ++j ) {
            for( int i = 0; i < m; ++i ) {
                if      ( magma_minproduct_s_isnan( *A(i,j) )) { c_nan++; }
                else if ( magma_minproduct_s_isinf( *A(i,j) )) { c_inf++; }
            }
        }
    }
    
    if ( cnt_nan != NULL ) { *cnt_nan = c_nan; }
    if ( cnt_inf != NULL ) { *cnt_inf = c_inf; }
    
    return (c_nan + c_inf);
}


/**
    Purpose
    -------

    magma_minproduct_snan_inf checks a matrix that is located on the CPU host
    for NAN (not-a-number) and INF (infinity) values.
    
    NAN is created by 0/0 and similar.
    INF is created by x/0 and similar, where x != 0.

    Arguments
    ---------
    @param[in]
    uplo    magma_minproduct_uplo_t
            Specifies what part of the matrix A to check.
      -     = Magma_minproductUpper:  Upper triangular part of A
      -     = Magma_minproductLower:  Lower triangular part of A
      -     = Magma_minproductFull:   All of A

    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in]
    dA      REAL array, dimension (LDDA,N), on the GPU device.
            The M-by-N matrix to be printed.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

    @param[out]
    cnt_nan INTEGER*
            If non-NULL, on exit contains the number of NAN values in A.

    @param[out]
    cnt_inf INTEGER*
            If non-NULL, on exit contains the number of INF values in A.
            
    @return
      -     >= 0:  Returns number of NAN + number of INF values.
      -     <  0:  If it returns -i, the i-th argument had an illegal value,
                   or another error occured, such as memory allocation failed.

    @ingroup magma_minproduct_saux2
    ********************************************************************/
extern "C"
magma_minproduct_int_t magma_minproduct_snan_inf_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproduct_int_t *cnt_nan,
    magma_minproduct_int_t *cnt_inf )
{
    magma_minproduct_int_t info = 0;
    if ( uplo != Magma_minproductLower && uplo != Magma_minproductUpper && uplo != Magma_minproductFull )
        info = -1;
    else if ( m < 0 )
        info = -2;
    else if ( n < 0 )
        info = -3;
    else if ( ldda < max(1,m) )
        info = -5;
    
    if (info != 0) {
        magma_minproduct_xerbla( __func__, -(info) );
        return info;
    }
    
    magma_minproduct_int_t lda = m;
    float* A;
    magma_minproduct_smalloc_cpu( &A, lda*n );
    magma_minproduct_sgetmatrix( m, n, dA, ldda, A, lda );
    
    magma_minproduct_int_t cnt = magma_minproduct_snan_inf( uplo, m, n, A, lda, cnt_nan, cnt_inf );
    
    magma_minproduct_free_cpu( A );
    return cnt;
}
