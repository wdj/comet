/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @author Mark Gates
       @precisions normal z -> s d c
*/
#include "common_magma_tally2.h"

#define PRECISION_z


//===========================================================================
// Set a matrix from CPU to multi-GPUs in 1D column block cyclic distribution.
// The dA arrays are pointers to the matrix data on the corresponding GPUs.
//===========================================================================
extern "C" void
magma_tally2_zsetmatrix_1D_col_bcyclic(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const magma_tally2DoubleComplex    *hA,   magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr       dA[], magma_tally2_int_t ldda,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb )
{
    magma_tally2_int_t info = 0;
    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( lda < m )
        info = -4;
    else if ( ldda < m )
        info = -6;
    else if ( ngpu < 1 )
        info = -7;
    else if ( nb < 1 )
        info = -8;
    
    if (info != 0) {
        magma_tally2_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    magma_tally2_int_t j, dev, jb;
    magma_tally2_device_t cdevice;

    magma_tally2_getdevice( &cdevice );

    for( j = 0; j < n; j += nb ) {
        dev = (j/nb) % ngpu;
        magma_tally2_setdevice( dev );
        jb = min(nb, n-j);
        magma_tally2_zsetmatrix_async( m, jb,
                                hA + j*lda, lda,
                                dA[dev] + j/(nb*ngpu)*nb*ldda, ldda, NULL );
    }

    magma_tally2_setdevice( cdevice );
}


//===========================================================================
// Get a matrix with 1D column block cyclic distribution from multi-GPUs to the CPU.
// The dA arrays are pointers to the matrix data on the corresponding GPUs.
//===========================================================================
extern "C" void
magma_tally2_zgetmatrix_1D_col_bcyclic(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr const dA[], magma_tally2_int_t ldda,
    magma_tally2DoubleComplex                *hA,   magma_tally2_int_t lda,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb )
{
    magma_tally2_int_t info = 0;
    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( ldda < m )
        info = -4;
    else if ( lda < m )
        info = -6;
    else if ( ngpu < 1 )
        info = -7;
    else if ( nb < 1 )
        info = -8;
    
    if (info != 0) {
        magma_tally2_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    magma_tally2_int_t j, dev, jb;
    magma_tally2_device_t cdevice;

    magma_tally2_getdevice( &cdevice );

    for( j = 0; j < n; j += nb ) {
        dev = (j/nb) % ngpu;
        magma_tally2_setdevice( dev );
        jb = min(nb, n-j);
        magma_tally2_zgetmatrix_async( m, jb,
                                dA[dev] + j/(nb*ngpu)*nb*ldda, ldda,
                                hA + j*lda, lda, NULL );
    }

    magma_tally2_setdevice( cdevice );
}


//===========================================================================
// Set a matrix from CPU to multi-GPUs in 1D row block cyclic distribution.
// The dA arrays are pointers to the matrix data on the corresponding GPUs.
//===========================================================================
extern "C" void
magma_tally2_zsetmatrix_1D_row_bcyclic(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const magma_tally2DoubleComplex    *hA,   magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr       dA[], magma_tally2_int_t ldda,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb )
{
    magma_tally2_int_t info = 0;
    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( lda < m )
        info = -4;
    else if ( ldda < (1+m/(nb*ngpu))*nb )
        info = -6;
    else if ( ngpu < 1 )
        info = -7;
    else if ( nb < 1 )
        info = -8;
    
    if (info != 0) {
        magma_tally2_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    magma_tally2_int_t i, dev, jb;
    magma_tally2_device_t cdevice;

    magma_tally2_getdevice( &cdevice );

    for( i = 0; i < m; i += nb ) {
        dev = (i/nb) % ngpu;
        magma_tally2_setdevice( dev );
        jb = min(nb, m-i);
        magma_tally2_zsetmatrix_async( jb, n,
                                hA + i, lda,
                                dA[dev] + i/(nb*ngpu)*nb, ldda, NULL );
    }

    magma_tally2_setdevice( cdevice );
}


//===========================================================================
// Get a matrix with 1D row block cyclic distribution from multi-GPUs to the CPU.
// The dA arrays are pointers to the matrix data for the corresponding GPUs.
//===========================================================================
extern "C" void
magma_tally2_zgetmatrix_1D_row_bcyclic(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_const_ptr const dA[], magma_tally2_int_t ldda,
    magma_tally2DoubleComplex                *hA,   magma_tally2_int_t lda,
    magma_tally2_int_t ngpu, magma_tally2_int_t nb )
{
    magma_tally2_int_t info = 0;
    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( ldda < (1+m/(nb*ngpu))*nb )
        info = -4;
    else if ( lda < m )
        info = -6;
    else if ( ngpu < 1 )
        info = -7;
    else if ( nb < 1 )
        info = -8;
    
    if (info != 0) {
        magma_tally2_xerbla( __func__, -(info) );
        return;  //info;
    }
    
    magma_tally2_int_t i, dev, jb;
    magma_tally2_device_t cdevice;

    magma_tally2_getdevice( &cdevice );

    for( i = 0; i < m; i += nb ) {
        dev = (i/nb) % ngpu;
        magma_tally2_setdevice( dev );
        jb = min(nb, m-i);
        magma_tally2_zgetmatrix_async( jb, n,
                                dA[dev] + i/(nb*ngpu)*nb, ldda,
                                hA + i, lda, NULL );
    }

    magma_tally2_setdevice( cdevice );
}