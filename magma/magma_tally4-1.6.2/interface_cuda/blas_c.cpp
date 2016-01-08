/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
 
       @author Mark Gates
       @generated from blas_z.cpp normal z -> c, Fri Jan 30 19:00:20 2015
*/

#include <stdlib.h>
#include <stdio.h>

#include "magma_tally4.h"
#include "error.h"

#define COMPLEX

#ifdef HAVE_CUBLAS

// ----------------------------------------
// Convert MAGMA_tally4 constants to CUBLAS v1 constants, which are the same as lapack.
// These must be static to avoid conflict with CUBLAS v2 translators.
#define cublas_trans_const( magma_tally4_const ) lapacke_trans_const( magma_tally4_const )
#define cublas_side_const(  magma_tally4_const ) lapacke_side_const(  magma_tally4_const )
#define cublas_diag_const(  magma_tally4_const ) lapacke_diag_const(  magma_tally4_const )
#define cublas_uplo_const(  magma_tally4_const ) lapacke_uplo_const(  magma_tally4_const )


// ========================================
// Level 1 BLAS

// --------------------
/** Returns index of element of vector x having max. absolute value;
    i.e., max (infinity) norm.
    
    @param[in]
    n       Number of elements in vector x. n >= 0.
            
    @param[in]
    dx      COMPLEX array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).
            
    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.
            
    @ingroup magma_tally4_cblas1
*/
extern "C" magma_tally4_int_t
magma_tally4_icamax(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dx, magma_tally4_int_t incx )
{
    return cublasIcamax( n, dx, incx );
}

// --------------------
/** Returns index of element of vector x having min. absolute value.

    @param[in]
    n       Number of elements in vector x. n >= 0.

    @param[in]
    dx      COMPLEX array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.

    @ingroup magma_tally4_cblas1
*/
extern "C" magma_tally4_int_t
magma_tally4_icamin(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dx, magma_tally4_int_t incx )
{
    return cublasIcamin( n, dx, incx );
}

// --------------------
/** Returns the sum of absolute values of vector x; i.e., one norm.

    @param[in]
    n       Number of elements in vector x. n >= 0.

    @param[in]
    dx      COMPLEX array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.

    @ingroup magma_tally4_cblas1
*/
extern "C" float
magma_tally4_scasum(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dx, magma_tally4_int_t incx )
{
    return cublasScasum( n, dx, incx );
}

// --------------------
/** Constant times a vector plus a vector; \f$ y = \alpha x + y \f$.

    @param[in]
    n       Number of elements in vectors x and y. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dx      COMPLEX array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dy      COMPLEX array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_tally4_cblas1
*/
extern "C" void
magma_tally4_caxpy(
    magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4FloatComplex_ptr       dy, magma_tally4_int_t incy )
{
    cublasCaxpy( n, alpha, dx, incx, dy, incy );
}

// --------------------
/** Copy vector x to vector y; \f$ y = x \f$.

    @param[in]
    n       Number of elements in vectors x and y. n >= 0.

    @param[in]
    dx      COMPLEX array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[out]
    dy      COMPLEX array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_tally4_cblas1
*/
extern "C" void
magma_tally4_ccopy(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4FloatComplex_ptr       dy, magma_tally4_int_t incy )
{
    cublasCcopy( n, dx, incx, dy, incy );
}

// --------------------
/** Returns dot product of vectors x and y; \f$ x^H y \f$.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in]
    dx      COMPLEX array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    dy      COMPLEX array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_tally4_cblas1
*/
extern "C"
magma_tally4FloatComplex magma_tally4_cdotc(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4FloatComplex_const_ptr dy, magma_tally4_int_t incy )
{
    return cublasCdotc( n, dx, incx, dy, incy );
}

#ifdef COMPLEX
// --------------------
/** Returns dot product (unconjugated) of vectors x and y; \f$ x^T y \f$.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in]
    dx      COMPLEX array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    dy      COMPLEX array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_tally4_cblas1
*/
extern "C"
magma_tally4FloatComplex magma_tally4_cdotu(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4FloatComplex_const_ptr dy, magma_tally4_int_t incy )
{
    return cublasCdotu( n, dx, incx, dy, incy );
}
#endif // COMPLEX

// --------------------
/** Returns 2-norm of vector x. Avoids unnecesary over/underflow.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in]
    dx      COMPLEX array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.

    @ingroup magma_tally4_cblas1
*/
extern "C" float
magma_tally4_scnrm2(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dx, magma_tally4_int_t incx )
{
    return cublasScnrm2( n, dx, incx );
}

// --------------------
/** Apply Givens plane rotation, where cos (c) is real and sin (s) is complex.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in,out]
    dx      COMPLEX array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).
            On output, overwritten with c*x + s*y.

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dy      COMPLEX array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).
            On output, overwritten with -conj(s)*x + c*y.

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @param[in]
    c       float. cosine.

    @param[in]
    s       COMPLEX. sine. c and s define a rotation
            [ c         s ]  where c*c + s*conj(s) = 1.
            [ -conj(s)  c ]

    @ingroup magma_tally4_cblas1
*/
extern "C" void
magma_tally4_crot(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dx, magma_tally4_int_t incx,
    magma_tally4FloatComplex_ptr dy, magma_tally4_int_t incy,
    float c, magma_tally4FloatComplex s )
{
    cublasCrot( n, dx, incx, dy, incy, c, s );
}

#ifdef COMPLEX
// --------------------
/** Apply Givens plane rotation, where cos (c) and sin (s) are real.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in,out]
    dx      COMPLEX array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).
            On output, overwritten with c*x + s*y.

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dy      COMPLEX array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).
            On output, overwritten with -conj(s)*x + c*y.

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @param[in]
    c       float. cosine.

    @param[in]
    s       float. sine. c and s define a rotation
            [  c  s ]  where c*c + s*s = 1.
            [ -s  c ]

    @ingroup magma_tally4_cblas1
*/
extern "C" void
magma_tally4_csrot(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dx, magma_tally4_int_t incx,
    magma_tally4FloatComplex_ptr dy, magma_tally4_int_t incy,
    float c, float s )
{
    cublasCsrot( n, dx, incx, dy, incy, c, s );
}
#endif // COMPLEX

#ifdef REAL
// --------------------
/** Apply modified plane rotation.

    @ingroup magma_tally4_cblas1
*/
extern "C" void
magma_tally4_crotm(
    magma_tally4_int_t n,
    float *dx, magma_tally4_int_t incx,
    float *dy, magma_tally4_int_t incy,
    const float *param )
{
    cublasCrotm( n, dx, incx, dy, incy, param );
}

// --------------------
/** Generate modified plane rotation.

    @ingroup magma_tally4_cblas1
*/
extern "C" void
magma_tally4_crotmg(
    float *d1, float       *d2,
    float *x1, const float *y1,
    float *param )
{
    cublasCrotmg( d1, d2, x1, y1, param );
}
#endif // REAL

// --------------------
/** Scales a vector by a constant; \f$ x = \alpha x \f$.

    @param[in]
    n       Number of elements in vector x. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in,out]
    dx      COMPLEX array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.

    @ingroup magma_tally4_cblas1
*/
extern "C" void
magma_tally4_cscal(
    magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_ptr dx, magma_tally4_int_t incx )
{
    cublasCscal( n, alpha, dx, incx );
}

#ifdef COMPLEX
// --------------------
/** Scales a vector by a real constant; \f$ x = \alpha x \f$.

    @param[in]
    n       Number of elements in vector x. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$ (real)

    @param[in,out]
    dx      COMPLEX array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.

    @ingroup magma_tally4_cblas1
*/
extern "C" void
magma_tally4_csscal(
    magma_tally4_int_t n,
    float alpha,
    magma_tally4FloatComplex_ptr dx, magma_tally4_int_t incx )
{
    cublasCsscal( n, alpha, dx, incx );
}
#endif // COMPLEX

// --------------------
/** Swap vector x and y; \f$ x <-> y \f$.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in,out]
    dx      COMPLEX array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dy      COMPLEX array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_tally4_cblas1
*/
extern "C" void
magma_tally4_cswap(
    magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dx, magma_tally4_int_t incx,
    magma_tally4FloatComplex_ptr dy, magma_tally4_int_t incy )
{
    cublasCswap( n, dx, incx, dy, incy );
}


// ========================================
// Level 2 BLAS

// --------------------
/** Perform matrix-vector product.
        \f$ y = \alpha A   x + \beta y \f$  (transA == Magma_tally4NoTrans), or \n
        \f$ y = \alpha A^T x + \beta y \f$  (transA == Magma_tally4Trans),   or \n
        \f$ y = \alpha A^H x + \beta y \f$  (transA == Magma_tally4ConjTrans).

    @param[in]
    transA  Operation to perform on A.

    @param[in]
    m       Number of rows of A. m >= 0.

    @param[in]
    n       Number of columns of A. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      COMPLEX array of dimension (ldda,n), ldda >= max(1,m).
            The m-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dx      COMPLEX array on GPU device.
            If transA == Magma_tally4NoTrans, the n element vector x of dimension (1 + (n-1)*incx); \n
            otherwise,                 the m element vector x of dimension (1 + (m-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dy      COMPLEX array on GPU device.
            If transA == Magma_tally4NoTrans, the m element vector y of dimension (1 + (m-1)*incy); \n
            otherwise,                 the n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_tally4_cblas2
*/
extern "C" void
magma_tally4_cgemv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr       dy, magma_tally4_int_t incy )
{
    cublasCgemv(
        cublas_trans_const( transA ),
        m, n,
        alpha, dA, ldda,
               dx, incx,
        beta,  dy, incy );
}

// --------------------
/** Perform rank-1 update, \f$ A = \alpha x y^H + A \f$.

    @param[in]
    m       Number of rows of A. m >= 0.

    @param[in]
    n       Number of columns of A. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dx      COMPLEX array on GPU device.
            The m element vector x of dimension (1 + (m-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    dy      COMPLEX array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @param[in,out]
    dA      COMPLEX array on GPU device.
            The m-by-n matrix A of dimension (ldda,n), ldda >= max(1,m).

    @param[in]
    ldda    Leading dimension of dA.

    @ingroup magma_tally4_cblas2
*/
extern "C" void
magma_tally4_cgerc(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4FloatComplex_const_ptr dy, magma_tally4_int_t incy,
    magma_tally4FloatComplex_ptr       dA, magma_tally4_int_t ldda )
{
    cublasCgerc(
        m, n,
        alpha, dx, incx,
               dy, incy,
               dA, ldda );
}

#ifdef COMPLEX
// --------------------
/** Perform rank-1 update (unconjugated), \f$ A = \alpha x y^H + A \f$.

    @param[in]
    m       Number of rows of A. m >= 0.

    @param[in]
    n       Number of columns of A. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dx      COMPLEX array on GPU device.
            The m element vector x of dimension (1 + (m-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    dy      COMPLEX array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @param[in,out]
    dA      COMPLEX array of dimension (ldda,n), ldda >= max(1,m).
            The m-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @ingroup magma_tally4_cblas2
*/
extern "C" void
magma_tally4_cgeru(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4FloatComplex_const_ptr dy, magma_tally4_int_t incy,
    magma_tally4FloatComplex_ptr       dA, magma_tally4_int_t ldda )
{
    cublasCgeru(
        m, n,
        alpha, dx, incx,
               dy, incy,
               dA, ldda );
}
#endif // COMPLEX

// --------------------
/** Perform Hermitian matrix-vector product, \f$ y = \alpha A x + \beta y \f$.

    @param[in]
    uplo    Whether the upper or lower triangle of A is referenced.

    @param[in]
    n       Number of rows and columns of A. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      COMPLEX array of dimension (ldda,n), ldda >= max(1,n).
            The n-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dx      COMPLEX array on GPU device.
            The m element vector x of dimension (1 + (m-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dy      COMPLEX array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_tally4_cblas2
*/
extern "C" void
magma_tally4_chemv(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr       dy, magma_tally4_int_t incy )
{
    cublasChemv(
        cublas_uplo_const( uplo ),
        n,
        alpha, dA, ldda,
               dx, incx,
        beta,  dy, incy );
}

// --------------------
/** Perform Hermitian rank-1 update, \f$ A = \alpha x x^H + A \f$.

    @param[in]
    uplo    Whether the upper or lower triangle of A is referenced. 

    @param[in]
    n       Number of rows and columns of A. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dx      COMPLEX array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dA      COMPLEX array of dimension (ldda,n), ldda >= max(1,n).
            The n-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @ingroup magma_tally4_cblas2
*/
extern "C" void
magma_tally4_cher(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    float alpha,
    magma_tally4FloatComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4FloatComplex_ptr       dA, magma_tally4_int_t ldda )
{
    cublasCher(
        cublas_uplo_const( uplo ),
        n,
        alpha, dx, incx,
               dA, ldda );
}

// --------------------
/** Perform Hermitian rank-2 update, \f$ A = \alpha x y^H + conj(\alpha) y x^H + A \f$.

    @param[in]
    uplo    Whether the upper or lower triangle of A is referenced. 

    @param[in]
    n       Number of rows and columns of A. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dx      COMPLEX array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    dy      COMPLEX array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @param[in,out]
    dA      COMPLEX array of dimension (ldda,n), ldda >= max(1,n).
            The n-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @ingroup magma_tally4_cblas2
*/
extern "C" void
magma_tally4_cher2(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dx, magma_tally4_int_t incx,
    magma_tally4FloatComplex_const_ptr dy, magma_tally4_int_t incy,
    magma_tally4FloatComplex_ptr       dA, magma_tally4_int_t ldda )
{
    cublasCher2(
        cublas_uplo_const( uplo ),
        n,
        alpha, dx, incx,
               dy, incy,
               dA, ldda );
}

// --------------------
/** Perform triangular matrix-vector product.
        \f$ x = A   x \f$  (trans == Magma_tally4NoTrans), or \n
        \f$ x = A^T x \f$  (trans == Magma_tally4Trans),   or \n
        \f$ x = A^H x \f$  (trans == Magma_tally4ConjTrans).

    @param[in]
    uplo    Whether the upper or lower triangle of A is referenced. 

    @param[in]
    trans   Operation to perform on A.

    @param[in]
    diag    Whether the diagonal of A is assumed to be unit or non-unit.

    @param[in]
    n       Number of rows and columns of A. n >= 0.

    @param[in]
    dA      COMPLEX array of dimension (ldda,n), ldda >= max(1,n).
            The n-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dx      COMPLEX array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @ingroup magma_tally4_cblas2
*/
extern "C" void
magma_tally4_ctrmv(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_diag_t diag,
    magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr       dx, magma_tally4_int_t incx )
{
    cublasCtrmv(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        cublas_diag_const( diag ),
        n,
        dA, ldda,
        dx, incx );
}

// --------------------
/** Solve triangular matrix-vector system (one right-hand side).
        \f$ A   x = b \f$  (trans == Magma_tally4NoTrans), or \n
        \f$ A^T x = b \f$  (trans == Magma_tally4Trans),   or \n
        \f$ A^H x = b \f$  (trans == Magma_tally4ConjTrans).

    @param[in]
    uplo    Whether the upper or lower triangle of A is referenced. 

    @param[in]
    trans   Operation to perform on A.

    @param[in]
    diag    Whether the diagonal of A is assumed to be unit or non-unit.

    @param[in]
    n       Number of rows and columns of A. n >= 0.

    @param[in]
    dA      COMPLEX array of dimension (ldda,n), ldda >= max(1,n).
            The n-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @param[in,out]
    dx      COMPLEX array on GPU device.
            On entry, the n element RHS vector b of dimension (1 + (n-1)*incx).
            On exit, overwritten with the solution vector x.

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @ingroup magma_tally4_cblas2
*/
extern "C" void
magma_tally4_ctrsv(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_diag_t diag,
    magma_tally4_int_t n,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr       dx, magma_tally4_int_t incx )
{
    cublasCtrsv(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        cublas_diag_const( diag ),
        n,
        dA, ldda,
        dx, incx );
}

// ========================================
// Level 3 BLAS

// --------------------
/** Perform matrix-matrix product, \f$ C = \alpha op(A) op(B) + \beta C \f$.

    @param[in]
    transA  Operation op(A) to perform on matrix A.

    @param[in]
    transB  Operation op(B) to perform on matrix B.

    @param[in]
    m       Number of rows of C and op(A). m >= 0.

    @param[in]
    n       Number of columns of C and op(B). n >= 0.

    @param[in]
    k       Number of columns of op(A) and rows of op(B). k >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      COMPLEX array on GPU device.
            If transA == Magma_tally4NoTrans, the m-by-k matrix A of dimension (ldda,k), ldda >= max(1,m); \n
            otherwise,                 the k-by-m matrix A of dimension (ldda,m), ldda >= max(1,k).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      COMPLEX array on GPU device.
            If transB == Magma_tally4NoTrans, the k-by-n matrix B of dimension (lddb,n), lddb >= max(1,k); \n
            otherwise,                 the n-by-k matrix B of dimension (lddb,k), lddb >= max(1,n).

    @param[in]
    lddb    Leading dimension of dB.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      COMPLEX array on GPU device.
            The m-by-n matrix C of dimension (lddc,n), lddc >= max(1,m).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_tally4_cblas3
*/
extern "C" void
magma_tally4_cgemm(
    magma_tally4_trans_t transA, magma_tally4_trans_t transB,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_const_ptr dB, magma_tally4_int_t lddb,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr       dC, magma_tally4_int_t lddc )
{
    cublasCgemm(
        cublas_trans_const( transA ),
        cublas_trans_const( transB ),
        m, n, k,
        alpha, dA, ldda,
               dB, lddb,
        beta,  dC, lddc );
}

// --------------------
/** Perform symmetric matrix-matrix product.
        \f$ C = \alpha A B + \beta C \f$ (side == Magma_tally4Left), or \n
        \f$ C = \alpha B A + \beta C \f$ (side == Magma_tally4Right),   \n
        where \f$ A \f$ is symmetric.

    @param[in]
    side    Whether A is on the left or right.

    @param[in]
    uplo    Whether the upper or lower triangle of A is referenced.

    @param[in]
    m       Number of rows of C. m >= 0.

    @param[in]
    n       Number of columns of C. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      COMPLEX array on GPU device.
            If side == Magma_tally4Left, the m-by-m symmetric matrix A of dimension (ldda,m), ldda >= max(1,m); \n
            otherwise,            the n-by-n symmetric matrix A of dimension (ldda,n), ldda >= max(1,n).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      COMPLEX array on GPU device.
            The m-by-n matrix B of dimension (lddb,n), lddb >= max(1,m).

    @param[in]
    lddb    Leading dimension of dB.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      COMPLEX array on GPU device.
            The m-by-n matrix C of dimension (lddc,n), lddc >= max(1,m).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_tally4_cblas3
*/
extern "C" void
magma_tally4_csymm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_const_ptr dB, magma_tally4_int_t lddb,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr       dC, magma_tally4_int_t lddc )
{
    cublasCsymm(
        cublas_side_const( side ),
        cublas_uplo_const( uplo ),
        m, n,
        alpha, dA, ldda,
               dB, lddb,
        beta,  dC, lddc );
}

// --------------------
/** Perform symmetric rank-k update.
        \f$ C = \alpha A A^T + \beta C \f$ (trans == Magma_tally4NoTrans), or \n
        \f$ C = \alpha A^T A + \beta C \f$ (trans == Magma_tally4Trans),      \n
        where \f$ C \f$ is symmetric.

    @param[in]
    uplo    Whether the upper or lower triangle of C is referenced.

    @param[in]
    trans   Operation to perform on A.

    @param[in]
    n       Number of rows and columns of C. n >= 0.

    @param[in]
    k       Number of columns of A (for Magma_tally4NoTrans) or rows of A (for Magma_tally4Trans). k >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      COMPLEX array on GPU device.
            If trans == Magma_tally4NoTrans, the n-by-k matrix A of dimension (ldda,k), ldda >= max(1,n); \n
            otherwise,                the k-by-n matrix A of dimension (ldda,n), ldda >= max(1,k).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      COMPLEX array on GPU device.
            The n-by-n symmetric matrix C of dimension (lddc,n), lddc >= max(1,n).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_tally4_cblas3
*/
extern "C" void
magma_tally4_csyrk(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr       dC, magma_tally4_int_t lddc )
{
    cublasCsyrk(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        n, k,
        alpha, dA, ldda,
        beta,  dC, lddc );
}

// --------------------
/** Perform symmetric rank-2k update.
        \f$ C = \alpha A B^T + \alpha B A^T \beta C \f$ (trans == Magma_tally4NoTrans), or \n
        \f$ C = \alpha A^T B + \alpha B^T A \beta C \f$ (trans == Magma_tally4Trans),      \n
        where \f$ C \f$ is symmetric.

    @param[in]
    uplo    Whether the upper or lower triangle of C is referenced.

    @param[in]
    trans   Operation to perform on A and B.

    @param[in]
    n       Number of rows and columns of C. n >= 0.

    @param[in]
    k       Number of columns of A and B (for Magma_tally4NoTrans) or rows of A and B (for Magma_tally4Trans). k >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      COMPLEX array on GPU device.
            If trans == Magma_tally4NoTrans, the n-by-k matrix A of dimension (ldda,k), ldda >= max(1,n); \n
            otherwise,                the k-by-n matrix A of dimension (ldda,n), ldda >= max(1,k).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      COMPLEX array on GPU device.
            If trans == Magma_tally4NoTrans, the n-by-k matrix B of dimension (lddb,k), lddb >= max(1,n); \n
            otherwise,                the k-by-n matrix B of dimension (lddb,n), lddb >= max(1,k).

    @param[in]
    lddb    Leading dimension of dB.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      COMPLEX array on GPU device.
            The n-by-n symmetric matrix C of dimension (lddc,n), lddc >= max(1,n).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_tally4_cblas3
*/
extern "C" void
magma_tally4_csyr2k(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_const_ptr dB, magma_tally4_int_t lddb,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr       dC, magma_tally4_int_t lddc )
{
    cublasCsyr2k(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        n, k,
        alpha, dA, ldda,
               dB, lddb,
        beta,  dC, lddc );
}

#ifdef COMPLEX
// --------------------
/** Perform Hermitian matrix-matrix product.
        \f$ C = \alpha A B + \beta C \f$ (side == Magma_tally4Left), or \n
        \f$ C = \alpha B A + \beta C \f$ (side == Magma_tally4Right),   \n
        where \f$ A \f$ is Hermitian.

    @param[in]
    side    Whether A is on the left or right.

    @param[in]
    uplo    Whether the upper or lower triangle of A is referenced.

    @param[in]
    m       Number of rows of C. m >= 0.

    @param[in]
    n       Number of columns of C. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      COMPLEX array on GPU device.
            If side == Magma_tally4Left, the m-by-m Hermitian matrix A of dimension (ldda,m), ldda >= max(1,m); \n
            otherwise,            the n-by-n Hermitian matrix A of dimension (ldda,n), ldda >= max(1,n).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      COMPLEX array on GPU device.
            The m-by-n matrix B of dimension (lddb,n), lddb >= max(1,m).

    @param[in]
    lddb    Leading dimension of dB.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      COMPLEX array on GPU device.
            The m-by-n matrix C of dimension (lddc,n), lddc >= max(1,m).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_tally4_cblas3
*/
extern "C" void
magma_tally4_chemm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_const_ptr dB, magma_tally4_int_t lddb,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr       dC, magma_tally4_int_t lddc )
{
    cublasChemm(
        cublas_side_const( side ),
        cublas_uplo_const( uplo ),
        m, n,
        alpha, dA, ldda,
               dB, lddb,
        beta,  dC, lddc );
}

// --------------------
/** Perform Hermitian rank-k update.
        \f$ C = \alpha A A^T + \beta C \f$ (trans == Magma_tally4NoTrans), or \n
        \f$ C = \alpha A^T A + \beta C \f$ (trans == Magma_tally4Trans),      \n
        where \f$ C \f$ is Hermitian.

    @param[in]
    uplo    Whether the upper or lower triangle of C is referenced.

    @param[in]
    trans   Operation to perform on A.

    @param[in]
    n       Number of rows and columns of C. n >= 0.

    @param[in]
    k       Number of columns of A (for Magma_tally4NoTrans) or rows of A (for Magma_tally4Trans). k >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      COMPLEX array on GPU device.
            If trans == Magma_tally4NoTrans, the n-by-k matrix A of dimension (ldda,k), ldda >= max(1,n); \n
            otherwise,                the k-by-n matrix A of dimension (ldda,n), ldda >= max(1,k).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      COMPLEX array on GPU device.
            The n-by-n Hermitian matrix C of dimension (lddc,n), lddc >= max(1,n).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_tally4_cblas3
*/
extern "C" void
magma_tally4_cherk(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    float alpha,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    float beta,
    magma_tally4FloatComplex_ptr       dC, magma_tally4_int_t lddc )
{
    cublasCherk(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        n, k,
        alpha, dA, ldda,
        beta,  dC, lddc );
}

// --------------------
/** Perform Hermitian rank-2k update.
        \f$ C = \alpha A B^T + \alpha B A^T \beta C \f$ (trans == Magma_tally4NoTrans), or \n
        \f$ C = \alpha A^T B + \alpha B^T A \beta C \f$ (trans == Magma_tally4Trans),      \n
        where \f$ C \f$ is Hermitian.

    @param[in]
    uplo    Whether the upper or lower triangle of C is referenced.

    @param[in]
    trans   Operation to perform on A and B.

    @param[in]
    n       Number of rows and columns of C. n >= 0.

    @param[in]
    k       Number of columns of A and B (for Magma_tally4NoTrans) or rows of A and B (for Magma_tally4Trans). k >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      COMPLEX array on GPU device.
            If trans == Magma_tally4NoTrans, the n-by-k matrix A of dimension (ldda,k), ldda >= max(1,n); \n
            otherwise,                the k-by-n matrix A of dimension (ldda,n), ldda >= max(1,k).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      COMPLEX array on GPU device.
            If trans == Magma_tally4NoTrans, the n-by-k matrix B of dimension (lddb,k), lddb >= max(1,n); \n
            otherwise,                the k-by-n matrix B of dimension (lddb,n), lddb >= max(1,k).

    @param[in]
    lddb    Leading dimension of dB.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      COMPLEX array on GPU device.
            The n-by-n Hermitian matrix C of dimension (lddc,n), lddc >= max(1,n).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_tally4_cblas3
*/
extern "C" void
magma_tally4_cher2k(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans,
    magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_const_ptr dB, magma_tally4_int_t lddb,
    float beta,
    magma_tally4FloatComplex_ptr       dC, magma_tally4_int_t lddc )
{
    cublasCher2k(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        n, k,
        alpha, dA, ldda,
               dB, lddb,
        beta,  dC, lddc );
}
#endif // COMPLEX

// --------------------
/** Perform triangular matrix-matrix product.
        \f$ B = \alpha op(A) B \f$ (side == Magma_tally4Left), or \n
        \f$ B = \alpha B op(A) \f$ (side == Magma_tally4Right),   \n
        where \f$ A \f$ is triangular.

    @param[in]
    side    Whether A is on the left or right.

    @param[in]
    uplo    Whether A is upper or lower triangular.

    @param[in]
    trans   Operation to perform on A.

    @param[in]
    diag    Whether the diagonal of A is assumed to be unit or non-unit.

    @param[in]
    m       Number of rows of B. m >= 0.

    @param[in]
    n       Number of columns of B. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      COMPLEX array on GPU device.
            If side == Magma_tally4Left, the n-by-n triangular matrix A of dimension (ldda,n), ldda >= max(1,n); \n
            otherwise,            the m-by-m triangular matrix A of dimension (ldda,m), ldda >= max(1,m).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      COMPLEX array on GPU device.
            The m-by-n matrix B of dimension (lddb,n), lddb >= max(1,m).

    @param[in]
    lddb    Leading dimension of dB.

    @ingroup magma_tally4_cblas3
*/
extern "C" void
magma_tally4_ctrmm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr       dB, magma_tally4_int_t lddb )
{
    cublasCtrmm(
        cublas_side_const( side ),
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        cublas_diag_const( diag ),
        m, n,
        alpha, dA, ldda,
               dB, lddb );
}

// --------------------
/** Solve triangular matrix-matrix system (multiple right-hand sides).
        \f$ op(A) X = \alpha B \f$ (side == Magma_tally4Left), or \n
        \f$ X op(A) = \alpha B \f$ (side == Magma_tally4Right),   \n
        where \f$ A \f$ is triangular.

    @param[in]
    side    Whether A is on the left or right.

    @param[in]
    uplo    Whether A is upper or lower triangular.

    @param[in]
    trans   Operation to perform on A.

    @param[in]
    diag    Whether the diagonal of A is assumed to be unit or non-unit.

    @param[in]
    m       Number of rows of B. m >= 0.

    @param[in]
    n       Number of columns of B. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      COMPLEX array on GPU device.
            If side == Magma_tally4Left, the m-by-m triangular matrix A of dimension (ldda,m), ldda >= max(1,m); \n
            otherwise,            the n-by-n triangular matrix A of dimension (ldda,n), ldda >= max(1,n).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in,out]
    dB      COMPLEX array on GPU device.
            On entry, m-by-n matrix B of dimension (lddb,n), lddb >= max(1,m).
            On exit, overwritten with the solution matrix X.

    @param[in]
    lddb    Leading dimension of dB.

    @ingroup magma_tally4_cblas3
*/
extern "C" void
magma_tally4_ctrsm(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_const_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr       dB, magma_tally4_int_t lddb )
{
    cublasCtrsm(
        cublas_side_const( side ),
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        cublas_diag_const( diag ),
        m, n,
        alpha, dA, ldda,
               dB, lddb );
}

#endif // HAVE_CUBLAS

#undef COMPLEX
