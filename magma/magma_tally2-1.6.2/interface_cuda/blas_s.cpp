/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
 
       @author Mark Gates
       @generated from blas_z.cpp normal z -> s, Fri Jan 30 19:00:20 2015
*/

#include <stdlib.h>
#include <stdio.h>

#include "magma_tally2.h"
#include "error.h"

#define REAL

#ifdef HAVE_CUBLAS

// ----------------------------------------
// Convert MAGMA_tally2 constants to CUBLAS v1 constants, which are the same as lapack.
// These must be static to avoid conflict with CUBLAS v2 translators.
#define cublas_trans_const_tally2( magma_tally2_const ) lapacke_trans_const_tally2( magma_tally2_const )
#define cublas_side_const_tally2(  magma_tally2_const ) lapacke_side_const_tally2(  magma_tally2_const )
#define cublas_diag_const_tally2(  magma_tally2_const ) lapacke_diag_const_tally2(  magma_tally2_const )
#define cublas_uplo_const_tally2(  magma_tally2_const ) lapacke_uplo_const_tally2(  magma_tally2_const )


// ========================================
// Level 1 BLAS

// --------------------
/** Returns index of element of vector x having max. absolute value;
    i.e., max (infinity) norm.
    
    @param[in]
    n       Number of elements in vector x. n >= 0.
            
    @param[in]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).
            
    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.
            
    @ingroup magma_tally2_sblas1
*/
extern "C" magma_tally2_int_t
magma_tally2_isamax(
    magma_tally2_int_t n,
    magma_tally2Float_const_ptr dx, magma_tally2_int_t incx )
{
    return cublasIsamax( n, dx, incx );
}

// --------------------
/** Returns index of element of vector x having min. absolute value.

    @param[in]
    n       Number of elements in vector x. n >= 0.

    @param[in]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.

    @ingroup magma_tally2_sblas1
*/
extern "C" magma_tally2_int_t
magma_tally2_isamin(
    magma_tally2_int_t n,
    magma_tally2Float_const_ptr dx, magma_tally2_int_t incx )
{
    return cublasIsamin( n, dx, incx );
}

// --------------------
/** Returns the sum of absolute values of vector x; i.e., one norm.

    @param[in]
    n       Number of elements in vector x. n >= 0.

    @param[in]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.

    @ingroup magma_tally2_sblas1
*/
extern "C" float
magma_tally2_sasum(
    magma_tally2_int_t n,
    magma_tally2Float_const_ptr dx, magma_tally2_int_t incx )
{
    return cublasSasum( n, dx, incx );
}

// --------------------
/** Constant times a vector plus a vector; \f$ y = \alpha x + y \f$.

    @param[in]
    n       Number of elements in vectors x and y. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dy      REAL array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_tally2_sblas1
*/
extern "C" void
magma_tally2_saxpy(
    magma_tally2_int_t n,
    float alpha,
    magma_tally2Float_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2Float_ptr       dy, magma_tally2_int_t incy )
{
    cublasSaxpy( n, alpha, dx, incx, dy, incy );
}

// --------------------
/** Copy vector x to vector y; \f$ y = x \f$.

    @param[in]
    n       Number of elements in vectors x and y. n >= 0.

    @param[in]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[out]
    dy      REAL array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_tally2_sblas1
*/
extern "C" void
magma_tally2_scopy(
    magma_tally2_int_t n,
    magma_tally2Float_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2Float_ptr       dy, magma_tally2_int_t incy )
{
    cublasScopy( n, dx, incx, dy, incy );
}

// --------------------
/** Returns dot product of vectors x and y; \f$ x^H y \f$.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    dy      REAL array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_tally2_sblas1
*/
extern "C"
float magma_tally2_sdot(
    magma_tally2_int_t n,
    magma_tally2Float_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2Float_const_ptr dy, magma_tally2_int_t incy )
{
    return cublasSdot( n, dx, incx, dy, incy );
}

#ifdef COMPLEX
// --------------------
/** Returns dot product (unconjugated) of vectors x and y; \f$ x^T y \f$.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    dy      REAL array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_tally2_sblas1
*/
extern "C"
float magma_tally2_sdot(
    magma_tally2_int_t n,
    magma_tally2Float_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2Float_const_ptr dy, magma_tally2_int_t incy )
{
    return cublasSdot( n, dx, incx, dy, incy );
}
#endif // COMPLEX

// --------------------
/** Returns 2-norm of vector x. Avoids unnecesary over/underflow.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.

    @ingroup magma_tally2_sblas1
*/
extern "C" float
magma_tally2_snrm2(
    magma_tally2_int_t n,
    magma_tally2Float_const_ptr dx, magma_tally2_int_t incx )
{
    return cublasSnrm2( n, dx, incx );
}

// --------------------
/** Apply Givens plane rotation, where cos (c) is real and sin (s) is real.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in,out]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).
            On output, overwritten with c*x + s*y.

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dy      REAL array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).
            On output, overwritten with -conj(s)*x + c*y.

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @param[in]
    c       float. cosine.

    @param[in]
    s       REAL. sine. c and s define a rotation
            [ c         s ]  where c*c + s*conj(s) = 1.
            [ -conj(s)  c ]

    @ingroup magma_tally2_sblas1
*/
extern "C" void
magma_tally2_srot(
    magma_tally2_int_t n,
    magma_tally2Float_ptr dx, magma_tally2_int_t incx,
    magma_tally2Float_ptr dy, magma_tally2_int_t incy,
    float c, float s )
{
    cublasSrot( n, dx, incx, dy, incy, c, s );
}

#ifdef COMPLEX
// --------------------
/** Apply Givens plane rotation, where cos (c) and sin (s) are real.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in,out]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).
            On output, overwritten with c*x + s*y.

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dy      REAL array on GPU device.
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

    @ingroup magma_tally2_sblas1
*/
extern "C" void
magma_tally2_srot(
    magma_tally2_int_t n,
    magma_tally2Float_ptr dx, magma_tally2_int_t incx,
    magma_tally2Float_ptr dy, magma_tally2_int_t incy,
    float c, float s )
{
    cublasSrot( n, dx, incx, dy, incy, c, s );
}
#endif // COMPLEX

#ifdef REAL
// --------------------
/** Apply modified plane rotation.

    @ingroup magma_tally2_sblas1
*/
extern "C" void
magma_tally2_srotm(
    magma_tally2_int_t n,
    float *dx, magma_tally2_int_t incx,
    float *dy, magma_tally2_int_t incy,
    const float *param )
{
    cublasSrotm( n, dx, incx, dy, incy, param );
}

// --------------------
/** Generate modified plane rotation.

    @ingroup magma_tally2_sblas1
*/
extern "C" void
magma_tally2_srotmg(
    float *d1, float       *d2,
    float *x1, const float *y1,
    float *param )
{
    cublasSrotmg( d1, d2, x1, y1, param );
}
#endif // REAL

// --------------------
/** Scales a vector by a constant; \f$ x = \alpha x \f$.

    @param[in]
    n       Number of elements in vector x. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in,out]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.

    @ingroup magma_tally2_sblas1
*/
extern "C" void
magma_tally2_sscal(
    magma_tally2_int_t n,
    float alpha,
    magma_tally2Float_ptr dx, magma_tally2_int_t incx )
{
    cublasSscal( n, alpha, dx, incx );
}

#ifdef COMPLEX
// --------------------
/** Scales a vector by a real constant; \f$ x = \alpha x \f$.

    @param[in]
    n       Number of elements in vector x. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$ (real)

    @param[in,out]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.

    @ingroup magma_tally2_sblas1
*/
extern "C" void
magma_tally2_sscal(
    magma_tally2_int_t n,
    float alpha,
    magma_tally2Float_ptr dx, magma_tally2_int_t incx )
{
    cublasSscal( n, alpha, dx, incx );
}
#endif // COMPLEX

// --------------------
/** Swap vector x and y; \f$ x <-> y \f$.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in,out]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dy      REAL array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_tally2_sblas1
*/
extern "C" void
magma_tally2_sswap(
    magma_tally2_int_t n,
    magma_tally2Float_ptr dx, magma_tally2_int_t incx,
    magma_tally2Float_ptr dy, magma_tally2_int_t incy )
{
    cublasSswap( n, dx, incx, dy, incy );
}


// ========================================
// Level 2 BLAS

// --------------------
/** Perform matrix-vector product.
        \f$ y = \alpha A   x + \beta y \f$  (transA == Magma_tally2NoTrans), or \n
        \f$ y = \alpha A^T x + \beta y \f$  (transA == Magma_tally2Trans),   or \n
        \f$ y = \alpha A^H x + \beta y \f$  (transA == Magma_tally2ConjTrans).

    @param[in]
    transA  Operation to perform on A.

    @param[in]
    m       Number of rows of A. m >= 0.

    @param[in]
    n       Number of columns of A. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      REAL array of dimension (ldda,n), ldda >= max(1,m).
            The m-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dx      REAL array on GPU device.
            If transA == Magma_tally2NoTrans, the n element vector x of dimension (1 + (n-1)*incx); \n
            otherwise,                 the m element vector x of dimension (1 + (m-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dy      REAL array on GPU device.
            If transA == Magma_tally2NoTrans, the m element vector y of dimension (1 + (m-1)*incy); \n
            otherwise,                 the n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_tally2_sblas2
*/
extern "C" void
magma_tally2_sgemv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    float alpha,
    magma_tally2Float_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_const_ptr dx, magma_tally2_int_t incx,
    float beta,
    magma_tally2Float_ptr       dy, magma_tally2_int_t incy )
{
    cublasSgemv(
        cublas_trans_const_tally2( transA ),
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
    dx      REAL array on GPU device.
            The m element vector x of dimension (1 + (m-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    dy      REAL array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @param[in,out]
    dA      REAL array on GPU device.
            The m-by-n matrix A of dimension (ldda,n), ldda >= max(1,m).

    @param[in]
    ldda    Leading dimension of dA.

    @ingroup magma_tally2_sblas2
*/
extern "C" void
magma_tally2_sger(
    magma_tally2_int_t m, magma_tally2_int_t n,
    float alpha,
    magma_tally2Float_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2Float_const_ptr dy, magma_tally2_int_t incy,
    magma_tally2Float_ptr       dA, magma_tally2_int_t ldda )
{
    cublasSger(
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
    dx      REAL array on GPU device.
            The m element vector x of dimension (1 + (m-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    dy      REAL array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @param[in,out]
    dA      REAL array of dimension (ldda,n), ldda >= max(1,m).
            The m-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @ingroup magma_tally2_sblas2
*/
extern "C" void
magma_tally2_sger(
    magma_tally2_int_t m, magma_tally2_int_t n,
    float alpha,
    magma_tally2Float_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2Float_const_ptr dy, magma_tally2_int_t incy,
    magma_tally2Float_ptr       dA, magma_tally2_int_t ldda )
{
    cublasSger(
        m, n,
        alpha, dx, incx,
               dy, incy,
               dA, ldda );
}
#endif // COMPLEX

// --------------------
/** Perform symmetric matrix-vector product, \f$ y = \alpha A x + \beta y \f$.

    @param[in]
    uplo    Whether the upper or lower triangle of A is referenced.

    @param[in]
    n       Number of rows and columns of A. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      REAL array of dimension (ldda,n), ldda >= max(1,n).
            The n-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dx      REAL array on GPU device.
            The m element vector x of dimension (1 + (m-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dy      REAL array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_tally2_sblas2
*/
extern "C" void
magma_tally2_ssymv(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    float alpha,
    magma_tally2Float_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_const_ptr dx, magma_tally2_int_t incx,
    float beta,
    magma_tally2Float_ptr       dy, magma_tally2_int_t incy )
{
    cublasSsymv(
        cublas_uplo_const_tally2( uplo ),
        n,
        alpha, dA, ldda,
               dx, incx,
        beta,  dy, incy );
}

// --------------------
/** Perform symmetric rank-1 update, \f$ A = \alpha x x^H + A \f$.

    @param[in]
    uplo    Whether the upper or lower triangle of A is referenced. 

    @param[in]
    n       Number of rows and columns of A. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dA      REAL array of dimension (ldda,n), ldda >= max(1,n).
            The n-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @ingroup magma_tally2_sblas2
*/
extern "C" void
magma_tally2_ssyr(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    float alpha,
    magma_tally2Float_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2Float_ptr       dA, magma_tally2_int_t ldda )
{
    cublasSsyr(
        cublas_uplo_const_tally2( uplo ),
        n,
        alpha, dx, incx,
               dA, ldda );
}

// --------------------
/** Perform symmetric rank-2 update, \f$ A = \alpha x y^H + conj(\alpha) y x^H + A \f$.

    @param[in]
    uplo    Whether the upper or lower triangle of A is referenced. 

    @param[in]
    n       Number of rows and columns of A. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    dy      REAL array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @param[in,out]
    dA      REAL array of dimension (ldda,n), ldda >= max(1,n).
            The n-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @ingroup magma_tally2_sblas2
*/
extern "C" void
magma_tally2_ssyr2(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t n,
    float alpha,
    magma_tally2Float_const_ptr dx, magma_tally2_int_t incx,
    magma_tally2Float_const_ptr dy, magma_tally2_int_t incy,
    magma_tally2Float_ptr       dA, magma_tally2_int_t ldda )
{
    cublasSsyr2(
        cublas_uplo_const_tally2( uplo ),
        n,
        alpha, dx, incx,
               dy, incy,
               dA, ldda );
}

// --------------------
/** Perform triangular matrix-vector product.
        \f$ x = A   x \f$  (trans == Magma_tally2NoTrans), or \n
        \f$ x = A^T x \f$  (trans == Magma_tally2Trans),   or \n
        \f$ x = A^H x \f$  (trans == Magma_tally2ConjTrans).

    @param[in]
    uplo    Whether the upper or lower triangle of A is referenced. 

    @param[in]
    trans   Operation to perform on A.

    @param[in]
    diag    Whether the diagonal of A is assumed to be unit or non-unit.

    @param[in]
    n       Number of rows and columns of A. n >= 0.

    @param[in]
    dA      REAL array of dimension (ldda,n), ldda >= max(1,n).
            The n-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dx      REAL array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @ingroup magma_tally2_sblas2
*/
extern "C" void
magma_tally2_strmv(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_diag_t diag,
    magma_tally2_int_t n,
    magma_tally2Float_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr       dx, magma_tally2_int_t incx )
{
    cublasStrmv(
        cublas_uplo_const_tally2( uplo ),
        cublas_trans_const_tally2( trans ),
        cublas_diag_const_tally2( diag ),
        n,
        dA, ldda,
        dx, incx );
}

// --------------------
/** Solve triangular matrix-vector system (one right-hand side).
        \f$ A   x = b \f$  (trans == Magma_tally2NoTrans), or \n
        \f$ A^T x = b \f$  (trans == Magma_tally2Trans),   or \n
        \f$ A^H x = b \f$  (trans == Magma_tally2ConjTrans).

    @param[in]
    uplo    Whether the upper or lower triangle of A is referenced. 

    @param[in]
    trans   Operation to perform on A.

    @param[in]
    diag    Whether the diagonal of A is assumed to be unit or non-unit.

    @param[in]
    n       Number of rows and columns of A. n >= 0.

    @param[in]
    dA      REAL array of dimension (ldda,n), ldda >= max(1,n).
            The n-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @param[in,out]
    dx      REAL array on GPU device.
            On entry, the n element RHS vector b of dimension (1 + (n-1)*incx).
            On exit, overwritten with the solution vector x.

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @ingroup magma_tally2_sblas2
*/
extern "C" void
magma_tally2_strsv(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_diag_t diag,
    magma_tally2_int_t n,
    magma_tally2Float_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr       dx, magma_tally2_int_t incx )
{
    cublasStrsv(
        cublas_uplo_const_tally2( uplo ),
        cublas_trans_const_tally2( trans ),
        cublas_diag_const_tally2( diag ),
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
    dA      REAL array on GPU device.
            If transA == Magma_tally2NoTrans, the m-by-k matrix A of dimension (ldda,k), ldda >= max(1,m); \n
            otherwise,                 the k-by-m matrix A of dimension (ldda,m), ldda >= max(1,k).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      REAL array on GPU device.
            If transB == Magma_tally2NoTrans, the k-by-n matrix B of dimension (lddb,n), lddb >= max(1,k); \n
            otherwise,                 the n-by-k matrix B of dimension (lddb,k), lddb >= max(1,n).

    @param[in]
    lddb    Leading dimension of dB.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      REAL array on GPU device.
            The m-by-n matrix C of dimension (lddc,n), lddc >= max(1,m).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_tally2_sblas3
*/
extern "C" void
magma_tally2_sgemm(
    magma_tally2_trans_t transA, magma_tally2_trans_t transB,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    float alpha,
    magma_tally2Float_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_const_ptr dB, magma_tally2_int_t lddb,
    float beta,
    magma_tally2Float_ptr       dC, magma_tally2_int_t lddc )
{
    cublasSgemm(
        cublas_trans_const_tally2( transA ),
        cublas_trans_const_tally2( transB ),
        m, n, k,
        alpha, dA, ldda,
               dB, lddb,
        beta,  dC, lddc );
}

// --------------------
/** Perform symmetric matrix-matrix product.
        \f$ C = \alpha A B + \beta C \f$ (side == Magma_tally2Left), or \n
        \f$ C = \alpha B A + \beta C \f$ (side == Magma_tally2Right),   \n
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
    dA      REAL array on GPU device.
            If side == Magma_tally2Left, the m-by-m symmetric matrix A of dimension (ldda,m), ldda >= max(1,m); \n
            otherwise,            the n-by-n symmetric matrix A of dimension (ldda,n), ldda >= max(1,n).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      REAL array on GPU device.
            The m-by-n matrix B of dimension (lddb,n), lddb >= max(1,m).

    @param[in]
    lddb    Leading dimension of dB.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      REAL array on GPU device.
            The m-by-n matrix C of dimension (lddc,n), lddc >= max(1,m).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_tally2_sblas3
*/
extern "C" void
magma_tally2_ssymm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo,
    magma_tally2_int_t m, magma_tally2_int_t n,
    float alpha,
    magma_tally2Float_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_const_ptr dB, magma_tally2_int_t lddb,
    float beta,
    magma_tally2Float_ptr       dC, magma_tally2_int_t lddc )
{
    cublasSsymm(
        cublas_side_const_tally2( side ),
        cublas_uplo_const_tally2( uplo ),
        m, n,
        alpha, dA, ldda,
               dB, lddb,
        beta,  dC, lddc );
}

// --------------------
/** Perform symmetric rank-k update.
        \f$ C = \alpha A A^T + \beta C \f$ (trans == Magma_tally2NoTrans), or \n
        \f$ C = \alpha A^T A + \beta C \f$ (trans == Magma_tally2Trans),      \n
        where \f$ C \f$ is symmetric.

    @param[in]
    uplo    Whether the upper or lower triangle of C is referenced.

    @param[in]
    trans   Operation to perform on A.

    @param[in]
    n       Number of rows and columns of C. n >= 0.

    @param[in]
    k       Number of columns of A (for Magma_tally2NoTrans) or rows of A (for Magma_tally2Trans). k >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      REAL array on GPU device.
            If trans == Magma_tally2NoTrans, the n-by-k matrix A of dimension (ldda,k), ldda >= max(1,n); \n
            otherwise,                the k-by-n matrix A of dimension (ldda,n), ldda >= max(1,k).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      REAL array on GPU device.
            The n-by-n symmetric matrix C of dimension (lddc,n), lddc >= max(1,n).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_tally2_sblas3
*/
extern "C" void
magma_tally2_ssyrk(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    float alpha,
    magma_tally2Float_const_ptr dA, magma_tally2_int_t ldda,
    float beta,
    magma_tally2Float_ptr       dC, magma_tally2_int_t lddc )
{
    cublasSsyrk(
        cublas_uplo_const_tally2( uplo ),
        cublas_trans_const_tally2( trans ),
        n, k,
        alpha, dA, ldda,
        beta,  dC, lddc );
}

// --------------------
/** Perform symmetric rank-2k update.
        \f$ C = \alpha A B^T + \alpha B A^T \beta C \f$ (trans == Magma_tally2NoTrans), or \n
        \f$ C = \alpha A^T B + \alpha B^T A \beta C \f$ (trans == Magma_tally2Trans),      \n
        where \f$ C \f$ is symmetric.

    @param[in]
    uplo    Whether the upper or lower triangle of C is referenced.

    @param[in]
    trans   Operation to perform on A and B.

    @param[in]
    n       Number of rows and columns of C. n >= 0.

    @param[in]
    k       Number of columns of A and B (for Magma_tally2NoTrans) or rows of A and B (for Magma_tally2Trans). k >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      REAL array on GPU device.
            If trans == Magma_tally2NoTrans, the n-by-k matrix A of dimension (ldda,k), ldda >= max(1,n); \n
            otherwise,                the k-by-n matrix A of dimension (ldda,n), ldda >= max(1,k).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      REAL array on GPU device.
            If trans == Magma_tally2NoTrans, the n-by-k matrix B of dimension (lddb,k), lddb >= max(1,n); \n
            otherwise,                the k-by-n matrix B of dimension (lddb,n), lddb >= max(1,k).

    @param[in]
    lddb    Leading dimension of dB.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      REAL array on GPU device.
            The n-by-n symmetric matrix C of dimension (lddc,n), lddc >= max(1,n).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_tally2_sblas3
*/
extern "C" void
magma_tally2_ssyr2k(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    float alpha,
    magma_tally2Float_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_const_ptr dB, magma_tally2_int_t lddb,
    float beta,
    magma_tally2Float_ptr       dC, magma_tally2_int_t lddc )
{
    cublasSsyr2k(
        cublas_uplo_const_tally2( uplo ),
        cublas_trans_const_tally2( trans ),
        n, k,
        alpha, dA, ldda,
               dB, lddb,
        beta,  dC, lddc );
}

#ifdef COMPLEX
// --------------------
/** Perform symmetric matrix-matrix product.
        \f$ C = \alpha A B + \beta C \f$ (side == Magma_tally2Left), or \n
        \f$ C = \alpha B A + \beta C \f$ (side == Magma_tally2Right),   \n
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
    dA      REAL array on GPU device.
            If side == Magma_tally2Left, the m-by-m symmetric matrix A of dimension (ldda,m), ldda >= max(1,m); \n
            otherwise,            the n-by-n symmetric matrix A of dimension (ldda,n), ldda >= max(1,n).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      REAL array on GPU device.
            The m-by-n matrix B of dimension (lddb,n), lddb >= max(1,m).

    @param[in]
    lddb    Leading dimension of dB.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      REAL array on GPU device.
            The m-by-n matrix C of dimension (lddc,n), lddc >= max(1,m).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_tally2_sblas3
*/
extern "C" void
magma_tally2_ssymm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo,
    magma_tally2_int_t m, magma_tally2_int_t n,
    float alpha,
    magma_tally2Float_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_const_ptr dB, magma_tally2_int_t lddb,
    float beta,
    magma_tally2Float_ptr       dC, magma_tally2_int_t lddc )
{
    cublasSsymm(
        cublas_side_const_tally2( side ),
        cublas_uplo_const_tally2( uplo ),
        m, n,
        alpha, dA, ldda,
               dB, lddb,
        beta,  dC, lddc );
}

// --------------------
/** Perform symmetric rank-k update.
        \f$ C = \alpha A A^T + \beta C \f$ (trans == Magma_tally2NoTrans), or \n
        \f$ C = \alpha A^T A + \beta C \f$ (trans == Magma_tally2Trans),      \n
        where \f$ C \f$ is symmetric.

    @param[in]
    uplo    Whether the upper or lower triangle of C is referenced.

    @param[in]
    trans   Operation to perform on A.

    @param[in]
    n       Number of rows and columns of C. n >= 0.

    @param[in]
    k       Number of columns of A (for Magma_tally2NoTrans) or rows of A (for Magma_tally2Trans). k >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      REAL array on GPU device.
            If trans == Magma_tally2NoTrans, the n-by-k matrix A of dimension (ldda,k), ldda >= max(1,n); \n
            otherwise,                the k-by-n matrix A of dimension (ldda,n), ldda >= max(1,k).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      REAL array on GPU device.
            The n-by-n symmetric matrix C of dimension (lddc,n), lddc >= max(1,n).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_tally2_sblas3
*/
extern "C" void
magma_tally2_ssyrk(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    float alpha,
    magma_tally2Float_const_ptr dA, magma_tally2_int_t ldda,
    float beta,
    magma_tally2Float_ptr       dC, magma_tally2_int_t lddc )
{
    cublasSsyrk(
        cublas_uplo_const_tally2( uplo ),
        cublas_trans_const_tally2( trans ),
        n, k,
        alpha, dA, ldda,
        beta,  dC, lddc );
}

// --------------------
/** Perform symmetric rank-2k update.
        \f$ C = \alpha A B^T + \alpha B A^T \beta C \f$ (trans == Magma_tally2NoTrans), or \n
        \f$ C = \alpha A^T B + \alpha B^T A \beta C \f$ (trans == Magma_tally2Trans),      \n
        where \f$ C \f$ is symmetric.

    @param[in]
    uplo    Whether the upper or lower triangle of C is referenced.

    @param[in]
    trans   Operation to perform on A and B.

    @param[in]
    n       Number of rows and columns of C. n >= 0.

    @param[in]
    k       Number of columns of A and B (for Magma_tally2NoTrans) or rows of A and B (for Magma_tally2Trans). k >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      REAL array on GPU device.
            If trans == Magma_tally2NoTrans, the n-by-k matrix A of dimension (ldda,k), ldda >= max(1,n); \n
            otherwise,                the k-by-n matrix A of dimension (ldda,n), ldda >= max(1,k).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      REAL array on GPU device.
            If trans == Magma_tally2NoTrans, the n-by-k matrix B of dimension (lddb,k), lddb >= max(1,n); \n
            otherwise,                the k-by-n matrix B of dimension (lddb,n), lddb >= max(1,k).

    @param[in]
    lddb    Leading dimension of dB.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      REAL array on GPU device.
            The n-by-n symmetric matrix C of dimension (lddc,n), lddc >= max(1,n).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_tally2_sblas3
*/
extern "C" void
magma_tally2_ssyr2k(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans,
    magma_tally2_int_t n, magma_tally2_int_t k,
    float alpha,
    magma_tally2Float_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_const_ptr dB, magma_tally2_int_t lddb,
    float beta,
    magma_tally2Float_ptr       dC, magma_tally2_int_t lddc )
{
    cublasSsyr2k(
        cublas_uplo_const_tally2( uplo ),
        cublas_trans_const_tally2( trans ),
        n, k,
        alpha, dA, ldda,
               dB, lddb,
        beta,  dC, lddc );
}
#endif // COMPLEX

// --------------------
/** Perform triangular matrix-matrix product.
        \f$ B = \alpha op(A) B \f$ (side == Magma_tally2Left), or \n
        \f$ B = \alpha B op(A) \f$ (side == Magma_tally2Right),   \n
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
    dA      REAL array on GPU device.
            If side == Magma_tally2Left, the n-by-n triangular matrix A of dimension (ldda,n), ldda >= max(1,n); \n
            otherwise,            the m-by-m triangular matrix A of dimension (ldda,m), ldda >= max(1,m).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      REAL array on GPU device.
            The m-by-n matrix B of dimension (lddb,n), lddb >= max(1,m).

    @param[in]
    lddb    Leading dimension of dB.

    @ingroup magma_tally2_sblas3
*/
extern "C" void
magma_tally2_strmm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n,
    float alpha,
    magma_tally2Float_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr       dB, magma_tally2_int_t lddb )
{
    cublasStrmm(
        cublas_side_const_tally2( side ),
        cublas_uplo_const_tally2( uplo ),
        cublas_trans_const_tally2( trans ),
        cublas_diag_const_tally2( diag ),
        m, n,
        alpha, dA, ldda,
               dB, lddb );
}

// --------------------
/** Solve triangular matrix-matrix system (multiple right-hand sides).
        \f$ op(A) X = \alpha B \f$ (side == Magma_tally2Left), or \n
        \f$ X op(A) = \alpha B \f$ (side == Magma_tally2Right),   \n
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
    dA      REAL array on GPU device.
            If side == Magma_tally2Left, the m-by-m triangular matrix A of dimension (ldda,m), ldda >= max(1,m); \n
            otherwise,            the n-by-n triangular matrix A of dimension (ldda,n), ldda >= max(1,n).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in,out]
    dB      REAL array on GPU device.
            On entry, m-by-n matrix B of dimension (lddb,n), lddb >= max(1,m).
            On exit, overwritten with the solution matrix X.

    @param[in]
    lddb    Leading dimension of dB.

    @ingroup magma_tally2_sblas3
*/
extern "C" void
magma_tally2_strsm(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n,
    float alpha,
    magma_tally2Float_const_ptr dA, magma_tally2_int_t ldda,
    magma_tally2Float_ptr       dB, magma_tally2_int_t lddb )
{
    cublasStrsm(
        cublas_side_const_tally2( side ),
        cublas_uplo_const_tally2( uplo ),
        cublas_trans_const_tally2( trans ),
        cublas_diag_const_tally2( diag ),
        m, n,
        alpha, dA, ldda,
               dB, lddb );
}

#endif // HAVE_CUBLAS

#undef REAL
